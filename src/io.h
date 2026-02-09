/*
  IO interface for the CXX4 netCDF routines
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
  Dependencies: netcdf library with cxx4 bindings (not cxx!)
*/
#ifndef IO_H
#define IO_H

#include <iostream>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <string>
#include "cmemt.h"

// Some functions outside the class
std::string file_exists(const std::string& name);
bool        bfile_exists(const std::string& name);

//
class io{
 private:
  netCDF::NcFile *ifile;
  std::string file;
 public:
  //
  std::vector<netCDF::NcDim> dims;
  std::vector<netCDF::NcVar> vars;
  // Constructor
  io(std::string filename, netCDF::NcFile::FileMode mode = netCDF::NcFile::write, bool verbose = true){
    initRead(filename, mode, verbose);
  }
 io():ifile(NULL){};
 void sync();
 
  // Destructor
  ~io(){
    dims.clear();
    vars.clear();
    delete ifile;
  }
  
  
  ////////////////////////////////////////
  // Methods implemented in the .cc file//
  ////////////////////////////////////////
  void initRead(std::string filename, netCDF::NcFile::FileMode mode = netCDF::NcFile::write, bool verbose = true);
  void initDim(std::string vname, int size = 0);
  void initDim(std::vector<std::string> vname, std::vector<int> vsize);
  void varAttr(std::string vname, std::string attr_n, std::string attr_v);

  std::vector<int> dimSize(std::string vname);

  //void initVar(std::string vname, std::vector<std::string> dnames){

  

  //////////////
  // TEMPLATES//
  //////////////

  // bool read_Tstep(std::string vname, mat<double> &res,int irec = 0);
  // Read variable assuming double
  template <class T> bool read_Tstep(std::string vname, mat<T> &res, int irec = 0, bool verbose = true){
  
    std::string inam = "io::read_Tstep: ";
    res.d.clear();
    
    //
    // Check which element corresponds to vname
    //
    bool found = false;
    int idx = 0;
    for(int ii = 0; ii<(int)vars.size(); ii++){

      if(vars[ii].getName().compare(vname) == 0){
	found = true;
	idx = ii;
	break;
      } 
    }
    if(!found){
      std::cerr << inam << "WARNING, variable "<<vname<<" not found in "<< file<<std::endl;
      return false;
    }
  
  
    //
    // Check dimensions, any unlimited
    // 
    int ndim = vars[idx].getDimCount();
    std::vector<netCDF::NcDim> vdims = vars[idx].getDims();
    //bool unlim[ndim];
    //int k = 0;
    int nun=0;
  
    std::vector<int> newdims;

    for(auto &it: vdims){
      //unlim[k++] = it.isUnlimited();
      if(!it.isUnlimited()) {
	newdims.push_back(it.getSize());
	nun++;
      }
    }

    //
    // allocate res with only the limited variables
    //
    // std::cout << newdims.size()<<std::endl;

    res.set(newdims);
    // std::cout << res.size(0)<<std::endl;

    if((ndim-newdims.size()) == 1){ // There is one unlimited dimention (time)
      std::vector<size_t> start, count;
      //
      for(auto &it: vdims){
	if(it.isUnlimited()){
	  irec = std::min(int(it.getSize()-1), irec);
	  start.push_back(irec);
	  count.push_back(1);
	}else{
	  start.push_back(0);
	  count.push_back(it.getSize());
	}
      }
      vars[idx].getVar(start, count, &res.d[0]);
      if(res.isNaN()){
	std::cerr << inam << "ERROR, variable ["<<vname<<"] contains NaNs, exiting!"<<std::endl;
	exit(0);
      }
      
      
      //
      if(verbose){
	std::cout << inam <<"read "<<vname<<" (t="<<irec<<") ["<<res.size(0);
	for(int tt=1;tt<res.ndims();tt++) std::cout << ", "<<res.size(tt);
	std::cout <<"]"<<std::endl;
      }

    }else if((ndim-newdims.size()) == 0){ // All dimensions are limited -> read all
      vars[idx].getVar(&res.d[0]);
      //
      if(verbose){
	std::cout << inam <<"read "<<vname<<" ["<<res.size(0);
	for(int tt=1;tt<res.ndims();tt++) std::cout << ", "<<res.size(tt);
	std::cout <<"]"<<std::endl;
      }
    }

    return true;
  }


  
  template <class T> bool write_Tstep(std::string vname, mat<T> &var, int irec = 0){
    std::string inam = "io::write_Tstep: ";
    
    // Check if the variable exists and copy it
    netCDF::NcVar ivar;
    //
    bool exists = false;
    for(auto &it: vars){
      if(it.getName().compare(vname) == 0){
	ivar = it;
	exists = true;
      }
    }
    if(!exists){
      std::cerr << inam <<"ERROR, "<<vname <<" has not been initialized before"<<std::endl;
      exit(0);
    }

    // Extract dimensions of the variable and check how many are unlimited
    std::vector<netCDF::NcDim> idims =  ivar.getDims();
    int nun = 0 ;
    for(auto &it: idims) if(it.isUnlimited()) nun++;
  


    // Check if it is a time var or not
    if(nun == 1){
      std::vector<size_t> start, count;
      for(auto &it: idims){
	if(it.isUnlimited()){
	  start.push_back(irec);
	  count.push_back(1);
	}else{
	  start.push_back(0);
	  count.push_back(it.getSize());
	}
      }
      ivar.putVar(start, count, &var.d[0]);
      std::cerr << inam << "writing ["<<vname << "] (t="<< irec <<") to "<<file<<std::endl;

    } else if(nun == 0){
      var.fillNaN(0.0);
      ivar.putVar(&var.d[0]);
      std::cerr << inam << "writing ["<<vname << "] to "<<file<<std::endl;

    } else{
      std::cerr << inam <<"ERROR, more than one time dimension is not implemented"<<std::endl;
      exit(0);
    }
    
    return true;
  }

// Storing only a subset of pixels for a given time step to avoid slow checkpoints and avoid rewriting the full grid.
//
// The subset is provided as a sequential pixel range [seq0, seq1] (inclusive),
// where seq = iy * nx + ix (row-major order).
//
// Assumptions:
// - The first two *non-unlimited* dimensions of the netCDF variable are (y, x).
// - `var` contains the FULL field (y, x, ...), not just the subset.
// - Supports var.ndims() in [2..4] (e.g., (y,x), (y,x,k), (y,x,k,l)).
template <class T> bool write_Tstep_pixels(std::string vname, mat<T> &var, int irec,
                        long seq0, long seq1, int nx){

  std::string inam = "io::write_Tstep_pixels: ";

  if(seq1 < seq0) return true; // nothing to do
  if(nx <= 0){
    std::cerr << inam << "ERROR, nx must be > 0" << std::endl;
    exit(0);
  }

  // Check if the variable exists
  netCDF::NcVar ivar;
  bool exists = false;
  for(auto &it: vars){
    if(it.getName().compare(vname) == 0){
      ivar = it;
      exists = true;
      break;
    }
  }
  if(!exists){
    std::cerr << inam << "ERROR, "<< vname <<" has not been initialized before" << std::endl;
    exit(0);
  }

  // Extract dimensions and locate the (first) unlimited dimension (time), if any.
  std::vector<netCDF::NcDim> idims = ivar.getDims();
  const int ndim = (int)idims.size();

  std::vector<int> nonunlim;
  nonunlim.reserve(ndim);
  for(int d=0; d<ndim; ++d){
    if(!idims[d].isUnlimited()) nonunlim.push_back(d);
  }

  if((int)nonunlim.size() < 2){
    std::cerr << inam << "ERROR, variable ["<<vname<<"] must have at least 2 non-unlimited dimensions" << std::endl;
    exit(0);
  }

  // Assume first two non-unlimited dimensions correspond to (y, x)
  const int ydim = nonunlim[0];
  const int xdim = nonunlim[1];

  const long ny_file = (long)idims[ydim].getSize();
  const long nx_file = (long)idims[xdim].getSize();
  const long ntot = ny_file * nx_file;

  // Clamp requested range to file domain
  if(seq0 < 0) seq0 = 0;
  if(seq1 >= ntot) seq1 = ntot - 1;
  if(seq1 < seq0) return true;

  // Sanity: ensure provided nx matches file nx
  if(nx_file != (long)nx){
    std::cerr << inam << "ERROR, nx mismatch: provided nx="<<nx<<" but file nx="<<nx_file<< std::endl;
    exit(0);
  }

  // Build baseline start/count (full extents) for all dims
  std::vector<size_t> start(ndim, 0), count(ndim, 1);
  for(int d=0; d<ndim; ++d){
    if(idims[d].isUnlimited()){
      start[d] = (size_t)irec;
      count[d] = 1;
    } else {
      start[d] = 0;
      count[d] = idims[d].getSize();
    }
  }

  const int var_nd = var.ndims();
  if(var_nd < 2){
    std::cerr << inam << "ERROR, input matrix must have at least 2 dimensions" << std::endl;
    exit(0);
  }
  if(var.size(0) != (int)ny_file || var.size(1) != (int)nx_file){
    std::cerr << inam << "ERROR, input matrix dims (y,x)=("<<var.size(0)<<","<<var.size(1)
              <<") do not match file dims ("<<ny_file<<","<<nx_file<<")" << std::endl;
    exit(0);
  }

  // Convert sequential range to (y,x) bounds
  long y0 = seq0 / (long)nx;
  long x0 = seq0 % (long)nx;
  long y1 = seq1 / (long)nx;
  long x1 = seq1 % (long)nx;

  // Pack one (y, x-range) slab into a contiguous buffer.
  auto pack_row_segment = [&](long yy, long xx0, long xcount, std::vector<T> &buf){
    size_t tail = 1;
    for(int k=2; k<var_nd; ++k) tail *= (size_t)var.size(k);
    buf.resize((size_t)xcount * tail);

    size_t p = 0;
    if(var_nd == 2){
      for(long xx=xx0; xx<xx0+xcount; ++xx) buf[p++] = var((int)yy,(int)xx);
    } else if(var_nd == 3){
      for(long xx=xx0; xx<xx0+xcount; ++xx)
        for(int a=0; a<var.size(2); ++a)
          buf[p++] = var((int)yy,(int)xx,a);
    } else if(var_nd == 4){
      for(long xx=xx0; xx<xx0+xcount; ++xx)
        for(int a=0; a<var.size(2); ++a)
          for(int b=0; b<var.size(3); ++b)
            buf[p++] = var((int)yy,(int)xx,a,b);
    } else {
      std::cerr << inam << "ERROR, write_Tstep_pixels only supports matrices with 2 to 4 dimensions" << std::endl;
      exit(0);
    }
  };

  std::vector<T> buf;

  // Write row-by-row so ranges that cross row boundaries work
  for(long yy=y0; yy<=y1; ++yy){
    long seg_x0 = (yy==y0) ? x0 : 0;
    long seg_x1 = (yy==y1) ? x1 : (long)nx-1;
    long xcount = seg_x1 - seg_x0 + 1;
    if(xcount <= 0) continue;

    start[ydim] = (size_t)yy;  count[ydim] = 1;
    start[xdim] = (size_t)seg_x0; count[xdim] = (size_t)xcount;

    pack_row_segment(yy, seg_x0, xcount, buf);
    ivar.putVar(start, count, &buf[0]);
  }

  std::cerr << inam << "writing partial ["<<vname<<"] (t="<< irec
            <<") seq=["<<seq0<<","<<seq1<<"] to "<<file<<std::endl;
  return true;
}

  
  template <class T> void initVar(std::string vname, std::vector<std::string> dnames){
    
    std::string inam = "io::initVar: ";
    
    // Check if variable already exists
    for(auto &it: vars){
      if(it.getName().compare(vname) == 0){
	std::cerr<<inam<<"WARNING, trying to initialize an existing variable ["<<vname<<
	  "] in "<<file<<std::endl;
	return;
      }
    }
    
    
    // Identify which dimensions are desired
    std::vector<netCDF::NcDim> idims;
    int nun = 0;
    for(int ii=0;ii<(int)dnames.size();ii++) {
      for(auto &it: dims){
	if(it.getName().compare(dnames[ii]) == 0) {
	  idims.push_back(it);
	  if(it.isUnlimited()) nun++;
	}
      }
    }
    
    if(typeid(T) == typeid(double)) vars.push_back(ifile->addVar(vname, netCDF::ncDouble,idims));
    if(typeid(T) == typeid(float))  vars.push_back(ifile->addVar(vname, netCDF::ncFloat,idims));
    if(typeid(T) == typeid(int))    vars.push_back(ifile->addVar(vname, netCDF::ncInt,   idims));
    if(typeid(T) == typeid(short))  vars.push_back(ifile->addVar(vname, netCDF::ncShort, idims));
    
  }
  
  
  bool is_var_defined(const std::string &vname){
    
    //
    // Check which element corresponds to vname
    //
    bool found = false;
    for(int ii = 0; ii<(int)vars.size(); ii++){
     
      if(vars[ii].getName().compare(vname) == 0){
	found = true;
	break;
      } 
    }
    if(!found) return false;
    else return true;
  }


  template <class T> bool writeOne(std::string vname, mat<T> &res){

    /* --- Get dimensions ---*/
    std::vector<int> dims = res.getdims();


    /* --- init dim names ---*/
    std::vector<std::string> dnames;
    dnames.resize(dims.size());

    int k = 0;
    for(auto &it: dnames) it = "d"+std::to_string(k++);

    initDim(dnames, dims);



    /* --- init var --- */
    initVar<T>(vname, dnames);
   

    /* --- Write var to disk ---*/
    write_Tstep<double>(vname, res, 0);


    return true;
  }

 
};

//

#endif
