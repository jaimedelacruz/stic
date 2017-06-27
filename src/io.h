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
