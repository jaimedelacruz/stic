/*
  IO interface for the CXX4 netCDF routines
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
  Dependencies: netcdf library with cxx4 bindings (not cxx!)
 */
#include <iostream>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>      // std::ifstream
#include <map>
#include "io.h"
#include "cmemt.h"
//
using namespace std;
using namespace netCDF;
//
void io::initRead(string filename, netCDF::NcFile::FileMode mode){
  //
  string inam = "io::initRead: ";
  file = filename;
  ifile = new NcFile(file, mode);
  //
  int nvar = ifile->getVarCount();
  //
  if(nvar > 0){
    multimap<string, NcVar> tmp = ifile->getVars();
    for ( auto &it: tmp) {
      vars.push_back(it.second);
    }
    
    // printout info
    cout << inam << "Found "<<nvar<<" variables"<< " in "<<file<<" ["<<vars[0].getName();
    for(int vv = 1;vv<nvar;vv++) cout <<", "<< vars[vv].getName();
    cout<<"]"<<endl;

  } //else{cout << inam <<"Warning, file is empty"<<endl;}

}
//
void io::initDim(string vname, int size){
  if(size == 0) dims.push_back(ifile->addDim(vname));
  else dims.push_back(ifile->addDim(vname,size));
  
}
//
void io::initDim(vector<string> vname, vector<int> vsize){
  string inam = "io::initDim: ";
  
  // Check that both arrays have the same size
  if(vname.size() != vsize.size()){
    cerr << "ERROR, arrays containing dimensions and sizes must have the same size"<<endl;
    exit(0);
  }

  for(int ii=0; ii<(int)vsize.size();ii++){
    if(vsize[ii] == 0) dims.push_back(ifile->addDim(vname[ii]));
    else dims.push_back(ifile->addDim(vname[ii], vsize[ii]));
  }
}
//
vector<int> io::dimSize(string vname){
  // bool exist = false;

  vector<int> idim;
  
  for(auto &it: vars){
    if(it.getName().compare(vname) == 0){
      idim.resize(it.getDimCount());
      int k = 0;
      for(auto &bla: it.getDims()) idim[k++] = int(bla.getSize());
    }
  }
  return idim;
}
//
std::string file_exists(const std::string& name) {
  std::ifstream f(name.c_str(),std::ifstream::in);
    if (f.good()) {
      f.close();
      return name;
    } else {
      f.close();
      std::cerr << "file_check: ERROR, file "<<name<<" does not exist!"<<std::endl;
      exit(0);
     return name;
    }
}
//
bool bfile_exists(const std::string& name) {
  std::ifstream f(name.c_str(),std::ifstream::in);
    if (f.good()) {
      f.close();
      return true;
    } else {
      f.close();
     return false;
    }
}
