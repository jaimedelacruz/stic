/*
  Wavelet class
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
  Dependencies: cmemt.h, GSL (library)
 */
#ifndef WAVELET_H
#define WAVELET_H
//
#include <vector>
#include <string>
#include <iostream>
#include <gsl/gsl_wavelet.h>
#include "cmemt.h"
//
enum wavelet_dir{
  dwt_forward,
  dwt_inverse
};
enum wavelet_type{
  dwt_haar,
  dwt_daub,
  dwt_bspline,
  dwt_haar_c,
  dwt_daub_c,
  dwt_bspline_c
};
//
wavelet_type string2wavelet(const std::string fam);

class wavelet{
 private:
  std::vector<int> dims;
  std::vector<gsl_wavelet*> w;
  std::vector<gsl_wavelet_workspace*> work, work2, work3;
  unsigned int nthreads, ndim;
  
 public:
  const gsl_wavelet_type *wtype;
  unsigned int order;
  std::string fam;
  wavelet_type otype;
  //
  wavelet(const std::vector<int> dims1, unsigned int nt = 1, 
	  const wavelet_type wtype1 = dwt_haar, unsigned int ord = 2){
    init(dims1, nt, wtype1, ord);
  }
  wavelet(){};
  //
  ~wavelet(){
    for(auto &it: w) gsl_wavelet_free (it);
    for(auto &it: work) gsl_wavelet_workspace_free (it);
    for(auto &it: work2) gsl_wavelet_workspace_free (it);
    for(auto &it: work3) gsl_wavelet_workspace_free (it);
    //
    work.resize(0);
    work2.resize(0);
    work3.resize(0);
    dims.resize(0);
  }
  //
  //  void init(std::vector<int> dims1, unsigned int nt = 1);
  void init(const std::vector<int> dims1, unsigned int nt = 1, 
	    const wavelet_type wtype1 = dwt_haar, unsigned int ord = 2);

  void transform(mat<double> &dat, wavelet_dir dir1 = dwt_forward);
  void transformSlices(mat<double> &dat, wavelet_dir idir = dwt_forward, int id1 = 0);
  
};


#endif
