/*
  Sparse inversion class
  Authors: Jaime de la Cruz Rodriguez (ISP-SU 2014) & Andres Asensio-Ramos (IAC-2014)
  Dependencies: cmemt.h, wavelet.{h,cc}, GSL (library), instrument.h
 */
#ifndef SPARSE_H
#define SPARSE_H
//
#include <mpi.h>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <gsl/gsl_sort_double.h>
#include "cmemt.h"
#include "wavelet.h"
#include "input.h"
#include "atmosphere.h"
//#include "cmilne.h"
#include "instruments.h"
#include "hinode.h"
#include "depthmodel.h"
#include "interpol.h"

//
enum spthres{ // Two different threshold methods
  spt_hard,
  spt_soft
};
//
class sparse2d{
 private: 
  wavelet wlt;
  unsigned nthreads;
  double thres, stime;
  spthres thres_mod;
  std::vector<int> dims;
  instrument *inst;
  std::vector<atmos*> atm;
  iput_t iput;
  mat<double> final, bsyn, isyn, iweight;
  int npix;
  unsigned long slsize;
  std::vector<double> sparsity, mmax, mmin, scalingParameters, chisq, LCommon;

  //
 public:
  int npar;
  int maxiter;
  double bestLogLike;


  //
  // METHODS
  //
  sparse2d(){};
  sparse2d(iput_t &input, std::vector<int> &dims1, double threshold, 
	   wavelet_type family, unsigned ord, spthres thresholdMode, unsigned nt1);
   
  ~sparse2d(){
    //   dsyn.clear();
    sparsity.clear();
  }


  /* --- prototypes --- */
  void init(iput_t &input, std::vector<int> &dims1, double threshold, 
	    wavelet_type family, unsigned ord, spthres thresholdMode, unsigned nt1);
  void threshold(mat<double> &x, double thr);
  void thresholdSingle(mat<double> &x, double thr);
  void SparseOptimization(mat<double> &obs, mat<double> &x, mat<double> &weights,  mdepthall_t &, mat<double> &pwe);
  double meritFunction(mat<double> &obs, mat<double> &x, mat<double> &dchi, mat<double> &noise, mdepthall_t &m, int compute_gradient = 1);
  void evalModel_mpi(mat<double> &x, mat<double> &syn, mat<double> &dsyn, mdepthall_t &m, int cgrad = 1);
  void evalModel_serial(mat<double> &x, mat<double> &syn, mat<double> &dsyn, mdepthall_t &m);
  void finalTransform(mat<double> &x);
  void transposePars(mat<double> &x, unsigned dir = 0);
  void getSlice(mat<double> &x, mat <double> &out, const int z);
  void putSlice(mat<double> &x, mat <double> &out, const int z);
  void checkParameters(mat<double> &xnew);
  inline double sqr(double var){return var * var;}
  double getStep(double &initalpha, mat<double> &obs, mat<double> &dchi,
		 mat<double> &pp, mat<double> &noise, mdepthall_t &m, double tol = 0.1);
  double my_getStep(double &initalpha, mat<double> &obs, double gy, mat<double> &dchi,
		 mat<double> &pp, mat<double> &noise, mdepthall_t &m, int nint = 1);
  double parabStepInt(double sa, double sb, double sc, double ga, double gb, double gc,
			      mat<double> &obs, mat<double> &dchi, mat<double> &y, mat<double> &noise,
		      mdepthall_t &m, int nint, double &bestchi);
  
  void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb,
	      double &fc, mat<double> &obs, mat<double> &dchi, mat<double> &pp,
	      mat<double> &noise, mdepthall_t &m);
  double brent(double ax, double bx, double cx,double tol,double &xmin,
	       mat<double> &obs, mat<double> &dchi, mat<double> &pp, mat<double> &noise, mdepthall_t &m);
  double evalStep(double alpha, mat<double> &obs, mat<double> &dchi, mat<double> &pp,
		  mat<double> &noise, mdepthall_t &m);
  double getNorm(int n, double *x, int incx = 1);
  mat<double> takeStep(mat<double> &par, mat<double> &grad, double step, bool do_thres=true);
  //
  double gettime(double t0 = -1.0){
    struct timeval dum;
    gettimeofday(&dum, NULL);
    if(t0 < 0.0) return dum.tv_sec + dum.tv_usec * 1.0E-6;
    else return (dum.tv_sec + dum.tv_usec * 1.0E-6) - t0;
  }
}; // class sparse

//
#endif
