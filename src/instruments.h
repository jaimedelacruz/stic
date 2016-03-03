#ifndef INSTRUMENTS_H
#define INSTRUMENTS_H

#include <complex>
#include <vector>
#include <iostream>
#include "cmemt.h"


class instrument{
 public:
  //
  instrument(){};
  virtual ~instrument(void){};
  virtual void degrade(mat<double> &syn, bool spectral = true, bool spatial = true, int ntt = -1){};
  virtual void degrade(double *syn, int ns){std::cerr<<"instrument::degrade: Dummy method!"<<std::endl;};

};


/* --- Generic coupling functions --- */
/*
void grad_one_thread(int nx, int ny, int nd, int npar, int npx, int npy, int ipix,
		     double *resi, double *rfi, double *gradi, double *psfi);
void compute_grad(int nx, int ny, int nd, int npar, int npx, int npy,
		  double *resi, double *rfi, double *gradi, double *psf, int nthreads);
void  degrade_one_thread(int nx, int ny, int nd, int npx, int npy, int ipix,
			 double *obsi, double *resi, double *psfi);
void degrade(int nx, int ny, int nd, int npx, int npy, double *obs, double *psf, int nthreads);
*/

/* --- Degrade spectrally --- */



#endif
