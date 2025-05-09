#include <vector>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <cstdio>
#include <algorithm>
//#include <omp.h>
#include <string>
#include "instruments.h"
#include "input.h"
#include "fpi.h"
#include "io.h"
#include "math_tools.h"
//
using namespace std;
using namespace netCDF;
//

/* --------------------------------------------------------------------------- */

sfpi::sfpi(region_t &in, int nthreads){

  /* --- Copy input --- */

  reg = in;
  nt = max(1,nthreads);
  ft.resize(nthreads);

  /* --- Open PSF data and prefilter pars --- */

  vector<double> ipsf, ppsf;
  vector<complex<double>> otf;
  
  {
    /* --- psf --- */
    
    mat<double> tmp;
    io ifil(file_exists(reg.ifile), NcFile::read, false);
    ifil.read_Tstep<double>("iprof", tmp, 0, false);
    ipsf = tmp.d;

    /* --- pref --- */

    ifil.read_Tstep<double>("pref", tmp, 0, false);
    pref.resize(reg.nw);
    tmp.d[0] = inv_convl(tmp.d[0]);
    
    for(int ii=0; ii<reg.nw; ii++)
      pref[ii] = 1.0 / (1.0 + pow(mth::sqr(2.0 * (reg.w0 + reg.dw * ii - tmp.d[0])/tmp.d[1]), tmp.d[2]));
    
  }


  
  /* --- Calculate padding and dims --- */
  
  int npsf = 0, npad = 0 , n = 0;
  n = reg.nw;
  //
  npsf = (int)ipsf.size();
  if((npsf/2)*2 == npsf) npsf--;
  npad = n + npsf;
  //
  ppsf.resize(npad, 0.0);
  otf.resize(npad/2 + 1);

  

  /* --- Prep. psf --- */

  double sum = 0.0;
  for(int kk = 0; kk < npsf; kk++) sum += ipsf[kk];
  for(int kk = 0; kk < npsf; kk++) ppsf[kk] = ipsf[kk] / (sum * npad);
  std::rotate(&ppsf[0], &ppsf[npsf/2], &ppsf[npad]); // shift half a domain


  /* --- Compute the transform of the PSF --- */
  
  fftw_plan pln = fftw_plan_dft_r2c_1d(npad, (double*)&ppsf[0],
					    (fftw_complex*)(&otf[0]), 
					    FFTW_ESTIMATE);


  /* --- Execute FFT and clean-up plan for PSF --- */
  
  fftw_execute(pln);
  fftw_destroy_plan(pln);

  
  
  /* --- allocate FFTW data --- */

  for(auto &it: ft){
    it.n = reg.nw;
    it.n1 = npsf;
    it.npad = npad;
    it.otf = otf;
    it.dat.resize(npad);
    it.ft.resize(npad/2 + 1);
    //
    it.fwd = fftw_plan_dft_r2c_1d(npad, (double*)&it.dat[0],
				  (fftw_complex*)(&it.ft[0]), 
				  FFTW_MEASURE);
    //
    it.rev = fftw_plan_dft_c2r_1d(npad, (fftw_complex*)(&it.ft[0]),
				  (double*)&it.dat[0],
				  FFTW_MEASURE);
  }


  
  //fprintf(stderr,"spectral::spectral: [%f] -> n=%d, npsf=%d, npad=%d\n", reg.w0, ft[0].n, ft[0].n1, ft[0].npad);

  
  
}


/* --------------------------------------------------------------------------- */

sfpi::~sfpi(){
  
  /* --- Deallocate FFTW plans --- */

  for(auto &it: ft){
    fftw_destroy_plan(it.fwd);
    fftw_destroy_plan(it.rev);
    it.dat.clear();
    it.otf.clear();
  }

  ft.clear();
  
}

/* --------------------------------------------------------------------------- */

void sfpi::degrade(mat<double> &syn, bool spectral, bool spatial, int ntt)
{
  
  /* --- Dims --- */
  
  int npix = syn.size(0), off = reg.off, ipix = 0, tid = 0, ns = syn.size(2);
  
  
  /* --- Start parallel block --- */
  
#pragma omp parallel default(shared) private(ipix, tid) num_threads(nt)
  {
    
    //tid = omp_get_thread_num();
#pragma omp for
    for(ipix = 0; ipix<npix; ipix++)
      degrade_one(ft[tid], &syn(ipix, off, 0), ns);
  }
}
void sfpi::degrade(double *syn, int ns)
{
    
  degrade_one(ft[0], &syn[reg.off*ns], ns);
      
}
  

void sfpi::degrade_one(spec_ft &ift, double *dat, int ns)
{
  
  int nw = ift.n;
  //double (&idat)[nw][ns] = *reinterpret_cast<double (*)[nw][ns]>(dat);

  
  
  for(int ss = 0; ss<ns; ss++){
    
    /* --- pad data array --- */
    
    for(int kk = 0; kk<ift.npad; kk++){
      if(kk < ift.n                               ) ift.dat[kk] = dat[kk*ns+ss] * pref[kk];
      else if(kk >= ift.n && kk < ift.n + ift.n1/2) ift.dat[kk] = dat[(ift.n-1)*ns+ss] * pref[ift.n-1];
      else ift.dat[kk] = dat[ss]*pref[0];//idat[0][ss];
    }

    
    /* --- FFT FORWARD--- */
    
    fftw_execute(ift.fwd);
    
    
    /* --- Convolve --- */
    
    for(int kk = 0; kk<(ift.npad/2+1); kk++)
      ift.ft[kk] *= ift.otf[kk];
    
    
    /* --- FFT BACKWARD --- */
    
    fftw_execute(ift.rev);
    
    
    
    /* --- Copy back in place --- */
    
    for(int kk = 0; kk < ift.n; kk++)
      dat[kk*ns+ss] = ift.dat[kk] / pref[kk];

  } // ss
  
  
}
