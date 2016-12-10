/* -------------------------------------
   
  Various math routines

  Coded by J. de la Cruz Rodriguez (ISP-SU 2016)

  -------------------------------------- */

#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <algorithm>
#include <complex>
#include <cstring>
#include <fftw3.h>
#include <cstdio>

namespace mth{

  
  /* --- sqr --- */
  
  template <class T> T sqr(T a){return a*a;};
  /* ------------------------------------------------------------------------------- */
  
  /* --- Sum, mean, stdev --- */
  
  template <class T> T sum(size_t n, T* a){
    T tot = 0;
    for(size_t ii=0;ii<n;ii++) tot += a[ii];
    return tot;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> T mean(size_t n, T* a){
    T tot = sum<T>(n, a);
    return tot/n;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> double stdev(size_t n, T* a){
    double res = 0.0;
    double a_mean = sum<T>(n, a) / (double)n;
    //
    for(size_t ii=0; ii<n; ii++) res += sqr<double>((double)a[ii]-a_mean);
    //
    return sqrt(res/(n-1.0));
  }
  /* ------------------------------------------------------------------------------- */
  
  /* --- more accurate Kahan addition --- */
  
  template <class T> double ksum(size_t n, T* arr){
    
    long double sum = 0.0L, c = 0.0L;
    
    for(size_t kk = 0; kk<n; kk++){
      long double y = arr[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */
  
  
  
  /* --- 1D FFTW convolution --- */
  
  template <class T> void convolve1D_fft(int n, double *d, int n1, double *psf)
    {
      
      int npad = n+n1;
      if((n1/2)*2 != n1) npad -= 1; // if odd, only need n1-1 points to pad the data
      
      
      /* --- allocate temp array to store padded data and PSF and copy data--- */
      
      double *padded = new double [npad], *ppsf = new double [npad]();
      //
      for(size_t ii = 0; ii<n; ii++) padded[ii] = (double)d[ii];
      for(size_t ii = n; ii<n+n1/2; ii++) padded[ii] = (double)d[n-1];
      for(size_t ii = n+n1/2; ii<npad; ii++) padded[ii] = (double)d[0];
      
      
      /* --- shift PSF 1/2 of the elements of the PSF cyclicly. Apply normalizations --- */
      
      double psf_tot = 1.0 / (sum<T>(n1, psf) * npad);
      for(size_t ii = 0; ii<n1; ii++) ppsf[ii] = (double)psf[ii] * psf_tot;
      std::rotate(&ppsf[0], &ppsf[n1/2], &ppsf[npad]);
      
      
      
      /* --- init FFT plans and execute FFTW --- */
      
      int nft = npad/2 + 1;
      std::complex<double> *ft  = new std::complex<double> [nft]();
      std::complex<double> *otf = new std::complex<double> [nft]();
      
      
      fftw_plan fplan = fftw_plan_dft_r2c_1d(npad, padded, (fftw_complex*)ft, FFTW_ESTIMATE);
      fftw_plan bplan = fftw_plan_dft_c2r_1d(npad, (fftw_complex*)ft, padded, FFTW_ESTIMATE);
      
      
      
      /* --- Compute fft of PSF and data --- */
      
      fftw_execute(fplan);
      fftw_execute_dft_r2c(fplan, ppsf, (fftw_complex*)otf); 
      
      
      /* --- Convolve, multiplication is overloaded for complex numbers --- */
      
      for(size_t ii=0; ii<nft; ii++) ft[ii] *= otf[ii];
      
      
      /* --- Transform back --- */
      
      fftw_execute(bplan);
      
      
      /* --- Destroy fftw plans --- */
      
      fftw_destroy_plan(fplan);
      fftw_destroy_plan(bplan);
      
      
      /* --- copy back result --- */
  
      for(size_t ii = 0; ii<n; ii++) d[ii] = (T)(padded[ii]);
      
      
      /* --- clean-up --- */
      
      delete [] padded;
      delete [] ppsf;
      delete [] ft;
      delete [] otf;
    }
  /* ------------------------------------------------------------------------------- */
  
} // namespace

#endif
