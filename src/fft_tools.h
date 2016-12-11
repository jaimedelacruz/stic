/* -------------------------------------
   
  Various fftw based routines in template format

  Coded by J. de la Cruz Rodriguez (ISP-SU 2016)

  -------------------------------------- */

#ifndef FFTTOOLS_H
#define FFTTOOLS_H

#include <cmath>
#include <algorithm>
#include <complex>
#include <cstring>
#include <fftw3.h>
#include <cstdio>

namespace mfft{
  
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

      double psf_tot = 0.0;
      for(size_t ii=0; ii<n1; ii++) psf_tot += psf[ii];
      psf_tot = 1.0 / (psf_tot * npad);
      //
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

  
  
  /* --- 
     1D FFTW convolution class, useful to perform many convolutions with the
     same PSF (e.g., inversions) because the PSF is only transformed once
     --- */
  
  template <class T> class fftconv1D {
  protected:
    size_t npad, n, n1, nft;
    std::complex<double> *otf, *ft;
    fftw_plan fplan, bplan;
    double *padded;
    bool started_plans;
  public:
  /* ------------------------------------------------------------------------------- */
    
  fftconv1D(size_t n_in, size_t n_psf, T *psf):
    npad(0),n(0),n1(0),nft(0),otf(NULL),ft(NULL), padded(NULL), started_plans(false){

      /* --- define dimensions --- */
      
      n = n_in, n1 = n_psf, npad = ((n1/2)*2 == n1) ? n1+n : n1+n-1;
      nft = npad/2 + 1;

      /* --- allocate arrays --- */
      
      double *ppsf   = new double [npad]();
      padded         = new double [npad]();
      //
      ft  = new std::complex<double> [nft]();
      otf = new std::complex<double> [nft]();

      
      /* --- shift PSF 1/2 of the elements of the PSF cyclicly. Apply normalizations --- */
      
      double psf_tot = 0.0;
      for(size_t ii=0; ii<n1; ii++) psf_tot += psf[ii];
      psf_tot = 1.0 / (psf_tot * npad);
      //
      for(size_t ii = 0; ii<n1; ii++) ppsf[ii] = (double)psf[ii] * psf_tot;
      std::rotate(&ppsf[0], &ppsf[n1/2], &ppsf[npad]);

      
      /* --- Init forward and backward plans --- */

      fplan = fftw_plan_dft_r2c_1d(npad, padded, (fftw_complex*)ft, FFTW_MEASURE);
      bplan = fftw_plan_dft_c2r_1d(npad, (fftw_complex*)ft, padded, FFTW_MEASURE);
      started_plans = true;
      

      /* --- transform psf --- */
      
      fftw_execute_dft_r2c(fplan, ppsf, (fftw_complex*)otf);


      /* --- clean-up --- */
      
      delete [] ppsf;
    }
    /* ------------------------------------------------------------------------------- */
    
    ~fftconv1D(){
      
      if(ft)  delete [] ft;
      if(otf) delete [] otf;
      if(padded) delete [] padded;

      if(started_plans){
	fftw_destroy_plan(fplan);
	fftw_destroy_plan(bplan);
      }

      ft = NULL, otf = NULL, padded = NULL, started_plans = false;
      n = 0, n1 = 0, npad = 0, nft = 0;
    }
  /* ------------------------------------------------------------------------------- */
    
    void convolve(size_t n_in, T *d){
      
      if(n_in != n){
	fprintf(stderr, "error: mth::fftconvol1D::convolve: n_in [%d] != n [%d], not convolving!\n", n_in, n);
	return;
      }

      
      /* --- copy data to padded array --- */
      
      for(size_t ii = 0; ii<n; ii++)         padded[ii] = (double)d[ii];
      for(size_t ii = n; ii<n+n1/2; ii++)    padded[ii] = (double)d[n-1];
      for(size_t ii = n+n1/2; ii<npad; ii++) padded[ii] = (double)d[0];

      
      /* --- Forward transform --- */

      fftw_execute_dft_r2c(fplan, (double*)padded, (fftw_complex*)ft);

      
      /* --- Convolve --- */
      
      for(size_t ii = 0; ii<nft; ii++) ft[ii] *= otf[ii];

      
      /* --- Backwards transform --- */

      fftw_execute(bplan);


      /* --- Copy back data (inplace) --- */

      for(size_t ii = 0; ii<n; ii++) d[ii] = (T)padded[ii];

    }
    
  }; // fftconvol1D class
  /* ------------------------------------------------------------------------------- */


}//namespace

#endif
