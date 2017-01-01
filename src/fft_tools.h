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
#include "math_tools.h"

namespace mfft{
  
  /* ------------------------------------------------------------------------------- */
  
  /* --- Internal memory routines --- */

  template <class T> T **mat2d(size_t nx1, size_t nx2, bool zero = false){
    T **p = new T* [nx1];
    //
    if(zero) p[0] = new T [nx1 * nx2]();
    else p[0] = new T [nx1 * nx2];
    //
    for(size_t x1=1;x1<nx1;++x1) p[x1] = p[x1-1] + nx2;
    return p;
  }
  
  template <class T> void del_mat(T **p){
    delete[] (p[0]);
    delete[] (p);
  }
  
  template <class T> T **var2dim(T *data, size_t nx1, size_t nx2){
    T **p = new T* [nx1];
    p[0] = data;
    for(size_t x1=1;x1<nx1;++x1) p[x1]=p[x1-1] + nx2;
    return p;
  }

  
  /* ------------------------------------------------------------------------------- */
  
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

  template <class T> double **reorder_psf_2D(size_t ny1, size_t nx1, T *psf_in,
					   size_t npy, size_t npx)
    {
      double **ppsf = mat2d<double>(npy, npx, true);
      double **tpsf = mat2d<double>(npx, npy, true);
      //
      T **psf = var2dim<T>(psf_in, ny1, nx1);
      
      double sum = 1.0 / ((double)mth::sum<T>(nx1*ny1, psf_in) * npx * npy);
            
      for(size_t yy=0;yy<ny1;yy++){
	for(size_t xx=0;xx<nx1;xx++) ppsf[yy][xx] = psf[yy][xx] * sum;
	std::rotate(ppsf[yy], &ppsf[yy][nx1/2], &ppsf[yy][npx]);
	for(size_t xx=0;xx<npx;xx++) tpsf[xx][yy] = ppsf[yy][xx];
      }

      for(size_t xx=0;xx<npx;xx++){
	std::rotate(tpsf[xx], &tpsf[xx][ny1/2], &tpsf[xx][npy]);
	for(size_t yy=0;yy<npy;yy++) ppsf[yy][xx] = tpsf[xx][yy];
      }

      
      del_mat<double>(tpsf);
      delete [] psf;

      return ppsf;
    }
  
  /* ------------------------------------------------------------------------------- */

  template <class T> void convolve2D_fft(size_t ny, size_t nx, T *img_in, size_t ny1,
				    size_t nx1, T *psf_in)
    {
      /* --- define dimensions --- */
      
      size_t npx = ((nx1/2) * 2 == nx1) ? nx+nx1 : nx+nx1-1,
	npy = ((ny1/2) * 2 == ny1) ? ny+ny1 : ny+ny1-1;
      size_t nftx = npx, nfty = npy / 2 + 1;
      size_t nft = nftx*nfty;


      /* --- compute median and stdev --- */
      
      double me = mth::sum<T>(nx*ny, img_in)/(double)(nx*ny), st = mth::stdev<T>(nx*ny, img_in);
      st = (st < 1.0e-5) ? st = 1.0 : 1.0/st;
      
      
      /* --- map input image and psf pointers to 2d arrays --- */

      T **img = var2dim<T>(img_in, ny, nx), **psf = var2dim<T>(psf_in, ny1, nx1);

      
      /* --- init temporary arrays --- */

      double **pad = mat2d<double>(npy, npx, true);

      
      /* --- Copy image to padded array. 
	 The padding is done my mirroring the image --- */

      size_t np = npx - nx;
      size_t np2= np-np/2+1;
      for(size_t yy = 0; yy<ny; yy++){
	
	for(size_t xx = 0; xx<nx; xx++)
	  pad[yy][xx] = ((double)img[yy][xx] - me) * st; // copy image
	  
	
	for(size_t xx = 0; xx<np/2; xx++)
	  pad[yy][xx+nx] = (double)pad[yy][nx-xx-1]; // bottom right 1/2
	  
	for(size_t xx = np/2; xx<np; xx++)
	  pad[yy][xx+nx] = (double)pad[yy][np2-xx]; // bottom right 2/2
	
      }

      np = npy-ny;
      np2= np-np/2+1;
      //
      for(size_t yy=0; yy<np/2; yy++){
	//	fprintf(stderr," %d %d %d %d\n",yy, yy+ny, ny-yy, npy);
	memcpy(pad[yy+ny], pad[ny-yy], npx*sizeof(double));	
      }
      for(size_t yy=np/2; yy<np; yy++) {
	//    	fprintf(stderr," %d %d %d %d\n",yy, yy+ny, np2-yy, npy);
	memcpy(pad[yy+ny], pad[np2-yy], npx*sizeof(double));
      }


      /* --- shift and pad psf --- */

      double **ppsf = reorder_psf_2D<T>(ny1, nx1, psf_in, npy, npx);

      

      /* --- start FFTW plans --- */

      std::complex<double> *ft  = new std::complex<double> [nft];
      std::complex<double> *otf = new std::complex<double> [nft];

      fftw_plan fplan = fftw_plan_dft_r2c_2d(npy, npx, &pad[0][0],
					     (fftw_complex*)ft, FFTW_ESTIMATE);
      fftw_plan bplan = fftw_plan_dft_c2r_2d(npy, npx, (fftw_complex*)ft,
					     &pad[0][0], FFTW_ESTIMATE);

      
      /* --- Perform FFTs --- */

      fftw_execute(fplan);
      fftw_execute_dft_r2c(fplan, &ppsf[0][0],  (fftw_complex*)otf);
      fftw_destroy_plan(fplan);


      
      /* --- Convolve data with PSF --- */
      
      for(size_t ii=0;ii<nft;ii++) ft[ii] *= otf[ii];
      delete [] otf;

      
      /* --- Convert back --- */

      fftw_execute(bplan);
      fftw_destroy_plan(bplan);
      delete [] ft;

      
      /* --- Copy in-place --- */

      for(size_t yy=0; yy<ny; yy++)
	for(size_t xx=0;xx<nx;xx++)
	  img[yy][xx] = (T)(pad[yy][xx] / st + me);
      

      
      /* --- Clean-up --- */
      
      del_mat<double>(pad);
      del_mat<double>(ppsf);
      delete [] psf;
      delete [] img;
    }
  
  /* ------------------------------------------------------------------------------- */



  
}//namespace

#endif
