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
#include "fpigen.h"
#include "io.h"
#include "math_tools.h"
#include "cmemt.h"
//
using namespace std;
using namespace netCDF;
//

/* --------------------------------------------------------------------------- */

void sfpigen::update(size_t ipix)
{

  this->ipix = ipix;
  
  /* --- Compute FPI transmission profile assuming 2 etalons and perpendicular incidence --- */
  
  dual_fpi((double)ech.d[ipix], (double)ecl.d[ipix], (double)erh.d[ipix], (double)erl.d[ipix]);

  
  /* --- Pad PSF and compute FFT --- */
  
  prepPSF();
  
}

/* --------------------------------------------------------------------------- */

void sfpigen::prepPSF()
{
  /* --- 
     Here we pad the PSF and shift it by 1/2 of the domain
     Remember that npsf in this case is set to the closest smaller
     odd number than the region number of elements (npsf)
     --- */

  
   /* --- Prep. psf --- */
  
  double sum = 0.0;
  for(int kk = 0; kk < npsf; kk++) sum += ipsf[kk];
  for(int kk = 0; kk < npsf; kk++) ppsf[kk] = ipsf[kk] / (sum * npad); // with FFTW3 we have to normalize by npad
  std::rotate(&ppsf[0], &ppsf[npsf/2], &ppsf[npad]); // shift half a domain


  /* --- Compute OTF --- */

  fftw_execute(ft.oplan);
  
  
}

/* --------------------------------------------------------------------------- */

void sfpigen::dual_fpi(double ech, double ecl, double erh, double erl)
{

  static const double ca = 6.28318530779; // * cos(3.1415926 / 2)
  
  /* --- Total reflectivity --- */
  
  double thr = hr + erh, tlr = lr + erl;

  
  /* --- Finesse --- */
  
  double fhr = 4.0 * thr / mth::sqr(1.0 - thr), flr = 4.0 * tlr / mth::sqr(1.0 - tlr);


  /* --- Phase --- */

  double phr = hc * ca, plr = lc * ca;

  
  /* --- Transmission profile --- */
  
  for(int ii=0; ii<npsf; ii++){
    ipsf[ii] = (1.0 / (1.0 + flr * mth::sqr(sin(plr / (tw[ii] + ecl + w0)))) * 
		1.0 / (1.0 + fhr * mth::sqr(sin(phr / (tw[ii] + ech + w0)))));
  }

 
  
}
  
/* --------------------------------------------------------------------------- */

sfpigen::sfpigen(region_t &in, int nthreads){

  /* --- Copy input --- */

  reg = in;
  nt = 1; // nthreads is dummy here

  
  /* --- Open calibration data and prefilter pars --- */
  
  {
    /* --- Read calibration data --- */
    
    mat<double> tmp;
    io ifil(file_exists(reg.ifile), NcFile::read, false);
    ifil.read_Tstep<double>("fpi_pars", tmp, 0, false);

    if(tmp.d.size() < 5){ // We need to read an array with 5 elements: hr, lr, hc, lc, w0
      fprintf(stderr,"sfpige::sfpigen: ERROR, calibration file is wronly formatted\n");
      exit(1);
    }
    hr = tmp.d[0], lr = tmp.d[1], hc = tmp.d[2], lc = tmp.d[3], w0 = tmp.d[4];

    ifil.read_Tstep<float>("fpi_erh", erh, 0, false);
    ifil.read_Tstep<float>("fpi_erl", erl, 0, false);
    ifil.read_Tstep<float>("fpi_ech", ech, 0, false);
    ifil.read_Tstep<float>("fpi_ecl", ecl, 0, false);

    
    
    /* --- prefilter --- */

    ifil.read_Tstep<double>("pref", tmp, 0, false);
    pref.resize(reg.nw);
    tmp.d[0] = inv_convl(tmp.d[0]);
    
    for(int ii=0; ii<reg.nw; ii++)
      pref[ii] = 1.0 / (1.0 + pow(mth::sqr(2.0 * (reg.w0 + reg.dw * ii - tmp.d[0])/tmp.d[1]), tmp.d[2]));
    
  }


  /* --- Align etalons at w0 --- */
  
  int nhr = int(0.5 + hc / (w0 * 0.5)), nlr = int(0.5 + lc / (w0 * 0.5));
  hc = nhr * w0 * 0.5, lc = nlr * w0 * 0.5;


  
  /* --- Init grid for PSF --- */

  npsf = ((reg.nw/2)*2 == reg.nw/2) ? reg.nw-1 : reg.nw; // Make it odd
  npad = reg.nw + npsf;

  ft.n = reg.nw, ft.n1 = npsf, ft.npad = npad;
  tw.resize(npsf);

  for(int ii=0; ii<npsf; ii++) tw[ii] = double(ii - npsf/2)*reg.dw; // Symmetric grid
  
  
  /* --- Resize arrays --- */
  
  ppsf.resize(npad, 0.0), ft.otf.resize(npad/2 + 1), ft.ft.resize(npad/2 + 1), tw.resize(npsf);
  ipsf.resize(npsf,0.0);


  

  
  /* --- Compute plans --- */
  
  ft.oplan = fftw_plan_dft_r2c_1d(npad, (double*)&ppsf[0],
				  (fftw_complex*)(&ft.otf[0]), 
				  FFTW_MEASURE);
  
  ft.fwd = fftw_plan_dft_r2c_1d(npad, (double*)&ft.dat[0],
				(fftw_complex*)(&ft.ft[0]), 
				FFTW_MEASURE);

  ft.rev = fftw_plan_dft_c2r_1d(npad, (fftw_complex*)(&ft.ft[0]),
				(double*)&ft.dat[0],
				FFTW_MEASURE);
  
  
}


/* --------------------------------------------------------------------------- */

sfpigen::~sfpigen(){
  
  /* --- Deallocate FFTW plans --- */
  
  fftw_destroy_plan(ft.fwd);
  fftw_destroy_plan(ft.rev);
  fftw_destroy_plan(ft.oplan);
  ft.dat.clear();
  ft.otf.clear();
  
}

/* --------------------------------------------------------------------------- */

void sfpigen::degrade(double *syn, int ns)
{
    
  degrade_one(ft, &syn[reg.off*ns], ns);
      
}
  

void sfpigen::degrade_one(spec_ft &ift, double *dat, int ns)
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
