#ifndef SPEC_H
#define SPEC_H

#include <vector>
#include <iostream>
#include <fftw3.h>
#include <string>
#include "instruments.h"
#include "input.h"

/* --- Struc to store FFTW plans --- */

typedef struct{
  int n, n1, npad;
  fftw_plan fwd, rev, oplan;
  std::vector<double> dat;
  std::vector<std::complex<double>> ft, otf;
  
} spec_ft;



/* --- Class definitions --- */

class spectral: public instrument{
 public:
  region_t reg;
  int nt;
  bool firsttime;

  std::vector<spec_ft> ft;
  

    
  /* --- constructor/Destructor --- */

 spectral():firsttime(true){};
  spectral(region_t &in, int nthreads = 1);
  ~spectral();
  
  /* --- Prototypes --- */
  void degrade_one(spec_ft &ft, double *dat, int ns);
  void degrade(mat<double> &syn, bool spectral = true, bool spatial = false, int ntt = 1);
  void degrade(double *syn, int ns);
  void init(int npsf, double const *psf);
  void update(size_t npsf, double *psf){init(int(npsf), psf);}
    
};


#endif
