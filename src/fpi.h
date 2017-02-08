#ifndef SFPI_H
#define SFPI_H

#include <vector>
#include <iostream>
#include <fftw3.h>
#include <string>
#include "instruments.h"
#include "input.h"
#include "spectral.h"


/* --- Class definitions --- */

class sfpi: public instrument{
 public:
  region_t reg;
  int nt;
  std::vector<spec_ft> ft;
  std::vector<double> pref;
    
  /* --- constructor/Destructor --- */

  sfpi(){};
  sfpi(region_t &in, int nthreads = 1);
  ~sfpi();
  
  /* --- Prototypes --- */
  void degrade_one(spec_ft &ft, double *dat, int ns);
  void degrade(mat<double> &syn, bool spectral = true, bool spatial = false, int ntt = 1);
  void degrade(double *syn, int ns);
  
};


#endif
