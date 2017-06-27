#ifndef SFPIGEN_H
#define SFPIGEN_H

#include <vector>
#include <iostream>
#include <fftw3.h>
#include <string>
#include "instruments.h"
#include "input.h"
#include "spectral.h"
#include "cmemt.h"


/* --- Class definitions --- */

class sfpigen: public instrument{
 private:
  double hr, lr, hc, lc, w0;
  mat<float> erh, erl, ech, ecl;
 public:
  region_t reg;
  int nt, npad, npsf;
  spec_ft ft;
  std::vector<double> pref, ipsf, ppsf, tw;
  
  
  
  /* --- constructor/Destructor --- */

  sfpigen(){};
  sfpigen(region_t &in, int nthreads = 1);
  ~sfpigen();
  
  /* --- Prototypes --- */
  void degrade_one(spec_ft &ft, double *dat, int ns);
  //void degrade(mat<double> &syn, bool spectral = true, bool spatial = false, int ntt = 1);
  void degrade(double *syn, int ns);
  void dual_fpi(double ech, double ecl, double erh, double erl);
  void update(size_t ipix);
  void prepPSF();

};


#endif
