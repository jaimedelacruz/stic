#ifndef HINODE_H
#define HINODE_H
//
#include <vector>
#include <complex>
#include <fftw3.h>
#include "instruments.h"
#include "spectral.h"
#include "cmemt.h"
//
class hinode: public spectral{  
 public:

  hinode(int nthreads, std::vector<int> dims){};
  void degrade(mat<double> &syn, bool spectral = true, bool spatial = true, int ntt = -1);

};



//
#endif
