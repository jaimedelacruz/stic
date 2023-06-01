#include <vector>
#include <iostream>
#include <cstdio>
#include <algorithm>
//#include <omp.h>
#include "instruments.h"
#include "input.h"
#include "specprefilter.h"

#include "io.h"
//
using namespace std;
using namespace netCDF;
//

/* --------------------------------------------------------------------------- */

void specprefilter::init(int npsf, double const *ipsf)
{

  if(npsf != reg.nw){
    fprintf(stderr,"specprefilter::init: error! the prefilter must have the same number of wavelengths than the region -> nw_region=%d != nprefilter=%d\n", reg.nw, npsf);
    exit(1);
  }
    
  iprof.resize(npsf);
  for(int ii=0; ii<npsf; ++ii)
    iprof[ii] = ipsf[ii];
  

}


specprefilter::specprefilter(region_t &in, int nthreads): reg(in){};

/* --------------------------------------------------------------------------- */

void degrade_one(double* const dat, int const ns, int const nw, std::vector<double> const& iprof)
{ 
  for(int ss = 0; ss<ns; ss++){
    double sum = 0.0, sum2=0.0;
    for(int ww=0; ww<nw; ++ww){
      sum += dat[ww*ns+ss]*iprof[ww];
      sum2 += iprof[ww];
    }

    sum /= sum2;

    for(int ww=0; ww<nw; ++ww)
      dat[ww*ns+ss] = sum;
    
  } // ss
}


void specprefilter::degrade(double *syn, int ns)
{
    
  degrade_one(&syn[reg.off*ns], ns, reg.nw, iprof);
      
}
  
/* --------------------------------------------------------------------------- */

specprefilter::~specprefilter(){};

/* --------------------------------------------------------------------------- */
