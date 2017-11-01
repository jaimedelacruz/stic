#ifndef SOPAC_H
#define SOPAC_H

#include "../atmos.h"

extern Atmosphere atmos;


inline double *scl_opac(int ndep, double *chi)
{
  register int i;
  double *rat = (double*)malloc(ndep * sizeof(double));
  for(i=0;i<ndep;i++) rat[i] =  chi[i] / atmos.rho[i];
  return rat;
}

#endif

