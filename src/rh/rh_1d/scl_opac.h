#ifndef SOPAC_H
#define SOPAC_H

#include "../atmos.h"

extern Atmosphere atmos;


inline double *scl_opac(int ndep, double *chi)
{
  register int i;
  //int ref_index = 0;
  //ActiveSet *as;
  double *rat = (double*)malloc(ndep * sizeof(double));


  /* --- Extract chi_cont --- */
  /*
  Locate(spectrum.Nspect, spectrum.lambda, atmos.lambda_ref, &ref_index);
  as = &spectrum.as[ref_index];
  alloc_as(ref_index, FALSE);
  
  readBackground_j(ref_index, 0, 0);
  // for(i=0;i<ndep;i++) rat[i] =  chi[i] / as->chi_c[i];

  free_as(ref_index, FALSE);
  */
  for(i=0;i<ndep;i++) rat[i] =  chi[i] / atmos.rho[i];

  return rat;
}

#endif

