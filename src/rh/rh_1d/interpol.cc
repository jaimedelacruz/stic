#include "interpol.h"
#include "rhf1d.h"
#include <math.h>

void hermitian_interpolation(int n, double *x, double *y, int nn, double *xp, double *yp, int lo)
{

  register int k;
  
  /* ---- Just a wrapper for the C++ routine --- */

  if(lo){
    double *bla = (double*)malloc(sizeof(double)*n);
    
    for(k=0;k<n;k++) bla[k] = log(y[k]);
    linpol<double,double>(n, x, bla, nn, xp, yp);
    for(k=0;k<n;k++) yp[k] = exp(yp[k]);

    free(bla);
  }else{
    
    linpol<double,double>(n, x, y, nn, xp, yp);
    
  }
}
