/* ------- file: -------------------------- spline.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:57:56 1999   --

       --------------------------                      ----------RH-- */

/* --- Cubic spline interpolation routines. --         -------------- */

 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

static bool_t  ascend;
static int     Ntable;
static double *xtable, xmin, xmax;

/* ------- begin -------------------------- splineCoef.c ------------ */

static double *M = NULL, *ytable;

void splineCoef(int N, double *x, double *y)
{
  register int j;
  static double *u = NULL;

  double  p, *q, hj, hj1, D, D1, mu;

  ascend = (x[1] > x[0]) ? TRUE : FALSE;
  xmin = (ascend) ? x[0] : x[N-1];
  xmax = (ascend) ? x[N-1] : x[0];

  q = M = (double *) realloc(M, N * sizeof(double));
  u = (double *) realloc(u, N * sizeof(double));
  hj = x[1] - x[0];
  D  = (y[1] - y[0]) / hj;

  q[0] = u[0] = 0.0;
  for (j = 1;   j < N-1;  j++) {
    hj1 = x[j+1] - x[j];
    mu  = hj / (hj + hj1);
    D1  = (y[j+1] - y[j]) / hj1;

    p = mu*q[j-1] + 2;
    q[j] = (mu - 1) / p;
    u[j] = ((D1 - D) * 6/(hj + hj1) - mu*u[j-1]) / p;

    hj = hj1;  D = D1;
  }

  M[N - 1] = 0.0;
  for (j = N-2;  j >= 0;  j--) {
    M[j] = q[j]*M[j+1] + u[j];
  }
  Ntable = N;
  xtable = x;  ytable = y;
}
/* ------- end ---------------------------- splineCoef.c ------------ */

/* ---------------------------------------- splineEval.c ------------ */

void splineEval(int N, double *x, double *y, bool_t hunt)
{
  register int n;

  int    j = 0;
  double hj, fx, fx1;

  for (n = 0;  n < N;  n++) {
    if (x[n] <= xmin)
      y[n] = (ascend) ? ytable[0] : ytable[Ntable-1];
    else if (x[n] >= xmax)
      y[n] = (ascend) ? ytable[Ntable-1] : ytable[0];
    else {
      if (hunt) 
	Hunt(Ntable, xtable, x[n], &j);
      else
	Locate(Ntable, xtable, x[n], &j);

      hj  = xtable[j+1] - xtable[j];
      fx  = (x[n] - xtable[j]) / hj;
      fx1 = 1 - fx;

      y[n] = fx1*ytable[j] + fx*ytable[j+1] +
	(fx1*(SQ(fx1) - 1) * M[j] + fx*(SQ(fx) - 1) * M[j+1]) * SQ(hj)/6.0;
    }
  }
}
/* ------- end ---------------------------- splineEval.c ------------ */


double signFortran2(const double val)
{
  return ((val >= 0.0)? 1.0 : -1.0);
}

double cent_deriv_steffen(double odx,double dx,
		   double yu,double y0, double yd)
{
  /* --- Derivatives from Steffen (1990) --- */
  
  const double S0 = (yd - y0) / dx;
  const double Su = (y0 - yu) / odx;
  const double P0 = fabs((Su*dx + S0*odx) / (odx+dx)) * 0.5;
  return (signFortran2(S0) + signFortran2(Su)) * fmin(fabs(Su),fmin(fabs(S0), P0));
}

void splineHermite(int const N, double* const x, double* const y, int const N1, double* const x1, double* const y1)
{
  // Coded by J. de la Cruz Rodriguez (ISP-SU, 2024) //
  
  register int n, j;
  int dn = 1, n0 = 0, n1 = N-1;
  int dj = 1, j0 = 0, j1 = N1-1;
  
  if((x[1]-x[0]) < 0){
    dn = -1, n0 = N-1, n1 = 0;
  }

  if((x1[1]-x1[0]) < 0){
    dj = -1, j0 = N1-1, j1 = 0;
  }
  
  
  // --- first calculate derivatives --- //

  double* const yp = (double*)calloc(N,sizeof(double));
  double odx = 0, dx = 0;
  
  
  yp[0] = (y[n0+dn]-y[n0]) / (x[n0+dn]-x[n0]);
  yp[n1] = (y[n1-dn]-y[n1]) / (x[n1-dn]-x[n1]);
  
  for(n=n0+dn; n != n1; n+=dn){ // avoid both outermost points
    odx = x[n]-x[n-dn];
    dx = x[n+dn]-x[n];
    yp[n] = cent_deriv_steffen(odx,dx,y[n-dn], y[n], y[n+dn]);
  }

  
  // --- Now calculate interpolated values --- //

  //int k=j0;
  double u = 0, u2 = 0, u3 = 0, ypu = 0, ypc = 0;
    
  for(n=n0; n != n1; n+= dn){
    dx = x[n+dn]-x[n];
    ypu = yp[n];
    ypc = yp[n+dn];
    
    for(j=j0; j != j1+dj; j += dj){
      if((x1[j] <= x[n+dn]) && (x1[j] > x[n])){
	u = (x1[j]-x[n])/dx, u2 = u*u, u3 = u2*u;
	y1[j] = (2.0*u3 - 3.0*u2 + 1.0)*y[n] + (u3-2.0*u2+u)*ypu + (3.0*u2-2.0*u3)*y[n+dn] + (u3-u2)*ypc;
	//k+=dj;
      };//else if(x1[j]>x[n+dn]) break; // Not working because, with the collisional rates, T is not ordered in x1.
    }
  } // intervals in the real data 


  
  // --- are there points outside the domain? --- //

  double const pmin = y[n0];
  double const pmax = y[n1];
  double const xmin1 = x[n0];
  double const xmax1 = x[n1];

  
  for(j=0; j<N1; ++j){
    if(x1[j] <= xmin1) y1[j] = pmin;
    else if(x1[j] >= xmax1) y1[j] = pmax;
  }
  
  
  free((void*)yp);
}

