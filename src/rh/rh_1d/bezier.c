/* ----------------------------------------------------------------------- 

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
           2017-03-12, JdlCR: Created!
	   2017-04-27, JdlCR: Removed linear approximation from the last interval,
	                      instead we set the last derivative using the linear
			      approximation. This ensures continuity and derivability
			      in the upwind point of the last interval.
	   2018-02-22, JdlCR: Stimulated emission can cause negative absortion, so removed the
	                      cap imposed in the control points (to be possitive). Also
			      introduced a better harmonic derivative implementation from 
			      Steffen (1990).

   ----------------------------------------------------------------------- */


#include <math.h>
#include <string.h>    // memcpy, memset
#include <stdlib.h>
//
#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "bezier.h"

/* --- Macros --                          -------------- */

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define swap(a,b,c) (c=a,a=b,b=c)

/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


/* --- Identity matrix ---- */

static const double ident[4][4] =
  {{1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0}};


/* ----------------------------------------------------------------------- */

void solveLinearFast(double A[4][4], double B[4])
{
  
  // --- Simple Gaussian elimination with partial pivoting --- //
  
  // A is the matrix with coeffs that multiply X. B is the right hand
  // term (on input). The result is overwritten into B. All operations
  // are done in-place
  // It uses the macro "swap" to swap the values of two rows when a pivot
  // is found

  // The algorithm is a pretty standard Gauss-Jordan algorithm, here optimized a
  // bit for performance.
  
  register int i,j,k;
  int maxrow,swapme;
  double maxel,tmp;
  
  
  for (i=0; i<4; i++) {

    // Find pivot

    maxel = fabs(A[i][i]);
    maxrow = i, swapme = 0;
    
    for (k=i+1; k<4; k++){
      tmp = fabs(A[k][i]);
      if(tmp > maxel){
	maxel = tmp;
	maxrow = k, swapme = 1;
      }
    }
    
    // swap
    if(swapme){
      for (k=i; k<4;k++) swap(A[maxrow][k],A[i][k],tmp);  
      swap(B[maxrow],B[i],tmp); 
    }

    // Set to zero relevant columns
    for (k=i+1; k<4; k++){
      tmp = -A[k][i]/A[i][i];
      A[k][i] = 0.0;
      for ( j=i+1; j<4; j++) {
	A[k][j] += tmp * A[i][j];
      }
      B[k] += tmp*B[i];
    }
  }

  
  // Solve upper triagonal system and store in B
  
  for (i=3; i>=0; i--) {
    B[i] /= A[i][i];
    for ( k=i-1;k>=0; k--) {
      B[k] -= A[k][i] * B[i];
    }
  }
  
}

/* ----------------------------------------------------------------------- */

inline void m4inv(double MI[4][4]){

  /* ---
     In-place Shipley-Coleman matrix inversion
     Fast, but ... how accurate?? Pivoting is always done in the diagonal.
     Copied here just in case the SIMD matrix inversion
     gives troubles.
     --- */
  
  register int k, i, j;
  for( k=0;k<4;k++){
    MI[k][k]=-1.0/MI[k][k];         // the pivot element
    for( i=0;i<4;++i) if(i!=k) MI[i][k]*=MI[k][k];//the pivot column
    for( i=0;i<4;++i)           //elements not in a pivot row or column
      if(i!=k)
	for( j=0;j<4;++j)
	  if(j!=k)
	    MI[i][j]+=MI[i][k]*MI[k][j];
    for( i=0;i<4;++i)           //elements in a pivot row
      if(i!=k)
	MI[k][i]*=MI[k][k];
  }
  for( i=0;i<4;++i)
    for( j=0;j<4;++j) MI[i][j]=-MI[i][j];
  return;
}
/* --------------------------------------------------------------- */

inline double sign(const double val)
{
  if(val<0) return -1.0;
  else      return  1.0;
}

/* --------------------------------------------------------------- */

inline double cent_deriv(double dsup,double dsdn,
			 double chiup,double chic, double chidn)
{
  /* --- Derivatives from Fritsch & Butland (1984) --- */

  double fim1, fi, alpha, wprime;
  
  fim1=(chic-chiup)/dsup;
  fi=(chidn-chic)/dsdn;

  if ((fim1*fi) > 0.0) {
    alpha = 0.3333333333333333333333333333333 * ( 1.0 + dsdn/(dsdn+dsup) );
    return (fim1*fi) / ( (1.0-alpha) * fim1 + alpha*fi );
   } else {
    return 0.0;
  }
}
/* --------------------------------------------------------------- */

/* inline double cent_deriv(double odx,double dx, */
/* 			 double chiup,double chic, double chidn) */
/* { */
/*   /\* --- Derivatives from Steffen (1990) --- *\/ */

/*   double ody=(chic-chiup)/odx, dy=(chidn-chic)/dx, wprime = 0.0, p; */
  
  
/*   if ((ody*dy) > 0) { */
/*     p = (dy * odx + ody * dx) / (dx + odx); */
/*     wprime = (2.0*sign(dy)) * min(min(fabs(ody), fabs(dy)), 0.5*fabs(p)); */
/*   } */
  
/*   return wprime; */
/* } */

/* --------------------------------------------------------------- */

inline void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
			   double chiup[4][4], double chic[4][4], double chidn[4][4])
{
  register int i,j;
  
  for(j=0;j<4;j++)
    for(i=0;i<4;i++)
      wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i], chidn[j][i]);
  
}

/* -------------------------------------------------------------------------- */

inline void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
		    double chiup[4], double chic[4], double chidn[4])
{
  register int i;
  
  for(i=0;i<4;i++)
    wprime[i] = cent_deriv(dsup, dsdn, chiup[i], chic[i], chidn[i]);
  
}

/* -------------------------------------------------------------------------- */

  /* --- matrix multiplication --- */

inline void m4m(double a[4][4], double b[4][4], double c[4][4])
{
  register int i, j, k;
  memset(&c[0][0],0,sizeof(double)*16);
  for(j = 0; j<4; j++)
    for(i = 0; i<4; i++)
      for(k = 0; k<4; k++)
	c[j][i] += a[k][i]*b[j][k]; 
}

/* -------------------------------------------------------------------------- */

/* --- matrix/vector multiplication, we use matrix as float to 
   be able to use Intel's matrix inversion --- */

inline void m4v(double a[4][4], double b[4], double c[4]){
  register int k, i;
  memset(&c[0],0,sizeof(double)*4);
  for(i = 0; i<4; i++)
    for(k = 0; k<4; k++)
      c[i] += ((double)a[i][k]) * b[k];
}

/* -------------------------------------------------------------------------- */

/* --- Extracts the Source vector at depth-points k --- */

inline void Svec(int k, double **S, double *Sf)
{
  Sf[0] = S[0][k], Sf[1] = S[1][k], Sf[2] = S[2][k], Sf[3] = S[3][k];
}

/* -------------------------------------------------------------------------- */

inline void Bezier3_coeffs(double dt, double *alpha, double *beta,
		    double *gamma, double *theta, double *eps)
{
  /* ---

     Integration coeffs. for cubic Bezier interpolants
     Use Taylor expansion if dtau is small
     
     --- */
  
  double dt2 = dt*dt, dt3 = dt2 * dt;//,  dt4;
  if(dt < 0.05){
    *eps =   1.0 - dt + 0.50 * dt2 - dt3 *  0.166666666666666666666666667;
    *alpha = 0.25 * dt - 0.20 * dt2 + dt3 * 0.083333333333333333333333333;// - dt4 / 840.0;
    *beta  = 0.25 * dt - 0.05 * dt2 + dt3 * 0.008333333333333333333333333;//  - dt4 / 42.0;
    *gamma = 0.25 * dt - 0.15 * dt2 + dt3 * 0.050;// - dt4 / 210.0;
    *theta = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025;// - dt4 / 84.0;
    return;
  }else{
    if((dt > 30.)){
      *eps = 0.0;
      *alpha = 6.0 / dt3;
      *beta = 1.0 - (*alpha)*(1.0-dt) - 3.0/dt;
      *gamma = *alpha*(dt-3.0);
      *theta = 3.0*(*alpha + (dt-4.0)/dt2);
      return;
    }else{
      *eps = exp(-dt);
      *alpha = -(-6.0+(6.0+6.0*dt+3.0*dt2+dt3)*(*eps))/dt3;
      *beta  = (-6.0 + dt*(6.0+dt*(dt-3.0)) +6.0*(*eps)) / dt3;
      *gamma = 3.0 * (2.0*dt-6.0 + *eps*(6.0+dt*(dt+4.0))) / dt3;
      *theta = 3.0 * (-2.0*(*eps)*(dt+3.0) +6.0+(dt-4.0)*dt) / dt3;
      return;
    }
  }
}
void PiecewiseStokesBezier3(int nspect, int mu, bool_t to_obs,
			    double *chi, double **S, double **I, double *Psi)
{
  /* ---
     Cubic DELO-Bezier solver for polarized light
     Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

     Reference(s):
     J. de la Cruz Rodriguez & N. Piskunov (2013)

     Comments: 
     
     --- */
  
  static const char routineName[] = "PiecewiseStokesBezier3";
  static const int siz_mat = 16*sizeof(double), siz_vec = 4*sizeof(double);

  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk, k_last, ref_index;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2, dchi_dn,
    I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt03;
  double dsup,dsdn,dt,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4];//, V1[4];
  double imu = 1.0 / geometry.muz[mu];
  double Md[4][4];
  ActiveSet *as;
  double *z, *chi1;
  
  if(FALSE){
    chi1 = (double*)malloc(Ndep * sizeof(double));//scl_opac(Ndep, chi);
    for(i=0;i<Ndep;i++) chi1[i] = chi[i] / atmos.rho[i];
    z = geometry.cmass;//
  }else{
    chi1 = chi;
    z = geometry.height;
  }
  
  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi1[k_start] + chi1[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  
  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_upw[0] = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Ibottom[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Itop[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    default:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[TOP]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  
  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;
  
  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi1[k] - chi1[k-dk])/dsup;

  
  /* ---  dchi/ds at central point --- */
  
  dchi_c = cent_deriv(dsup,dsdn,chi1[k-dk],chi1[k],chi1[k+dk]);
  
  
  /* --- upwind path_length (BEzier3 integration) --- */

  c2 = chi1[k]    - (dsup*0.3333333333333333333333333) * dchi_c;
  c1 = chi1[k-dk] + (dsup*0.3333333333333333333333333) * dchi_up;

  
  dtau_uw = 0.25 * dsup * (chi1[k] + chi1[k-dk] + c1 + c2);

  
  /* --- Ku, K0 and dKu, dSu --- */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);

  Svec(k_start,    S, Su);
  Svec(k_start+dk, S, S0);

  
  /* --- Assume side derivative in the first interval --- */
  
  for(n=0;n<4;n++){
    dSu[n] = (S0[n] - Su[n]) / dtau_uw;
    
    for(m=0;m<4;m++)
      dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
  }

  
  /* --- Solve transfer along ray --                   -------------- */
  
  k_last = k_end  + dk;
  //
  for (k = k_start+dk; k != k_last;  k += dk) {      
      
    /* ---  dchi/ds at downwind point --- */

    if(k != k_end){
      
      dsdn = fabs(z[k+dk] - z[k]) * imu;
      
      if(fabs(k-k_end)>1){
	dsdn2=fabs(z[k+2*dk] - z[k+dk]) * imu;
	dchi_dn = cent_deriv(dsdn,dsdn2,chi1[k],chi1[k+dk],chi1[k+2*dk]);       
      } else dchi_dn=(chi1[k+dk]-chi1[k])/dsdn;
     
      
    /* --- Make sure that c1 and c2 don't do below zero --- */

      c1 = (chi1[k+dk] - (dsdn*0.333333333333333333333333333) * dchi_dn);
      c2 = (chi1[k]    + (dsdn*0.333333333333333333333333333) * dchi_c);
      
      /* --- Bezier3 integrated dtau --- */
      
      dtau_dw = 0.25 * dsdn * (chi1[k] + chi1[k+dk] + c1 + c2);
      
      /* ---- get algebra in place --- */

      StokesK(nspect, k+dk, chi[k+dk], Kd);
      Svec(k+dk, S, Sd);
      
      cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
      cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);
      
    } else{

      /* --- Last interval, assume linear dependence for the derivatives 
	 at the central point. In that case all the info related to the 
	 downwind point can be ditched --- */
      
      for(n=0;n<4;n++){
	dS0[n] = (S0[n] - Su[n]) / dtau_uw;	
	for(m=0;m<4;m++) dK0[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
      }
    }
    
    dt = dtau_uw, dt03 = dt * 0.3333333333333333333333333333333;

    
    /* --- Bezier3 coeffs. ---- */
      
    Bezier3_coeffs(dt, &beta, &alpha, &theta, &gamma, &eps);

    
    /* --- Diagonal operator --- */
      
    if(Psi) Psi[k] = alpha + gamma;


    
    /* --- Build the linear system of equations to get the intensity --- */
    
    m4m(Ku, Ku, Ma); // Ku # Ku
    m4m(K0, K0, A ); // K0 # K0

    
    for(j=0;j<4;j++){
      for(i=0;i<4;i++){
	Md[j][i] = ident[j][i] + alpha * K0[j][i] - gamma *
	  -(dt03 * (A[j][i] + dK0[j][i] + K0[j][i]) + K0[j][i]);
	  
	Ma[j][i] = eps * ident[j][i] - beta * Ku[j][i] + theta *
	  (dt03 * (Ma[j][i] + dKu[j][i] + Ku[j][i]) - Ku[j][i]);
	  
	Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] - dt03 * Ku[j][i]);
	Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] + dt03 * K0[j][i]);
      }
    }
      
    /* ---
       Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
       (gam * dS0 - theta * dSu) * dtau / 3.0 to compute the 
       right-hand term
       --- */
        
    for(i = 0; i<4; i++){
      V0[i] = 0.0;
      for(j = 0; j<4; j++){
	V0[i] += Ma[i][j] * I[j][k-dk] + Mb[i][j] * Su[j] + Mc[i][j] * S0[j];
      }
      V0[i] += -dt03 * (gamma * dS0[i] - theta * dSu[i]);
    }
      
      
    /* --- Solve linear system to get the intensity --- */
     
    solveLinearFast(Md,V0);

    
    for(i=0;i<4;i++) I[i][k] = V0[i];

    
    /* --- Shift values for next depth --- */
      
    memcpy(Su,   S0, siz_vec);
    memcpy(S0,   Sd, siz_vec);
    memcpy(dSu, dS0, siz_vec);
      
    memcpy(Ku[0],   K0[0], siz_mat);
    memcpy(K0[0],   Kd[0], siz_mat);
    memcpy(dKu[0], dK0[0], siz_mat);
      
    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;
      
  }
  
  if(chi1 != chi)
    free(chi1);
}

/* --------------------------------------------------------------- */

void Piecewise_Bezier3(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  
  /* ---
     Cubic Bezier solver for unpolarized light
     Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

     Reference:
     J. de la Cruz Rodriguez & N. Piskunov (2013), Auer (2003)

     Comments: 
        JdlCR: We only check that the control points of the opacity
	       and source function are always above zero to avoid having
	       a negative interpolant. 
     
     --- */
  
  register int k,i;
  const char routineName[] = "Piecewise_Bezier3";

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, c1, c2, w[3],
         zmu, Bnu[2];
  double dsup,dsdn,dt,dt03,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double dS_up,dS_c,dchi_up,dchi_c,dchi_dn,dsdn2, tmp = 0.0;
  double *z, *chi1;
  
  if(FALSE){
    chi1 = (double*)malloc(Ndep * sizeof(double));//scl_opac(Ndep, chi);
    for(i=0;i<Ndep;i++) chi1[i] = chi[i] / atmos.rho[i];
    z = geometry.cmass;//
  }else{
    chi1 = chi;
    z = geometry.height;
  }
  
  zmu = 1.0 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --      -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * zmu * (chi1[k_start] + chi1[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_upw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Ibottom[nspect][mu];
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu);
      I_upw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Itop[nspect][mu];
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
  
  I[k_start] = I_upw;
  if (Psi) Psi[k_start] = 0.0;
  
  /* set variables for first iteration to allow simple 
     shift for all next iterations */

  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * zmu;
  dsdn = fabs(z[k+dk] - z[k]) * zmu;
  dchi_up= (chi1[k] - chi1[k-dk])/dsup;
  
  /*  dchi/ds at central point */

  dchi_c = cent_deriv(dsup,dsdn,chi1[k-dk],chi1[k],chi1[k+dk]);
  
  /* --- upwind path_length (Bezier3 integration) --- */

  c1 = (chi1[k]    - (dsup*0.333333333333333333) * dchi_c);
  c2 = (chi1[k-dk] + (dsup*0.333333333333333333) * dchi_up);
  
  dtau_uw =  dsup * (chi1[k] + chi1[k-dk] + c1 + c2) * 0.25;

  
  /* dS/dtau at upwind point */

  dS_up = (S[k]-S[k-dk]) / dtau_uw;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if(k != k_end){
      
      /* downwind path length */
      dsdn = fabs(z[k+dk] - z[k]   ) * zmu;

      
      /* dchi/ds at downwind point */
      
      if (fabs(k-k_end)>1) {
	dsdn2=fabs(z[k+2*dk] - z[k+dk]) * zmu;
	dchi_dn = cent_deriv(dsdn,dsdn2,chi1[k],chi1[k+dk],chi1[k+2*dk]);       
      } else {
	dchi_dn=(chi1[k+dk]-chi1[k])/dsdn;
      }


      /* --- Do not clip control points, it fails if there is much estimulated emission --- */
      
      c1 = (chi1[k]    + (dsdn*0.3333333333333333333) * dchi_c);
      c2 = (chi1[k+dk] - (dsdn*0.3333333333333333333) * dchi_dn);
      
      /* downwind optical path length */
      
      dtau_dw =  dsdn * (chi1[k] + chi1[k+dk] + c1 + c2) * 0.25;
      tmp += dtau_uw;

      /* ---  dS/dt at central point --- */
    
      dS_c = cent_deriv(dtau_uw,dtau_dw,S[k-dk],S[k],S[k+dk]);
      
    }else{

      /* --- Last interval, use linear approx for the derivative at 
	 the central point --- */
      
      dS_c = (S[k] -S[k-dk]) / dtau_uw;
    }
    
    
    /* --- compute interpolation parameters --- */
    
    dt03 = dtau_uw*0.33333333333333333333333;
    Bezier3_coeffs(dtau_uw, &beta, &alpha, &theta, &gamma, &eps);
    
    
    
    /* --- Source function control points --- */
       
    c1 = (S[k]    - dt03 * dS_c);
    c2 = (S[k-dk] + dt03 * dS_up);   
    
    /* --- Solve integral in this interval --- */
    
    I[k]= I[k-dk]*eps + alpha*S[k] + beta*S[k-dk] + gamma * c1 + theta * c2; 
    

    /* --- Diagonal operator --- */
    
    if (Psi) Psi[k] = alpha + gamma;
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    dS_up = dS_c;
  }

  if(chi1 != chi)
    free(chi1);
  
}

/* -------------------------------------------------------------------------- */




