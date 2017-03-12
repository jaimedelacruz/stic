/* ------- file: -------------------------- piecestokes.c -----------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Feb 23 16:18:46 2004 --

       --------------------------                      ----------RH-- */

/* --- Piecewise integration of the coupled Stokes transfer equations.
       Method is quasi-parabolic DELO method.

  See: - D. E. Rees, G. A. Murphy and C. J. Durrant 1989, ApJ 339,
         1093-1106.

       - H. Socas Navarro, J. Trujillo Bueno and B. Ruiz Cobo 2000,
         "Non-LTE Inversion of Stokes Profiles", ApJ 530, 977.

    -- For boundary condition THERMALIZED use the relation

          (I, Q, U, V) ~= (B - \mu dB / d\tau, 0, 0, 0)

       where dB/\tau is taken in the direction of the ray (i.e. NOT in
       the sense of optical depth, e.g. see Mihalas, 1978, p. 51).

       ------------------------------------------------------------------

       JdlCR: Added cubic Delo-Bezier solver for polarized light.
       See: J. de la Cruz Rodriguez & N. Piskunov (2013).
       
       ------------------------------------------------------------------

       --                                              -------------- */

#include <math.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "rhf1d.h"

/* --- Function prototypes --                          -------------- */

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];
extern MPI_t mpi;


typedef double mat[4][4];
typedef double vec[4];

static const double ident[4][4] = {{1.0, 0.0, 0.0, 0.0},
				   {0.0, 1.0, 0.0, 0.0},
				   {0.0, 0.0, 1.0, 0.0},
				   {0.0, 0.0, 0.0, 1.0}};

/* ------- begin -------------------------- PiecewiseStokes.c ------- */

void PiecewiseStokes(int nspect, int mu, bool_t to_obs,
		     double *chi_I, double **S, double **I, double *Psi)
{
  const char routineName[] = "PiecewiseStokes";
  register int k, n, m;

  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, dS_uw[4], dS_dw[4], c1, c2, w[3],
    I_upw[4], zmu, P[4], Bnu[2], Q[4][4], **R, K[4][4], K_upw[4][4];

  zmu = 0.5 / geometry.muz[mu];
  R = matrix_double(4, 4);

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = zmu * (chi_I[k_start] + chi_I[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);

  StokesK(nspect, k_start, chi_I[k_start], K_upw);

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
  for (n = 0;  n < 4;  n++)
    dS_uw[n] = (S[n][k_start] - S[n][k_start+dk]) / dtau_uw;

  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    w3(dtau_uw, w);
    StokesK(nspect, k, chi_I[k], K);

    if (k != k_end) {
      dtau_dw = zmu * (chi_I[k] + chi_I[k+dk]) *
	fabs(geometry.height[k] - geometry.height[k+dk]);

      for (n = 0;  n < 4;  n++) {
	dS_dw[n] = (S[n][k] - S[n][k+dk]) / dtau_dw;
	c1 = dS_uw[n]*dtau_dw + dS_dw[n]*dtau_uw;
	c2 = dS_uw[n] - dS_dw[n];
	P[n] = w[0]*S[n][k] + (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);
      }
      if (Psi) {
	c1 = dtau_uw - dtau_dw;
	Psi[k] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
      }
    } else {

      /* --- Piecewise linear integration at end of ray -- ---------- */

      for (n = 0;  n < 4;  n++)	P[n] = w[0]*S[n][k] + w[1]*dS_uw[n];
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }

    for (n = 0;  n < 4;  n++) {
      for (m = 0;  m < 4;  m++) {
	Q[n][m] = -w[1]/dtau_uw * K_upw[n][m];
	R[n][m] = (w[0] - w[1]/dtau_uw) * K[n][m];
      }
      Q[n][n] = 1.0 - w[0];
      R[n][n] = 1.0;
    }
    for (n = 0;  n < 4;  n++) {
      for (m = 0;  m < 4;  m++) 
	P[n] += Q[n][m] * I_upw[m];
    }
    /* --- Solve linear equations for I --             -------------- */

    SolveLinearEq(4, R, P, TRUE);
    if (mpi.stop) {
      freeMatrix((void **) R);
      return; /* Get out if there is a singular matrix */
    }
    
    /* --- Store results for Stokes vector --          -------------- */

    for (n = 0;  n < 4;  n++) I[n][k] = P[n];

    /* --- Re-use downwind quantities for next upwind position -- --- */

    dtau_uw = dtau_dw;
    for (n = 0;  n < 4;  n++) {
      I_upw[n] = I[n][k];
      dS_uw[n] = dS_dw[n];
      for (m = 0;  m < 4;  m++) K_upw[n][m] = K[n][m];
    }
  }
  freeMatrix((void **) R);
}
/* ------- end ---------------------------- PiecewiseStokes.c ------- */

inline double cent_deriv(double dsup,double dsdn, 
			 double chiup,double chic,double chidn)
{
  /* --- Derivative Fritsch & Butland (1984), reused J. Leenaarts
     implementation in this same code for unpolarized light --- */

  double fim1, fi, alpha, wprime;
  static const double onethird = 1.0/3.0;
  
  fim1=(chic-chiup)/dsup;
  fi=(chidn-chic)/dsdn;

  if (fim1*fi > 0) {
    alpha = onethird * ( 1.0 + dsdn/(dsdn+dsup) );
    wprime = (fim1*fi) / ( (1.0-alpha) * fim1 + alpha*fi );
  } else {
    wprime=0.0;
  }
  return wprime;
}
/* -------------------------------------------------------------------------- */

inline void cent_deriv_mat(mat wprime, double dsup,double dsdn,
		    mat chiup, mat chic, mat chidn)
{
  register int i,j;
  
  for(j=0;j<4;j++)
    for(i=0;i<4;i++)
      wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i], chidn[j][i]);
  
}

/* -------------------------------------------------------------------------- */

inline void cent_deriv_vec(vec wprime, double dsup,double dsdn,
		    vec chiup, vec chic, vec chidn)
{
  register int i;
  
  for(i=0;i<4;i++)
    wprime[i] = cent_deriv(dsup, dsdn, chiup[i], chic[i], chidn[i]);
  
}
/* -------------------------------------------------------------------------- */

  /* --- matrix multiplication --- */

inline void m4m(mat a, mat b, mat c){
    register int i, j, k;
    memset(&c[0][0],0,sizeof(double)*16);
    for(j = 0; j<4; j++)
      for(i = 0; i<4; i++)
	for(k = 0; k<4; k++)
	  c[j][i] += a[k][i]*b[j][k]; 
  }

/* -------------------------------------------------------------------------- */

/* --- matrix/vector multiplication --- */

inline void m4v(mat a, vec b, vec c){
  register int k, i;
  memset(&c[0],0,sizeof(double)*4);
  for(i = 0; i<4; i++)
    for(k = 0; k<4; k++)
      c[i] += a[i][k] * b[k];
}

/* -------------------------------------------------------------------------- */

inline void Svec(int k, double **S, double *Sf)
{
  Sf[0] = S[0][k], Sf[1] = S[1][k], Sf[2] = S[2][k], Sf[3] = S[3][k];
}

/* -------------------------------------------------------------------------- */

void PiecewiseStokesBezier3(int nspect, int mu, bool_t to_obs,
			    double *chi, double **S, double **I, double *Psi)
{
  /* ---
     Cubic DELO-Bezier solver for polarized light
     Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

     Reference:
     J. de la Cruz Rodriguez & N. Piskunov (2013)

     Comments: 
     
     --- */
  
  const char routineName[] = "PiecewiseStokesBezier3";
  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2,dchi_dn,
    I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt03;
  double dsup,dsdn,dt,dt2,dt3,dt4,eps,alpha,beta,gamma,theta;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4];
  double imu = 1.0 / geometry.muz[mu];
  double **Md = matrix_double(4, 4);
  double *z = geometry.height;

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start+dk]) *
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

  
  for (n = 0;  n < 4;  n++)
    dSu[n] = (S[n][k_start] - S[n][k_start+dk]) / dtau_uw;
  
  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;
  
  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;

  
  /* ---  dchi/ds at central point --- */
  
  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  
  /* --- upwind path_length (BEzier3 integration) --- */

  c2 = max(chi[k]    - (dsup/3.0) * dchi_c,  0.0);
  c1 = max(chi[k-dk] + (dsup/3.0) * dchi_up, 0.0);
  
  dtau_uw = 0.25 * dsup * (chi[k] + chi[k-dk] + c1 + c2);

  
  /* --- Ku, K0 and dKu, dSu --- */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);

  Svec(k_start,    S, Su);
  Svec(k_start+dk, S, S0);

  for(n=0;n<4;n++)
    for(m=0;m<4;m++) dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {      
      
    /* ---  dchi/ds at downwind point --- */
      
    dsdn = fabs(z[k+dk] - z[k]   ) * imu;
      
    if(fabs(k-k_end)>1){
      dsdn2=fabs(z[k+2*dk] - z[k+dk]) * imu;
      dchi_dn = cent_deriv(dsdn,dsdn2,chi[k],chi[k+dk],chi[k+2*dk]);       
    } else dchi_dn=(chi[k+dk]-chi[k])/dsdn;
      
      
    /* --- Make sure that c1 and c2 don't do below zero --- */
      
    c2 = max(chi[k]    + (dsdn/3.0) * dchi_c , 0.0);
    c1 = max(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);
      
      
    /* --- Bezier3 integrated dtau --- */
      
    dtau_dw = 0.25 * dsdn * (chi[k] + chi[k+dk] + c1 + c2);
    dt = dtau_uw;
      
      
    /* --- Bezier3 coeffs. ---- */
      
    dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt2 * dt2, dt03 = dt / 3.0;
      
    if(dt >= 1.e-3){
      eps = exp(-dt);
      alpha = (-6.0 + 6.0 * dt - 3.0 * dt2 + dt3 + 6.0 * eps)        / dt3;
      beta  = (6.0 + (-6.0 - dt * (6.0 + dt * (3.0 + dt))) * eps)    / dt3;
      gamma = 3.0 * (6.0 + (-4.0 + dt)*dt - 2.0 * (3.0 + dt) * eps)  / dt3;
      theta = 3.0 * ( eps * (6.0 + dt2 + 4.0 * dt) + 2.0 * dt - 6.0) / dt3;
    }else{
      eps = 1.0 - dt + 0.5 * dt2 - dt3 / 6.0 + dt4 / 24.0;
      //
      alpha = 0.25 * dt - 0.05 * dt2 + dt3 / 120.0 - dt4 / 840.0;
      beta  = 0.25 * dt - 0.20 * dt2 + dt3 / 12.0  - dt4 / 42.0; 
      gamma = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025 - dt4 / 210.0; 
      theta = 0.25 * dt - 0.15 * dt2 + dt3 * 0.05  - dt4 / 84.0; 
    }
      
    /* --- Diagonal operator --- */
      
    if(Psi) Psi[k] = alpha + gamma;
      
      
    /* ---- get algebra in place --- */
      
    StokesK(nspect, k+dk, chi[k+dk], Kd);
    Svec(k+dk, S, Sd);
    cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
    cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);
      
    m4m(Ku, Ku, Ma); // Ku # Ku
    m4m(K0, K0, A ); // K0 # K0
      
    for(j=0;j<4;j++){
      for(i=0;i<4;i++){
	A[j][i] = ident[j][i] + alpha * K0[j][i] - gamma *
	  -(dt03 * (A[j][i] + dK0[j][i] + K0[j][i]) + K0[j][i]);
	  
	Ma[j][i] = eps * ident[j][i] - beta * Ku[j][i] + theta *
	  (dt03 * (Ma[j][i] + dKu[j][i] + Ku[j][i]) - Ku[j][i]);
	  
	Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] - dt03 * Ku[j][i]);
	Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] + dt03 * K0[j][i]);
      }
    }
      
    /* ---
       Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
       (gam * dS0 - theta * dSu) * dtau / 3.0
       --- */
      
    for(i = 0; i<4; i++){
      V0[i] = 0.0;
      for(j = 0; j<4; j++){
	V0[i] += Ma[i][j] * I[j][k-dk] + Mb[i][j] * Su[j] + Mc[i][j] * S0[j];
      }
      V0[i] += dt03 * (gamma * dS0[i] - theta * dSu[i]);
    }
      
      
    /* --- Solve linear system to get the intensity --- */
      
    memcpy(Md[0], A[0], 16*sizeof(double));
    SolveLinearEq(4, Md, V0, TRUE);
      
    for(i=0;i<4;i++) I[i][k] = V0[i];
      
      
    /* --- Shift values for next depth --- */
      
    memcpy(Su,   S0, 4*sizeof(double));
    memcpy(S0,   Sd, 4*sizeof(double));
    memcpy(dSu, dS0, 4*sizeof(double));
    //memcpy(dS0, dSd, 4*sizeof(double));
      
    memcpy(Ku[0],   K0[0], 16*sizeof(double));
    memcpy(K0[0],   Kd[0], 16*sizeof(double));
    memcpy(dKu[0], dK0[0], 16*sizeof(double));
    //memcpy(dK0[0], dKd[0], 16*sizeof(double));
      
    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;
      
  }
      
  /* --- Linear integration in the last interval --- */
  
  k = k_end;
  dtau_uw = 0.5*imu * (chi[k] + chi[k-dk]) *
    fabs(geometry.height[k] - geometry.height[k-dk]);
  w3(dtau_uw, w);
      
  for (n = 0;  n < 4;  n++) V0[n] = w[0]*S[n][k] + w[1] * dSu[n];
  if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
      
  for (n = 0;  n < 4;  n++) {
    for (m = 0;  m < 4;  m++) {
      A[n][m] = -w[1]/dtau_uw * Ku[n][m];
      Md[n][m] = (w[0] - w[1]/dtau_uw) * K0[n][m];
    }
    A[n][n] = 1.0 - w[0];
    Md[n][n] = 1.0;
  }
      
  for (n = 0;  n < 4;  n++) 
    for (m = 0;  m < 4;  m++) 
      V0[n] += A[n][m] * I[m][k-dk];
      
  SolveLinearEq(4, Md, V0, TRUE);
      
  for (n = 0;  n < 4;  n++) I[n][k] = V0[n];
    
  
  freeMatrix((void **) Md);
}

/* ----------------------------------------------------------------------------- */


void PiecewiseStokesBezier2(int nspect, int mu, bool_t to_obs,
			    double *chi, double **S, double **I, double *Psi)
{
  /* ---
     Cuadratic DELO-Bezier solver for polarized light
     Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

     Reference:
     J. de la Cruz Rodriguez & N. Piskunov (2013)

     Comments: 
     
     --- */
  
  const char routineName[] = "PiecewiseStokesBezier2";
  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2,dchi_dn,
    I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt05;
  double dsup,dsdn,dt,dt2,dt3,dt4,eps,alpha,beta,gamma;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4];
  double imu = 1.0 / geometry.muz[mu];
  double **Md = matrix_double(4, 4);
  double *z = geometry.height;

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start+dk]) *
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

  
  for (n = 0;  n < 4;  n++)
    dSu[n] = (S[n][k_start] - S[n][k_start+dk]) / dtau_uw;
  
  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;
  
  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;

  
  /* ---  dchi/ds at central point --- */
  
  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  
  /* --- upwind path_length (BEzier2 integration) --- */

  c1 = max(chi[k] - (dsup*0.5) * dchi_c + chi[k-dk] + (dsup*0.5) * dchi_up,  0.0);
  dtau_uw =  dsup * (chi[k] + chi[k-dk] + c1*0.5) / 3.0;

  
  /* --- Ku, K0 and dKu, dSu --- */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);

  Svec(k_start,    S, Su);
  Svec(k_start+dk, S, S0);

  for(n=0;n<4;n++)
    for(m=0;m<4;m++) dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {

    /* ---  dchi/ds at downwind point --- */

    dsdn = fabs(z[k+dk] - z[k]   ) * imu;

    if(fabs(k-k_end)>1){
      dsdn2=fabs(z[k+2*dk] - z[k+dk]) * imu;
      dchi_dn = cent_deriv(dsdn,dsdn2,chi[k],chi[k+dk],chi[k+2*dk]);       
    } else dchi_dn=(chi[k+dk]-chi[k])/dsdn;
    
    
    /* --- Make sure that c1 and c2 don't do below zero --- */
    
    c1 = max(chi[k] + (dsdn*0.5) * dchi_c +chi[k+dk] - (dsdn*0.5) * dchi_dn, 0.0);

    
    /* --- Bezier3 integrated dtau --- */
    
    dtau_dw =  dsdn * (chi[k] + chi[k+dk] + c1 * 0.5) / 3.0;
    dt = dtau_uw;

    
    /* --- Bezier2 coeffs. ---- */

    dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt2 * dt2, dt05 = dt * 0.5;
    
    if(dt < 1.e-2){
      eps = 1.0 - dt + 0.5 * dt2 - dt3 / 6.0 + dt4 / 24.0;
      alpha = dt/3.0 - dt2/12.   + dt3/60.0;
      beta  = dt/3.0 - dt2*0.25  + dt3*0.1;
      gamma = dt/3.0 - dt2/6.0   + dt3*0.05;
    }else{
      eps = exp(-dt);
      alpha = (dt2 - 2.0*dt + 2.0 - 2.0*eps)    / dt2;
      beta  = (2.0 - (2.0 + 2.0*dt + dt2)*eps)  / dt2;
      gamma = (2.0*dt - 4.0 + (2.0*dt+4.0)*eps) / dt2;
    }

    
    /* --- We use the average of 2 control points so 
       include the 0.5 factor here --- */
    
    gamma *= 0.5;


    
    /* --- Diagonal operator --- */

    if(Psi) Psi[k] = alpha + gamma;
    
    
    /* ---- get algebra in place --- */
    
    StokesK(nspect, k+dk, chi[k+dk], Kd);
    Svec(k+dk, S, Sd);
    cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
    cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);
    
    m4m(Ku, Ku, Ma); // Ku # Ku
    m4m(K0, K0, A ); // K0 # K0
    
    for(j=0;j<4;j++){
      for(i=0;i<4;i++){
	A[j][i] = ident[j][i] + alpha * K0[j][i] - gamma *
	  (-(dt05 * (A[j][i] + dK0[j][i] + K0[j][i]) + K0[j][i]));

	Ma[j][i] = eps * ident[j][i] - beta * Ku[j][i] + gamma *
	  (dt05 * (Ma[j][i] + dKu[j][i] + Ku[j][i]) - Ku[j][i]);

	Mb[j][i] = beta * ident[j][i] + gamma * (ident[j][i] - dt05 * Ku[j][i]);
	Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] + dt05 * K0[j][i]);
      }
    }

    /* ---
       Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
       (gam * dS0 - gamma * dSu) * dtau / 5.0
       --- */
    
    for(i = 0; i<4; i++){
      V0[i] = 0.0;
      for(j = 0; j<4; j++){
	V0[i] += Ma[i][j] * I[j][k-dk] + Mb[i][j] * Su[j] + Mc[i][j] * S0[j];
      }
      V0[i] += dt05 * gamma * (dS0[i] - dSu[i]);
    }

    
    /* --- Solve linear system to get the intensity --- */

    memcpy(Md[0], A[0], 16*sizeof(double));
    SolveLinearEq(4, Md, V0, TRUE);
    
    for(i=0;i<4;i++) I[i][k] = V0[i];
    
    
    /* --- Shift values for next depth --- */
    
    memcpy(Su,   S0, 4*sizeof(double));
    memcpy(S0,   Sd, 4*sizeof(double));
    memcpy(dSu, dS0, 4*sizeof(double));
    
    memcpy(Ku[0],   K0[0], 16*sizeof(double));
    memcpy(K0[0],   Kd[0], 16*sizeof(double));
    memcpy(dKu[0], dK0[0], 16*sizeof(double));
    
    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;
  } // k loop
  
  /* --- Linear integration in the last interval --- */
  
  k = k_end;
  w3(dtau_uw, w);
  
  for (n = 0;  n < 4;  n++) V0[n] = w[0]*S[n][k] + w[1] * dSu[n];
  if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
  
  for (n = 0;  n < 4;  n++) {
    for (m = 0;  m < 4;  m++) {
      A[n][m] = -w[1]/dtau_uw * Ku[n][m];
      Md[n][m] = (w[0] - w[1]/dtau_uw) * K0[n][m];
    }
    A[n][n] = 1.0 - w[0];
    Md[n][n] = 1.0;
  }
  
  for (n = 0;  n < 4;  n++) 
    for (m = 0;  m < 4;  m++) 
      V0[n] += A[n][m] * I[m][k-dk];
  
  SolveLinearEq(4, Md, V0, TRUE);
  
  for (n = 0;  n < 4;  n++) I[n][k] = V0[n];
  
  freeMatrix((void **) Md);
}
