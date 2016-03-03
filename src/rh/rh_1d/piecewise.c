/* ------- file: -------------------------- piecewise1d.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Feb 24 08:33:59 2012 --

       --------------------------                      ----------RH-- */

/* --- Piecewise quadratic integration of transfer equation in 
       one dimension.

    -- For boundary condition THERMALIZED use the relation

          I ~= B - \mu dB / d\tau,

       where \tau is taken in the direction of the ray (i.e. NOT in the
       sense of optical depth, e.g. see Mihalas, 1978, p. 51).
       --                                              -------------- */


#include <math.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


void har_mean_deriv(double* wprime,double dsup,double dsdn, 
		    double chiup,double chic,double chidn)
{

  double fim1, fi, alpha;

  fim1=(chic-chiup)/dsup;
  fi=(chidn-chic)/dsdn;

  if (fim1*fi > 0) {
    alpha = 0.33333333 * ( 1.0 + dsdn/(dsdn+dsup) );
    *wprime = (fim1*fi) / ( (1.0-alpha) * fim1 + alpha*fi );
  } else {
    *wprime=0.0;
  }

}


/* ------- begin -------------------------- Piecewise_1D.c ---------- */

void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  register int k;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, dS_dw, c1, c2, w[3],
         zmu, Bnu[2];

  zmu = 0.5 / geometry.muz[mu];

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
  dtau_uw = zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);
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
    }
  }
  I[k_start] = I_upw;
  if (Psi) Psi[k_start] = 0.0;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    w3(dtau_uw, w);

    if (k != k_end) {

      /* --- Piecewise quadratic here --               -------------- */ 

      dtau_dw = zmu * (chi[k] + chi[k+dk]) *
	fabs(geometry.height[k] - geometry.height[k+dk]);
      dS_dw   = (S[k] - S[k+dk]) / dtau_dw;


      /* Tiago: commenting this out to always keep linear piecewise
      c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw);
      c2 = (dS_uw - dS_dw);
      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] +
	(w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);
      */

      /* --- Try piecewise linear if quadratic gives negative
             monochromatic intensity --                -------------- */ 
      /*
      if (I[k] < 0.0) { */
      c1   = dS_uw;
      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] + w[1]*c1;

      if (Psi) Psi[k] = w[0] - w[1]/dtau_uw;
      
      /*
      } else {
	if (Psi) {
	  c1 = dtau_uw - dtau_dw;
	  Psi[k] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
	}
      }
      */
    } else {
	
      /* --- Piecewise linear integration at end of ray -- ---------- */

      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] + w[1]*dS_uw;
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }
    I_upw = I[k];

    /* --- Re-use downwind quantities for next upwind position -- --- */

    dS_uw   = dS_dw;
    dtau_uw = dtau_dw;
  }
}
/* ------- end ---------------------------- Piecewise_1D.c ---------- */

/* ------- begin -------------------------- PieceBezier_1D.c -------- */

void PieceBezier_1D(int nspect, int mu, bool_t to_obs,
		    double *chi, double *S, double *I, double *Psi)
{
  const char routineName[] = "PieceBezier_1D";
  register int k;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_uw, dS_dw, U[3],
    zmu, Bnu[2], psi_uw, psi_0, psi_dw, Sc;

  zmu = 0.5 / geometry.muz[mu];

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

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_uw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_uw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_uw = geometry.Ibottom[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_uw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu);
      I_uw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_uw = geometry.Itop[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[TOP]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  dtau_uw = zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);
  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;
 
  I[k_start]            = I_uw;
  if (Psi) Psi[k_start] = 0.0;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {
    U3(dtau_uw, U);

    dtau_dw = zmu * (chi[k] + chi[k+dk]) *
      fabs(geometry.height[k] - geometry.height[k+dk]);
    dS_dw   = (S[k] - S[k+dk]) / dtau_dw;

    if (dS_uw * dS_dw < 0.0) {
      psi_uw = U[0] + (U[2] - 2.0*dtau_uw * U[1]) / (dtau_uw * dtau_uw);
      psi_0  = (2.0*dtau_uw * U[1] - U[2]) / SQ(dtau_uw);

      I[k] = (1.0 - U[0]) * I_uw + psi_uw*S[k-dk] + psi_0*S[k];
    } else {
      Sc = S[k] - dtau_uw/2.0 * (dtau_dw*dS_uw + dtau_uw*dS_dw) /
	(dtau_dw + dtau_uw);

      if (Sc < MIN(S[k-dk], S[k])  ||  Sc > MAX(S[k-dk], S[k])) {
        psi_uw = U[0] - U[2] / (dtau_uw * dtau_uw);
        psi_0  = U[2] / (dtau_uw * dtau_uw);

        I[k] = (1.0 - U[0]) * I_uw + psi_uw*S[k-dk] + psi_0*S[k];
      } else {
        psi_uw  = U[0] + (U[2] - (dtau_dw + 2.0*dtau_uw) * U[1]) /
	  (dtau_uw * (dtau_uw + dtau_dw));
        psi_0   = ((dtau_dw + dtau_uw) * U[1] - U[2]) / (dtau_dw * dtau_uw);
        psi_dw  = (U[2] - dtau_uw * U[1]) / (dtau_dw * (dtau_dw + dtau_uw));

        I[k] = (1.0 - U[0]) * I_uw + psi_uw*S[k-dk] + psi_0*S[k] + 
	  psi_dw*S[k+dk];
      }
    }
    /* --- Re-use downwind quantities for next upwind position -- --- */

    I_uw    = I[k];
    dtau_uw = dtau_dw;
    dS_uw   = dS_dw;

    /* --- Diagonal operator --                        -------------- */

    if (Psi) Psi[k] = psi_0;
  }
  /* --- Far boundary --                               -------------- */

  U3(dtau_uw, U);
  
  psi_uw   = U[0] + (U[2] - 2.0*dtau_uw * U[1]) / SQ(dtau_uw);
  psi_0    = (2.0*dtau_uw * U[1] - U[2]) / SQ(dtau_uw);
  I[k_end] = (1.0 - U[0]) * I_uw + psi_uw*S[k_end-dk] + psi_0*S[k_end];

  if (Psi) Psi[k_end] = psi_0;
}
/* ------- end ---------------------------- PieceBezier_1D.c -------- */

void Piecewise_Hermite_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  register int k;

  int    k_start, k_end, dk, Ndep = geometry.Ndep,i;
  double dtau_uw, dtau_dw, dS_uw, I_upw, dS_dw, c1, c2, w[3],
         zmu, Bnu[2];
  double dsup,dsdn,dt,dt2,dt3,ex,alpha,beta,alphaprim,betaprim,dsup2;
  double dS_up,dS_c,dchi_up,dchi_c,dchi_dn,dsdn2,dtau_up2;

  zmu = 0.5 / geometry.muz[mu];

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
  dtau_uw = zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);
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
    }
  }
  I[k_start] = I_upw;
  if (Psi) Psi[k_start] = 0.0;
  
  /* set variables for first iteration to allow simple 
     shift for all next iterations */
  k=k_start+dk;
  dsup = fabs(geometry.height[k] - geometry.height[k-dk]) / geometry.muz[mu];
  dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) / geometry.muz[mu];
  dchi_up= (chi[k] - chi[k-dk])/dsup;
  /*  dchi/ds at central point */
  har_mean_deriv(&dchi_c,dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  /* upwind path_length */
  dtau_uw = 0.5 * dsup * (chi[k-dk]+chi[k])  + 1./12. * dsup*dsup * (dchi_up - dchi_c );
  /* dS/dtau at upwind point */
  dS_up=(S[k]-S[k-dk])/dtau_uw;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {

      /* downwind path length */
       dsdn = fabs(geometry.height[k+dk] - geometry.height[k]   ) / geometry.muz[mu];
       
      /* dchi/ds at downwind point */
       if (fabs(k-k_end)>1) {
	 dsdn2=fabs(geometry.height[k+2*dk] - geometry.height[k+dk])/geometry.muz[mu];
	 har_mean_deriv(&dchi_dn,dsdn,dsdn2,chi[k],chi[k+dk],chi[k+2*dk]);       
       } else {
	 dchi_dn=(chi[k+dk]-chi[k])/dsdn;
       }

      /* downwind optical path length */
       dtau_dw = 0.5 * dsdn * (chi[k]+chi[k+dk]) + 1./12.* dsdn * dsdn * (dchi_c  - dchi_dn) ;

      dt=dtau_uw;
      dt2=dt*dt;
      dt3=dt2*dt; 

      /* compute interpolation parameters */
      if (dt > 0.05) {
	ex=exp(-dt);
	alpha     =  ( 6.0*(dt-2.0)       + (-dt3+ 6.0*(dt+2.0))*ex  ) / dt3;
	beta      =  ( (dt3-6.0*(dt-2.0)) - 6.0*(dt+2.0)*ex          ) / dt3;
	alphaprim =  ( (2.0*dt-6.0)       + (dt2+4.0*dt+6.0)*ex      ) / dt2;
	betaprim  =  ( (-dt2+4.0*dt-6.0)  + (2.0*dt+6.0)*ex          ) / dt2;
      } else {
	ex=1.0 - dt + 0.5*dt2 - 0.16666667*dt3;
	alpha     =   (2.0/15.0)*dt3 - (7.0/20.0)*dt2 + dt/2.0;
	beta      =   (1.0/30.0)*dt3 - (3.0/20.0)*dt2 + dt/2.0;
	alphaprim = - (1.0/20.0)*dt3 + (1.0/12.0)*dt2;
	betaprim  =   (1.0/30.0)*dt3 - (1.0/12.0)*dt2;
      } 
          
      /* dS/dt at central point */
      har_mean_deriv(&dS_c,dtau_uw,dtau_dw,S[k-dk],S[k],S[k+dk]);
      
      I[k]= I_upw*ex + alpha*S[k-dk] + beta*S[k] + alphaprim*dS_up + betaprim*dS_c;
    
      if (Psi) Psi[k] = beta;
	
  } else { 
	
    /* --- Piecewise linear integration at end of ray -- ---------- */
    
      dtau_uw = zmu * (chi[k] + chi[k-dk]) *
	fabs(geometry.height[k] - geometry.height[k-dk]);
      dS_uw = (S[k] - S[k-dk]) / dtau_uw;
      w3(dtau_uw, w);
      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] + w[1]*dS_uw;
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }
    I_upw = I[k];
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    dS_up = dS_c;

  }
}
/* ------- end -------------------- Piecewise_Hermite_1D.c ---------- */
