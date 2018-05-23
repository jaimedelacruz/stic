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
//#include "scl_opac.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];

/* --------------------------------------------------------------- */
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
/* --------------------------------------------------------------- */

void har_mean_deriv(double* wprime,double dsup,double dsdn, 
		    double chiup,double chic,double chidn)
{

  double fim1, fi, alpha;

  fim1=(chic-chiup)/dsup;
  fi=(chidn-chic)/dsdn;

  if (fim1*fi > 0) {
    alpha = 0.333333333333333 * ( 1.0 + dsdn/(dsdn+dsup) );
    *wprime = (fim1*fi) / ( (1.0-alpha) * fim1 + alpha*fi );
  } else {
    *wprime=0.0;
  }

}

/* inline double cent_deriv2(double dsup,double dsdn,  */
/* 			  double chiup,double chic,double chidn) */
/* { */

/*   double fim1, fi, alpha, wprime; */

/*   fim1=(chic-chiup)/dsup; */
/*   fi=(chidn-chic)/dsdn; */

/*   if (fim1*fi > 0) { */
/*     alpha = 0.3333333333333333333333 * ( 1.0 + dsdn/(dsdn+dsup) ); */
/*     wprime = (fim1*fi) / ( (1.0-alpha) * fim1 + alpha*fi ); */
/*   } else { */
/*     wprime=0.0; */
/*   } */
/*   return wprime; */
/* } */
int sign2(const double val)
{
  if(val >= 0.0) return 1.0;
  else           return -1.0;
}

double cent_deriv2(double odx,double dx, 
			 double chiup,double chic, double chidn)
{
  /* --- Derivatives from Steffen (1990) --- */

  double ody=(chic-chiup)/odx, dy=(chidn-chic)/dx, wprime = 0.0, p;
  
  
  if ((ody*dy) > 0) {
    p = (dy * odx + ody * dx) / (dx + odx);
    wprime = (2.0*sign2(dy)) * min(min(fabs(ody), fabs(dy)), 0.5*fabs(p));
  } 
  
  return wprime;
}

/* ------- begin -------------------------- Piecewise_1D.c ---------- */

void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  register int k, i;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, dS_dw, c1, c2, w[3],
         zmu, Bnu[2];

  double *z, *chi1;
  
  if(FALSE){
    chi1 = (double*)malloc(Ndep * sizeof(double));//scl_opac(Ndep, chi);
    for(i=0;i<Ndep;i++) chi1[i] = chi[i] / atmos.rho[i];
    z = geometry.cmass;//
  }else{
    chi1 = chi;
    z = geometry.height;
  }
  
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
  dtau_uw = zmu * (chi1[k_start] + chi1[k_start+dk]) *
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

      dtau_dw = zmu * (chi1[k] + chi1[k+dk]) *
	fabs(z[k] - z[k+dk]);
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

  if(chi1 != chi)
    free(chi1);
  
}
/* ------- end ---------------------------- Piecewise_1D.c ---------- */

/* ------- begin -------------------------- PieceBezier_1D.c -------- */

void PieceBezier_1D(int nspect, int mu, bool_t to_obs,
		    double *chi, double *S, double *I, double *Psi)
{
  const char routineName[] = "PieceBezier_1D";
  register int k, i;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_uw, dS_dw, U[3],
    zmu, Bnu[2], psi_uw, psi_0, psi_dw, Sc;

  zmu = 0.5 / geometry.muz[mu];
  
  double *z, *chi1;
  
  if(FALSE){
    chi1 = (double*)malloc(Ndep * sizeof(double));//scl_opac(Ndep, chi);
    for(i=0;i<Ndep;i++) chi1[i] = chi[i] / atmos.rho[i];
    z = geometry.cmass;//
  }else{
    chi1 = chi;
    z = geometry.height;
  }
  
  
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
  dtau_uw = zmu * (chi1[k_start] + chi1[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;
  
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

  dtau_uw = zmu * (chi1[k_start] + chi1[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;
 
  I[k_start]            = I_uw;
  if (Psi) Psi[k_start] = 0.0;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {
    U3(dtau_uw, U);

    dtau_dw = zmu * (chi1[k] + chi1[k+dk]) *
      fabs(z[k] - z[k+dk]);
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

  if(chi1 != chi)
    free(chi1);

  
}
/* ------- end ---------------------------- PieceBezier_1D.c -------- */

void Piecewise_Hermite_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  register int k, i;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, dS_dw, c1, c2, w[3],
         zmu, Bnu[2];
  double dsup,dsdn,dt,dt2,dt3,ex,alpha,beta,alphaprim,betaprim,dsup2;
  double dS_up,dS_c,dchi_up,dchi_c,dchi_dn,dsdn2,dtau_up2;
  //double *z = geometry.cmass;
  // double *chi1 = scl_opac(Ndep, chi);
  double *z, *chi1;
  
  if(FALSE){
    chi1 = (double*)malloc(Ndep * sizeof(double));//scl_opac(Ndep, chi);
    for(i=0;i<Ndep;i++) chi1[i] = chi[i] / atmos.rho[i];
    z = geometry.cmass;//
  }else{
    chi1 = chi;
    z = geometry.height;
  }
  
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
  dtau_uw = zmu * (chi1[k_start] + chi1[k_start+dk]) *
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
  dsup = fabs(z[k] - z[k-dk]) / geometry.muz[mu];
  dsdn = fabs(z[k+dk] - z[k]) / geometry.muz[mu];
  dchi_up= (chi1[k] - chi1[k-dk])/dsup;
  /*  dchi/ds at central point */
  har_mean_deriv(&dchi_c,dsup,dsdn,chi1[k-dk],chi1[k],chi1[k+dk]);
  /* upwind path_length */
  dtau_uw = 0.5 * dsup * (chi1[k-dk]+chi1[k])  + 1./12. * dsup*dsup * (dchi_up - dchi_c );
  /* dS/dtau at upwind point */
  dS_up=(S[k]-S[k-dk])/dtau_uw;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {

      /* downwind path length */
       dsdn = fabs(z[k+dk] - z[k]   ) / geometry.muz[mu];
       
      /* dchi/ds at downwind point */
       if (fabs(k-k_end)>1) {
	 dsdn2=fabs(z[k+2*dk] - z[k+dk])/geometry.muz[mu];
	 har_mean_deriv(&dchi_dn,dsdn,dsdn2,chi1[k],chi1[k+dk],chi1[k+2*dk]);       
       } else {
	 dchi_dn=(chi1[k+dk]-chi1[k])/dsdn;
       }

      /* downwind optical path length */
       dtau_dw = 0.5 * dsdn * (chi1[k]+chi1[k+dk]) + 1./12.* dsdn * dsdn * (dchi_c  - dchi_dn) ;

      dt=dtau_uw;
      dt2=dt*dt;
      dt3=dt2*dt; 

      /* compute interpolation parameters */
      if (dt > 0.05) {
	ex=((dt<12.0)?exp(-dt):0.0);
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
      
      //   fprintf(stderr,"[%3d] %e %e %e %e %e %e\n", k, ex, alpha, beta,alphaprim,betaprim, ex+alpha+beta+alphaprim+betaprim);

      
      I[k]= I_upw*ex + alpha*S[k-dk] + beta*S[k] + alphaprim*dS_up + betaprim*dS_c;
    
      if (Psi) Psi[k] = beta;
	
  } else { 
	
    /* --- Piecewise linear integration at end of ray -- ---------- */
    
      dtau_uw = zmu * (chi1[k] + chi1[k-dk]) *
	fabs(z[k] - z[k-dk]);
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
  if(chi1 != chi)
    free(chi1);

  //exit(0);
  
}
/* ------- end -------------------- Piecewise_Hermite_1D.c ---------- */

void Piecewise_lbrHermite_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  register int k, i;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, dS_dw, c1, c2, w[3],
         zmu, Bnu[2];
  double dsup,dsdn,dt,dt2,dt3,ex,alpha,beta,gamma,dsup2, obeta=0.0;
  double dS_up,dS_c,dchi_up,dchi_c,dchi_dn,dsdn2,dtau_up2;
  //double *z = geometry.cmass;
  // double *chi1 = scl_opac(Ndep, chi);
  double *z, *chi1, K0, K1, D, A, dir, E;
  double k_up, k_0, j_up, j_0, D_up, D_0;  
  
  if(FALSE){
    chi1 = (double*)malloc(Ndep * sizeof(double));//scl_opac(Ndep, chi);
    for(i=0;i<Ndep;i++) chi1[i] = chi[i] / atmos.rho[i];
    z = geometry.cmass;//
  }else{
    chi1 = chi;
    z = geometry.height;
  }
  
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
  dtau_uw = zmu * (chi1[k_start] + chi1[k_start+dk]) *
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
  dsup = -dk*(z[k] - z[k-dk]) / geometry.muz[mu];
  dsdn = -dk*(z[k+dk] - z[k]) / geometry.muz[mu];

  dchi_up = 0.0;
  dchi_c  = cent_deriv2(dsup,dsdn,chi1[k-dk],chi1[k],chi1[k+dk]);

  dtau_uw = 0.5 * dsup * (chi1[k-dk] + chi1[k])  +  dsup*dsup/12.0 * (dchi_up - dchi_c );
  dS_up=0.0;
  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {
      
      /* downwind path length */
      dsdn = -dk*(z[k+dk] - z[k]   ) / geometry.muz[mu];
      
      /* dchi/ds at downwind point */
      if (abs(k-k_end)>1){
	dsdn2=-dk*(z[k+2*dk] - z[k+dk]) / geometry.muz[mu];
	dchi_dn = cent_deriv2(dsdn,dsdn2,chi1[k],chi1[k+dk],chi1[k+2*dk]);       
      } else {
	dchi_dn = 0.0;
      }
      
      /* downwind optical path length */
      dtau_dw = 0.5 * dsdn * (chi1[k]+chi1[k+dk]) + dsdn *
	dsdn/12.0 * (dchi_c  - dchi_dn) ;
      dS_c    = cent_deriv2(dtau_uw,dtau_dw,S[k-dk],S[k],S[k+dk]);
      
    } else {
      dsdn = fabs(z[k] - z[k-dk]   ) / geometry.muz[mu];
      dchi_dn = 0.0, dS_c = 0.0;
      dsdn = dsup;
    }    
    
    dt = dtau_uw;
    K0 = 0.5 * dt;
    K1 = (dt*dt) * 0.08333333333333333333333333333333;
    
    
    D     =  1.0 + K0 + K1;
    A     = (1.0 - K0 + K1) / D;
    alpha =      (K0 - K1)  / D;
    beta  =      (K1 + K0)  / D;
    gamma =              K1 / D;
    
    E = alpha*S[k-dk] + beta*S[k] + gamma*(dS_up - dS_c);
    
    I[k]= I[k-dk]*A  + E;

    
    if (Psi) Psi[k] = beta + A*obeta; 

    I_upw = I[k];
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    dS_up = dS_c;
    obeta = beta;
  }
  
  if(chi1 != chi)
    free(chi1);
  //  exit(0);
}
