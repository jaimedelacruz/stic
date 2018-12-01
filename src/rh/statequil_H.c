/*
  
  The following routines are adapted from M3D to compute the statistical equilibrium
  equations and change conservation at the same time for H when it is an active species.
  The contribution to Ne of all other species is assumed to be LTE. H- is ignored in 
  this calculation as it has a very small impact and it makes the problem more non-linear.
  
  Coded in C by J. de la Cruz Rodriguez (ISP-SU 2018)
  Original routines in Fortran by J. Bjoergen (ISP-SU 2018)
  
  The downward radiative rates from the continuum are implicitly multiplied by the initial Ne.
  We remove that contribution and multiply it with the current value. The collisional terms
  must be recomputed too. Background opacities and profiles must be recomputed while Ne is being 
  iterated.

  In our tests the most stable solution is to start from LTE, not ZERO_RADIATION because
  Ne is usually computed in LTE in most simulations/inversions/EOS.

  NOTES: Do not use SVD to solve the Newton-Raphson linear system, it yields garbage.

 */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "accelerate.h"
#include "statistics.h"
#include "error.h"
#include "rh_1d/rhf1d.h"
#include "inputs.h"
#include "constant.h"
#include <string.h>
#include "statequil_H.h"

/* --- Function prototypes --                          -------------- */

void getfjk2(Element *element, double ne, int k, double *fjk, double *dfjk);
double getKuruczpf2(Element *element, int stage, int k);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];
extern MPI_t mpi;
extern InputData input;


inline double getKuruczpf2(Element *element, int stage, int k)
{
  bool_t hunt = TRUE;
  double Uk;
  Linear(atmos.Npf, atmos.Tpf, element->pf[stage], 
	 1, &atmos.T[k], &Uk, hunt);

  return Uk;
}
// -------------------------------------------------------------------//

void getfjk2(Element *element, double ne, int k, double *fjk, double *dfjk)
{
  register int i, j;

  double  sum1, sum2, CT_ne, Uk, Ukp1;
  Atom *atom;

  /* --- Get the fractional population f_j(ne, T) = N_j/N for element
         element and its partial derivative with ne. -- ------------- */



    /* --- Else use estimate from LTE from Kurucz partition
           functions --                                -------------- */

    static const double C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

    CT_ne   = 2.0 * pow(C1/atmos.T[k], -1.5) / ne;
    sum1    = 1.0;
    sum2    = 0.0;
    fjk[0]  = 1.0;
    dfjk[0] = 0.0;

    Uk = getKuruczpf2(element, 0, k);
    for (j = 1;  j < element->Nstage;  j++) {
      Ukp1 = getKuruczpf2(element, j, k);
      
      fjk[j]  = fjk[j-1] * CT_ne *
	exp(Ukp1 - Uk - element->ionpot[j-1]/(KBOLTZMANN*atmos.T[k]));
      dfjk[j] = -j * fjk[j] / ne;
      sum1   += fjk[j];
      sum2   += dfjk[j];
      Uk      = Ukp1;
    }

    for (j = 0;  j < element->Nstage;  j++) {
      fjk[j]  /= sum1;
      dfjk[j]  = (dfjk[j] - fjk[j] * sum2) / sum1;
    }
    
}


 double sum1(int n, int k, double **m){
  register int ii;
  double sum = 0.0;
  for(ii=0;ii<n;++ii) sum+= m[ii][k];
  return sum;
}

 double sum2(int n, int k, double **m){
  register int ii;
  double sum = 0.0;
  for(ii=0;ii<n;++ii) sum+= m[k][ii];
  return sum;
}
// -------------------------------------------------------------------//


void extractGamma(Atom *atom, double **G, int k)
{
  register int i, j, ij, kr;
  int nlev = atom->Nlevel;
  AtomicContinuum *cont;

  
  // --- Radiative contribution to gamma (preconditioned) --- //
  
  for (i = 0, ij = 0;  i < nlev;  i++) 
    for (j = 0;  j < nlev;  j++, ij++)
      G[i][j] = atom->Gamma[ij][k];
  
  
  // --- divide by "ne" the downward rates from the continuum ---//
  
  for(kr=0; kr<atom->Ncont; ++kr){
    cont = &atom->continuum[kr];
    i = cont->i;
    j = cont->j;
    
    G[i][j] /= atmos.ne[k];
  }
  
}


// -------------------------------------------------------------------//

void obtain_rates(double **w, double **c, int k, Atom *atom)
{
  register int kr, i, j, ij;
  int nlev = atom->Nlevel;
  AtomicLine *line;
  AtomicContinuum *cont;
  FixedTransition *fixed;
  
  // --- Set to Zero ---//
  
  memset(c[0], 0, (nlev*nlev)*sizeof(double));

  
  // --- Now loop through nlines and ncont and get rates --- //
  // --- In RH, in the rate matrix Gamma, the fast index (rightmost)
  // --- stores the first index of the Einstein coeffs, to the
  // --- notation is reversed: Gamma[i][j] -> Rji, Aji;
  // --- The collisional terms are already in the correct order
  
  for(kr=0; kr<atom->Ncont; ++kr){
    cont = &atom->continuum[kr];
    i = cont->i;
    j = cont->j;

    
    // --- Multiply the continuum downward transitions by
    // ---- the current estimate of Ne
    
    w[i][j] *= atmos.ne[k]; 
  }
  
  // --Recompute collisional rates for one depth-point --- //
  
  CollisionRateOne(atom, atom->fp_input, k);


  // --- read collisional rates --- //
  
  for (i = 0, ij = 0;  i < nlev;  i++) 
    for (j = 0;  j < nlev;  j++, ij++)
      c[i][j] = atom->C[ij][k];
}


// -------------------------------------------------------------------//

void rate_eq(double **dFF, double *F, double *var, double **ww, Atom *atom)
{
  register int ii, jj, kk;
  double sum;
  int nlev = atom->Nlevel;
  
  for(ii = 0; ii<nlev; ++ii){
    sum = sum1(nlev, ii, ww);
    
    F[ii] = var[ii] * sum;
    
    for(jj=0; jj<nlev; ++jj){
      F[ii] -= var[jj] * ww[ii][jj];
      
      if(ii == jj){
	dFF[ii][jj] = sum;
      }else{
	dFF[ii][jj] = -ww[ii][jj];
      } //else
      
    } //jj
  } // ii
}

// -------------------------------------------------------------------//

void charge_eq(double **dFF, double *F, double *var, double **w,  double **c,
	       Atom *atom, int k)
{

  register int  n, j, i, ij, ii, jj;
  int npar, anl, nn,  Nmaxstage = 0, anl1;

  double  CT_ne, fsum, dfsum,Uk,Ukp1,error,derror,nHtot,akj;
  Atom *H = &atmos.atoms[0];

  double *fjk = NULL, *dfjk = NULL,  PhiHmin, tg, sum, sum0;
  static const double C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);
  
  nHtot = atmos.H->ntotal[k];
  tg    = atmos.T[k];

  npar = atom->Nlevel + 1;
  anl  = atom->Nlevel, anl1 = anl-1;
  
  for(n=0; n < anl1; ++n) dFF[npar-1][n] = 0.0;
  
  dFF[npar-1][anl-1]  = -1.0/var[npar-1];
  dFF[npar-1][npar-1] = var[anl-1]/SQ(var[npar-1]);

  F[npar-1] = 1.0 - var[anl-1]/var[npar-1];

  Nmaxstage = 0;
  
  for (n = 1;  n < atmos.Nelem;  n++) 
    Nmaxstage = MAX(atmos.elements[n].Nstage, Nmaxstage);
  
  fjk  = (double *) calloc(Nmaxstage , sizeof(double));
  dfjk = (double *) calloc(Nmaxstage , sizeof(double));

  
  for (n = 1;  n < atmos.Nelem;  n++) {
    
    getfjk2(&atmos.elements[n], var[npar-1], k, fjk, dfjk);

    for(j=1; j<atmos.elements[n].Nstage; ++j){
      akj = atmos.elements[n].abund * j * nHtot;
      F[npar-1] -= akj*fjk[j] / var[npar-1];
      
      dFF[npar-1][npar-1] += akj*fjk[j]/SQ(var[npar-1]) - akj/var[npar-1]*dfjk[j];    
    }
    
  }

  free(fjk);
  free(dfjk);
    
  
  
  // --- Derivatives of Ne --- //

  for(ii=0; ii<anl; ++ii){
    
    if(ii == (anl-1)){

      sum = sum1(anl, anl-1, w)/var[npar-1];
      sum0= sum1(anl, anl-1, c)/var[npar-1];
      
      dFF[anl-1][npar-1] = var[anl-1]*(sum + 2.0*sum0);

      for(jj=0;jj<anl; ++jj)
	dFF[anl-1][npar-1] -= (var[jj]/var[npar-1])*c[anl-1][jj];
      
    }else{
      // sum = 0.0; for(jj=0;jj<anl; ++jj) sum += c[jj][ii];
      sum = sum1(anl, ii, c) / var[npar-1];
      dFF[ii][npar-1] = var[ii]*sum;

      for(jj=0; jj<anl1; ++jj)
	dFF[ii][npar-1] -= var[jj]*(c[ii][jj]/var[npar-1]);

      // --- recombination --- //

      dFF[ii][npar-1] += (var[ii] / var[npar-1]) * (c[anl-1][ii]) -
	(var[anl-1]/var[npar-1])*(w[ii][anl-1]+2.0*c[ii][anl-1]);      
    } // else
  } // ii

  //  freeMatrix((void**)c);
}

// -------------------------------------------------------------------//

void particle_cons(double **dFF, double *F, double *var, double tot, int nlev)
{
  register int kr;
  int pos1 = 0, npar = nlev+1;
  double maxvar = var[0];
  
  for(kr = 1; kr<nlev; ++kr)
    if(var[kr] > maxvar){
      maxvar = var[kr];
      pos1 = kr;
    }

  F[pos1] = 1.0;
  for( kr = 0; kr<nlev; ++kr) F[pos1] -= var[kr]/tot;
  
  dFF[pos1][npar-1] = 0.0;
  
  for(kr = 0; kr<nlev; ++kr) dFF[pos1][kr] = -1.0 / tot;
}


// -------------------------------------------------------------------//

void newtrap(double **dFF, double *F, double *var, int nlev)
{
  register int kr, ii;
  int npar = nlev+1;

  for(kr=0; kr<npar; ++kr){
   F[kr] *= -1.0; 
    for(ii=0; ii<npar; ++ii)
      dFF[kr][ii] *= var[ii]; // Multiply by var so we get a relative correction to n_k and ne
  }

  SolveLinearEq(npar, dFF, F, TRUE);
  
}

// -------------------------------------------------------------------//

void statEquil_H(Atom *atom, int isum, int mali_iter)
{
  static const double TOLX = 1.e-5;
  static const double vdamp = 10.0;
  static const double tiny_ne = 10.e6; // in M^{-3}, 10.in CGS.
  static const double tmax_update = 300000.0; 
  int max_iter;

  
  register int i, j, ij, k, kr, ii, jj;

  int    i_eliminate, Nlevel, niter, npar, Nlevel2, nlte_positive;
  double GamDiag, nmax_k, **Gamma_k, sum;
  double  *var, **Gamma_total, e, *npre,nepre, vdamp_var, tg, gradient_zero;
  double *RHS, **dFF, maxabs, tmp, max_error, **c, **Gamma_rad;
  int rst_vdamp = 0, lte_positive = 1, first = 1, *converged = NULL;
  
  
  getCPU(3, TIME_START, NULL);
  
  Nlevel = atom->Nlevel;
  Nlevel2=Nlevel*Nlevel;
  npar = Nlevel + 1;

  // --- Init storage --- //
  
  converged = (int*) calloc(atmos.Nspace,sizeof(int));
  
  Gamma_k = matrix_double(Nlevel, Nlevel);
  Gamma_rad   = matrix_double(Nlevel, Nlevel);

  Gamma_total = matrix_double(Nlevel, Nlevel);
  c  = matrix_double(Nlevel, Nlevel);
  npre = malloc(Nlevel*sizeof(double));
  
  var  = malloc(npar*sizeof(double));
  RHS  = malloc(npar*sizeof(double));
  dFF = matrix_double(npar, npar);


  
  for (k = 0;  k < atmos.Nspace;  k++) {
    
    if(atmos.T[k] > tmax_update){
      converged[k] = 1;
      continue;
    }

    extractGamma(atom, Gamma_rad, k);
    
    
    niter = 0;
    nepre = atmos.ne[k], vdamp_var = vdamp,
      tg = atmos.T[k], gradient_zero = 1.0;
    rst_vdamp = 0, lte_positive = 1, first = 1;
    max_iter = 1000;
    
    for(j=0; j<atom->Nlevel; ++j)
      npre[j] = atom->n[j][k];
    
    // --- Iterate --- //
    
    while(niter < max_iter){
      memset(RHS, 0, npar*sizeof(double));
      memset(&dFF[0][0], 0, npar*npar*sizeof(double));
      memcpy(Gamma_k[0], Gamma_rad[0], Nlevel2*sizeof(double));
      
      
      for(j=0; j<atom->Nlevel; ++j)
	var[j] = atom->n[j][k];
      
      var[npar-1] = atmos.ne[k];
      
      
      // --- Get rates --- //
      
      obtain_rates(Gamma_k, c,  k, atom);

      
      // --- add collisional rates to Gamma --- //

      for(i=0;i<Nlevel;++i)
	for(j=0;j<Nlevel;++j)
	  Gamma_total[i][j] = Gamma_k[i][j] + c[i][j];
      
      
      // --- Rate Eq. --- //
      
      rate_eq(dFF, RHS, var, Gamma_total, atom);


      // --- Charge Eq. --- //

      charge_eq(dFF, RHS, var, Gamma_k, c, atom, k);
      
      
      // --- Particle conservation --- //

      particle_cons(dFF, RHS, var, atom->ntotal[k], Nlevel);


      // --- Newton-Raphson --- //

      newtrap(dFF, RHS, var, Nlevel);
      
      
      // --- update n_k and ne_k with relative correction --- //
      
      for(kr=0; kr<Nlevel; ++kr)
	atom->n[kr][k] *= (1.0 + RHS[kr]/(1.0+vdamp_var*fabs(RHS[kr]))); 

      atmos.ne[k] *= (1.0 + RHS[npar-1]/(1.0+vdamp_var*fabs(RHS[npar-1])));


      
      // --- update parameters --- //
            
      nlte_positive = 1;
      for(j=0; j<Nlevel; ++j){
	if(atom->n[j][k]     < 0.0) nlte_positive = 0;
      }
    
      
      if((atmos.ne[k] < tiny_ne) || !nlte_positive){
	rst_vdamp = 1;
      }
      
      if(rst_vdamp && vdamp_var < (vdamp+1.0)){
	for(j=0; j<Nlevel; ++j)
	  var[j] = npre[j];
	var[npar-1] = nepre;
	
	vdamp_var = vdamp*3.0;
	rst_vdamp = 0;
	niter = 0;
	max_iter = 3000;
	lte_positive = 1;
      }else if(rst_vdamp){
	converged[k] = 0;
	break;
      }
      
      
      maxabs = fabs(RHS[0]);
      for(j=1; j<npar; ++j){
	tmp = fabs(RHS[j]);
	maxabs = ((tmp > maxabs) ? tmp : maxabs);
      }
      
      if((maxabs < TOLX) && (!rst_vdamp)){
	converged[k] = 1;
	break;
      }
      
      ++niter;
    } // while


  } // k
  
  //for(k=0;k<atmos.Nspace;k++){
    //if(!converged[k]) fprintf(stderr,"statEquil_H: NR iterations not converged [%d]\n", k);
    //}

  freeMatrix((void **) Gamma_k);
  freeMatrix((void **) dFF);
  freeMatrix((void **) Gamma_total);
  freeMatrix((void **) c);
  freeMatrix((void **) Gamma_rad);

  free(RHS);
  free(npre);
  free(var);
  free(converged);
  
  
  getCPU(3, TIME_POLL, "Stat Equil");
}
