

/* ------- file: -------------------------- rhf1d.c -----------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 16:40:14 2011 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with Feautrier difference scheme
       in static atmospheres, and with piecewise quadratic integration
       in moving atmospheres.

       --                                              -------------- */

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "math.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "rhf1d.h"
#include "dummyatmos.h"
#include "sortlambda_j.h"
#include "initial_j.h"
#include "scatter_j.h"
#include "iterate_j.h"

/* --- Function prototypes --                          -------------- */

//extern void hermitian_interpolation(int n, double *x, double *y, int nn, double *xp, double *yp);


/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];
rhinfo io;
BackgroundData bgdat;
rhbgmem *bmem; // To store background opac in mem
crhpop *save_popp;
MPI_t mpi;

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* ------- begin -------------------------- rhf1d.c ----------------- */

bool_t rhf1d(float muz, int rhs_ndep, double *rhs_T, double *rhs_rho, 
	     double *rhs_nne, double *rhs_vturb, double *rhs_v, 
	     double *rhs_B, double *rhs_inc, double *rhs_azi,
	     double *rhs_z, double *rhs_nhtot, double *rhs_tau ,
	     double *rhs_cmass, double gravity, bool_t stokes, ospec *sp,
	     crhpop *save_pop, int mynw, double *mylambda, int myrank, int savpop,
	     int iverbose, int *hydrostat)
{
  
  bool_t write_analyze_output, equilibria_only, quiet = ((iverbose == 0)? TRUE : FALSE);
  int    niter, nact, i, sNgperiod, sNgdelay, sPRDNITER;
  static int save_Nrays;
  static double save_muz, save_mux, save_muy, save_wmu;
  static enum StokesMode oldMode;

  double dpopmax;
  static Atom *atom;
  static Molecule *molecule;
 
  static bool_t firsttime = TRUE;
 
  int argc = 1;
  char *argv[] = {"rhf1d",NULL};

  if(firsttime){
    memset(&atmos, 0, sizeof(Atmosphere));
    memset(&geometry, 0, sizeof(Geometry));
    memset(&spectrum, 0, sizeof(Spectrum));
    atmos.cos_gamma = NULL;
    atmos.cos_2chi = NULL;
    atmos.sin_2chi = NULL;
  }
  
  atmos.moving = TRUE;
  atmos.Stokes = FALSE;
  atmos.Nspace = rhs_ndep;
  
  /* --- Read input data and initialize --             -------------- */

  
  setOptions(argc, argv, myrank, (int)quiet);
  if(iverbose == 0) commandline.quiet = quiet;
  
  
  if(firsttime){
    atmos.nH = NULL;
    getCPU(0, TIME_START, NULL);
    SetFPEtraps();
    geometry.Ndep = rhs_ndep;
    atmos.gravity = gravity;
    readInput();
    readAbundance(&atmos);
    DUMMYatmos(&atmos, &geometry, firsttime);
    oldMode = input.StokesMode;
    mpi.rank = myrank;
  }
  mpi.stop = false;


  /* --- Store Ng values --- */

  sNgperiod = input.Ngperiod;
  sNgdelay = input.Ngdelay;
  sPRDNITER = input.PRD_NmaxIter;
  
  /* --- Allocate space for arrays that define structure -- --------- */
  
  geometry.tau_ref = rhs_tau;//(double *) malloc(Ndep * sizeof(double));
  geometry.cmass   = rhs_cmass;//(double *) malloc(Ndep * sizeof(double));
  geometry.height  = rhs_z;//(double *) malloc(Ndep * sizeof(double));
  atmos.T      = rhs_T;//(double *) malloc(Ndep * sizeof(double));
  atmos.ne     = rhs_nne;//(double *) malloc(Ndep * sizeof(double));
  atmos.vturb  = rhs_vturb;//(double *) malloc(Ndep * sizeof(double));
  geometry.vel = rhs_v;//(double *) malloc(Ndep * sizeof(double));
  atmos.B = rhs_B;
  atmos.gamma_B = rhs_inc;
  atmos.chi_B = rhs_azi;
  atmos.H_LTE = TRUE;
  atmos.nHtot = rhs_nhtot;
  atmos.rho = rhs_rho;
  save_popp = save_pop;
  input.StokesMode = oldMode;
  spectrum.updateJ = TRUE;  
  getCPU(1, TIME_START, NULL);
  if (input.StokesMode > NO_STOKES)
    atmos.Stokes = TRUE;
  hydrostat[0] = (int)atmos.hydrostatic;


  /* --- Init atomic/molecular models, extra lambda 
     positions and background opac. Allocated only in the first call --- */
  
  if(firsttime){
    readAtomicModels();
    readMolecularModels();
    SortLambda_j(mynw, mylambda);

    init_Background_j();
    Initvarious();
    
    /* --- Save geometry values to change back after --    ------------ */
    
    save_Nrays = atmos.Nrays;   save_wmu = geometry.wmu[0];
    save_muz = geometry.muz[0]; save_mux = geometry.mux[0]; save_muy = geometry.muy[0];
  }
  
  firsttime = FALSE;

  
  /* --- Reallocate stuff and compute background opac --- */
  
  UpdateAtmosDep();
  Background_j(write_analyze_output=FALSE, equilibria_only=FALSE);
  //convertScales(&atmos, &geometry);
  
  //for(i=0;i<atmos.Nspace;i++) printf("[%3d] %e\n", i, geometry.tau_ref[i]);
  
  if(!mpi.stop){
    
    /* --- Init profiles, populations and scattering --- */
    
    getProfiles();
    initSolution_j( myrank, savpop);

    read_populations(save_pop,0);

    if((savpop == 0) && 1){
      input.Ngdelay = min(7,input.Ngdelay) ;
      input.Ngperiod = min(11,input.Ngperiod) ;
      //input.PRD_NmaxIter = min(3,input.PRD_NmaxIter) ;
    }

    // for(niter=0; niter<spectrum.nPRDlines; niter++)
    // fprintf(stderr,"prdline->frac[0][0]=%e\n", spectrum.PRDlines[niter]->frac[0][0]);
    
    initScatter();
    
    //getCPU(1, TIME_POLL, "Total Initialize");
    
    
    /* --- Solve radiative transfer for active ingredients -- --------- */
    
    Iterate_j(input.NmaxIter, input.iterLimit, &dpopmax);
    if(isnan(dpopmax) || isinf(dpopmax) || dpopmax < 0){
      mpi.stop = true;
    }

    input.Ngdelay=sNgdelay, input.Ngperiod = sNgperiod, input.PRD_NmaxIter = sPRDNITER;
    
    
    /* --- Adjust stokes mode in case we are running POLARIZATION_FREE --- */
    
    if(!mpi.stop){
      adjustStokesMode();
      niter = 0;
      
      while ((niter < input.NmaxScatter)) {
	if (solveSpectrum(FALSE, FALSE, 0, TRUE) <= input.iterLimit) break;
	niter++;
      }
    } else dpopmax = 1.0e13;
  } else dpopmax = 1.0e13;
  
  bool_t converged = dpopmax < 1.e-3;//input.iterLimit;

  /* --- Store populations if needed --- */

  if(savpop > 0 && converged)
    save_populations(save_pop);
  // exit(0);

  /* --- Compute output ray --- */
  if(converged){
    atmos.Nrays     = 1;
    geometry.Nrays  = 1;
    geometry.muz[0] = muz;
    geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
    geometry.muy[0] = 0.0;
    geometry.wmu[0] = 1.0;
    spectrum.updateJ = FALSE;
    
    calculateRay();
    
    
    
    /* --- Put back previous values for geometry  --- */
    
    atmos.Nrays     = save_Nrays;
    geometry.Nrays = save_Nrays;
    geometry.muz[0] = save_muz;
    geometry.mux[0] = save_mux;
    geometry.muy[0] = save_muy;
    geometry.wmu[0] = save_wmu;
    spectrum.updateJ = TRUE;


  }

  input.StokesMode = oldMode;

  
  /* --- Copy desired ray to output arrays---*/
  
  sp->nrays = atmos.Nrays;
  sp->nlambda = spectrum.Nspect;
  if(sp->I != NULL) free(sp->I);
  if(sp->Q != NULL) free(sp->Q);
  if(sp->U != NULL) free(sp->U);
  if(sp->V != NULL) free(sp->V);
  if(sp->lambda != NULL) free(sp->V);
  sp->I = calloc(spectrum.Nspect, sizeof(double));
  sp->Q = calloc(spectrum.Nspect, sizeof(double));
  sp->U = calloc(spectrum.Nspect, sizeof(double));
  sp->V = calloc(spectrum.Nspect, sizeof(double));
  sp->lambda = calloc(spectrum.Nspect, sizeof(double));

  memcpy(sp->lambda, &spectrum.lambda[0], spectrum.Nspect * sizeof(double));

  
  if(converged){
    memcpy(sp->I, spectrum.I[0], spectrum.Nspect * sizeof(double));
    if(atmos.Stokes && stokes){
      memcpy(sp->Q, spectrum.Stokes_Q[0], spectrum.Nspect * sizeof(double));
      memcpy(sp->U, spectrum.Stokes_U[0], spectrum.Nspect * sizeof(double));
      memcpy(sp->V, spectrum.Stokes_V[0], spectrum.Nspect * sizeof(double));	
    }
  }


  /* --- Deallocate n & nstar --- */
  
  for (nact = 0; nact < atmos.Natom; nact++) {
    atom = &atmos.atoms[nact];
    if(atom->n == atom->nstar){
      freeMatrix((void **) atom->nstar);
      atom->nstar = NULL;
      atom->n = NULL;
    }else{
      if (atom->nstar != NULL) freeMatrix((void **) atom->nstar);
      if (atom->n != NULL) freeMatrix((void **) atom->n);
      atom->nstar = NULL;
      atom->n = NULL;
    }
  }

  if(!quiet){
    fclose(commandline.logfile);
    //if(!mpi.stop)
    remove(commandline.logfileName);
    //else exit(0);
  }

  return converged;
}


void clean_saved_populations(crhpop *save_pop){

  int ii, kr;
  
  if(save_pop->nactive > 0){
    // fprintf(stderr,"clean_saved_populations: cleaning pops\n");
    for(ii = 0;ii < save_pop->nactive;ii++){
      
      free(save_pop->pop[ii].n);
      free( save_pop->pop[ii].ntotal);
      
      save_pop->pop[ii].ntotal = NULL;
      save_pop->pop[ii].n = NULL;
      
      if(save_pop->pop[ii].nprd > 0){
	
	for(kr=0;kr<save_pop->pop[ii].nprd; kr++){
	  // fprintf(stderr,"cleaning rho [%d] -> %p \n", kr, save_pop->pop[ii].line[kr]);
	  
	  free(save_pop->pop[ii].line[kr].rho);
	}

	free(save_pop->pop[ii].line);
	
      } // nprd > 0

      save_pop->pop[ii].nprd = 0;
      
    } // nactive
    
    
    free(save_pop->pop);
    free(save_pop->lambda);
    free(save_pop->J);
    if(input.backgr_pol) free(save_pop->J20);
    free(save_pop->tau_ref);
    
    save_pop->pop = NULL;
    save_pop->lambda = NULL;
    save_pop->J = NULL;
    save_pop->J20 = NULL;
    save_pop->tau_ref = NULL;
  }

  save_pop->nactive = 0;
}


void save_populations(crhpop *save_pop){

  Atom *atom;
  int    niter, nact, save_Nrays, nactotal, ii, kr, nprd;
  AtomicLine *line;
  register int k,j;

  
  
  // fprintf(stderr,"save_pop: nactive=%d\n", save_pop->nactive);
  
  if(save_pop->nactive > 0) clean_saved_populations(save_pop);

  save_pop->nactive = atmos.Nactiveatom;
  save_pop->ndep = atmos.Nspace;
  save_pop->nw = spectrum.Nspect;
  //
  save_pop->pop =     (crhatom*) malloc(atmos.Nactiveatom * sizeof(crhatom));
  save_pop->lambda =  (double*) malloc(spectrum.Nspect*sizeof(double));
  //
  save_pop->J =   (double*) malloc(spectrum.Nspect* atmos.Nspace*sizeof(double));
  save_pop->tau_ref = (double*) malloc(atmos.Nspace*sizeof(double));
  
  
  if(input.backgr_pol)
    save_pop->J20 = (double*) malloc(spectrum.Nspect* atmos.Nspace*sizeof(double));
  else
    save_pop->J20 = NULL;

  memcpy(save_pop->tau_ref, geometry.tau_ref, sizeof(double)*atmos.Nspace);
  
  
  /* --- copy J, J20 and lambda ---*/
  
  memcpy(&save_pop->J[0], &spectrum.J[0][0],
	 spectrum.Nspect*atmos.Nspace*sizeof(double));
  
  if(input.backgr_pol){
    memcpy(&save_pop->J20[0], &spectrum.J20[0][0],
	   spectrum.Nspect*atmos.Nspace*sizeof(double));
  }// else  fprintf(stderr, "NOT STORING J20!\n");
    
  memcpy(&save_pop->lambda[0], spectrum.lambda, spectrum.Nspect * sizeof(double));

  /* --- Copy populations for each active atom --- */
  
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    
    atom = atmos.activeatoms[nact];
    save_pop->pop[nact].nlevel = atom->Nlevel;
    save_pop->pop[nact].converged  = true;

    /* --- Allocate arrays ---*/
    save_pop->pop[nact].n = (double*) malloc( atom->Nlevel * atmos.Nspace*sizeof(double));
    save_pop->pop[nact].ntotal = (double*) malloc(atmos.Nspace*sizeof(double));


    /* --- copy populations and ntotal---*/
    // memcpy(&save_pop->pop[nact].n[0], &atom->n[0][0],
    //atom->Nlevel*atmos.Nspace*sizeof(double));
    memcpy(&save_pop->pop[nact].ntotal[0], &atom->ntotal[0], atmos.Nspace*sizeof(double));

    for(j=0;j<atom->Nlevel;j++)
      for(k=0;k<atmos.Nspace;k++)
	save_pop->pop[nact].n[j*atmos.Nspace+k] = atom->n[j][k] / atom->nstar[j][k];
    

    
    /* --- Check number of PRD lines in atom --- */
    
    nprd = 0;
    for (kr = 0;  kr < atom->Nline;  kr++)
      if (atom->line[kr].PRD) nprd++;
    

    /* --- Allocate arrays to store rho for each PRD line --- */
    
    save_pop->pop[nact].nprd = nprd;
    
    if(nprd >0){
      
      save_pop->pop[nact].line = (crhprd*) malloc(nprd * sizeof(crhprd));
      ii = 0;
      
      for (kr = 0;  kr < atom->Nline;  kr++){
	
	line = &atom->line[kr];
	
	if (line->PRD){
	  save_pop->pop[nact].line[ii].idx = kr;
	  save_pop->pop[nact].line[ii].nlambda = line->Nlambda;
	  save_pop->pop[nact].line[ii].rho = (double*)
	    malloc(line->Nlambda * atmos.Nspace * sizeof(double));
	  //
	  //  fprintf(stderr,"ii=%d, srho=%p, lrho=%p\n",
	  //	  ii, &save_pop->pop[nact].line[ii].rho[0], &line->rho_prd[0][0]);
	  //
	  memcpy(save_pop->pop[nact].line[ii].rho, &line->rho_prd[0][0],
		 line->Nlambda*atmos.Nspace*sizeof(double));
	  ii++;
	} // PRD line
      } // kr
    }// nprd > 0
    
  } // nact

  
}

void read_populations(crhpop *save_pop, int flag){
  Atom *atom;
  int    niter, nact, save_Nrays, nactotal, ii, nprd, kr, kkr, copied = 0;
  AtomicLine *line;
  register int k, j, la;
  double ratio, tmp;
  double *tmp1 = (double*)malloc(atmos.Nspace*sizeof(double));
  
  if(save_pop->pop == NULL || save_pop->nactive != atmos.Nactiveatom){
    // fprintf(stderr,"read_population: atmos->Nactiveatom[%d] != save_pop->nactive[%d]\n",
    //atmos.Nactiveatom, save_pop->nactive);
    return;
  }
  
  for(nact=0;nact < save_pop->nactive;nact++){
    atom = atmos.activeatoms[nact];
    
    /* --- Check dimensions --- */
    if(atom->Nlevel != save_pop->pop[nact].nlevel || atmos.Nspace != save_pop->ndep){
      continue;
    }
    copied += 1;
    
    
    for(j= 0; j < atom->Nlevel; j++ ){
      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->pop[nact].n[j*atmos.Nspace],
			      (int)atmos.Nspace, geometry.tau_ref, tmp1, 1);
      //memcpy(tmp1,&save_pop->pop[nact].n[j*atmos.Nspace],atmos.Nspace*sizeof(double));
      
      
      for(k = 0; k < atmos.Nspace; k++){
	atom->n[j][k] = tmp1[k] * atom->nstar[j][k];
      }
    }
    
    /* --- Make sure that the populations do not exceed ntotal --- */
    
    
    for(k = 0; k < atmos.Nspace; k++){
      ratio = 0.0;
      for(j= 0; j < atom->Nlevel; j++ )
	ratio += atom->n[j][k];
      
      ratio = atom->ntotal[k] / ratio;
      for(j= 0; j < atom->Nlevel; j++ )
	atom->n[j][k] *= ratio;
    }
    
    
    /* --- Copy rho for PRD lines? --- */
    
    if(save_pop->pop[nact].nprd > 0){
      
      /* --- Number of PRD lines in atom --- */
      nprd = 0;
      for (kr = 0;  kr < atom->Nline;  kr++){
	line = &atom->line[kr];
	if (line->PRD) nprd++;
      }

      if(nprd == save_pop->pop[nact].nprd){
	for(kr = 0; kr<save_pop->pop[nact].nprd; kr++){
	  
	  line = &atom->line[save_pop->pop[nact].line[kr].idx];
	  
	  if(line->PRD && (save_pop->pop[nact].line[kr].nlambda == line->Nlambda)){
	    // //   fprintf(stderr, "Copying rho array for nact=%d, line=%d\n", nact, kr);
	    //memcpy(&line->rho_prd[0][0], save_pop->pop[nact].line[kr].rho,
	    //		    line->Nlambda*atmos.Nspace*sizeof(double));
	    for(la=0;la<line->Nlambda;la++){
	      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->pop[nact].line[kr].rho[la*atmos.Nspace],
	    			      (int)atmos.Nspace, geometry.tau_ref, line->rho_prd[la],0);

	    }
	  }else{
	    fprintf(stderr,"read_populations: BAD BOOK-KEEPING, not a PRD line, not copying rho, idx=%d, kkr=%d \n", kr, kkr );

	  }
	  
	}
	
      }else{
	fprintf(commandline.logfile,"read_populations: nprd[%d] != save_pop.pop.nprd[%d], not copying PRD rho!\n", nprd,save_pop->pop[nact].nprd );
      }
      
    }//else{
     // fprintf(stderr,"read_populations, no PRD lines for ATOM=%d\n",nact);
    //}
  } // nact


  /* --- Copy radiation field --- */
  
  if(spectrum.Nspect == save_pop->nw || atmos.Nspace == save_pop->ndep){
    for(la=0;la<spectrum.Nspect;la++)
      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->J[la*atmos.Nspace],
      			      (int)atmos.Nspace, geometry.tau_ref, spectrum.J[la],0);
      //memcpy(spectrum.J[la], &save_pop->J[la*atmos.Nspace], atmos.Nspace*sizeof(double));
      
      
    if(input.backgr_pol){
      hermitian_interpolation((int)atmos.Nspace, save_pop->tau_ref, &save_pop->J20[la*atmos.Nspace],
      			      (int)atmos.Nspace, geometry.tau_ref, spectrum.J20[la],0);
      //memcpy(spectrum.J20[la], &save_pop->J20[la*atmos.Nspace], atmos.Nspace*sizeof(double));
	    
      //   memcpy(&spectrum.J20[0][0], &save_pop->J20[0],
      //	     spectrum.Nspect*atmos.Nspace*sizeof(double));
      
    }//else fprintf(stderr, "NOT READING J20!\n");
  }

  free(tmp1);
}




/* ------- end ---------------------------- rhf1d.c ----------------- */
