#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "accelerate.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "spectrum.h"
#include "rhf1d.h"

/* --- Global variables --                             -------------- */
extern Atmosphere atmos;
extern Geometry geometry;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
extern Spectrum spectrum;
extern rhinfo io;

void Initvarious(){
  int    nact;
  Atom  *atom;
  


  io.atom_file_pos = (long *) malloc(atmos.Nactiveatom * sizeof(long));
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    io.atom_file_pos[nact] = ftell(atom->fp_input);
  }
}



void UpdateAtmosDep(void) {
/* Updates the atmos-dependent factors for the atoms and molecules */
  const char routineName[] = "UpdateAtomsDep";
  int       ierror, nact, k, kr, la, Nlamu;
  double    vtherm;
  Atom     *atom;
  Molecule *molecule;
  AtomicLine    *line;
  MolecularLine *mrt;

  /* Put back initial Stokes mode */
  // input.StokesMode = mpi.StokesMode_save;

  //  mpi.zcut_hist[mpi.task] = mpi.zcut;
  /* Recalculate magnetic field projections */
  if (atmos.Stokes) {
    if (atmos.cos_gamma != NULL) {
      freeMatrix((void **) atmos.cos_gamma);
      atmos.cos_gamma = NULL;
    }
    if (atmos.cos_2chi != NULL) {
      freeMatrix((void **) atmos.cos_2chi);
      atmos.cos_2chi = NULL;
    }
    if (atmos.sin_2chi != NULL) {
      freeMatrix((void **) atmos.sin_2chi);
      atmos.sin_2chi = NULL;
    }
    Bproject();
  }
  
  /* Update atmos-dependent atomic  quantities --- --------------- */
  for (nact = 0; nact < atmos.Natom; nact++) {
    atom = &atmos.atoms[nact];
    
    /* Reallocate some stuff (because of varying Nspace) */
    atom->ntotal = (double *) realloc(atom->ntotal, atmos.Nspace * sizeof(double));
    atom->vbroad = (double *) realloc(atom->vbroad, atmos.Nspace * sizeof(double));
    
    if (atom->nstar != NULL) 
      freeMatrix((void **) atom->nstar);
    atom->nstar = matrix_double(atom->Nlevel, atmos.Nspace);
      
  
    /* When H is treated in LTE, n is just a pointer to nstar,
       so we don't need to free it */
    if (nact == 0) {
      
      if (!atmos.H_LTE) {
	//  freeMatrix((void **) atom->n);
	atom->n = matrix_double(atom->Nlevel, atmos.Nspace);
      }else  atom->n = atom->nstar;
    } else {
      /* Only allocate n again for active or atoms with read populations */
      if ((atom->active)){ //|| (atom->popsFile)) {
	if (atom->n != NULL)
	  freeMatrix((void **) atom->n);
	atom->n = matrix_double(atom->Nlevel, atmos.Nspace);
	//	  printf("Reallocating n[%d,%d]\n", atom->Nlevel, atmos.Nspace);
	
      } else {
	/* alias to nstar again, just in case */
	atom->n = atom->nstar; 
      }
    }
    
    
    for (k = 0;  k < atmos.Nspace;  k++)
      atom->ntotal[k] = atom->abundance * atmos.nHtot[k];
    
    if (atom->Nline > 0) {
      vtherm = 2.0*KBOLTZMANN/(AMU * atom->weight);
      for (k = 0;  k < atmos.Nspace;  k++)
	atom->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));
    }
  }
  
  /* Now only for active atoms */
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    
    
    /* Rewind atom files to point just before collisional data */
    if ((ierror = fseek(atom->fp_input, io.atom_file_pos[nact], SEEK_SET))) {
      sprintf(messageStr, "Unable to rewind atom file for %s", atom->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    
    /* Reallocate some stuff (because of varying Nspace) */

    /* Free collision rate array, will be reallocated by calls in Background_p */
    if (atom->C != NULL) {
      freeMatrix((void **) atom->C);
      atom->C = NULL;
    }

    /* Allocate Gamma, as iterate released the memory */
    atom->Gamma = matrix_double(SQ(atom->Nlevel), atmos.Nspace);

    
    /* Initialise some continuum quantities */
    for (kr = 0; kr < atom->Ncont; kr++) {
      atom->continuum[kr].Rij = (double *) realloc(atom->continuum[kr].Rij,
						   atmos.Nspace * sizeof(double));
      atom->continuum[kr].Rji = (double *) realloc(atom->continuum[kr].Rji,
						   atmos.Nspace * sizeof(double));
      for (k = 0;  k < atmos.Nspace;  k++) {
	atom->continuum[kr].Rij[k] = 0.0;
	atom->continuum[kr].Rji[k] = 0.0;
      }
    }
    
    
    /* Initialise some line quantities */
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];
      
      if (line->phi  != NULL) {
	freeMatrix((void **) line->phi);
	line->phi = NULL;
      }
      if (line->wphi != NULL) {
	free(line->wphi);
	line->wphi = NULL;
      }
      
      if (atmos.moving && line->polarizable && (input.StokesMode >= FIELD_FREE)) {
	
	if (line->phi_Q != NULL) {
	  freeMatrix((void **) line->phi_Q);
	  line->phi_Q = NULL;
	}
	if (line->phi_U != NULL) {
	  freeMatrix((void **) line->phi_U);
	  line->phi_U = NULL;
	}
	if (line->phi_V != NULL) {
	  freeMatrix((void **) line->phi_V);
	  line->phi_V = NULL;
	}
	

	if (input.magneto_optical) {
	  if (line->psi_Q != NULL) {
	    freeMatrix((void **) line->psi_Q);
	    line->psi_Q = NULL;
	  }
	  if (line->psi_U != NULL) {
	    freeMatrix((void **) line->psi_U);
	    line->psi_U = NULL;
	  }
	  if (line->psi_V != NULL) {
	    freeMatrix((void **) line->psi_V);
	    line->psi_V = NULL;
	  }
	}
      }
	

      /* realloc because of varying Nspace */
      line->Rij = (double *) realloc(line->Rij, atmos.Nspace * sizeof(double));
      line->Rji = (double *) realloc(line->Rji, atmos.Nspace * sizeof(double));

      for (k = 0;  k < atmos.Nspace;  k++) {
	line->Rij[k] = 0.0;
	line->Rji[k] = 0.0;
      }	
     
      
      if (line->PRD) {
	if (line->Ng_prd != NULL) {
	  NgFree(line->Ng_prd);
	  line->Ng_prd = NULL;
	}
	
	if (line->fp_GII != NULL) {
	  fclose(line->fp_GII);
	  line->fp_GII = NULL;
	}

	if (input.PRD_angle_dep == PRD_ANGLE_DEP)
	  Nlamu = 2*atmos.Nrays * line->Nlambda;
	else
	  Nlamu = line->Nlambda;
	
	/* Idea: instead of doing this, why not free and just set line->rho_prd = NULL,
	   (and also line->Qelast?), because profile.c will act on that and reallocate */
	if (line->rho_prd != NULL) freeMatrix((void **) line->rho_prd);
	line->rho_prd = matrix_double(Nlamu, atmos.Nspace);
	
	if (line->Qelast != NULL) {
	  line->Qelast = (double *) realloc(line->Qelast, atmos.Nspace * sizeof(double));
	} else {
	  line->Qelast = (double *) malloc(atmos.Nspace * sizeof(double));
	}

	/* Initialize the ratio of PRD to CRD profiles to 1.0 */
	for (la = 0;  la < Nlamu;  la++) {
	  for (k = 0;  k < atmos.Nspace;  k++)
	    line->rho_prd[la][k] = 1.0;
	}

	// reset interpolation weights 
	if (input.PRD_angle_dep == PRD_ANGLE_APPROX) {
	  Nlamu = 2*atmos.Nrays * line->Nlambda;	  
	  if (line->frac != NULL) {
	    freeMatrix((void **) line->frac);
	    line->frac = NULL;
	  }
	  if (line->id0 != NULL){
	    freeMatrix((void **) line->id0);
	    line->id0 = NULL;
	  }
	  if (line->id1 != NULL) {
	    freeMatrix((void **) line->id1);
	    line->id1 = NULL;
	  }
	}

      }
    }
  }

  //distribute_nH();

  
  /* Update atmos-dependent molecular  quantities --- --------------- */
  for (nact = 0; nact < atmos.Nmolecule; nact++) {
    molecule = &atmos.molecules[nact];

    /* Reallocate some stuff, because of varying Nspace */
    molecule->vbroad = (double *) realloc(molecule->vbroad, atmos.Nspace * sizeof(double));
    molecule->pf     = (double *) realloc(molecule->pf,     atmos.Nspace * sizeof(double));
    molecule->n      = (double *) realloc(molecule->n,      atmos.Nspace * sizeof(double));
    if (molecule->nv != NULL) {
      freeMatrix((void **) molecule->nv);
      molecule->nv = matrix_double(molecule->Nv, atmos.Nspace);
    }
    if (molecule->nvstar != NULL) {
      freeMatrix((void **) molecule->nvstar);
      molecule->nvstar = matrix_double(molecule->Nv, atmos.Nspace);
    }
    if (molecule->pfv != NULL) {
      freeMatrix((void **) molecule->pfv);
      molecule->pfv = matrix_double(molecule->Nv, atmos.Nspace);
    }
    

    vtherm = 2.0*KBOLTZMANN / (AMU * molecule->weight);
    for (k = 0;  k < atmos.Nspace;  k++)
      molecule->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));
    
    if (molecule->active) {
      /* Allocate Gamma, as iterate released the memory */
      molecule->Gamma = matrix_double(SQ(molecule->Nv), atmos.Nspace);
      
      LTEmolecule(molecule);

      /* Free CO collision rate array, will be reallocated in initSolution_p */
      if (strstr(molecule->ID, "CO")) {
	free(molecule->C_ul);
	molecule->C_ul = NULL;
      }

      /* Free some line quantities */
      for (kr = 0;  kr < molecule->Nrt;  kr++) {
	mrt = &molecule->mrt[kr];
	if (mrt->phi  != NULL) {
	  freeMatrix((void **) mrt->phi);
	  mrt->phi = NULL;
	}
	if (mrt->wphi != NULL) {
	  free(mrt->wphi);
	  mrt->wphi = NULL;
	}
      }

    } else {
      if (molecule->Npf > 0) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  molecule->pf[k] = partfunction(molecule, atmos.T[k]);
      }
    }

  }


  return;
}

