#include <stdlib.h>
#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "rhf1d.h"



extern enum Topology topology;

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern Geometry geometry;
extern char messageStr[];
//extern NCDF_Atmos_file infile;
//extern MPI_data mpi;
//extern IO_data io; 
extern rhinfo io;



void calculateRay(void) {
  /* performs necessary reallocations and inits, and solves for ray */
  
  int i, nact, ierror, mu,k;
  bool_t analyze_output, equilibria_only,prdh_limit_mem_save;
  Atom *atom;
  AtomicLine *line;
  
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
  
  /* Must calculate background opacity for new ray, need some prep first */
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    
    /* Rewind atom files to point just before collisional data */
    if ((ierror = fseek(atom->fp_input, io.atom_file_pos[nact], SEEK_SET))) {
      sprintf(messageStr, "Unable to rewind atom file for %s", atom->ID);
      Error(ERROR_LEVEL_2, "rhf1df", messageStr);
    }
	
    /* Free collisions matrix, going to be re-allocated in background */
    if (atom->C != NULL) {
      freeMatrix((void **) atom->C);
      atom->C = NULL;
    }
	
    /* Recalculate line profiles for new angle */
    for (i = 0; i < atom->Nline; i++) {
      line = &atom->line[i];
	  
      if (line->phi  != NULL) {
	freeMatrix((void **) line->phi);
	line->phi = NULL;
      }
      if (line->wphi != NULL) {
	free(line->wphi);
	line->wphi = NULL;
      }
      
      if (atmos.moving && line->polarizable && (input.StokesMode>FIELD_FREE)) {
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
	  	  
      Profile(line);  
    }
  }
      
  /* reallocate intensities for correct number of rays */
  if (spectrum.I != NULL) freeMatrix((void **) spectrum.I);
  spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
  if (atmos.Stokes || input.backgr_pol) {
    if (spectrum.Stokes_Q != NULL) freeMatrix((void **) spectrum.Stokes_Q);
    spectrum.Stokes_Q = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (spectrum.Stokes_U != NULL) freeMatrix((void **) spectrum.Stokes_U);
    spectrum.Stokes_U = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (spectrum.Stokes_V != NULL) freeMatrix((void **) spectrum.Stokes_V);
    spectrum.Stokes_V = matrix_double(spectrum.Nspect, atmos.Nrays);
  }

  if (input.PRD_angle_dep == PRD_ANGLE_APPROX &&  atmos.NPRDactive > 0) {

    // recalculate line-of-sight velocity
    
    if (spectrum.v_los != NULL) freeMatrix((void **) spectrum.v_los);
    spectrum.v_los = matrix_double( atmos.Nrays, atmos.Nspace);
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	spectrum.v_los[mu][k] = vproject(k, mu); 
      }
    }
    
    
    /* set switch so that shift of rho_prd is done with a fresh
       interpolation */
    prdh_limit_mem_save = FALSE;
    if (input.prdh_limit_mem) prdh_limit_mem_save = TRUE;
    input.prdh_limit_mem = TRUE;
    
  }
      
  Background_j(analyze_output=FALSE, equilibria_only=FALSE);

  /* --- Solve radiative transfer for ray --           -------------- */
  solveSpectrum(FALSE, FALSE, 0, TRUE);

  //set back PRD input option
  if (input.PRD_angle_dep == PRD_ANGLE_APPROX && atmos.NPRDactive > 0)  
    input.prdh_limit_mem = prdh_limit_mem_save ;
        
}
