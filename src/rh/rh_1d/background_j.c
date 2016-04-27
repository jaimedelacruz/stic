/* ------- file: -------------------------- background.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jul 24 12:52:46 2013 --

       --------------------------                      ----------RH-- */

/* Driving subroutine for background opacity sources.

 * Included at the moment:

  ++ Thomson scattering by free electrons

  ++ Hydrogen:
    -- Bound-free absorption and emission by neutral Hydrogen
    -- Free-free absorption and emission by neutral Hydrogen
    -- Rayleigh scattering by neutral Hydrogen and Helium
    -- Rayleigh scattering by molecular Hydrogen (H2)
    -- Bound-free absorption and emission by Hminus (H^-)
    -- Free-free absorption and emission by Hminus (H + e)
    -- Free-free absorption and emission by H2minus (H2 + e)
    -- Free-free absorption and emission by H2plus (H + H^+)

  ++ Metals:
    -- Bound-free absorption and emission from metals specified in
       file background.input.
    -- Bound-bound absorption and emission from metals specified in
       file background.input.
    -- LTE Bound-bound absorption and emission all elements from
       a Kurucz line list.

  ++ Molecules:
    -- Chemical equilibrium is calculated for molecules specified
       in file background.input and populations of constituent atoms
       are reduced accordingly.
    -- molecular opacities (LTE) may be taken into account by specifying
       data files with transition lists in the molecular input files.
    -- Bound-free absorption and emission by OH and CH molecules.

 * Atomic models are specified in atoms.input, molecules in 
   molecules.input

 * Entries for the atoms.input and molecules.input files should have
   the form, respectively:

    -------------------------------------------------------------------
   |                                                                  |
   |   Nmetal                                                         |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION   population file  |
   |                          .                                       |
   |                          .                                       |
   |                                                                  |
   |   Nmolecule                                                      |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION                    |
   |                 .                                                |
   |                 .                                                |
    -------------------------------------------------------------------


   Nmetal and Nmolecule are the number of metal and molecule entries.
   metalID is the two-character atomID, the next entry is either
   set to LTE or NLTE, model_file is the input file containing atomic
   data for this metal (generic atomic input data format), and
   population_file is the input file containing the NLTE population
   numbers from a previous calculation. This last entry is only read when
   the second entry is set to NLTE. 

   -- Units:
      Wavelengths are given in nm, densities in m^-3, opacities in m^2,
      and emissivities in J s^-1 Hz^-1 sr^-1.

 Note: The model atom file for hydrogen is specified in keyword.input.
       If H_LTE = TRUE is specified there LTE hydrogen populations are
       used. See: distribute_nH in the file hydrogen.c

 Note: Scattering opacity is added to total opacity after all
       contributions have been computed.

 Note: The quantities chi_ai and eta_ai store the angle-inpendent
       opacities and emissivities in case atmos.moving == TRUE.
       In static atmospheres these quantities are just mapped to
       atmos.chi_c and atmos.eta_c to save memory space.

 Note: If write_analyze_output == FALSE the auxiliary output files for
       the Background Record Structure (BRS), metals, and molecules
       are NOT written. This option is used when Background is called
       from programs like solveray (formal solution along one specific
       ray) in cases with moving atmospheres (angle-dependent opacity).

 Note: If equilibria_only is set to TRUE only the electron density,
       LTE populations and collisions, and chemical equilibria are
       evaluated.

 Note: Record numbers stored in atmos.backgrrecno refer to records
       of the size atmos.Nspace. If a polarized line is present 9
       (4 + 4 + 1, no magneto-optical effects), or 12 (7 + 4 +1, with
       magneto-optical effects) records are used, otherwise 3 (1 + 1 + 1)
       records.
       If the atmosphere is moving, or if a polarized line is present
       data is stored for each angle and wavelength, otherwise data
       is stored once for each wavelength only.
       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "geometry.h"
#include "rhf1d.h"
//#include "xdr.h"


#define COMMENT_CHAR  "#"
#define  FILE_EXT ".dat"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];
extern BackgroundData bgdat;
extern rhinfo io;
extern rhbgmem *bmem;
extern MPI_t mpi;

/* --- Routines to keep the background opacities in memory 
   Author: Jaime de la Cruz Rodriguez (ISP-SU 2015)
   ---*/
void allocateBack(int nspect){
  int nstokes = 1, recsize = atmos.Nrays * 2;  
    
	  
  /* --- Stokes Wavelength --- */
  if(atmos.backgrflags[nspect].ispolarized) nstokes = 4;
    
  
    
  /* --- Background opacity --- */
  bmem[nspect].chi_b = matrix_double(recsize, atmos.Nspace*nstokes);
  
  /* --- Background emissivity --- */
  bmem[nspect].eta_b = matrix_double(recsize, atmos.Nspace*nstokes);
  
  /* --- Background scattering opacity --- */
  bmem[nspect].sca_b = matrix_double(recsize, atmos.Nspace);
  
  /* --- off-diagonal elements of propagation matrix K ? --- */
  
  if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical)
    bmem[nspect].chip_b = matrix_double(recsize, 3 * atmos.Nspace);
    
    
  bmem[nspect].allocated = true;
}



int writeBackground_j(int la, int mu, bool_t to_obs,
		      double *chi_c, double *eta_c, double *sca_c,
		      double *chip_c){

  const char routineName[] = "writeBackground_j";
  long recnum =  2*mu + to_obs, reclen = atmos.Nspace*sizeof(double);
  int nstokes = 1;

  if(!bmem[la].allocated)
    allocateBack(la);

  
  if(atmos.backgrflags[la].ispolarized)
    nstokes = 4;

  
  /* --- Copy data to allocated arrays --- */
  
  memcpy(&bmem[la].chi_b[recnum][0], &chi_c[0], nstokes*reclen);
  memcpy(&bmem[la].eta_b[recnum][0], &eta_c[0], nstokes*reclen);
  memcpy(&bmem[la].sca_b[recnum][0], &sca_c[0], reclen);


  if (atmos.backgrflags[la].ispolarized && input.magneto_optical && chip_c != NULL)
    memcpy(&bmem[la].chip_b[recnum][0], chip_c, reclen*3);


}
int readBackground_j(int la, int mu, bool_t to_obs){

  const char routineName[] = "readBackground_j";
  long recnum, reclen;
  int nstokes = 1;
  ActiveSet *as = &spectrum.as[la];


  reclen = atmos.Nspace*sizeof(double);
  recnum = 2*mu + to_obs;
  
  if(atmos.backgrflags[la].ispolarized && input.StokesMode == FULL_STOKES) nstokes = 4;

  /* --- Copy data to allocated arrays --- */
  
  memcpy(as->chi_c, &bmem[la].chi_b[recnum][0], nstokes*reclen);
  memcpy(as->eta_c, &bmem[la].eta_b[recnum][0], nstokes*reclen);
  memcpy(as->sca_c, &bmem[la].sca_b[recnum][0], reclen);
  
  if (atmos.backgrflags[la].ispolarized &&
      input.magneto_optical && as->chip_c != NULL &&
      input.StokesMode == FULL_STOKES){
    memcpy(&as->chip_c[0], &bmem[la].chip_b[recnum][0], reclen*3);
  }
    
}


void init_Background_j(){
   const char routineName[] = "init_Background";
  char    inputLine[MAX_LINE_SIZE];
  bool_t  exit_on_EOF;
  int     n, nspect;
  FILE   *fp_fudge;
  char    file_background[MAX_MESSAGE_LENGTH], *fext = FILE_EXT;
  long nstokes, recsize;

  
  if (strcmp(input.fudgeData, "none")) {
    bgdat.do_fudge = TRUE;
    
    /* --- Read wavelength-dependent fudge factors to compensate for
       missing UV backround line haze --           -------------- */
    
    if ((fp_fudge = fopen(input.fudgeData, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", input.fudgeData);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    sprintf(messageStr,
	    "\n-Fudging background opacities with file\n  %s\n\n",
	    input.fudgeData);
    Error(MESSAGE, routineName, messageStr);
    
    getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    sscanf(inputLine, "%d", &bgdat.Nfudge);
    bgdat.lambda_fudge = (double *) malloc(bgdat.Nfudge * sizeof(double));
    bgdat.fudge = matrix_double(3, bgdat.Nfudge);
    for (n = 0;  n < bgdat.Nfudge;  n++) {
      getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      sscanf(inputLine, "%lf %lf %lf %lf", &bgdat.lambda_fudge[n],
	     &bgdat.fudge[0][n], &bgdat.fudge[1][n], &bgdat.fudge[2][n]);
    }
    for (n = 0;  n < 3*bgdat.Nfudge;  n++) bgdat.fudge[0][n] += 1.0;
    fclose(fp_fudge);
  } else
    bgdat.do_fudge = FALSE;


  
  /* --- Allocate memory for the boolean array that stores whether
         a wavelength overlaps with a Bound-Bound transition in the
         background, or whether it is polarized --     -------------- */
  // printf("%d\n", spectrum.Nspect);
  atmos.backgrflags = (flags *) malloc(spectrum.Nspect * sizeof(flags));
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    atmos.backgrflags[nspect].hasline = FALSE;
    atmos.backgrflags[nspect].ispolarized = FALSE;
  }
  /* --- Allocate memory for the list of record numbers that specifies
         for each wavelength where to find the background opacity,
         scattering opacity, and emissivity --         -------------- */

  if (atmos.moving || atmos.Stokes) {
    atmos.backgrrecno = 
      (int *) malloc(2*spectrum.Nspect*atmos.Nrays * sizeof(int)); 
  } else
    atmos.backgrrecno = (int *) malloc(spectrum.Nspect * sizeof(int));


   /* --- Read background files from Kurucz data file -- ------------- */
  atmos.Nrlk = 0;
  readKuruczLines(input.KuruczData);
  if (atmos.Nrlk > 0) {
    qsort(atmos.rlk_lines, atmos.Nrlk, sizeof(RLK_Line), rlk_ascend);
  }
  

  
  /* --- Init to null arrays to store data --- */
  bmem = (rhbgmem*) malloc(spectrum.Nspect * sizeof(rhbgmem));
  memset(bmem, 0, spectrum.Nspect * sizeof(rhbgmem));
  
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) bmem[nspect].allocated = false;


  return;
  
}


void Background_j(bool_t write_analyze_output, bool_t equilibria_only)
{
  const char routineName[] = "Background_p";
  register int k, nspect, n, mu, to_obs;
  
  static int ne_iter = 0;
  bool_t  do_fudge, fromscratch;
  int     index, Nfudge, NrecStokes;
  double *chi, *eta, *scatt, wavelength, *thomson, *chi_ai, *eta_ai, *sca_ai,
    Hmin_fudge, scatt_fudge, metal_fudge, *lambda_fudge, **fudge,
    *Bnu, *chi_c, *eta_c, *sca_c, *chip, *chip_c;
  Atom   *He;
  Element *element;
  flags   backgrflags;
  char    file_background[MAX_MESSAGE_LENGTH], *fext = FILE_EXT;


  getCPU(2, TIME_START, NULL);

  if (input.solve_ne == ONCE  || input.solve_ne == ITERATION ) {
    fromscratch = (input.solve_ne == ONCE  ||
		   (input.solve_ne == ITERATION  &&  ne_iter == 0)) ?
      TRUE : FALSE;
    Solve_ne(atmos.ne, fromscratch);
    ne_iter++;
  }
  SetLTEQuantities();

  if (input.NonICE)
    readMolecules(MOLECULAR_CONCENTRATION_FILE);
  else{
    ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);
    if(mpi.stop){
      return;
    }
  }
  if (equilibria_only) {

    /* --- If we only need ne, LTE populations and collisions, and
           chemical equilibrium leave here --          -------------- */

    getCPU(2, TIME_POLL, "Total Background");
    return;
  }

  getCPU(3, TIME_START, NULL);

  /* Get fudge data */
  lambda_fudge = bgdat.lambda_fudge;
  do_fudge     = bgdat.do_fudge;
  Nfudge       = bgdat.Nfudge;
  fudge        = bgdat.fudge;


  
  
  

  /* --- Allocate temporary storage space. The quantities are used
         for the following purposes:

       - chi, eta, scatt: Get contributions to opacity, emissivity,
         and scattering opacity, respectively, from a specific process
         for a given wavelength and possibly angle.

       - chi_c, eta_c, sca_c: Collect the total opacity, emissivity
         and scattering opacity for a given wavelength and possibly
         angle.

       - chi_ai, eta_ai: Collect the angle-independent part of
         opacity and emissivity for each wavelength so that these
         need not be recalculated in an angle-dependent case.
         When the atmosphere is not moving and has no magnetic fields
         these just point to the total quantities chi_c and eta_c.

   Note: In case of magnetic fields in the atmosphere chi, eta and 
         chip, and chi_c, eta_c and chip_c contain all four Stokes
         parameters, and should be allocated a 4 and 3 times larger
         storage space, respectively.
         --                                            -------------- */

  if (atmos.Stokes)
    NrecStokes = 4;
  else
    NrecStokes = 1;

  chi_c = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));
  eta_c = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));
  sca_c = (double *) malloc(atmos.Nspace * sizeof(double));

  chi   = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));
  eta   = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));
  scatt = (double *) malloc(atmos.Nspace * sizeof(double));
    
  if (atmos.Stokes && input.magneto_optical) {
    chip   = (double *) malloc(3*atmos.Nspace * sizeof(double));
    chip_c = (double *) malloc(3*atmos.Nspace * sizeof(double));
  } else {
    chip   = NULL;
    chip_c = NULL;
  }

  if (atmos.moving || atmos.Stokes) {
    chi_ai = (double *) malloc(atmos.Nspace * sizeof(double));
    eta_ai = (double *) malloc(atmos.Nspace * sizeof(double));
    sca_ai = (double *) malloc(atmos.Nspace * sizeof(double));
  } else {
    chi_ai = chi_c;
    eta_ai = eta_c;
    sca_ai = sca_c;
  }
  Bnu = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Thomson scattering by free electrons is wavelength independent
         in non-relativistic limit so we compute it only once -- ---- */

  thomson = (double *) malloc(atmos.Nspace * sizeof(double));
  Thomson(thomson);

  /* --- Check whether an atomic model is present for He -- --------- */

  He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;


  /* --- Go through the spectrum and add the different opacity and
         emissivity contributions. This is the main loop --  -------- */

  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    wavelength = spectrum.lambda[nspect];

    /* --- The Planck function at this wavelength --   -------------- */

    Planck(atmos.Nspace, atmos.T, wavelength, Bnu);

    /* --- Initialize the flags for this wavelength -- -------------- */

    atmos.backgrflags[nspect].hasline     = FALSE;
    atmos.backgrflags[nspect].ispolarized = FALSE;

    /* --- Initialize angle-independent quantities --  -------------- */

    for (k = 0;  k < atmos.Nspace;  k++) {
      chi_ai[k] = 0.0;
      eta_ai[k] = 0.0;
      sca_ai[k] = thomson[k];
    }
    /* --- Negative hydrogen ion, bound-free and free-free -- ------- */

    if (Hminus_bf(wavelength, chi, eta)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k];
	eta_ai[k] += eta[k];
      }
    }
    if (Hminus_ff(wavelength, chi)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k];
	eta_ai[k] += chi[k] * Bnu[k];
      }
    }
    /* --- Opacity fudge factors, applied to Hminus opacity -- ------ */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[0],
	     1, &wavelength, &Hmin_fudge, FALSE);
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] *= Hmin_fudge;
	eta_ai[k] *= Hmin_fudge;
      }
    }
    /* --- Opacities from bound-free transitions in OH and CH -- ---- */

    if (OH_bf_opac(wavelength, chi, eta)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k];
	eta_ai[k] += eta[k];
      }
    }
    if (CH_bf_opac(wavelength, chi, eta)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k];
	eta_ai[k] += eta[k];
      }
    }
    /* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

    if (Hydrogen_bf(wavelength, chi, eta)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k];
	eta_ai[k] += eta[k];
      }
    }
    Hydrogen_ff(wavelength, chi);
    for (k = 0;  k < atmos.Nspace;  k++) {
      chi_ai[k] += chi[k]; 
      eta_ai[k] += chi[k] * Bnu[k];
    }
    /* --- Rayleigh scattering by neutral hydrogen --  -------------- */

    if (Rayleigh(wavelength, atmos.H, scatt)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	sca_ai[k]  += scatt[k];
      }
    }
    /* --- Rayleigh scattering by neutral helium --    -------------- */

    if (He && Rayleigh(wavelength, He, scatt)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	sca_ai[k]  += scatt[k];
      }
    }
    /* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

    if (H2plus_ff(wavelength, chi)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k]; 
	eta_ai[k] += chi[k] * Bnu[k];
      }
    }
    /* --- Rayleigh scattering and free-free absorption by
           molecular hydrogen --                       -------------- */

    if (Rayleigh_H2(wavelength, scatt)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	sca_ai[k]  += scatt[k];
      }
    }
    if (H2minus_ff(wavelength, chi)) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi_ai[k] += chi[k]; 
	eta_ai[k] += chi[k] * Bnu[k];
      }
    }
    /* --- Bound-Free opacities due to ``metals'' --   -------------- */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[2],
	     1, &wavelength, &metal_fudge, FALSE);
    } else {
      metal_fudge = 1.0;
    }
    /* --- Note: Hydrogen bound-free opacities are calculated in
           routine Hydrogen_bf --                      -------------- */

    Metal_bf(wavelength, atmos.Natom-1, atmos.atoms+1, chi, eta);
    for (k = 0;  k < atmos.Nspace;  k++) {
      chi_ai[k] += chi[k] * metal_fudge;
      eta_ai[k] += eta[k] * metal_fudge;
    }
    /* --- Add the scattering opacity to the absorption part to store
           the total opacity --                        -------------- */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[1],
	     1, &wavelength, &scatt_fudge, FALSE);
    } else {
      scatt_fudge = 1.0;
    }
    for (k = 0;  k < atmos.Nspace;  k++) {
      sca_ai[k] *= scatt_fudge;
      chi_ai[k] += sca_ai[k];
    }
    /* --- Now the contributions that may be angle-dependent due to the
           presence of atomic or molecular lines --    -------------- */

    if (atmos.moving || atmos.Stokes) {
      for (mu = 0;  mu < atmos.Nrays;  mu++) {
        for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  index = 2*(nspect*atmos.Nrays + mu) + to_obs;

          /* --- First, copy the angle-independent parts -- --------- */

	  for (k = 0;  k < atmos.Nspace;  k++) {
	    chi_c[k] = chi_ai[k];
	    eta_c[k] = eta_ai[k];
            sca_c[k] = sca_ai[k];
	  }
	  
          /* --- Zero the polarized quantities, if necessary -- ----- */

	  if (atmos.Stokes) {
           for (k = atmos.Nspace;  k < 4*atmos.Nspace;  k++) {
              chi_c[k] = 0.0;
              eta_c[k] = 0.0;
            }
            if (input.magneto_optical)
              for (k = 0;  k < 3*atmos.Nspace;  k++) chip_c[k] = 0.0;
	  }
          /* --- Add opacity from passive atomic lines (including
                 hydrogen) --                          -------------- */

	  if (input.allow_passive_bb) {
	    backgrflags = passive_bb(wavelength, nspect, mu, to_obs,
				     chi, eta, chip);
	    if (backgrflags.hasline) {
	      atmos.backgrflags[nspect].hasline = TRUE;
	      if (backgrflags.ispolarized) {
                NrecStokes = 4;
                atmos.backgrflags[nspect].ispolarized = TRUE;
		if (input.magneto_optical) {
                  for (k = 0;  k < 3*atmos.Nspace;  k++)
                    chip_c[k] += chip[k];
                }
              } else
                NrecStokes = 1;

              for (k = 0;  k < NrecStokes*atmos.Nspace;  k++) {
                chi_c[k] += chi[k];
                eta_c[k] += eta[k];
              }

	    }
	  }
          /* --- Add opacity from Kurucz line list --  -------------- */

          if (atmos.Nrlk > 0) {
	    backgrflags = rlk_opacity(wavelength, nspect, mu, to_obs,
				      chi, eta, scatt, chip);
	    if (backgrflags.hasline) {
	      atmos.backgrflags[nspect].hasline = TRUE;
              if (backgrflags.ispolarized) {
                NrecStokes = 4;
                atmos.backgrflags[nspect].ispolarized = TRUE;
                if (input.magneto_optical) {
                  for (k = 0;  k < 3*atmos.Nspace;  k++)
                    chip_c[k] += chip[k];
                }
              } else
                NrecStokes = 1;

              for (k = 0;  k < NrecStokes*atmos.Nspace;  k++) {
                chi_c[k] += chi[k];
                eta_c[k] += eta[k];
              }
	      if (input.rlkscatter) {
		for (k = 0;  k < atmos.Nspace;  k++) {
		  sca_c[k] += scatt[k];
		  chi_c[k] += scatt[k];
		}
	      }
	    }
	  }
	  /* --- Add opacity from molecular lines --   -------------- */

	  backgrflags = MolecularOpacity(wavelength, nspect, mu, to_obs,
					 chi, eta, chip);
	  if (backgrflags.hasline) {
	    atmos.backgrflags[nspect].hasline = TRUE;
            if (backgrflags.ispolarized) {
              NrecStokes = 4;
              atmos.backgrflags[nspect].ispolarized = TRUE;
              if (input.magneto_optical) {
                for (k = 0;  k < 3*atmos.Nspace;  k++)
                  chip_c[k] += chip[k];
              }
            } else
              NrecStokes = 1;

            for (k = 0;  k < NrecStokes*atmos.Nspace;  k++) {
              chi_c[k] += chi[k];
              eta_c[k] += eta[k];
            }
	  }
	  /* --- Store angle-dependent results only if at least one
                 line was found at this wavelength --  -------------- */
	  
	  //	  atmos.backgrrecno[index] = backgrrecno;
	  if ((mu == atmos.Nrays-1 && to_obs) ||
	      (atmos.backgrflags[nspect].hasline && 
	       (atmos.moving || atmos.backgrflags[nspect].ispolarized))) {



	    if ((mu == atmos.Nrays-1 && to_obs) && !(atmos.backgrflags[nspect].hasline && 
						     (atmos.moving || atmos.backgrflags[nspect].ispolarized)) ){
	      writeBackground_j(nspect, 0, 0,
				chi_c, eta_c, sca_c, chip_c);
	    } else {
	      
	      
	      writeBackground_j(nspect, mu, to_obs,
				chi_c, eta_c, sca_c, chip_c);
	    }
	  }
	  
	}
      }    
    } else {
      /* --- Angle-independent case. First, add opacity from passive
	     atomic lines (including hydrogen) --      -------------- */

      if (input.allow_passive_bb) {
	backgrflags = passive_bb(wavelength, nspect, 0, TRUE,
				 chi, eta, NULL);
	if (backgrflags.hasline) {
	  atmos.backgrflags[nspect].hasline = TRUE;
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    chi_c[k] += chi[k];
	    eta_c[k] += eta[k];
	  }
	}
      }
      /* --- Add opacity from Kurucz line list --      -------------- */

      if (atmos.Nrlk > 0) {
	backgrflags = rlk_opacity(wavelength, nspect, 0, TRUE,
				  chi, eta, scatt, NULL);
	if (backgrflags.hasline) {
	  atmos.backgrflags[nspect].hasline = TRUE;
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    chi_c[k] += chi[k];
	    eta_c[k] += eta[k];
	  }
	  if (input.rlkscatter) {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      sca_c[k] += scatt[k];
	      chi_c[k] += scatt[k];
	    }
	  }
	}
      }
      /* --- Add opacity from molecular lines --       -------------- */

      backgrflags = MolecularOpacity(wavelength, nspect, 0, TRUE,
				     chi, eta, NULL);
      if (backgrflags.hasline) {
	atmos.backgrflags[nspect].hasline = TRUE;
	for (k = 0;  k < atmos.Nspace;  k++) {
	  chi_c[k] += chi[k];
	  eta_c[k] += eta[k];
	}
      }
      /* --- Store results --                          -------------- */

      //  atmos.backgrrecno[nspect] = backgrrecno;
      writeBackground_j(nspect, 0, 0,
		      chi_c, eta_c, sca_c, NULL);
    }
  }

  
 
  getCPU(3, TIME_POLL, "Background Opacity");

  /* --- Free the temporary space allocated in the ff routines -- --- */

  Hminus_ff(0.0, NULL);
  H2minus_ff(0.0, NULL);
  H2plus_ff(0.0, NULL);

  free(chi);    free(eta);  free(scatt);  free(Bnu);  free(thomson);
  free(chi_c);  free(eta_c);  free(sca_c);

  if (atmos.moving || atmos.Stokes) {
    free(chi_ai);
    free(eta_ai);
    free(sca_ai);
  }
  if (atmos.Stokes && input.magneto_optical) {
    free(chip);
    free(chip_c);
  }


  /* --- Free the element populations for species used in rlk_opacity --- */
  
  for (k = 0;  k < atmos.Nelem;  k++) {
    element = &atmos.elements[k];
    if (element->n != NULL) {
      freeMatrix((void **) element->n);
      element->n = NULL;
    }
    
  }
  
  getCPU(2, TIME_POLL, "Total Background");
}
/* ------- end ---------------------------- Background.c ------------ */
