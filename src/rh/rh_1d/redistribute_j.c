/* ------- file: -------------------------- redistribute.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 14:01:22 2009 --

       --------------------------                      ----------RH-- */

/* --- Administers iterations to redistribute intensity in PRD line while
       keeping the population number fixed. --         -------------- */
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "error.h"
#include "inputs.h"
#include "rhf1d.h"

/* --- Function prototypes --                          -------------- */
double MaxChange_j(struct Ng *Ngs, char *text, bool_t quiet);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
//extern MPI_data mpi;
extern char messageStr[];
extern MPI_t mpi;

/* ------- begin -------------------------- Redistribute.c ---------- */

void Redistribute(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Redistribute";
  register int kr, nact;

  bool_t  quiet, accel, eval_operator, redistribute;
  enum    Interpolation representation;
  int     niter, Nlamu;
  double  drho, drhomax, drhomaxa;
  Atom *atom;
  AtomicLine *line;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    
    /* --- Initialize structures for Ng acceleration PRD iteration -- */

    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];
      if (line->PRD && line->Ng_prd == NULL) {
	if (input.PRD_angle_dep == PRD_ANGLE_DEP)
	  Nlamu = 2*atmos.Nrays * line->Nlambda * atmos.Nspace;
	else
	  Nlamu = line->Nlambda*atmos.Nspace;
	
	line->Ng_prd = NgInit(Nlamu, input.PRD_Ngdelay,
			      input.PRD_Ngorder, input.PRD_Ngperiod,
			      line->rho_prd[0]);
      }
    }
  }
  /* --- Iterate over scattering integral while keeping populations
         fixed --                                      -------------- */

  niter = 1;
  while (niter <= NmaxIter) {

    drhomaxa = 0.0;
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      drhomax = 0.0;
      for (kr = 0;  kr < atom->Nline;  kr++) {
	line = &atom->line[kr];
	if (line->PRD) {
	   switch (input.PRD_angle_dep) {
	  case PRD_ANGLE_INDEP:
	    PRDScatter(line, representation=LINEAR);
	    break;

	  case PRD_ANGLE_APPROX:
	    PRDAngleApproxScatter(line, representation=LINEAR); 
	    break;
	    
	  case PRD_ANGLE_DEP:
	    PRDAngleScatter(line, representation=LINEAR);
	    break;
	  }

	  accel = Accelerate(line->Ng_prd, line->rho_prd[0]);
	  sprintf(messageStr, "  PRD: iter #%d, atom %s, line %d,",
		  line->Ng_prd->count-1, atom->ID, kr);
	  drho = MaxChange_j(line->Ng_prd, messageStr, quiet=FALSE);
	  if (mpi.stop) return; /* Get out if singular matrix, or NAN in dmax */
	  sprintf(messageStr, (accel) ? " (accelerated)\n" : "\n");
	  Error(MESSAGE, routineName, messageStr);

	  drhomax = MAX(drho, drhomax);
	}
	drhomaxa = MAX(drhomax, drhomaxa);
      }
    }

    /* --- Solve transfer equation with fixed populations -- -------- */

    solveSpectrum(eval_operator=FALSE, redistribute=TRUE);

    if (drhomaxa < iterLimit) break;
    niter++;
  }
}
/* ------- end ---------------------------- Redistribute.c ---------- */

/* ------- begin -------------------------- MaxChange_j.c ------------- */

double MaxChange_j(struct Ng *Ngs, char *text, bool_t quiet)
{
  register int k;

  double dmax = 0.0, *old, *new;

  if (Ngs->count < 2) return dmax;

  old = Ngs->previous[(Ngs->count - 2) % (Ngs->Norder + 2)];
  new = Ngs->previous[(Ngs->count - 1) % (Ngs->Norder + 2)];
  for (k = 0;  k < Ngs->N;  k++) {
    if (new[k])
      dmax = MAX(dmax, fabs((new[k] - old[k]) / new[k]));
  }
  if (!quiet)
    fprintf(commandline.logfile, "%s delta = %6.4E", text, dmax);

  if (isnan(dmax) || isinf(dmax)) mpi.stop = TRUE;

  return dmax;
}
/* ------- end ---------------------------- MaxChange_j.c ------------- */
