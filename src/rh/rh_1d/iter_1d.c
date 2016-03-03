/* ---------------------------------------- Iterate.c --------------- */
/* ------- One-D, plane-parallel version ---------------------------- */

/* Implements a simple single-line scattering problem with given
 * photon destruction probability epsilon and Planck function Bp.
 *
 * Han Uitenbroek
 * Last modified: Wed Jun 16 14:30:55 2010 --
 */
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern int Nlambda;
extern double *Bp, *epsilon, *phi, *wlamb, wphi, *chi, *Sny, **Iemerge;
extern char    messageStr[];
extern struct  Ng *NgS;

extern Geometry geometry;


/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate(int NmaxIter, double iterLimit)
{
  register int mu, la, l, nIter;

  int    Ndep = geometry.Ndep, Nrays = geometry.Nrays;
  bool_t P_only, quiet, accel;
  enum   FeautrierOrder F_order;
  double wqmu, *dJny, *chila, *diagonal, *P, *Psi;

  FILE *fp_iter;

  dJny     = (double *) malloc(Ndep * sizeof(double));
  chila    = (double *) malloc(Ndep * sizeof(double));
  diagonal = (double *) malloc(Ndep * sizeof(double));
  P        = (double *) malloc(Ndep * sizeof(double));
  Psi      = (double *) malloc(Ndep * sizeof(double));

  Iemerge  = matrix_double(Nrays, Nlambda);

  /* --- Write iterations  to file --                 --------------- */

  fp_iter = fopen("iter_1d_Snu", "w");

  /* --- Start main iteration loop --                 --------------- */

  nIter = 1;
  while (nIter <= NmaxIter) {
    getCPU(2, TIME_START, NULL);

    for (l = 0;  l < Ndep;  l++) {
      dJny[l] = diagonal[l] = 0.0;
    }
    for (mu = 0;  mu < Nrays;  mu++) {
      for (la = 0;  la < Nlambda;  la++) {
	for (l = 0;  l < Ndep;  l++) {
	  chila[l] = phi[la] * chi[l];
	}

        /* --- Formal solution and diagonal operator -- ------------- */

	Iemerge[mu][la] =
	  Feautrier(la, mu, chila, Sny, F_order=STANDARD, P, Psi);

        /* --- Angle and wavelength integration --      ------------- */

	wqmu = geometry.wmu[mu] * phi[la]*wlamb[la];
	for (l = 0;  l < Ndep;  l++) {
	  dJny[l]     += wqmu*P[l];
	  diagonal[l] += wqmu*Psi[l];
	}
      }
    }
    /* --- Accelerated lambda iteration --              ------------- */

    for (l = 0;  l < Ndep;  l++) {
      diagonal[l] *= wphi;
      dJny[l] = wphi*dJny[l] - diagonal[l] * Sny[l];
      Sny[l]  = ((1.0 - epsilon[l])*dJny[l] + epsilon[l]*Bp[l]) /
        (1.0 - (1.0 - epsilon[l])*diagonal[l]);
    }
    /* --- Ng acceleration if requested --              ------------- */

    accel = Accelerate(NgS, Sny);
    sprintf(messageStr, "-- Main iteration No: %d", nIter);
    if (MaxChange(NgS, messageStr, quiet=FALSE) <= iterLimit)  break;
    Error(MESSAGE, NULL, (accel == TRUE) ? " (accelerated)\n" : "\n");

    fwrite(Sny, sizeof(double), Ndep, fp_iter);

    getCPU(2, TIME_POLL, "Iterate");
    nIter++;
  }

  free(dJny);  free(chila);  free(diagonal);
  free(P);     free(Psi);
  NgFree(NgS);

  fclose(fp_iter);
}
/* ------- end ---------------------------- Iterate.c --------------- */
