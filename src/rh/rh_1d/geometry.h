/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Feb 24 08:21:55 2012 --

       --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

/* --- Define geometric quantities for 1-D plane-parallel version --  */


enum boundcond  {ZERO, THERMALIZED, IRRADIATED, REFLECTIVE};
enum mass_scale {GEOMETRIC, COLUMN_MASS, TAU500};
enum vertical   {TOP=0, BOTTOM};

typedef struct {
  enum     mass_scale  scale;
  enum     boundcond vboundary[2]; 
  int      Ndep, Nrays;
  double  *height, *cmass, *tau_ref, *mux, *muy, *muz, *wmu, *vel,
         **Itop, **Ibottom;
} Geometry;

typedef struct {
  char     *brs_fname;
  bool_t    do_fudge, write_BRS;
  int       Nfudge;
  int       brs_ncid,    brs_hl_var,    brs_ip_var,  brs_nrec_var;
  // for now here, in the future perhaps write things in io.h:
  int       j_ncid, jlambda_var, j20_var; 
  double   *lambda_fudge, **fudge;
  double **chi_c,**eta_c,*sca_c,**chip_c;
} BackgroundData;



/* --- Associated function prototypes --               -------------- */

void convertScales(Atmosphere *atmos, Geometry *geometry);
void getAngleQuad(Geometry *geometry);
void getBoundary(Geometry *geometry);
void MULTIatmos(Atmosphere *atmos, Geometry *geometry);
void writeGeometry(Geometry *geometry);


/* --- Formal solution related --                      -------------- */

double Feautrier(int nspect, int mu, double *chi, double *S,
		 enum FeautrierOrder order, double *P, double *Psi);
void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi);
void PieceBezier_1D(int nspect, int mu, bool_t to_obs,
		    double *chi, double *S, double *I, double *Psi);

void PiecewiseStokes(int nspect, int mu, bool_t to_obs,
		     double *chi_I, double **S, double **I, double *Psi);
void Piecewise_Hermite_1D(int nspect, int mu, bool_t to_obs,
			  double *chi, double *S, double *I, double *Psi);
void Piecewise_lbrHermite_1D(int nspect, int mu, bool_t to_obs,
			     double *chi, double *S, double *I, double *Psi);
#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
