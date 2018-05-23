/* ----------------------------------------------------------------------- 

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   ----------------------------------------------------------------------- */

#ifndef BEZIER_H
#define BEZIER_H


/* ----- Prototypes auxiliary functions --- */
double sign(const double val);
void SIMD_MatInv(float* src); // Matrix inversion Cramer method SSE instructions
void m4inv(double MI[4][4]);  // Matrix inversion Shipley-Coleman (not used)

void m4v(double a[4][4], double b[4], double c[4]);
void m4m(double a[4][4], double b[4][4], double c[4][4]);

void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
		    double chiup[4], double chic[4], double chidn[4]);

void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
		    double chiup[4][4], double chic[4][4], double chidn[4][4]);

double cent_deriv(double dsup,double dsdn, double chiup,double chic, double chidn);
void Svec(int k, double **S, double *Sf);

void Bezier3_coeffs(double dt, double *alpha, double *beta,
		    double *gamma, double *theta, double *eps);


/* --- Prototypes for cubic Bezier solvers --- */

void PiecewiseStokesBezier3(int nspect, int mu, bool_t to_obs,
			    double *chi, double **S, double **I, double *Psi);

void Piecewise_Bezier3(int nspect, int mu, bool_t to_obs,
		       double *chi, double *S, double *I, double *Psi);

/* -------------------------------------------------------------------------- */


#endif
