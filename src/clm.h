#ifndef CLM_H
#define CLM_H
/* 
   Implementation of the clm class (Levenberg-Marquardt least-squares-fit)
   author: J. de la Cruz Rodriguez (Department of Astronomy, Stockholm University 2016)
   
   DEPENDENCIES:
          EIGEN3, LaPack (routines commented out, but still working).

   PARAMETERS: (DEFAULT VALUES) 
      svd_thres  = 1.e-16 -> Relative threshold for small singular values
      ilambda    = 1.0    -> Initial damping parameter for the Hessian diag.
      lfac       = 10     -> Scale factor for lambda.
      lmax       = 1.e+4  -> Max lambda parameter
      lmin       = 1.e-4  -> Min lambda parameter
      chi2_thres = 1.0    -> Inversion will be finished if Chi2 <= chi2_thres
      xtol       = 1.e-5  -> Min relative change in Chi2
      maxreject  = 6      -> Allowed number of consecutive 
                             changes to lambda when iterations fail.

    USER FUNCT: The user must provide a function that will compute RF and the residues.
                It will take the following pars:

	    void funct(int npar, int nd, double *x, double *res, double **rf, void *mydata);

	    where:
	        npar = number of free parameters
		nd   = number of data points
		x    = current estimate of parameters
		res  = [nd], array to store the residues
		rf   = [npar][nd], if != NULL, storage of the Jacobian for each parameter
		mydata = A pointer to data that the user needs to use uin funct 
		        (it can be a struct, object, variable...whatever, just cast 
			it into the right type inside funct).


	    
		   
*/
#include <vector>
#include <sys/time.h>


/* --- 
   Struct to control each of the parameters 
   --- */

struct clmf{
  double limit[2];
  bool capped;
  double maxchange[2];
  double scl;
  bool cyclic;
  bool bouncing;
  bool userscaled;
  bool relchange;
};




/* --- 
   Type for the function that will compute chi2 and the response function 
   --- */

typedef int (*clm_func)(int npar, int nd, double *x, double *res, double **rf, void *mydata);
double sumarr(double *arr, size_t n);
double sumarr2(double *arr, int n);


/* --- Class definitions --- */

class clm{
 private:
  int nd, npar, nzero, lwork;
  bool error;
 public:
  std::vector<clmf> fcnt;
  std::vector<double> diag, tmp;
  bool verb;
  double xtol, chi2_thres, svd_thres, lfac, lmax, lmin, ilambda;
  int maxreject, proc;

  
  /* --- Constructor / Destructor --- */
  
  clm(int ind, int inpar);
  ~clm();


  
  /* --- prototypes --- */

  void checkParameters(double *x);
  void checkMaxChange(double *dx, double *x);
  void normalizeParameters(double *x);
  void scaleParameters(double *x);
  double fitdata(clm_func ifx, double *x, void *mydat, int maxiter = 50);
  double compute_chi2(double *res);
  double getChi2Pars(double *res, double **rf, double lambda,
		     double *x, double *xnew, void *mydat, clm_func fx);
  //  void compute_trial3(double *res, double **rf, double lambda,
  //	      double *x, double *xnew);
  void compute_trial2(double *res, double **rf, double lambda,
		      double *x, double *xnew);
  void backsub(double **u, double *w, double **v, int n,
	       double *b, double *x);
  void scaleRF(double **rf);
  double checkLambda(double lamb);
  void zero(double *res, double **rf);

  double getTime(double t0 = -1.0){
    struct timeval dum;
    gettimeofday(&dum, NULL);
    if(t0 < 0.0) return dum.tv_sec + dum.tv_usec * 1.0E-6;
    else return (dum.tv_sec + dum.tv_usec * 1.0E-6) - t0;
  }
};

#endif
