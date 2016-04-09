#ifndef CLM_H
#define CLM_H
/* 
   Implementation of the clm class (Levenberg-Marquardt least-squares-fit)
   author: J. de la Cruz Rodriguez (Department of Astronomy, Stockholm University 2016)
   
   DEPENDENCIES:
         There are two implementations of the routine that gets the corrections to the model. one based on GSL and the other one based on LAPACK (DEFAULT). A third implementation is commented out in the implementation file, which makes use of EIGEN3.

   PARAMETERS: (DEFAULT VALUES) 
      svd_thres  = 1.e-16 -> Relative threshold for small singular values
      ilambda    = 1.e-2  -> Initial damping parameter for the Hessian diag.
      lfac       = 5.0    -> Scale factor for lambda.
      lmax       = 1.e+5  -> Max lambda parameter
      lmin       = 1.e-10 -> Min lambda parameter
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


/* --- 
   Struct to control each of the parameters 
   --- */

struct clmf{
  double limit[2];
  bool capped;
  double maxchange;
  double scl;
  bool cyclic;
  bool bouncing;
  bool userscaled;
};




/* --- 
   Type for the function that will compute chi2 and the response function 
   --- */

typedef int (*clm_func)(int npar, int nd, double *x, double *res, double **rf, void *mydata);
double sumarr(double *arr, int n);
double sumarr2(double *arr, int n);


/* --- Class definitions --- */

class clm{
 private:
  int nd, npar, nzero, lwork;
 public:
  std::vector<clmf> fcnt;
  std::vector<double> diag;
  bool verb;
  double xtol, chi2_thres, svd_thres, lfac, lmax, lmin, ilambda;
  int maxreject;

  
  /* --- Constructor / Destructor --- */
  
  clm(int ind, int inpar);
  ~clm();


  
  /* --- prototypes --- */

  void checkParameters(double *x);
  void checkMaxChange(double *dx);
  void normalizeParameters(double *x);
  void scaleParameters(double *x);
  double fitdata(clm_func ifx, double *x, void *mydat, int maxiter = 50);
  double compute_chi2(double *res);
  double getChi2Pars(double *res, double **rf, double lambda,
		     double *x, double *xnew, void *mydat, clm_func fx);
  void compute_trial3(double *res, double **rf, double lambda,
		      double *x, double *xnew);
  void backsub(double **u, double *w, double **v, int n,
	       double *b, double *x);
  void scaleRF(double **rf);
  double checkLambda(double lamb);
  void zero(double *res, double **rf);
};

#endif
