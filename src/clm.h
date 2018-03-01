#ifndef CLM_H
#define CLM_H
/* 
   Implementation of the clm class (Levenberg-Marquardt least-squares-fit)
   author: J. de la Cruz Rodriguez (Department of Astronomy, Stockholm University 2016)
   
   DEPENDENCIES:
          EIGEN3, for SVD decomposition.

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
#include <eigen3/Eigen/Dense>


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

/* --- Struct to store the regularization terms --- */

struct reg_t{
  bool to_reg;
  int npar, nreg;
  double scl;
  double **dreg;
  double *reg;

reg_t():to_reg(false), npar(0), nreg(0), scl(1.0), dreg(NULL), reg(NULL){};
  reg_t(int npar_in, int nreg_in, double scl_in);
  reg_t(const reg_t &in);
  ~reg_t();
  
  void set(int npar_in, int nreg_in, double scl_in);
  void zero();
  double getReg();
  void del();
  void copyReg(double *reg_in);
  reg_t &operator=(const reg_t &in);
  
  
};



/* --- 
   Type for the function that will compute chi2 and the response function 
   --- */

typedef int (*clm_func)(int npar, int nd, double *x, double *syn_in, double *res,
			double **rf, void *mydata, reg_t &regul, bool store);

double sumarr(double *arr, size_t n);
double sumarr2(double *arr, int n);


/* --- Class definitions --- */

class clm{
 private:
  int nd, npar, nzero, lwork;
  bool error;
 public:
  std::vector<clmf> fcnt;
  std::vector<double> diag, tmp, bestSyn, iSyn;
  std::vector<unsigned> ptype, ntype;
  std::vector<std::vector<unsigned>> pidx;
  bool verb, regularize, first;
  double xtol, chi2_thres, svd_thres, lfac, lmax, lmin, ilambda, regul_scal, regul_scal_in, reset_par;
  int maxreject, proc, nvar, use_geo_accel, delay_bracket;

  
  /* --- Constructor / Destructor --- */
  
  clm(int ind, int inpar);
  ~clm();


  
  /* --- prototypes --- */

  void checkParameters(double *x);
  void checkMaxChange(double *dx, double *x);
  void normalizeParameters(double *x);
  void scaleParameters(double *x);
  double fitdata(clm_func ifx, double *x, void *mydat, int maxiter, reg_t &reg);
  double compute_chi2(double *res, double penalty);
  double getChi2Pars(double *res, double **rf, double lambda,
		     double *x, double *xnew, void *mydat, clm_func fx, reg_t &regul);
  double getChi2ParsLineSearch(double *res, double **rf, double &lambda,
			       double *x, double *xnew, void *mydat,
			       clm_func fx, reg_t &regul, double rchi2, bool braket = true);
  //  void compute_trial3(double *res, double **rf, double lambda,
  //	      double *x, double *xnew);
  void compute_trial3(double *res, double **rf, double lambda,
		      double *x, double *xnew, reg_t &regul, void *mydat, clm_func fx);
  //void backsub(double **u, double *w, double **v, int n,
  //	       double *b, double *x);
  void scaleRF(double **rf);
  double checkLambda(double lamb);
  void zero(double *res, double **rf);
  void getParTypes();
  void backSub(int n, Eigen::MatrixXd &u, Eigen::VectorXd &w, Eigen::MatrixXd &v,  double *b);
  void geoAcceleration(double *x, double *dx, double h, Eigen::MatrixXd &A, Eigen::MatrixXd &LL,
		       double *res, void *mydat, reg_t &dregul, double **rf, clm_func fx, double lam);


  double getTime(double t0 = -1.0){
    struct timeval dum;
    gettimeofday(&dum, NULL);
    if(t0 < 0.0) return dum.tv_sec + dum.tv_usec * 1.0E-6;
    else return (dum.tv_sec + dum.tv_usec * 1.0E-6) - t0;
  }
};

#endif
