/* 
   Implementation of the clm class (Levenberg-Marquardt least-squares-fit)
   author: J. de la Cruz Rodriguez (Stockholm University 2016)

   Documentation included in the header file clm.h.

   Dependencies: Eigen3 for SVD decomp.

   Modifications: 
           2016-04-08, JdlCR: Added Kahan summation to perform the matrix multiplication,
	                      it improves convergence significantly!
	   
	   2016-06-26, JdlCR: Changed to Eigen3 routines for SVD, 
	                      the implementation is much cleaner and faster than LaPack.

	   2016-12-13, JdlCR: added the possibility to add l-2 regulatization. The user 
	                      must provide a function that computes the regularization 
			      term for all parameter.

*/

#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include "clm.h"

using namespace std;
using namespace Eigen;

/* -------------------------------------------------------------------------------- */

/* --- Internal routines --- */


double **mat2d(int nx1, int nx2){
  double **p;
  p = new double* [nx1];
  p[0] = new double [nx1 * nx2]();
  for(int x1=1;x1<nx1;++x1) p[x1] = p[x1-1] + nx2;
  return p;
}

void del_mat(double **p){
  delete[] (p[0]);
  delete[] (p);
}

double **var2dim(double *data, int nx1, int nx2){
  double **p;
  p=new double* [nx1];
  p[0] = data;
  for(int x1=1;x1<nx1;++x1) p[x1]=p[x1-1] + nx2;
  return p;
}

/* --- Kahan summation for the elements of a vector--- */

inline double sumarr(double *arr, size_t n){

  long double sum = 0.0L, c = 0.0L;
  
  for(size_t kk = 0; kk<n; kk++){
    long double y = arr[kk] - c;
    long double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
    
  return (double)sum;
}

/* --- Kahan summation for the square of the elements of a vector --- */

inline double sumarr2(double *arr, int n){

  long double sum = 0.0, c = 0.0;
  
  for(int kk = 0; kk<n; kk++){
    long double y = arr[kk]*arr[kk] - c;
    long double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
    
  return (double)sum;
}


/* -------------------------------------------------------------------------------- */

clm::clm(int ind, int inpar){

  /* --- copy dimensions --- */
  
  nd = ind;
  npar = inpar;
  lwork = -1; // To be used on the first call by LaPack's SVD
  tmp.resize(nd);

  
  /* --- resize fit control arrays --- */

  fcnt.resize(npar);
  diag.resize(npar);
  memset(&fcnt[0], 0, npar*sizeof(clmf));

  
  /* --- 
     Default values of the fit control, 
     can be changed externally 
     --- */
  
  xtol = 1.e-5;       // Controls the minimum relative change to Chi2
  chi2_thres = 1.0;   // Exit inversion if Chi2 is lower than this
  svd_thres = 1.e-13; // Cut-off "relative" thres. for small singular values
  lmax = 1.0e4;       // Maximum lambda value
  lmin = 1.0e-5;      // Minimum lambda value
  lfac = 10.0;        // Change lambda by this amount
  ilambda = 1.0;      // Initial damping parameter for the Hessian giag.
  maxreject = 7;      // Max failed evaluations of lambda.
  error = false;
  regularize = false; // Compute regularization terms in fx ?
  regul_scal = 1.0;   // scale factor for regularization terms
  proc = 0;           // print-out processor number
}

/* -------------------------------------------------------------------------------- */

clm::~clm()
{
  fcnt.clear();
  diag.clear();
}

/* -------------------------------------------------------------------------------- */

void clm::checkParameters(double *x)
{

  /* --- check that parameter values are within limits --- */

  for(int ii = 0; ii<npar; ii++){

    double range = fcnt[ii].limit[1] - fcnt[ii].limit[0];  
    double par = x[ii];
    
    if(fcnt[ii].cyclic){
      if(par < fcnt[ii].limit[0]) par += range;
      if(par > fcnt[ii].limit[1]) par -= range;
    }

    if(fcnt[ii].bouncing){
      if(par < fcnt[ii].limit[0]) par = -par;
      if(par > fcnt[ii].limit[1]) par = fcnt[ii].limit[1] - (par-fcnt[ii].limit[1]);
    }

    if(par <  fcnt[ii].limit[0]) par =  fcnt[ii].limit[0];
    if(par >  fcnt[ii].limit[1]) par =  fcnt[ii].limit[1];

    x[ii] = par;
  }
  
}

/* -------------------------------------------------------------------------------- */

double clm::checkLambda(double lamb)
{
  return max(min(lamb, lmax), lmin);
}

/* -------------------------------------------------------------------------------- */

void clm::checkMaxChange(double *dx, double *x)
{

  /* --- check that corrections to the parameter 
     are within reasonable limits --- */
  
  double maxcha = 0.0;
  for(int ii = 0; ii<npar; ii++){
    if(!fcnt[ii].capped) continue;
    
    if(dx[ii] >= 0) maxcha = fcnt[ii].maxchange[1];
    else maxcha = fcnt[ii].maxchange[0];

    if(fcnt[ii].relchange) maxcha *= fabs(x[ii]);
    

    if(fabs(dx[ii]) > fabs(maxcha)) dx[ii] *= fabs(maxcha / dx[ii]);

  }

}

/* -------------------------------------------------------------------------------- */

void clm::normalizeParameters(double *x)
{
  /* --- Normalize parameters with the given norm --- */

  for(int ii = 0; ii<npar; ii++) x[ii] /= fcnt[ii].scl; 
}

/* -------------------------------------------------------------------------------- */

void clm::scaleParameters(double *x)
{
  /* --- Scale parameters with the given norm --- */

  for(int ii = 0; ii<npar; ii++) x[ii] *= fcnt[ii].scl; 
}

/* -------------------------------------------------------------------------------- */

void clm::zero(double *res, double **rf)
{
  if(res) memset(&res[0], 0, nd*sizeof(double));
  if(rf)  memset(&rf[0][0], 0, nd*npar*sizeof(double));
}
  
/* -------------------------------------------------------------------------------- */

double clm::compute_chi2(double *res, double penalty)
{
  // std::cerr<< (double)sumarr2(res,nd)/double(nd)<<" "<<penalty<<std::endl;
  return (double)sumarr2(res,nd)/double(nd) + penalty*regul_scal;				 
}

/* -------------------------------------------------------------------------------- */

double clm::getChi2Pars(double *res, double **rf, double lambda,
			double *x, double *xnew, void *mydat, clm_func fx, double *dregul)
{

  double newchi2 = 1.e13;
  
  /* --- Get new estimate of the model for the current lambda value --- */
  
  memset(xnew, 0, sizeof(double)*npar);
  compute_trial2(res, rf, lambda, x, xnew, dregul);


  /* --- Evaluate the new model, no response function is needed --- */

  double *new_res = new double [nd]();
  double *new_dregul = NULL;
  if(dregul) new_dregul = new double [nd]();
  
  int status = fx(npar, nd, xnew, new_res, NULL, mydat, new_dregul);
  if(status){
    if(verb)
      fprintf(stderr, "clm::fitdata: [p:%4d] ERROR in the evaluation of FX, aborting inversion\n", proc);
    error = true;
  }else newchi2 = compute_chi2(new_res, (dregul)?new_dregul[npar]:0.0);

  delete [] new_res;
    
    
  /* --- compute Chi2 --- */
  if(dregul) dregul[npar] = new_dregul[npar]; // The derivatives should not be updated here, only the penalty...
  if(new_dregul) delete [] new_dregul;
  return (double)newchi2;
}

/* -------------------------------------------------------------------------------- */

void clm::scaleRF(double **rf)
{
  for(int pp = 0; pp<npar;pp++){
    double scl = fcnt[pp].scl;
    for(int ww=0; ww<nd; ww++)
      rf[pp][ww] *= scl;
  }

}

/* -------------------------------------------------------------------------------- */

double clm::fitdata(clm_func fx, double *x, void *mydat, int maxiter)
{

  /* --- Init variables --- */
  
  double chi2 = 1.e13, ochi2 = 1.e13, bestchi2 = 1.e13, olambda = 0.0, t0 = 0, t1 = 0;
  int iter = 0, nretry = 0;
  bool exitme = false, toolittle = false;
  string rej = "";
  memset(&diag[0],0,npar*sizeof(double));
  error = false;

  double *dregul = NULL;
  if(regularize) dregul = new double [npar+1](); // To store derivatives of regularization terms
  
  
  /* --- Init array for residues and response function --- */
  
  double **rf = mat2d(npar,nd);
  double *res = new double [nd]();
  double *bestpars = new double [nd]();
  double *xnew = new double [nd]();

  
  /* --- check parameters --- */
  
  checkParameters(x);


  /* --- adjust lambda parameter --- */

  double lambda = checkLambda(ilambda);
  t0 = getTime();

  /* --- Evaluate residues and RF, init chi2 --- */

  int status = fx(npar, nd, x, res, rf, mydat, dregul);
  if(status){
    if(verb)
      fprintf(stderr, "clm::fitdata: [p:%4d] ERROR in the evaluation of FX, aborting inversion\n", proc);
    ochi2 = 1.e13;
    error = true;
    iter = maxiter+1;
  }else{
    scaleRF(rf);
    memcpy(&bestpars[0], &x[0], npar*sizeof(double));
    //
    ochi2 = compute_chi2(res, (dregul)?dregul[npar]*regul_scal:0.0);
    bestchi2 = ochi2;
    
    if(verb)
      fprintf(stdout, "[p:%4d, Init] chi2=%f (%f), lambda=%e\n", proc, ochi2,(dregul)?dregul[npar]*regul_scal:0.0 , lambda);
  }

  /* --- Main iterations --- */

  while(iter <= maxiter){

    /* --- Get new estimate of the model and Chi2 for a given lambda parameter.
       The new parameters are already checked for limits and too large corrections 
       inside getChi2Pars 
       --- */
    if(error) break;
    chi2 = getChi2Pars(res, rf, lambda, x, xnew, mydat, fx, dregul);
    if(chi2 != chi2) error = true;
    if(error) break;
    
   

    /* --- Check if we have improved or not --- */

    double reldchi = 2.0 * (chi2 - bestchi2) / (chi2 + bestchi2);
    olambda = lambda;

    if(chi2 < bestchi2){
      
      /* --- Prep lambda for next iter. If we have been increasing lambda,
	 do not decrease it again until next iteration 
	 --- */
      
      if(nretry == 0) lambda = checkLambda(lambda / lfac);
      

      /* --- is the improvements below our threshold? --- */

      if(fabs(reldchi) < xtol){
	if(toolittle) exitme = true;
	else toolittle = true;
      }else toolittle = false;
   
      
      /* --- Store new best guessed model --- */
      
      bestchi2 = chi2;
      memcpy(&bestpars[0], xnew, npar*sizeof(double));
      memcpy(&x[0],        xnew, npar*sizeof(double));
      nretry = 0;

      rej = " ";
      
    }else{
      /* --- Prep lambda for next trial --- */
      
      lambda = checkLambda(lambda * lfac);
      nretry++;
      rej = " *";
      
      if(nretry < maxreject){
	if(verb)
	  fprintf(stderr,"[p:%4d,i:%4d]  ->  chi2=%f (%f), increasing lambda [%e -> %e]\n",
		  proc,iter, chi2, (dregul)?dregul[npar]*regul_scal:0.0,olambda, lambda);
	continue;
      }
	
    }

    /* --- printout --- */
    
    t1 = getTime();

    if(verb)
      fprintf(stderr,"[p:%4d,i:%4d] chi2=%14.5f (%f), dchi2=%e, lambda=%e, elapsed=%5.3fs %s\n",
	      proc,iter, chi2, (dregul)?dregul[npar]*regul_scal:0.0,chi2-ochi2, olambda, t1-t0,rej.c_str());
    
    ochi2 = chi2;

    
    /* --- Check if we are breaking the loop --- */
    
    if(chi2 < chi2_thres){
      if(verb)
	fprintf(stderr,"clm::fitdata: [p:%4d] Chi2 threshold reached [%f] -> chi2=%f\n",
		proc,chi2_thres, bestchi2);
      break;
    }
    
    if(nretry >= maxreject){
      if(verb) 
	fprintf(stderr,"clm::fitdata: [p:%4d] Too many failed attempts, finalizing inversion\n", proc);
      break;
    }

    if(exitme){
      if(verb) 
	fprintf(stderr,"clm::fitdata:  [p:%4d] relative change in chi2 is too low, inversion finished\n", proc);
      break;
    }
    

    /* --- prepare Jacobian for the next iteration --- */
    
    t0 = t1;
    zero(res, rf);
    if(dregul) memset(dregul, 0, npar*sizeof(double));
    //
    status = fx(npar, nd, x, res, rf, mydat, dregul);
    //
    if(status){
      if(verb)
	fprintf(stderr, "clm::fitdata: [p:%4d] ERROR in the evaluation of FX, aborting inversion\n", proc);
      error = true;
      break;
    }
    
    scaleRF(rf);
    iter++;
  }
  
  /* --- Copy results to output array --- */

  memcpy(&x[0], &bestpars[0], npar*sizeof(double));
  
    
  /* --- clean-up --- */

  del_mat(rf);
  delete [] res;
  delete [] bestpars;
  delete [] xnew;
  delete [] dregul;
  
  return (double)bestchi2;
}

/* -------------------------------------------------------------------------------- */

void clm::compute_trial2(double *res, double **rf, double lambda,
			 double *x, double *xnew, double *dregul)
{

  
  /* --- Init Eigen-3 arrays --- */
  
  MatrixXd A(npar, npar);
  VectorXd B(npar);
  Map<VectorXd> RES(xnew, npar);

  
  /* --- 
     compute the curvature matrix and the right-hand side of eq.: 
     A*(dx) = B
     where A = J.T # J
           B = J.T # res
	   and "dx" is the correction to the current model
     --- */

  for(int yy = 0; yy<npar; yy++){
    
    /* --- Compute the Hessian matrix --- */

    for(int xx = 0; xx<=yy; xx++){
      for(int ww = 0; ww<nd; ww++) tmp[ww] = rf[yy][ww] * rf[xx][ww];
      
      A(yy,xx) = sumarr(&tmp[0], nd);

      /* --- Add regularization terms --- */
      
      if(dregul) A(yy,xx) += dregul[xx]*dregul[yy]*regul_scal*regul_scal;
      
      A(xx,yy) = A(yy,xx); // Remember that A is symmetric!
      
    }

    
    /* --- It works better to store the largest diagonal terms
       in this cycle and multiply lambda by this value than the 
       current estimate. I am setting a cap on how much larger
       it can be relative to the current value 
       --- */
    
    diag[yy] = max(A(yy,yy), diag[yy]);
    if(diag[yy] == 0.0) diag[yy] = 1.0;
    //
    double idia = diag[yy], thrfac = 1.0e4;
    if(idia > A(yy,yy)*thrfac) idia = A(yy,yy)*thrfac;
    
    
    /* --- Damp the diagonal of A --- */
    
    A(yy,yy) += lambda * idia;
    //A(yy,yy) *= (1.0 + lambda);
    
    /* --- Compute J^t * Residue --- */
    
    for(int ww = 0; ww<nd; ww++) tmp[ww] = rf[yy][ww] * res[ww];
    B[yy] = sumarr(&tmp[0], nd);
    
  } // yy


  
  /* --- 
     Solve linear system with SVD decomposition and singular value thresholding.
     The SVD is computed with Eigen3 instead of LaPack.
     --- */
  
  JacobiSVD<MatrixXd,ColPivHouseholderQRPreconditioner> svd(A, ComputeThinU | ComputeThinV);
  svd.setThreshold(svd_thres);
  RES = svd.solve(B);

  
  /* --- 
     New estimate of the parameters, xnew = x + dx.
     Check for maximum change, add to current pars and normalize.
     --- */
  
  scaleParameters(xnew);
  checkMaxChange(xnew, x);
  
  for(int ii = 0; ii<npar; ii++) xnew[ii] += x[ii];


  
  /* --- Check that new parameters are within limits --- */
  
  checkParameters(xnew);
  

 
}
