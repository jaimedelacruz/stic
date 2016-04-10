/* 
   Implementation of the clm class (Levenberg-Marquardt least-squares-fit)
   author: J. de la Cruz Rodriguez (Stockholm University 2016)
*/

#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <cstring>
#include "clm.h"

using namespace std;

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

inline double sumarr(double *arr, int n){

  double sum = 0.0, c = 0.0;
  
  for(int kk = 0; kk<n; kk++){
    double y = arr[kk] - c;
    double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
    
  return (double)sum;
}

/* --- Kahan summation for the square of the elements of a vector --- */

inline double sumarr2(double *arr, int n){

  double sum = 0.0, c = 0.0;
  
  for(int kk = 0; kk<n; kk++){
    double y = arr[kk]*arr[kk] - c;
    double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
    
  return sum;
}


/* -------------------------------------------------------------------------------- */

clm::clm(int ind, int inpar){

  /* --- copy dimensions --- */
  
  nd = ind;
  npar = inpar;
  lwork = -1; // To be used on the first call by LaPack's SVD

  
  /* --- resize fit control arrays --- */

  fcnt.resize(npar);
  diag.resize(npar);

  
  /* --- 
     Default values of the fit control, 
     can be changed externally 
     --- */
  
  xtol = 1.e-5;       // Controls the minimum relative change to Chi2
  chi2_thres = 1.0;   // Exit inversion if Chi2 is lower than this
  svd_thres = 1.e-16; // Cut-off "relative" thres. for small singular values
  lmax = 1.0e6;       // Maximum lambda value
  lmin = 1.0e-12;     // Minimum lambda value
  lfac = 10.0;        // Change lambda by this amount
  ilambda = 0.1;      // Initial damping parameter for the Hessian giag.
  maxreject = 6;      // Max failed evaluations of lambda.
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

void clm::checkMaxChange(double *dx)
{

  /* --- check that corrections to the parameter 
     are within reasonable limits --- */
  
  for(int ii = 0; ii<npar; ii++){
    if(!fcnt[ii].capped) continue;
    
    if(fabs(dx[ii]) > fcnt[ii].maxchange)
      dx[ii] *= fcnt[ii].maxchange / fabs(dx[ii]);
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

double clm::compute_chi2(double *res)
{
  return (double)sumarr2(res,nd)/double(nd);				 
}

/* -------------------------------------------------------------------------------- */

double clm::getChi2Pars(double *res, double **rf, double lambda,
			double *x, double *xnew, void *mydat, clm_func fx)
{

  /* --- Get new estimate of the model for the current lambda value --- */
  
  memset(xnew, 0, sizeof(double)*npar);
  compute_trial3(res, rf, lambda, x, xnew);


  /* --- Evaluate the new model, no response function is needed --- */

  double *new_res = new double [nd]();
  
  fx(npar, nd, xnew, new_res, NULL, mydat);
  double newchi2 = compute_chi2(new_res);

  delete [] new_res;
    
    
  /* --- compute Chi2 --- */
    
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
  
  double chi2 = 1.e13, ochi2 = 1.e13, bestchi2 = 1.e13, olambda = 0.0;
  int iter = 0, nretry = 0;
  bool exitme = false;
  string rej = "";
  memset(&diag[0],0,npar*sizeof(double));
  
  
  /* --- Init array for residues and response function --- */
  
  double **rf = mat2d(npar,nd);
  double *res = new double [nd]();
  double *bestpars = new double [nd]();
  double *xnew = new double [nd]();

  
  /* --- check parameters --- */
  
  checkParameters(x);


  /* --- adjust lambda parameter --- */

  double lambda = checkLambda(ilambda);


  /* --- Evaluate residues and RF, init chi2 --- */

  fx(npar, nd, x, res, rf, mydat);
  scaleRF(rf);
  memcpy(&bestpars[0], &x[0], npar*sizeof(double));
  //
  ochi2 = compute_chi2(res);
  bestchi2 = ochi2;

  if(verb)
    fprintf(stdout, "[Init] chi2=%f, lambda=%e\n", ochi2, lambda);


  /* --- Main iterations --- */

  while(iter <= maxiter){

    /* --- Get new estimate of the model and Chi2 for a given lambda parameter.
       The new parameters are already checked for limits and too large corrections 
       inside getChi2Pars 
       --- */
    
    chi2 = getChi2Pars(res, rf, lambda, x, xnew, mydat, fx);


    /* --- Check if we have improved or not --- */

    double reldchi = 2.0 * (chi2 - bestchi2) / (chi2+bestchi2);
    olambda = lambda;

    if(chi2 < bestchi2){
      
      /* --- Prep lambda for next iter --- */
      
      lambda = checkLambda(lambda / lfac);

      
      /* --- If the relative change was too small, exit --- */
      
      if((fabs(reldchi) < xtol)) exitme = true;

      
      /* --- Store new best guessed model --- */
      
      bestchi2 = chi2;
      memcpy(bestpars, xnew, npar*sizeof(double));
      memcpy(x, xnew, npar*sizeof(double));
      nretry = 0;
      rej = " ";
      
    }else{
      /* --- Prep lambda for next trial --- */
      lambda = checkLambda(lambda * lfac);
      nretry++;
      rej = " *";
      
      if(nretry < maxreject){
	if(verb)
	  fprintf(stderr,"[%4d]  ->  chi2=%f, increasing lambda [%e -> %e]\n",
		  iter, chi2, olambda, lambda);
	continue;
      }
	
    }

    /* --- printout --- */

    if(verb)
      fprintf(stderr,"[%4d] chi2=%14.5f, dchi2=%e, lambda=%e %s\n",
	      iter, chi2, chi2-ochi2, olambda, rej.c_str());
    
    ochi2 = chi2;
    

    
    /* --- Check if we are breaking the loop --- */
    
    if(chi2 < chi2_thres){
      if(verb)
	fprintf(stderr,"clm::fitdata: Chi2 threshold reached [%f] -> chi2=%f\n",
		chi2_thres, bestchi2);
      break;
    }
    
    if(nretry >= maxreject){
      if(verb) 
	fprintf(stderr,"clm::fitdata: Too many failed attempts, finalizing inversion\n");
      break;
    }

    if(exitme){
      if(verb) 
	fprintf(stderr,"clm::fitdata: relative change in chi2 is too low, inversion finished\n");
      break;
    }
    

    /* --- prepare Jacobian for the next iteration --- */
    
    zero(res, rf);
    fx(npar, nd, x, res, rf, mydat);

    scaleRF(rf);
    iter++;
  }
  
  /* --- Copy results to output array --- */

  memcpy(x, bestpars, npar*sizeof(double));
  
    
  /* --- clean-up --- */

  del_mat(rf);
  delete [] res;
  delete [] bestpars;
  delete [] xnew;

  return (double)bestchi2;
}


/* -------------------------------------------------------------------------------- */

extern "C"{

  /* --- Prototypes for the LAPACK routines --- */
  
  void dgesdd_( char *JOBZ, int &M, int &N, double *A, int &LDA, double *S,
  		double *U, int &LDU, double *VT, int &LDVT, double *WORK,
  		int &lwork, int* iwork, int &info);
  
  void dgelsd_( int &M, int &N, int &NRHS, double *A, int &LDA, double *B,
		int &LDB, double *S, double &RCOND, int &RANK,
                double *work, int &lwork, int *iwork, int &info );
}

/* -------------------------------------------------------------------------------- */

void clm::compute_trial3(double *res, double **rf, double lambda,
			double *x, double *xnew)
{

  
  /* --- Init arrays --- */
  
  double **A = mat2d(npar, npar);
  double **V = mat2d(npar, npar);

  double *B = new double [npar];
  double *w = new double [npar];
  int   *iw = new int  [8*npar]();
  
  double *work = NULL;
  double  *tmp = new double [nd];

  
  
  /* --- 
     compute the curvature matrix and the right-hand side of eq.: 
     A*(dx) = B
     where A = J.T # J
           B = J.T # res
	   and "dx" is the correction to the current model
     --- */

  for(int yy = 0; yy<npar; yy++){
    for(int xx = 0; xx<=yy; xx++){
      
      memset(tmp, 0, nd*sizeof(double));
      for(int ww = 0; ww<nd; ww++)
	tmp[ww] = rf[yy][ww] * rf[xx][ww];
      
      A[yy][xx] = sumarr(tmp, nd);
    }
      
    
    /* --- Compute gradient of chi2 --- */
    
    memset(tmp, 0, nd*sizeof(double));
    for(int ww = 0; ww<nd; ww++) tmp[ww] = rf[yy][ww] * res[ww];
    B[yy] = sumarr(tmp, nd);
    
  }

  
  /* Fill in the other half of A (symmetric)--- */

  for(int yy = 0; yy<npar; yy++)
    for(int xx = 0; xx<yy; xx++)
      A[xx][yy] = A[yy][xx];


  /* --- Damp the diagonal of A --- */
  
  for(int yy = 0; yy<npar; yy++){

    diag[yy] = max(A[yy][yy], diag[yy]);
    if(diag[yy] == 0.0) diag[yy] = 1.0;
    
    A[yy][yy] += lambda * (diag[yy]*0.5 + A[yy][yy]*0.5);

  }

  
  /* --- Solve a least square fit to ||B - AX||^2 (internally using SVD) --- */
  
  char bla[1] = {'O'};
  int dummy = 0;
  int npp = npar;
  
  
  /* --- Request size of work array --- */
  
  if(lwork == -1){
    lwork = -1;
    double wrkopt=0.0;
    dgesdd_(&bla[0], npp, npp, &A[0][0], npp, w, NULL, npp, &V[0][0] ,
	    npp, &wrkopt, lwork, iw, dummy);
    lwork = (int)wrkopt;
  }
  
  
  /* --- Call DGESDD --- */
  
  work = new double [lwork];
  dummy = 0;
  dgesdd_(&bla[0], npp, npp, &A[0][0], npp, w, NULL, npp, &V[0][0] ,
	  npp, work, lwork, iw, dummy);
  
  
  /* --- Check for too small singular values --- */
  
  double minval = w[0] * svd_thres;
  for(int ii = 1; ii < npar; ii++) if(w[ii] < minval) w[ii] = 0.0;
  
  
  
  /* --- Solve the system, xnew at this point is the correction to x, 
     not the new parameters. We should do xnew += x. --- */
  
  backsub(A,w,V,npar,B,xnew);
  

  
  /* --- check for maximum change, add to current pars and normalize --- */
  
  
  scaleParameters(xnew);
  checkMaxChange(xnew);

  for(int ii = 0; ii<npar; ii++) xnew[ii] += x[ii];
  
  checkParameters(xnew);

  
  /* --- Clean-up mem --- */
  
  del_mat(A);
  del_mat(V);
  delete [] B;
  delete [] w;
  delete [] tmp;
  delete [] work;
  delete [] iw;
}

/* -------------------------------------------------------------------------------- */

void clm::backsub(double **u, double *w, double **v, int n, double *b, double *x)
{
  //
  // Ax = B -> A^-1 * B = VT W^-1 U *B
  // matrices indexes are a bit weird so it works with LAPACK's routines
  // JdlCR MODIFICATION: Added more accurate routine to add the product of rows/columns.
  //
  
  double *tmp = new double [n]();
  double *tmp1 = new double [n]();
  
  for(int jj=0;jj<n;jj++){    
    if(fabs(w[jj]) == 0.0) continue;
    
    double sum = 0.0;
    for(int ii = 0; ii<n; ii++)
      tmp1[ii]= u[jj][ii]*b[ii];
    
    tmp[jj] = sumarr(tmp1,n)/w[jj];
  }
  
  for(int j=0;j<n;j++){    
    double sum = 0.0;
    for(int jj=0;jj<n;jj++)
      tmp1[jj]= tmp[jj]*v[j][jj];
    x[j] = sumarr(tmp1,n);
  }  

  delete [] tmp;
  delete [] tmp1;
}
