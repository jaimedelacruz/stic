#include <math.h>
//#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* --- 
   Solve Linear systems with SVD for rank defficient matrices
   
   Author: J. de la Cruz Rodriguez (ISP-SU 2016)

   Dependencies: LaPack

   --- */


/* -------------------------------------------------------------------------------- */

/* --- Prototypes for the LAPACK routines --- */

extern void dgesdd_( char *JOBZ, int *M, int *N, double *A, int *LDA, double *S,
  		double *U, int *LDU, double *VT, int *LDVT, double *WORK,
  		int *lwork, int* iwork, int *info);

/* -------------------------------------------------------------------------------- */

double **mat2d(int nx1, int nx2){
  double **p;
  int x1 = 0;
  
  p = (double**)malloc(nx1*sizeof(double*));
  p[0] = (double*)calloc(nx1 * nx2, sizeof(double));
  
  for(x1=1;x1<nx1;++x1) p[x1] = p[x1-1] + nx2;
  
  return p;
}

/* -------------------------------------------------------------------------------- */

void del_mat(double **p){
  free(p[0]);
  free(p);
}

/* -------------------------------------------------------------------------------- */

inline double sumarr(double *arr, int n){

  double sum = 0.0, c = 0.0, y = 0.0, t = 0.0;
  int kk = 0;
  
  for(kk = 0; kk<n; kk++){
    y = arr[kk] - c;
    t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
    
  return (double)sum;
}

/* -------------------------------------------------------------------------------- */

void backsub(double **u, double *w, double **v, int n, double *b)
{
  /* --- 
     Ax = B -> A^-1 * B = VT W^-1 U *B
     JdlCR: Matrices indexes are a bit weird so it works with LAPACK's routines.
            In LaPack both U and V are transposed compared to standard C notation.
            Added more accurate routine to add the product of rows/columns.

     The result is returned in array b

  --- */
  
  int jj, ii, j;
  double sum = 0.0;
  double *tmp  = (double*)calloc(n, sizeof(double));
  double *tmp1 = (double*)calloc(n, sizeof(double));

  for(jj=0;jj<n;jj++){    
    if(fabs(w[jj]) == 0.0) continue;
    
    sum = 0.0;
    for(ii = 0; ii<n; ii++)
      tmp1[ii]= u[jj][ii]*b[ii];
    
    tmp[jj] = sumarr(tmp1,n)/w[jj];
  }
  
  for(j=0;j<n;j++){    
    sum = 0.0;
    for(jj=0;jj<n;jj++)
      tmp1[jj]= tmp[jj]*v[j][jj];
    b[j] = sumarr(tmp1,n);
  }  


  free(tmp);
  free(tmp1);
}

/* -------------------------------------------------------------------------------- */


void SolveLinearSvd(int n, double **A_in, double *b)
{

  /* --- Allocate temporary arrays --- */
  
  double **V = (double**)mat2d(n,n);
  double **A = (double**)mat2d(n,n);
  double *w = (double*)malloc(n*sizeof(double));
  int *iw = (int*)malloc(8*n*sizeof(int));
  double *work = NULL;
  int lwork = -1;
  double wrkopt = 0.0, thres = 1.e-15;
  int npp = 0, dummy = 0, i = 0, j=0;
  char bla[1] = {'O'};


  /* --- Copy A --- */
  
  //memcpy(&A[0][0], &A_in[0][0], n*n*sizeof(double));
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)
      A[j][i] = A_in[i][j];

  
  /* --- Request work space from lapack --- */

  npp = n;
  dgesdd_(&bla[0], &npp, &npp, &A[0][0], &npp, w, NULL, &npp, &V[0][0],
	  &npp, &wrkopt, &lwork, iw, &dummy);
  lwork = (int)wrkopt;

  

  /* --- Allocate work space --- */

  work = (double*)malloc(((int)lwork)*sizeof(double));
  dummy = 0;


  
  /* --- Perform SVD, "divide and conquer" algorithm from LaPack --- */

  dgesdd_(&bla[0], &npp, &npp, &A[0][0], &npp, w, NULL, &npp, &V[0][0] ,
	  &npp, work, &lwork, iw, &dummy);
  


  /* --- look for very small singular values --- */

  thres *= w[0];
  for(i = 1; i<n; i++){
    if(w[i] <= thres) w[i] = 0.0;    
  }
  
  
  /* --- Back-substitution --- */

  backsub(A,w,V,n,b);

  

  /* --- clean-up --- */

  free(w);
  free(iw);
  free(work);
  del_mat(V);
  del_mat(A);
}

/* -------------------------------------------------------------------------------- */
