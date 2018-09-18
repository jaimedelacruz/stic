#include <cmath>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <cstdio>
//
#include "solveLinearCXX.h"

using namespace std;
using namespace Eigen;

void solveLinearCXX(int n, double **a, double *b, bool_t improve)
{
  
  /* --- map pointers with Eigen-3 containers --- */
  
  Map<Matrix<double, Dynamic,1>> B(b, n);
  Matrix<double, Dynamic, 1> C = B;
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>> A(&a[0][0], n, n);
  

  /* --- Solve system --- */

  BDCSVD<Matrix<double,Dynamic, Dynamic, RowMajor>> svd(A, ComputeThinU | ComputeThinV);
  C = svd.solve(B);


  /* ---- Improve solution! --- */
  
  if(improve){
    Matrix<double, Dynamic, 1> residue = B - A*C;
    C += svd.solve(residue);
    residue = B-A*C;
  } // if
  
  B = C;
  
}
/*
void solveLinearCXX(int n, double **a, double *b, bool_t improve)
{
  // --- map pointers with Eigen-3 containers --- //
  
  static const double thres = 1.e-15; static const int ITMAX = 10;
  
  Map<Matrix<double, Dynamic, 1>> B(b, n);
  Matrix<double, Dynamic, 1> C = B;
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>> A(&a[0][0], n, n);
  
  
  // --- Solve system --- //

  ColPivHouseholderQR<Matrix<double, Dynamic, Dynamic, RowMajor>> qr(A);
  C = qr.solve(B);
  
  if(thres > 0.0){
    Matrix<double, Dynamic, 1> residue = B - A*C;
    int it = 0;

    while((double(residue.norm()) > thres) && (it++ < ITMAX)){
      C += qr.solve(residue);
      residue = B-A*C;
    } // while
  }

  B = C;

}
*/
