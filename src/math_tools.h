/* -------------------------------------
   
  Various math routines in template format

  Coded by J. de la Cruz Rodriguez (ISP-SU 2016)

  -------------------------------------- */

#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>

namespace mth{

  
  /* --- sqr --- */
  
  template <class T> T sqr(T a){return a*a;};
  /* ------------------------------------------------------------------------------- */

  /* --- Sum, mean, stdev --- */
  
  template <class T> T sum(size_t n, T* a){
    T tot = 0;
    for(size_t ii=0;ii<n;ii++) tot += a[ii];
    return tot;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> T mean(size_t n, T* a){
    T tot = sum<T>(n, a);
    return tot/n;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> double stdev(size_t n, T* a){
    double res = 0.0;
    double a_mean = sum<T>(n, a) / (double)n;
    //
    for(size_t ii=0; ii<n; ii++) res += sqr<double>((double)a[ii]-a_mean);
    //
    return sqrt(res/(n-1.0));
  }
  /* ------------------------------------------------------------------------------- */
  
  /* --- more accurate Kahan addition --- */
  
  template <class T> double ksum(size_t n, T* arr){
    
    long double sum = 0.0L, c = 0.0L;
    
    for(size_t kk = 0; kk<n; kk++){
      long double y = arr[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */

  template <class T> double ksum2(size_t n, T* arr){
    
    long double sum = 0.0L, c = 0.0L;
    
    for(size_t kk = 0; kk<n; kk++){
      long double y = arr[kk]*arr[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */

  /* --- Dot product of vectors --- */

  template <class T> double dot(size_t n, T *a, T *b){
    double res = 0.0;
    for(size_t ii=0; ii<n; ii++) res += a[ii]*b[ii];
    return res;
  }
  /* ------------------------------------------------------------------------------- */

  template <class T> double kdot(size_t n, T *a, T* b){
    long double sum = 0.0L, c = 0.0L;
    for(size_t kk = 0; kk<n; kk++){
      long double y = a[kk]*b[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */
  

  
} // namespace

#endif
