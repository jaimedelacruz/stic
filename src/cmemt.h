#ifndef CMEMT_H
#define CMEMT_H
/*
  Class mem.
  Allows for "virtual" multidimensional arrays and array operations.
  Implemented up to 5 dimensions (at least the indexing part, the rest is ready).
  There is no implementation file because templates must be in a header file.

  AUTHOR: Jaime de la Cruz Rodriguez (ISP-SU 2014).

  CHANGES:
          2014-10-03, JdlCR: Changed interface to use C++ vectors instead of
	                     the old C vectors. This should help to avoid 
			     memory leaks.
		     
	  2014-10-22, JdlCR: Added [=] operator to make copies of the class.
	  
	  2014-10-25, JdlCR: Added [+,-,*,/] operators (not +=, -=, *=, /=),
	                     also made the [=] operator more efficient.

	  2015-04-10, JdlCR: Changed nc to a 'stack' 6 element C array, 
	                     hopefully much easier to optimize for the compiler
			     than a C++ dynamic vector class.
			   
	  2016-11-23, JdlCR: Check for NaNs.

*/

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <cstring>
#include <cmath>
//
template <class T> std::string formatVect(std::vector<T> &in){
  std::stringstream res;
  if(in.size() == 0) return std::string("[ ]");

  res << "["<<in[0];
  for(unsigned int ii = 1; ii<in.size(); ii++) res << ", "<<in[ii];
  res<<"]";
  return res.str();

}
//
template <class T> class mat {
 protected:
  //
  std::vector<int> n;
  //std::vector<long unsigned> nc;
  long unsigned nc[6];
  long unsigned nel;
  //
 public:
  std::vector<T> d;
  bool verbose;



  //
  // Constructor (overloaded)
  //
  mat(std::vector<int> dim) {
    set(dim);
  }
  mat(int n1, int n2){ 
    set({n1,n2});
  }
  mat(int n1, int n2, int n3){
    set({n1,n2,n3}); 
  }
  mat(int n1, int n2, int n3,int n4){
    set( {n1, n2,n3,n4});
  }
  mat(int n1, int n2, int n3,int n4, int n5){
    set({n1,n2,n3,n4,n5});
  }
  mat(){};



  //
  // Destructor
  //
  ~mat() {
    d.clear();
    n.clear();
    // nc.clear();
  }
  

  
  //
  // Init the class 
  //
  void set(int dim){
    set({dim});
  }
  void set(std::vector<int> dim1){
    //
    if((int)dim1.size() > 6){
      std::cerr<<"mat::set: Error, too many dimensions, maximum number is 6"<<std::endl;
      exit(0);
    }
    //
    memset(&nc[0],0,sizeof(unsigned long)*6);
    
    nel = 1;
    int ndim = dim1.size();
    n.resize(ndim);
    //
    for(int dd = 0;dd<ndim;dd++) {
      nel *=  dim1[ndim-dd-1];
      n[ndim-dd-1] = (long int)dim1[ndim-dd-1];
      if(dd == 0) nc[0] = n[ndim-dd-1];
      else nc[dd] = nc[dd-1] * n[ndim-dd-1]; 
    }

    // Resize vector to nel = number of elements
    d.resize(nel, T(0L));
  }


  //
  // Other methods
  //
  void reform(std::vector<int> newdims){
    
    std::string inam = "mat::reform: ";
    int ndim = newdims.size();

    // get new number of elements
    long unsigned nel1 = 1;
    for(auto &it: newdims) nel1 *= it;
    
    // Compare
    if(nel != nel1){
      std::cerr<<inam<<"ERROR, the total number of elements ["<<nel<<
	"] cannot change ["<<nel1<<"]"<<std::endl;
      return;
    }

    // Change n and nc to accomodate new dimensions
    n.resize(ndim, int(0));
    //nc.resize(ndim, (long unsigned)0);
    memset(&nc[0],0,sizeof(unsigned long)*6);

    
    
    for(int dd=0; dd<ndim;dd++){
      n[ndim-dd-1] = (long int)newdims[ndim-dd-1];
      if(dd == 0) nc[0] = n[ndim-dd-1];
      else nc[dd] = nc[dd-1] * n[ndim-dd-1]; 
    }
    
  }
  //
  int size(int idx){
    int ndim = n.size();
    if(abs(idx) < ndim){
      if(idx < 0) return n[ndim+idx];
      else return n[idx];
    } else return 0;
  }
  //
  inline std::vector<int> &getdims(){
    return n;
  }

  // Retrieve number of dimensions
  inline int ndims(){return n.size();}

  // Retrieve total number of elements
  inline long long int n_elements(){return nel;}
  
  // Zero all elements of the array.
  inline void zero(){std::fill(d.begin(), d.end(), 0L);}
  
  // indexing of memory (overloaded)
  inline T &operator()(const int x) {
    return d[x];
  }
  inline T &operator()(const int y, const int x) {
    return d[ y * nc[0] + x];
  }
  inline T &operator()(const int z, const int y, const int x) {
    return d[z * nc[1] + y * nc[0] + x];
  }
  inline T &operator()(const int t, const int z, const int y, const int x) {
    return d[t * nc[2]  + z * nc[1] + y * nc[0] + x];
  }
  inline T &operator()(const int r, const int t, const int z, const int y, const int x) {
    return d[r * nc[3] + t * nc[2]  + z * nc[1] + y * nc[0] + x];
  }
  inline T &operator()(const int w, const int r, const int t, const int z, const int y, const int x) {
    return d[w * nc[4] + r * nc[3] + t * nc[2]  + z * nc[1] + y * nc[0] + x];
  }


  // Equal operator, makes a copy
  /*
  mat<T> operator= ( mat<T> &orig){
    
    mat<T> res;
    res.set(orig.getdims());
    //
    memcpy(&d[0], &orig.d[0], orig.d.size() * sizeof(T));
    
    return res;
  }
  */
    void operator= ( mat<T> &orig){
      set(orig.getdims());
      std::vector<int> dims = orig.getdims();
      set(dims);
      memcpy(&d[0], &orig.d[0], sizeof(T) * orig.d.size());
    }

  
  // Add
  mat<T> operator+ (const mat<T> &b){
    mat<T> res;
    if(b.d.size() != d.size()) return res;
    
    res.set(n);
    for(long unsigned ii = 0; ii<d.size(); ii++) 
      res.d[ii] = d[ii] + b.d[ii];

    return res;
  }
  //Subtract
  mat<T> operator- (const mat<T> &b){
    mat<T> res;
    if(b.d.size() != d.size()) return res;
    
    res.set(n);
    for(long unsigned ii = 0; ii<d.size(); ii++) 
      res.d[ii] = d[ii] - b.d[ii];

    return res;
  }
  //Multiplication
  mat<T> operator* (const mat<T> &b){
    mat<T> res;
    if(b.d.size() != d.size()) return res;
    
    res.set(n);
    for(long unsigned ii = 0; ii<d.size(); ii++) 
      res.d[ii] = d[ii] * b.d[ii];

    return res;
  }
  //Division
  mat<T> operator/ (const mat<T> &b){
    mat<T> res;
    if(b.d.size() != d.size()) return res;
    
    res.set(n);
    for(long unsigned ii = 0; ii<d.size(); ii++) 
      res.d[ii] = d[ii] / b.d[ii];

    return res;
  }


  // Retrieve minimum value in the array
  T min(){
    T mi = d[0];
    for(long unsigned ii = 1; ii < nel;ii++) mi = std::min(mi,d[ii]);
    return mi;
  }

  // Retrieve maximum value in the array
  T max(){
    T ma = d[0];
    for(long long int ii = 1; ii < nel;ii++) ma = std::max(ma,d[ii]);
    return ma;
  }

  // Retrieve the location of the maximum value in the array
  T maxloc(){
    long long int p = 0;
    T ma = d[0];
    for(long long int ii = 1; ii < nel;ii++){
      if(d[ii] > ma) {
	p = ii;
	ma = d[ii];
      }
    }
    return p;
  }
  
  // Retrieve the location of the minimum value in the array
  long long int minloc(){
    long long int p = 0;
    T mi = d[0];
    for(long unsigned ii = 1; ii < nel;ii++){
      if(d[ii] < mi) {
	p = ii;
	mi = d[ii];
      }    
    }
    return p;
  }

  // Perform the sum of all the elements of the array
  inline double sum(){
    double s = 0.0L;
    for(auto &it: d) s += it;
    return s;
  }

  bool isNaN() const{
    for(size_t ii=0; ii<nel; ii++) if(std::isnan(d[ii])) return true;
    return false;
  }

  void fillNaN(const T val = 0){
    for(size_t ii=0; ii<nel; ii++) if(std::isnan(d[ii])) d[ii] = 0.0;
  }

};
 
#endif
