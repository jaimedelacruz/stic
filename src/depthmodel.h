#ifndef MDEPTH_H
#define MDEPTH_H
//
#include <vector>
#include <string>
#include <iostream>
#include "io.h"
#include "input.h"
#include "interpol.h"
#include "ceos.h"
#include "cmemt.h"
//#include "atmosphere.h"

//-------------------------------------------
// simple class containig model parameters
//-------------------------------------------s
class mdepth{
 public:
  //-------------------------------------------
  double *temp, *nne, *rho, *pgas, *b,
    *inc, *azi, *v, *vturb, *ltau,
    *z, *pel, *tau, *cmass;
  int ndep;
  mat<double> cub;
  static const double boundary_pgas_default;
  //-------------------------------------------

  
  //-----------------------------------------
  // constructor/destructor
  //-----------------------------------------
  mdepth(){};
  mdepth(int n){
    setsize(n);
  }
  ~mdepth(){
    //  setsize(0);
  }; 

  /* --- prototypes ---*/
  void setsize(int n);
  void nodes2depth(int n, double *x, double *y, int nn, double *xx, double *yy, int interpol = 0);
  void expand(nodes_t &nodes, double *p, int interpol = 0);
  void fixBoundary(int boundary, ceos &eos);
  void getPressureScale(int boundary, ceos &eos);
  void zero(void);
  void fill_densities(ceos &eos);
  mdepth& operator= ( mdepth &m);
//-------------------------------------------
};
typedef mdepth mdepth_t;



/* --- Class MDEPTHALL --- */
class mdepthall{
 public:
  mat<double> temp, nne, rho, pgas, b, inc,
    azi, v, vturb, ltau, z, pel, boundary;
  int ndep, btype;
  ceos eos;
  mat<double> cub;

  
  
 mdepthall(int iny, int inx, int indep): eos(4.44){
    setsize(iny, inx, indep);
  }
  mdepthall(){};
  ~mdepthall(){setsize(0,0,0);};

  void model_parameters (mat<double> &tmp, nodes_t &n, int nt = 1);
  void model_parameters2(mat<double> &tmp, nodes_t &n, int nt = 1);
  int  read_model(std::string &filename, int tstep = 0,  bool require_tau = false);
  int  read_model2(std::string &filename,int tstep = 0,  bool require_tau = false);
  void compress(int n, double *x, double *y, int nn, double *xx, double *yy);
  void compress(int n, float *x, float *y, int nn, double *xx, double *yy);

  void setsize(int ny, int nx, int ndep, bool verbose = true);
  void convertBoundary(int bound, bool verbose = true);
  void expandAtmos(nodes_t &nodes, mat<double> &pars, int interpolation = 0);
  void expand(int n, double *x, double *y, int nn, double *xx, double *yy, int interpolation = 0);
  void write_model(std::string &filename, int tstep = 0);
  void write_model2(std::string &filename, int tstep = 0);

};
typedef mdepthall mdepthall_t;



#endif
