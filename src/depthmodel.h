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
#include "eoswrap.h"
#include "cmemt.h"

//-------------------------------------------
// simple class containig model parameters
//-------------------------------------------s
class mdepth{
 public:
  //-------------------------------------------
  double *temp, *nne, *rho, *pgas, *bl,
    *bh, *azi, *v, *vturb, *ltau,
    *z, *pel, *tau, *cmass;
  double tr_loc, tr_amp;
  int ndep;
  mat<double> cub;
  double bound_val;
  static const double boundary_pgas_default;
  mdepth *ref_m;
  //-------------------------------------------

  
  //-----------------------------------------
  // constructor/destructor
  //-----------------------------------------
  mdepth(){ref_m = NULL; ndep = 0;};
  mdepth(int n){
    ref_m = NULL;
    setsize(n);
  }
  mdepth(const mdepth &in);
  
  
  ~mdepth(){
    //  setsize(0);
  }; 

  /* --- prototypes ---*/
  void setsize(int n);
  void nodes2depth(int n, double *x, double *y, int nn, double *xx, double *yy, int interpol, bool extrapolate);
  void expand(nodes_t &nodes, double *p, int interpol, int mtype = 0);
  void fixBoundary(int boundary, eoswrap &eos);
  void nne_enhance(nodes_t &nodes, int n, double *pars, eoswrap &eos);
  void getPressureScale(int depth_t, int boundary, eoswrap &eos);
  void getScales(eoswrap &eos, int bound);
  void zero(void);
  void fill_densities(eoswrap &eos, int keep_nne, int k0, int k1);
  void hydrostatic(eoswrap &eos, int depth_t);
  void optimize_depth(eoswrap &eos, float tcut, int smooth = 0);
  void optimize_depth_ltau(eoswrap &eos, float tcut);
  void to_txt(const std::string &fname);
  double* compress(nodes_t &nodes)const;

  
  mdepth& operator= ( mdepth &m);

  std::vector<double> model2vector();
  void vector2model(std::vector<double> &vec);
  
//-------------------------------------------
};
typedef mdepth mdepth_t;



/* --- Class MDEPTHALL --- */
class mdepthall{
 public:
  mat<double> temp, nne, rho, pgas, bl, bh,
    azi, v, vturb, ltau, z, pel, boundary, tr_amp, tr_loc;
  mat<int> tr_N;
  int ndep, btype;
  // ceos eos;
  mat<double> cub;

   
  
  mdepthall(int iny, int inx, int indep){//: eos(4.44){
    setsize(iny, inx, indep);
  }
  mdepthall(){};
  ~mdepthall(){setsize(0,0,0);};

  void model_parameters (mat<double> &tmp, nodes_t &n, int nt = 1);
  void model_parameters2( mat<double> &tmp, nodes_t &n, int nt = 1);
  int  read_model2(iput_t const& input, std::string &filename,int tstep = 0,  bool require_tau = false);
  void compress(int n, double *x, double *y, int nn, double *xx, double *yy);
  void compress(int n, float *x, float *y, int nn, double *xx, double *yy);

  void setsize(int ny, int nx, int ndep, bool verbose = true);
  void convertBoundary(int bound, bool verbose = true);
  void expandAtmos(nodes_t &nodes, mat<double> &pars, int interpolation = 0);
  void expand(int n, double *x, double *y, int nn, double *xx, double *yy, int interpolation = 0);
  void write_model(std::string &filename, int tstep = 0);
  void write_model2(iput_t const& input, std::string &filename, int tstep = 0);

};
typedef mdepthall mdepthall_t;



#endif
