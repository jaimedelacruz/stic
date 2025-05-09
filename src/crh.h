
#ifndef CRH_H
#define CRH_H

#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include "cmemt.h"
#include "input.h"
#include "depthmodel.h"
#include "cprofiles2.h"
#include "atmosphere.h"
#include "rhf1d.h"

class crh: public atmos{
 public:
  static const double pmax[9];
  static const double pmin[9];
  static const double pscal[9];
  static const double pstep[9];

  // iput_t input;
  int nlambda, nlines, nregions;
  std::vector<double> lambda, cmass, nhtot;
  crhpop save_pop;
  
  /* --- Prototypes --- */
  std::vector<double> get_max_limits(nodes_t &n, int mode =1);
  std::vector<double> get_min_limits(nodes_t &n, int mode = 1);
  std::vector<double> get_scaling(nodes_t &n, int mode =1);
  std::vector<double> get_steps(nodes_t &n, int mode = 1);
  bool synth(mdepth &m, double *syn, int computing_derivatives = 0, cprof_solver sol = bez_ltau, bool store_pops = true);
  void cleanup();
  void lambdaIDX(int nw, double *lambda);
  void checkBounds(mdepth_t &m);


  
  /* --- Contructor / Destructor --- */
  crh(iput_t &iput, double grav = 4.44);
  ~crh();
  

  
  
  
  
};
#endif
