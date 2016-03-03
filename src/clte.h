/*
  CLASS CLTE
  
  PURPOSE: Compute synthetic profiles assuming LTE populations
  
  AUTHOR(S): J. de la Cruz Rodriguez (ISP-SU)
  
  DEPENDENCIES: cprofiles, ceos, input
  
  MODIFICATIONS: -

*/
#ifndef CLTE_H
#define CLTE_H
//
#include <vector>
#include <string>
#include "ceos.h"
#include "cprofiles.h"
#include "input.h"
#include "cmemt.h"
#include "atmosphere.h"
//
class clte: public atmos{
 public:
  static const double lte_const;
  static const double pmax[6];
  static const double pmin[6];
  static const double pscal[6];
  static const double pstep[6];
  std::vector<double> step, isyn;

  std::vector<line_t> lines;
  std::vector<double> lambda;
  int nlambda, nlines, nregions;
  //ceos eos;
  
  /* --- Other objects included --- */
  //ceos eos; // Now ncluded in atmos base class
  cprofiles prof;
  
  /* --- Constructor/Destructor --- */
  // clte(){};
  clte(iput_t &input, double grav = 4.44);
  ~clte(){};

  /* --- methods --- */
 inline double lte_opac(double temp, double n_u, double gf, double elow, double nu0);
 //void synth(mdepth_t &m, mat<double> &syn, cprof_solver sol = bez_z);
 void synth(mdepth &m, double *syn, cprof_solver sol = bez_ltau, bool store_pops = true);
 std::vector<double> get_max_limits(nodes_t &n);
 std::vector<double> get_min_limits(nodes_t &n);
 std::vector<double> get_scaling(nodes_t &n);
 std::vector<double> get_steps(nodes_t &n);
 void cleanup(void){};

 
};


#endif
