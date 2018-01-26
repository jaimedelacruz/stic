#ifndef INPUT_H
#define INPUT_H
#include <string>
#include <vector>
#include <iostream>
#include "cmemt.h"
//
/* --- line list struct --- */
struct line{
  char elem[8], label[15]; // remember that last char is for the NULL termination
  double Jup, Jlow, Gup, Glow;
  double w0, nu0, width, eion;
  double gf, e_low, e_up, amass;
  double g_rad, g_str, g_vdw;
  double b_sig, b_alp, b_vbar, b_gvw;
  int anum, ion, firsttime, off;
  // Zeeman splitting
  int nZ;
  std::vector<double> strength, splitting;
  std::vector<int> iL;
};
typedef line line_t;


/* --- region struct --- */
struct region{
  double w0, dw, cscal;
  int nw, off;
  std::vector<double> wav, nu;
  std::vector<int> idx;
  std::string inst, ifile;
};
typedef region region_t;


/* --- node types ---*/
enum nodes_type_t{
  none_node,
  temp_node,
  v_node,
  vturb_node,
  bl_node,
  bh_node,
  azi_node,
  pgas_node
};


/* --- nodes struct --- */
struct nodes{
  int nnodes, temp_off, v_off, vturb_off, bl_off, bh_off, azi_off, pgas_off, tosend, bound,depth_t;
  std::vector<double> temp;
  std::vector<double> v;
  std::vector<double> vturb;
  std::vector<double> bl;
  std::vector<double> bh;
  std::vector<double> azi;
  std::vector<nodes_type_t> ntype;
  int toinv[7];
  int regul_type[7];
};
typedef nodes nodes_t;



/* --- input structure --- */
struct iput{
  unsigned long buffer_size, buffer_size1;
  int nt, ny, nx, ns, npar, npack, mode, nInv, inst_len, atmos_len, ab_len,
    nw_tot, boundary, ndep, solver, centder, thydro, dint, keep_nne, svd_split, random_first, depth_model,
    use_geo_accel, nresp, getResponse[8];
  double mu, chi2_thres, sparse_threshold, dpar, init_step, marquardt_damping, svd_thres, regularize, tcut;
  std::string imodel, omodel, iprof, oprof, myid, instrument,
    atmos_type, wavelet_type, oatmos, abfile;
  int xx, yy, ipix, nPacked;
  std::vector<double> chi;
  int myrank, nprocs, cgrad;
  unsigned max_inv_iter, master_threads, wavelet_order;
  std::vector<unsigned long> ntosend;
  std::vector<std::string> ilines;
  std::vector<region_t> regions;
  std::vector<line_t> lines;
  nodes_t nodes;
  bool verbose;
};
typedef iput iput_t; 


/* --- Functions --- */
std::vector<std::string> strsplit(std::string &var, std::string token, bool rmspaces = true);
std::string removeSpaces(std::string input);
iput_t read_input(std::string filename, bool verbose = 1);
void read_lines(std::string filename, iput_t &input, bool verbose = 1);
std::vector<double> fill_lambdas(iput_t &input, bool air = false);


void equidist(std::vector<double> &var, double min, double max);
void equidist(std::vector<double> &var, std::vector<double> &itau, bool cent_grid = false);
double nodeLocation(std::vector<double> &itau, double iloc);
//int set_nodes(nodes_t &n, double min, double max, bool verbose = false);
int set_nodes(nodes_t &n, std::vector<double> &itau, int dint, bool verbose = false);

/*--- Convert from air -> vacuum and viceversa --- */
double convl(double lambda);
double inv_convl(double lambda);


//
#endif
