/* ---
   This is a base class for the RT type. In the end
   a lot of functions ended up here, and probably they should be moved elsewhere.

   The 1D pixel to pixel fitting routines could be part of the depthmodel base class
   instead of the RT base class.

   Regularization could also be moved.

   Coded by J. de la Cruz Rodriguez (ISP-SU 2016/2017)

   --- */
#include <vector>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <random>
#include "atmosphere.h"
#include "spectral.h"
#include "ceos.h"
#include "input.h"
#include "physical_consts.h"
#include "interpol.h"
#include "clm.h"
#include "math_tools.h"
#include "mmem.h"
//
using namespace std;
//
const double atmos::maxchange[7] = {7000., 7.0e5, 5.0e5, 600., 600., phyc::PI/5, 0.6};

inline double SQ(const double a){return a*a;};
inline double CUB(const double a){return a*a*a;};

vector<double> atmos::get_max_change(nodes_t &n, int mode){
  
  int nnodes = (int)n.nnodes, ntype = n.ntype.size();
  if(nnodes != ntype) return maxc;
  if(nnodes > 0){
    maxc.resize(nnodes);
    
    for(int k = 0; k<nnodes; k++){
      if     (n.ntype[k] == temp_node ) maxc[k] = maxchange[0];
      else if(n.ntype[k] == v_node    ) maxc[k] = maxchange[1];
      else if(n.ntype[k] == vturb_node) maxc[k] = maxchange[2];
      else if(n.ntype[k] == bl_node    ) maxc[k] = maxchange[3];
      else if(n.ntype[k] == bh_node  ) maxc[k] = maxchange[4];
      else if(n.ntype[k] == azi_node  ) maxc[k] = maxchange[5];
      else if(n.ntype[k] == pgas_node ) maxc[k] = maxchange[6];
      else                              maxc[k] = 0;
    }
  }
  
  return maxc;
}



void atmos::responseFunction(int npar, mdepth_t &m_in, double *pars, int nd, double *out, int pp, double *syn){
  
  if(npar != input.nodes.nnodes){
    fprintf(stderr,"atmos::responseFunction: %d (npar) != %d (nodes.nnodes)\n", npar, input.nodes.nnodes);
    return;
  }

  mdepth m = m_in;
  // mdepth m(m_in.ndep);
  //m.cub.d = m_in.cub.d;
  //m.bound_val = m_in.bound_val;
  
  bool store_pops = false;
  int centder = input.centder; 
  //if(input.nodes.ntype[pp] == v_node) centder = 1;

  
  /* --- Init perturbation --- */
  
  double *ipars = new double [npar];
  memcpy(&ipars[0], &pars[0], npar*sizeof(double));
  memset(&out[0], 0, nd*sizeof(double));
  //
  double pval = ipars[pp];
  double pertu = 0.0;
    
  //
  if(input.depth_model == 0){
    if((input.nodes.ntype[pp] == temp_node) && (pval > 10000.)){
      pertu = input.dpar * pval * 0.25;
    }else
      pertu = input.dpar * scal[pp] * 0.3;
  }else  pertu = input.dpar * scal[pp];
    
  
  /* --- Centered derivatives ? --- */
  
  if(centder > 0){
    
    /* --- Up and down perturbations --- */
    double up = pertu * 0.5;
    double down = -0.5 * pertu;

    
    /* --- Init temporary vars for storing spectra --- */
    double *spec = new double [nd]();
    
    
    /* --- Compute spectra with perturbations on both sides --- */

    /* --- Upper side perturb --- */
    
    if((pval + up) > mmax[pp]){
      up = 0.0;
      ipars[pp] = pval;
      memcpy(&out[0], &syn[0], nd*sizeof(double));
    } else{
      ipars[pp] = pval + up;
      m.expand(input.nodes, &ipars[0], input.dint, input.depth_model);
      checkBounds(m);

      /* -- recompute Hydro Eq. ? --- */
      
      //if(input.nodes.ntype[pp] == temp_node && input.thydro == 1)
	m.getPressureScale(input.nodes.depth_t, input.boundary, eos);
	//m.nne_enhance(input.nodes, npar, &ipars[0], eos);

	synth(m, &out[0], 1, (cprof_solver)input.solver, store_pops);
    }

    
    /* --- Lower side perturb --- */

    if((pval + down) < mmin[pp]){
      down = 0.0;
      ipars[pp] = pval;
      memcpy(&spec[0], &syn[0], nd*sizeof(double));
    } else{
      ipars[pp] = pval + down;
      m.expand(input.nodes, &ipars[0], input.dint, input.depth_model);
      checkBounds(m);

      // if(input.nodes.ntype[pp] == temp_node && input.thydro == 1)
      m.getPressureScale(input.nodes.depth_t, input.boundary, eos);
	//m.nne_enhance(input.nodes, npar, &ipars[0], eos);

      synth(m, &spec[0], 1, (cprof_solver)input.solver, store_pops);
    }

    /* --- Compute finite difference --- */
    
    ipars[pp] = pval;
    pertu = (up - down);
    if(pertu != 0.0) pertu = 1.0 / pertu;
    else pertu = 0.0;

    
    /* --- Finite differences --- */
    
    for(int ii = 0; ii<nd; ii++)
      out[ii] =  (out[ii] - spec[ii]) * pertu;
    
    
    /* --- Clean-up --- */
    
    delete [] spec;
    
  } else{ // Sided-derivative
    
    /* --- One sided derivative --- */
    
    if((ipars[pp]+pertu) < mmax[pp]){
      
      /* --- Up --- */

      ipars[pp] += pertu;
    }else{
      
      /* --- Down --- */

      ipars[pp] -= pertu;
    }

    pertu = (ipars[pp] - pval);
    //fprintf(stderr,"pval=%e, pertu=%e, ipars=%e\n", pval, pertu, ipars[pp]);
    
    /* --- Synth --- */
    
    m.expand(input.nodes, &ipars[0],  input.dint, input.depth_model);
    checkBounds(m);

    //  if((input.nodes.ntype[pp] == temp_node) && (input.thydro == 1))
    m.getPressureScale(input.nodes.depth_t, input.boundary, eos);
    //m.nne_enhance(input.nodes, npar, &ipars[0], eos);
      
    synth(m, &out[0], 1, (cprof_solver)input.solver, store_pops);
    
    /* --- Finite differences ---*/
    
    for(int ii = 0; ii<nd; ii++)
      out[ii]  =  (out[ii] - syn[ii]) /  pertu;
  }
  
  
  /* --- clean-up --- */
  
  delete [] ipars;
    
}

void atmos::randomizeParameters(const nodes_t &n, int npar, double *pars, const int rvel){

  static const double vtau[4] = {-7.0, -5.0, -3.0, 1.0};
  static const double vvel[4] = {10, 6.0, 3.0, 0.0};
  std::vector<double> res(n.v.size(), 0.0);
  static bool firsttime = true;
  
  static std::mt19937 rng;
  if(firsttime){
    rng.seed(std::random_device()());
  }
  static std::uniform_int_distribution<std::mt19937::result_type> rand_dist(0,1.e4);
  firsttime = false;
  
  
  int nvel = (int)n.v.size();
  if(nvel > 1){
    hermpol<double>(4, vtau, vvel, nvel, &n.v[0], &res[0]);
  }

  
  //srand (time(NULL));
  int ninit = (int)scal.size();
  
  if(npar != ninit){
    fprintf(stderr,"atmos::randomizeParameters: atmos pars[%d] != scal[%d] -> not randomizing!\n", npar, ninit);
    return;
  }

  nodes_type_t last = none_node;
  double pertu, rnum=0.3;
  
  for(int pp = 0; pp<npar; pp++){

    /* --- Only apply a constant perturbation --- */
    if(last != n.ntype[pp]){
      rnum = rand_dist(rng)*1.e-4;//(double)rand() / RAND_MAX;
      last = n.ntype[pp];
    }

    if(n.ntype[pp] == temp_node){
      pertu = 2.0 * (rnum - 0.5) * scal[pp] * 0.2;
    }else if(n.ntype[pp] == v_node)
      if((rvel == 0)&&(nvel>1)){
	pars[pp]= res[pp-n.v_off]*1.e5; pertu = 0.0;
      }else if((rvel == 1)&&(nvel>1)){
	pars[pp]= res[pp-n.v_off]*-1.e5; pertu = 0.0;
      }else
	pertu =  (rnum - 0.5) * scal[pp];
    else if(n.ntype[pp] == vturb_node){
      pertu =  (rnum - 0.75) * 2.e5;
    }else if(n.ntype[pp] == bl_node)
      pertu =  (rnum-0.5) * scal[pp];
    else if(n.ntype[pp] == bh_node)
      pertu =  (rnum-0.5) * scal[pp];
    else if(n.ntype[pp] == azi_node)
      pertu =  (rnum-0.5) * phyc::PI;   
    else if(n.ntype[pp] == pgas_node){
      pertu = 0.0;//1.0 + (rnum-0.35);
      //pars[pp] = 0.0;
    }
    if(input.depth_model == 0)
      pars[pp] = checkParameter(pars[pp] + pertu, pp);
    else
      pars[pp] = pars[pp] + pertu;
  }

  
}

/* --------------------------------------------------------------------------------------------------- */

int const_dregul(int n, double *ltau, double *var, double weight, reg_t &dregul, double m, int off, int roff, int typ)
{
  
  /* --- 
     Penalize deviations from me, removed the factor x2 from the derivative
     to be consistent with our LM definition 
     
     Penalty: pe = we * (xi - c)**2
     d/dxi pe  = 2*we*(xi-x)

     I am actually dividing everything by npar,
     so that the penalty remains relatively constant regardless
     of the number of nodes. It is a mean value when all the penalties
     are summed.
     --- */  

  double *reg = dregul.reg, **LL = dregul.LL, **dreg = dregul.dreg; 
  
  
  double c = sqrt(weight/n);
  double dtau = 1.0;// ttau = ((n>1)?fabs(ltau[n-1] - ltau[0]):1.0);
  double dtmax = 0.0;
  
  for(int ii=1; ii<n; ii++){
    double dt = fabs(ltau[ii]-ltau[ii-1]);
    if(dt > dtmax) dtmax = dt;
  }
  if(dtmax == 0.0) dtmax = 1.0;


  
  for(int yy=0;yy<n; yy++){

    
    if((yy == 0) && (n > 1)) dtau = fabs(ltau[0] - ltau[1])/dtmax;
    else if((yy > 0) && (yy < (n-1)) && (n > 1)) dtau = fabs(ltau[yy+1] - ltau[yy-1])*0.5/dtmax;
    else if((yy == (n-1)) && (n > 1)) dtau = fabs(ltau[yy] - ltau[yy-1])/dtmax;
    
    //dtau = 1.0;
    double tmp = c/ dtau;
    reg[roff+yy] = tmp*(var[yy] - m);
    dreg[roff+yy][off+yy] = tmp;
    dregul.rt[yy+roff] = typ;

    tmp *= tmp;
    LL[off+yy][off+yy] += c/dtau;
    
  }
  
  
  return n;
}

/* --------------------------------------------------------------------------------------------------- */

int mean_dregul(int n, double *ltau, double *var, double weight, reg_t &dregul, int off, int roff, int typ)
{
  
  /* --- 
     Penalize deviations from me, removed the factor x2 from the derivative
     to be consistent with our LM definition. Very similar to const_dregul, 
     but it is non-diagonal because all variables contribute to the mean
     and therefore we have all the block full with values.

     Penalty: pe = we**0.5 * (xi - sum xj/N)

     Derivative:
     d/dxi pe = 2 * we*(xi - sum xj/N) * (dij - 1/N)

     where dij is a Kronecker delta that takes value 1 for j==i terms 
     and zero for non-diagonal terms. Obviously (sum xj/N) is the mean value
     of the physical parameter as a function of height.

     --- */

  double *reg = dregul.reg, **LL = dregul.LL, **dreg = dregul.dreg; 

  
  if(n == 1){
    reg[roff] = 0.0;
    dreg[roff][off] = 0.0;
    return 1;
  }

  double dn = double(n), idn = 1.0/dn;
  double c = sqrt(weight/n);
  double me = mth::mean(n, var), dtau = 1.0;

  double dtmax = 0.0;
  for(int ii=1; ii<n; ii++){
    double dt = fabs(ltau[ii]-ltau[ii-1]);
    if(dt > dtmax) dtmax = dt;
  }
  if(dtmax == 0.0) dtmax = 1.0;

  
  for(int yy=0;yy<n; yy++){

    
    if((yy == 0) && (n > 1)) dtau = fabs(ltau[0] - ltau[1]) / dtmax;
    else if((yy > 0) && (yy < (n-1)) && (n > 1)) dtau = fabs((ltau[yy+1] - ltau[yy-1]) / dtmax) * 0.5;
    else if((yy == (n-1)) && (n > 1)) dtau = fabs(ltau[yy] - ltau[yy-1]) / dtmax;
    
    
    reg[roff+yy] = c*(var[yy] - me) / dtau;
    dregul.rt[yy+roff] = typ;

    for(int xx=0; xx<n; xx++){

      if(xx == yy) dreg[roff+yy][off+xx] = c * (1.0 - idn) / dtau;
      else         dreg[roff+yy][off+xx] = c * (    - idn) / dtau;
      
    }
  }

  
  for(int yy=0; yy<n; yy++){
    for(int xx=0; xx<n; xx++){
      double sum = 0.0;
      for(int kk=0;kk<n;kk++)
	sum += dreg[roff+kk][off+xx] * dreg[roff+kk][off+yy];
      LL[yy+roff][xx+off] += sum;
    }
  }
      
  


  
  return n;
}

/* --------------------------------------------------------------------------------------------------- */

int tikhonov1_dregul(int n, double *ltau, double *var, double weight, reg_t &dregul, int off, int roff, int typ)
{

  double *reg = dregul.reg, **LL = dregul.LL, **dreg = dregul.dreg; 
  double c = weight / (double)(n-1);
  double c_sqrt = sqrt(c);
  
  /* --- Derivative in the first and last points:

     Penalty: pe =  c * ((var[k+1] - var[k])/|dltau|)
     d/x_i pe     = c / |dltau|
     d/x_i-1 pe  =  c * (-1) / |dltau| = -d/dx_i pe

     In the LM we have divided the x2 in all terms so we do it here too in 
     the implementation.
     
     --- */

  
  //double ttau = fabs(ltau[n-1] - ltau[0]);
  double dtmax = 0.0;
  for(int ii=1; ii<n; ii++){
    double dt = fabs(ltau[ii]-ltau[ii-1]);
    if(dt > dtmax) dtmax = dt;
  }
  if(dtmax == 0.0) dtmax = 1.0;
  
  
  for(int yy = 1; yy<n; yy++){

    /* --- 2 terms around the diagonal of J --- */
    
    double dtau = (abs(ltau[yy] - ltau[yy-1]) / dtmax);
    dregul.rt[yy-1+roff] = typ;	
    reg[yy-1+roff] = c_sqrt * (var[yy] - var[yy-1]) / dtau;
    //
    double tmp = c_sqrt / dtau;
    dreg[yy-1+roff][yy+off] = tmp;
    dreg[yy-1+roff][yy-1+off] = -tmp;
    //
    tmp *= tmp;
    LL[yy+off-1][yy+off-1] += tmp;
    LL[yy+off  ][yy+off  ] += tmp;
    
    LL[yy+off][yy+off-1] += - tmp;
    LL[yy+off-1][yy+off]   += - tmp;
  }
  
  return n-1;
}

/* --------------------------------------------------------------------------------------------------- */


int scalGrad_dregul(int n, double *ltau, double *var, double weight, reg_t &dregul, int off, int roff, int typ)
{

  double *reg = dregul.reg, **LL = dregul.LL, **dreg = dregul.dreg; 
  double c = weight/ (double)(n-1);
  double c_sqrt = sqrt(c);

  
  /* --- Derivative in the first and last points:

     Penalty: pe =  c * ((var[k+1] - var[k])/|dltau|)
     d/x_i pe     = c / |dltau|
     d/x_i-1 pe  =  c * (-1) / |dltau| = -d/dx_i pe

     In the LM we have divided the x2 in all terms so we do it here too in 
     the implementation.
     
     --- */

  
  //double ttau = fabs(ltau[n-1] - ltau[0]);
  double dtmax = 0.0;
  for(int ii=1; ii<n; ii++){
    double dt = fabs(ltau[ii]-ltau[ii-1]);
    if(dt > dtmax) dtmax = dt;
  }
  if(dtmax == 0.0) dtmax = 1.0;
  
  
  for(int yy = 1; yy<n; yy++){

    /* --- 2 terms around the diagonal of J --- */

    double ya = var[yy-1], yb = var[yy];
    
    double dtau = (abs(ltau[yy] - ltau[yy-1]) / dtmax);
    double ysum = sqrt(1.0 + SQ(ya) + SQ(yb));
    double tmp = c_sqrt / dtau;

    dregul.rt[yy-1+roff] = typ;	
    reg[yy-1+roff] = tmp * 2.0 * (yb-ya) / ysum;
    //
    dreg[yy-1+roff][yy+off] = tmp * (2.0*(1+ya*ya + ya*yb)) / CUB(ysum);
    dreg[yy-1+roff][yy-1+off] = -tmp * (2.0*(1.0+ya*yb + yb*yb)) / CUB(ysum);
    //
    tmp *= tmp*0.5;
    ysum = CUB(1.0 + SQ(ya) + SQ(yb));


    LL[yy+off-1][yy+off-1] += tmp * (8.*(1.-3.*SQ(ya) - 2.*CUB(ya)*yb + SQ(yb) + 6*ya*(yb+CUB(yb)))) / ysum;
    LL[yy+off  ][yy+off  ] += tmp * (8.*(1.+SQ(ya) + 6.*ya*yb*(1.+SQ(ya)) - 3.*SQ(yb) - 2.*ya*CUB(yb))) / ysum;
    
    double tmp1 = 8.*(-1. + SQ(SQ(ya)) -4.*ya*yb - 6.*SQ(ya)*SQ(yb) + SQ(SQ(yb))) / ysum;
    LL[yy+off][yy+off-1] +=  tmp * tmp1 ;
    LL[yy+off-1][yy+off]   +=  tmp * tmp1;
  }
  
  return n-1;
}
/* --------------------------------------------------------------------------------------------------- */

int secDer_dregul(int n, double *ltau, double *var, double weight, reg_t &dregul, int off, int roff, int typ)
{

  double *reg = dregul.reg, **LL = dregul.LL, **dreg = dregul.dreg; 

  double c = weight / double(n-2);
  double c_sqrt = sqrt(c);
  
  /* --- non-equidistant second derivative:
     ypp = A*y[ii+1] + B*y[ii] + C*y[ii-1]

     with:
     A = 2  / (ddx*(ddx+udx))
     B = -2 / (ddx*udx)
     C = 2  / (udx*(ddx+udx))

     The derivatives are trivially A,B and C
     
     --- */

  double dtmax = 0.0;
  for(int ii=1; ii<n; ii++){
    double dt = fabs(ltau[ii]-ltau[ii-1]);
    if(dt > dtmax) dtmax = dt;
  }
  if(dtmax == 0.0) dtmax = 1.0;

  
  
  for(int yy = 1; yy<(n-1); yy++){

    /* --- 3 terms around the diagonal of J --- */
    
    double dt0 = fabs(ltau[yy] - ltau[yy-1]) / dtmax;
    double dt1 = fabs(ltau[yy+1] - ltau[yy]) / dtmax;
    double tmp = 2.0/(dt0+dt1);
    double A = tmp / dt1, C = tmp / dt0;
    double B = -2.0 / (dt0*dt1);
    
    dregul.rt[yy-1+roff] = typ;
    reg[yy-1+roff] = c_sqrt * (B*var[yy] + C*var[yy-1] + A*var[yy+1]);
    dreg[yy-1+roff][yy-1+off] += c_sqrt * C;
    dreg[yy-1+roff][yy+off]   += c_sqrt * B;
    dreg[yy-1+roff][yy+1+off] += c_sqrt * A;
 
    LL[yy-1+off][yy-1+off] += c * C*C;
    LL[yy+off  ][yy+off  ] += c * B*B;
    LL[yy+1+off][yy+1+off] += c * A*A;

    tmp = c*C*B;
    LL[yy-1+off][yy+0+off] += tmp; LL[yy+0+off][yy-1+off] += tmp;

    tmp = c*A*C;
    LL[yy-1+off][yy+1+off] += tmp; LL[yy+1+off][yy-1+off] += tmp;

    tmp = c*A*B;
    LL[yy+1+off][yy+0+off] += tmp; LL[yy+0+off][yy+1+off] += tmp;
  }
  
  return n-2;
}



/* --------------------------------------------------------------------------------------------------- */

void getDregul2(double *m, int npar, reg_t &dregul, nodes_t &n)
{
  //static const double weights[7] = {0.006, 0.06, 0.1, 0.1, 0.1, 0.1,20.0};
  //  const double weights[7] = {0.0002, 0.6, 0.6, 0.1, 0.1, 0.1,20.0};
  double *weights = n.rewe;
  double  *ltau = NULL, we = 0.0;
  nodes_type_t ntype = none_node;
  int off = 0, roff = 0;

  
  dregul.zero();
  
  /* --- Penalize temp ? Only allow Tikhonov, the rest don't make sense for temp --- */

  if(n.toinv[0] && (n.regul_type[0] > 0)){ // temperature
    ntype = temp_node;
    off = (int)n.temp_off;
    ltau =  &n.temp[0];
    we = weights[0]*dregul.scl;
    
    int nn = (int)n.temp.size();
    switch(n.regul_type[0]){
    case(1):
      if((nn) >= 2){
	roff += tikhonov1_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff,0);
      }
      break;
    case(2):
      if(nn >= 2)
	roff += scalGrad_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff,0);
      break;
    case(4):
      if(nn >= 3)
	roff += secDer_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff, 0);
      break;
    case(5):
      if(nn >= 3){
	//	roff += tikhonov1_dregul(nn-1, &ltau[1], &m[off+1], we*0.6, dregul.dreg, dregul.reg, off+1, roff);
	//	roff +=    line3dregul(nn, ltau, &m[off], we, dregul, off, roff);
	roff+= secDer_dregul(nn, &ltau[0], &m[off],    we*0.3, dregul , off, roff,0);
	roff+= tikhonov1_dregul(nn, &ltau[0], &m[off], we*0.7, dregul , off, roff,0);
      }
      break;
    default:
      break;
    }
  }

  /* --- Penalize vlos? --- */
  
  if(n.toinv[1] && (n.regul_type[1] > 0)){ // vlos
    int nn = (int)n.v.size();
    ntype = v_node;
    off = (int)n.v_off;
    ltau =  &n.v[0];
    we = weights[1]*dregul.scl;
    
    switch(n.regul_type[1]){
    case(1):
      if(nn >= 2)
	roff += tikhonov1_dregul(nn, ltau, &m[off], we, dregul , off, roff,1);
      break;
    case(2):
      roff += mean_dregul(nn, ltau, &m[off], we, dregul , off, roff,1);
      break;
    case(3):
      roff += const_dregul(nn, ltau, &m[off], we, dregul, 0.0, off, roff,1);
      break;
    case(4):
      if(nn >= 3)
	roff += secDer_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff, 0);
      break;
    default:
      break;
    }
  }

  
   /* --- Penalize vturb? --- */
  
  if(n.toinv[2] && (n.regul_type[2] > 0)){ // vturb
    int nn = (int)n.vturb.size();
    ntype = vturb_node;
    off = (int)n.vturb_off;
    ltau =  &n.vturb[0];
    we = weights[2]*dregul.scl;
    switch(n.regul_type[2]){
    case(1):
      if(nn >= 2)
	roff += tikhonov1_dregul(nn, ltau, &m[off], we, dregul, off, roff,2);
      break;
    case(2):
      if(nn >= 2)
	roff += mean_dregul(nn, ltau, &m[off], we, dregul, off, roff,2);
      break;
    case(3):
      roff += const_dregul(nn, ltau, &m[off], we, dregul, 0.0, off, roff,2);
      break;
    case(4):
      if(nn >= 3)
	roff += secDer_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff, 0);
      break;
    case(5):
      if(nn >= 3){
	roff += secDer_dregul(nn, &ltau[0], &m[off], we*0.4, dregul , off, roff, 0);
	roff += const_dregul(nn, ltau, &m[off], we*0.6, dregul, 0.0, off, roff,2);
      }
      break;
    default:
      break;
    }
  }

  
  /* --- Penalize Blong? --- */
  
  if(n.toinv[3] && (n.regul_type[3] > 0)){ // B
    int nn = (int)n.bl.size();
    ntype = bl_node;
    off = (int)n.bl_off;
    ltau =  &n.bl[0];
    we = weights[3]*dregul.scl;

    switch(n.regul_type[3]){
    case(1):
      if(nn >= 2)
	roff += tikhonov1_dregul(nn, ltau, &m[off], we, dregul, off, roff,3);
      break;
    case(2):
      roff += mean_dregul(nn, ltau, &m[off], we, dregul, off, roff,3);
      break;
    case(3):
      roff += const_dregul(nn, ltau, &m[off], we, dregul, 0.0, off, roff,3);
      break;
    case(4):
      if(nn >= 3)
	roff += secDer_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff, 0);
      break;
    default:
      break;
    }
  }

  /* --- Penalize Bhor? --- */
  
  if(n.toinv[4] && (n.regul_type[4] > 0)){ // inc
    int nn = (int)n.bh.size();
    ntype = bh_node;
    off = (int)n.bh_off;
    ltau =  &n.bh[0];
    we = weights[4]*dregul.scl;

    switch(n.regul_type[4]){
    case(1):
      if(nn >= 2)
	roff += tikhonov1_dregul(nn, ltau, &m[off], we, dregul, off, roff,4);
      break;
    case(2):
      roff += mean_dregul(nn, ltau, &m[off], we, dregul, off, roff,4);
      break;
    case(3):
      roff += const_dregul(nn, ltau, &m[off], we, dregul, 0.0, off, roff,4);
      break;
    case(4):
      if(nn >= 3)
	roff += secDer_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff, 0);
      break;
    default:
      break;
    }
  }

  /* --- Penalize azi? --- */
  
  if(n.toinv[5] && (n.regul_type[5] > 0)){ // azi
    int nn = (int)n.azi.size();
    ntype = azi_node;
    off = (int)n.azi_off;
    ltau =  &n.azi[0];
    we = weights[5]*dregul.scl;

    switch(n.regul_type[5]){
    case(1):
      if(nn >= 2)
	roff += tikhonov1_dregul(nn, ltau, &m[off], we, dregul, off, roff,5);
      break;
    case(2):
      roff += mean_dregul(nn, ltau, &m[off], we, dregul, off, roff,5);
      break;
    case(4):
      if(nn >= 3)
	roff += secDer_dregul(nn, &ltau[0], &m[off], we, dregul , off, roff, 0);
      break;
    default:
      break;
    }
  }
  

  
  /* --- Penalize pgas enhancement factor ? --- */
  
  if(n.toinv[6] && (n.regul_type[6] > 0)){ // pgas
    
    
    /* --- Penalize deviations from 1.0 --- */
    int nn = 1; double dum = 0.1;
    //ltau =  &n.[0];
    we =  weights[6]*dregul.scl;
    off = (int)n.pgas_off;
    roff += const_dregul(nn, &dum, &m[off], we, dregul, 1.0, off, roff, 6);
  }

  if(0){
    for(int ii=0; ii<dregul.nreg; ii++){
      for(int jj=0; jj<dregul.npar; jj++) fprintf(stderr,"%e ",dregul.dreg[ii][jj]);
      cerr<<endl;
    }
    cerr<<endl;
    
    for(int ii=0; ii<dregul.nreg; ii++) fprintf(stderr,"%e ",dregul.reg[ii]);
    cerr<<endl;
  }

  
  // for(int ii=0; ii<dregul.nreg; ii++) fprintf(stderr,"%e %d\n", dregul.reg[ii], dregul.nreg);
					  
  //fprintf(stderr,"%e\n", dregul.getReg() );
  //memset(&dregul.reg[0], 0,dregul.nreg*sizeof(double)); 
  
}

/* --------------------------------------------------------------------------------------------------- */

int getChi2(int npar1, int nd, double *pars1, double *syn_in, double *dev, double **derivs, void *tmp1,  reg_t &dregul, bool store){

  
  /* --- Cast tmp1 into a double --- */
  
  atmos &atm = *((atmos*)tmp1); 
  double *ipars = new double [npar1]();
  mdepth &m = *atm.imodel;
  double nd1 = 1.0 / sqrt((double)nd);

  
    /* --- Expand atmosphere ---*/
  
  for(int pp = 0; pp<npar1; pp++)
    ipars[pp] = pars1[pp] * atm.scal[pp];

  
  if(store){
    mdepth &m1 = *atm.imodel->ref_m;
    m1.expand(atm.input.nodes, &ipars[0], atm.input.dint, atm.input.depth_model);
    atm.checkBounds(m1);
    m1.getPressureScale(atm.input.nodes.depth_t, atm.input.boundary, atm.eos);
    delete [] ipars;
    
    return 0;
  }
  
  

  
  m.expand(atm.input.nodes, &ipars[0], atm.input.dint, atm.input.depth_model);
  atm.checkBounds(m);
  m.getPressureScale(atm.input.nodes.depth_t, atm.input.boundary, atm.eos);
  
  
  
  /* --- Compute synthetic spetra --- */
  
  memset(&atm.isyn[0], 0, nd*sizeof(double));
  //  for(int ii=0; ii<m.ndep;ii++) fprintf(stderr,"%e %e %e %e %e\n", m.cmass[ii], m.temp[ii], m.v[ii], m.vturb[ii], m.pgas[ii]);
  bool conv = atm.synth( m , &atm.isyn[0], 0, (cprof_solver)atm.input.solver, true);  
  
  
  if(!conv){
    atm.cleanup();
    return 1;
  }
  
 /* --- Compute derivatives? --- */
  
  if(derivs){
    
    for(int pp = npar1-1; pp >= 0; pp--){
      
      if(derivs[pp]){
	//atm.cleanup();
	/* --- Compute response function ---*/
	
	memset(&derivs[pp][0], 0, nd*sizeof(double));
	atm.responseFunction(npar1, m, &ipars[0], nd,
			      &derivs[pp][0], pp, &atm.isyn[0]);

	
	/* --- Degrade response function --- */
	
	atm.spectralDegrade(atm.input.ns, (int)1, nd, &derivs[pp][0]);

	
	/* --- renormalize the response function by the 
	   scaling factor and divide by the noise --- */
		
	for(int ii = 0; ii<nd; ii++)
	  derivs[pp][ii] *= (atm.scal[pp] / atm.w[ii]) * nd1;
	  //derivs[pp][ii]  /= atm.w[ii];
	
      }
    }    
  }


  /* --- Degrade synthetic spectra --- */
  
  atm.spectralDegrade(atm.input.ns, (int)1, nd, &atm.isyn[0]);
  memcpy(syn_in, &atm.isyn[0], nd*sizeof(double));


  /* --- Compute residue --- */

  for(int ww = 0; ww < nd; ww++){
    dev[ww] = (atm.obs[ww] - atm.isyn[ww]) / atm.w[ww] * nd1;
  }


  /* ---  compute regularization --- */ 

  if(dregul.to_reg)
    getDregul2(pars1, npar1, dregul, atm.input.nodes);
      
  /* --- clean up --- */
  
  delete [] ipars;

  return 0;
}

reg_t init_dregul(int npar, nodes_t &n, double scl, double scl1, int ntrans)
{
  
  /* --- detect number of penalty functions --- */

  int npen = 0, nn = 0;
  if(n.toinv[0] && (n.regul_type[0] > 0)){
    nn = n.temp.size();
    if     ((n.regul_type[0] == 1) && (nn >= 2)) npen += nn-1;
    else if((n.regul_type[0] == 2) && (nn >= 2)) npen += nn-1;
    else if((n.regul_type[0] == 4) && (nn >= 3)) npen += nn-2;
    else if((n.regul_type[0] == 5) && (nn >= 3)) npen += nn-2 + nn-1;
  }

  if(n.toinv[1] && (n.regul_type[1] > 0)){
    nn = n.v.size();
    if((n.regul_type[1] == 1) && (nn >= 2)) npen += nn-1;
    else if(n.regul_type[1] == 2)           npen += nn;
    else if((n.regul_type[1] == 4) && (nn >= 3)) npen += nn-2;
  }


  if(n.toinv[2] && (n.regul_type[2] > 0)){
    nn = n.vturb.size();
    if((n.regul_type[2] == 1) && (nn >= 2)) npen += nn-1;
    else if((n.regul_type[2] == 2) && (nn >=2)) npen += nn;
    else if((n.regul_type[2] == 3))             npen += nn;
    else if((n.regul_type[2] == 4) && (nn >= 3)) npen += nn-2;
    else if((n.regul_type[2] == 5) && (nn >= 3)) npen += nn-2 + nn;
    else                     npen += nn;
  }

  if(n.toinv[3] && (n.regul_type[3] > 0)){
    nn = n.bl.size();
    if     ((n.regul_type[3] == 1) && (nn >= 2)) npen += nn-1;
    else if((n.regul_type[3] == 4) && (nn >= 3)) npen += nn-2;
    else                                         npen += nn;
  }

  if(n.toinv[4] && (n.regul_type[4] > 0)){
    nn = n.bh.size();
    if((n.regul_type[4] == 1) && (nn >= 2)) npen += nn-1;
    else if((n.regul_type[4] == 4) && (nn >= 3)) npen += nn-2;
    else                     npen += nn;
  }
  
  if(n.toinv[5] && (n.regul_type[5] > 0)){
    nn = n.azi.size();
    if((n.regul_type[5] == 1) && (nn >= 2)) npen += nn-1;
    else if((n.regul_type[5] == 4) && (nn >= 3)) npen += nn-2;
    else if((n.regul_type[5] == 2)) npen += nn;
  }

  if(n.toinv[6] && (n.regul_type[6] > 0)){
    npen += 1;
  }
  
  reg_t res;
  if(npen > 0)
    res.set(npar, npen, scl, scl1, ntrans);
  return res;
}


double atmos::fitModel2(mdepth_t &m, int npar, double *pars, int nobs, double *o, mat<double> &weights){


    
  /* --- compute tau scale ---*/
  
  for(int k = 0; k < (int)m.ndep; k++) m.tau[k] = pow(10.0, m.ltau[k]);
  
  
  /* --- point to obs and weights --- */
  mdepth_t m1 = m, m2 = m, best_m = m;
  bool depth_per = ((this->input.depth_model > 0) ? true : false);

  obs = &o[0];
  w = &weights.d[0];
  imodel = ((depth_per)? &m1 : &m);

  if(depth_per){
    imodel->ref_m = &m;
    m.ref_m = &m;
    checkBounds(m);
    checkBounds(m1);
    memset(pars, 0, npar*sizeof(double));
  }else          imodel->ref_m = NULL;


  /* --- get regularization --- */
  reg_t regul;
  if(input.nodes.regularize[0] >= 1.e-8)
    regul = init_dregul(npar, input.nodes, input.nodes.regularize[0], input.nodes.regularize[1], input.nodes.nregul);

  
  /* --- Invert pixel randomizing parameters if iter > 0 --- */
  
  int ndata = input.nw_tot * input.ns;
  double ipars[npar], bestPars[npar], ichi, bestChi = 1.e10;
  isyn.resize(ndata);
  vector<double> bestSyn;
  bestSyn.resize(ndata);

  
  /* --- Init clm --- */

  clm lm = clm(ndata, npar);
  lm.xtol = 1.e-2;
  lm.verb = input.verbose;
  lm.use_geo_accel = input.use_geo_accel;
  
  if(input.marquardt_damping > 0.0) lm.ilambda = input.marquardt_damping;
  else                              lm.ilambda = 1.0;
  lm.maxreject = 8;
  lm.svd_thres = max(input.svd_thres, 1.e-16);
  lm.chi2_thres = input.chi2_thres;
  lm.lmax = 1.e5;
  lm.lmin = 1.e-4;
  lm.lfac = sqrt(10.);
  lm.proc = input.myrank;
  lm.delay_bracket = input.delay_bracket;
  
  if(input.nodes.regularize[0] >= 1.e-5){
    lm.regularize = true;
    lm.regul_scal_in = input.nodes.regularize[0];
  } else lm.regularize = false;
  
  /* ---  Set parameter limits --- */
  
  for(int pp = 0; pp<npar; pp++){

    if(!depth_per){
      lm.fcnt[pp].limit[0] = mmin[pp]/scal[pp];
      lm.fcnt[pp].limit[1] = mmax[pp]/scal[pp];
      lm.reset_par = false;
    }else{
      lm.fcnt[pp].limit[0] = -maxc[pp]/scal[pp];
      lm.fcnt[pp].limit[1] =  maxc[pp]/scal[pp];
      lm.reset_par = true;
    }
    
    lm.fcnt[pp].scl = 1.0;//scal[pp];
    lm.ptype[pp] = (input.svd_split) ? (unsigned)input.nodes.ntype[pp] : (unsigned)0;
    if(input.nodes.ntype[pp] ==  pgas_node && input.svd_split > 0)  lm.ptype[pp] = (unsigned)temp_node; // Treat Pgas as a temperature variable

      
    if(input.nodes.ntype[pp] == azi_node) lm.fcnt[pp].cyclic = true;
    else                                  lm.fcnt[pp].cyclic = false;
    
    lm.fcnt[pp].bouncing = false;
    lm.fcnt[pp].capped = 1;
    
    if(input.nodes.ntype[pp] == temp_node){
      lm.fcnt[pp].relchange = ((depth_per) ? false : true);
      lm.fcnt[pp].maxchange[0] = 0.5;
      lm.fcnt[pp].maxchange[1] = 2.0;
    }else{
      lm.fcnt[pp].relchange = false;
      lm.fcnt[pp].maxchange[0] = maxc[pp]/scal[pp];
      lm.fcnt[pp].maxchange[1] = maxc[pp]/scal[pp];
    }
    
    if(!depth_per)
      pars[pp] = checkParameter(pars[pp], pp);
  }
  

  
  /* --- Loop iters --- */
  int do_vel_grad = -1;
  
  for(int iter = 0; iter < input.nInv; iter++){

    cleanup();

    
    
    /* --- init parameters for this inversion --- */

    if(!depth_per)
      memcpy(&ipars[0], &pars[0], npar * sizeof(double));
    else
      memset(&ipars[0], 0, npar * sizeof(double));

    
    if(iter > 0 || input.random_first){
      if(input.vgrad) do_vel_grad+= 1;
      randomizeParameters(input.nodes , npar, &ipars[0], do_vel_grad);
      if(depth_per){
	imodel->ref_m->expand(input.nodes, &ipars[0], input.dint, input.depth_model);
	checkBounds(*imodel->ref_m);
	imodel->ref_m->getPressureScale(input.nodes.depth_t, input.boundary, eos);
      }
    }
    
    /* --- Work with normalized parameters --- */
    
    for(int pp = 0; pp<npar; pp++){
      ipars[pp] /= scal[pp];      
    }
    
    
    /* --- Call clm --- */

    double chi2 = lm.fitdata(getChi2, &ipars[0], (void*)this, input.max_inv_iter, regul);

    
    
    /* --- Re-start populations --- */
    
    cleanup();
    

    /* --- Store result? ---*/
    
    if(chi2 < bestChi){
      bestChi = chi2;
      memcpy(&bestPars[0], &ipars[0], npar*sizeof(double));
      memcpy(&best_m.cub.d[0], &m.cub.d[0], m.ndep*14*sizeof(double));
      memcpy(&bestSyn[0], &lm.bestSyn[0], ndata*sizeof(double));
    }

    if(depth_per){
      memcpy(&m.cub.d[0], &m2.cub.d[0], m.ndep*14*sizeof(double));
      memcpy(&m1.cub.d[0], &m2.cub.d[0], m.ndep*14*sizeof(double));
    }
    
    /* --- reached thresthold? ---*/
    
    if(bestChi <= input.chi2_thres) break;
    
  } // iter

  
  double sum = 0.0;
  for(int ww = 0; ww<ndata;ww++){
    double tmp = (obs[ww] - bestSyn[ww]) / weights.d[ww];
    sum += (tmp*tmp);
  }
  //fprintf(stderr,"Recomp chi2=%13.5f\n", sum/ndata);
  memcpy(&obs[0], &bestSyn[0], ndata*sizeof(double));
  memcpy(&m.cub.d[0], &best_m.cub.d[0], m.ndep*14*sizeof(double));
  
  /* --- Clean-up --- */
  
  isyn.clear();
  bestSyn.clear();

  return bestChi;
}


void atmos::spectralDegrade(int ns, int npix, int ndata, double *obs){

  /* --- loop and degrade --- */
  
  if(inst == NULL) return;

  int nreg = (int)input.regions.size();
  for(int ipix = 0; ipix<npix; ipix++){
    for(int ireg = 0; ireg<nreg; ireg++){
      inst[ireg]->degrade(&obs[ipix * ndata], ns);
    }
  }

  
}


//-------------------------------------------------------------------------- //


void atmos::responseFunctionFull(mdepth_t m, int nd, double *out_in, double *syn, int pp)
{
  if(pp >= 6){
    pp++;
  }
  
  double dpar = ((input.dpar >= 1.e-3)?input.dpar : 0.01);
  double pertu = scal[pp] * dpar;
  int ndep = (int)m.ndep, centder = input.centder;
  bool store_pops = false;
  //double (&out)[ndep][nd] = *reinterpret_cast< double (*)[ndep][nd]>(out_in);
  double **out = mmem::var2dim(out_in, ndep, nd);
  

  if(centder){
    vector<double> syup(nd, 0.0), sydow(nd,0.0);

    for(int kk=0;kk<ndep;kk++){
      double  up=0, down=0;
      double pval = m.cub(pp, kk);

      if(pp<6){
	up = pertu, down = -pertu;
      }else{
	up = pval * dpar, down = -up;
      }
	
      /* --- upper perturbation --- */
      
      if((pval + pertu) > mmax[kk]){
	up = 0.0;
	memcpy(&syup[0],syn, nd*sizeof(double));
      }else{
	m.cub(pp,kk) = pval + up;
	/*
	if(pp == 7){
	  eos.nne_from_T_rho_nne(m.temp[kk], m.pgas[kk],  m.rho[kk], m.nne[kk]);
	  eos.store_partial_pressures(m.ndep, kk, eos.xna, m.nne[kk]);
	}
	*/
	synth(m, &syup[0], 1, (cprof_solver)input.solver, store_pops);
	
      }


      /* --- Lower perturbation --- */

      if((pval - pertu) < mmin[kk]){
	down = 0.0;
	memcpy(&sydow[0],syn, nd*sizeof(double));
      }else{
	m.cub(pp,kk) = pval + down;

	/*
	if(pp == 7){
	  eos.nne_from_T_rho_nne(m.temp[kk], m.pgas[kk],  m.rho[kk], m.nne[kk]);
	  eos.store_partial_pressures(m.ndep, kk, eos.xna, m.nne[kk]);
	}
	*/
	synth(m, &sydow[0], 1, (cprof_solver)input.solver, store_pops);
      }

      
      // --- restore value --- //
      
      m.cub(pp,kk) = pval;
      
      if(pp == 7){
	eos.nne_from_T_rho_nne(m.temp[kk], m.pgas[kk],  m.rho[kk], m.nne[kk]);
	eos.store_partial_pressures(m.ndep, kk, eos.xna, m.nne[kk]);
      }
      
      // --- compute response ---//

      up += down;
      if(up > 1.e-6) up = 1.0/up;
      else           up = 0.0;
      
      for(int ww=0;ww<nd;ww++) out[kk][ww] = (syup[ww]-sydow[ww]) * up;
    }// kk
    
    syup.clear();
    sydow.clear();
    
  }else{ // side derivative
    
    for(int kk=0;kk<ndep;kk++){
      double pval = m.cub(pp, kk);
      double per=0;
      
      if(pp < 6){
	per = (((pval+pertu) < mmax[pp])? pertu : -pertu);
      }else{
	per = (((pval*(1.0+dpar)) < mmax[pp])? dpar*pval : -dpar*pval);
      }
      m.cub(pp,kk) += per;

      /*
      if(pp == 7){
	eos.nne_from_T_rho_nne(m.temp[kk], m.pgas[kk],  m.rho[kk], m.nne[kk]);
	eos.store_partial_pressures(m.ndep, kk, eos.xna, m.nne[kk]);
      }
      */
	
      synth(m, &out[kk][0], 1, (cprof_solver)input.solver, store_pops);

      m.cub(pp, kk) = pval;

      /*
      if(pp == 7){
	eos.nne_from_T_rho_nne(m.temp[kk], m.pgas[kk],  m.rho[kk], m.nne[kk]);
	eos.store_partial_pressures(m.ndep, kk, eos.xna, m.nne[kk]);
      }
      */
      
      per = 1.0/per;
      for(int ww=0;ww<nd;ww++) out[kk][ww] = (out[kk][ww]-syn[ww]) * per;
    } // kk

  }//else


  delete [] out;
  
}



