#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "crh.h"
#include "cmemt.h"
#include "input.h"
#include "physical_consts.h"
#include "interpol.h"
extern "C" {
#include "rhf1d.h"
#include "rh.h"
}
//
using namespace std;
using namespace phyc;
//
const double crh::pmax[7]  = {80000., 50.e5, 12e5, 5000.0, 5000., PI, 15.0};
const double crh::pmin[7]  = {2700. ,-50.e5,  +0.0,-5000.0,  +0.0,  +0.0, 0.1};
const double crh::pscal[7] = {1200. , 6.0e5, 6.0e5, 1000.0, 1000.0, PI, 10.0};
const double crh::pstep[7] = {1.e-1 , 1.e-1, 1.0e-1, 2.0e-1, 2.0e-1, 1.0e-1, 1.0e-1};

/* ----------------------------------------------------------------*/

vector<double> crh::get_max_limits(nodes_t &n){

  int nnodes = (int)n.nnodes;
  mmax.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) mmax[k] = pmax[0];
    else if(n.ntype[k] == v_node    ) mmax[k] = pmax[1];
    else if(n.ntype[k] == vturb_node) mmax[k] = pmax[2];
    else if(n.ntype[k] == bl_node    ) mmax[k] = pmax[3];
    else if(n.ntype[k] == bh_node  ) mmax[k] = pmax[4];
    else if(n.ntype[k] == azi_node  ) mmax[k] = pmax[5];
    else if(n.ntype[k] == pgas_node ) mmax[k] = pmax[6];
    else                              mmax[k] = 0;
  }
  return mmax;
}

/* ----------------------------------------------------------------*/


vector<double> crh::get_min_limits(nodes_t &n){

  int nnodes = (int)n.nnodes;
  mmin.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) mmin[k] = pmin[0];
    else if(n.ntype[k] == v_node    ) mmin[k] = pmin[1];
    else if(n.ntype[k] == vturb_node) mmin[k] = pmin[2];
    else if(n.ntype[k] == bl_node    ) mmin[k] = pmin[3];
    else if(n.ntype[k] == bh_node  ) mmin[k] = pmin[4];
    else if(n.ntype[k] == azi_node  ) mmin[k] = pmin[5];
    else if(n.ntype[k] == pgas_node ) mmin[k] = pmin[6];
    else                              mmin[k] = 0;
  }
  return mmin;
}

/* ----------------------------------------------------------------*/

vector<double> crh::get_scaling(nodes_t &n){

  int nnodes = (int)n.nnodes;
  scal.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) scal[k] = pscal[0];
    else if(n.ntype[k] == v_node    ) scal[k] = pscal[1];
    else if(n.ntype[k] == vturb_node) scal[k] = pscal[2];
    else if(n.ntype[k] == bl_node    ) scal[k] = pscal[3];
    else if(n.ntype[k] == bh_node  ) scal[k] = pscal[4];
    else if(n.ntype[k] == azi_node  ) scal[k] = pscal[5];
    else if(n.ntype[k] == pgas_node ) scal[k] = pscal[6];
    else                              scal[k] = 1.0;
  }
  return scal;
}

/* ----------------------------------------------------------------*/

vector<double> crh::get_steps(nodes_t &n){

  int nnodes = (int)n.nnodes;
  step.resize(nnodes);

  double ipstep[7];
  double sum = 0.0;
  
  for(int ii=0;ii<7;ii++){
    sum += pstep[ii];
    ipstep[ii] = pstep[ii];
  }
  sum /= 7.0;
  
  for(int ii=0;ii<7;ii++){
    ipstep[ii]/= sum;
    // cerr<<ipstep[ii]<<endl;
  }
  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) step[k] = ipstep[0];
    else if(n.ntype[k] == v_node    ) step[k] = ipstep[1];
    else if(n.ntype[k] == vturb_node) step[k] = ipstep[2];
    else if(n.ntype[k] == bl_node    ) step[k] = ipstep[3];
    else if(n.ntype[k] == bh_node  ) step[k] = ipstep[4];
    else if(n.ntype[k] == azi_node  ) step[k] = ipstep[5];
    else if(n.ntype[k] == pgas_node ) step[k] = ipstep[6];
    else                              step[k] = 1.0;
  }
  return step;
}

/* ----------------------------------------------------------------*/

crh::crh(iput_t &inpt, double grav): atmos(inpt, grav){
  /*
  input.lines = inpt.lines;
  input.regions = inpt.regions;
  input.solver = inpt.solver;
  input.mu = inpt.mu;
  */
  input = inpt;
  
  /* --- Copy lines array and --- */
  nlines = input.regions.size();
  
 
  /* --- Copy wavelength array and check which lines are computed at each wav --- */
  nlambda = 0;
  for(auto &it: input.regions){
    
    it.off = nlambda;
    nlambda += it.nw;	
    it.wav.resize(it.nw);
    it.nu.resize(it.nw);
    for(int ii = 0; ii<it.nw; ii++){
      // Compute wavelength array in vacuum
      it.wav[ii] = inv_convl( convl(it.w0) + it.dw * double(ii)); // lambda
      it.nu[ii] = CC  / (it.wav[ii] * 1.e-8);  // nu (Freq. in s^-1)
    }

  }

  /* --- Create array of all lambdas --- */
  lambda.resize(nlambda);
  int kk = 0;
  for(auto &it: input.regions){
    for(int ii = 0; ii<it.nw; ii++){
      lambda[kk++] = it.wav[ii]*0.1;
    }
  }


  /* --- add lambda = 5000 A at the end of the array for reference --- */
  // lambda.push_back(inv_convl(5000.0));

  /* --- Init saved pop --- */
  save_pop.pop = NULL;
  save_pop.nactive = 0;


  /* --- Init limits for inversion if nodes are present --- */

  if(input.nodes.nnodes > 0){
    vector<double> dummy;
    dummy = this->get_scaling(input.nodes);
    dummy = this->get_max_limits(input.nodes);
    dummy = this->get_min_limits(input.nodes);
    dummy = this->get_max_change(input.nodes);
  }
  
  

  /* --- Set-up output file --- */

  // sprintf()
  
  
}

/* ----------------------------------------------------------------*/

bool crh::synth(mdepth_t &m_in, double *syn, cprof_solver sol, bool save_pops){

  static int ncall = 0;
  ncall++;

  /* --- Copy model, RH seems to tamper with the model --- */
  
  mdepth m(m_in.ndep);
  m.cub.d = m_in.cub.d;
  
  /* --- Init vectors --- */
  
  float xa=0.0, xe=0.0;
  vector<float> frac, part;
  frac.resize(m.ndep,0.0);
  part.resize(m.ndep,0.0);
  nhtot.resize(m.ndep,0.0);
  double *B = new double [m.ndep], *inc = new double [m.ndep];

  
  
  /* --- Init storage for RH spectra --- */
  
  ospec sp = {};
  sp.I = NULL;
  sp.Q = NULL;
  sp.U = NULL;
  sp.V = NULL;
  sp.lambda = NULL;

  
  /* --- Restore nHtot and convert model to SI units --- */
  
  for(int kk = 0; kk < m.ndep; kk++){


    /* --- 
       nHtot 
       --- */
    
    //eos.read_partial_pressures(kk, frac, part, xa, xe);
    
    double xna = m.pgas[kk] / (phyc::BK * m.temp[kk]) - m.nne[kk];
    nhtot[kk] = xna * eos.ABUND[0] / eos.totalAbund * 1.e6; // nHtot based on abundance
    //nhtot[kk] = m.rho[kk] / (phyc::AMU * eos.avmol) *  eos.ABUND[0] / eos.totalAbund * 1.e6;
    
    /* --- Convert units to SI --- */

    if(fabs(m.cmass[kk]) >= 1.e-20)
      m.cmass[kk] = pow(10., m.cmass[kk]); // CMass should be in log10 scale
    m.cmass[kk] *= 10.0; // *= G_TO_KG / CM_TO_M**2
    m.rho[kk] *= 1000.; // G_TO_KG / CM_TO_M**3
    m.v[kk] *= -1.e-2; // CM_TO_M and change the sign so upflows are positive
    m.vturb[kk] *= 1.e-2; // CM_TO_M
    m.z[kk] *= 1.0e-2;    // CM_TO_M
    m.nne[kk] *= 1.e6;    // 1 / CM_TO_M**3
    m.tau[kk] = pow(10.0, m.ltau[kk]);

    B[kk] = sqrt(m.bl[kk]*m.bl[kk] + m.bh[kk]*m.bh[kk]) * 1.0e-4 ; // B in tesla
    if(B[kk] > 0.0) 
      inc[kk] = acos(m.bl[kk] * 1.0e-4 / B[kk]);
    else inc[kk] = 0.0;
    if(std::isnan(inc[kk])) inc[kk] = 0.0;
    
  }

  int savep = 0, hydrostat = 0;
  if(save_pops) savep = 1;
  
  
  /* --- Call RH --- */
  
  bool conv = rhf1d(input.mu, m.ndep, &m.temp[0], &m.rho[0], &m.nne[0], &m.vturb[0], &m.v[0],
		     &B[0], &inc[0], &m.azi[0], &m.z[0], &nhtot[0], &m.tau[0],
		    &m.cmass[0], 4.44, (bool_t)true, &sp, &save_pop, nlambda, &lambda[0],
		    input.myrank, savep, (int)input.verbose, &hydrostat);
  
  delete [] B;
  delete [] inc;

  /* --- Did not converge? --- */

  if(!conv){

    for(int ww = 0; ww<nlambda*4; ww++)
      syn[ww] = 1.e13;

    if(sp.I != NULL) delete [] sp.I;
    if(sp.Q != NULL) delete [] sp.Q;
    if(sp.U != NULL) delete [] sp.U;
    if(sp.V != NULL) delete [] sp.V;
    if(sp.lambda != NULL) delete [] sp.lambda;
    
    return (bool)false;
  }

  
  /* --- Retrieve spectra at the observed grid --- */
  
  for(auto &reg: input.regions){

    /* --- Get the indexes where the spectra are stored if firsttime --- */
    
    if((int)reg.idx.size() != sp.nlambda) lambdaIDX(sp.nlambda, sp.lambda);

    
    /* --- copy to output array and convert to CGS units --- */
    
    double scl = 1.0e3 / reg.cscal;
    
    int ele = reg.off*4;
    for(int ww=0;ww<reg.nw;ww++){
      syn[ele++] = sp.I[reg.idx[ww]] * scl;
      syn[ele++] = sp.Q[reg.idx[ww]] * scl;
      syn[ele++] = sp.U[reg.idx[ww]] * scl;
      syn[ele++] = sp.V[reg.idx[ww]] * scl;
    }
  }


  
  /* --- convert nHtot to Pgas using the electron density, the H abundance and temperature --- */
  
  if(hydrostat > 0){
    for(int kk = 0; kk < m.ndep; kk++){
      
      m.pgas[kk] = (nhtot[kk] * eos.totalAbund / eos.ABUND[0] + m.nne[kk]) *
	phyc::BK * m_in.temp[kk] * 1.e-6;
      
      if(kk > 0) m_in.pgas[kk] = m.pgas[kk]; // preserve pgas at the boundary otherwise the inversion can be unstable
      m_in.nne[kk] = m.nne[kk] * 1.0e-6;
    }
  }

  /* --- Deallocate sp --- */
  
  if(sp.I != NULL) delete [] sp.I;
  if(sp.Q != NULL) delete [] sp.Q;
  if(sp.U != NULL) delete [] sp.U;
  if(sp.V != NULL) delete [] sp.V;
  if(sp.lambda != NULL) delete [] sp.lambda;

  return conv;
}

/* ----------------------------------------------------------------*/

void crh::cleanup(void){
  clean_saved_populations(&save_pop);
}

/* ----------------------------------------------------------------*/

crh::~crh(void){
  cleanup();
}

/* ----------------------------------------------------------------*/

void crh::lambdaIDX(int nw, double *lamb){


  /* --- loop regions --- */
  for(auto &it: input.regions){

    it.idx.resize(it.nw,0.0);
    for(int ww=0; ww<it.nw;ww++){
      bool found = false;
      for(int ss=0;ss<nw;ss++){
	if(abs(it.wav[ww]*0.1 - lamb[ss]) < 1.e-5){
	  it.idx[ww] = ss;
	  found = true;
	  continue;
	} // if
      } // ss
      if(!found) cerr<<"crh::lambdaIDX: ERROR, could not found idx!"<<endl;
    } // ww
    
  } // it
  
}

/* ------------------------------------------------------------------- */

void crh::checkBounds(mdepth_t &m)
{
  
    for(size_t ii=0; ii<m.ndep; ii++){
      m.temp[ii] =  std::max(pmin[0], std::min(m.temp[ii],  pmax[0]));
      m.v[ii] =     std::max(pmin[1], std::min(m.v[ii],     pmax[1]));
      m.vturb[ii] = std::max(pmin[2], std::min(m.vturb[ii], pmax[2]));
      m.bl[ii] =     std::max(pmin[3], std::min(m.bl[ii],     pmax[3]));
      m.bh[ii] =   std::max(pmin[4], std::min(m.bh[ii],   pmax[4]));
      m.azi[ii] =   std::max(pmin[5], std::min(m.azi[ii],   pmax[5]));
    }
  
  
}
