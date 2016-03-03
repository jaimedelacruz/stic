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
const double crh::pmax[6]  = {50000., 20.e5, 9.0e5, 5000.0, PI, 2.0*PI};
const double crh::pmin[6]  = {2400. ,-20.e5,  +0.0,   +0.0,  +0.0,  +0.0};
const double crh::pscal[6] = {2000. , 4.0e5, 3.0e5, 1000.0, 2*PI, 2*PI};
const double crh::pstep[6] = {1.e-1 , 1.e-1, 1.0e-1, 2.0e-1, 1.0e-1, 1.0e-1};

/* ----------------------------------------------------------------*/

vector<double> crh::get_max_limits(nodes_t &n){

  int nnodes = (int)n.nnodes;
  mmax.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) mmax[k] = pmax[0];
    else if(n.ntype[k] == v_node    ) mmax[k] = pmax[1];
    else if(n.ntype[k] == vturb_node) mmax[k] = pmax[2];
    else if(n.ntype[k] == b_node    ) mmax[k] = pmax[3];
    else if(n.ntype[k] == inc_node  ) mmax[k] = pmax[4];
    else if(n.ntype[k] == azi_node  ) mmax[k] = pmax[5];
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
    else if(n.ntype[k] == b_node    ) mmin[k] = pmin[3];
    else if(n.ntype[k] == inc_node  ) mmin[k] = pmin[4];
    else if(n.ntype[k] == azi_node  ) mmin[k] = pmin[5];
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
    else if(n.ntype[k] == b_node    ) scal[k] = pscal[3];
    else if(n.ntype[k] == inc_node  ) scal[k] = pscal[4];
    else if(n.ntype[k] == azi_node  ) scal[k] = pscal[5];
    else                              scal[k] = 1.0;
  }
  return scal;
}

/* ----------------------------------------------------------------*/

vector<double> crh::get_steps(nodes_t &n){

  int nnodes = (int)n.nnodes;
  step.resize(nnodes);

  double ipstep[6];
  double sum = 0.0;
  
  for(int ii=0;ii<6;ii++){
    sum += pstep[ii];
    ipstep[ii] = pstep[ii];
  }
  sum /= 6.0;
  
  for(int ii=0;ii<6;ii++){
    ipstep[ii]/= sum;
    // cerr<<ipstep[ii]<<endl;
  }
  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) step[k] = ipstep[0];
    else if(n.ntype[k] == v_node    ) step[k] = ipstep[1];
    else if(n.ntype[k] == vturb_node) step[k] = ipstep[2];
    else if(n.ntype[k] == b_node    ) step[k] = ipstep[3];
    else if(n.ntype[k] == inc_node  ) step[k] = ipstep[4];
    else if(n.ntype[k] == azi_node  ) step[k] = ipstep[5];
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

void crh::synth(mdepth_t &m_in, double *syn, cprof_solver sol, bool save_pops){

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

  
  /* --- Init storage for RH spectra --- */
  
  ospec sp = {};
  sp.I = NULL;
  sp.Q = NULL;
  sp.U = NULL;
  sp.V = NULL;
  sp.lambda = NULL;

  
  /* --- Restore nHtot and convert model to SI units --- */
  
  for(int kk = 0; kk < m.ndep; kk++){

    eos.read_partial_pressures(kk, frac, part, xa, xe);
    nhtot[kk] = frac[eos.IXH1-1] * part[eos.IXH1-1] * 1.e6;
    // fprintf(stderr, "%e %e %e %e %e %e %e \n",
    //    m.ltau[kk], m.z[kk]*1.e-5, m.temp[kk], m.nne[kk], nhtot[kk], m.v[kk],
    //	    m.vturb[kk]);

    m.cmass[kk] *= 10.0; // *= G_TO_KG / CM_TO_M**2
    m.rho[kk] *= 1000.; // G_TO_KG / CM_TO_M**3
    m.v[kk] *= 1.e-2; // CM_TO_M
    m.vturb[kk] *= 1.e-2; // CM_TO_M
    m.z[kk] *= 1.0e-2;    // CM_TO_M
    m.nne[kk] *= 1.e6;    // 1 / CM_TO_M**3
    m.tau[kk] = pow(10.0, m.ltau[kk]);
    m.b[kk] *= 1.0e-4; // B in tesla
  }


  int savep = 0;
  if(save_pops) savep = 1;
  
  
  /* --- Call RH --- */
  
  bool conv = rhf1d(input.mu, m.ndep, &m.temp[0], &m.rho[0], &m.nne[0], &m.vturb[0], &m.v[0],
		     &m.b[0], &m.inc[0], &m.azi[0], &m.z[0], &nhtot[0], &m.tau[0],
		    &cmass[0], 4.44, (bool_t)true, &sp, &save_pop, nlambda, &lambda[0],
		    input.myrank, savep);
  
    
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


  
  /* --- Restore units in model --- */
  /*
  for(int kk = 0; kk < m.ndep; kk++){
    m.cmass[kk] *= 0.1;
    m.rho[kk] *= 0.001;
    m.v[kk] *= 1.0e2;
    m.vturb[kk] *= 1.0e2;
    m.z[kk] *= 1.0e2;    
    m.nne[kk] *= 1.0e-6;
    m.b[kk] *= 1.e4;
  }
  */
  

  /* --- Deallocate sp --- */
  
  if(sp.I != NULL) delete [] sp.I;
  if(sp.Q != NULL) delete [] sp.Q;
  if(sp.U != NULL) delete [] sp.U;
  if(sp.V != NULL) delete [] sp.V;
  if(sp.lambda != NULL) delete [] sp.lambda;
  
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
