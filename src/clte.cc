/*
  CLASS CLTE

  AUTHOR: Jaime de la Cruz Rodriguez (ISP-SU 2015)

  

*/
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "clte.h"
#include "cmemt.h"
#include "input.h"
#include "physical_consts.h"
#include "cprofiles2.h"
//
using namespace std;
using namespace phyc;
//
const double clte::pmax[7]  = {50000., 20.e5, 7.0e5, 5000.0, PI, PI, 100.0};
const double clte::pmin[7]  = {2050. ,-20.e5,  +0.0,   +0.0,  +0.0,  +0.0, 0.01};
const double clte::pscal[7] = {500. , 1.0e5, 1.0e5, 500.0, PI, PI, 1.0};
const double clte::pstep[7] = {1.e0 , 1.e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0};
//
const double clte::lte_const = PI*EE*EE/(ME*CC); // In units of freq.


// -------------------------------------------------------------------------
// Get scaling of parameters (only for inversions)
// -------------------------------------------------------------------------
vector<double> clte::get_max_limits(nodes_t &n){

  int nnodes = (int)n.nnodes;
  mmax.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) mmax[k] = pmax[0];
    else if(n.ntype[k] == v_node    ) mmax[k] = pmax[1];
    else if(n.ntype[k] == vturb_node) mmax[k] = pmax[2];
    else if(n.ntype[k] == b_node    ) mmax[k] = pmax[3];
    else if(n.ntype[k] == inc_node  ) mmax[k] = pmax[4];
    else if(n.ntype[k] == azi_node  ) mmax[k] = pmax[5];
    else if(n.ntype[k] == pgas_node ) mmax[k] = pmax[6];
    else                              mmax[k] = 0;
  }
  return mmax;
}
vector<double> clte::get_min_limits(nodes_t &n){

  int nnodes = (int)n.nnodes;
  mmin.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) mmin[k] = pmin[0];
    else if(n.ntype[k] == v_node    ) mmin[k] = pmin[1];
    else if(n.ntype[k] == vturb_node) mmin[k] = pmin[2];
    else if(n.ntype[k] == b_node    ) mmin[k] = pmin[3];
    else if(n.ntype[k] == inc_node  ) mmin[k] = pmin[4];
    else if(n.ntype[k] == azi_node  ) mmin[k] = pmin[5];
    else if(n.ntype[k] == pgas_node ) mmin[k] = pmin[6];
    else                              mmin[k] = 0;
  }
  return mmin;
}
vector<double> clte::get_scaling(nodes_t &n){

  int nnodes = (int)n.nnodes;
  scal.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) scal[k] = pscal[0];
    else if(n.ntype[k] == v_node    ) scal[k] = pscal[1];
    else if(n.ntype[k] == vturb_node) scal[k] = pscal[2];
    else if(n.ntype[k] == b_node    ) scal[k] = pscal[3];
    else if(n.ntype[k] == inc_node  ) scal[k] = pscal[4];
    else if(n.ntype[k] == azi_node  ) scal[k] = pscal[5];
    else if(n.ntype[k] == pgas_node ) scal[k] = pscal[6];
    else                              scal[k] = 1.0;
  }
  return scal;
}
vector<double> clte::get_steps(nodes_t &n){

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
    else if(n.ntype[k] == b_node    ) step[k] = ipstep[3];
    else if(n.ntype[k] == inc_node  ) step[k] = ipstep[4];
    else if(n.ntype[k] == azi_node  ) step[k] = ipstep[5];
    else if(n.ntype[k] == pgas_node ) step[k] = ipstep[6];
    else                              step[k] = 1.0;
  }
  return step;
}
// -------------------------------------------------------------------------
// Constructor
// -------------------------------------------------------------------------
clte::clte(iput_t &inpt, double grav): atmos(inpt, grav){
  /*
  input.lines = inpt.lines;
  input.regions = inpt.regions;
  input.solver = inpt.solver;
  input.mu = inpt.mu;
  */
  input = inpt;
    
  /* --- Init Zeeman splitting and strength--- */
  
  nlines = input.lines.size();
  for(auto &it: input.lines)
    prof.init_zeeman_components(it);
 
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
      lambda[kk++] = convl(it.wav[ii]);
    }
  }


  /* --- add lambda = 5000 A at the end of the array for reference --- */
  lambda.push_back(inv_convl(5000.0));


  /* --- Init limits for inversion if nodes are present --- */

  if(input.nodes.nnodes > 0){
    vector<double> dummy;
    dummy = this->get_scaling(input.nodes);
    dummy = this->get_max_limits(input.nodes);
    dummy = this->get_min_limits(input.nodes);
    dummy = this->get_max_change(input.nodes);
  }
  
  
}

// -------------------------------------------------------------------------
// LTE opacity, combination of Mihalas (1971), pag. 68 - Eq. 3.4 &
// Rutten (2003) eq. 2.98, pag. 31
// -------------------------------------------------------------------------
double clte::lte_opac(double temp, double n_u, double gf, double elow, double nu0){
  return lte_const * gf * n_u * exp(-elow / (BK * temp)) *
    (1.0 - exp( -(HH * nu0) / (BK * temp)));
}

// -------------------------------------------------------------------------
// Synthesize profiles given a depth-stratified model
// -------------------------------------------------------------------------
bool clte::synth(mdepth &m, double *syn, cprof_solver sol, bool store_pops){
  string inam = "clte::synth: ";
  
  int ndep =   m.ndep;
  int nw =     (int)lambda.size(); 

  /* --- Init arrays in class cprofiles ---*/
  
  prof.init(nw, ndep);
  //prof.set_zero_abmat();

  
  /* --- Init sizes --- */
  vector<double> scatt;
  scatt.resize(nw);
  
  double lineop = 0;
  double damping = 0;

  vector<float> part, frac;
  float na, ne;
  
  prof.sf.resize(ndep);
  memset(&prof.sf[0], 0, ndep*sizeof(double));
  memset(&scatt[0], 0, nw*sizeof(double));
  
  /* --- Loop in height and get things that depend on the EOS --- */
  for(int k = 0; k<ndep; k++){

    eos.read_partial_pressures(k, frac, part, na, ne);
    
    /* --- Campute contop. for all lambdas --- */
    eos.contOpacity(m.temp[k], nw,  &lambda[0], &prof.mki[k][0], &scatt[0], frac, na, ne);

    
    /* --- Store output for later, remember that eos.fract is in fact 
       partial number density / partition function (must multiply by pf) to get the number
       density. We must run with mode 0 or the routines in contop.f90 do not work 
       --- */
    double nh = frac[eos.IXH1-1] * part[eos.IXH1-1];
    double nhe = frac[eos.IXHE1-1] * part[eos.IXHE1-1];
    
    /* --- Loop regions and compute profiles for each wavelength--- */
    for(auto &it: input.regions){

      for(int w = 0; w< it.nw; w++){ // Loop lambda
	  for(auto &li: input.lines){ // Loop lines
	    if(fabs(it.wav[w] - li.w0) > li.width) continue;

	    /* --- get absopt. coeff in LTE--- */
	    lineop = lte_opac(m.temp[k], (double)frac[li.off], li.gf, li.e_low, li.nu0);

	    /* --- damping --- */
	    double dlnu = prof.get_doppler_factor(m.temp[k], m.vturb[k], li.amass) * li.nu0; //doppler_width
	    damping = prof.damp(li, m.temp[k], m.vturb[k], m.nne[k], nh, nhe, dlnu);

	    /* --- Compute Voigt-Faraday Profiles, stored in variables of the prof.class --- */
	    prof.zeeman_profile(it.nu[w], li, m.v[k], m.b[k], dlnu, damping);

	    /* --- Now get the terms of the ABS. Matrix, stored internally in the cprofile class --- */
	    prof.zeeman_opacity( m.inc[k], m.azi[k], lineop, k, w + it.off);
      
	  } // lines
      } // w
    } // regions
  } // k

  
  /* --- Loop regions and compute profiles for each wavelength--- */
  for(auto &it: input.regions){
    for(int w = 0; w< it.nw; w++){ // Loop lambda
      prof.set_zero();
      
      for(int k = 0; k<ndep; k++){
	
	/* --- LTE source function --- */
	prof.sf[k] = prof.plank_nu(it.nu[w], m.temp[k]);

	/* --- Normalize elements of the abs. Matrix by ki and store in the vector version of the matrix elements --- */
	double iki =  prof.mki[k][w+it.off];
	prof.ki[k] = iki;
	prof.kq[k] = prof.mkq[k][w+it.off] / iki;
	prof.ku[k] = prof.mku[k][w+it.off] / iki;
	prof.kv[k] = prof.mkv[k][w+it.off] / iki;
	prof.fq[k] = prof.mfq[k][w+it.off] / iki;
	prof.fu[k] = prof.mfu[k][w+it.off] / iki;
	prof.fv[k] = prof.mfv[k][w+it.off] / iki;

	/* --- If the height scale is tau_500 then normalize ki by the opacity at 500 nm --- */
	if((sol == bez_ltau || (sol == lin_ltau))) prof.ki[k] /= prof.mki[k][nw-1];
      } // k

      /* --- Compute formal solution at this wavelength, select method --- */
      double *iprof = &syn[4 * (w + it.off)];
      if(     sol == bez_z)    prof.delobez3(ndep, m.z,   iprof, input.mu);
      else if(sol == bez_ltau) prof.delobez3(ndep, m.tau, iprof, input.mu);
      else if(sol == lin_z)    prof.delolin( ndep, m.z,   iprof, input.mu);
      else if(sol == lin_ltau) prof.delolin( ndep, m.tau, iprof, input.mu);
      else{
	cerr << inam << "ERROR, solver [" << sol << "] not implemented, exiting" << endl;
	exit(0);
      }

      for(int ss = 0; ss<4; ss++ ) iprof[ss] /= it.cscal;
      
    } // w
  } // regions

  
  /* --- Deallocate profiles --- */
  
  prof.cleanup();

  return true;
}
