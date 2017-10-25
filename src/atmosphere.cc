#include <vector>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <functional>
#include "atmosphere.h"
#include "spectral.h"
#include "ceos.h"
#include "input.h"
#include "physical_consts.h"
#include "clm.h"
#include "math_tools.h"
//
using namespace std;
//
const double atmos::maxchange[7] = {3000., 5.0e5, 1.5e5, 600., 600., phyc::PI/5, 0.6};



vector<double> atmos::get_max_change(nodes_t &n){

  int nnodes = (int)n.nnodes;
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
    if(input.nodes.ntype[pp] == temp_node) pertu = input.dpar * pval * 0.25;
    else                                   pertu = input.dpar * scal[pp];
  }else  pertu = input.dpar * scal[pp];
    
  
  /* --- Centered derivatives ? --- */
  
  if(centder > 0){
    
    /* --- Up and down perturbations --- */
    double up = pertu * 0.5;
    double down = - pertu * 0.5;

    
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

      synth(m, &out[0], (cprof_solver)input.solver, store_pops);
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

      synth(m, &spec[0], (cprof_solver)input.solver, store_pops);
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

    
    /* --- Synth --- */
    
    m.expand(input.nodes, &ipars[0],  input.dint, input.depth_model);
    checkBounds(m);

    //  if((input.nodes.ntype[pp] == temp_node) && (input.thydro == 1))
    m.getPressureScale(input.nodes.depth_t, input.boundary, eos);
    //m.nne_enhance(input.nodes, npar, &ipars[0], eos);
      
    synth(m, &out[0], (cprof_solver)input.solver, store_pops);
    
    
    /* --- Finite differences ---*/
    
    for(int ii = 0; ii<nd; ii++)
      out[ii]  =  (out[ii] - syn[ii]) / pertu;
  }
  
  
  /* --- clean-up --- */
  
  delete [] ipars;
    
}

void atmos::randomizeParameters(nodes_t &n, int npar, double *pars){
  srand (time(NULL));
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
      rnum = (double)rand() / RAND_MAX;
      last = n.ntype[pp];
    }
    
    if(n.ntype[pp] == temp_node)
      pertu = 2.0 * (rnum - 0.5) * scal[pp] * 0.2;
    else if(n.ntype[pp] == v_node)
      pertu =  (rnum - 0.5) * scal[pp];
    else if(n.ntype[pp] == vturb_node){
      pertu =  (rnum - 0.75) * 2.e5;
    }else if(n.ntype[pp] == bl_node)
      pertu =  (rnum-0.35) * scal[pp];
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
  }

  
}

/* --------------------------------------------------------------------------------------------------- */

double const_dregul(int n, double *ltau, double *var, double weight, double *res, double me, bool divide_by_n)
{
  
  /* --- 
     Penalize deviations from me, removed the factor x2 from the derivative
     to be consistent with our LM definition 
     --- */
  double penalty = 0.0, nn = (double)n;
  
  double c1 = (divide_by_n && (n > 0)) ? ((nn-1.0) / (nn)) : 1.0;
  
  
  if(n > 1){
    double c =  weight / fabs(ltau[n-1] - ltau[0]);
    
    for(int ii = 0; ii<n; ii++){
      //
      int i0=std::max(ii-1, 0);
      int i1=std::min(ii+1, n-1);
      //
      double trange = fabs(ltau[i0]-ltau[i1]) * 0.5;
      double tmp = var[ii] - me;
      //
      res[ii] = c * trange * tmp;
      penalty += c1 * c * trange * tmp * tmp;
    }
  }else{
    double tmp = var[0]-me;
    penalty = weight / 7.0 * tmp*tmp;
    res[0] = weight  / 7.0 * tmp;
  }

  return penalty;
}

/* --------------------------------------------------------------------------------------------------- */

double tikhonov1_dregul(int n, double *ltau, double *var, double weight, double *res)
{

  if(n == 1){
    res[0] = 0.0;
    return 0.0;
  }
  
  double penalty = 0.0, c = weight / (n-1.0);
  
  
  /* --- Derivative in the first and last points:

     Penalty = c * ((var[k+1] - var[k])/dltau)**2

     der_0 = 2c * (var[0] - var[1]) / dltau**2
     der_i = 2c * ((var[i] -var[i-1]) / dltau[i]**2 + (var[i] - var[i+1])/dltau[i+1]**2) 
     der_n = 2c * (var[n] - var[n-1]) / dltau[n]**2

     In the LM we have divided the x2 in all terms so we do it here too in the implementation.

     --- */
  
  res[0]   =  c * (var[0] - var[1]) / mth::sqr((ltau[0] - ltau[1]));// ltau_range;
  res[n-1] =  c * (var[n-1] - var[n-2]) / mth::sqr((ltau[n-1] - ltau[n-2])) ;/// ltau_range;
  
  
  /* --- Derivative in intermediate points --- */

  if(n > 2)
    for(int kk=1;kk<n-1;kk++ ){
      double dt0 = (ltau[kk]-ltau[kk-1]), dt1 = (ltau[kk+1]-ltau[kk]);//, dtau = dt0+dt1;
      res[kk] =  c * ( (var[kk] - var[kk-1]) / (dt0*dt0) + (var[kk] - var[kk+1]) / (dt1*dt1 )) ;//* (dtau / ltau_range);
    }
  

  /* --- Penalty term --- */

  for(int ii = 1; ii<n; ii++){
    penalty += mth::sqr((var[ii] - var[ii-1]) / (ltau[ii] - ltau[ii-1]));
  }


return penalty*c;///ltau_range;
}

/* --------------------------------------------------------------------------------------------------- */

void getDregul2(double *m, int npar, double *dregul, nodes_t &n)
{
  const double weights[7] = {0.0001, 7.0, 2.0, 1.0, 1.0, 1.0, 1.5};
  double penalty = 0.0, *ltau = NULL, we = 0.0;
  nodes_type_t ntype = none_node;
  int off = 0;
  memset(dregul, 0, (npar+1)*sizeof(double)); // set derivatives to zero
  
  
  /* --- Penalize temp ? Only allow Tikhonov, the rest don't make sense for temp --- */

  if(n.toinv[0] && (n.regul_type[0] > 0)){ // temperature
    ntype = temp_node;
    off = (int)n.temp_off;
    ltau =  &n.temp[0];
    we = weights[0];
    
    int nn = (int)n.temp.size();
    switch(n.regul_type[0]){
    case(1):
      penalty += tikhonov1_dregul(nn, ltau, &m[off], we, &dregul[off]);
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
    we = weights[1];
    
    switch(n.regul_type[1]){
    case(1):
      penalty += tikhonov1_dregul(nn, ltau, &m[off], we, &dregul[off]);
      break;
    case(2):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], mth::mean(nn, &m[off]), true);
      break;
    case(3):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], 0.0, false);	    
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
    we = weights[2];

    switch(n.regul_type[2]){
    case(1):
      penalty += tikhonov1_dregul(nn, ltau, &m[off], we, &dregul[off]);
      break;
    case(2):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], mth::mean(nn, &m[off]), true);
      break;
    case(3):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], 0.0, false);	    
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
    we = weights[3];

    switch(n.regul_type[3]){
    case(1):
      penalty += tikhonov1_dregul(nn, ltau, &m[off], we, &dregul[off]);
      break;
    case(2):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], mth::mean(nn, &m[off]), true);
      break;
    case(3):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], 0.0, false);	    
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
    we = weights[4];

    switch(n.regul_type[4]){
    case(1):
      penalty += tikhonov1_dregul(nn, ltau, &m[off], we, &dregul[off]);
      break;
    case(2):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], mth::mean(nn, &m[off]), true);
      break;
    case(3):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], 0.0, false);	    
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
    we = weights[5];

    switch(n.regul_type[5]){
    case(1):
      penalty += tikhonov1_dregul(nn, ltau, &m[off], we, &dregul[off]);
      break;
    case(2):
      penalty += const_dregul(nn, ltau, &m[off], we, &dregul[off], mth::mean(nn, &m[off]), true);
      break;
    default:
      break;
    }
  }
  

  
  /* --- Penalize pgas enhancement factor ? --- */
  
  if(n.toinv[6] && (n.regul_type[6] > 0)){ // pgas

    
    /* --- Penalize deviations from 1.0 --- */

    off = (int)n.pgas_off;
    double tmp = (m[off] - 1.0);
    penalty += tmp * tmp * weights[6];
    dregul[off] = weights[6] * tmp;
  }
  
  /* --- Store it in the last element of the array --- */
  
  dregul[npar] = penalty;
  
}

/* --------------------------------------------------------------------------------------------------- */

int getChi2(int npar1, int nd, double *pars1, double *dev, double **derivs, void *tmp1, double *dregul, bool store){

  
  /* --- Cast tmp1 into a double --- */
  
  atmos &atm = *((atmos*)tmp1); 
  double *ipars = new double [npar1]();
  mdepth &m = *atm.imodel;
  
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
  bool conv = atm.synth( m , &atm.isyn[0], (cprof_solver)atm.input.solver, true);

  
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
	  derivs[pp][ii] *= (atm.scal[pp] / atm.w[ii]);
	
      }
    }    
  }


  /* --- Degrade synthetic spectra --- */
  
  atm.spectralDegrade(atm.input.ns, (int)1, nd, &atm.isyn[0]);
  


  /* --- Compute residue --- */

  for(int ww = 0; ww < nd; ww++){
    dev[ww] = (atm.obs[ww] - atm.isyn[ww]) / atm.w[ww];
  }


  /* ---  compute regularization --- */ 

  if(dregul)
    getDregul2(pars1, npar1, dregul, atm.input.nodes);
      
  /* --- clean up --- */
  
  delete [] ipars;

  return 0;
}

double atmos::fitModel2(mdepth_t &m, int npar, double *pars, int nobs, double *o, mat<double> &weights){
  
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

  
  /* --- compute tau scale ---*/
  
  for(int k = 0; k < (int)m.cub.size(1); k++) m.tau[k] = pow(10.0, m.ltau[k]);


  
  /* --- Invert pixel randomizing parameters if iter > 0 --- */
  
  int ndata = input.nw_tot * input.ns;
  double ipars[npar], bestPars[npar], ichi, bestChi = 1.e10;
  isyn.resize(ndata);
  vector<double> bestSyn;
  bestSyn.resize(ndata);

  
  /* --- Init clm --- */

  clm lm = clm(ndata, npar);
  lm.xtol = 3.e-3;
  lm.verb = input.verbose;
  if(input.marquardt_damping > 0.0) lm.ilambda = input.marquardt_damping;
  else                              lm.ilambda = 1.0;
  lm.maxreject = 6;
  lm.svd_thres = max(input.svd_thres, 1.e-16);
  lm.chi2_thres = input.chi2_thres;
  lm.lmax = 1.e4;
  lm.lmin = 1.e-4;
  lm.lfac = sqrt(10.0);
  lm.proc = input.myrank;
  if(input.regularize >= 1.e-5){
    lm.regularize = true;
    lm.regul_scal_in = input.regularize;
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
      lm.fcnt[pp].maxchange[0] = 0.25;
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
  
  for(int iter = 0; iter < input.nInv; iter++){

    cleanup();

    
    
    /* --- init parameters for this inversion --- */
    
    memcpy(&ipars[0], &pars[0], npar * sizeof(double));
    if(iter > 0 || input.random_first){
      randomizeParameters(input.nodes , npar, &ipars[0]);

      if(depth_per){
	imodel->expand(input.nodes, &ipars[0], input.dint, input.depth_model);
	checkBounds(*imodel);
	imodel->getPressureScale(input.nodes.depth_t, input.boundary, eos);
      }
    }
    
    /* --- Work with normalized parameters --- */
    
    for(int pp = 0; pp<npar; pp++) 
      ipars[pp] /= scal[pp];
    
    
    
    /* --- Call clm --- */

    double chi2 = lm.fitdata(getChi2, &ipars[0], (void*)this, input.max_inv_iter);

    
    
    /* --- Re-start populations --- */
    
    cleanup();
    

    /* --- Store result? ---*/
    
    if(chi2 < bestChi){
      bestChi = chi2;
      memcpy(&bestPars[0], &ipars[0], npar*sizeof(double));
      memcpy(&best_m.cub.d[0], &m.cub.d[0], m.ndep*14*sizeof(double));
    }

    if(depth_per){
      memcpy(&m.cub.d[0], &m2.cub.d[0], m.ndep*14*sizeof(double));
      memcpy(&m1.cub.d[0], &m2.cub.d[0], m.ndep*14*sizeof(double));
    }
    
    /* --- reached thresthold? ---*/
    
    if(bestChi <= input.chi2_thres) break;
    
  } // iter

  

  
  /* --- re-scale fitted parameters and expand atmos ---*/
  
  for(int pp = 0; pp<npar; pp++)
    pars[pp] = bestPars[pp] * scal[pp];

  memset(&isyn[0],0,ndata*sizeof(double));
  if(!depth_per){
    m.expand(input.nodes, &pars[0], input.dint, input.depth_model);
    checkBounds(m);
    m.getPressureScale(input.nodes.depth_t , input.boundary, eos);
  } else  memcpy(&m.cub.d[0], &best_m.cub.d[0], m.ndep*14*sizeof(double));

  synth( m , &isyn[0], (cprof_solver)input.solver, false);
  spectralDegrade(input.ns, (int)1, input.nw_tot, &isyn[0]);
  
  
  double sum = 0.0;
  for(int ww = 0; ww<ndata;ww++){
    double tmp = (obs[ww] - isyn[ww]) / weights.d[ww];
    sum += (tmp*tmp);
  }
  fprintf(stderr,"Recomp chi2=%13.5f\n", sum/ndata);
  memcpy(&obs[0], &isyn[0], ndata*sizeof(double));
  
  /* --- Clean-up --- */
  
  isyn.clear();
  bestSyn.clear();
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
