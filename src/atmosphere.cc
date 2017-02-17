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
const double atmos::maxchange[7] = {2500., 5.0e5, 2.0e5, 500., phyc::PI/5, phyc::PI/5, 0.6};



vector<double> atmos::get_max_change(nodes_t &n){

  int nnodes = (int)n.nnodes;
  maxc.resize(nnodes);

  for(int k = 0; k<nnodes; k++){
    if     (n.ntype[k] == temp_node ) maxc[k] = maxchange[0];
    else if(n.ntype[k] == v_node    ) maxc[k] = maxchange[1];
    else if(n.ntype[k] == vturb_node) maxc[k] = maxchange[2];
    else if(n.ntype[k] == b_node    ) maxc[k] = maxchange[3];
    else if(n.ntype[k] == inc_node  ) maxc[k] = maxchange[4];
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

  mdepth m(m_in.ndep);
  m.cub.d = m_in.cub.d;
  m.bound_val = m_in.bound_val;
  
  bool store_pops = false;
  int centder = input.centder; 
  if(input.nodes.ntype[pp] == v_node) centder = 1;

  
  /* --- Init perturbation --- */
  
  double *ipars = new double [npar];
  memcpy(&ipars[0], &pars[0], npar*sizeof(double));
  memset(&out[0], 0, nd*sizeof(double));
  //
  double pval = ipars[pp];
  double pertu = 0.0;
    
  //
  if(input.nodes.ntype[pp] == temp_node) pertu = input.dpar * pval * 0.4;
  else                                   pertu = input.dpar * scal[pp];

  
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
      m.expand(input.nodes, &ipars[0], input.dint);
      checkBounds(m);

      /* -- recompute Hydro Eq. ? --- */
      
      //if(input.nodes.ntype[pp] == temp_node && input.thydro == 1)
	m.getPressureScale(input.boundary, eos);
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
      m.expand(input.nodes, &ipars[0], input.dint);
      checkBounds(m);

      // if(input.nodes.ntype[pp] == temp_node && input.thydro == 1)
	m.getPressureScale(input.boundary, eos);
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
    
    m.expand(input.nodes, &ipars[0],  input.dint);
    checkBounds(m);

    //  if((input.nodes.ntype[pp] == temp_node) && (input.thydro == 1))
    m.getPressureScale(input.boundary, eos);
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
      pertu = 2.0 * (rnum - 0.5) * scal[pp] * 0.1;
    else if(n.ntype[pp] == v_node)
      pertu =  (rnum - 0.5) * scal[pp];
    else if(n.ntype[pp] == vturb_node)
      pertu =  (rnum - 0.5) * scal[pp];
    else if(n.ntype[pp] == b_node)
      pertu =  rnum * scal[pp];
    else if(n.ntype[pp] == inc_node)
      pertu =  rnum * scal[pp]*phyc::PI;
    else if(n.ntype[pp] == azi_node)
      pertu =  rnum * scal[pp]*phyc::PI;   
    else if(n.ntype[pp] == pgas_node)
      pertu = rnum;
      
    pars[pp] = checkParameter(pars[pp] + pertu, pp);
  }

  
}
void getDregul2(double *m, int npar, double *dregul, nodes_t &n)
{

  std::vector<double> tmp, tmp1;
  double penalty = 0.0, weights[3] = {0.15,10.0,5.0};

  
  /* --- Tikhonov's regularization on first derivative for Temp --- */
  
  if(n.toinv[0] && true){
    int nn = (int)n.temp.size();
    
    if(nn > 1){
      tmp.resize(nn, 0.0);
      tmp1.resize(nn, 0.0);

      mth::cent_der<double>(nn,&n.temp[0], &m[n.temp_off], &tmp[0]);
      mth::cent_der<double>(nn,&n.temp[0], &tmp[0], &tmp1[0]);
      mth::cent_der<double>(nn,&n.temp[0], &tmp1[0], &dregul[n.temp_off]);
      mth::cmul<double>(nn, &dregul[n.temp_off], 2*weights[0]/nn);
      penalty += weights[0] * mth::ksum2(nn, &tmp1[0]) / nn;
    }
  }

  
  /* --- Tikhonov's regularization on first derivative for vlos --- */
  
  if(n.toinv[1] && true){
    int nn = (int)n.v.size();
    
    if(nn > 1){
      tmp.resize(nn, 0.0);
      mth::cent_der<double>(nn,&n.v[0], &m[n.v_off], &tmp[0]);
      mth::cent_der<double>(nn,&n.v[0],        &tmp[0], &dregul[n.v_off]);
      mth::cmul<double>(nn, &dregul[n.v_off], 2*weights[1]/nn);
      penalty +=weights[1] * mth::ksum2(nn, &tmp[0]) / nn;
    }
  }
  
  /* --- Penalize deviations from zero for Vturb --- */

  if(n.toinv[2] && true){
    int nn = (int)n.vturb.size();          
    double sum2 = 0.0;
    for(size_t ii=0; ii<nn; ii++) {
      dregul[n.vturb_off+ii] = 2.0 * weights[2] * m[n.vturb_off+ii] / nn;
      sum2 +=  mth::sqr( m[n.vturb_off+ii]);
    }
    penalty +=  weights[2] * sum2 / nn;
  }
  
  
  dregul[npar] = penalty;
  
}

void getDregul(mdepth &m, int npar, double *dregul, nodes_t &n)
{
  const double scale[2] = {0.05, 0.05};
  double penalty = 0.0, sum = 0.0, total = 0.0;
  double bla[2] = {0.0,0.0};
  
  if(n.toinv[2]){
    
    /* --- Vturb regularization: penalize deviations from zero --- */
    
    for(size_t ii=0; ii<m.ndep; ii++) penalty += mth::sqr((m.vturb[ii]-sum)*1.e-5);
    penalty *= scale[0]/m.ndep;
    
    
    /* --- Now derivative of vturb --- */
    
    for(size_t ii=0; ii<m.ndep; ii++) bla[0] += 2.0e-5 * (m.vturb[ii]-sum);
    bla[0] *= scale[0]/m.ndep;
    
  } else bla[0] = 0.0;
  
  
  if(n.toinv[1]){
    
    
    /* --- Vlos regularization: penalize deviations from the mean value --- */
    
    sum = 0.0;
    total = mth::sum(m.ndep, m.v) / m.ndep;
    for(size_t ii=0; ii<m.ndep; ii++) sum += mth::sqr((m.v[ii]-total)*1.e-5);
    penalty +=  sum * scale[1] / m.ndep;
    
    
    /* --- Now derivative of Vlos --- */
    
    for(size_t ii=0; ii<m.ndep; ii++) bla[1] += 2.0e-5 * (m.v[ii]-total);
    bla[1] *= scale[1]/m.ndep;
    
  }else bla[1] = 0.0;
    
    for(size_t ii=0; ii<npar; ii++) {
      if     (n.ntype[ii] == vturb_node) dregul[ii] = bla[0];
      else if(n.ntype[ii] == v_node    ) dregul[ii] = bla[1];
      else                               dregul[ii] = 0.0;
    }
  

  /* --- Store the Chi2 penalty in the last element of dregul --- */
  
  dregul[npar] = penalty;
}

int getChi2(int npar1, int nd, double *pars1, double *dev, double **derivs, void *tmp1, double *dregul){

  
  /* --- Cast tmp1 into a double --- */
  
  atmos &atm = *((atmos*)tmp1); 
  double *ipars = new double [npar1]();
  mdepth &m = *atm.imodel;

  
  /* --- Expand atmosphere ---*/
  
  for(int pp = 0; pp<npar1; pp++)
    ipars[pp] = pars1[pp] * atm.scal[pp]; 
  
  m.expand(atm.input.nodes, &ipars[0], atm.input.dint);
  atm.checkBounds(m);
  m.getPressureScale(atm.input.boundary, atm.eos);

  
  /* --- Compute synthetic spetra --- */
  
  memset(&atm.isyn[0], 0, nd*sizeof(double));
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
  
  obs = &o[0];
  w = &weights.d[0];
  imodel = &m;

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
  lm.lmin = 1.e-5;
  lm.proc = input.myrank;
  if(input.regularize >= 1.e-5){
    lm.regularize = true;
    lm.regul_scal = input.regularize;
  } else lm.regularize = false;
  
  /* ---  Set parameter limits --- */
  
  for(int pp = 0; pp<npar; pp++){
    lm.fcnt[pp].limit[0] = mmin[pp]/scal[pp];
    lm.fcnt[pp].limit[1] = mmax[pp]/scal[pp];
    lm.fcnt[pp].scl = 1.0;//scal[pp];
    
    if(input.nodes.ntype[pp] == azi_node) lm.fcnt[pp].cyclic = true;
    else                                  lm.fcnt[pp].cyclic = false;
    
    lm.fcnt[pp].bouncing = false;
    lm.fcnt[pp].capped = 1;
    
    if(input.nodes.ntype[pp] == temp_node){
      lm.fcnt[pp].relchange = true;
      lm.fcnt[pp].maxchange[0] = 0.25;
      lm.fcnt[pp].maxchange[1] = 2.0;
    }else{
      lm.fcnt[pp].relchange = false;
      lm.fcnt[pp].maxchange[0] = maxc[pp]/scal[pp];
      lm.fcnt[pp].maxchange[1] = maxc[pp]/scal[pp];
    }
    pars[pp] = checkParameter(pars[pp], pp);
  }
  

  
  /* --- Loop iters --- */
  
  for(int iter = 0; iter < input.nInv; iter++){

    cleanup();
    
    /* --- init parameters for this inversion --- */
    
    memcpy(&ipars[0], &pars[0], npar * sizeof(double));
    if(iter > 0) randomizeParameters(input.nodes , npar, &ipars[0]);

    
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
    }

    /* --- reached thresthold? ---*/
    
    if(bestChi <= input.chi2_thres) break;
    
  } // iter


  /* --- re-scale fitted parameters and expand atmos ---*/
  
  for(int pp = 0; pp<npar; pp++)
    pars[pp] = bestPars[pp] * scal[pp];
  
  memset(&isyn[0],0,ndata*sizeof(double));
  m.expand(input.nodes, &pars[0], input.dint);
  checkBounds(m);
  m.getPressureScale(input.boundary, eos); 
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
