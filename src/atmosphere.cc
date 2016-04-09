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
//
using namespace std;
//
const double atmos::maxchange[6] = {1500., 5.0e5, 2.0e5, 800., phyc::PI/5, phyc::PI/5};




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
  
  bool store_pops = false; 
  

  /* --- Init perturbation --- */
  double *ipars = new double [npar];
  memcpy(&ipars[0], &pars[0], npar*sizeof(double));
  memset(&out[0], 0, nd*sizeof(double));
  //
  double pval = ipars[pp];
  double pertu = 0.0;
    
  //
  //if(input.nodes.ntype[pp] == temp_node) pertu = input.dpar * pval * 0.5;
  //else
  pertu = input.dpar * scal[pp];
  
  /* --- Centered derivatives ? --- */
  
  if(input.centder > 0){
    
    /* --- Up and down perturbations --- */
    double up = pertu * 0.75;
    double down = - pertu*0.75;

    
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

      /* -- recompute Hydro Eq. ? --- */
      
      //if(input.nodes.ntype[pp] == temp_node && input.thydro == 1)
	m.getPressureScale(input.boundary, eos);
      
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
      
      // if(input.nodes.ntype[pp] == temp_node && input.thydro == 1)
	m.getPressureScale(input.boundary, eos);
      
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
    //  if((input.nodes.ntype[pp] == temp_node) && (input.thydro == 1))
      m.getPressureScale(input.boundary, eos);
    
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
      pertu =  2.0 * (rnum - 0.5) * scal[pp];
    else if(n.ntype[pp] == vturb_node)
      pertu =  2.0 * (rnum - 0.5) * scal[pp];
    else if(n.ntype[pp] == b_node)
      pertu =  rnum * scal[pp];
    else if(n.ntype[pp] == inc_node)
      pertu =  rnum * scal[pp];
    else if(n.ntype[pp] == azi_node)
      pertu =  rnum * scal[pp];   
      
    pars[pp] = checkParameter(pars[pp] + pertu, pp);
  }

  
}


int getChi2(int npar1, int nd, double *pars1, double *dev, double **derivs, void *tmp1){

  
  /* --- Cast tmp1 into a double --- */
  
  atmos &atm = *((atmos*)tmp1); 
  double *ipars = new double [npar1]();
  mdepth &m = *atm.imodel;

  
  /* --- Expand atmosphere ---*/
  
  for(int pp = 0; pp<npar1; pp++)
    ipars[pp] = pars1[pp] * atm.scal[pp]; 
  
  m.expand(atm.input.nodes, &ipars[0], atm.input.dint);
  m.getPressureScale(atm.input.boundary, atm.eos);

  
  /* --- Compute synthetic spetra --- */
  
  memset(&atm.isyn[0], 0, nd*sizeof(double));
  atm.synth( m , &atm.isyn[0], (cprof_solver)atm.input.solver, true);


  
 /* --- Compute derivatives? --- */
  if(derivs){
    
    for(int pp = npar1-1; pp >= 0; pp--)
      if(derivs[pp]){
	//	cerr<<"COMPUTING RF for PAR="<<pp<<" / "<<npar1-1<<endl;
	
	/* --- Compute response function ---*/
	memset(&derivs[pp][0], 0, nd*sizeof(double));
	atm.responseFunction(npar1, m, &ipars[0], nd,
			      &derivs[pp][0], pp, &atm.isyn[0]);
	
	atm.spectralDegrade(atm.input.ns, (int)1, nd, &derivs[pp][0]);

	/* --- renormalize the response function by the 
	   scaling factor and divide by the noise --- */
	
	for(int ii = 0; ii<nd; ii++)
	  derivs[pp][ii] *= (atm.scal[pp] / atm.w[ii]);
	
      }
    
  }


  /* --- Degrade synthetic spectra --- */
  
  atm.spectralDegrade(atm.input.ns, (int)1, nd, &atm.isyn[0]);


   
  /* --- Compute difference ---*/
  
  double ichi = 0.0;
  for(unsigned ii=0; ii<nd; ii++) {
    dev[ii] = - (atm.isyn[ii] - atm.obs[ii]) / atm.w[ii]; 
    ichi += dev[ii] * dev[ii];
  }
  ichi /= (double)nd;


  
  /* --- DEBUG --- */
  static FILE *chif;
  static bool firsttime = true;
  if(!derivs){
    if(firsttime){
      string dnam = string("chi2_")+to_string(atm.input.myrank)+string(".txt");
      chif = fopen(dnam.c_str(), "w");
      firsttime = false;
    }
    fprintf(chif,"CHI2=%f\n", ichi);
    fflush(chif);
    /*
    //    mdepth &m = *atm->imodel;
    string dnam = string("imodel_")+to_string(atm.input.myrank)+string(".txt");
    FILE *mid = fopen(dnam.c_str(),"w");

    for(int kk=0; kk<m.ndep;kk++){
      fprintf(mid, "%e %e %e %e %e %e %e %e %e %e",
	      m.ltau[kk], m.z[kk], m.temp[kk], m.pgas[kk], m.rho[kk], m.v[kk], m.vturb[kk], m.b[kk], m.inc[kk], m.azi[kk]);
      //for(int bb=0;bb<npar1;bb++) fprintf(mid, "%e ", derivs[bb][kk]);
      fprintf(mid,"\n");
    }
    fclose(mid);
    */
  }
  /*else{
    //  mdepth &m = *atm->imodel;
    string dnam1 = string("rf_")+to_string(atm.input.myrank)+string(".txt");

    FILE *mid1 = fopen(dnam1.c_str(),"w");


    for(int kk=0; kk<npar1;kk++){
      for(int dd=0; dd<nd; dd++)
	fprintf(mid1, "%e ", derivs[kk][dd]);
      fprintf(mid1," \n");
    }
    fclose(mid1);
  }
  */
  /* --- END DEBUG --- */

  
  

  /* --- clean up --- */
  
  delete [] ipars;
  if(derivs)  atm.cleanup();

  
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
  lm.xtol = 1.e-3;
  lm.verb = true;
  lm.ilambda = 0.1;
  lm.maxreject = 6;
  lm.svd_thres = 1.e-16;
  
  /* ---  Set parameter limits --- */
  for(int pp = 0; pp<npar; pp++){
    lm.fcnt[pp].limit[0] = mmin[pp]/scal[pp];
    lm.fcnt[pp].limit[1] = mmax[pp]/scal[pp];
    lm.fcnt[pp].scl = 1.0;//scal[pp];
    
    if(input.nodes.ntype[pp] == azi_node){
      lm.fcnt[pp].cyclic = true;
    }else{
      lm.fcnt[pp].cyclic = true;
    }
 
    lm.fcnt[pp].bouncing = false;
    lm.fcnt[pp].capped = 1;
    lm.fcnt[pp].maxchange = maxc[pp]/scal[pp];
    
    pars[pp] = checkParameter(pars[pp], pp);
  }
  

  
  
  /* --- Loop iters --- */
  
  for(int iter = 0; iter < input.nInv; iter++){
    
    /* --- init parameters for this inversion --- */
    
    memcpy(&ipars[0], &pars[0], npar * sizeof(double));
    if(iter > 0) randomizeParameters(input.nodes , npar, &ipars[0]);

    
    /* --- Work with normalized parameters --- */
    
    for(int pp = 0; pp<npar; pp++) 
      ipars[pp] /= scal[pp];
    
    
    
    /* --- Call mpfit --- */

    double chi2 = lm.fitdata(getChi2, &ipars[0], (void*)this, input.max_inv_iter);

    /* --- Re-start populations --- */
    
    cleanup();
    

    /* --- Store result? ---*/
    
    if(chi2 < bestChi){
      bestChi = chi2;
      memcpy(&bestPars[0], &ipars[0], npar*sizeof(double));
      memcpy(&bestSyn[0], &isyn[0], ndata*sizeof(double));
    }

    /* --- reached thresthold? ---*/
    
    if(bestChi <= input.chi2_thres) break;
    
  } // iter


  /* --- re-scale fitted parameters and expand atmos ---*/
  
  for(int pp = 0; pp<npar; pp++) pars[pp] = bestPars[pp] * scal[pp]; 
  m.expand(input.nodes, &pars[0], input.dint);
  m.getPressureScale(input.boundary, eos);  
  synth( m , &obs[0], (cprof_solver)input.solver, false);
  spectralDegrade(input.ns, (int)1, input.nw_tot, &obs[0]);


  
  
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
