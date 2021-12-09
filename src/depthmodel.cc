#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <cmath>
//
#include "depthmodel.h"
#include "interpol.h"
#include "ceos.h"
#include "physical_consts.h"
#include "math_tools.h"

using namespace std;

/* --- Default boundary condition (if none provided) --- */

const double mdepth::boundary_pgas_default = 7.e-1;

/* ------------------------------------------------------ */

void mdepth::to_txt(const std::string &fname)
{
  FILE *fout = fopen(fname.c_str(), "w");
  for(int ii=0; ii<ndep; ++ii){
    for(int ss = 0; ss< 14; ++ss)
      fprintf(fout, "%e  ", cub(ss,ii));
    fprintf(fout, "\n");
  }
  fclose(fout);
}

/* ------------------------------------------------------ */

std::vector<double> mdepth::model2vector()
{
  vector<double> res(ndep*12+2, 0.0);
  
  if(ndep >0)
    memcpy(&res[0], &cub(0,0), ndep*12*sizeof(double));

  int  off = ndep*12;
  res[off++] = tr_loc;
  res[off++] = tr_amp;

  return res; 
}

/* ------------------------------------------------------ */

void mdepth::vector2model(std::vector<double> &vec)
{

  
  int nn = (int)vec.size()/12;
  if(ndep != nn) setsize(nn);

  memcpy(&cub(0,0), &vec[0], 13*ndep*sizeof(double));
  int  off = ndep*12;
  tr_loc = vec[off++];
  tr_amp = vec[off++];
  
}

/* ------------------------------------------------------ */

void var_int(int n, double *x, double *y, double *xx, int k0, int k1, bool use_log) // inplace
{


  double *tmp = new double [n];
  int nn = k1 - k0 + 1;
  
  if(!use_log){
    hermpol(nn, &x[k0], &y[k0], n, xx, tmp, false);
  }else{
    int idmi = mth::argmin(n, y);
    double imi = y[idmi];
    
    use_log = ((imi > 0.0) ? true : false);
    
    if(use_log)
      for(int ii=k0;ii<=k1;ii++) y[ii] = log(y[ii]);
    
    hermpol(nn, &x[k0], &y[k0], n, xx, tmp, false);
    
    if(use_log)
      for(int ii=0; ii<n; ii++) tmp[ii] = exp(tmp[ii]);
  }
  
  memcpy(y, tmp, n*sizeof(double));
  
  
  delete [] tmp;
}

/* ------------------------------------------------------ */

void mdepth::optimize_depth(eoswrap &eos, float tcut, int smooth)
{
  static const double grph=2.26e-24, crhmbf=2.9256e-17;

  bool reset_tau = false;
  vector<double> aind(ndep, 0.0), xind(ndep, 0.0), tmp(ndep, 0.0);
  double tdiv = 0., rdiv = 0.0, taudiv = 0.0, lg11 = log10(1.1); 
  int k0 = 0, k1 = ndep-1;
  double *dep = NULL;

  /* --- Get temperature cut --- */
  
  for(int k=0;k<=k1;k++){
    if(temp[k] > tcut){
      k0 = k;
      continue;
    }else break;
  }

  double range = double(k1-k0)/(ndep-1.0);
  for(int k=0; k<ndep; k++) xind[k] = double(k)*range +k0;

  this->fill_densities(eos, true, 0, k1);
  

  /* --- compute approximate tau scale --- */
  
  if((ltau[0] == 0.) && (nne[0] != 0.) && rho[0] !=0.){
    double okap, kap, otau, tau;
    reset_tau = true;
    
    tau = 1.e-9;
    kap =  1.03526e-16*nne[0]*crhmbf/pow(temp[0], 1.5)*
      exp(0.754*phyc::EV/phyc::BK/temp[0])*rho[0]/grph;
    ltau[0] = tau;
    
    for(int ii=1;ii<ndep;ii++){
      otau = tau, okap = kap;
      kap =  1.03526e-16*nne[ii]*crhmbf/pow(temp[ii], 1.5)*
	exp(0.754*phyc::EV/phyc::BK/temp[ii])*rho[ii]/grph;

      ltau[ii] = ltau[ii-1] + 0.5 * fabs(z[ii]-z[ii-1]) * (kap + okap);
    }
    
    for(int k=0; k<ndep; k++){
      ltau[k] = log10(ltau[k]);
      //cerr<<ltau[k]<<endl;
    }
    // exit(0);
  }
  
  if(ltau[0] != 0)       dep = ltau;
  else if(cmass[0] != 0) dep = cmass;



  
  /* --- Get location of strong gradients --- */
  
  for(int k=k0+1;k<=k1;k++){
    tdiv = fabs(log10(temp[k]) - log10(temp[k-1])) / lg11;
    rdiv = fabs(log10(rho[k])  - log10(rho[k-1]))  / lg11;
    taudiv = fabs(dep[k] - dep[k-1]) * 10.;/// 0.1;
    aind[k] = aind[k-1] + std::max(std::max(tdiv,rdiv), taudiv);
  }

  /* --- Smooth gradients a bit --- */
  
  float rat = double(k1-k0) / aind[k1];
  for(int k=0; k<=k1; k++) aind[k] = aind[k] * rat + k0;
  if(smooth) mth::smooth(k1-k0+1, &aind[k0], smooth);
  
  
  /* --- Check that all points are separated --- */

  for(int ii=k0+1; ii< k1; ii++)
    if(fabs(aind[ii] - aind[ii-1]) < 1.e-3)
      aind[ii] = 0.5 * (aind[ii-1] + aind[ii+1]);
  
  
  /* --- Get new grid --- */

  int ndep1 = k1 - k0 + 1;
  var_int(ndep, &aind[0], temp, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],    v, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],vturb, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],   bl, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],   bh, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],  azi, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],  rho, &xind[0], k0, k1, true);
  var_int(ndep, &aind[0], pgas, &xind[0], k0, k1, true);
  var_int(ndep, &aind[0],  nne, &xind[0], k0, k1, true);
  var_int(ndep, &aind[0],  pel, &xind[0], k0, k1, true);
  var_int(ndep, &aind[0],    z, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0], ltau, &xind[0], k0, k1, false);
  var_int(ndep, &aind[0],cmass, &xind[0], k0, k1, false);

  for(int ii=0; ii<ndep; ii++) tau[ii] = pow(10.0, ltau[ii]);
  

  if(reset_tau) memset(ltau, 0, ndep*sizeof(ndep));

  // FILE *bla = fopen("bla.txt", "w");
  // for(int ii=0; ii<ndep; ii++)
  //   fprintf(bla, "%3d %8.4f %8.4f %e %e %f %e\n", ii, z[ii]*1.e-5, ltau[ii], temp[ii], rho[ii], xind[ii], aind[ii]);
  // fclose(bla);
  
  //cerr<<endl;
  // exit(0);
  
}

/* ------------------------------------------------------ */

void mdepth::optimize_depth_ltau(eoswrap &eos, float tcut)
{
  static const double grph=2.26e-24, crhmbf=2.9256e-17;

  bool reset_tau = false;
  vector<double> aind(ndep, 0.0), xind(ndep, 0.0), tmp(ndep, 0.0);
  //double tdiv = 0., rdiv = 0.0, taudiv = 0.0, lg11 = log10(1.1); 
  int k0 = 0, k1 = ndep-1;
  // double *dep = NULL;

  /* --- Get temperature cut --- */
  
  for(int k=0;k<=k1;k++){
    if(temp[k] > tcut){
      k0 = k;
      continue;
    }else break;
  }

  for(int k=0; k<=k1; k++) xind[k] = (double)k;

  if(ltau[0] == 0){
    double wav = 5000., kap= 0.0, scat=0., okap=0, itau=0.0;
    ltau[k0] = 1.e-9;

    bool use_pgas = ((pgas[k0] > 1.e-8)? true : false);
    if(use_pgas)
      eos.contOpacity_TPg(temp[k0], pgas[k0], 1, &wav, &kap, &scat, 1.e-4);
    else
      eos.contOpacity_TRho(temp[k0], rho[k0], 1, &wav, &kap, &scat, nne[k0]*temp[k0]*eos.bk);

    for(int k=k0+1; k<=k1; k++){
      okap = kap;
      if(use_pgas)
	eos.contOpacity_TPg(temp[k], pgas[k], 1, &wav, &kap, &scat, 1.e-4);
      else
	eos.contOpacity_TRho(temp[k], rho[k], 1, &wav, &kap, &scat, nne[k]*temp[k]*eos.bk);

      
      
      double dz = fabs(z[k]-z[k-1]);
      ltau[k] = ltau[k-1] + 0.5 * dz * (kap + okap);
      
      if(ltau[k] > 200.){
	k1 = k;
	break;
      }
      
    }
    
    //    ltau[k0] = exp(2.0 * log(ltau[k0+1]) - log(ltau[k0+2]));
    for(int k=k0;k<=k1;k++) ltau[k] = log10(ltau[k]);
  }
  double ma = ltau[k1], mi = ltau[k0];
  for(int ii=0; ii<ndep;ii++) aind[ii] = double(ii) / (ndep-1.0) * (ma-mi) + mi;
  
  
  var_int(ndep, ltau, temp,  &aind[0], k0, k1, false);
  var_int(ndep, ltau, v,     &aind[0], k0, k1, false);
  var_int(ndep, ltau, vturb, &aind[0], k0, k1, false);
  var_int(ndep, ltau, bl,    &aind[0], k0, k1, false);
  var_int(ndep, ltau, bh,    &aind[0], k0, k1, false);
  var_int(ndep, ltau, azi,   &aind[0], k0, k1, false);
  var_int(ndep, ltau, pgas,  &aind[0], k0, k1, true);
  var_int(ndep, ltau, rho,   &aind[0], k0, k1, true);
  var_int(ndep, ltau, nne,   &aind[0], k0, k1, true);
  var_int(ndep, ltau, pel,   &aind[0], k0, k1, true);
  var_int(ndep, ltau, cmass, &aind[0], k0, k1, false);
  var_int(ndep, ltau, z,     &aind[0], k0, k1, false);
  
  memcpy(ltau, &aind[0], ndep*sizeof(double));
  
  
  //for(int ii=0; ii<ndep;ii++) fprintf(stderr,"[%3d] %f %f\n", ii, aind[ii], temp[ii]);
    // exit(0);
  
  
  
}

/* ------------------------------------------------------ */

void mdepth::setsize(int n){

  /* --- resize --- */
  cub.set({14, n});
  ndep = n;
  
  
  /* --- Assign pointers (keep this order so we can 
     copy the buffer directly ---*/
  
  temp =  &cub(0,0);
  v =     &cub(1,0);
  vturb = &cub(2,0);
  bl =    &cub(3,0);
  bh  =   &cub(4,0);
  azi =   &cub(5,0);
  pgas =  &cub(6,0);
  rho  =  &cub(7,0);
  nne  =  &cub(8,0);
  ltau =  &cub(9,0);
  z =     &cub(10,0);
  cmass = &cub(11,0);
  pel =   &cub(12,0);
  tau  =  &cub(13,0);

  tr_loc = 0;
  tr_amp = 1.0;
  
}

void mdepth::zero(){
  cub.zero();
}
mdepth& mdepth::operator= (  mdepth &m)
{

  this->setsize(m.ndep);
  memcpy(&(this->cub(0,0)), &m.cub(0,0), 14*m.ndep*sizeof(double));

  this->tr_loc = m.tr_loc;
  this->tr_amp = m.tr_amp;
  
  this->bound_val = m.bound_val;
  this->ref_m = m.ref_m;
  
  return *this;
}

/* ------------------------------------------------------------------ */

mdepth::mdepth(const  mdepth &m)
{

  this->setsize(m.ndep);
  memcpy(&(this->cub.d[0]), &m.cub.d[0], 14*m.ndep*sizeof(double));
  
  this->tr_loc = m.tr_loc;
  this->tr_amp = m.tr_amp;
  
  this->bound_val = m.bound_val;
  this->ref_m = m.ref_m;
  
}


/* ------------------------------------------------------------------ */

//void eoswrap::hydrostatic_cmass(int ndep, double *tau, double *t, double *Pg, double *rho, double *nel,
//			     double *z, double *cmass, double *ltau

void mdepth::hydrostatic(eoswrap &eos, int depth_m)
{

  if(depth_m == 0)   
    eos.hydrostatic((int)ndep, &tau[0], &temp[0], &pgas[0], &rho[0], &nne[0], &pel[0],
		    &z[0], &cmass[0], pgas[0], (float)1.e-5);
  else
    eos.hydrostatic_cmass((int)ndep, tau, temp, pgas, rho, nne, z, cmass, ltau, pgas[0]);
  
  
}

/* ------------------------------------------------------------------ */


void mdepth::nodes2depth(int n, double *x, double *y, int nn, double *xx, double *yy, int interpol, bool extrapolate)
{
  if     (n <  1)           return;
  else if(n == 1)           for(int kk=0; kk<nn; kk++) yy[kk] = y[0];
  else if(n == 2)            linpol(n, x, y, nn, xx, yy, extrapolate);
  else if(n >= 3){
    if     (interpol == 0)   linpol(n, x, y, nn, xx, yy, extrapolate);
    else if(interpol == 1)  bezpol2(n, x, y, nn, xx, yy, extrapolate);
    else if(interpol == 2)  hermpol(n, x, y, nn, xx, yy, extrapolate);
    else if(interpol == 3){
      if(n > 3) vlint(n, x, y, nn, xx, yy);
      else      linpol(n, x, y, nn, xx, yy, extrapolate);
    }
  }
}

void mdepth::nne_enhance(nodes_t &nodes, int n, double *pars, eoswrap &eos){

  /* --- are we inverting the nne enhancement? --- */
  
  bool doit = false;
  double mult = 1.0;
  //
  for(int ii = 0; ii < nodes.nnodes; ii++)
    if(nodes.ntype[ii] == pgas_node){
      doit = true;
      mult = pars[ii];
    }
  if(!doit) return;

  //fprintf(stderr,"   [mult=%f]\n", mult);
  
  /* --- Enhance electron pressure from ltau_500 <= -4.5 --- */

  const double dx = 0.2, x0 = -4.8;
  
  for(int ii = 0; ii<ndep; ii++){
    double at = -tanh((ltau[ii]-x0)*phyc::PI/(dx)) * 0.5 + 0.5;
    double corr = mult * at + (1.0 - at);
    nne[ii] *= corr;

    eos.nne_from_T_Pg_nne (temp[ii], pgas[ii],  rho[ii], nne[ii]);
    eos.store_partial_pressures(ndep, ii, eos.xna, eos.xne);
  }

}

std::vector<double> expand_transition_region2(int const N, double tr_loc, double const tr_amp, int const dtr_in)
{

  //constexpr static const int dtr_in = 7;

  // --- Init output array with ones --- //
  std::vector<double> tr(N, 1.0);
  if(tr_loc == 0) return tr;

  // --- check that TR is within limits --- //
  tr_loc = std::min<double>(tr_loc, N-1);
  int    const dtr = std::min<int>(dtr_in, int(tr_loc));
  double const i0  = std::max(0.0, tr_loc-dtr);
  double const i1  = tr_loc;

  // --- define linear integration limits in index numbers --- //
  int const ii0 = int(floor(i0)+0.1);
  int const ii1 = int(ceil(i1)+0.1);
  int const nn  = ii1-ii0;

  // --- Array with the desired output x-values --- //
  std::vector<double> xx(nn, 0.0);
  for(int ii=0; ii<nn; ++ii)
    xx[ii] = ii+ii0;

  // --- Fill the left-hand side of the atmos with TR value --- //
  for(int ii=0; ii<ii0; ++ii)
    tr[ii] = tr_amp;

  // --- Interpolate gradient --- //
  double const x[2] = {i0, i1};
  double const y[2] = {tr_amp, 1.0};
  linpol(2, x, y, nn, &xx[0], &tr[ii0], false);
  
  return tr;
}

std::vector<double> expand_transition_region(int const N, double tr_loc, double const tr_amp, int const dtr_in)
{

  // --- Init output array with ones --- //
  std::vector<double> tr(N, 1.0);
  if((tr_loc == 0) || (tr_amp <= 1.0)) return tr;

  // --- check that TR is within limits --- //
  tr_loc = std::min<double>(tr_loc, N-1);
  int    const dtr = std::min<int>(dtr_in, int(tr_loc));
  double const i0  = std::max(0.0, tr_loc-dtr);
  double const i1  = tr_loc;

  // --- define calculation limits in index numbers --- //
  int const ii0 = int(floor(i0)+0.1);
  int const ii1 = int(ceil(i1)+0.1);
  double const dnn = i1-i0;

  // --- Fill the left-hand side of the atmos with TR value --- //
  for(int ii=0; ii<ii0; ++ii)
    tr[ii] = tr_amp;

  // --- Take the tanh decay --- //
  double imax = 0, imin = 0;
  
  for(int ii=ii0; ii<=ii1; ++ii){
    double const x = ((i1-double(ii)) / dnn - 0.5) * 1.4*  phyc::PI;
    double const y =  tanh(x);
    tr[ii] = y;
    
    imax = std::max(imax, y);
    imin = std::min(imin, y);
  }

  // --- And scale it properly --- //
  double const scl = (tr_amp - 1.0) / (imax-imin);
  for(int ii=ii0; ii<=ii1; ++ii){
    tr[ii] = (tr[ii] - imin) * scl + 1.0;
  }
  
  return tr;
}



void mdepth::expand(nodes_t &n, double *p, int interpol, int mtype){

  double *dep = NULL;
  if(n.depth_t == 0) dep = ltau;
  else               dep = cmass;

  
  if(mtype == 0){ // Nodes are the actual value of the model
    
    if(n.toinv[0]){
      // --- expand nodes --- //
      int len = (int)n.temp.size();
      nodes2depth(len, &n.temp[0], &p[n.temp_off], ndep, dep, &temp[0], interpol, false);

    }
    
    if(n.toinv[1]){
      int len = (int)n.v.size();
      nodes2depth(len, &n.v[0], &p[n.v_off], ndep, dep, &v[0], interpol, false);
    }
    
    if(n.toinv[2]){
      int len = (int)n.vturb.size();
      nodes2depth(len, &n.vturb[0], &p[n.vturb_off], ndep, dep, &vturb[0], interpol, false);
    }
    
    if(n.toinv[3]){
      int len = (int)n.bl.size();
      nodes2depth(len, &n.bl[0], &p[n.bl_off], ndep, dep, &bl[0], interpol, false);
    }
    
    if(n.toinv[4]){
      int len = (int)n.bh.size();
      nodes2depth(len, &n.bh[0], &p[n.bh_off], ndep, dep, &bh[0], interpol, false);
    }
    
    if(n.toinv[5]){
      int len = (int)n.azi.size();
      nodes2depth(len, &n.azi[0], &p[n.azi_off], ndep, dep, &azi[0], interpol, false);
    }

    
    if(n.toinv[6]){
      if     (n.bound == 1) pgas[0] = bound_val*p[n.pgas_off];
      else if(n.bound == 2) rho[0]  = bound_val*p[n.pgas_off];
      else if(n.bound == 3) nne[0]  = bound_val*p[n.pgas_off];
      else pgas[0] =  boundary_pgas_default*p[n.pgas_off];
    }



    // --- apply TR parameterization --- //
    if(n.toinv[7] > 0){
      double const tr_loc_ext = std::max(p[n.tr_off], 0.0);
      double const tr_amp_ext = std::max<double>(p[n.tr_off+1], 1.0);
      fprintf(stderr, "loc=%f, amp=%f\n", tr_loc_ext, tr_amp_ext);
      
      std::vector<double> const tr(std::move(expand_transition_region2(ndep, tr_loc_ext, tr_amp_ext, n.fit_tr)));
      
      for(int ii = 0; ii<ndep; ++ii)
	temp[ii] *= tr[ii];


      tr_amp = tr_amp_ext;
      tr_loc = tr_loc_ext;
    }else{
      tr_amp = 1.0;
      tr_loc = 0.0;
    }
    
    
    
  }else{ // NICOLE/SIR behaviour: node is a correction

    double *tmp = new double [ndep];
    mdepth *r = this->ref_m;
    
    if(n.toinv[0]){
      int len = (int)n.temp.size();
      nodes2depth(len, &n.temp[0], &p[n.temp_off], ndep, dep, tmp, interpol, true);
      for(int ii=0;ii<ndep;ii++){
	temp[ii] = r->temp[ii] + tmp[ii];
      }
    }

    
    if(n.toinv[1]){
      int len = (int)n.v.size();
      nodes2depth(len, &n.v[0], &p[n.v_off], ndep, dep, tmp, interpol, false);
      for(int ii=0;ii<ndep;ii++) v[ii] = r->v[ii] + tmp[ii];
    }

    if(n.toinv[2]){
      int len = (int)n.vturb.size();
      nodes2depth(len, &n.vturb[0], &p[n.vturb_off], ndep, dep, tmp, interpol, false);
      for(int ii=0;ii<ndep;ii++) vturb[ii] = r->vturb[ii] + tmp[ii];
    }
    
    if(n.toinv[3]){
      int len = (int)n.bl.size();
      nodes2depth(len, &n.bl[0], &p[n.bl_off], ndep, dep, tmp, interpol, false);
      for(int ii=0;ii<ndep;ii++) bl[ii] = r->bl[ii]  + tmp[ii];
    }     

    if(n.toinv[4]){
      int len = (int)n.bh.size();
      nodes2depth(len, &n.bh[0], &p[n.bh_off], ndep, dep, tmp, interpol, false);
      for(int ii=0;ii<ndep;ii++) bh[ii] = r->bh[ii] + tmp[ii];
    }

    if(n.toinv[5]){
      int len = (int)n.azi.size();
      nodes2depth(len, &n.azi[0], &p[n.azi_off], ndep, dep, tmp, interpol, false);
      for(int ii=0;ii<ndep;ii++) azi[ii] = r->azi[ii] + tmp[ii];
    }


    if(n.toinv[6]){
      if     (n.bound == 1) pgas[0] = r->pgas[0]*(1.0+p[n.pgas_off]);
      else if(n.bound == 2) rho[0]  = r->rho[0]*(1.0+p[n.pgas_off]);
      else if(n.bound == 3) nne[0]  = r->nne[0]*(1.0+p[n.pgas_off]);
      else pgas[0] =  boundary_pgas_default*p[n.pgas_off];
    }
    //fprintf(stderr,"Pgas=%e %e\n", pgas[0], bound_val);
    
    r = NULL;
    delete [] tmp;
  }
  


  
  return;
}

void mdepth::fill_densities(eoswrap &eos, int keep_nne, int k0, int k1){

  /* --- which scale do we have? --- */
  
  int bound = -1;
  double spgas = 0.0, srho = 0.0, snne = 0.0;

  for(int kk = 0; kk<ndep; kk++){
    spgas += pgas[kk];
    srho  += rho[kk];
    snne  += nne[kk];
  }

  if(spgas > 0.0) bound = 0;
  else if(srho > 0.0) bound = 1;
  else if(snne > 0.0) bound = 3;
  
  /* --- Fill the density arrays, the partial pressures are stored internally
     inside fill_densities --- */

  int ndep1 = k1-k0+1;
  eos.fill_densities(ndep1, &temp[k0], &pgas[k0], &rho[k0], &pel[k0], &nne[k0], bound,  keep_nne, 1.e-5);
  
  

  /*

  double sz = 0.0, stau = 0.0, scm = 0.0;
  for(int kk = 0; kk<ndep; kk++){
    sz += fabs(z[kk]);
    stau += fabs(ltau[kk]);
    scm += fabs(cmass[kk]);
  }

  bound = -1;
  if(stau > 0.0)     bound = 0;
  else if(sz > 0.0)  bound = 1;
  else if(scm > 0.0) bound = 2;


  
  getScales(eos, bound);
  
  */
}
void mdepth::getScales(eoswrap &eos, int bound){

  vector<double> kappa;
  kappa.resize(ndep);
  int nw = 1;
  double wav = 5000.0, scat = 0.0;
  vector<float> frac, part;
  float na=0, ne=0;

  
  /* --- get cont opac --- */
  
  for(int k=0; k<ndep; k++){
    eos.read_partial_pressures(k, frac, part, na, ne);
    eos.contOpacity(temp[k], nw,  &wav, &kappa[k],
			      &scat, frac, na, ne);
  }
  
  
  if(bound == 0){  /* --- if we know tau --- */


    tau[0] = pow(10.0, ltau[0]);
    cmass[0] = 0.0;//(tau[0] / kappa[0]) * rho[0];
    z[0] = 0.0;
    
    for(int k = 1; k < ndep; k++){
      tau[k] = pow(10.0, ltau[k]);
      z[k] = z[k-1] - 2.0 * (tau[k] - tau[k-1]) / (kappa[k] + kappa[k-1]);
      cmass[k] = cmass[k-1] + 0.5*(rho[k-1] + rho[k]) * (z[k-1] - z[k]);
    } // k
    double cmass_off = exp(2.0 * log(cmass[1]) - log(cmass[2]));
    for(int k = 0; k < ndep; k++) cmass[k] = log10(cmass[k]+cmass_off);

  }else if(bound == 1){  /* --- if we know Z --- */

    eos.read_partial_pressures(0, frac, part, na, ne);
    tau[0] = 0;//0.5 * kappa[0] * (z[0] - z[1]);
    cmass[0] = 0;//(na + ne) * (phyc::BK * temp[0] / eos.gravity);
    ltau[0] = log10(tau[0]);

    for(int k = 1; k < ndep; k++){
      tau[k] = tau[k-1] + 0.5 * (kappa[k-1] + kappa[k]) * (z[k-1] - z[k]);
      cmass[k] = cmass[k-1] + 0.5 * (rho[k-1] + rho[k]) * (z[k-1] - z[k]);
    }
    
    double cmass_off = exp(2.0 * log(cmass[1]) - log(cmass[2]));
    double tau_off   =  exp(2.0 * log(tau[1]) - log(tau[2]));
    
    for(int k = 0; k < ndep; k++){
      tau[k] += tau_off;
      cmass[k] = log10(cmass[k]+cmass_off);
      ltau[k] = log10(tau[k]);
    }

    
  }else if(bound ==2){ /* --- If we know cmass --- */
    
    z[0] = 0.0;
    tau[0] = 0.0; //kappa[0]/rho[0] * cmass[0];
    
    
    for(int k = 1; k < ndep; k++){
      z[k] = z[k-1] - 2.0 * (cmass[k] - cmass[k-1]) / (rho[k-1] + rho[k]);
      tau[k] = tau[k-1] + 0.5 * (kappa[k-1] + kappa[k]) * (z[k-1] - z[k]);
    }
    
    /* --- Extrapolate tau at the top --- */
    
    double toff = exp(2.0 * log(tau[1]) - log(tau[2]));
    
    for(int k = 0; k < ndep; k++){
      tau[k] = log10(tau[k]+toff);
    } 
  }

  
}


void mdepth::fixBoundary(int boundary, eoswrap &eos){
  
  if(boundary == 0) pgas[0] = boundary_pgas_default;
  else if(boundary == 1) nne[0] = eos.nne_from_T_Pg(temp[0], pgas[0],  rho[0]);
  else if(boundary == 2) nne[0] = eos.nne_from_T_rho(temp[0], pgas[0], rho[0]);
  else if(boundary == 3) rho[0] = eos.rho_from_T_nne(temp[0], pgas[0], nne[0]);
  else                   rho[0] = eos.rho_from_T_pel(temp[0], pgas[0], pel[0]);
  
}

void mdepth::getPressureScale(int depth_t, int boundary, eoswrap &eos){

  /* --- If pgas was not given, convert rho or nne or pel to pgas --- */
  
  fixBoundary(boundary, eos);

  
  /* --- Solve hydrostatic eq. --- */
  
  hydrostatic(eos, depth_t);



  // --- testing --- //

}


void mdepthall::setsize(int ny, int nx, int ndep, bool verbose){
  std::string inam = "mdepthall::setsize: ";

  std::vector<int> dims = {ny, nx, ndep};
  if(verbose) std::cout << inam << "["<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"]"<<std::endl;
		 

  cub.set({ny, nx, 12, ndep});
}


void mdepthall::compress(int n, double *x, double *y, int nn, double *xx, double *yy){
  if(nn == 0) return;
  else if(nn == 1){
    double tmp = 0.0;
    for(int zz = 0; zz<n;zz++) tmp += y[zz];
    yy[0] = tmp / (double)n;
    return;
  } else {
    linpol<double, double>(n, x, y, nn, xx, yy, false);
    return;
  }
}
void mdepthall::compress(int n, float *x, float *y, int nn, double *xx, double *yy){
  if(nn == 0) return;
  else if(nn == 1){
    double tmp = 0.0;
    for(int zz = 0; zz<n;zz++) tmp += y[zz];
    yy[0] = tmp / (double)n;
    return;
  } else {
    linpol<float, double>(n, x, y, nn, xx, yy, true);
    return;
  }
}

void mdepthall::model_parameters(mat<double> &tmp, nodes_t &n, int nt){

  int nnodes = n.nnodes;
  int nx = temp.size(1);
  int ny = temp.size(0);
  int ndep = temp.size(2);
    
  tmp.set({ny, nx, nnodes});
    
  for(int yy = 0; yy < ny; yy++)
    for(int xx = 0; xx < nx; xx++){

      int k = 0;
      int nn = 0;
	
      // Temp
      nn = (int)n.temp.size();
      compress(ndep, &ltau(yy,xx,0), &temp(yy,xx,0), nn, &n.temp[0], &tmp(yy,xx,k));
      k += nn;

      // v_los
      nn = (int)n.v.size();
      compress(ndep, &ltau(yy,xx,0), &v(yy,xx,0), nn, &n.v[0], &tmp(yy,xx,k));
      k += nn;

      // vturb
      nn = (int)n.vturb.size();
      compress(ndep, &ltau(yy,xx,0), &vturb(yy,xx,0), nn, &n.vturb[0], &tmp(yy,xx,k));
      k += nn;

      // Bl
      nn = (int)n.bl.size();
      compress(ndep, &ltau(yy,xx,0), &bl(yy,xx,0), nn, &n.bl[0], &tmp(yy,xx,k));
      k += nn;
	
      // Bh
      nn = (int)n.bh.size();
      compress(ndep, &ltau(yy,xx,0), &bh(yy,xx,0), nn, &n.bh[0], &tmp(yy,xx,k));
      k += nn;
      
      // azi
      nn = (int)n.azi.size();
      compress(ndep, &ltau(yy,xx,0), &azi(yy,xx,0), nn, &n.azi[0], &tmp(yy,xx,k));
      k += nn;

      //	for(int zz = 0; zz < ndep; zz++) ivar()
    }

}

void mdepthall::model_parameters2( mat<double> &tmp, nodes_t &n, int nt){

  int const nnodes = n.nnodes;
  int const nx = cub.size(1);
  int const ny = cub.size(0);
  int const ndep = cub.size(3);
    
  tmp.set({ny, nx, nnodes});

  int idx = ((n.depth_t == 0)?9:11);
  
  
  for(int yy = 0; yy < ny; yy++)
    for(int xx = 0; xx < nx; xx++){

      int k = 0;
      int nn = 0;


      // Transition region amplification
      if(n.toinv[7]>0){
	tmp(yy,xx,n.tr_off) = std::max<double>(tr_loc(yy,xx),4); // Init at the fourth grid cell
	tmp(yy,xx,n.tr_off+1) = std::max<double>(tr_amp(yy,xx),1.0);

	// --- correct the temperature for the old tr amp before inferring the nodes --- //
	std::vector<double> const tr = expand_transition_region(ndep, tr_loc(yy,xx), tr_amp(yy,xx), tr_N(yy,xx));

	
	double* const __restrict__ Tg = &cub(yy,xx,0,0);
	for(int ii=0; ii<ndep; ++ii)
	  Tg[ii] /= tr[ii];

	// -- Take the new grid ---//
	tr_N(yy,xx) = n.fit_tr;
      }
      
      // Temp
      if(n.toinv[0]){
	nn = (int)n.temp.size();
	compress(ndep, &cub(yy,xx,idx,0), &cub(yy,xx,0,0), nn, &n.temp[0], &tmp(yy,xx,k));
	k += nn;
      }

      // v_los
      if(n.toinv[1]){
	nn = (int)n.v.size();
	compress(ndep, &cub(yy,xx,idx,0), &cub(yy,xx,1,0), nn, &n.v[0], &tmp(yy,xx,k));
	k += nn;
      }
      
      // vturb
      if(n.toinv[2]){
	nn = (int)n.vturb.size();
	compress(ndep, &cub(yy,xx,idx,0), &cub(yy,xx,2,0), nn, &n.vturb[0], &tmp(yy,xx,k));
	k += nn;
      }

      // Blong
      if(n.toinv[3]){
	nn = (int)n.bl.size();
	compress(ndep, &cub(yy,xx,idx,0), &cub(yy,xx,3,0), nn, &n.bl[0], &tmp(yy,xx,k));
	k += nn;
      }
      
      // Bhor
      if(n.toinv[4]){
	nn = (int)n.bh.size();
	compress(ndep, &cub(yy,xx,idx,0), &cub(yy,xx,4,0), nn, &n.bh[0], &tmp(yy,xx,k));
	k += nn;
      }
      
      // azi
      if(n.toinv[5]){
	nn = (int)n.azi.size();
	compress(ndep, &cub(yy,xx,idx,0), &cub(yy,xx,5,0), nn, &n.azi[0], &tmp(yy,xx,k));
	k += nn;
      }

      // Pgas boundary
      
      if(n.toinv[6]>0){
	tmp(yy,xx,k++) = 1.0;
      }


      
      
    }

}


int mdepthall::read_model2(iput_t const& input, std::string &filename, int tstep, bool require_tau){
  io ifile(filename, netCDF::NcFile::read);
  std::string inam = "mdepthall::read_model: ";
  int idep, bound = 0;
  mat<double> tmp;
    
  /* --- get dimensions --- */
  std::vector<int> dims = ifile.dimSize("temp");
  int ndims = (int)dims.size();

  
  if(ndims == 4) dims.erase(dims.begin()); // if time, remove first element
  else if(ndims != 3) {
    std::cout << inam << "ERROR, ndims must be 3 or 4: [(nt), ny, nx, ndep], but is "<<ndims<<std::endl;
  }     

  /* --- Allocate cube --- */
  cub.set({dims[0], dims[1], 13, dims[2]});
  ndep = dims[2];


  
  /* --- read vars, assuming they exists --- */
  if(ifile.is_var_defined("temp")){
    ifile.read_Tstep<double>("temp", tmp, tstep);

    /* --- Copy to consecutive array --- */
    for(int yy = 0; yy< dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++)
	memcpy(&cub(yy,xx,0,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
  }else {
    std::cerr << inam << "ERROR, "<<filename <<" does not contain a temperature array, exiting"<<std::endl;
    exit(0);
  }

  boundary.set({dims[0], dims[1]});


  /* --- Read Vlos --- */
  if(ifile.is_var_defined("vlos")){
    ifile.read_Tstep<double>("vlos", tmp, tstep);
     for(int yy = 0; yy< dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++)
	memcpy(&cub(yy,xx,1,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
  }
  

  
  /* --- Read Vmic --- */
   if(ifile.is_var_defined("vturb")){
    ifile.read_Tstep<double>("vturb", tmp, tstep);
    for(int yy = 0; yy< dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++)
	memcpy(&cub(yy,xx,2,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
   }


   
   /* --- Read blong --- */
   if(ifile.is_var_defined("blong")){
     ifile.read_Tstep<double>("blong", tmp, tstep);
     for(int yy = 0; yy< dims[0]; yy++)
       for(int xx = 0; xx < dims[1]; xx++)
	 memcpy(&cub(yy,xx,3,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
     
   
   


     /* --- Read bhor --- */
     if(ifile.is_var_defined("bhor")){
       ifile.read_Tstep<double>("bhor", tmp, tstep);
       for(int yy = 0; yy< dims[0]; yy++)
	 for(int xx = 0; xx < dims[1]; xx++)
	   memcpy(&cub(yy,xx,4,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
     }


     
   }else{

     /* --- Are B and inc defined? -> Backwards compatibility  --- */
     
     mat<double> tmp1;

     if(ifile.is_var_defined("b") && ifile.is_var_defined("inc")){
       
       ifile.read_Tstep<double>("b", tmp, tstep);
       ifile.read_Tstep<double>("inc", tmp1, tstep);
       
       for(int yy = 0; yy< dims[0]; yy++)
	 for(int xx = 0; xx < dims[1]; xx++)
	   for(int zz = 0; zz< dims[2]; zz++){
	     cub(yy,xx,3,zz) = tmp(yy,xx,zz) * cos(tmp1(yy,xx,zz));
	     cub(yy,xx,4,zz) = tmp(yy,xx,zz) * sin(tmp1(yy,xx,zz));
	   }
     }
   }



   /* --- Read azi --- */
   if(ifile.is_var_defined("azi")){
     ifile.read_Tstep<double>("azi", tmp, tstep);
     for(int yy = 0; yy< dims[0]; yy++)
       for(int xx = 0; xx < dims[1]; xx++)
	 memcpy(&cub(yy,xx,5,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
   }
   



   /* --- Read Pgas --- */
   if(ifile.is_var_defined("pgas")){
     ifile.read_Tstep<double>("pgas", tmp, tstep);
     for(int yy=0; yy<dims[0]; yy++) for(int xx = 0;xx<dims[1]; xx++){
	 boundary(yy,xx) = tmp(yy,xx,0);
	 memcpy(&cub(yy,xx,6,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
       }
   }


 
   
   /* --- Read Rho --- */
   if(ifile.is_var_defined("rho")){
     ifile.read_Tstep<double>("rho", tmp, tstep);
     if(bound == 0) {
       for(int yy=0; yy<dims[0]; yy++)
	 for(int xx = 0;xx<dims[1]; xx++){
	   boundary(yy,xx) = tmp(yy,xx,0);
	 }
     }
     
     for(int yy=0; yy<dims[0]; yy++)
       for(int xx = 0;xx<dims[1]; xx++)
	 memcpy(&cub(yy,xx,7,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
   }



   
   /* --- Read nne --- */
   if(ifile.is_var_defined("nne")){
     ifile.read_Tstep<double>("nne", tmp, tstep);
     if(bound == 0) {
       for(int yy=0; yy<dims[0]; yy++)
	 for(int xx = 0;xx<dims[1]; xx++){
	   boundary(yy,xx) = tmp(yy,xx,0);
	 }
     }
     
     for(int yy=0; yy<dims[0]; yy++)
       for(int xx = 0;xx<dims[1]; xx++)
	 memcpy(&cub(yy,xx,8,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
   }


   /* --- Init boundary --- */
   if(fabs(cub(0,0,6,0)) > 0.0) bound = 1;
   else if(fabs(cub(0,0,7,0)) > 0.0) bound = 2;
   else if(fabs(cub(0,0,8,0)) > 0.0) bound = 3;

   if(bound > 0)
     for(int yy=0; yy<dims[0]; yy++)
       for(int xx = 0;xx<dims[1]; xx++)
	 boundary(yy,xx) = cub(yy,xx,5+bound, 0);
   cerr<<"mdepthall::read_model2: Bound -> "<<bound<<endl;
   


   /* --- Read LTAU500 --- */
   bool set_ltau = false;
   if(ifile.is_var_defined("ltau500")){
    ifile.read_Tstep<double>("ltau500", tmp, tstep);
    set_ltau = true;
    for(int yy=0; yy<dims[0]; yy++)
      for(int xx = 0;xx<dims[1]; xx++)
	memcpy(&cub(yy,xx,9,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
   }



   /* --- Read Z --- */
   bool set_z = false;
   if(ifile.is_var_defined("z")){
     ifile.read_Tstep<double>("z", tmp, tstep);
     set_z = true;
     if(tmp.ndims() == 1){
       cerr<< inam <<"replicating z-scale in all pixels"<<endl;
       for(int yy = 0; yy<dims[0];yy++)
	 for(int xx = 0; xx< dims[1]; xx++)
	   for(int zz = 0; zz<dims[2]; zz++)
	     cub(yy,xx,10,zz) = tmp.d[zz];
     }else{
       for(int yy=0; yy<dims[0]; yy++)
	 for(int xx = 0;xx<dims[1]; xx++)
	   memcpy(&cub(yy,xx,10,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
       
     }
     
   }
   
   /* -- Read cmass? --- */
   
   bool set_cmass = false;
   if(ifile.is_var_defined("cmass")){
     set_cmass = true;
     ifile.read_Tstep<double>("cmass", tmp, tstep);
     for(int yy=0; yy<dims[0]; yy++)
       for(int xx = 0;xx<dims[1]; xx++)
	 memcpy(&cub(yy,xx,11,0), &tmp(yy,xx,0), dims[2]*sizeof(double));
   }

   /* --- Tr amplification factor --- */
   if(ifile.is_var_defined("transition_region_scale")){
     ifile.read_Tstep<double>("transition_region_scale", tr_amp, tstep);
   }else{
     tr_amp.set({dims[0], dims[1]});
     long const nTot = long(dims[0]) * long(dims[1]);
     for(long ii=0; ii<nTot; ++ii)
       tr_amp(ii) = 1.0;
   }
   
   /* --- Tr location --- */
   if(ifile.is_var_defined("transition_region_loc")){
     ifile.read_Tstep<double>("transition_region_loc", tr_loc, tstep);
   }else{
     tr_loc.set({dims[0], dims[1]});
     long const nTot = long(dims[0]) * long(dims[1]);
     for(long ii=0; ii<nTot; ++ii)
       tr_loc(ii) = 30; 
   }
   
   /* --- Tr N --- */
   if(ifile.is_var_defined("transition_region_nGrid")){
     ifile.read_Tstep<int>("transition_region_nGrid", tr_N, tstep);
   }else{
     tr_N.set({dims[0], dims[1]});
     long const nTot = long(dims[0]) * long(dims[1]);
     for(long ii=0; ii<nTot; ++ii)
       tr_N(ii) = input.fit_tr;
   }

   
   if(set_ltau == false && set_cmass == false){
     std::cerr << inam << "ERROR, "<<filename
	      <<" does not contain a depth-scale [ltau500] or [cmass], exiting"
	      <<std::endl;
    exit(0);
  }
   

   return bound;
}
/*
void mdepthall::convertBoundary(int bound, bool verbose){

  string inam = "depthall::convertBoundary: ";
  
  vector<int> dims = boundary.getdims();
  
  if(verbose) cout << inam << "converting boundary to Pgas if needed ... ";
  switch(bound)
    {
    case(2):
      for(int yy=0;yy<dims[0];yy++){
	cerr << yy <<" ";
	for(int xx=0;xx<dims[1];xx++){
	  boundary(yy,xx) = eos.nne_from_T_rho(temp(yy,xx,0), pgas(yy,xx,0), rho(yy,xx,0));
	}
      }
      break;
    case(3):
      for(int yy=0;yy<dims[0];yy++)
	for(int xx=0;xx<dims[1];xx++)
	  boundary(yy,xx) = eos.nne_from_T_rho(temp(yy,xx,0), pgas(yy,xx,0), nne(yy,xx,0));
      break;
    case(4):
      for(int yy=0;yy<dims[0];yy++)
	for(int xx=0;xx<dims[1];xx++)
	  boundary(yy,xx) = eos.nne_from_T_rho(temp(yy,xx,0), pgas(yy,xx,0), pel(yy,xx,0));
      break;
    defaul:
      break;
    }
  
  if(verbose) cout << endl;
  
}
*/
void mdepthall::expand(int n, double *x, double *y, int nn, double *xx, double *yy, int interpolation){

  if     (n == 1)                for(int kk=0;kk<nn;kk++) yy[kk] = y[0];
  else if((n == 2))              linpol(n, x, y, nn, xx, yy, true);
  else if(n >= 3){
    if(interpolation == 0)       linpol(n, x, y, nn, xx, yy, true);
    else if(interpolation == 1) bezpol2(n, x, y, nn, xx, yy, true);
    else if(interpolation == 2) hermpol(n, x, y, nn, xx, yy, true);
    else if(interpolation == 3){
      if(n >= 3) vlint(n, x, y, nn, xx, yy);
      else     linpol(n, x, y, nn, xx, yy, true);
    }
  }
  else return;
}


void mdepthall::expandAtmos(nodes_t &n, mat<double> &pars, int interpolation){

  int ny = pars.size(0);
  int nx = pars.size(1);

  int idx = ((n.depth_t == 0)?9:11);
  
  for(int yy = 0; yy< ny; yy++) for(int xx=0;xx<nx;xx++){

      /* --- temperature --- */
      if(n.toinv[0]){
	int len = (int)n.temp.size();
	expand(len, &n.temp[0], &pars(yy,xx,n.temp_off), ndep, &cub(yy,xx,idx,0), &cub(yy,xx,0,0), interpolation);

	
	
      }


      /* --- vlos --- */
      if(n.toinv[1]){
	int len = (int)n.v.size();
	expand(len, &n.v[0], &pars(yy,xx,n.v_off), ndep, &cub(yy,xx,idx,0), &cub(yy,xx,1,0), interpolation);
      }


      /* --- vturb --- */
      if(n.toinv[2]){
	int len = (int)n.vturb.size();
	expand(len, &n.vturb[0], &pars(yy,xx,n.vturb_off), ndep, &cub(yy,xx,idx,0), &cub(yy,xx,2,0), interpolation);
      }

      
      /* --- Blong --- */
      if(n.toinv[3]){
	int len = (int)n.bl.size();
	expand(len, &n.bl[0], &pars(yy,xx,n.bl_off), ndep, &cub(yy,xx,idx,0), &cub(yy,xx,3,0), interpolation);
      }


      /* --- Bhor --- */
      if(n.toinv[4]){
	int len =(int)n.bh.size();
	expand(len, &n.bh[0], &pars(yy,xx,n.bh_off), ndep, &cub(yy,xx,idx,0), &cub(yy,xx,4,0), interpolation);
      }



      /* --- Azi --- */
      if(n.toinv[5]){
	int len = (int)n.azi.size();
	expand(len, &n.azi[0], &pars(yy,xx,n.azi_off), ndep, &cub(yy,xx,idx,0), &cub(yy,xx,5,0), interpolation);
      }
      
    } // xx & yy
  
  
}


void mdepthall::write_model(string &filename, int tstep){

  static bool firsttime = true;

  /* --- init output file ont he first call --- */
  static io ofile(filename, netCDF::NcFile::replace);


  /* --- Init vars & dims if firsttime --- */
  if(firsttime){
    /* --- Dims --- */
    vector<int> dims = temp.getdims();
    dims.insert(dims.begin(), 0);
    ofile.initDim({"time","y", "x", "ndep"}, dims);

    /* --- vars -- */
     ofile.initVar<float>(string("temp"),    {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("vlos"),    {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("vturb"),   {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("blong"),       {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("bhor"),     {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("azi"),     {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("ltau500"), {"time","y", "x", "ndep"});
     ofile.initVar<float>(string("pgas"),    {"time","y", "x", "ndep"});

     firsttime = false;
    
  }


  /* --- write time stamp --- */
  ofile.write_Tstep(string("temp"),    temp,  tstep);
  ofile.write_Tstep(string("vlos"),    v,     tstep);
  ofile.write_Tstep(string("vturb"),   vturb, tstep);
  ofile.write_Tstep(string("blong"),    bl,     tstep);
  ofile.write_Tstep(string("bhor"),     bh,   tstep);
  ofile.write_Tstep(string("azi"),     azi,   tstep);
  ofile.write_Tstep(string("ltau500"), ltau,  tstep);
  ofile.write_Tstep(string("pgas"),    pgas,  tstep);
}



void mdepthall::write_model2(iput_t const& input, string &filename, int tstep){

  static bool firsttime = true;

  /* --- init output file ont he first call --- */
  
  static io ofile(filename, netCDF::NcFile::replace);

  
  /* --- Dims --- */
  vector<int> cdims = cub.getdims();
  vector<int> dims = {0, cdims[0], cdims[1], cdims[3]};


  
  /* --- Init vars & dims if firsttime --- */
  if(firsttime){

    ofile.initDim({"time","y", "x", "ndep"}, dims);

    
    /* --- vars -- */
    
    ofile.initVar<float>(string("vlos"),    {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("temp"),    {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("vturb"),   {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("blong"),       {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("bhor"),     {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("azi"),     {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("ltau500"), {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("z"),       {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("pgas"),    {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("rho"),     {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("nne"),     {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("cmass"),     {"time","y", "x", "ndep"});
    ofile.initVar<float>(string("transition_region_loc"),     {"time","y", "x"});
    ofile.initVar<float>(string("transition_region_scale"),     {"time","y", "x"});
    ofile.initVar<int>(string("transition_region_nGrid"),     {"time","y", "x"});

    firsttime = false;
    
  }
  
  {
  mat<double> tmp((vector<int>){dims[1], dims[2], dims[3]});
  
  /* --- write time step --- */
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,0,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("temp"),    tmp,  tstep);

  
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,1,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("vlos"),   tmp,  tstep);

  
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,2,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("vturb"),   tmp, tstep);

  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,3,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("blong"),       tmp,     tstep);

  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,4,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("bhor"),     tmp,   tstep);

  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,5,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("azi"),     tmp,   tstep);

  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,9,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("ltau500"), tmp,  tstep);

  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,10,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("z"), tmp,  tstep);

  
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,6,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("pgas"),    tmp,  tstep);
  
  
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,7,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("rho"),    tmp,  tstep);
  
  
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,8,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("nne"),    tmp,  tstep);

  
  for(int yy=0; yy<dims[1]; yy++)
    for(int xx = 0;xx<dims[2];xx++)
      memcpy(&tmp(yy,xx,0), &cub(yy,xx,11,0), dims[3]*sizeof(double));
  ofile.write_Tstep<double>(string("cmass"),    tmp,  tstep);
  }

  {
    ofile.write_Tstep<double>(string("transition_region_loc"),   tr_loc,  tstep);
    ofile.write_Tstep<double>(string("transition_region_scale"),   tr_amp,  tstep);
    ofile.write_Tstep<int>(string("transition_region_nGrid"),   tr_N,  tstep);
  }

  
}


void Ncompress(int const n,  double*  x,  double*  y, int const nn,  double*  xx, double*  yy)
{
  if(nn == 0) return;
  else if(nn == 1){
    double tmp = 0.0;
    for(int zz = 0; zz<n;zz++) tmp += y[zz];
    yy[0] = tmp / (double)n;
    return;
  } else {
    linpol<double, double>(n, x, y, nn, xx, yy, false);
    return;
  }
}

double* mdepth::compress(nodes_t &n)const
{
  int const nnodes = n.nnodes;
  //int const ndep = temp.size(2);
  double *tmp = new double [nnodes]();
  
  int k = 0;
  int nn = 0;
	
  // Temp
  nn = (int)n.temp.size();
  Ncompress(ndep, ltau, temp, nn, &n.temp[0], &tmp[k]);
  k += nn;
  
  // v_los
  nn = (int)n.v.size();
  Ncompress(ndep, ltau, v, nn, &n.v[0], &tmp[k]);
  k += nn;
  
  // vturb
  nn = (int)n.vturb.size();
  Ncompress(ndep, ltau, vturb, nn, &n.vturb[0], &tmp[k]);
  k += nn;
  
  // Bl
  nn = (int)n.bl.size();
  Ncompress(ndep, ltau, bl, nn, &n.bl[0], &tmp[k]);
  k += nn;
  
  // Bh
  nn = (int)n.bh.size();
  Ncompress(ndep, ltau, bh, nn, &n.bh[0], &tmp[k]);
  k += nn;
  
  // azi
  nn = (int)n.azi.size();
  Ncompress(ndep, ltau, azi, nn, &n.azi[0], &tmp[k]);
  k += nn;
  
  
  return tmp;
}
