/*
  Sparse inversion class
  Authors: Jaime de la Cruz Rodriguez (ISP-SU 2014) & Andres Asensio-Ramos (IAC-2014)
  Dependencies: cmemt.h, wavelet.{h,cc}, GSL (library), instrument.h
 */
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <cstring>
#include <cstdio>
#include <netcdf>
//#include <omp.h>
#include <gsl/gsl_sort_double.h>
#include <sys/time.h>
#include "sparse.h"
#include "wavelet.h"
#include "comm.h"
#include "io.h"
#include "clte.h"
#include "interpol.h"
//
using namespace std;
using namespace netCDF;

/* ---- Some aux functions for the line search functions --- */
void mysort4(double &a, double &b, double &c, double &d, double &sa, double &sb, double &sc, double &sd){

  double tmp[4] = {a, b, c, d};
  double tmp1[4] = {sa, sb, sc, sd};
  
  unsigned long idx[4] = {0,1,2,3};
  gsl_sort_index(&idx[0], &tmp[0], 1, 4);

  a = tmp[idx[0]];
  b = tmp[idx[1]];
  c = tmp[idx[2]];
  d = tmp[idx[3]];

  sa = tmp1[idx[0]];
  sb = tmp1[idx[1]];
  sc = tmp1[idx[2]];
  sd = tmp1[idx[3]];

}

void mysort3(double &a, double &b, double &c, double &sa, double &sb, double &sc){

  double tmp[3] = {a, b, c};
  double tmp1[3] = {sa, sb, sc};
  
  unsigned long idx[3] = {0,1,2};
  gsl_sort_index(&idx[0], &tmp[0], 1, 3);

  a = tmp[idx[0]];
  b = tmp[idx[1]];
  c = tmp[idx[2]];

  sa = tmp1[idx[0]];
  sb = tmp1[idx[1]];
  sc = tmp1[idx[2]];

}


void MeSt(double d, double e, double f, double &me, double &st){
  me = (d+e+f)/3.0;
  st = sqrt(0.5*(sqr<double>(d-me) + sqr<double>(e-me) + sqr<double>(f-me)));
  
}



//
void sparse2d::evalModel_mpi(mat<double> &x, mat<double> &syn, mat<double> &dsyn, mdepthall_t &m, int compute_gradient){
  

  unsigned long ntot = dims[0] * dims[1];
  int ncom = int(std::floor(ntot / double(iput.npack)));
  if((unsigned long)(ncom * iput.npack) != ntot) ncom++;
  //
  int nprocs = iput.nprocs;
  int iproc = 0;
  unsigned long ipix = 0;
  int tocom = ncom;
  mat<double> dum; // dummy chi parameter


  // Init slaves
  for(int ss = 1; ss<=min(nprocs-1,ncom); ss++) {
    
    comm_master_pack_data(iput, syn, x, ipix, ss, m, compute_gradient);
    //tocom--;
    }



  // While loop
  int per  = 0;
  int oper  = -1;
  float pno =  100.0 / (float(ntot) - 1.0);
  unsigned long irec = 0;
  //
  if(compute_gradient == 1) cerr << "\rSynth + derivatives -> "<<per<<"% ";
  else  cerr << "\rSynth -> "<<per<<"%                    ";

  
  while(irec < ntot){
    // Receive processed data from any slave (iproc)
    comm_master_unpack_data(iproc, iput, syn, x, dum, irec, dsyn, compute_gradient, m);
    per = irec * pno;

    // Send more data to that same slave (iproc)
    if(ipix < ntot) comm_master_pack_data(iput, syn, x, ipix, iproc, m, compute_gradient);
    
    // Keep count of communications left
    tocom--;
    
    // Printout
    if(per > oper){
      oper = per;
      if(compute_gradient == 1) cerr << "\rSynth + derivatives -> "<<per<<"% ";
      else  cerr << "\rSynth -> "<<per<<"%                    ";
    }
  }

}
void sparse2d::evalModel_serial(mat<double> &x, mat<double> &syn, mat<double> &dsyn,  mdepthall_t &m){
  
  int id = 0;
  int ncopy = dims[2]*dims[3]*sizeof(double);
  int per=0, oper=-1;

  for(int yy = 0; yy < dims[0]; yy++){
    for(int xx = 0; xx < dims[1]; xx++){
      // atm[id]->synth_grad(&x(yy,xx,0), &syn(yy,xx,0,0), &dsyn(yy,xx,0,0,0), iput.dpar);
      //for(int pp = 0;pp<npar;pp++) memcpy(&dsyn[pp](yy,xx,0,0), &tmp[pp][0], ncopy);
    }
    per = yy * 100./double(dims[0]-1);
    if(oper < per){
      oper = per;
      cerr << "\rSynth -> "<<per<<"% ";
    }
  }
}

double sparse2d::meritFunction(mat<double> &obs, mat<double> &x,  mat<double> &dchi, mat<double> &noise,  mdepthall_t &m, int compute_gradient){
  
  string inam = "sparse2d::meritFunction: ";
  bool gradient = false;
  if(compute_gradient > 0) gradient = true;
  {
    // Convert parameters to physical space
    wlt.transformSlices(x, dwt_inverse);

    for(int pp = 0; pp<npar; pp++)  
       for(int yy = 0; yy < dims[0]; yy++)
	 for(int xx = 0; xx < dims[1]; xx++)
	   x(pp,yy,xx) *= scalingParameters[pp];
    
  }
  double chi2 = 0.0;
  {
    //
    // Init arrays to store synthetic spectra and derivatives
    //
    mat<double> syn;
    syn.set(dims);
    //
    mat<double> dsyn;
    double norm = dims[0] * dims[1] * dims[2] * dims[3];

    if(gradient){
      vector<int> ndims = {dims[0], dims[1], npar, dims[2], dims[3]};
      dsyn.set(ndims);
      stime = 0.0;
    }
    
    sparsity.resize(npar,0.0L);

    //
    // Get profiles and numerical response-functions
    //
    transposePars(x, unsigned(1));
    double dum = gettime();
    if(iput.nprocs > 1) evalModel_mpi(x, syn, dsyn, m, compute_gradient);
    else evalModel_serial(x, syn, dsyn, m);
    stime += gettime() - dum;
    transposePars(x, unsigned(0));
    

    //
    // Degrade synthetic profiles and derivatives
    //
    isyn = syn;
    inst->degrade(syn, true, false);
    
    
    if(gradient){
      //cerr << ", degrading spectra and derivatives ... ";
      mat<double> kk(dims);
      double ttt = gettime();
      for(int pp = 0; pp<npar;pp++) {
	
	for(int yy = 0; yy<dims[0];yy++)
	  for(int xx = 0; xx<dims[1]; xx++)
	    memcpy(&kk(yy,xx,0,0), &dsyn(yy,xx,pp,0,0), dims[2]*dims[3]*sizeof(double));
	
	inst->degrade(kk, true, false);
	
	for(int yy = 0; yy<dims[0];yy++)
	  for(int xx = 0; xx<dims[1]; xx++)
	    memcpy(&dsyn(yy,xx,pp,0,0), &kk(yy,xx,0,0), dims[2]*dims[3]*sizeof(double));
      }
      
      cerr << gettime()-ttt<<"s ";
    }
   
    
    //
    // Get chi
    //
    int yy,xx,ww,ss,id,pp;
    // mat<double> slice(dims[0], dims[1]);
    
    //
    for(yy = 0; yy < dims[0]; yy++)
      for( xx = 0; xx < dims[1]; xx++)
	for( ww = 0; ww < dims[2]; ww++)
	  for( ss = 0; ss < dims[3]; ss++){			       
	    syn(yy,xx,ww,ss) -= obs(yy,xx,ww,ss); // REMEMBER THIS!!! ESPECIALLY IN THE NEXT NESTED LOOP!
	    syn(yy,xx,ww,ss) /= noise(ww,ss);     // REMEMBER THIS!!!
	    chi2 += iweight(yy,xx) * sqr(syn(yy,xx,ww,ss)) * 0.5;
	    syn(yy,xx,ww,ss) *= iweight(yy,xx) /  noise(ww,ss); // REMEMBER THIS BELOW!!!
	  }
    
    /* --- Normalize by the number of elements --- */
    chi2 /= norm;
    
    //
    // dchi2 -> gradient of the merit function
    // Changed order of the transform compared to Andres' original implementation:
    // First add all parameters to dchi and then transform to Wavelet space (much faster!)
    // Note that [syn] contains a lot of factors that would be applyed redundantly 
    //
    if(gradient){
      dchi.zero();
#pragma omp parallel default(shared) private(pp,xx,yy,ww,ss,id) num_threads(nthreads)
      {
	// id = omp_get_thread_num();
#pragma omp for 
	for( yy = 0; yy < dims[0]; yy++)
	  for( xx = 0; xx < dims[1]; xx++)
	    for( pp = 0; pp<npar;pp++)
	      for( ww = 0; ww < dims[2]; ww++)
		for( ss = 0; ss < dims[3]; ss++)
		  dchi(pp,yy,xx) += scalingParameters[pp] * dsyn(yy,xx,pp,ww,ss) * syn(yy,xx,ww,ss);
    } // parallel block
      
      /* --- Normalize by the number of elements --- */
      for(auto &it: dchi.d) it /= norm;
      
      wlt.transformSlices(dchi, dwt_forward);
      
    } // if gradient
  } //Chi2 block
    
  

  // Convert back the parameters
  
  wlt.transformSlices(x, dwt_forward);
  
  for(int pp = 0; pp<npar; pp++)      
    for(int yy = 0; yy < dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++)
	x(pp,yy,xx) /= scalingParameters[pp];
      
  
  
  return chi2;
}
//

//
void sparse2d::SparseOptimization(mat<double> &obs, mat<double> &x, mat<double> &noise,  mdepthall_t &m, mat<double> &pwe){
  /*	
    Parameters must be an array with shape (ny, nx, npar)!
  */
  string inam = "sparse2d::SparseOptimization: ";
  


  // Init pixel dependent weights and normalize weights
  
  if(pwe.d.size() == 0){
    std::vector<int> ndims = {dims[0], dims[1]};
    iweight.set(ndims);
    for(int yy=0; yy<dims[0]; yy++) 
      for(int xx=0; xx<dims[1]; xx++)
	for(int ww = 0; ww<dims[2]; ww++) iweight(yy,xx) += obs(yy,xx,ww,0);
    
    //double tmp =  double(iweight.n_elements()) / iweight.sum();
    //for(auto &it: iweight.d) it *= tmp;
    double tmp =   iweight.sum()/double(iweight.n_elements());
    for(auto &it: iweight.d) if(it > 0.0) it = tmp/it;
  } else iweight = pwe;
  
  // Transpose parameters so dimensions are (Npar, Ny, Nx)
  transposePars(x, unsigned(0));


  { // Transform all parameters to Wavelet space
    
    // Wavelet transform
    for(int pp = 0; pp < npar; pp++) 
      for(int yy = 0; yy < dims[0]; yy++)
	for(int xx = 0; xx < dims[1]; xx++)
	  x(pp,yy,xx) /= scalingParameters[pp];
      
      // Wavelet transform
      wlt.transformSlices(x, dwt_forward);
      
      
  } // Wavelet block


  /* --- Init merit function --- */
  double logLike = 1e10, logLikeOld = 1e10, deltaLogLike = 0.0;
  mat<double> dlogLike(x.getdims());
  mat<double> xOld(x.getdims());
  

  /* --- Init Fista algorithm --- */
  mat<double> y = x;
  double t = 1.0;
  mat<double> xbest = x, v = x;

  double t0=0, t1=0, tsum=0;
  struct timeval dum;
  bestLogLike = 1.e11;

  double L = 2.0, tol = 2.e-2;
  bool lastreject = true;
  int naccept = 0;
  int nrej = 0;
  double delta_f = 0.0;
  
  
  /* --- Main loop --- */
  {
    mat<double> dlogLikeOld(dlogLike.getdims());
    cerr<<endl;
    //
    for(unsigned i = 1; i<maxiter; i++){

      /* --- get time --- */
      gettimeofday(&dum, NULL);
      t0 = dum.tv_sec + dum.tv_usec * 1.0E-6;

      xOld = x;
      dlogLikeOld = dlogLike;
      double theta = 2.0 / ((double)i + 1.0);

      for(int ii=0;ii<(int)y.d.size();ii++)
	y.d[ii] = (1.0-theta)*x.d[ii] + theta*v.d[ii];

      /* --- Check that parameters are within valid rages --- */
      checkParameters(y);
      
      /* --- Compute the gradient of y and Chi2[y] ---*/
      logLikeOld = logLike;
      logLike = meritFunction(obs, y, dlogLike, noise, m, 1);
      if(i == 1) logLikeOld = logLike;
      
      /* --- Get norm of the gradient --- */
      double gnorm = getNorm((int)dlogLike.d.size(),&dlogLike.d[0]);
      if(delta_f > 0 || i == 1) L = fabs(iput.init_step/gnorm);
      
      
      /* --- Get step size ---*/
      logLike =  my_getStep(L, obs,  logLike, dlogLike, y, noise, m, 1);
      delta_f = logLike - logLikeOld;



      /* --- Take step and threshold parameter values --- */
      mat<double> xnew = takeStep(y, dlogLike, L, true);

      
      /* --- Update v for next step --- */
      for(int kk=0;kk<(int)v.d.size();kk++)
	v.d[kk] = x.d[kk] + (xnew.d[kk] - x.d[kk]) / theta;

      /* --- Update t --- */
      //double tnew = 0.5 * (1.0 + sqrt(1.0 + 4.0 * t * t));
      
      
      /* --- Test for restarting following O'Donogue & Candes --- */
      double restart = 0.0; //dot product of (y-xnew) * (new - x)
      for(long unsigned ii = 0; ii<x.d.size(); ii++)
	restart += (y.d[ii] - xnew.d[ii]) * (xnew.d[ii] - x.d[ii]);
      string rej;


      deltaLogLike = logLike - logLikeOld;
      logLikeOld = logLike;
      
      /* --- Best Chi^2 or reject? --- */
      if(logLike < bestLogLike) {
	xbest = x;
	bsyn = isyn;
	rej = "  ";
	bestLogLike = logLike;
	lastreject = false;
	naccept++;
	nrej = 0;
      } else{
	rej=" *";
	lastreject = true;
	naccept = 0;
	nrej++;
      }

      x = xnew;

      
      if(restart > 0.0){
	//	t = 1.0;
	y = x;
	v = x;
	cerr <<endl<<"Restarting direction!"<<endl;
	logLikeOld = logLike;
	deltaLogLike = logLike - logLikeOld;
      }
            

      /* --- Printout info --- */
      gettimeofday(&dum, NULL);
      t1 = dum.tv_sec + dum.tv_usec * 1.0E-6;
      tsum += t1-t0;
      double tmean = tsum / double(i+1);
      fprintf(stderr,
	      "\r[%4d]: Chi2=%12.6f, dChi2=%13.6f, restart=%11.4f, step=%e, total_time=%6.1fs, synth_time=%6.1fs, estimatedLeft=%7lds %s\n", 
	      i, logLike, deltaLogLike, 
	      restart, L, t1-t0, stime, (long unsigned)((maxiter-1-i)*tmean + 0.5), rej.c_str() );
      
      //
      chisq.push_back(logLike / double(npix * dims[2] * dims[3]));
      if(nrej >= 5) break;
    } // i < maxiter
    
  } // main loop block


  /* --- Prepare results --- */
  cerr <<endl<< inam<<"Best Chi2="<< bestLogLike <<endl; 
  x = xbest;
  finalTransform(x);
  transposePars(x, unsigned(1));
  obs = bsyn;
}
void sparse2d::finalTransform(mat<double> &x){
  
  //mat<double> temp(dims[0], dims[1]);
  
  wlt.transformSlices(x, dwt_inverse);

  // Convert parameters to physical space
  // for(int pp = 0; pp<npar; pp++){
    
  // getSlice(x, temp, pp);
      
  // wlt.transform(temp, dwt_inverse);
  for(int pp = 0; pp<npar; pp++)
    for(int yy = 0; yy < dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++)
	x(pp,yy,xx) *= scalingParameters[pp];
  //}
// 
}

void sparse2d::thresholdSingle(mat<double> &x, double thr){
  
  if(thres_mod == spt_hard){
    // Init indexes array
    vector<long unsigned> idx;
    idx.resize(x.d.size());
    
    vector<double> tmp = x.d;
    for(auto &it: tmp) it = fabs(it);


    // sort values and get index numbers
    gsl_sort_index(&idx[0], &tmp[0], 1, x.d.size());

    // Set to zero, forcing a % of the values to be zero
    long unsigned idx_max = (long unsigned)((1.0-thr) * x.d.size() + 0.5);
    for(long unsigned ii = 0; ii<idx_max; ii++) x.d[idx[ii]] = 0.0L;
    
    return;

  } else{ // soft-thresholding
    for(auto &it: x.d) max(0.0, (1.0 - thr) / max(fabs(it), 1e-10)) * it;
  }
  
}
//
void sparse2d::threshold(mat<double> &x, double thr){
  vector<mat<double>> tmp;
  tmp.resize(nthreads);
  for(auto &it: tmp) it.set({dims[0], dims[1]});
  

  //
  // Check limits of the model in physical space
  //
  //checkParameters(x);
    
  // For all parameters in the model					
  int par=0, id=0;
#pragma omp parallel default(shared) private(par,id) num_threads(nthreads)
  //id = omp_get_thread_num(); // get thread number
  
#pragma omp for
  for( par = 0; par<npar; par++){
    //getSlice(x, tmp[id], par);
    memcpy(&tmp[id].d[0], &x(par,0,0), (unsigned long)(dims[0]*dims[1])*sizeof(double));
    
    // Threshold values for this parameters
    thresholdSingle(tmp[id], thres);
    
    
    // Get sparsity
    for(auto &it: tmp[id].d) if(fabs(it) < 1.e-15) sparsity[par]++;
    sparsity[par] /= npix;
    
    
    // Copy back
    memcpy( &x(par,0,0), &tmp[id].d[0], (unsigned long)(dims[0]*dims[1])*sizeof(double));

    //    putSlice(x, tmp[id], par);
    
  } // par
  
}
sparse2d::sparse2d(iput_t &input, std::vector<int> &dims1, double threshold, 
		   wavelet_type family, unsigned ord, spthres thresholdMode, unsigned nt1){

  init(input, dims1,  threshold, family,  ord,  thresholdMode,  nt1);
}

void sparse2d::init(iput_t &input, std::vector<int> &dims1, double threshold, 
		   wavelet_type family, unsigned ord, spthres thresholdMode, unsigned nt1){
    
  std::string inam = "sparse::sparse: ";
  
  //
  // Init some variables
  // 
  nthreads = nt1;
  thres_mod = thresholdMode;
  thres = threshold;
  dims = dims1;
  npix = dims[0] * dims[1];
  maxiter = input.max_inv_iter;
  npar = (int)input.nodes.nnodes;
  slsize = dims[0] * dims[1] * sizeof(double);
  
  //
  // Init wavelet class
  //
  {
    std::vector<int> wdims = {dims[0], dims[1]};
    wlt.init(wdims, nthreads, family, ord);
  }
  
  
  //
  // Init instrument
  // 
  if(!(input.instrument.compare(string("hinode")))) inst = new hinode(nthreads, dims);
  else {
    std::cerr << inam << "instrument ["<<input.instrument << "] not implemented yet!" <<std::endl;
    inst = new instrument();
  }

  
  
  //
  // Init atmos and get max/min
  //
  //int kk=0;
  atm.resize(nthreads);
  for(auto &it: atm) {
    
    // if(input.atmos_type == string("cmilne")) it = new cmilne(line, input);
    if(input.atmos_type == string("lte")) it = new clte(input);
    else {
      std::cerr << inam << "ERROR ["<<input.atmos_type << "] not implemented yet!" <<std::endl;
      exit(0);
    }
  }
  
  mmax              = atm[0]->get_max_limits(input.nodes);
  mmin              = atm[0]->get_min_limits(input.nodes);
  scalingParameters = atm[0]->get_scaling(input.nodes);
  LCommon           = atm[0]->get_steps(input.nodes);

  //
  // Init arrays to store parameters and transformed parameters
  //
  {
    std::vector<int> ndims = {npar, npix};
    //  parametersPbP.set(ndims);
    // parametersPbPThresholded.set(ndims);
  }
  {
    std::vector<int> ndims = {npar, dims[0], dims[1]};
    // final.set(ndims);
  }
  
  iput = input;
}

void sparse2d::transposePars(mat<double> &x, unsigned dir){
  
  mat<double> tmp = x;
  
  if(dir == 0){
    x.set({npar, dims[0], dims[1]});
    for(int pp = 0; pp < npar; pp++)
      for(int yy = 0; yy < dims[0]; yy++)
	for(int xx = 0; xx < dims[1]; xx++)
	  x(pp,yy,xx) = tmp(yy,xx,pp);
  } else{  
    x.set({ dims[0], dims[1], npar});
    for(int pp = 0; pp < npar; pp++)
      for(int yy = 0; yy < dims[0]; yy++)
	for(int xx = 0; xx < dims[1]; xx++)
	  x(yy,xx,pp) = tmp(pp,yy,xx);
  }
   
}
void sparse2d::getSlice(mat<double> &x, mat <double> &out, const int z){
  memcpy(&out.d[0], &x(z,0,0), dims[0]*dims[1]*sizeof(double));
}
void sparse2d::putSlice(mat<double> &x, mat <double> &out, const int z){
    memcpy( &x(z,0,0), &out.d[0], dims[0]*dims[1]*sizeof(double));
}

/* ------------------------------------------- */
double sparse2d::my_getStep(double &step, mat<double> &obs, double gy, mat<double> &dchi, mat<double> &y,
			    mat<double> &noise, mdepthall_t &m, int nint){

  /* --- beta is a parameter that will be used to get the step size --- */
  double beta = 0.75;
  double bestchi = 1.e13, beststep = step;

  /* --- Seems to improve convergence --- */
  step *= 1./beta;
  
  /* --- Get chi2 at the current step size --- */
  vector<double> gx, istep;
  gx.resize(1);
  gx[0] = evalStep(step, obs, dchi, y, noise, m);
  istep.resize(1);
  istep[0] = step;
  bestchi = gx[0];
  beststep = step;
  
  /* --- init loop to braket the step size --- */
  int it = 0, idx=0;
  double gxold = 1.e12;

  /* --- braket step size assuming that gy > step > gx[0] ---*/
  mat<double> x = takeStep(y, dchi, step, true);


  /* --- Compute term to check progress --- */
  double norm2 = 0, dprod = 0;
  for(int ii=0;ii<(int)y.d.size();ii++){
    norm2 += sqr(x.d[ii] - y.d[ii]);
    dprod += dchi.d[ii] * (x.d[ii]-y.d[ii]);
  }
  double term = gy + dprod + 0.5/step * norm2;
  fprintf(stderr,"BCHI[%e] = %f, term=%e\n", istep[idx], gx[idx], term);

  /* --- iterate --- */
  while(gx[idx] < gxold){
  //while(gx[idx] > term && (it++ < 10)){
    gxold = gx[idx];
    istep.push_back(istep[idx++]*beta);
    gx.push_back(evalStep(istep[idx], obs, dchi, y, noise, m));
    
    if(gx[idx] < bestchi){
      bestchi = gx[idx];
      beststep = istep[idx];
    }
    
    mat<double> x = takeStep(y, dchi, istep[idx], true);
    
    norm2 = 0, dprod = 0;
    for(int ii=0;ii<(int)y.d.size();ii++){
      norm2 += sqr(x.d[ii] - y.d[ii]);
      dprod += dchi.d[ii] * (x.d[ii]-y.d[ii]);
    }
    term = gy + dprod + 0.5/step * norm2;
    
    fprintf(stderr,"BCHI[%e] = %f, term=%e\n", istep[idx], gx[idx], term);
    it++;
  }
  
  
  bestchi = gx[max(idx-1,0)];
  step = istep[max(idx-1,0)];


  return bestchi;
}

/* ------------------------------------------- */

double sparse2d::parabStepInt(double sa, double sb, double sc, double ga, double gb, double gc,
			      mat<double> &obs, mat<double> &dchi, mat<double> &y, mat<double> &noise,
			      mdepthall_t &m, int nint, double &bestchi){

  /* --- Order according to step size --- */
  mysort3(sa, sb, sc, ga, gb, gc);
  
  double sd = 0.0, gd = 0.0;
  bestchi = min(sa,min(sb,sc));

  /* --- Parabola fits can overshoot, get a symmetric step --- */
  if(ga > gc) sd = sb - (sc - sb);
  else        sd = sb + (sb - sa);
  gd = evalStep(sd, obs, dchi, y, noise, m);
  fprintf(stderr,"ECHI[%e] = %f\n", sd, gd);

  
  mysort4(ga, gb, gc, gd, sa, sb, sc, sd);

  
  for(int ii = 0; ii<nint;ii++){

    /* --- Order according to step size --- */
    mysort3(sa, sb, sc, ga, gb, gc);

    /* --- Get parab coeff. --- */
    double me = 0.0, st = 0.0;
    MeSt(sa, sb, sc, me, st);
    vector<double> cc = parab_fit<double>((sa-me)/st, (sb-me)/st, (sc-me)/st, ga, gb, gc);
    sd = (-0.5 * cc[1]/cc[2])*st + me;

    /* --- Eval Chi at the minimum of the parabola ---*/
    gd = evalStep(sd, obs, dchi, y, noise, m);
    fprintf(stderr,"BCHI[%e] = %f\n", sd, gd);

    /* --- Is it smaller than the others? -> improve interval for next iter ---*/
    bool better = false;
    if((gd<ga) || (gd<gb) || (gd<gc)) better = true;
    mysort4(ga,gb,gc,gd,sa,sb,sc,sd);

    if(!better) break;
    
  } // ii


  /* --- return best step, which has been ordered to the "a" location ---*/
  fprintf(stderr,"FCHI[%e] = %f\n", sa, ga);
      
  bestchi = ga;
  return sa;
}
  

/* ------------------------------------------- */


double sparse2d::getStep(double &alpha, mat<double> &obs, mat<double> &dchi, mat<double> &pp,
			 mat<double> &noise, mdepthall_t &m, double tol){

  /* --- Init parameters to bracket the optimized step and call mnbrak--- */
  double ax = 0.0, bx = alpha, cx,fa,fb,fc, xm=0.0;
  
  mnbrak(ax,bx,cx,fa,fb,fc,obs,dchi,pp,noise, m);

  /* --- Refine the search with the brent method --- */
  double cs = brent(ax,bx,cx,tol,xm,obs,dchi,pp,noise,m);
  
  alpha = xm;

  return cs;
}

/* ------------------------------------------- */

double sparse2d::evalStep(double alpha,  mat<double> &obs, mat<double> &dchi,
			  mat<double> &pp, mat<double> &noise, mdepthall_t &m){

  mat<double> ipp = takeStep(pp, dchi, alpha, true);
  double chi = meritFunction(obs, ipp, dchi, noise, m, 0);
  return chi;
}

/* ------------------------------------------- */
void sparse2d::checkParameters(mat<double> &xnew){
  
  /* --- Convert to parameter space and check limits --- */
  wlt.transformSlices(xnew, dwt_inverse);


  for(int pp = 0; pp<npar; pp++)  
    for(int yy = 0; yy < dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++)
	xnew(pp,yy,xx) = atm[0]->checkParameter(xnew(pp,yy,xx)*scalingParameters[pp], pp) / scalingParameters[pp];
	  
  /* --- Convert back to projected space--- */
  wlt.transformSlices(xnew, dwt_forward);
  
}

/* ------------------------------------------- */


mat<double> sparse2d::takeStep(mat<double> &par, mat<double> &grad, double step,
			       bool do_thres){

  /* --- Define variables --- */
  mat<double> xnew = par;

  /* --- take step --- */
  for(int pp = 0; pp < npar; pp++)
    for(int yy = 0; yy < dims[0]; yy++)
      for(int xx = 0; xx < dims[1]; xx++){
	xnew(pp,yy,xx) = ( par(pp,yy,xx) - grad(pp,yy,xx) * step *
			    LCommon[pp] );
	
      }


  //checkParameters(xnew);

  
  
  /* --- Do thresholding? --- */
  if(do_thres)
    threshold(xnew, thres);

  checkParameters(xnew);

  
  return xnew;
}
/* ------------------------------------------- */


double FMAX(double a, double b)
{
  return (a>b)?a:b;
}

/* ------------------------------------------- */


#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void sparse2d::mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb,
	      double &fc, mat<double> &obs, mat<double> &dchi, mat<double> &pp,
		      mat<double> &noise, mdepthall_t &m){

    double ulim,u,r,q,fu,dum;

  fa=evalStep(ax, obs, dchi, pp, noise, m);
  fb=evalStep(bx, obs, dchi, pp, noise, m);
  if (fb > fa) {
    SHFT(dum,ax,bx,dum)
    SHFT(dum,fb,fa,dum)
  }
  cx=(bx)+GOLD*(bx-ax);
  fc=evalStep(fc, obs, dchi, pp, noise, m);
  while (fb > fc) {
    r=(bx-ax)*(fb-fc);
    q=(bx-cx)*(fb-fa);
    u=(bx)-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(bx)+GLIMIT*(cx-bx);
    // fprintf(stderr,"[%14.7E,%14.7E,%14.7E]\n",u,ulim,(bx-u)*(u-cx));
    if ((bx-u)*(u-cx) > 0.0) {
      fu=evalStep(u, obs, dchi, pp, noise, m);;
      if (fu < fc) {
        ax=(bx);
        bx=u;
        fa=(fb);
        fb=fu;
        return;
      } else if (fu > fb) {
        cx=u;
        fc=fu;
        return;
      }
      u=(cx)+GOLD*(cx-bx);
      fu=evalStep(u, obs, dchi, pp, noise, m);
    } else if ((cx-u)*(u-ulim) > 0.0) {
      fu=evalStep(u, obs, dchi, pp, noise, m);
      if (fu < fc) {
        SHFT(bx,cx,u,cx+GOLD*(cx-bx))
        SHFT(fb,fc,fu,evalStep(u, obs, dchi, pp, noise, m))
      }
    } else if ((u-ulim)*(ulim-cx) >= 0.0) {
      u=ulim;
      fu=evalStep(u, obs, dchi, pp, noise, m);
    } else {
      u=(cx)+GOLD*(cx-bx);
      fu=evalStep(u, obs, dchi, pp, noise, m);
    }
    SHFT(ax,bx,cx,u)
    SHFT(fa,fb,fc,fu)
  }
  
}

#undef GOLD
#undef GLIMIT
#undef TINY

/* ------------------------------------------- */

#define ITMAX 25
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

double sparse2d::brent(double ax, double bx, double cx,double tol,double &xmin,
	       mat<double> &obs, mat<double> &dchi, mat<double> &pp, mat<double> &noise, mdepthall_t &m)
{
  int iter;
  double a,b,dd,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=evalStep(x, obs, dchi, pp, noise, m);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=dd;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        dd=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
        dd=p/q;
        u=x+dd;
        if (u-a < tol2 || b-u < tol2) dd=SIGN(tol1,xm-x);
      }
    } else {
      dd=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(dd) >= tol1 ? x+dd : x+SIGN(tol1,dd));
    fu=evalStep(u, obs, dchi, pp, noise, m);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }
  //exit(fprintf(stderr,"Too many iterations in brent"));
  xmin=x;
  return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

/* ------------------------------------------- */


double sparse2d::getNorm(int n, double *x, int incx){
  /*
    Ported to C++ directly from DNRM2 from LAPACK
   */
  
  /* --- Init vars --- */
  double absxi=0, norm=0, scale=0, ssq=0;
  int ix = 0;
  
  /* --- Check what to do based on the number of elements ---*/
  if(n < 1 || incx < 1) norm = 0;
  else if(n == 1) norm = abs(x[0]);
  else{
    scale = 0.0;
    ssq = 1.0;
    for(ix = 0; ix< (n*incx); ix+=incx){
      if(x[ix] != 0){
	absxi = abs(x[ix]);
	if(scale < absxi){
	  ssq = 1.0 + ssq*sqr(scale/absxi);
	  scale = absxi;
	}else{
	  ssq += sqr(absxi/scale);
	}//else
      }//if
    }//for ix
    norm = scale * sqrt(ssq);
  }//else

  return norm;
}



