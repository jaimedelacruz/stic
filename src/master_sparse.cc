#include <iostream>
#include <string>
#include <netcdf>
#include <omp.h>
#include "io.h"
#include "cmemt.h"
#include "atmosphere.h"
#include "sparse.h"
#include <mpi.h>
#include "comm.h"
#include "depthmodel.h"
#include "clte.h"
#include "crh.h"
//
using namespace netCDF;
using namespace std;
//

void slaveInversion(iput_t &iput, mdepthall_t &m, mat<double> &obs, mat<double> &x, mat<double> &chi2){

  /* --- Init dimensions --- */
  unsigned long ntot = (unsigned long)(x.size(0) * x.size(1));
  int ncom = (int)(std::floor(ntot / (double)iput.npack));
  if((unsigned long)(ncom * iput.npack) != ntot) ncom++;
  int nprocs = iput.nprocs;
  int iproc = 0;
  unsigned long ipix = 0;
  int tocom = ncom;

  mat<double>  dsyn; // dummy 
  int compute_gradient = 0; // dummy parameter here
  chi2.set({x.size(0), x.size(1)});
  

  if(nprocs > 1){

  
    /* --- Init slaves with first package --- */
    for(int ss = 1; ss<=min(nprocs-1,ncom); ss++)
      comm_master_pack_data(iput, obs, x, ipix, ss, m, compute_gradient);

    int per  = 0;
    int oper  = -1;
    float pno =  100.0 / (float(ntot) - 1.0);
    unsigned long irec = 0;
    cerr << "\rProcessed -> "<<per<<"% ";


    /* --- manage packages as long as needed --- */
    while(irec < ntot){
    
      // Receive processed data from any slave (iproc)
      comm_master_unpack_data(iproc, iput, obs, x, chi2, irec, dsyn, compute_gradient, m);
      per = irec * pno;
      //cerr << ipix << " " << irec<<" "<<ntot << endl;
      // Send more data to that same slave (iproc)
      if(ipix < ntot) comm_master_pack_data(iput, obs, x, ipix, iproc, m, compute_gradient);
    
      // Keep count of communications left
      tocom--;
    
      // Printout
      if(per > oper){
	oper = per;
	cerr << "\rProcessed -> "<<per<<"% ";
      }
    }
  
    cerr << " "<<endl;

  }
  
}

void master_inverter(mdepthall_t &model, mat<double> &pars, mat<double> &obs, mat<double> &w, iput_t &input)
{

  int ndep = (int)model.ndep, nx = input.nx, ny = input.ny;
  mdepth_t m(ndep);
  atmos *atm;
    
  /* --- Init atmos --- */

  if(input.atmos_type == string("lte")) atm = new clte(input, 4.44);
  else if(input.atmos_type == string("rh")) atm = new crh(input, 4.44);
  else{
    cerr<<"master_invert: ERROR, atmos unknown -> "<<input.atmos_type<<endl;
    exit(0);
  }


  /* --- (TO-DO, change this!) --- */
  
  vector<instrument*> inst;
  int nreg = atm->input.regions.size();
  inst.resize(nreg);
  
  for(int kk = 0; kk<nreg; kk++){
    if(atm->input.regions[kk].inst == "spectral")
      inst[kk] = new spectral(atm->input.regions[kk], 1);
    else inst[kk] = new instrument();
  }
  atm->inst = &inst[0];


  /* --- Loop and invert --- */
  
  int per = 0, oper = -1, kk = 0;
  
  for(int yy = 0; yy<input.ny; yy++)
    for(int xx = 0; xx<input.nx; xx++){

      /* --- Copy data to single pixel model --- */
      
      memcpy(&m.cub.d[0], &model.cub(yy,xx,0,0), 11*ndep*sizeof(double));


      /* --- invert --- */
      
      atm->fitModel( m, input.npar, &pars(yy,xx,0),
		    (int)(input.nw_tot*input.ns), &obs(yy,xx,0,0), w);

      
      per = ++kk / (nx*ny);
      if(per > oper){
	oper = per;
	fprintf(stderr,"\rmaster_invert: %d%s",per,"%");
      }
    }
  
  fprintf(stderr,"master_invert: %d%s",per,"%");


  /* --- Clean-up --- */
  delete atm;
  for(auto &it: inst)
    delete it;
}


void do_master_sparse(int myrank, int nprocs,  char hostname[]){

  /* --- Printout number of processes --- */
  cerr << "SparseInv: Initialized with "<<nprocs <<" processes"<<endl;
  mat<double> model, obs, wav, w, syn, chi2;;
  mdepthall_t im;

  
  /* --- Read input files and fill in variables---*/
  iput_t input = read_input("input.cfg", true);
  read_lines("lines.cfg", input, input.verbose);
  wav.set(vector<int>{input.nw_tot});
  wav.d = fill_lambdas(input, false);

  
  input.nprocs = nprocs;
  input.myid = (string)"master, "; 

  /* --- Set OpenMP threads in the master --- */
  
  omp_set_num_threads(input.master_threads);
  
  
  /* --- Open input files with the models and profiles --- */
  
  io ipfile, omfile, imfile(file_exists(input.imodel), NcFile::read);
  
  bool inversion = false;
  if(input.mode == 1 || input.mode == 3) inversion = true;

  if(inversion)
    ipfile.initRead(file_exists(input.iprof), NcFile::read);
  
  
  
  /* --- Get dimensions and fill input struct --- */
  
  input.nt = 1;
  if(inversion){
    int off = 0;
    vector<int> odims = ipfile.dimSize(string("profiles"));
    if(odims.size() == 5) {
      off = 1;
      input.nt = odims[0];
    }
    input.ny = odims[0+off];
    input.nx = odims[1+off];
    input.ns = odims[3+off];
    if(input.nw_tot != odims[2+off]){
      cerr << input.myid <<"ERROR, number of wavelength points in regions != "<<input.iprof<<endl;
      exit(0);
    }
  }else{
    int off = 0;
    vector<int> odims = imfile.dimSize(string("temp"));
    if(odims.size() == 4) {
      off = 1;
      input.nt = odims[0];
    }
    input.ny = odims[0+off];
    input.nx = odims[1+off];
    input.ns = 4;

    obs.set(vector<int>{input.ny, input.nx, input.nw_tot, input.ns});
  }
  vector<int> dims = {input.ny, input.nx, input.nw_tot, input.ns};
  
 
  /* --- Read time-independent variables --- */
  
  if(inversion){
    ipfile.read_Tstep<double>(string("weights"), w);
    ipfile.read_Tstep<double>(string("wav"), wav);

    if(w.d.size() == 0) {
      w.set({input.nw_tot, input.ns});
      fill(w.d.begin(), w.d.end(), 3.0e-3); // Default noise?
    }
  }

  // Read model for tstep = 0
  input.boundary = im.read_model2(input.imodel, 0, true);
  input.ndep = im.ndep;
  
  /* ---
     Open output files and init vars to store results 
     (dimension = 0 means unlimited)
     We decide here if the variables are going to be
     stored as floats or doubles, regardless of their 
     type in memory
     --- */
  io opfile(input.oprof,  NcFile::replace);
  opfile.initDim({"time","y", "x", "wav", "stokes"},{0,input.ny, input.nx, input.nw_tot, input.ns});
  opfile.initVar<double>(string("profiles"), {"time","y", "x", "wav", "stokes"});
  opfile.initVar<double>(string("wav"), {"wav"});
  opfile.write_Tstep<double>(string("wav"), wav);


  
  if(inversion){
    opfile.initVar<float>(string("weights"), {"wav", "stokes"});
    opfile.write_Tstep<double>(string("weights"), w);

    omfile.initRead(input.omodel, NcFile::replace);

    vector<int> dims = ipfile.dimSize("profiles");
    if(dims.size() == 5) dims.erase(dims.begin(), dims.begin()+1); // remove time
 
    //
    // Propagate information to slaves via MPI
    //
    double mmax = im.cub(0,0,9,0), mmin = im.cub(0,0,9,0);
    for(int zz = 1; zz < input.ndep; zz++){
      if(im.cub(0,0,9,zz) > mmax) mmax = im.cub(0,0,9,zz);
      if(im.cub(0,0,9,zz) < mmin) mmin = im.cub(0,0,9,zz);
    }
    
    /* --- Place Nodes --- */
    //input.npar = set_nodes(input.nodes,  mmin, mmax, input.verbose);
    
    {
      vector<double> itau;
      itau.resize(input.ndep);
      for(int ii=0;ii<input.ndep;ii++) itau[ii] = im.cub(0,0,9,ii);
      input.npar = set_nodes(input.nodes, itau, input.verbose);
    }
    
  } // if inversion



  
  /* --- MPI bussiness --- */
  
  if(input.npack < 0){
    input.npack = (int)((double)(input.ny*input.nx) / nprocs) + 1;
    if(input.verbose) cerr<<input.myid<<"Using NPACK="<<input.npack<<endl;
  }
  
  comm_get_buffer_size(input);
  MPI_Barrier(MPI_COMM_WORLD); // Wait until all processors reach this point
  comm_send_parameters(input);
  
  if(inversion) comm_send_weights(input, w); //  

  
  
  
  /* --- Init sparse class --- */
  
  sparse2d inv;
  if(input.mode == 3)
    inv.init(input,  dims, input.sparse_threshold, 
		 string2wavelet(input.wavelet_type), input.wavelet_order, 
		 spt_hard, input.master_threads);
  
  //
  // Main loop
  //
  for(int tt = 0; tt<input.nt; tt++){ // Loop in time

    /* --- Read Tstep data --- */
    
    mat<double> pweight;
    if(inversion){
      ipfile.read_Tstep(string("profiles"), obs, tt);
      if(ipfile.is_var_defined((string)"pixel_weights")){
	ipfile.read_Tstep<double>(string("pixel_weights"), pweight, tt);
      }
    }
    
    if(tt > 0) im.read_model2(input.imodel, tt, true);

    /* --- Check dimensions in inversion mode --- */
    
    if(inversion){
      if(obs.size(0) != im.cub.size(0) && obs.size(1) != im.cub.size(1)){
	cerr << input.myid <<"ERROR, the input model and the observations do not have the same dimensions in X,Y axes:"<<endl;
	cerr << "   -> "<<input.imodel<<" "<<formatVect<int>(im.temp.getdims())<<endl;
	cerr << "   -> "<<input.iprof<<" "<<formatVect<int>(obs.getdims())<<endl;
	
	comm_kill_slaves(input, nprocs);
	exit(0);
      }

      if(obs.size(2) != input.nw_tot){
	cerr << input.myid <<"ERROR, number of wavelenghts in input.cfg ["<<input.nw_tot<< "] does not match the observations ["<< obs.size(2) << "]"<<endl;
	comm_kill_slaves(input, nprocs);
	exit(0);
	
      }

      /* --- get free-parameters from input model --- */
      
      
    } // inversion mode

    im.model_parameters2(model, input.nodes);


    
    /* --- Invert data --- */
    
    if     (input.mode == 1){

      if(nprocs == 1)
	master_inverter(im, model, obs, w, input);
      else
	slaveInversion(input, im, obs, model, chi2); // implemented above!
      
    }else if(input.mode == 2) slaveInversion(input, im, obs, model, chi2); // it won't invert if mode == 2
    else if(input.mode == 3) inv.SparseOptimization(obs, model, w, im, pweight);
    
    
    if(inversion){

      /* --- Expand fitted parameters into depth-stratified atmos --- */

      im.expandAtmos(input.nodes, model, input.dint);
      
      /* --- init output for inverted model --- */
      
      if(tt == 0){
	vector<int> odim = model.getdims();
	odim.insert(odim.begin(), 0);
	omfile.initDim({"time","y", "x", "par"},odim);
	omfile.initVar<float>(string("model"), {"time","y", "x", "par"});
      }
    

      /* --- Write model parameters, profiles and depth-stratified atmos --- */
      
      omfile.write_Tstep(string("model"), model, tt);
      im.write_model2(input.oatmos, tt);
    }

    opfile.write_Tstep(string("profiles"), obs, tt);

  }

  
  
  /* --- Tell slaves to exit while(1) loop --- */
  
  comm_kill_slaves(input, nprocs);
  
}

