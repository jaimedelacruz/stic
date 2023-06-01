/*
  Wavelet class
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
  Dependencies: cmemt.h, GSL (library)
 */
#include <iostream>
#include <cstring>
//#include <omp.h>
#include "wavelet.h"
#include "cmemt.h"
//
using namespace std;
//

wavelet_type string2wavelet(const std::string fam){
  
  if(!fam.compare(std::string("haar"))){
    return dwt_haar;
  } else if(!fam.compare(std::string("haar_c"))){
    return dwt_haar_c;
  } else if(!fam.compare(std::string("daub"))){
    return dwt_daub;
  } else if(!fam.compare(std::string("daub_c"))){
    return dwt_daub_c;
  } else if(!fam.compare(std::string("bspline"))){
    return dwt_bspline;
  } else if(!fam.compare(std::string("bspline_c"))){
    return dwt_bspline_c;
  } else{
    std::cerr <<"string2wavelet: Error, spline type not implemented, falling back to Daubechies"<<std::endl;
    return dwt_daub;
  }
}
//

void wavelet::init(const vector<int> dims1, unsigned int nt, const wavelet_type wtype1, unsigned int ord){
  
  string inam = "wavelet::init: ";

  //
  // Wavelet type
  //
  if(wtype1 == dwt_haar){
    fam = "Haar";
    wtype = gsl_wavelet_haar;
    otype = wtype1;
  } else if(wtype1 == dwt_haar_c){
    fam = "Haar_c";
    wtype = gsl_wavelet_haar_centered;
    otype = wtype1;
  } else if(wtype1 == dwt_daub){
    fam = "Daubechies";
    wtype = gsl_wavelet_daubechies;
    otype = wtype1;
  } else if(wtype1 == dwt_daub_c){
    fam = "Daubechies_c";
    wtype = gsl_wavelet_daubechies_centered;
    otype = wtype1;
  } else if(wtype1 == dwt_bspline){
    fam = "Bspline";
    wtype = gsl_wavelet_bspline;
    otype = wtype1;
  }else if(wtype1 == dwt_bspline_c){
    fam = "Bspline_c";
    wtype = gsl_wavelet_bspline_centered;
    otype = wtype1;
  } else{
    cerr << inam <<"Error, spline type not implemented, falling back to Daubechies [4]"<<endl;
    wtype = gsl_wavelet_daubechies;
    fam = "Daubechies";
    otype = dwt_daub;
    ord = 4;
  }
  order = ord;


  //
  // Check number of dimensions
  //
  ndim = dims1.size();
  if(ndim > 3){
    cerr << inam << "Error, number of dimensions must be 3 >= ndim > 0 -> ndim = " << ndim << endl;
    exit(0);
  }
  dims = dims1;

  //
  // allocate arrays
  //
  nthreads = nt;
  
  // Work space and wavelets, one per thread (array of pointers)
  w.resize(nthreads);
  work.resize(nthreads);
  if(ndim > 1) work2.resize(nthreads);
  if(ndim > 2) work3.resize(nthreads);
  //
  for(unsigned int nn = 0; nn < nthreads; nn++){
    w[nn] = gsl_wavelet_alloc(wtype, order);
    work[nn] = gsl_wavelet_workspace_alloc (dims[0]);
    if(ndim > 1) work2[nn] = gsl_wavelet_workspace_alloc (dims[1]);
    if(ndim > 2) work3[nn] = gsl_wavelet_workspace_alloc (dims[2]);
  }
  //
  cerr << inam << "Using "<<fam<<"["<<order<<"], nthreads="<<nthreads<<", dims"<<formatVect<int>(dims)<<endl;
}
//
void wavelet::transform(mat<double> &dat, wavelet_dir idir){
  string inam;
  bool forw;
  gsl_wavelet_direction dir;
  
  if(idir == dwt_forward) {
    inam = "wavelet::transform [forward]: ";
    forw = true;
    dir = gsl_wavelet_forward;
  }else {
    inam = "wavelet::transform [inverse]: ";
    forw = false;
    dir  = gsl_wavelet_backward;
  }
 
 
  //
  // Perform wavelet transform, the order of loops depends on the direction of the transforms
  // in multi-dimensional DWT
  //
  if(ndim == 1) {
    //   cerr << inam <<"computing 1d transform ... ";
    gsl_wavelet_transform(w[0], &dat.d[0], 1, dat.size(0), dir, work[0]);
  }else if(ndim == 2){

    // cerr << inam <<"computing 2d transform ... ";
    
    int xx, yy, id;
    unsigned long dx = 1;
    unsigned long dy = dims[1];
    if(forw){ // forward
      
      //
      // 2D FORWARD
      //      
#pragma omp parallel default(shared) private(xx,yy,id) num_threads(nthreads)
      {
	//	id = omp_get_thread_num(); // get thread number
	
	// Fist do all rows
#pragma omp for 
	for(yy = 0; yy < dims[0]; yy++) 
	  gsl_wavelet_transform(w[id], &dat(yy,0), dx, dat.size(1), dir, work[id]);
	
	// Now columns
#pragma omp for
	for(xx = 0; xx < dat.size(1); xx++)
	  gsl_wavelet_transform(w[id], &dat(0,xx), dy, dat.size(0), dir, work2[id]);	
	
	
      } // parallel block
    } else{ // Inverse 

      //
      // 2D INVERSE
      //
#pragma omp parallel default(shared) private(xx,yy,id) num_threads(nthreads)
      {
	//id = omp_get_thread_num(); // get thread number
	// Now columns
#pragma omp for // Schedule 1 horizontal slice at the time
	for(xx = 0; xx < dat.size(1); xx++)
	  gsl_wavelet_transform(w[id], &dat(0,xx), dy, dat.size(0), dir, work2[id]);
	
	// Now rows
#pragma omp for  // Schedule 1 horizontal slice at the time
	for(yy = 0; yy < dims[0]; yy++) 
	  gsl_wavelet_transform(w[id], &dat(yy,0), dx, dat.size(1), dir, work[id]);
      } // parallel block
    } // inverse
    
    } // ndim == 2
    else if(ndim == 3){
    
      //  cerr << inam <<"computing 3d transform ... ";
    int xx, yy, zz, id;
    unsigned long dx = 1;
    unsigned long dy = dims[2];
    unsigned long dz = dims[1]*dims[2];

    if(forw){

      //
      // 3D FORWARD
      //
#pragma omp parallel default(shared) private(xx,zz,yy,id) num_threads(nthreads)
      {
	//id = omp_get_thread_num(); // get thread number
	// Transform along "X"
#pragma omp for  // Schedule 1 horizontal slice at the time
	for(zz = 0; zz < dims[0]; zz++)
	  for(yy = 0; yy < dims[1]; yy++){
	    //
	    gsl_wavelet_transform(w[id], &dat(zz,yy,0), dx, dat.size(2), dir, work3[id]);
	    //
	  }
	
	// Now "Y"
#pragma omp for // Schedule 1 horizontal slice at the time
	for(zz = 0; zz < dat.size(0); zz++)
	  for(xx = 0; xx < dat.size(2); xx++){
	    //
	    gsl_wavelet_transform(w[id], &dat(zz,0,xx), dy, dat.size(1), dir,work2[id]);
	    //
	  }

	// Now "Z"
#pragma omp for // Schedule 1 horizontal slice at the time		
	for(yy = 0; yy<dat.size(1);yy++) 
	  for(xx = 0; xx < dat.size(2); xx++){
	    //
	    gsl_wavelet_transform(w[id], &dat(0,yy,xx), dz, dat.size(0), dir,work[id]);
	    //
	  }

      } // parallel block
    } else {

      //
      // 3D INVERSE
      // 
#pragma omp parallel default(shared) private(xx,zz,yy,id) num_threads(nthreads)
      {
	//	id = omp_get_thread_num(); // get thread number

	// Now "Z"
#pragma omp for
	for(yy = 0; yy<dat.size(1);yy++) 
	  for(xx = 0; xx < dat.size(2); xx++){
	    //
	    gsl_wavelet_transform(w[id], &dat(0,yy,xx), dz, dat.size(0), dir, work[id]);
	    //
	  }	
	
	// Now "Y"
#pragma omp for
	for(zz = 0; zz < dat.size(0); zz++)
	  for(xx = 0; xx < dat.size(2); xx++){
	    //
	    gsl_wavelet_transform(w[id], &dat(zz,0,xx), dy, dat.size(1), dir, work2[id]);
	    //
	  }

	
	// Transform along "X"
#pragma omp for
	for(zz = 0; zz < dims[0]; zz++)
	  for(yy = 0; yy < dims[1]; yy++){
	    //
	    gsl_wavelet_transform(w[id], &dat(zz,yy,0), dx, dat.size(2), dir, work3[id]);
	    //
	  }
      } // parallel block
    } // Inverse
  } // ndim == 3
}

void wavelet::transformSlices(mat<double> &dat, wavelet_dir idir, int id1){
  string inam = "wavelet::transformSlices: ";
  
  vector<int> nn = dat.getdims();
  if(nn.size() != 3){
    cerr << inam << "ERROR, data must be 3D "<<formatVect<int>(nn)<<endl;
    return;
  }

  bool forw;
  gsl_wavelet_direction dir;
  
  if(idir == dwt_forward) {
    forw = true;
    dir = gsl_wavelet_forward;
  }else {
    forw = false;
    dir  = gsl_wavelet_backward;
  }

  unsigned long dx = 1;
  unsigned long dy = nn[2];
  int xx, yy, zz, id;

  
  if(forw){ // forward
    
    //
    // 2D FORWARD
    //      
#pragma omp parallel default(shared) private(zz,xx,yy,id) num_threads(nthreads)
    {
      //id = omp_get_thread_num(); // get thread number
      
      // Fist do all rows
#pragma omp for
      for(zz = 0; zz < nn[0]; zz++)
	for(yy = 0; yy < dims[0]; yy++) 
	  gsl_wavelet_transform(w[id], &dat(zz,yy,0), dx, nn[2], dir, work[id]);
      
      // Now columns
#pragma omp for 
      for(zz = 0; zz < nn[0]; zz++)
	for(xx = 0; xx < dat.size(1); xx++)
	  gsl_wavelet_transform(w[id], &dat(zz,0,xx), dy, nn[1], dir, work2[id]);	
      
      
    } // parallel block
  } else{ // Inverse 
    
    //
    // 2D INVERSE
    //
#pragma omp parallel default(shared) private(xx,zz,yy,id) num_threads(nthreads)
    {
      // id = omp_get_thread_num(); // get thread number
      // Now columns
#pragma omp for
      for(zz = 0; zz < nn[0]; zz++)
	for(xx = 0; xx < dat.size(1); xx++)
	gsl_wavelet_transform(w[id], &dat(zz,0,xx), dy, nn[1], dir, work2[id]);
      
      // Now rows
#pragma omp for
            for(zz = 0; zz < nn[0]; zz++)
	      for(yy = 0; yy < dims[0]; yy++) 
		gsl_wavelet_transform(w[id], &dat(zz,yy,0), dx, nn[2], dir, work[id]);
    } // parallel block
  } // inverse
  
 
}

  
