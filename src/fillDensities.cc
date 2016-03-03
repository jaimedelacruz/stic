#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include "io.h"
#include "ceos.h"
#include "depthmodel.h"
#include "netcdf.h"
#include "cmemt.h"

using namespace std;
using namespace netCDF;

int main(int narg, char *argv[])
{
  
  /* --- Define I/O --- */
  
  string ifile, ofile, scal;
  if(narg != 5){
    cerr<<"usage: ./fillDensities filein.nc fileout.nc [z | tau] [0 | 1]"<<endl;
    exit(0);
  }
  
  ifile = string(argv[1]);
  ofile = string(argv[2]);
  scal = string(argv[3]);
  int hydrostat = atoi(argv[4]);
  cerr<<"Compute hydrostatic eq. -> "<<hydrostat <<endl;
  //
  if(scal != string("z") && scal != string("tau")){
    cerr<<"ERROR, Unknown depth scale ["<<scal<<"] -> Must be 'z' or 'tau'"<<endl;
    exit(0);
  }
  bool do_z = false;
  if(scal == string("z")) do_z = true;
  
    
  /* --- Get dimensions --- */

  vector<int> dim;
  {
    io ifil(ifile, NcFile::read);
    dim = ifil.dimSize("temp");
  }
  int ndep = dim[3];

  /* --- Read 1st time-step and check densities/pressures --- */
  
  mdepthall_t m;
  int bound = m.read_model2(ifile,(int)0, false);
  //
  int touse = -1;
  if     (fabs(m.cub(0,0,6,1) - m.cub(0,0,6,0)) > 0.0) touse = 0;
  else if(fabs(m.cub(0,0,7,1) - m.cub(0,0,7,0)) > 0.0) touse = 1;
  else if(fabs(m.cub(0,0,8,1) - m.cub(0,0,8,0)) > 0.0) touse = 2;
  //
  if(touse == -1){
    cerr<<"ERROR, you must provide at least one non-zero density/pressure scale"<<endl;
    exit(0);
  }

  vector<string> names = {"Pgas","Dens", "Nelect"};
  cerr<<"EOS: using ["<<names[touse]<<"]"<<endl;
  


  
  /* --- Some vars --- */

  vector<double> kappa;
  kappa.resize(ndep);
  int nw = 1;
  double wav = 5000.0, scat = 0.0;
  vector<float> frac, part;
  float na=0, ne=0;
  mdepth mm(m.ndep);
  
  /* --- Loop t-steps --- */
  
  for(int tt = 0; tt<dim[0]; tt++){

    int npix = dim[1] * dim[2];
    int per = 0, oper = -1, ipix = 0;

      
    /* --- Read Time step --- */
    
    if(tt > 0) m.read_model2(ifile, tt, false);

    
    
    /* --- loop pixels --- */
    
    for(int yy = 0; yy < dim[1]; yy++){
      for(int xx = 0; xx < dim[2]; xx++){

	memcpy(&mm.cub(0,0), &m.cub(yy,xx,0,0), m.ndep*11*sizeof(double));
	
	double *temp = mm.temp;
	double *pgas = mm.pgas;
	double *rho = mm.rho;
	double *nne = mm.nne;
	double *tau = mm.ltau;
	double *z   = mm.z;
	double *pel = mm.pel;
	double *cmass=mm.cmass;
	

	
	/* --- EOS: Fill pressure / desities --- */
	if(hydrostat == 0 ){
	  for(int kk = 0; kk < ndep; kk++){
	    if(touse == 0)
	      nne[kk] = m.eos.nne_from_T_Pg( temp[kk], pgas[kk], rho[kk]);
	    else if(touse == 1)
	      nne[kk] = m.eos.nne_from_T_rho(temp[kk], pgas[kk], rho[kk]);
	    else if(touse == 2)
	      rho[kk] =  m.eos.nne_from_T_rho(temp[kk], pgas[kk], nne[kk]);
	    
	    m.eos.store_partial_pressures(ndep, kk, m.eos.xna, m.eos.xne);
	    m.eos.read_partial_pressures(kk, frac, part, na, ne);
	    
	    
	    /* --- Get kappa_5000 --- */
	    
	    m.eos.contOpacity(temp[kk], nw,  &wav, &kappa[kk],
			      &scat, frac, na, ne);
	    
	  }
	  /* --- compute the depth-scale --- */
	  
	  if(do_z){
	    //
	    tau[0] = 0.5 * kappa[0] * (z[0] - z[1]);
	    for(int k = 1; k < ndep; k++)
	      tau[k] = tau[k-1] + 0.5 * (kappa[k-1] + kappa[k]) * (z[k-1] - z[k]);
	    
	    for(int k = 0; k < ndep; k++){
	      tau[k] = log10(tau[k]);
	    }
	  }else{
	    
	    z[0] = 0.0;
	    double otau = pow(10.0, tau[0]);
	    
	    for(int k = 1; k < ndep; k++){
	      double itau = pow(10.0, tau[k]);
	      z[k] = z[k-1] - 2.0 * (itau - otau) / (kappa[k] + kappa[k-1]);
	      otau = itau;
	    } // k
	  } // tau
	}else{
	  
	  /* --- Do hydrostatic eq. --- */
	  for(int kk=0; kk<m.ndep; kk++) mm.tau[kk] = pow(10.0, tau[kk]);
	  mm.getPressureScale(bound, m.eos);
	  //cerr<<"BOUND="<<bound<<endl;


	}

	memcpy(&m.cub(yy,xx,0,0), &mm.cub(0,0), 11*m.ndep*sizeof(double));
	
	per = ++ipix * (100. / npix);
	fprintf(stderr,"\r[t=%d, pix=%d/%d] Processing -> %d%s", tt,ipix-1,npix-1,per,"%");
	
      } // xx
    }//yy
    
    cerr<<endl;
    
    
    /* --- Write to file --- */

    m.write_model2(ofile, tt);
    
  } // tt
  
}
