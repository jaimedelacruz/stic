/* ------- file: -------------------------- initial_xdr.c -----------
1;95;0c
       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jan 16 17:29:31 2012 --

       --------------------------                      ----------RH-- */

/* --- Reads and/or computes the initial solution (populations and/or
       mean intensity J).

       XDR (external data representation) version.

       Possible options:

         LTE_POPULATIONS    -- Assume LTE populations initially.
         ZERO_RADIATION     -- Solve statistical equilibrium with
                               zero radiation field
         OLD_POPULATIONS    -- Read old populations from file
         OLD_POPS_AND_J     -- Read both old populations and J from file
         ESCAPE_PROBABILITY -- Not yet implemented
         OLD_J              -- Use mean intensities from previous solution
                               (Only implemented for wavelength_table).

       --                                              -------------- */


#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "initial_j.h"
#include "pesc.h"
//#include "mtime.h"

#define IMU_FILE_TEMPLATE "scratch/Imu_p%d.dat"

/* --- Function prototypes --                          -------------- */
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];

extern enum Topology topology;


void zeroRadiation(Atom *atom, int nact)
{

  static const double hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
  static const double twoc   = 2.0*CLIGHT / CUBE(NM_TO_M);
  static const double fourPI = 4.0 * PI;
  register int nspect, la, i, j, ij, k, n, kr;
  double wla, twohnu3_c2, gijk;
  AtomicLine *line;
  AtomicContinuum *continuum;
  ActiveSet *as;

  
  initGammaAtom(atom, 1.0);
  
  /* --- Then add radiative contributions of active transitions --  */

  for(kr=0;kr<atom->Nline; kr++){
    line = &atom->line[kr];
    i  = line->i;
    j  = line->j;
    ij = i*atom->Nlevel + j;
    
    for (k = 0;  k < atmos.Nspace;  k++)
      atom->Gamma[ij][k] += line->Aji;
  }
  
  
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    as = spectrum.as + nspect;
    
    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {    
      if(ATOMIC_CONTINUUM == as->art[nact][n].type){
	continuum = as->art[nact][n].ptype.continuum;
	//if(continuum->atom == atom){
	la = nspect - continuum->Nblue;
	i  = continuum->i;
	j  = continuum->j;
	ij = i*atom->Nlevel + j;
	
	wla = fourPI * getwlambda_cont(continuum, la) /
	  continuum->lambda[la];
	twohnu3_c2 = twoc / CUBE(continuum->lambda[la]);
	for (k = 0;  k < atmos.Nspace;  k++) {
	  gijk = atom->nstar[i][k]/atom->nstar[j][k] *
	    exp(-hc_k/(continuum->lambda[la] * atmos.T[k]));
	  atom->Gamma[ij][k] += gijk * twohnu3_c2 *
	    continuum->alpha[la]*wla;
	}
	//	}
      }
    }
  }
  /* --- Solve statistical equilibrium equations --  ------------ */
  
  statEquil(atom, (input.isum == -1) ? 0 : input.isum);
}

void initSolution_alloc2(int myrank) {
  const char routineName[] = "initSolution_p";
  register int nspect, nact,k,kr,mu,la;
  char    permission[3], file_imu[MAX_MESSAGE_LENGTH];
  int     Nsr, Nplane, index, oflag;
  int to_obs,sign,ncoef,ilow,Nlamu,lamu,Np;
  Molecule   *molecule;
  Atom       *atom;
  AtomicLine *line;
  long int idx, lc, lak, onc, omin, lamuk, Nlam;
  double *lambda,fac,lambda_prv,lambda_gas,lambda_nxt,dl,dl1,frac,lag;
  double q0,q_emit,qN, waveratio=1.0, t0=0, t1=0;
  static bool_t firsttime = true, firstl, lastl;
  static const double vsign[2] = {-1.0, 1.0};

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    atom->converged = FALSE;
    atom->mxchange = 1.e44;
  }
   

  /* --- Allocate space for angle-averaged mean intensity -- -------- */

  if (!input.limit_memory) {
    if (spectrum.J != NULL) freeMatrix((void **) spectrum.J);
    spectrum.J = matrix_double(spectrum.Nspect, atmos.Nspace);

    /* --- If we do background polarization we need space for the
           anisotropy --                               -------------- */

    if (input.backgr_pol){
      if (spectrum.J20 != NULL) freeMatrix((void **) spectrum.J20);
      spectrum.J20 = matrix_double(spectrum.Nspect, atmos.Nspace);
    }
  }

 /* --- For the PRD angle approximation  we need to store J in
     the gas frame,                                   -------- */
  if (input.PRD_angle_dep == PRD_ANGLE_APPROX &&  atmos.NPRDactive > 0) {
    
    //if (spectrum.Jgas != NULL) freeMatrix((void **) spectrum.Jgas);
    //spectrum.Jgas  = matrix_double(spectrum.Nspect, atmos.Nspace);
    if (spectrum.v_los != NULL) freeMatrix((void **) spectrum.v_los);
    spectrum.v_los = matrix_double(    atmos.Nrays, atmos.Nspace);     
    
    /* Calculate line of sight velocity */
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	spectrum.v_los[mu][k] = vproject(k, mu); // / vbroad[k];
      }
    }

    /* --- allocate space for gII redistribution function in PRD 
       lines                                                     --- */

    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      
      atom = atmos.activeatoms[nact];
      
      for (kr = 0;  kr < atom->Nline;  kr++) {
	
	line = &atom->line[kr];
	
	if (line->PRD) {
	  /*
	  
	  //line->gII  = (double **) malloc(atmos.Nrays * sizeof(double *));
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    for (la = 0 ;  la < line->Nlambda;  la++) {
	      
	      // second index in line.gII array
	      lak= k * line->Nlambda + line->Nlambda;
	      
	      q_emit = (line->lambda[la] - line->lambda0) * CLIGHT /
		(line->lambda0 * atom->vbroad[k]);
	      
	      
	      if (fabs(q_emit) < PRD_QCORE) {
		q0 = -PRD_QWING;
		qN =  PRD_QWING;
	      } else {
		if (fabs(q_emit) < PRD_QWING) {
		  if (q_emit > 0.0) {
		    q0 = -PRD_QWING;
		    qN = waveratio * (q_emit + PRD_QSPREAD);
		  } else {
		    q0 = waveratio * (q_emit - PRD_QSPREAD);
		    qN = PRD_QWING;
		  }
		} else {
		  q0 = q_emit - PRD_QSPREAD;
		  qN = q_emit + PRD_QSPREAD;
		}
	      }
	      Nsr = MAX((int) ((qN - q0) / PRD_DQ) + 1, Np);
	      
	      
	      //line->gII[lak]= (double *) malloc(Np*sizeof(double));
	      
	    }//la
	  }//k
	*/
	  
	  Nsr=3*MAX(3*PRD_QWING,2*PRD_QSPREAD)/PRD_DQ+1;
	  if (line->gII != NULL) freeMatrix((void **) line->gII);
	  line->gII  = matrix_double(atmos.Nspace*line->Nlambda,Nsr);
	  line->gII[0][0] = -1.0;
	}
      }
    }

  
   /* precompute prd_rho interpolation coefficients if requested */
    
    if (!input.prdh_limit_mem) {

      for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {

	atom = atmos.activeatoms[nact];

	for (kr = 0;  kr < atom->Nline;  kr++) {

	  line = &atom->line[kr];
	  
	  if (line->PRD) {
	      
	    Nlam = 2*atmos.Nrays * line->Nlambda;

	    if (line->frac != NULL) freeMatrix((void **) line->frac);
	    line->frac = matrix_float(Nlam, atmos.Nspace);

	    if (line->id0 != NULL) freeMatrix((void **) line->id0);
	    line->id0  = matrix_int(Nlam, atmos.Nspace);

	    if (line->id1 != NULL) freeMatrix((void **) line->id1);
	    line->id1  = matrix_int(Nlam, atmos.Nspace);

	    //if(line->Jgas != NULL) del_d2dim(line->Jgas, 0, 0);
	    //line->Jgas = d2dim(0, line->Nlambda-1, 0, atmos.Nspace-1);
	    
	    
	    for (la = 0;  la < line->Nlambda;  la++) {
	      for (mu = 0;  mu < atmos.Nrays;  mu++) {
		for (to_obs = 0;  to_obs <= 1;  to_obs++) {
		  sign = vsign[to_obs];//(to_obs) ? 1.0 : -1.0;
		  lamu = 2*(atmos.Nrays*la + mu) + to_obs;
		  
		  for (k = 0;  k < atmos.Nspace;  k++) {
		    
		    // wavelength in local rest frame 
		    lag=line->lambda[la] * (1.+spectrum.v_los[mu][k]*sign/CLIGHT);
		    		    
		    if (lag <= line->lambda[0]) {
		      // out of the lambda table, constant extrapolation
		      line->frac[lamu][k]=0.0;
		      line->id0[lamu][k]=0;
		      line->id1[lamu][k]=1;
		    } else if (lag >= line->lambda[line->Nlambda-1] ) {
		      // out of the lambda table, constant extrapolation
		      line->frac[lamu][k]=1.0;
		      line->id0[lamu][k]=line->Nlambda-2;
		      line->id1[lamu][k]=line->Nlambda-1;
		    } else {
		      // Locate index of line->lambda of point directly to the left of lag
		      Locate(line->Nlambda,line->lambda,lag,&ilow);
		      line->frac[lamu][k] = (float)(lag-line->lambda[ilow])/ (line->lambda[ilow+1]-line->lambda[ilow]);
		      line->id0[lamu][k]=ilow;
		      line->id1[lamu][k]=ilow+1;
		    }
		  }
		}
	      }
	    }
	  } // PRD line
	} // kr
      } // atoms
      
      // --- allocate coefficients for Jgas --- //
      
      Nlamu=2*2*atmos.Nrays*(spectrum.nJlam)*atmos.Nspace;
      
      if(spectrum.iprdh != NULL) free((void*)spectrum.iprdh);
      spectrum.iprdh =  ((unsigned int *)calloc( Nlamu , sizeof(unsigned int   )));
      
      if(spectrum.cprdh != NULL) free((void*)spectrum.cprdh);
      spectrum.cprdh =  ((unsigned int *)calloc( Nlamu , sizeof(unsigned int)));

      Nlamu=2*2*atmos.Nrays*(spectrum.nJlam+2)*atmos.Nspace;
      if (spectrum.nc != NULL) free((void*)(spectrum.nc-1));
      spectrum.nc=  ((unsigned int*)calloc(Nlamu, sizeof(unsigned int)))+1;
      
      lc = 0, onc = 0, omin = 0;
      lambda = spectrum.Jlam;
      
      for (la = 0;  la < spectrum.nJlam;  la++) {
	for (mu = 0;  mu < atmos.Nrays;  mu++) {
	  firstl = true;
	  
	  for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	    sign = vsign[to_obs];//(to_obs) ? 1.0 : -1.0;
	    //lamu = 2*(atmos.Nrays*la + mu) + to_obs;
	    
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      lastl = true;
	      lamuk = la * (atmos.Nrays*2*atmos.Nspace) 
		+ mu     * (2*atmos.Nspace)
		+ to_obs * (atmos.Nspace)
		+ k ;
	      
	      ncoef=0;
	      lc = spectrum.nc[lamuk-1];
	      
	      // previous, current and next wavelength shifted to gas rest frame 
	      fac = (1.+spectrum.v_los[mu][k]*sign/CLIGHT);
	      lambda_prv = lambda[ la-1                 ]*fac;
	      lambda_gas = lambda[ la                   ]*fac;
	      lambda_nxt = lambda[ la+1                 ]*fac;
	      
	      //omin = 0;		
		
	      // do lambda_prv and lambda_gas bracket lambda points?
	      dl = (lambda_gas - lambda_prv); // It contains the conversion factor to ushort (65535)
	      dl1= (lambda_nxt - lambda_gas); // It contains the conversion factor to ushort (65535)

	      lastl = true;
	      for (idx = 0; idx < spectrum.nJlam ; idx++) {
		
		if ((lambda[idx] > lambda_prv) && (lambda[idx] <= lambda_gas)){
		  
		  ncoef += 1;
		  spectrum.iprdh[lc]=(unsigned int)idx;
		  spectrum.cprdh[lc++]=(unsigned int)(round(65535.*(lambda[idx]-lambda_prv)/dl));
		  /* if(firstl){ */
		  /*   omin = min(idx, omin); */
		  /*   firstl = false; */
		  /* } */
		  lastl = false; 
		  
		  continue;
		  //omin = min(omin,idx);
		  //lc++;
		}

		if(!lastl) break;

	      }
	      
	      lastl = true;
	      for (idx = 0; idx < spectrum.nJlam ; idx++) {
		
		if ((lambda[idx] > lambda_gas) && (lambda[idx] < lambda_nxt)){
		  ncoef += 1;
		  spectrum.iprdh[lc]=(unsigned int)idx;
		  spectrum.cprdh[lc++]=(unsigned int)(round(65535.*(1. - (lambda[idx]-lambda_gas)/dl1))); // The 255 must be 1.0 if using floats
		  
		  /* if(firstl){ */
		  /*   omin = min(idx, omin); */
		  /*   firstl = false; */
		  /* } */
		  lastl = false;
		  
		  
		  continue;
		  //lc++;
		}
		
		if(!lastl) break;
		
	      } // idx
	      
	      spectrum.nc[lamuk]=onc+ncoef;
	      onc = spectrum.nc[lamuk];
	      //      }
	      
	      // lc = spectrum.nc[lamuk-1];
	      
	    } // k
	  } // to_obs
	} // mu
	//	omin = max(omin-1, 0);
	
      }// la
      
      
    }
    
    /* /\* precompute Jgas interpolation coefficients if requested *\/ */
    /* if (!input.prdh_limit_mem) { */

    /*   lambda = spectrum.lambda; */
    /*   idx=5*atmos.Nrays*(spectrum.Nspect+2)*atmos.Nspace; */
      
    /*   /\* --- keeps track of where to get indices and interpolation */
    /*          coefficients in spectrum.iprhh and spectrum.cprdh --- *\/ */
    /*   if (spectrum.nc != NULL) free((void*)(spectrum.nc-1)); */
    /*   spectrum.nc=  ((int*)calloc(idx, sizeof(int)))+1; */

    /*   /\* --- now we know the number of interpolation coefficients, */
    /* 	 it's stored in the last element of spectrum.nc, */
    /* 	 so allocate space                                     --- *\/ */
      
    /*   //  idx=spectrum.nc[2*atmos.Nrays*spectrum.Nspect*atmos.Nspace-1]; */
    /*   // On average the total number of coeffs must be 2 per wavelength */

    /*   // --- we do this trick to have two points at infinity at the extremes --- // */

      
      
      
    /*   if (spectrum.iprdh != NULL) free((void*)(spectrum.iprdh-1)); */
    /*   spectrum.iprdh= ((int *)calloc( idx , sizeof(int   )))+1;  */

    /*   if (spectrum.cprdh != NULL) free((void*)(spectrum.cprdh-1)); */
    /*   spectrum.cprdh= ((unsigned char*)calloc( idx, sizeof(unsigned char)))+1; */
      

    /*   //t0 = gettime(); */

    /*   lc = 0, onc = 0, omin = 0; */
    /*   for (la = 0;  la < spectrum.Nspect;  la++) { */
    /* 	for (mu = 0;  mu < atmos.Nrays;  mu++) { */
    /* 	  firstl = true; */

    /* 	  for (to_obs = 0;  to_obs <= 1;  to_obs++) { */
    /* 	    sign = (to_obs) ? 1.0 : -1.0; */
    /* 	    for (k = 0;  k < atmos.Nspace;  k++) { */
    /* 	      lastl = true; */

	      
    /* 	      lamuk = la * (atmos.Nrays*2*atmos.Nspace)  */
    /* 		+ mu     * (2*atmos.Nspace) */
    /* 		+ to_obs * (atmos.Nspace) */
    /* 		+ k ; */
	      
    /* 	      ncoef=0; */
	      
    /* 	      // previous, current and next wavelength shifted to gas rest frame  */
    /* 	      fac = (1.+spectrum.v_los[mu][k]*sign/CLIGHT); */
    /* 	      lambda_prv = lambda[ la-1                 ]*fac; */
    /* 	      lambda_gas = lambda[ la                   ]*fac; */
    /* 	      lambda_nxt = lambda[ la+1                 ]*fac; */
	      
    /* 	      dl = 255./(lambda_gas - lambda_prv); // It contains the conversion factor to ushort (65535) */
    /* 	      dl1= 255./(lambda_nxt - lambda_gas); // It contains the conversion factor to ushort (65535) */
	      
    /* 	      for (idx = omin; idx < spectrum.Nspect ; idx++) { */
		
		
    /* 		// do lambda_prv and lambda_gas bracket lambda points? */
    /* 		if ((lambda[idx] > lambda_prv) && (lambda[idx] <= lambda_gas)){ */
    /* 		  ncoef += 1; */
    /* 		  spectrum.iprdh[lc]=idx; */
    /* 		  spectrum.cprdh[lc++]=(unsigned char)(round((lambda[idx]-lambda_prv)*dl)); */
		  
    /* 		  if(firstl){ */
    /* 		    omin = min(idx, omin); */
    /* 		    firstl = false; */
    /* 		  } */
    /* 		  lastl = false; */
		  
    /* 		  continue; */
    /* 		  //omin = min(omin,idx); */
    /* 		  //lc++; */
    /* 		}else if ((lambda[idx] > lambda_gas) && (lambda[idx] < lambda_nxt)){ */
    /* 		  ncoef += 1; */
    /* 		  spectrum.iprdh[lc]=idx; */
    /* 		  spectrum.cprdh[lc++]=(unsigned char)(round(255. - (lambda[idx]-lambda_gas)*dl1)); // The 255 must be 1.0 if using floats */
    /* 		  if(firstl){ */
    /* 		    omin = min(idx, omin); */
    /* 		    firstl = false; */
    /* 		  } */
    /* 		  lastl = false; */

		  
    /* 		  continue; */
    /* 		  //lc++; */
    /* 		} */

    /* 		if(!lastl) break; */
		
    /* 	      } // idx */

	      
      /* 	      /\* // do lambda_prv and lambda_gas bracket lambda points? *\/ */
      /* 	      /\* if (lambda_prv !=  lambda_gas) { *\/ */
      /* 	      /\* 	dl= lambda_gas - lambda_prv; *\/ */
      /* 	      /\* 	for (idx = 0; idx < spectrum.Nspect ; idx++) {      *\/ */
      /* 	      /\* 	  if (lambda[idx] > lambda_prv && lambda[idx] <= lambda_gas) ncoef=ncoef+1; *\/ */
      /* 	      /\* 	} *\/ */
      /* 	      /\* } else { *\/ */
      /* 	      /\* 	// edge case, use constant extrapolation for lambda[idx]<lambda gas *\/ */
      /* 	      /\* 	for (idx = 0; idx < spectrum.Nspect ; idx++) { *\/ */
      /* 	      /\* 	  if (lambda[idx] <=  lambda_gas) ncoef=ncoef+1; *\/ */
      /* 	      /\* 	} *\/ */
      /* 	      /\* }  *\/ */
	      
      /* 	      /\* // do lambda_gas and lambda_nxt bracket lambda points? *\/ */
      /* 	      /\* if (lambda_gas != lambda_nxt) { *\/ */
      /* 	      /\* 	dl= lambda_nxt - lambda_gas; *\/ */
      /* 	      /\* 	for (idx = 0; idx < spectrum.Nspect ; idx++) {      *\/ */
      /* 	      /\* 	  if (lambda[idx] > lambda_gas && lambda[idx] < lambda_nxt) ncoef=ncoef+1; *\/ */
      /* 	      /\* 	} *\/ */
      /* 	      /\* } else { *\/ */
      /* 	      /\* 	// edge case, use constant extrapolation for lambda[idx]>lambda gas *\/ */
      /* 	      /\* 	for (idx = 0; idx < spectrum.Nspect ; idx++) { *\/ */
      /* 	      /\* 	  if (lambda[idx] >=  lambda_gas) ncoef=ncoef+1; *\/ */
      /* 	      /\* 	} *\/ */
      /* 	      /\* }  *\/ */

      /* 	      /\* --- number of point this lambda contributes to is */
      /* 	 	 computed as a difference --- *\/ */
      /* 	      //fprintf(stderr,"lamuc=%ld\n", lamuk); */
      /* 	      //      if (lamuk == 0) { */
      /* 	      //	spectrum.nc[lamuk] = ncoef; */
      /* 		//} else { */
      /* 	      //	      spectrum.nc[lamuk]=spectrum.nc[lamuk-1]+ncoef; */
      /* 	      spectrum.nc[lamuk]=onc+ncoef; */
      /* 	      onc = spectrum.nc[lamuk]; */
      /* 		//      } */
	      
      /* 	      lc = spectrum.nc[lamuk-1]; */

      /* 	    } // k */
      /* 	  } // to_obs */
      /* 	} // mu */
      /* 	omin = max(omin-1, 0); */
      /* 	//	fprintf(stderr, "la=%ld, omin=%ld, max=%ld\n", la, omin,spectrum.Nspect ); */
      /* } // la */
      //t1 = gettime();
      //fprintf(stderr,"DT third loop=%fs\n", t1-t0);
      //t0=t1;
    //  idx =  5*atmos.Nrays*(spectrum.Nspect+2)*atmos.Nspace;
    //   if((idx-spectrum.nc[2*atmos.Nrays*spectrum.Nspect*atmos.Nspace-1])<0){
       
    //	 fprintf(stderr,"allocated=%ld, needed=%ld, difference=%ld\n",idx, spectrum.nc[2*atmos.Nrays*spectrum.Nspect*atmos.Nspace-1], idx-spectrum.nc[2*atmos.Nrays*spectrum.Nspect*atmos.Nspace-1]);
    //	 exit(0);
       

     /* --- Run through all lamuk points again, and now store indices
            to lambda array and the corresponding interpolation
            coefficients                                          --- */
      /* for (la = 0;  la < spectrum.Nspect;  la++) { */
      /* 	for (mu = 0;  mu < atmos.Nrays;  mu++) { */
      /* 	  for (to_obs = 0;  to_obs <= 1;  to_obs++) { */

      /* 	    sign = (to_obs) ? 1.0 : -1.0; */

      /* 	    for (k = 0;  k < atmos.Nspace;  k++) { */

      /* 	      lamuk = la * (atmos.Nrays*2*atmos.Nspace)  */
      /* 		+ mu     * (2*atmos.Nspace) */
      /* 		+ to_obs * (atmos.Nspace) */
      /* 		+ k ; */

      /* 	      // starting index for storage for this lamuk point */
      /* 	      lc = (lamuk==0) ? 0 : spectrum.nc[lamuk-1]; */
	 	      
      /* 	      // previous, current and next wavelength shifted to gas rest frame  */
      /* 	      fac = (1.+spectrum.v_los[mu][k]*sign/CLIGHT); */
      /* 	      lambda_prv = lambda[ MAX(la-1,0)                 ]*fac; */
      /* 	      lambda_gas = lambda[ la                          ]*fac; */
      /* 	      lambda_nxt = lambda[ MIN(la+1,spectrum.Nspect-1) ]*fac; */
	      
      /* 	      // do lambda_prv and lambda_gas bracket lambda points? */
      /* 	      if (lambda_prv !=  lambda_gas) { */
      /* 		dl= lambda_gas - lambda_prv; */
      /* 		for (idx = 0; idx < spectrum.Nspect ; idx++) {      */
      /* 		  if (lambda[idx] > lambda_prv && lambda[idx] <= lambda_gas) { */
      /* 		    // bracketed point found */
      /* 		    spectrum.iprdh[lc]=idx; */
      /* 		    spectrum.cprdh[lc]=(lambda[idx]-lambda_prv)/dl; */
      /* 		    lc++; */
      /* 		  } */
      /* 		} */
      /* 	      } else { */
      /* 		// edge case, use constant extrapolation for lambda[idx]<lambda gas */
      /* 		for (idx = 0; idx < spectrum.Nspect ; idx++) { */
      /* 		  if (lambda[idx] <=  lambda_gas)  {   */
      /* 		    spectrum.iprdh[lc]=idx; */
      /* 		    spectrum.cprdh[lc]=1.0; */
      /* 		    lc++; */
      /* 		  } */
      /* 		} */
      /* 	      }  */
		
      /* 	      // do lambda_gas and lambda_nxt bracket lambda points? */
      /* 	      if (lambda_gas != lambda_nxt) { */
      /* 		dl= lambda_nxt - lambda_gas; */
      /* 		for (idx = 0; idx < spectrum.Nspect ; idx++) {      */
      /* 		  if (lambda[idx] > lambda_gas && lambda[idx] < lambda_nxt) { */
      /* 		    // bracketed point found */
      /* 		    spectrum.iprdh[lc]=idx; */
      /* 		    spectrum.cprdh[lc]=1.0 - (lambda[idx]-lambda_gas)/dl; */
      /* 		    lc++; */
      /* 		  } */
      /* 		} */
      /* 	      } else { */
      /* 		// edge case, use constant extrapolation for lambda[idx]>lambda gas */
      /* 		for (idx = 0; idx < spectrum.Nspect ; idx++) { */
      /* 		  if (lambda[idx] >=  lambda_gas)  { */
      /* 		    spectrum.iprdh[lc]=idx; */
      /* 		    spectrum.cprdh[lc]=1.0; */
      /* 		    lc++; */
      /* 		  } */
      /* 		} */
      /* 	      }  */
	   
      /* 	    } // k */
      /* 	  } // to_obs */
      /* 	} // mu */
      /* } // la */

    //  }  //input.prdh_limit_mem if switch      
  } // PRD_ANGLE_APPROX if switch

  // t1 = gettime();
  // fprintf(stderr,"DT fourth loop=%fs\n", t1-t0);
  // t0=t1;
  /* --- Allocate space for the emergent intensity --  -------------- */

  if (atmos.Stokes || input.backgr_pol) {
    if (spectrum.Stokes_Q != NULL) {
      freeMatrix((void **) spectrum.Stokes_Q);
      spectrum.Stokes_Q = NULL;
    }
    if (spectrum.Stokes_U != NULL) {
      freeMatrix((void **) spectrum.Stokes_U);
      spectrum.Stokes_U = NULL;
    }
    if (spectrum.Stokes_V != NULL) {
      freeMatrix((void **) spectrum.Stokes_V);
      spectrum.Stokes_V = NULL;
    }
  }
  
  switch (topology) {
  case ONE_D_PLANE:
    if (spectrum.I != NULL) freeMatrix((void **) spectrum.I);
    spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (atmos.Stokes || input.backgr_pol) {
      spectrum.Stokes_Q = matrix_double(spectrum.Nspect, atmos.Nrays);
      spectrum.Stokes_U = matrix_double(spectrum.Nspect, atmos.Nrays);
      spectrum.Stokes_V = matrix_double(spectrum.Nspect, atmos.Nrays);
    }
    break;
  case TWO_D_PLANE:
    Nsr = spectrum.Nspect * atmos.Nrays;
    spectrum.I = matrix_double(Nsr, atmos.N[0]);
    if (atmos.Stokes || input.backgr_pol) {
      spectrum.Stokes_Q = matrix_double(Nsr, atmos.N[0]);
      spectrum.Stokes_U = matrix_double(Nsr, atmos.N[0]);
      spectrum.Stokes_V = matrix_double(Nsr, atmos.N[0]);
    }
    break;
  case THREE_D_PLANE:
    spectrum.I = matrix_double(spectrum.Nspect * atmos.Nrays, 
			       atmos.N[0] * atmos.N[1]);
    if (atmos.Stokes || input.backgr_pol) {
      Nsr    = spectrum.Nspect * atmos.Nrays;
      Nplane = atmos.N[0] * atmos.N[1];

      spectrum.I = matrix_double(Nsr, Nplane);
      if (atmos.Stokes || input.backgr_pol) {
	spectrum.Stokes_Q = matrix_double(Nsr, Nplane);
	spectrum.Stokes_U = matrix_double(Nsr, Nplane);
	spectrum.Stokes_V = matrix_double(Nsr, Nplane);
      }
    }
    break;
  case SPHERICAL_SYMMETRIC:
    spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (atmos.Stokes) {
      Error(ERROR_LEVEL_2, routineName,
	    "Cannot do a full Stokes solution in spherical geometry");
    }    
    break;
  default:
    sprintf(messageStr, "Unknown topology (%d)", topology);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }



  /* Things to be done only for the first task */
  if (firsttime) {
    /* --- Need storage for angle-dependent specific intensities for
       angle-dependent PRD --                        -------------- */

    if (atmos.NPRDactive > 0 && input.PRD_angle_dep == PRD_ANGLE_DEP) {
      oflag = 0;
      if (input.startJ == OLD_J) {
	if (spectrum.updateJ) {
	  strcpy(permission, "r+");
	  oflag |= O_RDWR;
	} else {
	  strcpy(permission, "r");
	  oflag |= O_RDONLY;
	}
      } else {
	strcpy(permission, "w+");
	oflag |= (O_RDWR | O_CREAT);
      }
      /* Imu file name, this may not work very well in the mpi version... */
      sprintf(file_imu, IMU_FILE_TEMPLATE, myrank);
      
      if ((spectrum.fd_Imu = open(file_imu, oflag, PERMISSIONS)) == -1) {
	sprintf(messageStr, "Unable to open %s file %s with permission %s",
		(spectrum.updateJ) ? "update" : "input",
		file_imu, permission);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      /* --- Fill the index list that keeps track of the location
	 of intensity Imu in file spectrum.fd_Imu at wavelength
	 corresponding to nspect. --                 -------------- */
      
      spectrum.PRDindex = (int *) malloc(spectrum.Nspect * sizeof(int));
      index = 0;
      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
	if (containsPRDline(&spectrum.as[nspect])) {
	  spectrum.PRDindex[nspect] = index++;
	}
      }
    }


    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];
      
      /* --- Allocate memory for the rate equation matrix -- ---------- */
      atom->Gamma = matrix_double(SQ(atom->Nlevel), atmos.Nspace);
    }


    for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
      molecule = atmos.activemols[nact];
      
      /* --- Allocate memory for the rate equation matrix -- ---------- */
      molecule->Gamma = matrix_double(SQ(molecule->Nv), atmos.Nspace);
    }

    firsttime = false;
  } /* End of first task condition */



  return;
}



void initSolution_j( int myrank, int savepop)
{
  const char routineName[] = "initSolution_j";
  register int k, i, ij, nspect, n, nact;
  int     la, j, status;
  double  gijk, wla, twohnu3_c2, hc_k, twoc, fourPI;
  ActiveSet  *as;
  Molecule   *molecule;
  Atom       *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  int conv = 0;
  getCPU(2, TIME_START, NULL);
  

  /* --- Allocate arrays, taken fro T. Pereira's 
     version to not mess up with the PRD business --- */
  
  initSolution_alloc2(myrank);
  

  /* --- Fill matrix J with old values from previous run ----- -- */
  if (input.startJ == OLD_J){
    fprintf(stderr,"initSolution_j: OLD_J not supported for inversions, exiting!");
    exit(0);
  }

  //if(!savepop) return;// Reusing departure coeffs
  
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    /* --- Initialize the mutex lock for the operator Gamma if there
           are more than one threads --                -------------- */

    if (input.Nthreads > 0) {
      if ((status = pthread_mutex_init(&atom->Gamma_lock, NULL))) {
	sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		status);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    switch(atom->initial_solution) {
    case LTE_POPULATIONS:
      for (i = 0;  i < atom->Nlevel;  i++) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  atom->n[i][k] = atom->nstar[i][k];
      }
      break;

    case ZERO_RADIATION:
      zeroRadiation(atom, nact);
	    
      break;

    case PESC:
      //readPopulations(atom);
      zeroRadiation(atom, nact);
      for(i=0; i<100; i++){
	conv = pesc(atom, i, 1.e-2);
	//if(i>2)exit(0);
	if(conv) break;
      }
      if(!conv) zeroRadiation(atom, nact);
      break;

    default:;
    break;
    }
  }
  /* --- Now the molecules that are active --          -------------- */
  
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    /* --- Calculate the LTE vibration level populations here. They
           cannot be calculated yet in readMolecule since chemical
           equilibrium has to be established first --  -------------- */

    for (i = 0;  i < molecule->Nv;  i++) {
      for (k = 0;  k < atmos.Nspace;  k++)
	molecule->nvstar[i][k] = molecule->n[k] *
	  molecule->pfv[i][k] / molecule->pf[k];
    }


    /* --- Initialize the mutex lock for the operator Gamma if there
           are more than one thread --                 -------------- */

    if (input.Nthreads > 0) {
      if ((status = pthread_mutex_init(&molecule->Gamma_lock, NULL))) {
	sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		status);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    switch(molecule->initial_solution) {

    case LTE_POPULATIONS:
      for (i = 0;  i < molecule->Nv;  i++) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  molecule->nv[i][k] = molecule->nvstar[i][k];
      }
      break;
      
    case OLD_POPULATIONS:
      //readMolPops(molecule);
      break;

    default:;
    }

    /* --- Calculate collisions for molecule (must be done here because
           rotation-vibration transitions are dominated by hydrogen and
           H2 collisions for which chemical equilibrium needs to be
           established first --                        -------------- */

    if (strstr(molecule->ID, "CO"))
      COcollisions(molecule);
    else if (strstr(molecule->ID, "H2"))
      H2collisions(molecule);
    else {
      sprintf(messageStr, "Collisions for molecule %s not implemented\n",
	      molecule->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
}
/* ------- end ---------------------------- initSolution.c ---------- */


