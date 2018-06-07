
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
#include "geometry.h"
#include "rhf1d.h"
#include "spectrum.h"
#include "pesc.h"

/* --- Function prototypes --                          -------------- */
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
extern Geometry geometry;

extern enum Topology topology;
extern rhbgmem *bmem;
int readBackground_j(int la, int mu, bool_t to_obs);

/* ------- file: -------------------------- voigt.c -----------------


--- Voigt profile generators:

       -- Armstrong 1967, JQSRT 7, pp. 61-88
          (slow for damping parameters larger than 1.5, accurate to
          6 significant figures).

       --                                              -------------- */

int searchLam(int n, double lamref, double *lam)
{
  register int la;
  for (la=0; la<n; la++){
    if(fabs(lamref-lam[la]) < 1.e-8){
      return la;
    }
  }
  return 0;
}


double phinorm(Atom *atom, int kr, int k, double lambda, double xmu)
{
  /* Line profile phi for line transition with 
     index kr 
     cos theta = xmu 
     wavelength lambda 
     depth point k
     Normalization is int phinorm  dfrequency [Hz] = 1.0
*/
  double x,shift,thermal,width,lab, bla;
  AtomicLine *trans = &atom->line[kr];

  
  lab=trans->lambda0;
  shift= vproject(k, 0);
  thermal = sqrt(2*atmos.T[k] * KBOLTZMANN / atom->weight/AMU) ;
  width=atmos.vturb[k];
  width = sqrt (thermal*thermal + width*width) ;  /* width in cm/s */

  x= lambda*(1. - xmu*shift/CLIGHT); /* wavelength in lab frame */
  x=(x/lab - 1.) / (width/CLIGHT)  ;    /* Dopp units */

  /* normalized, Can replace exp(-x*x) with Voigt */
  //bla = exp(-x*x)/SQRTPI/(width/lab);
  // fprintf(stderr,"normphi[%d]=%e -> %f %f\n", k,bla,lab, lambda);

  bla =  VoigtArmstrong(5.e-4,x)/SQRTPI/(width/lab) ;

  return bla;
}

double f1(double x)
{
  const char routineName[] = "f1";
  double f,gmma=0.5772156649;

  /* f1(x) = expint_1(x)*exp(x) */
  
  if(x <= 0.02) 
    {
      f=(-log(x)-gmma+x)*exp(x);
    }
  if(x >= 0.02 & x <= 1.5) 
    {
      f=log((x+1.)/x) - (0.36+0.03/sqrt(x+0.01))/(x+1.)/(x+1.) ;
    }
  if(x > 1.5 & x <= 10.) 
    {
      f=log((x+1.)/x)- (0.36+0.03*sqrt(x+0.01))/(x+1)/(x+1);
    }
  if(x > 10.) 
    {
      f = 1./x -1./x/x + 2./x/x/x ; 
    }
  
  return f ;
}

void uv_line( AtomicLine *trans, int kr, int k, int la, double xmu,
	 double *udn, double *vdn, double *vup)
{
  /* returns ull' and vll' of RH92  for lines, equations 2.2
     Inputs;
     kr = index of transition 
     k = atmospheric bin index
     lambda = wavelength at which calculation is needed 
     xmu = cosine angle
     output:  udn vdn vup as a function of spatial index k (uup=0)

     Note: udn, vdn refers to the downward transition of kr
     vup refers to the upward transition of kr
  */
  int i,indx,j,krp ;
  double HC4PI, phi;
  double beta, C2, Phi, xs;
  double lambda = trans->lambda[la];

  //fprintf(stderr,"line=%f, kr=%d\n", trans->lambda0, kr);
  
  HC4PI=HPLANCK*CLIGHT/4./PI/lambda ;
  phi=phinorm(trans->atom, kr, k, lambda, xmu) ;
  *udn = HC4PI  * phi * trans->Aji ;
  *vdn = HC4PI  * phi * trans->Bji ;
  *vup = HC4PI  * phi * trans->Bij ;
}


void uv_cont( AtomicContinuum *trans, int kr, int k, int la, double xmu,
	      double *udn, double *vdn, double *vup)
{
  int i,indx,j,krp ;
  double HC4PI, phi;
  double beta, C2, Phi, xs, lambda = trans->lambda[la];
  
  krp=kr;
  C2=2*HPLANCK*CLIGHT/lambda/lambda/lambda ;
  Linear(trans->Nlambda,trans->lambda, trans->alpha,1, &lambda, &xs, TRUE) ;
  j=trans->j;  /* j is upper level */
  i=trans->i;
  Phi=trans->atom->nstar[i][k]/trans->atom->nstar[j][k];   /* Phi here is Phi(RH92)/n_e */ 
  beta= HPLANCK*CLIGHT/KBOLTZMANN/lambda/atmos.T[k] ;
  *udn = Phi * C2 * exp(-beta) *xs ;
  *vdn = Phi * exp(-beta) *xs ;
  *vup = xs; 
  
}  

double pe(int bb, double t, double tc, double alpha, double *dq)
{
  double beta,b3, q=0.,tl,tmp,etc;
  *dq=0.;
  if(tc > 50.) return 0. ;
  etc=exp(-tc) ;
  if(bb == 1) { 
    beta=2.*PI;
    q=etc/(2.+ beta*t) ;
    *dq = -(tc*beta +2.*tc/t + beta)*q/(beta*t+2.) ;
  }
  else {
    beta = fmax(3.*(t+tc)/alpha,1.);
    b3=beta*beta*beta;
    q= exp(- b3*(t+tc) - alpha*(beta-1.))/2/beta;
    *dq = - b3 * q ;
  }
  return q;
} 

double escape(double lambda, double *chi, double *chib, double *S, double *P, double *Q, double *Psi, int bb)
{
  
  /*
    
    Given the input atmosphere, wavelength (line center or bf edge), opacity chi,
    chib,   S  source function, and the flag bb=1 for lines, bb=0 for continua, 
    
    this returns
    
    P  I+ + I- equiv J
    Q  I+ - I- equiv F
    Psi  approximate lambda operator
    output intensity returned upon call. 
    
  */
  
  const char routineName[] = "escape";
  register int k;
  int     Nspace = atmos.Nspace;
  double tau0 = 0.0, Bnu[2], zmu, Iplus,*tau,*tauc;
  double t,tc,dt; 
  double alpha;
  double sum;
  double zz,ss;
  double *z = geometry.height;
  
  ss=3.;  /* spectral index of continuum, nu^-ss */
  
  
  tau  = (double *) calloc(Nspace, sizeof(double));
  tauc  = (double *) calloc(Nspace, sizeof(double));
  
  /* frequency averaged approximations to escape probabilities  */  
  
  t=0.; tc=0. ; dt=0.; 
  zmu=1. ;
  
  for (k = 1;  k < Nspace-1;  k++) {
    zz=0.5*(z[k-1] - z[k+1])/zmu;
    t+=chi[k] * zz ;
    tc+=chib[k] * zz ;
    tau[k]=t+tc;
    tauc[k]=tc;
  }
  
  tau[0]=tau[1]*0.5;
  tauc[0]=tauc[1]*0.5;
  tau[Nspace-1]=2.*tau[Nspace-2] ;
  tauc[Nspace-1]=2.*tauc[Nspace-2];
  
  double h=0., ep=0., dp=0., factor=2.,dx=0.;
  k=Nspace-1 ;
  P[k]=S[k] ;
  Q[k]=0;
  Psi[k] = 1. ;
  sum=0.;

  
  for (k = Nspace-2;  k > 1 ;  k--) {
    
    t=tau[k] ; tc=tauc[k] ;
    
    alpha=HPLANCK*(CLIGHT/KBOLTZMANN)/lambda/atmos.T[k];
    ep=pe(bb,t,tc,alpha,&dp);
    
    Psi[k] = 1. - factor*ep ;    /* first order solution */ 

    /* orig linear 
       dt = (tau[k+1]+tauc[k+1] -tau[k-1]-tauc[k-1])/2.;  
       dt = t+tc - tau[k-1]-tauc[k-1];  
       dt = tau[k+1]+tauc[k+1] -t -tc;  
       h = -S[k] * dp *dt ; 
    */
    
    /*    dx = log( (tau[k+1]+tauc[k+1]) / (tau[k]+tauc[k]) ) ; */
    dx = log( (tau[k+1]+tauc[k+1]) / (tau[k-1]+tauc[k-1]) )/2. ;
    h = -S[k] * dp * (tau[k] * dx) ;  
    
    sum+=h;
    
    P[k]= S[k]*(1. - factor*ep) + sum ;   /* Formal Solution */
    Q[k]= -S[k]*factor*ep + sum ; 
    
  }
  
  P[0]=P[1] ;
  Psi[0]=Psi[1] ;
  Q[0]=Q[1] ;
  Iplus=P[0] ;
  
  
  free(tau), free(tauc);

  return Iplus;
}

void gamma_esc(Atom *atom)
{

  // Atom *atom;
  int i,j,k,l,la,kl,kr,krp,ii,ij,ji,indx,imu,nat,Nlevel=atom->Nlevel,Nspace,Nspectrum,bb, recnum;
  double lambda, integrand,xmu=1.0,wt,x;
  double chihm,dfreq,udnk,vdnk,vupk,ieff,HC=HPLANCK*CLIGHT;
  double Bnu[1],**accum_chi, **accum_u, **accum_v;
  double *chi,*s, *p,*psistar,*q, *chic,*etac,*tot_chi,*accum_eta, *z = geometry.height;
  double ss,*tau,*tauc,*chib,*sl ;
  double zmu=1.,phi;
  double zz=0.,t=0.,tc=0.,ie=0.;
  ActiveSet *as;

  AtomicLine *line = NULL;
  AtomicContinuum *cont = NULL;
 
  
  Nspace=atmos.Nspace;
  Nspectrum=spectrum.Nspect;
  //atmos->divf=(double *) malloc(Nspace * sizeof(double));
  chi=(double *) calloc(Nspace , sizeof(double));
  s=(double *) calloc(Nspace , sizeof(double));
  sl=(double *) calloc(Nspace , sizeof(double));
  p=(double *) calloc(Nspace , sizeof(double));
  q=(double *) calloc(Nspace , sizeof(double));
  psistar=(double *) calloc(Nspace , sizeof(double));
  chic=(double *) calloc(Nspace , sizeof(double));
  tau=(double *) calloc(Nspace , sizeof(double));
  tauc=(double *) calloc(Nspace , sizeof(double));
  chib=(double *) calloc(Nspace , sizeof(double));
  etac=(double *) calloc(Nspace , sizeof(double));
  tot_chi=(double *) calloc(Nspace , sizeof(double));
  accum_eta=(double *) calloc(Nspace , sizeof(double));
  
  accum_chi=matrix_double(Nlevel,Nspace);
  accum_u=matrix_double(Nlevel,Nspace);
  accum_v=matrix_double(Nlevel,Nspace);
  

  Nlevel = atom->Nlevel;
  for(kr=0;kr<atom->Nline; kr++){
    memset(atom->line[kr].Rij, 0, Nspace*sizeof(double));
    memset(atom->line[kr].Rji, 0, Nspace*sizeof(double));
  }
  for (k = 0;  k < Nspace;  k++) {
    for (i=1; i< Nlevel; i++) {
      for (j=0; j< i; j++) {
	ij = i*Nlevel + j;
	ji = j*Nlevel + i;
	atom->Gamma[ij][k] = 0.;
	atom->Gamma[ji][k] = 0.;
      }
    }
  }
  
  //--- Init gamma with collisional rate --- //
  
  initGammaAtom(atom, 1.0);


  // --- Main loop --- //
  
  for(kr=0;kr<atom->Nline; kr++){
    xmu = 1.0;
    zmu=1.,phi=0.0;
    zz=0.,t=0.,tc=0.,ie=0.;
    line = &atom->line[kr];
    bb = 1;
    
    i = line->i, j = line->j;
    kl=searchLam(line->Nlambda, line->lambda0, line->lambda);
    lambda = line->lambda[kl];
    ij = i*Nlevel + j;   
    ji = j*Nlevel + i;

    l = line->Nblue+kl;
    //
    as = spectrum.as + l;
    alloc_as(l, FALSE);
    
    readBackground_j(l, 0, 1);
    Opacity(l, 0, 1, 1);

    memcpy(chic, as->chi_c, Nspace*sizeof(double));
    memcpy(etac, as->eta_c, Nspace*sizeof(double));
    memcpy(chi, as->chi, Nspace*sizeof(double));

    
    for (k=0; k < Nspace ; k++){
      //void uv_line( AtomicLine *trans, int kr, int k, int la, double xmu,
      // double *udn, double *vdn, double *vup)
      uv_line(line, kr, k, kl, xmu, &udnk,&vdnk,&vupk);
      x=atom->n[i][k]*vupk - atom->n[j][k]*vdnk;

      zz=(z[k-1] - z[k+1])/zmu/2.;
      t+=chi[k] * zz ;
      tc+=chib[k] * zz ;
      tau[k]=t+tc;
      tauc[k]=tc;
      
      accum_chi[i][k]=x; /* absorption in the transition */
      accum_chi[j][k]=-x; /* negative absorption in the transition */
      accum_u[j][k]=udnk; /* emission */ 
      chi[k]=x;
      accum_v[j][k]=vdnk ; 
      accum_v[i][k]=vupk ; 
      accum_eta[k]=atom->n[j][k]*udnk;
      tot_chi[k]=chi[k]+chic[k] ; 
      s[k]=(accum_eta[k]+etac[k])/(chi[k] +chic[k]); 
      sl[k]=(accum_eta[k])/(chi[k]); 
      
    } // k

    spectrum.I[l][0] = escape( lambda, chi, chic, s, p, q, psistar, bb);
    
    for (k=0; k < Nspace ; k++){
      /*  lower to upper */
      ie=p[k] - s[k] * psistar[k]; 
      /*	    ie=fmax(ie,0);   */
      atom->Gamma[ji][k] +=line->Bij*ie ;
      /*   upper to lower */
      atom->Gamma[ij][k] += line->Aji*(1. - psistar[k])  + line->Bji*ie ;
      line->Rji[k]+=atom->Gamma[ji][k];
      line->Rij[k]+=atom->Gamma[ij][k];
      //phi=phinorm(atom, kr, k, lambda, xmu) ;
      //atmos->divf[k]+= q[k]*(chi[k]+chic[k])/phi;
    }
    
    free_as(l, FALSE);

  }// kr

  line = NULL;
  
  for(kr=0;kr<atom->Ncont; kr++){
    xmu = 1.0;
    zmu=1.,phi=0.0;
    zz=0.,t=0.,tc=0.,ie=0.;
    cont = &atom->continuum[kr];
    
    i = cont->i, j = cont->j;
    kl=0, bb=0;
    lambda = cont->lambda[kl];
    ij = i*Nlevel + j;   
    ji = j*Nlevel + i;

    l = cont->Nblue;
    
    as = spectrum.as + l;
    alloc_as(l, FALSE);
    
    readBackground_j(l, 0, 1);
    // Opacity(l, 0, 1, 1);

    memcpy(chic, as->chi_c, Nspace*sizeof(double));
    memcpy(etac, as->eta_c, Nspace*sizeof(double));
    //memcpy(chi, as->chi, Nspace*sizeof(double));


    
    for (k=0; k < Nspace ; k++){
      //void uv_line( AtomicLine *trans, int kr, int k, int la, double xmu,
      // double *udn, double *vdn, double *vup)
      uv_cont(cont, kr, k, kl, xmu, &udnk,&vdnk,&vupk);
      x=atom->n[i][k]*vupk - atom->n[j][k]*vdnk;

      zz=(z[k-1] - z[k+1])/zmu/2.;
      t+=chi[k] * zz ;
      tc+=chib[k] * zz ;
      tau[k]=t+tc;
      tauc[k]=tc;
      
      accum_chi[i][k]=x; /* absorption in the transition */
      accum_chi[j][k]=-x; /* negative absorption in the transition */
      accum_u[j][k]=udnk; /* emission */ 
      chi[k]=x;
      accum_v[j][k]=vdnk ; 
      accum_v[i][k]=vupk ; 
      accum_eta[k]=atom->n[j][k]*udnk;
      tot_chi[k]=chi[k]+chic[k] ; 
      s[k]=(accum_eta[k]+etac[k])/(chi[k] +chic[k]); 
      sl[k]=(accum_eta[k])/(chi[k]); 
      
    } // k

    spectrum.I[l][0] = escape( lambda, chi, chic, s, p, q, psistar, bb);
    
    
    double ss=3., ny0=0., alpha=0.,Blc=0.,gc=0.,Acl=0.;

      
    Blc=4*PI*cont->alpha0/HPLANCK/(ss-1.) ;
    ny0=(atom->E[j]-atom->E[i]) / HPLANCK ;
    for (k=0; k < Nspace ; k++){
      alpha=HPLANCK*ny0/KBOLTZMANN/atmos.T[k];
      gc=2.*PI*(M_ELECTRON/HPLANCK)*(KBOLTZMANN/HPLANCK)*atmos.T[k] ;
      gc*=sqrt(gc) *2.*atom->g[j]/atmos.ne[k];
      
      /*  lower to upper
	  /(ss+0) here means assumes J=J_0 shortward of continuum, 
	  s+k means J = J_0 nu^-k 
	  CRITICAL FACTOR PHOTOIONIZATION HERE ss+2, .8 is a "fit"
      */ 
      
      ie=(p[k] - s[k] * psistar[k])/(ss+.8);   
      
      atom->Gamma[ji][k] += Blc*ie ; 
      /* if(DEBUG > 2) {printf("%3d Blc  %9.2e %9.2e \n",k, Blc,ie); }*/
      
      /*   upper to lower  */
      Acl=8.*PI*atom->g[i] * (ny0*cont->alpha0) * 
	(ny0/CLIGHT)*(ny0/CLIGHT)*f1(alpha)/ gc ;
      atom->Gamma[ij][k] += Acl*(1. - psistar[k])  + Blc*ie/atom->g[i]/gc ; 
      cont->Rji[k]+=atom->Gamma[ji][k];
      cont->Rij[k]+=atom->Gamma[ij][k];
      
      //atmos->divf[k]+= q[k]*(chi[k]+chic[k])*ny0*(1. + 1./(ss+.8));
    }
    free_as(l, FALSE);

  }// kr





  // -- cleanup ---//
  free(chi), free(s), free(sl), free(p), free(q), free(psistar), free(chic),
    free(tau), free(tauc), free(chib), free(etac), free(tot_chi), free(accum_eta);

  freeMatrix((void**)accum_chi), freeMatrix((void**)accum_u), freeMatrix((void**)accum_v);
}

int pesc(Atom *atom, int it, double EMAX)
{
  register int i,j,k,ij,ji,lit;
  int Nlevel=atom->Nlevel, Nspace = atmos.Nspace;
  double **nold,**A,*B,change,maxchange,tt,xmu,wmu, cmult = 1.0;
  int imx,kmx=0;
  bool_t accel;
  
  nold = matrix_double(Nlevel, Nspace);
  A = matrix_double(Nlevel, Nlevel);
  B = (double *) calloc(Nlevel, sizeof(double));
  
  /*  GAMMA MATRIX FOR THIS ONE ATOM */
  
  gamma_esc(atom);
  
    
  maxchange=0. ;
  for(k=0;k<Nspace; k++){ /* main loop over space */
    /* set nold and initialize A to zero */
    for(i=0;i<Nlevel; i++){
      nold[i][k]=atom->n[i][k] ;
      for(j=0;j<Nlevel; j++) A[i][j]=0.;
    }
    
    for(i=0;i<Nlevel; i++){
      B[i]=0. ;
      for(j=0;j<Nlevel; j++){
	ij = i*Nlevel + j;
	ji = j*Nlevel + i;
	A[i][j]=atom->Gamma[ij][k] ;
	A[j][i]=atom->Gamma[ji][k] ;
      }
    }
    
    for(i=0;i<Nlevel; i++){
      A[i][i]=0.;
      for(j=0;j<Nlevel; j++){
	ji = j*Nlevel + i;
	A[i][i]-=atom->Gamma[ji][k] ;
      }
    }
    
    /* compute imax */
    int imax;
    double nmax;
    nmax=atom->n[0][k];
    imax=0;
    for (i=1;i<Nlevel;i++){
      if( nmax < atom->n[i][k]) { 
	nmax=atom->n[i][k] ; 
	imax=i;
      }
    }

    for(i=0;i<Nlevel; i++){
      A[imax][i]=1.0 ;  
    } 
    
    B[imax]=atom->ntotal[k];

    SolveLinearEq(Nlevel,A,B,TRUE);

    tt=0. ;
    for(i=0;i<Nlevel; i++) {
      if(B[i] <= 0.) { 
	printf("NEGATIVE POPULATION %d %d %e \n",k,i, B[i]);
      }
      atom->n[i][k]=B[i] ; /*update populations */
      tt+=atom->n[i][k];
      change=atom->n[i][k]/nold[i][k]-1 ;
      B[i]=0.;
      if(fabs(change) >= fabs(maxchange)) {
	maxchange=change ;
	imx=i;
	kmx=k;
      }
    }
    /*    if(DEBUG >= 1) printf("tt    %d %d %e  \n",i,j,tt); */
  }

  free(B), freeMatrix((void**)A), freeMatrix((void**)nold);
  
  //printf("[%3d] iteration max[%2d][%2d] %9.2e\n",it,kmx,imx,fabs(maxchange));
  if(fabs(maxchange) <= EMAX && cmult <= (1.0+EMAX/100.)){
    /* printf(" \n maxchange %e kmx %d imax %d\n", */
    /* 	   fabs(maxchange),kmx,imx);  */
    /* printf(" \n converged %s \n",&atom->ID);  */
    return 1;
  }
  
  return 0;
}


