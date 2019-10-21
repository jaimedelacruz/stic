/* ----
   WITTMANN EOS routines, adapted from the SIR code implementation to C++ templates
   J. de la Cruz Rodriguez (ISP-SU). Original reference Mihalas (1970).
   Added proper treatment of rho compared to the original routines.

   Added interface routines to the background opacity code.
   
   2017-02-28, JdlCR: Added partition functions from Uppsala and a saha solver
                      that includes more than 3 levels if the partition functions
		      are known
		     
   2017-08-18, JdlCR: Added proper treatment of rho: it accounts that H2 molecules 
                      reduce the amount of neutral and ionized Hydrogen and therefore
		      it must be accounted for in rho.

   ---- */

#include <cmath>
#include <string>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>      // std::ifstream

//
#include "physical_consts.h"
#include "witt.h"
#include "partition.h"
#include "eoswrap.h"
#include "input.h"

using namespace std;
using namespace eos;
//

/* --------------------------------------------------------------------------------------- */

template <class T> T pow10(T x){return (T)exp(2.3025850929940459*x);};

/* --------------------------------------------------------------------------------------- */

const float witt::prec = 1.e-5; 
const int   witt::ncontr = 28;
const double witt::saha_fac =  pow( (2.0 * phyc::PI * phyc::ME * phyc::BK) / (phyc::HH*phyc::HH), 1.5);

/* --------------------------------------------------------------------------------------- */

template <class T> void witt::acota(T &x, T x0, T x1)
{
  if(x < x0) x = x0;
  if(x > x1) x = x1;
}

template void witt::acota(float &, float, float);
template void witt::acota(double &, double, double);

/* --------------------------------------------------------------------------------------- */

//
template <class T> void witt::acotasig(T &x, T x0, T x1){
  if(x < 0){
    x = -x;
    acota<T>(x, x0, x1);
    x = -x;
  }else acota<T>(x, x0, x1);
}

template void witt::acotasig(float &, float, float);
template void witt::acotasig(double &, double, double);


/* --------------------------------------------------------------------------------------- */

eos::witt::witt(std::vector<line_t> &lines, std::string &abFile, double grav):
  avweight(0.0), gravity(pow10<double>(grav)), eoswrap()
{

  /* --- Fill array with default ABUNDance --- */
  
  //ABUND.resize(Nelements,0.0);
  memcpy(&ABUND[0], &eos::ABUND_default[0], 99*sizeof(float));
  
  
  /* --- if n >0 then copy given ABUNDances to array --- */
  
  //if(n > 0) memcpy(&ABUND[0], ab, n*sizeof(float));


  /* --- Compute sums --- */
  
  avweight = 0.0, tABUND = 0.0, Ab_others = 0.0, muH = 0.0, rho_from_H = 0.0;
  
  for(int ii=0; ii<Nelements; ii++){
    ABUND[ii] = pow10<double>(ABUND[ii]);
    ABUND[ii] /= ABUND[0]; // Make X_H = 1.0
    
    if(ii > 0) Ab_others += ABUND[ii]; // Used in gasc and pe_pg
    
    avweight += ABUND[ii] * phyc::AMASS[ii];
    tABUND += ABUND[ii];
  }

  // mu'= sum amass_i * ABUND_i / amass_0 
  muH = avweight / phyc::AMASS[0] / ABUND[0]; //avweight is here the total sum, until we normalize
  rho_from_H = muH * phyc::AMASS[0]*phyc::AMU / phyc::BK;
  
  /* --- Now compute average particle weight --- */
  
  avweight = avweight / tABUND * phyc::AMU;


  /* --- Init species table --- */
  
  // init_Species_table(lines);
  readAbund(abFile);
  
  
}

/* --------------------------------------------------------------------------------------- */ 

template <class T> T eos::witt::init_pe_from_T_Pg(T t, T pg)
{
  T nu,saha,aaa,bbb,ccc,ybh,pe;
    	
  nu=0.9091;//       ! assume that only Hydrogen is ionized
  saha=-0.4771+2.5*log10(t)-log10(pg)-(13.6*5040./t);
  saha=pow10<T>(saha);
  aaa=1.0+saha;
  bbb=-(nu-1.)*saha;
  ccc=-saha*nu;
  ybh=(-bbb+sqrt(bbb*bbb-4.*aaa*ccc))/(2.*aaa); //! ionization fraction
  pe=pg*ybh/(1.+ybh);
    
  return pe;
}

template float eos::witt::init_pe_from_T_Pg(float, float);
template double eos::witt::init_pe_from_T_Pg(double, double);

/* --------------------------------------------------------------------------------------- */

template <class T> void witt::getBackgroundPartials(T t, T Pgas, T Pe, double *n)
{

  T  pp[5], u = 0;
  int nLev = 0;
  T xpa[6], pf[6], ein[6];
  
  /* --- Use more acurate PF from Atlas-9, including up to 6 levels --- */

  /* --- He/He+/He++ --- */
  
  nLev = getXpart<T>(1, t, Pgas, Pe, xpa, pf, ein, true);
  n[3] = xpa[0], n[4] = xpa[1], n[5] = xpa[2];


  /* --- C --- */
  
  nLev = getXpart<T>(5, t, Pgas, Pe, xpa, pf, ein, true);
  n[6] = xpa[0];

  
  /* --- Al --- */

  nLev = getXpart<T>(12, t, Pgas, Pe, xpa, pf, ein, true);
  n[7] = xpa[0];


  /* --- Si/Si+ --- */

  nLev = getXpart<T>(13, t, Pgas, Pe, xpa, pf, ein, true);
  n[8] = xpa[0], n[9] = xpa[1];
  

  /* --- Ca/Ca+ --- */

  nLev = getXpart<T>(19, t, Pgas, Pe, xpa, pf, ein, true);
  n[10] = xpa[0], n[11] = xpa[1];

  /* --- Mg/Mg+ --- */

  nLev = getXpart<T>(11, t, Pgas, Pe, xpa, pf, ein, true);
  n[12] = xpa[0], n[13] = xpa[1];


  /* --- Fe --- */

  nLev = getXpart<T>(25, t, Pgas, Pe, xpa, pf, ein, true);
  n[14] = xpa[0];


  /* --- N --- */
  
  nLev = getXpart<T>(6, t, Pgas, Pe, xpa, pf, ein, true);
  n[15] = xpa[0];


  /* --- O --- */

  nLev = getXpart<T>(7, t, Pgas, Pe, xpa, pf, ein, true);
  n[16] = xpa[0];
  
  
  /* --- Now get H- --- */
  
  T pg = Pgas;
  gasc<T>(t, Pe, pg, pp); 
  
  
  n[0]  = pp[0]*pp[4]/(t*phyc::BK) * 0.5; // H / pf[H]
  n[1]  = pp[1]*pp[4]/(t*phyc::BK);     // H+ / 1.0
  n[2]  = pp[3]*pp[4]/(t*phyc::BK);     // H- / 1.0
  
  
}

template  void witt::getBackgroundPartials(float, float, float, double*);
template  void witt::getBackgroundPartials(double, double, double, double*);

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::nsaha(T t, T xne, T u0, T u1, T eion)
{
  double fexparg = exp(- double(eion) / double(t*phyc::BK));
  double res = 2 * saha_fac * double(u1) / double(u0) * sqrt(double(t)) * double(t)*
    fexparg / double(xne);
  return (T)res;
}

template float eos::witt::nsaha(float, float, float, float, float);
template double eos::witt::nsaha(double, double, double, double, double);


/* --------------------------------------------------------------------------------------- */

template <class T> int eos::witt::getXpart(int iatom, T t, T Pgas, T Pe, T *xpa, T *pf, T *ein, bool divide_by_u)
{

  
  /* --- pre-compute xna and xne --- */

  T TBK = t * phyc::BK, xna = (Pgas - Pe) / TBK, xne = Pe / TBK;
  T n_tot = xna * ABUND[iatom] / tABUND; // partial density of iatom
  
  
  /* --- Get partition function and scaled ionization potential --- */
  
  int nLev = pfn::partition_f<T>(iatom, t, xne, xna, ein, pf, true);


  
  /* --- Solve Saha to get the energy of each stage for this atom --- */

  xpa[0] = 1.0;
  
  for(int ii=1; ii<nLev; ii++) xpa[ii] = nsaha<T>(t, xne, pf[ii-1], pf[ii], ein[ii-1]*phyc::EV);
  for(int ii=nLev-1;ii>0;ii--) xpa[0] = 1.0 + xpa[ii]*xpa[0];
  xpa[0] = 1.0/xpa[0];

  for(int ii=1; ii<nLev;ii++) xpa[ii] *= xpa[ii-1]; 

  if(divide_by_u)
    for(int ii=0; ii<nLev;ii++) xpa[ii] *= n_tot / pf[ii];
  else
    for(int ii=0; ii<nLev;ii++) xpa[ii] *= n_tot;
  
  return nLev;
}

template int eos::witt::getXpart(int, float, float, float, float*, float *, float *,bool);
template int eos::witt::getXpart(int, double, double, double, double*, double *, double *,bool);

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::getN_and_U(int iatom, int istage, T t, T Pgas,
					   T Pe, T &u_ion, bool divide_by_u)
{
   
  T tk = (phyc::BK * t);
  T n_tot = (Pgas -Pe) / tk  * ABUND[iatom] / tABUND ;
  T n_e   = Pe / tk;

  T u0, u1, u2, du0, du1, du2;
  partition_f<T>(iatom, t, u0, u1, u2, du0, du1, du2);

  T Eion1 = phyc::EION1[iatom] * phyc::EV, Eion2 = phyc::EION2[iatom] * phyc::EV;


  
  /* --- Get ion ratios --- */
  
  T n1_n0 = nsaha<T>(t, n_e, u0, u1, Eion1);
  T n2_n1 = (iatom > 0)? nsaha<T>(t, n_e, u1, u2, Eion2) : 0.0;
  

  /* --- get N0 --- */
  
  T n_ion = n_tot / (1.0 + n1_n0 + n2_n1*n1_n0);
  u_ion = u0;

  /* --- Now scale depending on the ion stage --- */
  
  if(istage == 2) n_ion *= n1_n0, u_ion = u1;
  else if(istage >= 3) n_ion *= n2_n1*n1_n0, u_ion = u2;

  if(divide_by_u) n_ion /= u_ion;
  return n_ion;
}

template float eos::witt::getN_and_U(int, int, float, float, float, float&, bool);
template double eos::witt::getN_and_U(int, int, double, double, double, double&, bool);

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::pe_from_pg(T temp, T Pgas, T *fe_out)
{
  
  float dif = 1.1;
  T Pe = 1.2 * init_pe_from_T_Pg<T>(temp, Pgas), oPe = Pe;
  int it = 0;
  
  while((dif > prec) && (it++ < 250)){
    if(dif > 0.5) Pe = (oPe + Pe)*0.5;
    oPe = Pe;
    Pe = pe_pg<T>(temp, Pe, Pgas, fe_out);
    dif = 2.0 * fabs(float(Pe - oPe)/(Pe + oPe));
  }

  return Pe;
}

template float  eos::witt::pe_from_pg (float, float, float*  );
template double eos::witt::pe_from_pg (double, double, double*);

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::rho_from_pg(T temp, T Pgas, T &Pe)
{

  /* ---- 
     we can derive rho relative to H, but we need Pe and N_e / N_H 
     (Mihalas 1970, page 78 - Eq. 3-52): 
     rho = (pe * mu'*mH) / (f_e * BK * T) 
     --- */
  
  T fe_out = 0.0;
  Pe = pe_from_pg<T>(temp, Pgas, &fe_out);
  T rho = Pe * rho_from_H / (fe_out * temp);

  
  return rho;
}

template double  eos::witt::rho_from_pg (double, double, double&);
template float   eos::witt::rho_from_pg (float, float, float&  );


/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::rho_from_pe(T temp, T Pe, T &pg)
{

  /* ---- 
     we can derive rho relative to H, but we need Pe and N_e / N_H 
     (Mihalas 1970, page 78 - Eq. 3-52): 
     rho = (pe * mu'*mH) / (f_e * BK * T) 
     --- */
  
  T fe_out = 0.0;
  pg =  pg_from_pe<T>(temp, Pe, &fe_out);
  T rho = Pe * rho_from_H / (fe_out * temp);
  
  return rho;
}

template double  eos::witt::rho_from_pe (double, double, double&);
template float   eos::witt::rho_from_pe (float, float, float&  );

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::pg_from_rho(T temp, T rho, T &Pe)
{


  /* --- Init Pgas assuming no molecules --- */

  T pgas = rho / avweight * phyc::BK * temp;
  Pe = pe_from_pg<T>(temp, pgas);
  T irho = rho_from_pe<T>(temp, Pe, pgas);
  

  /* --- Iterate until you get a Pe that is consistent with rho --- */

  double dif =  1.0;
  
  int it = 0;
  while ((dif >= 1.e-5) && (it++ < 200)){
    Pe *= (1.0 + rho / irho)*0.5 ;
    irho = rho_from_pe<T>(temp, Pe, pgas);
    dif = fabs((irho - rho) / (rho));
  }
  return pgas;
}

template double  eos::witt::pg_from_rho (double, double, double&);
template float   eos::witt::pg_from_rho (float, float, float&);

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::pe_pg(T temp, T Pe, T Pgas, T *fe_out)
{

  T cmol[2], alfai,
    u0, u1, u2, du0=0, du1=0, du2=0,dcmol[2],
    g1 = 0, g2 = 0, g3 = 0, g4=0, g5=0, theta = 5040./temp,
    a=0, b=0, c=0, d=0, e=0, f1=0, f2=0, f3=0, f4=0, f5=0, fe = 0,
    c3 = 0, c2=0, c1=0, const6, const7;

   
  if(Pe < 0.0) Pe = 1.e-15, g4 = 0.0, g5 = 0.0;
  else{
    molecb<T>(theta, cmol, dcmol);
    acota<T>(cmol[0], -30., -30.); acota<T>(cmol[1], -30., 30.);

    g4 = Pe * pow10<T>( cmol[0]);
    g5 = Pe * pow10<T>( cmol[1]);
  }
  
  /* --- get partition functions and solve saha --- */
  
  partition_f<T>(0, temp, u0, u1, u2, du0, du1, du2);  // First H and H-
  g2 = saha<T>(theta, phyc::EION1[0], u0, u1, Pe); // p(h+)/p(h)
  g3 = saha<T>(theta, 0.754, 1.0, u0, Pe);         // p(h)/p(h-)

  acota<T>(g3, 1.e-30, 1.e30);
  g3 = 1.0 / g3, g1 = 0.0;
  

  /* --- Get number of electrons from each contributing species --- */
  
  for(int i=1; i<ncontr; i++){
    partition_f<T>(i, temp, u0, u1, u2, du0, du1, du2);
    a = saha<T>(theta, phyc::EION1[i], u0, u1, Pe);
    b = saha<T>(theta, phyc::EION2[i], u1, u2, Pe); 
    c = 1.+a*(1.+b);
    alfai = ABUND[i] / ABUND[0];
    g1 += alfai/c*a*(1.+2.*b);
  }
  a=1.+g2+g3;
  e=g2/g5*g4 ;
  b=2.*(1.+e);
  c=g5; 
  d=g2-g3 ;

  acotasig<T>(a, 1.e-15, 1.e15);
  acotasig<T>(d, 1.e-15, 1.e15);

  c1=c*b*b+a*d*b-e*a*a;
  c2=2.*a*e-d*b+a*b*g1;                                                         
  c3=-(e+b*g1) ;                                                        
  f1=0.5*c2/c1 ;                                                              
  f1=-f1+sign<T>(1.,c1)*sqrt(f1*f1-c3/c1) ;
  f5=(1.-a*f1)/b ;
  f4=e*f5 ;
  f3=g3*f1 ;
  f2=g2*f1 ;
  fe=f2-f3+f4+g1 ;

  acota<T>(fe, 1.e-30, 1.e30);
  T phtot=Pe/fe;

  if(f5 < 1.e-4){
    const6 = g5/Pe*f1*f1;
    const7=f2-f3+g1; 
    for(int ii =0; ii<5; ii++)
      f5=phtot*const6, f4=e*f5, fe=const7+f4, phtot=Pe/fe;
  }

                              
  Pe=Pgas/(1.+(f1+f2+f3+f4+f5+Ab_others)/fe);
  if(Pe <= 0) Pe=1.e-15;

  if(fe_out) fe_out[0] = fe;
  
  return Pe;
}

template float  witt::pe_pg (float, float, float, float*);
template double witt::pe_pg (double, double, double, double*);

/* --------------------------------------------------------------------------------------- */

template <class T> T witt::pg_from_pe(T t, T Pe, T *fe_out)
{
  //T p[99], dpp[99], dp[99];
  //T theta = 5040.0 / t;
  //gasb<T>(theta, Pe, p, dp, dpp);
  //return p[83];

  T pp[5], fe = 0.0, Pg = 0.0;
  gasc<T>(t, Pe, Pg, pp, &fe);


  if(fe_out) fe_out[0] = fe;
  
  return Pg;
}

template float witt::pg_from_pe(float, float, float*);
template double witt::pg_from_pe(double, double, double*);

/* --------------------------------------------------------------------------------------- */

template <class T> void witt::molecb(T X, T *Y, T *dy)
{
  Y[0]=-11.206998+X*(2.7942767+X*(7.9196803E-2-X*2.4790744E-2)); // H2
  Y[1]=-12.533505+X*(4.9251644+X*(-5.6191273E-2+X*3.2687661E-3)); // H
  T dx=(-X*X)/5040.;
                                                                        
  dy[0]=dx*(2.7942767+X*(2*7.9196803e-2-X*3*2.4790744e-2)) ;
  dy[1]=dx*(4.9251644+X*(-2*5.6191273e-2-X*3*3.2687661e-3)) ;

  /*!     
    Y(3)=-11.89+X*(3.8084-X*2.4402E-2) ! CH                           
     Y(4)=-11.85+X*(4.1411+X*(-6.7847E-2+4.9178E-3*X)) ! NH            
!     Y(5)=-12.199+X*(4.884+X*(-7.2794E-2+5.1747E-3*X)) ! OH            
!     Y(6)=-11.205+2.0112*X ! SIH                                       
!     Y(7)=-10.514+X*(2.2206-1.3654E-2*X) ! MGH                         
!     Y(8)=-10.581+X*(2.1713+X*(-7.7446E-2+6.5014E-3*X)) ! CAH          
!     Y(9)=-11.711+X*(3.1571-1.5205E-2*X) ! ALH                         
!     Y(10)=-12.24+X*(4.8009-2.6693E-2*X) ! HCL                         
!     Y(11)=-11.849+X*(4.1621-4.1213E-2*X) ! HS                         
!     Y(12)=-12.6+X*(6.3336-1.2019E-2*X) ! C2                           
!     Y(13)=-12.667+X*(7.779-1.1674E-2*X) ! CN                          
!     Y(14)=-13.641+X*(11.591+X*(-8.8848E-2+X*7.3465E-3)) ! CO          
!     Y(15)=-12.07+X*(4.7321-1.4647E-2*X) ! SIC                         
!     Y(16)=-13.122+X*(8.1528-2.5302E-2*X) ! CS                         
!     Y(17)=-13.435+X*(10.541+X*(-2.8061E-1+X*(5.883E-2+X*4.6609E-3)))  
!     Y(18)=-12.606+X*(6.9634+X*(-8.2021E-2+X*6.8654E-3)) ! NO          
!     Y(19)=-12.43+X*(5.409-2.341E-2*X) ! SIN                           
!     Y(20)=-12.172+X*(5.2755-2.3355E-2*X) ! SN                         
!     Y(21)=-13.087+X*(5.3673-1.4064E-2*X) ! O2                         
!     Y(22)=-13.034+8.3616*X ! SIO                                      
!     Y(23)=-11.39+X*(4.267-1.7738E-2*X) ! MGO                          
!     Y(24)=-12.493+X*(5.3382+X*(-5.8504E-2+4.8035E-3*X)) ! ALO         
!     Y(25)=-13.367+X*(8.869+X*(-7.1389E-1+X*(1.5259E-1-1.1909E-2*X)))  
!     Y(26)=-12.876+X*(6.8208+X*(-7.5817E-2+6.3995E-3*X)) ! VO          
!     Y(27)=-13.368+X*(9.0764+X*(-2.8354E-1+2.5784E-2*X)) ! ZRO         
!     Y(28)=-12.645+X*(5.6644-2.2882E-2*X) ! SO                         
!     Y(29)=-11.336+X*(4.4639-1.6552E-2*X) ! NACL                       
!     Y(30)=-12.515+6.7906*X ! SIS                                      
!     Y(31)=-10.638+4.2139*X ! CACL                                     
!     Y(32)=-12.06+X*(5.3309-1.6459E-2*X) ! ALCL                        
!     Y(33)=-12.508+X*(2.727-1.7697E-2*X) ! CL2                         
!     Y(34)=-12.651+X*(4.697-2.5267E-2*X) ! S2                          
!     Y(35)=-24.883+X*(8.2225-2.6757E-1*X) ! CH2                        
!     Y(36)=-24.82+X*(8.4594+X*(-1.2208E-1+8.6774E-3*X)) ! NH2          
!     Y(37)=-25.206+X*(10.311+X*(-9.0573E-2+5.3739E-3*X)) ! H2O         
!     Y(38)=-24.314+X*(8.1063-3.4079E-2*X) ! H2S                        
!     Y(39)=-25.168+13.401*X ! HCN                                      
!     Y(40)=-25.103+X*(12.87-3.8336E-2*X) ! HCO                         
!     Y(41)=-25.078+X*(9.3435+X*(-1.06E-1+8.2469E-3*X)) ! HNO           
!     Y(42)=-25.161+X*(7.8472+X*(-1.0399E-1+7.8209E-3*X)) ! HO2         
!     Y(43)=-27.038+X*(14.376+X*(6.8899E-2+6.0371E-2*X)) ! C3           
!     Y(44)=-25.889+13.317*X ! SIC2                                     
!     Y(45)=-27.261+X*(16.866-1.0144E-2*X) ! CO2                        
!     Y(46)=-26.25+X*(11.83+X*(-4.2021E-2+3.4397E-3*X)) ! N2O           
!     Y(47)=-26.098+X*(10.26+X*(-1.0101E-1+8.4813E-3*X)) ! NO2          
!     Y(48)=-26.115+X*(6.5385-1.9332E-2*X) ! O3                         
!     Y(49)=-27.496+X*(13.549+X*(-2.205E-2+2.1407E-2*X)) ! TIO2         
!     Y(50)=-27.494+X*(15.99+X*(-2.3605E-1+2.1644E-2*X)) ! ZRO2         
!     Y(51)=-25.244+X*(11.065-3.7924E-2*X) ! AL2O                       
!     Y(52)=-23.748+X*(9.5556-2.4867E-2*X) ! ALCL2                      
!     Y(53)=-37.194+X*(13.371-3.4229E-2*X) ! CH3                        
!     Y(54)=-37.544+X*(12.895-4.9012E-2*X) ! NH3                        
!     Y(55)=-37.931+17.216*X ! C2H2                                     
!     Y(56)=-38.274+X*(16.264-3.2379E-2*X) ! HCOH                       
!     Y(57)=-38.841+X*(18.907-3.5705E-2*X) ! HCNO                       
!     Y(58)=-50.807+X*(17.956-3.6861E-2*X) ! CH4                        
!     Y(59)=-11.4575+X*(3.1080922+X*(-3.3159806E-1+4.314945E-2*X)) ! NAH
!     Y(60)=-10.964723+X*(2.270225+X*(-7.66888E-2+6.519213E-3*X)) ! KH  
!     Y(61)=-10.807839+X*(2.744854+X*(5.758024E-2+3.315373E-3*X)) ! BEH 
!     Y(62)=-10.491008+X*(2.051217+X*(-7.643729E-2+6.425358E-3*X)) ! SRH
!     Y(63)=-11.01929+X*(3.13829+X*(1.214975-1.77077E-1*X)) ! SRO       
!     Y(64)=-10.446909+X*(2.024548+X*(-7.680739E-2+6.471443E-3*X)) ! BAH
!     Y(65)=-10.921254+X*(3.847116+X*(1.189653-1.662815E-1*X)) ! BAO    
!     Y(66)=-13.561415+X*(7.528355+X*(-5.031809E-1+6.787392E-2*X)) ! SCO
!     Y(67)=-14.107593+X*(12.226315+X*(-1.019148+1.02046E-1*X)) ! YO    
!     Y(68)=-14.231303+X*(11.28907+X*(-1.108545+1.274056E-1*X)) ! LAO   
!     Y(69)=-12.03193+X*(3.012432+X*(1.798894E-1-1.79236E-2*X)) ! SI2   
!     Y(70)=-11.344906+X*(2.836732+X*(-1.134115E-1+1.99661E-2*X)) ! LIH 
!     Y(71)=-25.913244+X*(11.856324+X*(1.05407-1.541093E-1*X)) ! VO2    
!     Y(72)=-26.934577+X*(13.421189+X*(2.671335E-1-3.475224E-2*X)) ! SIO
!     Y(73)=-23.225499+X*(4.820973+X*(6.722119E-1-5.903448E-2*X)) ! SIH2
!     Y(74)=-25.079417+X*(7.196588+X*(1.196713E-1+1.0484E-2*X)) ! SI3   
!     Y(75)=-26.331315+X*(12.500392+X*(6.531014E-1-1.162598E-1*X)) ! C2H
!     Y(76)=-11.673776+X*(3.245147+X*(1.334288E-1-1.524113E-2*X)) ! BH  
!     Y(77)=-12.972588+X*(7.80983-X*(6.263376E-2-4.763338E-3*X)) ! BO   
!     Y(78)=-12.654+X*(6.2558-3.0868E-2*X) ! HF                         
!     Y(79)=-11.639+X*(3.9742-X*(1.7229E-1-1.687E-2*X)) ! LIO           
!     Y(80)=-11.92+X*(6.1426-1.8981E-2*X) ! LIF                         
!     Y(81)=-11.576+X*(5.1447+X*(1.4625E-2-8.9675E-3*X)) ! LICL         
!     Y(82)=-12.579+X*(4.8824+X*(-1.3848E-1+1.25E-2*X)) ! FEO           
!     Y(83)=-11.666+X*(5.1748-1.9602E-2*X) ! NAF                        
!     Y(84)=-11.292+X*(4.851-1.8104E-2*X) ! MGF                         
!     Y(85)=-12.453+X*(7.1023-1.6086E-2*X) ! ALF                        
!     Y(86)=-11.8+X*(5.7274-X*(1.9201E-1-2.1306E-2*X)) ! KF             
!     Y(87)=-10.798+X*(2.6664+X*(1.2037E-1-1.778E-2*X)) ! MGCL          
!     Y(88)=-11.453+X*(4.9299+X*(-1.4643E-1+1.4143E-2*X)) ! KCL         
!     Y(89)=-11.484+X*(3.2399-2.2356E-2*X) ! FECL                       
!     Y(90)=-24.304+X*(9.5257-4.2841E-2*X) ! LIOH                       
!     Y(91)=0.                      
  */
  
}
template void witt::molecb(float, float*, float*);
template void witt::molecb(double, double*, double*);

/* --------------------------------------------------------------------------------------- */

template <class T> void eos::witt::contOpacity(T temp, T iPgas, T iPe, int nw, double *w, double *opac)
{

  /* --- Some definitions --- */
  
  double iT = (double)temp,  TK  = phyc::BK*iT, TKEV = TK / phyc::EV;
  double HTK  = phyc::HH/TK, TLOG = log(iT), Pgas = (double)iPgas, Pe = (double)iPe;
  double xna = (Pgas - Pe) / TK, xne = Pe / TK;
  
  
  /* --- Get the partial densities of background absorvers --- */
  
  double n[17] = {}, *scat = new double [nw]();
  getBackgroundPartials<double>(iT, Pgas, Pe, n);

  
  /* --- Get background opacities, the cop routines work with partial densities / partition function! --- */
  
  cop(iT,  TKEV, TK, HTK, TLOG, xna, xne, w, opac, scat, n[0], n[1], n[2], n[3],
      n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16],
      nw, 0, 0); 
  
  delete [] scat;
}

template void witt::contOpacity(float, float, float, int, double*, double*);
template void witt::contOpacity(double, double, double, int, double*, double*);


/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::saha(T theta, T chi, T u1, T u2, T Pe)
{
  double res = u2 * exp(2.3025850929940459 * (9.0804625434325867-(double)(theta*chi))) / ((double)u1*Pe*pow((double)theta, 2.5));
  return (T)res;
}

template float  witt::saha (float, float, float, float, float);
template double witt::saha (double, double, double, double, double);

/* --------------------------------------------------------------------------------------- */

template <class T> T eos::witt::dsaha(T theta, T chi, T du1, T du2)
{
  double res = du2-du1+(theta/5040.)*(2.5+(double)chi*theta*log(10.0));
  return (T)res;
}

template float  witt::dsaha (float, float, float, float);
template double witt::dsaha (double, double, double, double);

/* --------------------------------------------------------------------------------------- */

 template <class T> void witt::gasc(T t, T pe, T &pg, T *pp, T *fe_out)
{
  T cmol[2], chi1, chi2, u0, u1,
    u2, du0, du1, du2, dcmol[2], theta=5040./t, alfai;// alfai[ncontr]

  molecb<T>(theta, cmol, dcmol);
  T g4=pe*pow10<T>( cmol[0]);
  T g5=pe*pow10<T>( cmol[1]); 

  chi1 = phyc::EION1[0], chi2 = phyc::EION2[0];
  partition_f<T>(0, t, u0, u1, u2, du0, du1, du2);
  
  T g2 = saha<T>(theta,chi1,u0,u1,pe);
  T g3 = 1.0/saha<T>(theta,0.754,1.0,u0,pe); // Ph- / Ph
  T g1=0.0, a, b, c, d, e, c1, c2, c3, f1, f2,f3,f4,f5, fe, phtot;
  
  for(int ii=1; ii<ncontr; ii++){
    alfai = ABUND[ii]/ABUND[0], chi1 = phyc::EION1[ii], chi2 = phyc::EION2[ii];
    partition_f<T>(ii, t, u0, u1, u2, du0, du1, du2);
    
    a=saha<T>(theta,chi1,u0,u1,pe); // n1 / n0
    b=saha<T>(theta,chi2,u1,u2,pe); // n2 / n1
    c=1.+a*(1.+b); // fraction of n0/ntot
    //pp[ii]=alfai[ii]/c; // ABUND /ntot (partial pressure of neutral species ii)
    
    /* --- 
       1*n1 + 2*n2 + ... j*n_j so we count how many electrons
       come from each ionized species 
       ---*/
    g1 += alfai/c *a*(1.+2.*b);
  }
                                                                          
  a=1.+g2+g3;
  e=g2/g5*g4;
  b=2.0*(1.0 + e);
  c=g5;
  d=g2-g3;
  c1=c*b*b+a*d*b-e*a*a;
  c2=2.*a*e-d*b+a*b*g1;
  c3=-(e+b*g1);
  f1=0.5*c2/c1;
  f1=-f1+sign<T>(1.0,c1)*sqrt(f1*f1-c3/c1);
  f5=(1.0-a*f1)/b;
  f4=e*f5;
  f3=g3*f1;
  f2=g2*f1;
  fe=f2-f3+f4+g1; // N_e / N_H
  phtot=pe/fe; 

                                                                 
  if(f5 <= 1.e-4){
    T diff = 1.0, const6=g5/pe*f1*f1, const7=f2-f3+g1; 

    int it = 0;
    while((diff > 1.e-5) && (it++ < 5)){
      T of5 = f5;
      f5=phtot*const6;
      f4=e*f5;
      fe=const7+f4;
      phtot=pe/fe;
      diff = 0.5 * fabs(f5-of5) / (f5 + of5);
    }
  }
  
  //pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe);
  pg=pe*(1.0+(f1+f2+f3+f4+f5+Ab_others)/fe);                                       

  if(fe_out) fe_out[0] = fe;
  
  pp[0]=f1;//! p(h)/p(h')                      
  pp[1]=f2;//! p(h+)/p(h')                     
  pp[2]=f5;//! p(h2)/p(h') 
  pp[3]=f3;//! p(h-)/p(h')
  pp[4] = phtot; // p(h')
}

template void witt::gasc(float, float, float&, float*, float*);
template void witt::gasc(double, double, double&, double*, double*);

/* --------------------------------------------------------------------------------------- */

// template <class T> void witt::gasb(T theta, T pe, T *p, T *dp, T *ddp)
// {
//   T cmol[91], alfai[ncontr], chi1[ncontr], chi2[ncontr], u0[ncontr], u1[ncontr],
//     u2[ncontr], du0[ncontr], du1[ncontr], du2[ncontr], dcmol[91], g4,g5,
//     dg4,dg5,ddg4,ddg5, t = 5040.0 / theta;
  
//   if(pe <= 0.0) pe = 1.e-10, g4 = 0, g5=0,dg4=0,dg5=0,ddg4=0,ddg5=0;
//   else{
//     molecb<T>(theta,cmol,dcmol);
//     acota<T>(cmol[0], -30.0, 30.0); acota<T>(cmol[1], -30.0, 30.0);
    
//     g4=pe*pow10<T>(cmol[0]); 
//     dg4=dcmol[0]*log(10.0);
//     ddg4=1./pe;
//     g5=pe*pow10<T>( cmol[1]); 
//     dg5=dcmol[1]*log(10.0); 
//     ddg5=1./pe;
//   }

//   for(int ii=0;ii<ncontr; ii++){
//     alfai[ii] = ABUND[ii]/ABUND[0];
//     chi1[ii] = phyc::EION1[ii];
//     chi2[ii] = phyc::EION2[ii];
//     partition_f<T>(ii, t, u0[ii], u1[ii], u2[ii], du0[ii], du1[ii], du2[ii]);
//   }

//   T g2 = saha<T>(theta,chi1[0],u0[0],u1[0],pe);
//   T dg2 = dsaha<T>(theta,chi1[0],du0[0],du1[0]);
//   T ddg2=-1./pe;
//   //!derivada log(g2) con pe                     
//   p[91]=g2; 
//   dp[91]=dg2;
//   ddp[91]=ddg2;
//   //! p(h)/p(h-)                   
//   T g3=saha<T>(theta,0.754,1.,u0[0],pe);
//   acota<T>(g3,1.e-30,1.e30);
  
//   g3=1.0/g3;
//   T dg3=-1.0 * dsaha<T>(theta,0.754,0.,du0[0]);
//   T ddg3=1./pe;
//   p[92]=g3; 
//   dp[92]=dg3;
//   ddp[92]=ddg3;
//   T g1=0.0;
//   // !las dl son derivadas no de log(g) sino de g        
//   T dlg1=0.;
//   T ddlg1=0.,a,da,dda,b,dlb,ddlb,c,ss1,ss,dss,ddss;

  
//   for(int ii=1;ii<ncontr;ii++){
//     a=saha<T>(theta,chi1[ii],u0[ii],u1[ii],pe);
//     da=dsaha<T>(theta,chi1[ii],du0[ii],du1[ii]);
//     dda=-1./pe;
//     b=saha<T>(theta,chi2[ii],u1[ii],u2[ii],pe);
//     dlb=b*dsaha<T>(theta,chi2[ii],du1[ii],du2[ii]);
//     ddlb=-b/pe;
//     c=1.+a*(1.+b); 

//     acotasig<T>(c,1.e-20,1.e20);
//     //  ! p/ph' for neutral he,li, ...                  
//     p[ii]=alfai[ii]/c;
//     dp[ii]=-(a*da*(1.+b)+a*dlb)/c;
//     ddp[ii]=-(a*dda*(1.+b)+a*ddlb)/c;
//     ss1=(1.+2.*b);
//     acotasig<T>(ss1,1.e-20,1.e20);
//     ss=p[ii]*a*ss1;
//     dss=dp[ii]+da+2.*dlb/ss1;
//     ddss=ddp[ii]+dda+2.*ddlb/ss1;
                                                                        
// //!      g1=g1+p(i)*a*ss1                                                 
//     g1=g1+ss; 
//     //                               !ojo estas no son derivadas del log  
//     dlg1=dlg1+ss*dss;
//     ddlg1=ddlg1+ss*ddss;
//   }
  
//   a=1.+g2+g3; 
//   T dla=g2*dg2+g3*dg3, ddla=g2*ddg2+g3*ddg3;

//   if(g5 < 1.e-35){
//     fprintf(stderr,"error: witt::gasb: error, Pel is too small [Pe=%e]\n", g5);
//     exit(0);
//   }
  
//   acotasig<T>(g5,1.e-20,1.e20);
//   T dlc, ddlc, d, dld, ddld, e, de, ddle, c1, dlc1, c2, dlc2, ddlc2, c3, dlc3, ddlc3, f1,
//     dc1, dc2, dle, dde, ddlc1, dlf1;
  
//   b=2.*(1.+g2/g5*g4);
//   dlb=(b-2.)*(dg2-dg5+dg4);
//   ddlb=(b-2.)*(ddg2-ddg5+ddg4);
//   c=g5;
//   dlc=dg5*g5;
//   ddlc=ddg5*g5;
//   d=g2-g3;
//   dld=g2*dg2-g3*dg3;
//   ddld=g2*ddg2-g3*ddg3;
//   e=g2/g5*g4;
//   de=dg2-dg5+dg4;
//   dde=ddg2-ddg5+ddg4;
//   dle=e*de;
//   ddle=e*dde;
                                                                        
//   acotasig<T>(a,1.e-15,1.e15);
//   acotasig<T>(b,1.e-15,1.e15);
//   acotasig<T>(c,1.e-15,1.e15);
//   acotasig<T>(d,1.e-15,1.e15); 
//   acotasig<T>(e,1.e-15,1.e15); 
                                                                        
                                                                        
//   c1=c*b*b+a*d*b-e*a*a;
//   dlc1=dlc*b*b+(c*2.*b+a*d)*dlb+dla*(d*b-2.*e*a)+dld*a*b-dle*a*a;
//   ddlc1=ddlc*b*b+(c*2.*b+a*d)*ddlb+ddla*(d*b-2.*e*a)+ddld*a*b-ddle*a*a;
//   c2=2.*a*e-d*b+a*b*g1;
//   dlc2=dla*(2.*e+b*g1)+dlb*(a*g1-d)-dld*b+dle*2.*a+a*b*dlg1;
//   ddlc2=ddla*(2.*e+b*g1)+ddlb*(a*g1-d)-ddld*b+ddle*2.*a+a*b*ddlg1;
//   c3=-(e+b*g1);
//   dlc3=-dle-dlb*g1-b*dlg1;
//   ddlc3=-ddle-ddlb*g1-b*ddlg1;
  
//   acotasig<T>(c1,1.e-15,1.e15);
                                                                        
//   f1=0.5*c2/c1;
//   //!      dlf1=.5*(dlc2/c1-(dlc1*c2)/(c1*c1))                              
//   dc1=dlc1/c1;
//   dc2=dlc2/c2;
//   dlf1=f1*(dc2-dc1);
//   T ddlf1=.5*(ddlc2/c1-(ddlc1*c2)/(c1*c1));                        
//   T ddc1=ddlc1/c1;
//   T ddc2=ddlc2/c2;
//   ddlf1=f1*(ddc2-ddc1);
                                                                     
//   dlf1=-dlf1+sign<T>(1.,c1)*(2.*f1*dlf1-dlc3/c1+dlc1*c3/(c1*c1))/(2.*sqrt(f1*f1-c3/c1));
//   ddlf1=-ddlf1+sign<T>(1.,c1)*(2.*f1*ddlf1-ddlc3/c1+ddlc1*c3/(c1*c1)) /(2.*sqrt(f1*f1-c3/c1));
//   f1=-f1+sign<T>(1.,c1)*sqrt(f1*f1-c3/c1);
                                                                        
//   T f5=(1.-a*f1)/b, dlf5, ddlf5, df5, ddf5, f4, df4, ddf4,  ddlf4,
//     ddlf3, f2, dlf2, ddlf2, divf1, ddivf1, dlf200, ddlf200, dlf4, f3, dlf3;
//   if(fabs(f5) < 1.e-30){
//     dlf5=0.;
//     ddlf5=0.; 
//     df5=0.; 
//     ddf5=0.; 
//   }else{
//     dlf5=(-dla*f1-a*dlf1)/b-((1.-a*f1)*dlb)/(b*b);
//     ddlf5=(-ddla*f1-a*ddlf1)/b-((1.-a*f1)*ddlb)/(b*b);
//     df5=dlf5/f5;
//     ddf5=ddlf5/f5;
//   }
//   f4=e*f5;
  
//   if(fabs(f4) < 1.e-30){
//     dlf4=0.;
//     ddlf4=0.;
//     df4=0.;
//     ddf4=0.;
//   }else{
//     dlf4=f5*dle+e*dlf5;
//     ddlf4=f5*ddle+e*ddlf5;
//     df4=dlf4/f4 ;
//     ddf4=ddlf4/f4 ;
//   }
                                                                        
  
//   f3=g3*f1;
//   dlf3=f3*dg3+g3*dlf1;
//   ddlf3=f3*ddg3+g3*ddlf1;
//   f2=g2*f1;
//   dlf2=f2*dg2+g2*dlf1;
//   ddlf2=f2*ddg2+g2*ddlf1;
                                                                        
//   if(abs(f1) < 1.e-30){ 
//     divf1=0.; 
//     ddivf1=0. ;
//   }else{
//     divf1=dlf1/f1;
//     ddivf1=ddlf1/f1;
//   }
                                                                        
//   dlf200=dg2+divf1;
//   ddlf200=ddg2+ddivf1;
                                                                        
//   T fe=f2-f3+f4+g1, dlfe, ddlfe, dfe, ddfe;
//   if(abs(fe) < 1.e-30){
//     dlfe=0.;
//     ddlfe=0.;
//     dfe=0.;
//     ddfe=0.;
//   }else{
//     dlfe=dlf2-dlf3+dlf4+dlg1;
//     ddlfe=ddlf2-ddlf3+ddlf4+ddlg1;
//     dfe=dlfe/fe;
//     ddfe=ddlfe/fe;
//   }
                                                                        
                                                                        
//   acotasig<T>(fe,1.e-15,1.e15);
//   T phtot=pe/fe;
//   //T dlphtot=-phtot*dfe;
//   //T ddlphtot=(1.-ddlfe*phtot)/fe;
//   T dphtot=-dfe;
//   T ddphtot=1./pe-ddfe;
//   T const6, dconst6, ddconst6, const7, dlconst7,  ddlconst7, //dconst7,ddconst7 
//     pg, dlpg, ddlpg, divf3, ddivf3;
                                                                        
//   if(f5 <= 1.e-4){
//     const6=g5/pe*f1*f1;
//     //!dlf1/f1                    
//     dconst6=dg5+2.*divf1;
//     //!ddlf1/f1                
//     ddconst6=ddg5-1./pe+2.*ddivf1;
//     const7=f2-f3+g1;
//     acotasig<T>(const7,1.e-15,1.e15);
//     dlconst7=dlf2-dlf3+dlg1;
//     ddlconst7=ddlf2-ddlf3+ddlg1;
//     //dconst7=dlconst7/const7;
//     //ddconst7=ddlconst7/const7;
    
//     for(int ii=0; ii<5; ii++){ 
//       f5=phtot*const6;
//       df5=dphtot+dconst6;
//       ddf5=ddphtot+ddconst6;
//       f4=e*f5;
//       df4=de+df5;
//       ddf4=dde+ddf5;
//       fe=const7+f4;
//       acotasig<T>(fe,1.e-15,1.e15);
//       dfe=(dlconst7+df4*f4)/fe;
//       ddfe=(ddlconst7+ddf4*f4)/fe;
//       phtot=pe/fe;
//       dphtot=-dfe;
//       ddphtot=1./pe-ddfe;
//     }
                                                                        
//     dlf5=df5*f5;
//     ddlf5=ddf5*f5;
//     dlf4=df4*f4;
//     ddlf4=ddf4*f4;
//     dlfe=dfe*fe;
//     ddlfe=ddfe*fe;
//     //dlphtot=dphtot*phtot;
//     //ddlphtot=ddphtot*phtot;
//   }
  
//   pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe);
//   dlpg=pe*(dlf1+dlf2+dlf3+dlf4+dlf5)/fe-(pg-pe)*(dlfe/fe);
//   ddlpg=pe*(ddlf1+ddlf2+ddlf3+ddlf4+ddlf5)/fe-(pg-pe)*(ddlfe/fe);
//   ddlpg=ddlpg+pg/pe;
//   //! p(h)/p(h')                                            
//   p[0]=f1; 
//   //!dlf1/f1                                   
//   dp[0]=divf1;
//   //!ddlf1/f1                                
//   ddp[0]=ddivf1;
  
//   acotasig<T>(pg,1.e-20,1.e20);
//   //! gas pressure                                         
//   p[83]=pg;
//   dp[83]=dlpg/pg;
//   ddp[83]=ddlpg/pg;
//   //! p(h')                                             
//   p[84]=phtot;
//   //!dlphtot/phtot                     
//   dp[84]=dphtot;
//   //!ddlphtot/phtot                  
//   ddp[84]=ddphtot;
//   //! p(h+)/p(h')                                          
//   p[85]=f2;
//   dp[85]=dlf200;
//   ddp[85]=ddlf200;
//   //! p(h-)/p(h')                                          
//   p[86]=f3;
//   if(fabs(f3) < 1.e-30){
//     divf3=0.;
//     ddivf3=0.;
//   }else{
//     divf3=dlf3/f3;
//     ddivf3=ddlf3/f3;
//   }
                                                                        
//   dp[86]=divf3;
//   ddp[86]=ddivf3;
//   //! p(h2+)/p(h')                                         
//   p[87]=f4;
//   //!dlf4/f4                                    
//   dp[87]=df4;
//   //!ddlf4/f4                                 
//   ddp[87]=ddf4;
//   //! p(h2)/p(h')                                          
//   p[88]=f5;
//   //!dlf5/f5                                    
//   dp[88]=df5;
//   //!ddlf5/f5                                 
//   ddp[88]=ddf5;
//   //! pe/p(h')                                             
//   p[89]=fe;
//   //!dlfe/fe                                    
//   dp[89]=dfe;
//   //!ddlfe/fe                                 
//   ddp[89]=ddfe;
  
//   //! n(e)=pe/kt                             
//   p[90]=pe/(1.38054e-16*t);
//   dp[90]=-1./t;
//   ddp[90]=1./pe;
  
// }
// template void witt::gasb(float, float, float*, float*, float*);
// template void witt::gasb(double, double, double*, double*, double*);


/* --------------------------------------------------------------------------------------- */

  template <class T> void eos::witt::partition_f(int nel_in, T t_in, T &u1_in, T &u2_in, T &u3_in, T &du1_in, T &du2_in, T &du3_in)
{
  
/* ----
   Converted from Fortran using f2c (quick-n-dirty!).
   Removed all dependencies on f2c.h and made use of std::max.


   J. de la Cruz Rodriguez (ISP-SU 2017)

   ---- */
  
  /* System generated locals */
  double r__1, r__2, x, y, dx, t = (double)t_in, u1=0.0, u2=0.0 , u3=0.0, du1=0.0, du2=0.0, du3=0.0;
  int nel = nel_in + 1; // To get same indexing as fortran;
  
  
  /*     This subroutine computes partition functions. NEL is element number and */
  /*     T is temperature. */
  
  x = log(5040.0 / t);
  dx = -1.0 / t;
  y = t * .001;
  
  if (nel == 1) {
    goto L1;
  } else if (nel == 2) {
    goto L2;
  } else if (nel == 3) {
    goto L3;
  } else if (nel == 4) {
    goto L4;
  } else if (nel == 5) {
    goto L5;
  } else if (nel == 6) {
    goto L6;
  } else if (nel == 7) {
    goto L7;
  } else if (nel == 8) {
    goto L8;
  } else if (nel == 9) {
    goto L9;
  } else if (nel == 10) {
    goto L10;
  } else if (nel == 11) {
    goto L11;
  } else if (nel == 12) {
    goto L12;
  } else if (nel == 13) {
    goto L13;
  } else if (nel == 14) {
    goto L14;
  } else if (nel == 15) {
    goto L15;
  } else if (nel == 16) {
    goto L16;
  } else if (nel == 17) {
    goto L17;
  } else if (nel == 18) {
    goto L18;
  } else if (nel == 19) {
    goto L19;
  } else if (nel == 20) {
    goto L20;
  } else if (nel == 21) {
    goto L21;
  } else if (nel == 22) {
    goto L22;
  } else if (nel == 23) {
    goto L23;
  } else if (nel == 24) {
    goto L24;
  } else if (nel == 25) {
    goto L25;
  } else if (nel == 26) {
    goto L26;
  } else if (nel == 27) {
    goto L27;
  } else if (nel == 28) {
    goto L28;
  } else if (nel == 29) {
    goto L29;
  } else if (nel == 30) {
    goto L30;
  } else if (nel == 31) {
    goto L31;
  } else if (nel == 32) {
    goto L32;
  } else if (nel == 33) {
    goto L33;
  } else if (nel == 34) {
    goto L34;
  } else if (nel == 35) {
    goto L35;
  } else if (nel == 36) {
    goto L36;
  } else if (nel == 37) {
    goto L37;
  } else if (nel == 38) {
    goto L38;
  } else if (nel == 39) {
    goto L39;
  } else if (nel == 40) {
    goto L40;
  } else if (nel == 41) {
    goto L41;
  } else if (nel == 42) {
    goto L42;
  } else if (nel == 43) {
    goto L43;
  } else if (nel == 44) {
    goto L44;
  } else if (nel == 45) {
    goto L45;
  } else if (nel == 46) {
    goto L46;
  } else if (nel == 47) {
    goto L47;
  } else if (nel == 48) {
    goto L48;
  } else if (nel == 49) {
    goto L49;
  } else if (nel == 50) {
    goto L50;
  } else if (nel == 51) {
    goto L51;
  } else if (nel == 52) {
    goto L52;
  } else if (nel == 53) {
    goto L53;
  } else if (nel == 54) {
    goto L54;
  } else if (nel == 55) {
    goto L55;
  } else if (nel == 56) {
    goto L56;
  } else if (nel == 57) {
    goto L57;
  } else if (nel == 58) {
    goto L58;
  } else if (nel == 59) {
    goto L59;
  } else if (nel == 60) {
    goto L60;
  } else if (nel == 61) {
    goto L61;
  } else if (nel == 62) {
    goto L62;
  } else if (nel == 63) {
    goto L63;
  } else if (nel == 64) {
    goto L64;
  } else if (nel == 65) {
    goto L65;
  } else if (nel == 66) {
    goto L66;
  } else if (nel == 67) {
    goto L67;
  } else if (nel == 68) {
    goto L68;
  } else if (nel == 69) {
    goto L69;
  } else if (nel == 70) {
    goto L70;
  } else if (nel == 71) {
    goto L71;
  } else if (nel == 72) {
    goto L72;
  } else if (nel == 73) {
    goto L73;
  } else if (nel == 74) {
    goto L74;
  } else if (nel == 75) {
    goto L75;
  } else if (nel == 76) {
    goto L76;
  } else if (nel == 77) {
    goto L77;
  } else if (nel == 78) {
    goto L78;
  } else if (nel == 79) {
    goto L79;
  } else if (nel == 80) {
    goto L80;
  } else if (nel == 81) {
    goto L81;
  } else if (nel == 82) {
    goto L82;
  } else if (nel == 83) {
    goto L83;
  } else if (nel == 84) {
    goto L84;
  } else if (nel == 85) {
    goto L85;
  } else if (nel == 86) {
    goto L86;
  } else if (nel == 87) {
    goto L87;
  } else if (nel == 88) {
    goto L88;
  } else if (nel == 89) {
    goto L89;
  } else if (nel == 90) {
    goto L90;
  } else if (nel == 91) {
    goto L91;
  } else if (nel == 92) {
    goto L92;
  }
 L1:
  u1 = 2;
  /* H */
  if (t > 1.3e4) {
    u1 = t * 3.8e-5 + 1.51;
    du1 = 3.8e-5 / u1;
  }
  if (t > 16200) {
    u1 = t * (t * 3.52e-8 - .0011428) + 11.41;
    du1 = (t * 2 * 3.52e-8 - .0011428) / u1;
  }
  u2 = 1;
  u3 = 0;
  goto LEND;
 L2:
  u1 = 1;
  /* HE */
  if (t > 3e4) {
    u1 = t * (t * 1.6095e-8 - 9.4103e-4) + 14.8;
    du1 = (t * 2 * 1.6095e-8 - 9.4103e-4) / u1;
  }
  u2 = 2;
  u3 = 1;
  goto LEND;
 L3:
  u1 = 2.081 - y * (.068926 - y * .014081);
  /* LI */
  du1 = (y * 2 * .014081 - .068926) * .001 / u1;
  if (t > 6e3) {
    u1 = t * (t * 8.5586e-8 - 7.3292e-4) + 3.4864;
    du1 = (t * 2 * 8.5586e-8 - 7.3292e-4) / u1;
  }
  u2 = 1;
  u3 = 2;
  goto LEND;
 L4:
  /* Computing MAX */
  u1 = std::max(1.0,t * 7.032e-5 + .631);
  /* BE */
  if (u1 != 1) {
    du1 = 7.032e-5 / u1;
  }
  u2 = 2;
  u3 = 1;
  goto LEND;
 L5:
  u1 = y * .010438 + 5.9351;
  /* B */
  du1 = 1.0437999999999999e-5 / u1;
  u2 = 1;
  u3 = 2;
  goto LEND;
 L6:
  u1 = y * (y * (.017629 - y * 3.9091e-4) + .020485) + 8.6985;
  /* C */
  du1 = (y * 2 * .017629 + .020485 - y * 3 * y * 3.9091e-4) / u1 * 
    .001;
  if (t > 1.2e4) {
    u1 = t * (t * 9.0844e-8 - .0013907) + 13.97;
    du1 = (t * 2 * 9.0844e-8 - .0013907) / u1;
  }
  u2 = t * 1.6833e-5 + 5.838;
  du2 = 1.6833e-5 / u2;
  if (t > 2.4e4) {
    u2 = t * (t * 2.0861e-8 - 6.9347e-4) + 10.989;
    du2 = (t * 2 * 2.0861e-8 - 6.9347e-4) / u2;
  }
  u3 = 1;
  if (t > 19500) {
    u3 = t * 8e-5 - .555;
    du3 = 8e-5 / u3;
  }
  goto LEND;
 L7:
  u1 = y * (.017491 - y * (.010148 - y * .0017138)) + 3.9914;
  /* N */
  du1 = (.017491 - y * 2 * .010148 + y * 3 * y * .0017138) / u1 * 
    .001;
  if (t > 8800) {
    u1 = t * 2.54e-4 + 2.171;
    du1 = 2.54e-4 / u1;
  }
  if (t > 1.8e4) {
    u1 = t * (t * 8.633e-8 - .0017139) + 11.396;
    du1 = (t * 2 * 8.633e-8 - .0017139) / u1;
  }
  u2 = t * 1.42e-4 + 8.06;
  du2 = 1.42e-4 / u2;
  if (t > 3.3e4) {
    u2 = t * (t * 4.4612e-8 - .0018931) + 26.793;
    du2 = (t * 2 * 4.4612e-8 - .0018931) / u2;
  }
  u3 = t * (t * 1.8228e-9 - 2.6651e-5) + 5.9835;
  du3 = (t * 2 * 1.8228e-9 - 2.6651e-5) / u3;
  if (t < 7310.5) {
    u3 = 5.89;
    du3 = 0;
  }
  goto LEND;
 L8:
  u1 = t * 1.1e-4 + 8.29;
  /* O */
  du1 = 1.1e-4 / u1;
  if (t > 1.9e4) {
    u1 = t * (t * 1.657e-7 - .006019) + 66.81;
    du1 = (t * 2 * 1.657e-7 - .006019) / u1;
  }
  /* Computing MAX */
  u2 = std::max(4.0,t * 8e-5 + 3.51);
  if (u2 != 4) {
    du2 = 8e-5 / u2;
  }
  if (t > 36400) {
    u2 = t * (t * 6.885e-8 - .004216) + 68.7;
    du2 = (t * 2 * 6.885e-8 - .004216) / u2;
  }
  u3 = t * 1.1348e-4 + 7.865;
  du3 = 1.1348e-4 / u3;
  goto LEND;
 L9:
  u1 = y * (y * (y * (.026771 - y * .0013035) - .20884) + .77683) + 
    4.5832;
  /* F */
  du1 = .77683 - y * 2 * .20884 + y * 3 * y * .026771 - y * 4 * y * y * 
    .0013035;
  du1 = du1 / u1 * .001;
  if (t > 8750) {
    u1 = 5.9;
    du1 = 0;
  }
  if (t > 2e4) {
    u1 = t * (t * 2.312e-8 - 9.229e-4) + 15.16;
    du1 = (t * 2 * 2.312e-8 - 9.229e-4) / u1;
  }
  u2 = t * 8.9e-5 + 8.15;
  du2 = 8.9e-5 / u2;
  u3 = t * 1.38e-4 + 2.315;
  du3 = 1.38e-4 / u3;
  goto LEND;
 L10:
  u1 = 1;
  /* NE */
  if (t > 26900) {
    u1 = t * (t * 4.359e-8 - .002113) + 26.3;
    du1 = (t * 2 * 4.359e-8 - .002113) / u1;
  }
  u2 = t * 4e-5 + 5.4;
  du2 = 4e-5 / u2;
  u3 = t * 7.956e-5 + 7.973;
  du3 = 7.956e-5 / u3;
  goto LEND;
 L11:
  /* Computing MAX */
  u1 = std::max(2.0, t * 9.3e-5 + 1.72);
  /* NA */
  if (u1 != 2.0) {
    du1 = 9.3e-5 / u1;
  }
  if (t > 5400) {
    u1 = t * 5.66e-4 - .83;
    du1 = 5.66e-4 / u1;
  }
  if (t > 8500) {
    u1 = t * (t * 1.3861e-7 - .0012415) + 4.5568;
    du1 = (t * 2 * 1.3861e-7 - .0012415) / u1;
  }
  u2 = 1;
  u3 = t * 5.69e-6 + 5.69;
  du3 = 5.69e-6 / u3;
  goto LEND;
 L12:
  u1 = exp(-4.027262 - x * (x * (x * (x * .784131 + 2.393895) + 
				  2.889176) + 6.173172)) + 1;
  du1 = -dx * (x * (x * (x * 4 * .784131 + 7.1816849999999999) + 
		    5.7783519999999999) + 6.173172);
  du1 = du1 * (u1 - 1) / u1;
  if (t > 8e3) {
    u1 = t * (t * 7.4531e-8 - 7.8909e-4) + 2.757;
    du1 = (t * 2 * 7.4531e-8 - 7.8909e-4) / u1;
  }
  u2 = exp(-7.721172 - x * (x * (x * .212417 + 1.966097) + 7.600678)) 
    + 2;
  du2 = -dx * (x * (x * 3 * .212417 + 3.932194) + 7.0600678);
  du2 = du2 * (u2 - 2) / u2;
  if (t > 2e4) {
    u2 = t * (t * 4.7841e-8 - .0010817) + 7.1041;
    du2 = (t * 2 * 4.7841e-8 - .0010817) / u2;
  }
  u3 = 1;
  goto LEND;
 L13:
  u1 = y * (.27833 - y * (.047529 - y * .0030199)) + 5.2955;
  /* AL */
  du1 = (.27833 - y * (.095058000000000004 - y * 3 * .0030199)) * .001 
    / u1;
  /* Computing MAX */
  u2 = std::max(1.0,t * 3.245e-5 + .725);
  if (u2 != 1.0) {
    du2 = 3.245e-5 / u2;
  }
  if (t > 22400) {
    u2 = t * (t * 1.485e-7 - .005987) + 61.06;
    du2 = (t * 2 * 1.485e-7 - .005987) / u2;
  }
  /* Computing MAX */
  u3 = std::max(2.0, t * 3.43e-6 + 1.976);
  if (u3 != 2) {
    du3 = 3.43e-6 / u3;
  }
  if (t > 18140) {
    u3 = t * (t * 4.382e-9 - 1.59e-4) + 3.522;
    du3 = (t * 2 * 4.382e-9 - 1.59e-4) / u3;
  }
  goto LEND;
 L14:
  u1 = y * (y * (y * (.013109 - y * 6.2013e-4) - .11622) + .86319) + 
    6.7868;
  /* SI */
  du1 = (y * (y * (.039327000000000001 - y * 4 * 6.2013e-4) - 
	      .23244000000000001) + .86319) * .001 / u1;
  if (t > 10400) {
    u1 = t * (t * 7.282e-7 - .01465) + 86.01;
    du1 = (t * 2 * 7.282e-7 - .01465) / u1;
  }
  u2 = t * 4e-5 + 5.47;
  du2 = 4e-5 / u2;
  if (t > 1.8e4) {
    u2 = t * (t * 6.188e-8 - .00222) + 26.44;
    du2 = (t * 2 * 6.188e-8 - .00222) / u2;
  }
  /* Computing MAX */
  u3 = std::max(1.0,t * 1.1e-5 + .911);
  if (u3 != 1.0) {
    du3 = 1.1e-5 / u3;
  }
  if (t > 33300) {
    u3 = t * (t * 2.617e-8 - .001408) + 19.14;
    du3 = (t * 2 * 2.617e-8 - .001408) / u3;
  }
  goto LEND;
 L15:
  u1 = y * (y * (.057306 - y * .0010381) - .22476) + 4.2251;
  /* P */
  du1 = (y * (.11461200000000001 - y * 3 * .0010381) - .22476) * .001 /
    u1;
  if (t > 6e3) {
    u1 = t * 5.2e-4 + 1.56;
    du1 = 5.2e-4 / u1;
  }
  u2 = y * (y * (y * (.071913 - y * .0035156) - .55371) + 2.2494) + 
    4.4151;
  du2 = (y * (y * (.21573900000000001 - y * 4 * .0035156) - 
	      1.1074200000000001) + 2.2494) * .001 / u2;
  if (t > 7250) {
    u2 = t * 5.38e-4 + 4.62;
    du2 = 5.38e-4 / u2;
  }
  u3 = t * 3.4e-5 + 5.595;
  du3 = 3.4e-5 / u3;
  goto LEND;
 L16:
  u1 = t * 2.15e-4 + 7.5;
  /* S */
  du1 = 2.1e-4 / u1;
  if (t > 11600) {
    u1 = t * (t * 2.125e-7 - .004906) + 38.76;
    du1 = (t * 2 * 2.125e-7 - .004906) / u1;
  }
  u2 = t * 2.43e-4 + 2.845;
  du2 = 2.43e-4 / u2;
  if (t > 10500) {
    u2 = t * (t * 1.323e-8 - 1.68e-4) + 6.406;
    du2 = (t * 2 * 1.323e-8 - 1.68e-4) / u2;
  }
  u3 = t * 1.88e-4 + 7.38;
  du3 = 1.88e-4 / u3;
  goto LEND;
 L17:
  u1 = t * 6e-5 + 5.2;
  /* CL */
  du1 = 6e-5 / u1;
  if (t > 18400) {
    u1 = t * .0048 - 81.6;
    du1 = .0048 / u1;
  }
  u2 = t * 2.43e-4 + 7;
  du2 = 2.43e-4 / u2;
  u3 = t * 2.62e-4 + 2.2;
  du3 = 2.62e-4 / u3;
  goto LEND;
 L18:
  u1 = 1;
  /* AR */
  u2 = t * 3.8e-5 + 5.2;
  du2 = 3.8e-5 / u2;
  u3 = t * 1.554e-4 + 7.474;
  du3 = 1.554e-4 / u3;
  goto LEND;
 L19:
  u1 = y * (.023169 - y * (.017432 - y * .0040938)) + 1.9909;
  /* K */
  du1 = (.023169 - y * (.034863999999999999 - y * 3 * .0040938)) * 
    .001 / u1;
  if (t > 5800) {
    u1 = t * .002124 - 9.93;
    du1 = .002124 / u1;
  }
  u2 = 1;
  u3 = t * 1.93e-5 + 5.304;
  du3 = 1.93e-5 / u3;
  goto LEND;
 L20:
  u1 = exp(-1.731273 - x * (x * (x * (x * .508553 + 1.326861) + 
				  1.645456) + 5.004556)) + 1;
  du1 = -dx * (x * (x * (x * 2.0342120000000001 + 3.9805830000000002) + 
		    3.2909120000000001) + 5.004556);
  du1 = du1 * (u1 - 1) / u1;
  u2 = exp(-1.582112 - x * (x * (x * .539672 + 1.890737) + 3.996089)) 
    + 2;
  du2 = -dx * (x * (x * 3 * .539672 + 3.7814739999999998) + 3.996089);
  du2 = du2 * (u2 - 2) / u2;
  u3 = 1;
  goto LEND;
 L21:
  u1 = exp(x * (x * (x * .517796 + 1.173504) - 1.2392) + 2.071563) + 
    4;
  /* SC */
  du1 = dx * (x * (x * 3 * .517796 + 2.3470080000000002) - 1.2392);
  du1 = du1 * (u1 - 4) / u1;
  u2 = exp(x * (x * .054658 - .596238) + 2.988362) + 3;
  du2 = dx * (x * 2 * .054658 - .596238);
  du2 = du2 * (u2 - 3) / u2;
  u3 = 10;
  /* APPROXIMATELY */
  goto LEND;
 L22:
  u1 = exp(x * (x * (x * .278963 + .799613) - 1.227798) + 3.200453) + 
    5;
  /* TI */
  du1 = dx * (x * (x * 3 * .278963 + 1.599226) - 1.227798);
  du1 = du1 * (u1 - 5) / u1;
  if (t < 5500) {
    u1 = t * (t * 5.819e-7 - 2.838e-4) + 16.37;
    du1 = (t * 2 * 5.819e-7 - 2.838e-4) / u1;
  }
  u2 = exp(x * (x * .115693 - .551431) + 3.94529) + 4;
  du2 = dx * (x * 2 * .115693 - .551431);
  du2 = du2 * (u2 - 4) / u2;
  u3 = t * 8.5e-4 + 16.4;
  du3 = 8.5e-4 / u3;
  goto LEND;
 L23:
  u1 = exp(x * (x * (x * .1622 + .724694) - .906352) + 3.769611) + 4;
  /* V */
  du1 = dx * (x * (x * 3 * .1622 + 1.4493879999999999) - .906352);
  du1 = du1 * (u1 - 4) / u1;
  u2 = exp(x * (x * .21043 - .757371) + 3.755917) + 1;
  du2 = dx * (x * 2 * .21043 - .757371);
  du2 = du2 * (u2 - 1) / u2;
  u3 = t * .0103 - 18;
  du3 = .0103 / u3;
  if (t < 2250) {
    u3 = t * .0024;
    du3 = .0024 / u3;
  }
  goto LEND;
 L24:
  u1 = exp(x * (x * (x * .09527 + .154709) - 2.923459) + 1.225042) + 
    7;
  /* CR */
  du1 = dx * (x * (x * 3 * .09527 + .30941800000000003) - 2.923459);
  du1 = du1 * (u1 - 7) / u1;
  u2 = exp(.128752 - x * (x * (x * .230073 + 1.096548) + 4.143973)) + 
    6;
  du2 = -dx * (x * (x * 3 * .230073 + 2.1930960000000002) + 4.143973);
  du2 = du2 * (u2 - 6) / u2;
  u3 = t * .0021 + 10.4;
  du3 = .0021 / u3;
  goto LEND;
 L25:
  u1 = exp(-.86963 - x * (x * (x * (x * .265557 + 1.061055) + 2.13632) 
			   + 5.531252)) + 6;
  /* MN */
  du1 = -dx * (x * (x * (x * 1.062228 + 3.1831650000000002) + 4.27264) 
	       + 5.531252);
  du1 = du1 * (u1 - 6) / u1;
  u2 = exp(-.282961 - x * (x * (x * .159822 + .814675) + 3.77279)) + 
    7;
  du2 = -dx * (x * (x * 3 * .159822 + 1.6293500000000001) + 3.77279);
  du2 = du2 * (u2 - 7) / u2;
  u3 = 10;
  /* APPROXIMATELY */
  goto LEND;
 L26:
  u1 = exp(x * (x * (x * .118218 + .76027) - .979745) + 2.930047) + 
    9;
  /* FE */
  du1 = dx * (x * (x * 3 * .118218 + 1.52054) - .979745);
  du1 = du1 * (u1 - 9) / u1;
  if (t < 4e3) {
    u1 = t * (t * 2.04e-7 + .001306) + 15.85;
    du1 = (t * 2 * 2.04e-7 + .001306) / u1;
  }
  if (t > 9e3) {
    u1 = t * (t * 1.2477e-6 - .0095922) + 39.149;
    du1 = (t * 2 * 1.2477e-6 - .0095922) / u1;
  }
  u2 = exp(x * (x * .280982 - .612094) + 3.501597) + 10;
  du2 = dx * (x * 2 * .280982 - .612094);
  du2 = du2 * (u2 - 10) / u2;
  if (t > 1.8e4) {
    u2 = t * (t * 5.1567e-7 - .0061104) + 68.356;
    du2 = (t * 2 * 5.1567e-7 - .0061104) / u2;
  }
  u3 = t * (t * 5.7514e-8 + 5.5048e-4) + 17.336;
  du3 = (t * 2 * 5.7514e-8 + 5.5048e-4) / u3;
  goto LEND;
 L27:
  u1 = t * .0049 + 8.65;
  /* CO */
  du1 = .0049 / u1;
  u2 = t * .00358 + 11.2;
  du2 = .00358 / u2;
  u3 = t * .00142 + 15;
  du3 = .00142 / u3;
  goto LEND;
 L28:
  u1 = exp(x * (x * (.077498 - x * .278468) - .401323) + 3.084552) + 
    9;
  /* NI */
  du1 = dx * (x * (.15499599999999999 - x * 3 * .278468) - .401323);
  du1 = du1 * (u1 - 9) / u1;
  u2 = exp(1.593047 - x * (x * .115654 + 1.528966)) + 6;
  du2 = -dx * (x * 2 * .115654 + 1.528966);
  du2 = du2 * (u2 - 6) / u2;
  u3 = t * 6.9e-4 + 13.3;
  du3 = 6.9e-4 / u3;
  goto LEND;
 L29:
  /* Computing MAX */
  u1 = std::max(2.0,t * 1.51e-4 + 1.5);
  /* CU */
  if (u1 != 2) {
    du1 = 1.51e-4 / u1;
  }
  if (t > 6250) {
    u1 = t * 4.58e-4 - .3;
    du1 = 4.58e-4 / u1;
  }
  /* Computing MAX */
  r__1 = 1, r__2 = t * 1.49e-4 + .22;
  u2 = std::max(1.0,t * 1.49e-4 + .22);
  if (u2 != 1.0) {
    du2 = 1.49e-4 / u2;
  }
  u3 = t * 9.4e-5 + 8.025;
  du3 = 9.4e-5 / u3;
  goto LEND;
 L30:
  /* Computing MAX */
  u1 = std::max(1.0,t * 5.11e-5 + .632);
  /* ZN */
  if (u1 != 1) {
    du1 = 5.11e-5 / u1;
  }
  u2 = 2;
  u3 = 1;
  goto LEND;
 L31:
  u1 = y * (y * (y * (.054876 - y * .0025054) - .4643) + 1.9338) + 
    1.7931;
  /* GA */
  du1 = (y * (y * (.164628 - y * 4 * .0025054) - .92859999999999998) + 
	 1.9338) * .001 / u1;
  if (t > 6e3) {
    u1 = t * 2.03e-4 + 4.18;
    du1 = 2.03e-4 / u1;
  }
  u2 = 1;
  u3 = 2;
  goto LEND;
 L32:
  u1 = t * 4.08e-4 + 6.12;
  /* GE */
  du1 = 4.08e-4 / u1;
  u2 = t * 1.78e-4 + 3.445;
  du2 = 1.78e-4 / u2;
  u3 = 1.1;
  /* APPROXIMATELY */
  goto LEND;
 L33:
  u1 = t * 3.65e-4 + 2.65;
  /* AS */
  du1 = 3.65e-4 / u1;
  u2 = y * (y * (y * (.030408 - y * .0011609) - .33383) + 2.284) - 
    .25384;
  du2 = (y * (y * (.091224 - y * 4 * .0011609) - .66766000000000003) + 
	 2.284) * .001 / u2;
  if (t > 1.2e4) {
    u2 = 8;
    du2 = 0;
  }
  u3 = 8;
  /* APPROXIMATELY */
  goto LEND;
 L34:
  u1 = t * 1.71e-4 + 6.34;
  /* SE */
  du1 = 1.71e-4 / u1;
  u2 = y * (y * .032053 - .15392) + 4.1786;
  du2 = (y * .032053 * 2 - .15392) * .001 / u2;
  u3 = 8;
  /* APPROXIMATELY */
  goto LEND;
 L35:
  u1 = t * 1.12e-4 + 4.12;
  /* BR */
  du1 = 1.12e-4 / u1;
  u2 = t * 3.08e-4 + 5.22;
  du2 = 3.08e-4 / u2;
  u3 = t * 2.86e-4 + 2.3;
  du3 = 2.86e-4 / u3;
  goto LEND;
 L36:
  u1 = 1;
  /* KR */
  u2 = t * 7.4e-5 + 4.11;
  du2 = 7.4e-5 / u2;
  u3 = t * 2.23e-4 + 5.35;
  du3 = 2.23e-4 / u3;
  goto LEND;
 L37:
  /* Computing MAX */
  u1 = std::max(2.0,t * 1.94e-4 + 1.38);
  /* RB */
  if (u1 != 2) {
    du1 = 1.94e-4 / u1;
  }
  if (t > 6250) {
    u1 = t * .00279 - 14.9;
    du1 = .00279 / u1;
  }
  u2 = 1;
  u3 = t * 4.85e-5 + 4.207;
  du3 = 4.85e-5 / u3;
  goto LEND;
 L38:
  u1 = y * (y * (y * (.021424 - y * .0010231) - .10746) + .20148) + 
    .87127;
  /* SR */
  du1 = (y * (y * (.064271999999999996 - y * .0010231 * 4) - .21492) 
	 + .20148) * .001;
  du1 /= u1;
  if (t > 6500) {
    u1 = t * .001224 - 6.12;
    du1 = .001224 / u1;
  }
  /* Computing MAX */
  u2 = std::max(2.0,t * 2.6e-4 + .84);
  if (u2 != 2.0) {
    du2 = 2.6e-4 / u2;
  }
  u3 = 1;
  /* APPROXIMATELY */
  goto LEND;
 L39:
  u1 = t * .00258 + .2;
  /* Y */
  u2 = t * .001855 + 7.15;
  u3 = t * 9.9e-5 + 9.71;
  goto LEND;
 L40:
  u1 = t * (t * 2.199e-6 - .01866) + 76.31;
  /* ZR */
  if (t < 6236) {
    u1 = t * (t * 5.386e-7 + .002806) + 6.8;
  }
  u2 = exp(3.721329 - x * .906502) + 4;
  u3 = t * .001385 + 12.3;
  goto LEND;
 L41:
  /* Computing MAX */
  u1 = std::max(1.0,t * .0143 - 19);
  /* NB */
  u2 = t * .01015 - 4;
  u3 = 25;
  /* APPROXIMATELY */
  goto LEND;
 L42:
  /* Computing MAX */
  u1 = std::max(7.0,t * .0015 + 2.);
  /* MO */
  if (t > 7e3) {
    u1 = t * .00728 - 38.1;
  }
  u2 = t * .00117 + 1.25;
  if (t > 6900) {
    u2 = t * .00548 - 28.5;
  }
  u3 = t * 1.464e-4 + 24.04;
  goto LEND;
 L43:
  u1 = y * (y * (y * (y * (.048401 - y * .0021538) - .4078) + 1.6525) 
	    + .30648) + 4.439;
  /* TC */
  if (t > 6e3) {
    u1 = 24;
  }
  u2 = y * (y * (y * (y * (.049656 - y * .0019087) - .502) + 2.369) - 
	    2.963) + 8.1096;
  if (t > 6e3) {
    u2 = 17;
  }
  u3 = 220;
  /* APPROXIMATELY */
  goto LEND;
 L44:
  u1 = t * .00717 - 3;
  /* RU */
  u2 = t * .00426 + 3;
  u3 = 22;
  goto LEND;
 L45:
  u1 = y * (y * (.043125 - y * (.0087907 - y * 5.9589e-4)) + 3.8468) + 
    6.9164;
  /* RH */
  u2 = y * (y * (y * (y * 2.1218e-4 + .002014) - .038257) + 1.7476) + 
    7.2902;
  u3 = 30;
  /* APPROXIMATELY */
  goto LEND;
 L46:
  /* Computing MAX */
  u1 = std::max(1.0,t * 9.86e-4 - 1.75);
  /* PD */
  u2 = t * 3.62e-4 + 5.6;
  u3 = 20;
  /* APPROXIMATELY */
  goto LEND;
 L47:
  /* Computing MAX */
  u1 = std::max(2.0,t * 7.88e-5 + 1.537);
  /* AG */
  /* Computing MAX */
  u2 = std::max(1.0,t * 3.4e-5 + .73);
  u3 = t * 1.248e-4 + 6.773;
  goto LEND;
 L48:
  /* Computing MAX */
  u1 = std::max(1.0,t * 7.6e-5 + .43);
  /* CD */
  u2 = 2;
  u3 = 1;
  goto LEND;
 L49:
  u1 = t * 3.92e-4 + 2.16;
  /* IN */
  u2 = 1;
  u3 = 2;
  goto LEND;
 L50:
  u1 = t * 6.16e-4 + 2.14;
  /* SN */
  u2 = t * 2.27e-4 + 2.06;
  u3 = 1.05;
  /* APPROXIMATELY */
  goto LEND;
 L51:
  u1 = t * 4.86e-4 + 2.34;
  /* SB */
  u2 = t * 5.36e-4 + .69;
  u3 = 3.5;
  /* APPROXIMATELY */
  goto LEND;
 L52:
  u1 = t * 4.56e-4 + 3.948;
  /* TE */
  u2 = y * (y * (.06939 - y * .0024271) - .25894) + 4.2555;
  if (t > 1.2e4) {
    u2 = 7;
  }
  u3 = 5;
  /* APPROXIMATELY */
  goto LEND;
 L53:
  /* Computing MAX */
  u1 = std::max(4.0,t * 9.5e-5 + 3.8);
  /* I */
  u2 = t * 3e-4 + 4.12;
  u3 = 7;
  /* APPROXIMATELY */
  goto LEND;
 L54:
  u1 = 1;
  /* XE */
  u2 = t * 6.876e-5 + 3.75;
  u3 = t * 2.323e-4 + 4.121;
  goto LEND;
 L55:
  /* Computing MAX */
  r__1 = 2, r__2 = t * 1.67e-4 + 1.56;
  u1 = std::max(r__1,r__2);
  /* CS */
  if (t > 4850) {
    u1 = t * .00104 - 2.68;
  }
  u2 = 1;
  u3 = t * 4.971e-5 + 3.769;
  goto LEND;
 L56:
  /* Computing MAX */
  r__1 = 1, r__2 = t * 9.85e-4 - 1.8;
  u1 = std::max(r__1,r__2);
  /* BA */
  if (u1 != 1) {
    du1 = 9.85e-4 / u1;
  }
  if (t > 6850) {
    u1 = t * .00308 - 16.2;
    du1 = .00308 / u1;
  }
  u2 = t * 5.94e-4 + 1.11;
  du2 = 5.94e-4 / u2;
  u3 = 1;
  goto LEND;
 L57:
  u1 = t * 9.5e-4 + 15.42;
  /* LA */
  if (t > 5060) {
    u1 = t * .0038 + 1;
  }
  u2 = t * .00356 + 13.2;
  u3 = 12;
  /* APPROXIMATELY */
  goto LEND;
 L58:
  u1 = exp(x * (x * (x * .179675 + .119673) - 1.98399) + 5.202903) + 
    9;
  /* CE */
  u2 = exp(5.634882 - x * (x * (x * .052221 + .310515) + 1.459196)) + 
    8;
  u3 = exp(3.629123 - x * (x * (x * (.03186 - x * .014676) + .372409) 
			    + 1.340945)) + 9;
  goto LEND;
 L59:
  u2 = exp(4.32396 - x * (x * (x * .028999 + .149498) + 1.191467)) + 
    9;
  /* PR */
  u1 = u2;
  /* APPROXIMATELY */
  u3 = exp(x * (x * (x * .277916 + .489574) - 1.614554) + 3.206855) + 
    10;
  goto LEND;
 L60:
  u1 = exp(x * (x * (x * (x * .127326 + .50666) + .082258) - 2.779176) 
	   + 4.456882) + 9;
  /* ND */
  u2 = exp(x * (x * (x * (x * .038225 + .26392) + .17193) - 2.039946) 
	   + 4.689643) + 8;
  u3 = u2;
  /* APPROXIMATELY */
  goto LEND;
 L61:
  u1 = 20;
  /* PM APPROXIMATELY */
  u2 = 25;
  /* APPROXIMATELY */
  u3 = 100;
  /* APPROXIMATELY */
  goto LEND;
 L62:
  u1 = exp(x * (x * (x * .566263 + .9964) - 1.851549) + 3.549595) + 
    1;
  /* SM */
  u2 = exp(x * (x * (x * .161944 + .358695) - 1.418222) + 4.052404) + 
    2;
  u3 = exp(3.222807 - x * (x * (x * (x * .251011 + .533833) - .056205) 
			    + .699473)) + 1;
  goto LEND;
 L63:
  u1 = exp(1.024374 - x * (x * (x * (x * .286737 + .827789) + 1.540805)
			    + 4.533653)) + 8;
  /* EU */
  u2 = exp(x * (x * (x * .05684 + .379584) - 1.50646) + 1.92776) + 9;
  u3 = 8;
  /* APPROXIMATELY */
  goto LEND;
 L64:
  u1 = exp(x * (x * (x * .388845 + .800411) - 1.583513) + 4.009587) + 
    5;
  /* GD */
  u2 = exp(4.362107 - x * (x * (x * (x * .055475 + .076453) - .074813) 
			    + 1.208124)) + 6;
  u3 = exp(3.412951 - x * (x * (.042489 - x * .004017) + .50271)) + 
    5;
  goto LEND;
 L65:
  u1 = exp(x * (x * (x * .240203 + .570094) - 1.249355) + 4.791661) + 
    16;
  /* TB */
  u2 = exp(4.472549 - x * (x * (x * .131631 + .00588) + .295965)) + 
    15;
  u3 = u2;
  /* APPROXIMATELY */
  goto LEND;
 L66:
  u1 = exp(3.029646 - x * (x * (.086671 - x * .216214) + 3.121036)) + 
    17;
  /* DY */
  u2 = exp(3.465323 - x * (x * (x * (x * .303575 + .431447) - .382265) 
			    + 1.27062)) + 18;
  u3 = u2;
  /* APPROXIMATELY */
  goto LEND;
 L67:
  u3 = exp(1.610084 - x * (x * (.133139 - x * .071196) + 2.373926)) + 
    16;
  /* HO */
  u1 = u3;
  u2 = u3;
  /* APPROX. */
  goto LEND;
 L68:
  u1 = exp(2.895648 - x * (x * (x * (x * .095813 + .215267) + .561515) 
			    + 2.968603)) + 13;
  /* ER */
  u2 = exp(3.202542 - x * (x * (x * (x * .186042 + .343738) - .226622) 
			    + .852209)) + 14;
  u3 = u2;
  /* APPROX. */
  goto LEND;
 L69:
  u1 = exp(1.021172 - x * (x * (x * .034811 + 1.081603) + 4.94757)) + 
    8;
  /* TM */
  u2 = exp(x * (x * (x * .813303 + 1.940395) - 1.295327) + 2.173152) + 
    9;
  u3 = exp(x * (x * (x * .554397 + .799911) - 3.383369) - .567398) + 
    8;
  goto LEND;
 L70:
  u1 = exp(-2.350549 - x * (x * (x * .269237 + 1.93869) + 6.688837)) + 
    1;
  /* YB */
  u2 = exp(-3.047465 - x * (x * (x * .44757 + 2.355267) + 7.390444)) + 
    2;
  u3 = exp(-6.192056 - x * (x * (x * .940171 + 4.579385) + 10.560552)) 
    + 1;
  goto LEND;
 L71:
  u1 = exp(x * (x * (x * .193362 + .608536) - 1.140264) + 1.537094) + 
    4;
  /* LU */
  /* Computing MAX */
  r__1 = 1, r__2 = t * 1.52e-4 + .66;
  u2 = std::max(r__1,r__2);
  if (t > 5250) {
    u2 = t * 4.86e-4 - 1.09;
  }
  u3 = 5;
  /* APPROXIMATELY */
  goto LEND;
 L72:
  u1 = y * (y * (.57862 - y * (.072887 - y * .0036848)) + .407) + 
    4.1758;
  /* HF */
  u2 = t * .003095 - 2.979;
  u3 = 30;
  /* APPROXIMATELY */
  goto LEND;
 L73:
  u1 = y * (y * (y * (y * 3.0739e-4 + .0074861) + .34936) + .81776) + 
    3.0679;
  /* TA */
  u2 = y * (y * (.56443 - y * (.031036 - y * 8.9565e-4)) + 2.0103) + 
    1.6834;
  u3 = 15;
  goto LEND;
 L74:
  u1 = y * (y * (y * (y * (.041924 - y * .00184) - .34373) + 1.4433) - 
	    .25057) + .3951;
  /* W */
  if (t > 1.2e4) {
    u1 = 23;
  }
  u2 = y * (y * (.3303 - y * (.0084971 - y * 5.5794e-4)) + 1.0396) + 
    1.055;
  u3 = 20;
  goto LEND;
 L75:
  u1 = y * (y * (y * (.09075 - y * .0039331) - .42096) + .72721) + 
    5.5671;
  /* RE */
  if (t > 1.2e4) {
    u1 = 29;
  }
  u2 = y * (y * (y * (.050724 - y * .0018544) - .28532) + .59999) + 
    6.5699;
  if (t > 1.2e4) {
    u2 = 22;
  }
  u3 = 20;
  goto LEND;
 L76:
  u1 = y * (y * (.68181 - y * (.044252 - y * .0019975)) - .32516) + 
    8.6643;
  /* OS */
  u2 = y * (y * (.65292 - y * (.064984 - y * .0028792)) - .3814) + 
    9.7086;
  u3 = 10;
  goto LEND;
 L77:
  u1 = y * (y * (y * (y * (.033511 - y * .0013376) - .34389) + 1.9388) 
	    - 2.412) + 11.07;
  /* IR */
  if (t > 1.2e4) {
    u1 = 30;
  }
  u2 = 15;
  u3 = 20;
  goto LEND;
 L78:
  u1 = t * .00127 + 16.4;
  /* PT */
  du1 = .00127 / u1;
  u2 = y * (y * (.57234 - y * (.061219 - y * .0026878)) - 1.0363) + 
    6.5712;
  du2 = (y * (1.1446799999999999 - y * (.18365700000000001 - y * 4 * 
					 .0026878)) - 1.0363) * .001;
  du2 /= u2;
  u3 = 15;
  goto LEND;
 L79:
  u1 = t * 2.79e-4 + 1.24;
  /* AU */
  du1 = 2.79e-4 / u1;
  u2 = y * (y * (y * .0016586 + .0028439) - .040809) + 1.0546;
  du2 = (y * (y * .0016586 * 3 + .0056877999999999998) - .040809) * 
    .001 / u2;
  u3 = 7;
  goto LEND;
 L80:
  u1 = 1;
  /* HG */
  u2 = 2;
  /* Computing MAX */
  r__1 = 1, r__2 = t * 3.976e-5 + .669;
  u3 = std::max(r__1,r__2);
  if (u3 != 1) {
    du3 = 3.976e-5 / u3;
  }
  goto LEND;
 L81:
  /* Computing MAX */
  r__1 = 2, r__2 = t * 3.35e-4 + .63;
  u1 = std::max(r__1,r__2);
  /* TL */
  if (u1 != 2) {
    du1 = 3.35e-4 / u1;
  }
  u2 = 1;
  u3 = 2;
  goto LEND;
 L82:
  /* Computing MAX */
  r__1 = 1, r__2 = t * 2.35e-4 + .42;
  u1 = std::max(r__1,r__2);
  /* PB */
  if (u1 != 1) {
    du1 = 2.35e-4 / u1;
  }
  if (t > 6125) {
    u1 = t * 5e-4 - 1.2;
    du1 = 5e-4 / u1;
  }
  /* Computing MAX */
  r__1 = 2, r__2 = t * 7.9e-5 + 1.72;
  u2 = std::max(r__1,r__2);
  if (u2 != 2) {
    du2 = 7.9e-5 / u2;
  }
  u3 = 1;
  goto LEND;
 L83:
  u1 = t * 2.87e-4 + 2.78;
  /* BI */
  du1 = 2.87e-4 / u1;
  /* Computing MAX */
  r__1 = 1, r__2 = t * 1.41e-4 + .37;
  u2 = std::max(r__1,r__2);
  if (u2 != 1) {
    du2 = 1.41e-4 / u2;
  }
  u3 = 2.5;
  /* APPROXIMATELY */
  goto LEND;
 L84:
  u1 = 5;
  /* PO */
  u2 = 5;
  u3 = 4;
  goto LEND;
 L85:
  u1 = 4;
  /* AT */
  u2 = 6;
  u3 = 6;
  goto LEND;
 L86:
  u1 = 1;
  /* RN */
  u2 = 4;
  u3 = 6;
  goto LEND;
 L87:
  u1 = 2;
  /* FR */
  u2 = 1;
  u3 = 4.5;
  goto LEND;
 L88:
  u1 = 1;
  /* RA */
  u2 = 2;
  u3 = 1;
  goto LEND;
 L89:
  u1 = 6;
  /* AC */
  u2 = 3;
  u3 = 7;
  goto LEND;
 L90:
  u1 = 8;
  /* TH */
  u2 = 8;
  u3 = 8;
  goto LEND;
 L91:
  u1 = 50;
  /* PA */
  u2 = 50;
  u3 = 50;
  goto LEND;
 L92:
  u1 = 25;
  /* U */
  u2 = 25;
  u3 = 25;
  goto LEND;

 LEND:
  u1_in = (T)u1, u2_in = (T)u2, u3_in = (T)u3;
  du1_in = (T)du1, du2_in = (T)du3, du3_in = (T)du3;
  return;  
} /* partition_f__ */

template void  witt::partition_f(int, double, double &, double &, double &, double &, double &, double & );
template void   witt::partition_f(int, float, float &, float &, float &, float &, float &, float & );


void witt::witt::hydrostatic(int ndep, double *tau, double *t, double *Pg, double *rho, double *nel,
		       double *pel, double *z, double *cmass, double pgas_bound, float tol){


  string inam = "witt::hydrostatic: ";
  int maxiter = 50, imax = 0;

  
  // Init vars
  int nw = 1;
  double wav = 5000.0;
  double kappa = 0.0, scat = 0.0, kappa_old = 0.0;

  // Init boundary
  Pg[0] = pgas_bound;
  z[0] = 0.0;
  
  rho[0] = rho_from_pg<double>(t[0], Pg[0], pel[0]);
  contOpacity<double>(t[0], Pg[0], pel[0], nw, &wav, &kappa_old);
  kappa_old /= rho[0];
  cmass[0] = Pg[0] / gravity; //(xna + xne) * (kb * t[0] / gravity);
  nel[0] = pel[0]/(phyc::BK*t[0]);

  // Loop height
  for(int k = 1; k < ndep; k++){
    
    double dtau = tau[k] - tau[k-1];
    double dif = 1e10;
    int iter = 0;
    kappa = kappa_old;
    
    /* --- Iterate because kappa 
       depends on Pg and vice-versa 
       --- */
    while((dif > tol) && (iter < maxiter)){
      
      dif = Pg[k]; // Store old value

      /* --- Integrate kappa assuming linear dependence with tau --- */
      
      // if(iter == 0){
      //	Pg[k] = Pg[k-1] + gravity * dtau  / (kappa_old);
      //}else{
	///Pg[k] = Pg[k-1] + gravity * dtau / (kappa - kappa_old) * log(kappa/kappa_old);
      	Pg[k] = Pg[k-1] + gravity * 2.0 * dtau / (kappa + kappa_old) ;//* log(kappa/kappa_old);

	//}
      
      /* ---  Get opacity --- */
      rho[k] = rho_from_pg<double>(t[k], Pg[k], pel[k]);
      contOpacity<double>(t[k], Pg[k], pel[k], nw, &wav, &kappa);
      
      //contOpacity_TPg(t[k], Pg[k], nw, &wav, &kappa, &scat);
      kappa /= rho[k]; // Convert to [cm^2 / g]
      

      /* --- Get relative change --- */
      
      dif = abs(2.0 * (dif-Pg[k]) / ( dif+Pg[k]));
      iter++;
    }
    if((iter-1)> imax) imax = iter-1;
    
    /* --- Fill in arrays --- */
    
    //nel[k] = xne; 
    //rho[k] = RHOest;
    // pel[k] = bk * xne * t[k];

    nel[k] = pel[k]/(phyc::BK*t[k]);

    
    /* --- Compute z-scale and cmass--- */
    
    z[k] = z[k-1] - 2.0 * dtau / (kappa * rho[k] + kappa_old * rho[k-1]);
    cmass[k] = cmass[k-1] + 2.0 * dtau / (kappa + kappa_old);

    /* --- Store partial pressures and partition function 
       for the ratiative transfer
       --- */
    //store_partial_pressures(ndep, k, xna, xne);
    kappa_old = kappa; // Store for next iteration


  }
  //double toff = exp(2.0 * log(cmass[1]) - log(cmass[2]));
  for(int kk=0; kk<ndep; kk++) cmass[kk] = log10(cmass[kk]);
  
  // cerr<<"ceos::hydrostatic: ITMAX="<<imax<<endl;
}


void witt::witt::fill_densities(int ndep, double *t, double *pgas, double *rho, double *pel,
			  double *nne, int touse,  int keep_nne, float tol){

  /* --- touse switches which of the following to use to compute the rest:
     (0) pgas, (1) rho, (2) pel, (3) nne.
     
     Store partial pressures and partition functions for the RT calculations
     --- */
  if((keep_nne == 0) && ((touse == 0) || (touse == 1))){
    if(touse == 0){
      for(int k = 0; k<ndep; k++){
	nne[k] = nne_from_T_Pg(t[k], pgas[k], rho[k]);
	pel[k] = phyc::BK*nne[k]*t[k];
	//store_partial_pressures(ndep, k, xna, xne);
      }
    }else if(touse == 1){
      for(int k = 0; k<ndep; k++){
	nne[k] = nne_from_T_rho(t[k], pgas[k], rho[k]);
	pel[k] = phyc::BK*nne[k]*t[k];
	//store_partial_pressures(ndep, k, xna, xne);
      }
    }  
    else{
      cerr << "ceos::fill_densities: ERROR, touse["<<touse<<"] takes values between 0-3"<<endl;
    }
  }else{
    if(touse == 0){
      for(int k = 0; k<ndep; k++){
	nne_from_T_Pg_nne(t[k], pgas[k], rho[k], nne[k]);
	pel[k] = phyc::BK*nne[k]*t[k];
	//store_partial_pressures(ndep, k, xna, xne);
      }
    }else if(touse == 1){
      for(int k = 0; k<ndep; k++){
	nne_from_T_rho_nne(t[k], pgas[k], rho[k], nne[k]);
	pel[k] = phyc::BK*nne[k]*t[k];
	//	store_partial_pressures(ndep, k, xna, xne);	    
      }
    } else if(touse == 2){
      for(int k = 0; k<ndep; k++){
	rho[k] = rho_from_T_pel(t[k], pgas[k], pel[k], tol);
	nne[k] = pel[k]/(t[k]*phyc::BK);
	//store_partial_pressures(ndep, k, xna, xne);
      }
    } else if(touse == 3){
      for(int k = 0; k<ndep; k++){
	rho[k] = rho_from_T_nne(t[k], pgas[k], nne[k], tol);
	pel[k] = phyc::BK*nne[k]*t[k];
	//store_partial_pressures(ndep, k, xna, xne);
      }
    } else{
      cerr << "ceos::fill_densities: ERROR, touse["<<touse<<"] takes values between 0-3"<<endl;
    }
  }
}

void witt::witt::initAbundances(vector<iabund> &ab, bool verbose)
{

  /* --- Init default abundances --- */
  
  for(int ii=0; ii<MAX_ELEM; ii++){
    ABUND[ii] = pow(10., ABUND_default[ii]);
    ABUND[ii] /= ABUND[0];
  }
  
  /* --- replace abunds --- */
  
  for(auto &it: ab){
    
    for(int ii=0; ii<MAX_ELEM; ii++){
      if(!strcmp(phyc::ELEMEN[ii], it.elem)){
	ABUND[ii] = pow(10., it.abund);
	if(verbose) fprintf(stderr, "witt::witt::initAbundances: Changed [%s] -> %7.3f\n", it.elem, it.abund);
	break;
      }
    } // ii
  } // it


  double sum = 0.0;
  for(int ii = 0; ii<MAX_ELEM; ii++) sum += ABUND[ii];
  
  //  for(int ii = 0; ii<MAX_ELEM; ii++) ABUND[ii] /= ABUND[];
  tABUND = sum;//sum;///ABUND[0];

   
  //
  // Estimate average molecular weight, excluding electrons?
  //
  double wsum = 0;
  for(int ii = 0; ii<MAX_ELEM; ii++){
    wsum += phyc::AMASS[ii] * ABUND[ii];
  }

  avweight = wsum / sum;
}

void witt::witt::readAbund(std::string file)
{

  /* --- Open file and read abundances --- */

  std::vector<iabund> ab;
  iabund ia;
  int Nread;
  
  std::ifstream in(file, std::ios::in | std::ios::binary);
  std::string iline;
  if(in){
    while (std::getline(in, iline)) {
      
      if(iline[0] == '#' || iline == "") continue;
      iline.erase(0,2); // Remove dummy spaces
      int n = iline.find(" ");      
      Nread = sscanf(iline.substr(0, n).c_str(), "%2s", ia.elem);
      
      if(ia.elem[1] != '\0' && ia.elem[1] != ' ') ia.elem[1] = tolower(ia.elem[1]);
      else ia.elem[1] = ' ';
      ia.elem[2] = '\0';
      
      ia.abund = atof(iline.substr(n+1).c_str()) - 12.0;
      ab.push_back(ia);    
    }//while
  }else
    fprintf(stderr,"ceos::readAbund: WARNING: file [%s] not found, using default abundance values\n",
	    file.c_str());

  /* --- Update abundances --- */
  
  initAbundances(ab, false);
  
}


void eos::witt::hydrostatic_cmass(int ndep, double *tau, double *t, double *Pg, double *rho, double *nel,
				  double *z, double *cmass, double *ltau, double &pgas_bound){

  
  const string inam = "witt::hydrostatic_cmass: ";
  double wav = 5000.0, pel; 
  double dum, kappa_old, kappa, cm, ocm;
  int  nw = 1;
  
  
  /* --- Init vars --- */
  cm = pow(10.0, cmass[0]);
  Pg[0] = gravity * cm;
  
  rho[0] = rho_from_pg<double>(t[0], Pg[0], pel);
  contOpacity<double>(t[0], Pg[0], pel, nw, &wav, &kappa_old);
  kappa_old /= rho[0];
  cmass[0] =  Pg[0] / gravity; //(xna + xne) * (kb * t[0] / gravity);
  nel[0] = pel/(phyc::BK*t[0]);
  tau[0] = Pg[0] * kappa / (rho[0]*gravity);//0.0; //kappa / rho[0] * cm, ltau[0] = log10(tau[0]);
  z[0] = 0.0;


  for(int k = 1; k<ndep; k++){
    kappa_old = kappa;
    ocm = cm;

    
    /* --- Solve Gas pressure --- */
    
    cm = pow(10.0, cmass[k]);
    Pg[k] = gravity * cm;


    /* --- Fill other scales and variables --- */
    rho[k] = rho_from_pg<double>(t[k], Pg[k], pel);
    contOpacity(t[k], Pg[k], pel, nw, &wav, &kappa);
    //rho[k] = RHOest, nel[k] = xne;
    nel[k] = pel/(t[k]*phyc::BK);
    
    double dz =  2.0 * (cm - ocm) / (rho[k-1] + rho[k]);
    z[k] = z[k-1] - dz;
    tau[k] = tau[k-1] + dz * 0.5 * (kappa + kappa_old), ltau[k] = log10(tau[k]);

    // store_partial_pressures(ndep, k, xna, xne);

  }

  double toff = exp(2.0 * log(tau[1]) - log(tau[2]));
	    
	    
  for(int k = 0; k < ndep; k++){
    //tau[k] += toff;
    ltau[k] = log10(tau[k]);
    //  fprintf(stderr,"%d %e %e %e\n", k, ltau[k], cmass[k], z[k]*1.e-5);
  }
  
  //   for(int k = 0; k<ndep; k++)
  //  fprintf(stderr,"[%3d] %f  %f  %f  %e  %e\n", k, cmass[k], ltau[k], z[k]*1.e-5, t[k], Pg[k]);

  //exit(0);

}
