/*
  CLASS CPROFILES
  
  PURPOSE: Compute line profiles given atom level populations
  
  AUTHOR(S): J. de la Cruz Rodriguez (ISP-SU 2015)
  
  MODIFICATIONS: 
        2015-12-27, JdlCR: Compute the Zeeman splitting only once and re-use them.
	                   Bug-fix, was using F were should use L = 2*F in the calculation
			   of the disperssion profile.

 */
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "cprofiles.h"
#include "input.h"
#include "physical_consts.h"
#include "depthmodel.h"
//
using namespace std;
using namespace phyc;
//

/* --- Factorials --- */
const double cprofiles::fact[31] =
  {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0,
   479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0, 20922789888000.0, 355687428096000.0,
   6402373705728000.0, 121645100408832000.0, 2432902008176640000.0, 51090942171709440000.0,
   1124000727777607680000.0, 25852016738884978212864.0, 620448401733239409999872.0,
   15511210043330986055303168.0, 403291461126605650322784256.0, 10888869450418351940239884288.0,
   304888344611713836734530715648.0, 8841761993739700772720181510144.0, 265252859812191032188804700045312.0};

/* --- Zeeman splitting constant --- */

const double cprofiles::LARMOR= phyc::EE / (4.0 * phyc::PI * phyc::ME * phyc::CC);

/* --- Identity matrix --- */

const double cprofiles::ident[4][4] = {{1.0,0.0,0.0,0.0}, {0.0,1.0,0.0,0.0},
				       {0.0,0.0,1.0,0.0}, {0.0,0.0,0.0,1.0}};


//----------------------------------------------------------------
// Constructor
//----------------------------------------------------------------
cprofiles::cprofiles(int nw, int ndep){
  nnw = 0;
  nndep = 0;
  init(nw, ndep);
}


void cprofiles::init(int nw, int ndep){

  /* --- do not reallocate if already allocated --- */
  if(nnw == nw && nndep == ndep)
    return;

  /* --- Otherwise allocate everything --- */
  nndep = ndep;
  nnw = nw;
  
  vector<int> dims = {ndep, nw};
  mki.set(dims);
  mkq.set(dims);
  mku.set(dims);
  mkv.set(dims);
  mfq.set(dims);
  mfu.set(dims);
  mfv.set(dims);
  
  ki.resize(ndep);
  kq.resize(ndep);
  ku.resize(ndep);
  kv.resize(ndep);
  fq.resize(ndep);
  fu.resize(ndep);
  fv.resize(ndep);
  
} 

//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
cprofiles::~cprofiles(){
  ki.clear();
  kq.clear();
  ku.clear();
  kv.clear();
  fq.clear();
  fu.clear();
  fv.clear();
}

//-------------------------------------------------------------------------
// Set to zero the terms of the Abs. Matrix at all depths
//-------------------------------------------------------------------------
void cprofiles::set_zero(){ 
  std::fill(ki.begin(), ki.end(), 0);
  std::fill(kq.begin(), kq.end(), 0);
  std::fill(ku.begin(), ku.end(), 0);
  std::fill(kv.begin(), kv.end(), 0);
  std::fill(fq.begin(), fq.end(), 0);
  std::fill(fu.begin(), fu.end(), 0);
  std::fill(fv.begin(), fv.end(), 0);
}


//-------------------------------------------------------------------------
// Zero opacity matrices
//-------------------------------------------------------------------------
void cprofiles::set_zero_abmat(){
  mki.zero();
  mkq.zero();
  mku.zero();
  mkv.zero();
  mfq.zero();
  mfu.zero();
  mfv.zero();
}

//-------------------------------------------------------------------------
// Compute Zeeman splitting
//-------------------------------------------------------------------------
void cprofiles::init_zeeman_components(line_t &line)
{
  /* --- 
     
     Compute the Zeeman splitting and strength for a line, 
     the Zeeman components are already scaled so the 
     sum{Red}=1, sum{Blue}=1 and sum{Par}=1
     
     --- */
  
  
  line.nZ = 0;
  int nup = int(2 * line.Jup) + 1; 
  int delta_j = (int)(line.Jup - line.Jlow);

  
  for(int iup = 1; iup <= nup; iup++){
    float Mup = line.Jup + 1 - iup;
    for(int ilow = 1; ilow <= 3; ilow++){
      float Mlow = Mup - 2 + ilow;
      if(abs(Mlow) <= line.Jlow){
	
	/* --- Compute relative Zeeman strength, 
	   Landi Degl'innocenti & Landolfi (2004), 
	   table 3.1 - pag. 81 --- */
	
	double strength = 0.0; // If abs(delta_j) > 1 then don't compute Zeeman splitting. 
	//
	if(delta_j == 1){
	  
	  if(ilow == 1)      strength = 1.5 * (line.Jlow + Mlow + 1.0) * (line.Jlow + Mlow + 2.0)
			       / ((line.Jlow+1.0)*(2.0*line.Jlow + 1.0) * (2.0 * line.Jlow + 3.0));
	  else if(ilow == 2) strength = 3.0 * (line.Jlow - Mlow + 1.0) * (line.Jlow + Mlow + 1.0)
			       / ((line.Jlow+1.0)*(2.0*line.Jlow + 1.0) * (2.0 * line.Jlow + 3.0));
	  else               strength = 1.5 * (line.Jlow - Mlow + 1.0) * (line.Jlow - Mlow + 2.0)
			       / ((line.Jlow+1.0)*(2.0*line.Jlow + 1.0) * (2.0 * line.Jlow + 3.0));
	  
	} else if(delta_j == 0){
	  
	  if(ilow == 1)      strength = 1.5 * (line.Jlow - Mlow) * (line.Jlow + Mlow + 1.0)
			       / (line.Jlow * (line.Jlow + 1.0) * (2.0 * line.Jlow + 1.0));
	  else if(ilow == 2) strength = 3.0 * Mlow * Mlow
			       / (line.Jlow * (line.Jlow + 1.0) * (2.0 * line.Jlow + 1.0));
	  else               strength = 1.5 * (line.Jlow + Mlow) * (line.Jlow - Mlow + 1.0)
			       / (line.Jlow * (line.Jlow + 1.0) * (2.0 * line.Jlow + 1.0));
	  
	} else if(delta_j == -1){
	  
	  if(ilow == 1)      strength = 1.5 * (line.Jlow - Mlow) * (line.Jlow - Mlow - 1.0)
			       / (line.Jlow * (2.0 * line.Jlow - 1.0) * (2.0 * line.Jlow + 1.0));
	  else if(ilow == 2) strength = 3.0 * (line.Jlow - Mlow) * (line.Jlow + Mlow)
			       / (line.Jlow * (2.0 * line.Jlow - 1.0) * (2.0 * line.Jlow + 1.0));
	  else               strength = 1.5 * (line.Jlow + Mlow) * (line.Jlow + Mlow - 1.0)
			       / (line.Jlow * (2.0 * line.Jlow - 1.0) * (2.0 * line.Jlow + 1.0));
	  
	}
	
	
	/* --- Zeeman splitting and strength ---*/
	
	double splitting = line.Gup*Mup - line.Glow*Mlow;
	line.splitting.push_back(splitting);
	line.strength.push_back(strength);
	line.iL.push_back(ilow-1);
	line.nZ++;
      }
    }
  }
}

//----------------------------------------------------------------
// Plank function at a given frequency fora given temperature
//----------------------------------------------------------------
double cprofiles::plank_nu(const double nu, const double temp){
  
  double c1 = (2.0 * HH * nu*nu*nu) / (CC*CC) ;
  double x = HH * nu / (BK * temp);

  if(x < 80.0) return c1 / (exp(x) - 1.0);
  else return         c1  * exp(-x);
}

//----------------------------------------------------------------
// Get Doppler factor, without the freq
//----------------------------------------------------------------
double cprofiles::get_doppler_factor(double temp, double vturb, double amass){
  return sqrt(2.0 * BK * temp / (amass * AMU)  + vturb * vturb) / CC;
}

//----------------------------------------------------------------
// Damping, assuming VALD3 input
//----------------------------------------------------------------
double cprofiles::damp(line_t &line, double temp, double vturb, double nne, double nh, double nhe, double dlnu){

   /* --- radiative damping --- */
  double adamp = line.g_rad; // in input.cc, if g_rad is 0, then the value is filled with an approx. formula.
  adamp +=  vanderWaals(line, temp, nh, nhe);  
  adamp += stark(line.g_str, temp, nne);
  adamp /= 4.0 * PI* dlnu;

  return adamp;
}

//----------------------------------------------------------------
// Stark effect (collisions with charged particles),
// Assuming VALD input, we get gamma_4 * nne (at 10000 K) to correct
// for temperature, use x(temp/10000.)^(1./6.)
//----------------------------------------------------------------
double cprofiles::stark(double g_str, double temp, double nne){
  return g_str * nne * pow(temp/10000.,1.0/6.0); // Assuming Vald input!
}

//----------------------------------------------------------------
// van der Waals broadening (collisions with neutral particles, typically H)
// Includes Barklem formulation
//----------------------------------------------------------------
double cprofiles::vanderWaals(line_t &line,  double temp, double nh,  double nhe){
  if(line.g_vdw >= 20.0){
    
    /* --- Barklem (constants initialized in input.cc) --- */
    //     http://www.astro.uu.se/~barklem/howto.html //
    
    if(line.firsttime == 1){
      line.b_sig *= 2.80028e-21;
      double gx = (2.0 - line.b_alp*0.5) - 1.0;
      
      double gammaf =
	1.0 + (-0.5748646 + (0.9512363 + (-0.6998588 + (0.4245549 - 0.1010678 * gx) * gx) * gx) * gx) * gx;
      
      line.b_gvw = pow((4.0 / PI), (line.b_alp * 0.5)) * gammaf * line.b_sig * 1.0e4;
      line.b_vbar =  21172.6 * ( 1.0 / 1.008 + 1.0 / line.amass);
      line.firsttime = 0;
    }
    
    /* --- Actual calculation of vdW broadening --- */
    
    double vbar = sqrt(temp * line.b_vbar);
    return pow((vbar / 1.E4), (1.0 - line.b_alp)) *  line.b_gvw * (nh + 0.42*nhe) * 2.e6;
    
  }else if(line.g_vdw != 0.0){
    
    /* --- Unsoeld 1955 from Nikolai Piskunov's routines, assuming VALD input --- */
    
    return pow(10.0,line.g_vdw) * pow(temp / 10000., 0.3) * (nh + 0.42*nhe);
    
  }else return 0.0; 
}

//----------------------------------------------------------------
// This function calculates the 3-j symbol
// J_i and M_i have to be twice the actual value of J and M
// Adapted from Polarization in spectral lines, Landi Degl'innocenti & Landolfi (2004) 
//----------------------------------------------------------------
double cprofiles::w3js(const int j1, const int j2, const int j3, const int m1, const int m2, const int m3){

  
  if ((m1+m2+m3) != 0) return 0.0;

  double w3js = 0.0;
  int ia = j1 + j2;
  if (j3 > ia) return 0.0;
  int ib = j1 - j2;
  if (j3 < abs(ib)) return 0.0;
  if (abs(m1) > j1) return 0.0;
  if (abs(m2) > j2) return 0.0;
  if (abs(m3) > j3) return 0.0;

  int jsum = j3 + ia;
  int ic = j1 - m1;
  int id = j2 - m2;

  int ie = j3 - j2 + m1;
  int im = j3 - j1 - m2;
  int zmin = max(max(0,-ie),-im);
  int ig = ia - j3;
  int ih = j2 + m2;
  int zmax = min(min(ig,ih),ic);
  double cc = 0.0;

  for(int z = zmin; z <= zmax; z += 2){
    double denon = fact[z/2] * fact[(ig-z)/2] * fact[(ic-z)/2] *
      fact[(ih-z)/2] * fact[(ie+z)/2] * fact[(im+z)/2];
    if((z % 4) != 0) denon = -denon;
    cc += 1.0 / denon;
  }

  double cc1 = fact[ig/2] * fact[(j3+ib)/2] * fact[(j3-ib)/2] / fact[(jsum+2)/2];
  double cc2 = fact[(j1+m1)/2] * fact[ic/2] * fact[ih/2] * fact[id/2] * fact[(j3-m3)/2]*fact[(j3+m3)/2];
  cc *= sqrt(cc1*cc2);

  if(((ib-m3) % 4) != 0) cc = -cc;
  w3js = cc;
  if(abs(w3js) < 1.e-8) w3js = 0.0;
  return w3js;
}

// -------------------------------------------------------------------------
// Compute the Voigt-Faraday function, for a given wavelength and height
// Adapted to C++ from A. Asensio-Ramos' Fortran Milne-Eddington routines 
// -------------------------------------------------------------------------
void cprofiles::zeeman_profile(double nu, const line_t &line, double vel, double bfield, double dlnu, double damping){
  
  /* --- Init vars --- */
  memset(&voigt[0],   0, sizeof(double)*3);
  memset(&faraday[0], 0, sizeof(double)*3);

  
  /* --- From Landi Degl'innocenti 2004, pag. 385, eq. 9.23-9.24 --- */
  
  double v  = (line.nu0 - nu)  / dlnu;
  double va = line.nu0 * vel   / (CC * dlnu);
  double vb = LARMOR * bfield  / dlnu;

  
  for(int ii=0;ii<line.nZ; ii++){      
    
    /* --- Voigt-Faraday profile ---*/
    
    double voi=0, far=0;
    voigt_complex(damping, v - va + vb * line.splitting[ii], voi, far);
    
    voigt[line.iL[ii]] += voi * line.strength[ii];
    faraday[line.iL[ii]] += 2.0 * far * line.strength[ii]; // L = 2 * F
    
  }
  
  /* --- Normalize the profile: prof / (sqrt(pi) * dldop)
     In Gray 2005, the sqrt(pi) is compensated with the pi 
     factor in the absorption coeff., but we keep it for readability
     --- */
  dlnu = 1.0 / (dlnu*sqrt(PI));
  //
  for(int ii = 0; ii<3; ii++){
    voigt[ii] *= dlnu;
    faraday[ii] *= dlnu;
  }
  //
}

// -------------------------------------------------------------------------
// Compute the terms of the absorption matrix for a given height
// -------------------------------------------------------------------------
void cprofiles::zeeman_opacity(double inc, double azi,  double lineop, int idep, int wav){

  /* --- Compute sin/cos of the inclination and azimuth --- */
  azi *= 2.0;//* DTOR;
  //
  double sinin = sin(inc);
  double cosin = cos(inc);
  double sinin2 = sinin * sinin;
  double cosin2 = cosin * cosin;
  double sin2az = sin(azi);
  double cos2az = cos(azi);
  lineop *= 0.5; // Includes the 1/2 multiplying each term of the equation
  
  /* --- Compute elements of the ABS matrix (Landi Degl'innocenti 2004, eq. 9.32)--- */
  mki(idep, wav) += lineop * (voigt[1] * sinin2 + 0.5 * (voigt[0] + voigt[2]) * (1.0 + cosin2));

  
  double tmp =  lineop * (  voigt[1] - 0.5 * (  voigt[0] +   voigt[2])) * sinin2;
  double tmp1 = lineop * (faraday[1] - 0.5 * (faraday[0] + faraday[2])) * sinin2;

  mkq(idep, wav) += tmp  * cos2az;
  mfq(idep, wav) += tmp1 * cos2az;
  
  mku(idep, wav) += tmp  * sin2az;
  mfu(idep, wav) += tmp1 * sin2az;
  
  mkv(idep, wav) += lineop * (voigt[2]   -   voigt[0]) * cosin;
  mfv(idep, wav) += lineop * (faraday[2] - faraday[0]) * cosin;
}

//-------------------------------------------------------------------------
// Compute the effective Lande factor with jlow, jup, glow, gup,
// Eq. 3.44 - Landi Degl'Innocenti & Landolfi (2004)
//-------------------------------------------------------------------------
double cprofiles::lande_factor(double j1, double j2, double g1, double g2){
  return 0.5 * (g1 + g2) + 0.25 * (g1 - g2) * (j1 * (j1 + 1) - j2 * (j2 + 1));
}


//-------------------------------------------------------------------------
// Delo-Lin formal solver, using z as input
//-------------------------------------------------------------------------
void cprofiles::delolin(int ndep, double *z, double *stokes, double mu){
  
  //  int ndep = (int)z.size();

  
  /* --- Init arrays --- */
  double *dtau = new double [ndep];
  double stk[4];
  memset(&dtau[0], 0, sizeof(double) * ndep);
  memset(&stk[0],  0, sizeof(double) * 4);

  
  int k0 = 0, k1 = ndep-1;
  double itau = 0.0;

  /* --- Compute dtau_nu scale using linear approx --- */
  double imu = fabs(1.0 / mu);
  for(int k = 1; k<ndep; k++){
    dtau[k] = 0.5 * (ki[k-1] + ki[k]) * fabs(z[k-1] - z[k]) * imu;
    itau += dtau[k];
    
    /* --- get integration limits --- */
    if(itau <= 1.E-4) k0 = k;
    if(itau <= 100) k1 = k;
  }

  /* --- Init integration at the lower boundary I = SF --- */
  stk[0] = sf[k1];
  double Ku[4][4], K0[4][4], Su[4], S0[4];

  
  /* --- Init source vector & Abs. matrix at upwind point --- */
  abmat(k1, Ku);
  source_vector(k1, Su);
  
  for(int k = k1-1; k >= k0; k--){

    int k1 = k + 1;

    /* --- Init source vector and Abs. matrix for depth "k" --- */
    abmat(k, K0);
    source_vector(k, S0);

    
    /* --- dtau, exponentials and integration coeffs. --- */
    double dt = dtau[k1];
    double eps = 0.0;
    double cu = 0, c0 = 0;
    //
    if(dt >= 1.e-5){

      eps = exp(-dt);
      double u0 = 1.0 - eps;
      double u1 = dt - 1.0 + eps;
      
      c0 = u1 / dt;
      cu = u0 - c0;
      
    }else{
      eps = 1.0 - dt + 0.5 * dt*dt - dt*dt*dt/6.0;

      double dt2 = dt * dt;
      c0 = dt * 0.5 - dt2 / 6.0;
      cu = dt * 0.5 - dt2 / 3.0;
      
    }

    /* --- Compute terms for the integration --- */
    double mat1[4][4];
    double mat2[4][4];
    double vec1[4] = {};

    //
    for(int j = 0; j<4; j++){
      for(int i = 0; i<4; i++){
	mat1[j][i] = ident[j][i] * eps - cu * Ku[j][i];
	mat2[j][i] = ident[j][i]       + c0 * K0[j][i];
      }
    }
    
    m4inv(mat2); // Invert matrix
    m4v(mat1, stk, vec1); // Matrix x vector

    
    for(int i = 0; i<4; i++){
      vec1[i] += cu * Su[i] + c0 * S0[i];
    }
    
    m4v(mat2, vec1, stk); // Matrix x vector
    
    /* --- Copy variables to upwind arrays for next height ---*/
    memcpy(&Su[0],    &S0[0],     4*sizeof(double));
    memcpy(&Ku[0][0], &K0[0][0], 16*sizeof(double));
  }
  
  memcpy(&stokes[0], &stk[0], 4*sizeof(double));
  delete [] dtau;
}

//-------------------------------------------------------------------------
// Delo-Bezier formal solver, using z as input
// REFERENCE: de la Cruz Rodriguez & Piskunov (2013)
// Coded by J. de la Cruz Rodriguez
//
// NOTE: "u" stands for upwind and "d" for downwind point,
//       relative to the direction of the integration.
//
// NOTE2: The fast index is assumed to loop columns in a 4x4 matrix.
//        Unlike fortran and IDL, this is the rightmost index.
//
// NOTE3: It assumes that kq,ku,kv,fq,fu,fv are normalized by ki.
//-------------------------------------------------------------------------
void cprofiles::delobez3(int ndep, double *z, double *stokes, double mu){


  
  /* --- Init arrays --- */
  double *dtau = new double [ndep];
  double stk[4];
  memset(&dtau[0], 0, sizeof(double) * ndep);
  memset(&stk[0],  0, sizeof(double) * 4);

  
  int k0 = 0, k1 = ndep-1;
  double itau = 0.0;

  /* --- Compute dtau_nu scale using a bezier interpolant approx --- */
  //
  double dzu = z[1]  - z[0];
  double deu  = (ki[1] - ki[0]) / dzu;
  double odki = deu;
  double dki, dzd, ded;
  //
  double imu = fabs(1.0 / mu);
  for(int k = 1; k<ndep-1; k++){

    
    int kd = k + 1; // index of the downwind point
    
    dzd = z[kd] - z[k];
    ded = (ki[kd] - ki[k]) / dzd;

    
    /* --- Derivative of the opacity following Fritsch & Butland (1984) --- */
    if(deu*ded > 0.0){
      double lambda = (1.0 + dzd / (dzd + dzu)) / 3.0;
      dki = (deu / (lambda * ded + (1.0 - lambda) * deu)) * ded;
    } else dki = 0.0;

    
    /* --- integrate opacity using cubic bezier splines
       The Bezier3 integral should be:
       dz * (ki_0 + ki_u + cntrl1 + cntr2) / 4 
       --- */
    
    dtau[k] = fabs(dzu) * ((ki[k] - dki/3.0 * dzu) + (ki[k-1] + odki/3.0 * dzu) + ki[k] + ki[k-1]) * 0.25 * imu;
    itau += dtau[k];
    
    /* --- Store values for next interval --- */
    
    dzu  = dzd;
    deu  = ded;
    odki = dki;

    
    /* --- Integration limits --- */

    if(itau <= 1.E-4) k0 = k;
    if(itau <= 100.0) k1 = k;
    else break; // If tau is already larger than 100, then stop integrating the opacity
  }
  
  /* --- 
     Reached the boundary? Use quadratic Bezier, and compute 
     eC with the derivative from the upwind point 
     --- */
  if(k1 == ndep-2) {
    dtau[ndep-1] = fabs(dzu) * ( (ki[ndep-2] + odki/3.0 * dzu) + ki[ndep-1] + ki[ndep-2]) / 3.0 * imu;
    itau += dtau[ndep-1];
    if(itau <= 1.5E1) k1 = ndep-1;
  }
   
  /* --- Define some vars and at the lower boundary make I = SF --- */
  
  stk[0] = sf[k1];
  //
  double Ku[4][4], K0[4][4], Su[4], S0[4], tmpa[4][4], tmpb[4][4], A[4][4], tmpc[4][4]; 
  double dkq[ndep], dku[ndep], dkv[ndep], dfq[ndep], dfu[ndep], dfv[ndep],
    dSv[ndep][4], vtemp[4][ndep], vtemp1[4][ndep];
  double dK_u[4][4], dK_0[4][4], v0[4];

  
  /* -- centered derivatives of all independent terms in the Abs matrix --- */
  
  cent_deriv(ndep, dtau, &kq[0], &dkq[0], k0, k1);
  cent_deriv(ndep, dtau, &ku[0], &dku[0], k0, k1);
  cent_deriv(ndep, dtau, &kv[0], &dkv[0], k0, k1);
  cent_deriv(ndep, dtau, &fq[0], &dfq[0], k0, k1);
  cent_deriv(ndep, dtau, &fu[0], &dfu[0], k0, k1);
  cent_deriv(ndep, dtau, &fv[0], &dfv[0], k0, k1);

  /* --- Derivative of the elements of the source vector --- */
  
  for(int k = max(0,k0-1); k<min(ndep,k1+2); k++){
    vtemp[0][k] = sf[k];
    vtemp[1][k] = sf[k] * kq[k];
    vtemp[2][k] = sf[k] * ku[k];
    vtemp[3][k] = sf[k] * kv[k];
  }
  
  cent_deriv(ndep, dtau, &vtemp[0][0], &vtemp1[0][0], k0, k1);
  cent_deriv(ndep, dtau, &vtemp[1][0], &vtemp1[1][0], k0, k1);
  cent_deriv(ndep, dtau, &vtemp[2][0], &vtemp1[2][0], k0, k1);
  cent_deriv(ndep, dtau, &vtemp[3][0], &vtemp1[3][0], k0, k1);


  for(int k = k0; k <= k1; k++){ // Transpose Source vector array
    dSv[k][0] = vtemp1[0][k];
    dSv[k][1] = vtemp1[1][k];
    dSv[k][2] = vtemp1[2][k];
    dSv[k][3] = vtemp1[3][k];
  }
  
  /* --- Init source vector, K and dK/dtau at upwind point --- */
  abmat(k1, Ku);
  source_vector(k1, Su);
  abmat(dkq[k1], dku[k1], dkv[k1], dfq[k1], dfu[k1], dfv[k1], dK_u);
  
  
  /* --- Integrate over all relevant depths --- */
  for(int k = k1-1; k >= k0; k--){
    
    int ku = k + 1; // index of the upwind point

    /* --- Init source vector and Abs. matrix for depth "k" --- */
    abmat(k, K0);
    source_vector(k, S0);
    abmat(dkq[k], dku[k], dkv[k], dfq[k], dfu[k], dfv[k], dK_0);


    
    /* --- Integration coeffs. and exponential --- */
    double dt = dtau[ku];
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt03 = dt / 3.0;
    double eps, alp, bet, gam, mu;
    //
    if(dt >= 1.e-2){
      //
      eps = exp(-dt);
      //
      alp = (-6.0 + 6.0 * dt - 3.0 * dt2 + dt3 + 6.0 * eps)        / dt3;
      bet = (6.0 + (-6.0 - dt * (6.0 + dt * (3.0 + dt))) * eps)    / dt3;
      gam = 3.0 * (6.0 + (-4.0 + dt)*dt - 2.0 * (3.0 + dt) * eps)  / dt3;
      mu  = 3.0 * ( eps * (6.0 + dt2 + 4.0 * dt) + 2.0 * dt - 6.0) / dt3;
    }else{ // Taylor expansion of the exponential
      //
      double dt4 = dt2 * dt2;
      eps = 1.0 - dt + 0.5 * dt2 - dt3 / 6.0 + dt4 / 24.0;
      //
      alp = 0.25 * dt - 0.05 * dt2 + dt3 / 120.0 - dt4 / 840.0;
      bet = 0.25 * dt - 0.20 * dt2 + dt3 / 12.0  - dt4 / 42.0; 
      gam = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025 - dt4 / 210.0; 
      mu  = 0.25 * dt - 0.15 * dt2 + dt3 * 0.05  - dt4 / 84.0; 
    }    
    
    m4m(Ku, Ku, tmpa); //  Ku^2
    m4m(K0, K0, A);    //  K0^2

    /* --- Compute temporary matrices --- */
    for(int j = 0; j<4; j++){
      for(int i = 0; i<4; i++){
	A[j][i] = ident[j][i] + alp * K0[j][i] - gam * -(dt03 * (A[j][i] + dK_0[j][i] + K0[j][i]) + K0[j][i]);
	tmpa[j][i] = eps * ident[j][i] - bet * Ku[j][i] + mu * (dt03 * (tmpa[j][i] + dK_u[j][i] + Ku[j][i]) - Ku[j][i]);
	tmpb[j][i] = bet * ident[j][i] + mu * (ident[j][i] - dt03*Ku[j][i]);
	tmpc[j][i] = alp * ident[j][i] + gam* (ident[j][i] + dt03*K0[j][i]);
      }
    }

    // Here I am doing tmpa*stk + tmpb * Su + tmpc * S0 + (gam * dS0 - mu * dSu) * dtau / 3.0
    for(int i = 0; i<4; i++){
      v0[i] = 0.0;
      for(int j = 0; j<4; j++){
	v0[i] += tmpa[i][j] * stk[j] + tmpb[i][j] * Su[j] + tmpc[i][j] * S0[j];
      }
      v0[i] += dt03 * (gam * dSv[k][i] - mu * dSv[ku][i]);
    }
    m4inv(A); // Invert "A"

    
    /* --- Get new intensity at depth k --- */
    m4v(A, v0, stk); 

    
    /* --- Copy variables for next interval --- */
    memcpy(&Su[0], &S0[0],      4 * sizeof(double)); // central point -> upwind
    memcpy(&Ku[0][0], &K0[0][0],     16 * sizeof(double)); // central point -> upwind
    memcpy(&dK_u[0][0], &dK_0[0][0], 16 * sizeof(double)); // central point -> upwind
  }


  //  fprintf(stdout, "%E %E %E %E \n",stk[0],stk[1],stk[2],stk[3]);
  memcpy(&stokes[0], &stk[0], 4*sizeof(double));
  delete [] dtau;
}

void cprofiles::cent_deriv(int n, const double *dx, const double *y, double *yp, int k0, int k1){
  // Assumes that yp has been allocated: yp[n]
  int kinit = max(1,k0);
  int kend = min(k1, n-2);
  
  double oder = (y[kinit] - y[kinit-1]) / dx[kinit];
  if(k0 == 0) yp[0] = oder;
  
  for(int k = kinit; k <= kend; k++){
    double der = (y[k+1] - y[k]) / dx[k+1];
    
    if(der*oder > 0.0){
      double lambda = (1.0 + dx[k+1] / (dx[k+1] + dx[k])) / 3.0;
      yp[k] = (oder / (lambda * der + (1.0 - lambda) * oder)) * der;
    } else yp[k] = 0.0; // Set der to zero at extrema;
    
    oder = der;
  }
  if(k1 == (n-1)) yp[n-1] = oder;			  
}
