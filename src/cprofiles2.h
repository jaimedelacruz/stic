#ifndef CPROF2_H
#define CPROF2_H
//
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
#include <cmath>
#include "physical_consts.h"
#include "cmemt.h"
#include "input.h"
//
enum cprof_solver{
  bez_ltau,
  lin_ltau,
  bez_z,
  lin_z,
};
//
class cprofiles{
 public:
  static constexpr double A[7] = {122.607931777104326, 214.382388694706425, 181.928533092181549,
				  93.155580458138441, 30.180142196210589, 5.912626209773153,
				  0.564189583562615};
  static constexpr double B[8] = {122.60793177387535, 352.730625110963558, 457.334478783897737, 
				  348.703917719495792, 170.354001821091472, 53.992906912940207, 
				  10.479857114260399, 1.0};
   
  static constexpr double LARMOR= phyc::EE / (4.0 * phyc::PI * phyc::ME * phyc::CC);

  int nndep, nnw;
  std::vector<double> ki, kq, ku, kv, fq, fu, fv, sf;
  double **mki, **mkq, **mku, **mkv, **mfq, **mfu, **mfv;
  
  double voigt[3], faraday[3];
  typedef double mat4[4][4]; // define a type to pass a matrix as an argument and keep the shape
  typedef double vect4[4];
  

  /* --- methods ---*/
  
  cprofiles(){};

  /* 
     --- Voigt-Faraday profiles implemented here to make inlining easier for the compiler --- 
  */
  inline void voigt_complex(const double a, const double v, double &vgt, double &far){
    double sav = fabs(v) + a;
    std::complex<double> tav(a, -v);
    std::complex<double> uav = tav*tav;
    std::complex<double> w4;
  
    /* --- HUMLICEK'S APPROXIMATION --- */
    if(sav >= 15.0){
      w4 = tav * 0.5641896 / (0.5 + uav);
    } else if(sav >= 5.5){
      w4 = tav * (1.410474 + uav * 0.5641896) / (0.75+uav * (3.0 + uav));
    } else if(a >= (0.195 * fabs(v) - 0.176)){
      w4 = (16.4955 + tav * (20.20933 + tav * (11.96482 + tav * (3.778987 +
								 tav * 0.5642236)))) / (16.4955 + tav * (38.82363 + tav * (39.27121 +
															   tav * (21.69274 + tav * (6.699398 + tav)))));
    } else{
      w4 = tav * (36183.31 - uav * (3321.9905 - uav * (1540.787 -
						       uav * (219.0313 - uav * (35.76683 - uav * (1.320522 -
												  uav * 0.56419))))));
      
      std::complex<double> v4 = (32066.6 - uav * (24322.84 - uav * (9022.228 - 
								    uav * (2186.181 - uav * (364.2191 - uav * (61.57037 -
													       uav * (1.841439 - uav)))))));
      w4 = exp(uav) - w4 / v4;
    }
    
    /* ---  Note that FVGT below is in fact 2 * (Faradey-Voigt function) ---*/
    vgt = w4.real();
    far = 0.5 * w4.imag();
  }
  
  inline void voigtf(double damp, double vv, double &H, double &F){
    
    std::complex<double> Z(damp,-fabs(vv));
    
    Z = ((((((A[6]*Z+A[5])*Z+A[4])*Z+A[3])*Z+A[2])*Z+A[1])*Z+A[0]) /
      (((((((Z+B[6])*Z+B[5])*Z+B[4])*Z+B[3])*Z+B[2])*Z+B[1])*Z+B[0]);
    
    
    H = Z.real();
    F = ((vv<0.0)?-0.5:0.5) * Z.imag();
  
  }

  
  /* 
     --- Construct abs matrix ---
     JdlCR: mat4 is a double[4][4] defined above
     needed this format to pass a reference 
     to a matrix and keep the variable shape 
     within the function.
  */
  inline void abmat(int idep, mat4 &k){
    k[0][0] = 0.0;//ki[idep];
    k[0][1] = kq[idep];///ki[idep];
    k[0][2] = ku[idep];///ki[idep];
    k[0][3] = kv[idep];///ki[idep];
    //
    k[1][0] = kq[idep];///ki[idep];
    k[1][1] = 0.0;//ki[idep];
    k[1][2] = fv[idep];///ki[idep];
    k[1][3] =-fu[idep];///ki[idep];
    //
    k[2][0] = ku[idep];///ki[idep];
    k[2][1] =-fv[idep];///ki[idep];
    k[2][2] = 0.0;//ki[idep];
    k[2][3] = fq[idep];///ki[idep];
    //
    k[3][0] = kv[idep];///ki[idep];
    k[3][1] = fu[idep];///ki[idep];
    k[3][2] =-fq[idep];///ki[idep];
    k[3][3] = 0.0;//ki[idep];
  }

  void abmat(double ikq, double iku, double ikv, double ifq, double ifu, double ifv, mat4 &k){
    k[0][0] = 0.0;//ki[idep];
    k[0][1] = ikq;//[idep];///ki[idep];
    k[0][2] = iku;//[idep];///ki[idep];
    k[0][3] = ikv;//[idep];///ki[idep];
    //
    k[1][0] = ikq;//[idep];///ki[idep];
    k[1][1] = 0.0;//ki[idep];
    k[1][2] = ifv;//[idep];///ki[idep];
    k[1][3] =-ifu;//[idep];///ki[idep];
    //
    k[2][0] = iku;//[idep];///ki[idep];
    k[2][1] =-ifv;//[idep];///ki[idep];
    k[2][2] = 0.0;//ki[idep];
    k[2][3] = ifq;//[idep];///ki[idep];
    //
    k[3][0] = ikv;//[idep];///ki[idep];
    k[3][1] = ifu;//[idep];///ki[idep];
    k[3][2] =-ifq;//[idep];///ki[idep];
    k[3][3] = 0.0;//ki[idep];
  }
  

  inline void source_vector(int idep, vect4 &S){
    S[0] =            sf[idep];
    S[1] = kq[idep] * sf[idep]; // / ki[idep];
    S[2] = ku[idep] * sf[idep]; // / ki[idep];
    S[3] = kv[idep] * sf[idep]; // / ki[idep];
  }
  
  /* --- matrix multiplication --- */
  void m4m(mat4 &a, mat4 &b, mat4 &c){
    memset(&c[0][0],0,sizeof(double)*16);
    
    for(int j = 0; j<4; j++)
      for(int i = 0; i<4; i++)
	for(int k = 0; k<4; k++)
	  c[j][i] += a[k][i]*b[j][k]; 
  }

  /* --- matrix/vector multiplication --- */
  void m4v(mat4 &a, vect4 &b, vect4 &c){
    memset(&c[0],0,sizeof(double)*4);
    
    for(int i = 0; i<4; i++){
      for(int k = 0; k<4; k++){
	c[i] += a[i][k] * b[k];
      }
    }
    
  }

  void m4inv(mat4 &MI){
    for(int k=0;k<4;k++){
      MI[k][k]=-1.0/MI[k][k];         // the pivot element 
      for(int i=0;i<4;++i) if(i!=k) MI[i][k]*=MI[k][k];//the pivot column 
      for(int i=0;i<4;++i)           //elements not in a pivot row or column
	if(i!=k)
	  for(int j=0;j<4;++j)
	    if(j!=k)
	      MI[i][j]+=MI[i][k]*MI[k][j];
      for(int i=0;i<4;++i)           //elements in a pivot row
	if(i!=k)
	  MI[k][i]*=MI[k][k];
    }
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j) MI[i][j]=-MI[i][j];
    return;
  }

  //----------------------------------------------------------------
  // Constructor
  //----------------------------------------------------------------
  cprofiles(int nw, int ndep){
    nnw = 0;
    nndep = 0;
    init(nw, ndep);
  }


  void init(int nw, int ndep){

    /* --- do not reallocate if already allocated --- */
    
    //if(nnw == nw && nndep == ndep)
    //  return;

    /* --- Otherwise allocate everything --- */
    
    nndep = ndep;
    nnw = nw;
  
    mki = (double**)mat2d(ndep, nw);
    mkq = (double**)mat2d(ndep, nw);
    mku = (double**)mat2d(ndep, nw);
    mkv = (double**)mat2d(ndep, nw);
    mfq = (double**)mat2d(ndep, nw);
    mfu = (double**)mat2d(ndep, nw);
    mfv = (double**)mat2d(ndep, nw);

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
  ~cprofiles(){
    //cleanup();
  }

  //----------------------------------------------------------------
  // Cleanup
  //----------------------------------------------------------------
  void cleanup(){
    ki.clear();
    kq.clear();
    ku.clear();
    kv.clear();
    fq.clear();
    fu.clear();
    fv.clear();
  
    if(mki != NULL) del_mat(mki);
    if(mkq != NULL) del_mat(mkq);
    if(mku != NULL) del_mat(mku);
    if(mkv != NULL) del_mat(mkv);
    if(mfq != NULL) del_mat(mfq);
    if(mfu != NULL) del_mat(mfu);
    if(mfv != NULL) del_mat(mfv);
  }


  //-------------------------------------------------------------------------
  // Set to zero the terms of the Abs. Matrix at all depths
  //-------------------------------------------------------------------------
  void set_zero(){ 
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
  void set_zero_abmat(){
    memset(&mki[0][0], 0, nndep*nnw*sizeof(double));
    memset(&mkq[0][0], 0, nndep*nnw*sizeof(double));
    memset(&mku[0][0], 0, nndep*nnw*sizeof(double));
    memset(&mkv[0][0], 0, nndep*nnw*sizeof(double));
    memset(&mfq[0][0], 0, nndep*nnw*sizeof(double));
    memset(&mfu[0][0], 0, nndep*nnw*sizeof(double));
    memset(&mfv[0][0], 0, nndep*nnw*sizeof(double));
  }

  //-------------------------------------------------------------------------
  // Compute Zeeman splitting
  //-------------------------------------------------------------------------
  void init_zeeman_components(line_t &line)
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
	if(fabs(Mlow) <= line.Jlow){
	
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
  double plank_nu(const double nu, const double temp){
  
    double c1 = (2.0 * phyc::HH * nu*nu*nu) / (phyc::CC*phyc::CC) ;
    double x = phyc::HH * nu / (phyc::BK * temp);

    if(x < 80.0) return c1 / (exp(x) - 1.0);
    else return         c1  * exp(-x);
  }

  //----------------------------------------------------------------
  // Get Doppler factor, without the freq
  //----------------------------------------------------------------
  double get_doppler_factor(double temp, double vturb, double amass){
    return sqrt(2.0 * phyc::BK * temp / (amass * phyc::AMU)  + vturb * vturb) / phyc::CC;
  }

  //----------------------------------------------------------------
  // Damping, assuming VALD3 input
  //----------------------------------------------------------------
  double damp(line_t &line, double temp, double vturb, double nne, double nh, double nhe, double dlnu){

    /* --- radiative damping --- */
    double adamp = line.g_rad; // in input.cc, if g_rad is 0, then the value is filled with an approx. formula.
    adamp +=  vanderWaals(line, temp, nh, nhe);  
    adamp += stark(line.g_str, temp, nne);
    adamp /= 4.0 * phyc::PI* dlnu;

    return adamp;
  }

  //----------------------------------------------------------------
  // Stark effect (collisions with charged particles),
  // Assuming VALD input, we get gamma_4 * nne (at 10000 K) to correct
  // for temperature, use x(temp/10000.)^(1./6.)
  //----------------------------------------------------------------
  double stark(double g_str, double temp, double nne){
    return g_str * nne * pow(temp/10000.,1.0/6.0); // Assuming Vald input!
  }

  //----------------------------------------------------------------
  // van der Waals broadening (collisions with neutral particles, typically H)
  // Includes Barklem formulation
  //----------------------------------------------------------------
  double vanderWaals(line_t &line,  double temp, double nh,  double nhe){
    if(line.g_vdw >= 20.0){
    
      /* --- Barklem (constants initialized in input.cc) --- */
      //     http://www.astro.uu.se/~barklem/howto.html //
    
      if(line.firsttime == 1){
	line.b_sig *= 2.80028e-21;
	double gx = (2.0 - line.b_alp*0.5) - 1.0;
      
	double gammaf =
	  1.0 + (-0.5748646 + (0.9512363 + (-0.6998588 + (0.4245549 - 0.1010678 * gx) * gx) * gx) * gx) * gx;
      
	line.b_gvw = pow((4.0 / phyc::PI), (line.b_alp * 0.5)) * gammaf * line.b_sig * 1.0e4;
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
  double w3js(const int j1, const int j2, const int j3, const int m1, const int m2, const int m3){

  
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
    int zmin = std::max(std::max(0,-ie),-im);
    int ig = ia - j3;
    int ih = j2 + m2;
    int zmax = std::min(std::min(ig,ih),ic);
    double cc = 0.0;

    for(int z = zmin; z <= zmax; z += 2){
      double denon = phyc::fact[z/2] * phyc::fact[(ig-z)/2] * phyc::fact[(ic-z)/2] *
	phyc::fact[(ih-z)/2] * phyc::fact[(ie+z)/2] * phyc::fact[(im+z)/2];
      if((z % 4) != 0) denon = -denon;
      cc += 1.0 / denon;
    }

    double cc1 = phyc::fact[ig/2] * phyc::fact[(j3+ib)/2] * phyc::fact[(j3-ib)/2] / phyc::fact[(jsum+2)/2];
    double cc2 = phyc::fact[(j1+m1)/2] * phyc::fact[ic/2] * phyc::fact[ih/2] * phyc::fact[id/2] * phyc::fact[(j3-m3)/2]*phyc::fact[(j3+m3)/2];
    cc *= sqrt(cc1*cc2);

    if(((ib-m3) % 4) != 0) cc = -cc;
    w3js = cc;
    if(fabs(w3js) < 1.e-8) w3js = 0.0;
    return w3js;
  }

  // -------------------------------------------------------------------------
  // Compute the Voigt-Faraday function, for a given wavelength and height
  // Adapted to C++ from A. Asensio-Ramos' Fortran Milne-Eddington routines 
  // -------------------------------------------------------------------------
  void zeeman_profile(double nu, const line_t &line, double vel, double bfield, double dlnu, double damping){
  
    /* --- Init vars --- */
    memset(&voigt[0],   0, sizeof(double)*3);
    memset(&faraday[0], 0, sizeof(double)*3);

  
    /* --- From Landi Degl'innocenti 2004, pag. 385, eq. 9.23-9.24 --- */
  
    double v  = (line.nu0 - nu)  / dlnu;
    double va = line.nu0 * vel   / (phyc::CC * dlnu);
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
    dlnu = 1.0 / (dlnu*phyc::SQPI);
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
  void zeeman_opacity(double inc, double azi,  double lineop, int idep, int wav){

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
    mki[idep][wav] += lineop * (voigt[1] * sinin2 + 0.5 * (voigt[0] + voigt[2]) * (1.0 + cosin2));

  
    double tmp =  lineop * (  voigt[1] - 0.5 * (  voigt[0] +   voigt[2])) * sinin2;
    double tmp1 = lineop * (faraday[1] - 0.5 * (faraday[0] + faraday[2])) * sinin2;

    mkq[idep][wav] += tmp  * cos2az;
    mfq[idep][wav] += tmp1 * cos2az;
  
    mku[idep][wav] += tmp  * sin2az;
    mfu[idep][wav] += tmp1 * sin2az;
  
    mkv[idep][wav] += lineop * (voigt[2]   -   voigt[0]) * cosin;
    mfv[idep][wav] += lineop * (faraday[2] - faraday[0]) * cosin;
  }

  //-------------------------------------------------------------------------
  // Compute the effective Lande factor with jlow, jup, glow, gup,
  // Eq. 3.44 - Landi Degl'Innocenti & Landolfi (2004)
  //-------------------------------------------------------------------------
  double lande_factor(double j1, double j2, double g1, double g2){
    return 0.5 * (g1 + g2) + 0.25 * (g1 - g2) * (j1 * (j1 + 1) - j2 * (j2 + 1));
  }


  //-------------------------------------------------------------------------
  // Delo-Lin formal solver, using z as input
  //-------------------------------------------------------------------------
  void delolin(int ndep, double *z, double *stokes, double mu){
  
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
	  mat1[j][i] = phyc::ident[j][i] * eps - cu * Ku[j][i];
	  mat2[j][i] = phyc::ident[j][i]       + c0 * K0[j][i];
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
  void delobez3(int ndep, double *z, double *stokes, double mu){


  
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
  
    for(int k = std::max(0,k0-1); k<std::min(ndep,k1+2); k++){
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
	  A[j][i] = phyc::ident[j][i] + alp * K0[j][i] - gam * -(dt03 * (A[j][i] + dK_0[j][i] + K0[j][i]) + K0[j][i]);
	  tmpa[j][i] = eps * phyc::ident[j][i] - bet * Ku[j][i] + mu * (dt03 * (tmpa[j][i] + dK_u[j][i] + Ku[j][i]) - Ku[j][i]);
	  tmpb[j][i] = bet * phyc::ident[j][i] + mu * (phyc::ident[j][i] - dt03*Ku[j][i]);
	  tmpc[j][i] = alp * phyc::ident[j][i] + gam* (phyc::ident[j][i] + dt03*K0[j][i]);
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


    memcpy(&stokes[0], &stk[0], 4*sizeof(double));
    delete [] dtau;
  }

  void cent_deriv(int n, const double *dx, const double *y, double *yp, int k0, int k1){
    // Assumes that yp has been allocated: yp[n]
    int kinit = std::max(1,k0);
    int kend = std::min(k1, n-2);
  
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

  double **mat2d(int nx1, int nx2){
    double **p;
    p = new double* [nx1];
    p[0] = new double [nx1 * nx2]();
    for(int x1=1;x1<nx1;++x1) p[x1] = p[x1-1] + nx2;
    return p;
  }

  void del_mat(double **p){
    delete[] (p[0]);
    delete[] (p);
    p = NULL;
  }
};
#endif
