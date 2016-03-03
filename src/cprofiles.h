#ifndef CPROF_H
#define CPROF_H
//
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include "depthmodel.h"
#include "input.h"
#include "ceos.h"
#include "physical_consts.h"
#include "cmemt.h"
//
enum cprof_solver{
  bez_ltau,
  lin_ltau,
  bez_z,
  lin_z,
};
//
class cprofiles{
  
 private:

 public:
  int nndep, nnw;
  std::vector<double> ki, kq, ku, kv, fq, fu, fv, sf;
  mat<double> mki, mkq, mku, mkv, mfq, mfu, mfv;
  
  double voigt[3], faraday[3];
  static const double fact[31];
  static const double LARMOR;
  static const double ident[4][4];
  typedef double mat4[4][4]; // define a type to pass a matrix as an argument and keep the shape
  typedef double vect4[4];
  
/* --- For the calculation of the Voigt-Faraday function --- */

  static constexpr double A[7] = {122.607931777104326, 214.382388694706425, 181.928533092181549,
			       93.155580458138441, 30.180142196210589, 5.912626209773153,
			       0.564189583562615};
  static constexpr double B[8] = {122.60793177387535, 352.730625110963558, 457.334478783897737, 
			       348.703917719495792, 170.354001821091472, 53.992906912940207, 
			       10.479857114260399, 1.0};

  
  /* --- methods ---*/
  cprofiles(int nw, int ndep);
  cprofiles(){};
  ~cprofiles();

  void init(int nw, int ndep);
  double w3js( int j1, int j2, int j3, int m1, int m2, int m3);
  void zeeman_profile(double nu, const line_t &line, double vel,
		      double bfield, double doppler_factor, double damping);
  double damp(line_t &line, double temp, double vturb, double nne,  double nh, double nhe, double dlnu);
  double vanderWaals(line_t &line,  double temp, double nh, double nhe);
  double stark(double g_str, double temp, double nne);
  double plank_nu(const double nu, const double temp);
  void init_zeeman_components(line_t &line);
  void zeeman_opacity(double inc,double azi,  double lineop, int idep, int wav);
  double lande_factor(double jlow, double jup, double glow, double gup);
  void set_zero();
  void set_zero_abmat();
  void delolin(int ndep, double *z, double *stokes, double mu = 1.0);
  void delobez3(int ndep, double *z, double *stokes, double mu = 1.0);
  double get_doppler_factor(double temp, double vturb, double amass);
  void cent_deriv(int n, const double *dx, const double *y, double *yp, int k0, int k1);

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
  } else if(a >= (0.195 * abs(v) - 0.176)){
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

  

  
};
#endif
