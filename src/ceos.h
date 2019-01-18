/*
  CEOS class
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
  Dependencies: Fortran EOS (N. Piskunov, UU)
 */
#ifndef CEOS_H
#define CEOS_H
//
#include <iostream>
#include <vector>
#include <string>
#include "input.h"
#include "cmemt.h"
#include "eoswrap.h"

/* Interface for Fortran routines */
extern "C" {
  int eqcount_(const char [][3], const char [][8], int *, int &, int &, int &, int, int);

  int eqlist_(const float *, const char [][3], const char [][8], int *, int *, char *, int &,
	      int &, int &,  int&, int, int, int);

  void eqstat_(int &, float &, float &, float &, float *, const char [][3],
	       const float *, int &, int *, char *, float *, float *, float *,
	       float *, int &, int &, float &, float &, float &, int &, int, int);
  void eqstat_rho_(int &, float &, float &, float &, float *, const char [][3],
	       const float *, int &, int *, char *, float *, float *, float *,
		   float *, int &, int &, float &, float &, float &, int &, int, int);
  //     subroutine eqstat_rho(mode,temp,Pg,Pe,abund,elemen,amass,
  //   &                  ELESIZ,spindx,splist,xfract,poti,atwght,
  //   &                  nlines,nlist,xne,xna,rho,niter)
  //     subroutine     eqstat(mode,temp,Pg,Pe,abund,elemen,amass,
  //   &                  ELESIZ,spindx,splist,xfract,pfunc,poti,atwght,
  //   &                  nlines,nlist,xne,xna,rho,niter)
  
  void check_(const char [][3], char [][8], int *, int &, int &, int &);

  /*
  void contop_(float &T, float &TKEV, float &TK, float &HKT, float &TLOG,
	       float &XNA, float &XNE, double *WLGRID, double *OPACITY,
	       double *SCATTER,  
	       float &H1, float &H2, float &HMIN, float &HE1, float &HE2,
	       float &HE3, float &C1, float &AL1, float &SI1, float &SI2,
	       float &CA1, float &CA2, float &MG1, float &MG2, float &FE1,
	       float &N1, float &O1, int &nWLGRID, int &NLINES, int &NTOTALLIST);
  */
}



/* class definition */
class ceos: public eoswrap{
 private:
  
  
 public:
  // Atomic data (see implementation file)
  //static const int   MAX_ELEM = 99;

  static const float AMASS[MAX_ELEM];
  static const char  ELEMEN[MAX_ELEM][3];
  static const float ABUND_default[MAX_ELEM];
  //float ABUND[MAX_ELEM];
  int NELEM, NLIST, NLINES;

  
  // Physical constants
  static constexpr double bk = 1.3806488E-16;
  static constexpr double cc = 2.99792458E10;
  static constexpr double mp = 1.672621777E-24;
  static constexpr double me = 9.10938215E-28;

  // EOS book keeping
  double avmol, wsum, asum, gravity, totalAbund;
  int    IXH1,IXH2,IXHMIN,IXHE1,IXHE2,IXHE3,IXC1,IXAL1,IXSI1, 
         IXSI2,IXCA1,IXCA2,IXMG1,IXMG2,IXFE1,IXN1,IXO1;

  std::vector<float> fract, pf, potion, xamass;
  std::vector<std::string> species;
  std::vector<int>   idxspec, ion;
  std::vector<int> idx, idx1;
  std::vector<std::string> uspec;

  float xne, xna, RHOest;
  std::vector<char> totallist;
  int ntotallist;

  mat<float> buf;
  
  // Functions
  
  ceos(double grav = 4.44);
  ceos(std::vector<line_t> &lines, double grav = 4.44);
  ceos(std::vector<line_t> &lines, std::vector<iabund> &ab, double grav = 4.44);
  ceos(std::vector<line_t> &lines, std::string &abfile, double grav = 4.44);

  void initAbundances(std::vector<iabund> &ab, bool verbose = false);
  void initEOS(std::vector<line_t> &lines);
  
  
  ~ceos(){
    totallist.clear();
  }

  double nne_from_T_Pg    (double T,  double Pg,  double &rho, double iPe = -1.0);
  double nne_from_T_rho   (double T, double &iPg,  double rho);
  double rho_from_T_pel   (double T, double &iPg,   double Pe, float tol = 1.0e-5);
  double rho_from_T_nne   (double T, double &iPg,  double nne, float tol = 1.0e-5);

  double nne_from_T_Pg_nne (double T,  double Pg,  double &rho, double nne);
  double nne_from_T_rho_nne(double T, double &iPg,  double rho, double nne);
  double nne_from_T_rho_nne_old(double T, double &iPg,  double rho, double nne, float tol=1.e-5);
  float nne_from_T_rho_old   (float  T, float  &iPg,   float rho, float tol = 1.e-5);
  double nne_from_T_rho_old  (double T, double &iPg,  double rho, float tol = 1.e-5);

  
  float nne_from_T_Pg    (float T, float Pg,  float &rho, float Pe = -1.0);
  float nne_from_T_rho   (float T, float &Pg,  float rho);
  float rho_from_T_pel   (float T, float &Pg,   float Pe, float tol = 1.0e-5);
  float rho_from_T_nne   (float T, float &Pg,  float nne, float tol = 1.0e-5);

  float nne_from_T_Pg_nne (float T,  float Pg,  float &rho, float nne);
  float nne_from_T_rho_nne(float T, float &iPg,  float rho, float nne);

  
  void  contOpacity_TPg  (double T, double Pg, int nw, double *w, double *opac,
			  double *scattering, double Pe=-1.0);
  
  void  contOpacity_TRho  (double T, double rho, int nw, double *w, double *opac,
			  double *scattering, double Pe=-1.0);

  void contOpacity        (double T, int nw, double *w,
			   double *opac, double *scattering, std::vector<float> &frac, float na, float ne);

  void  hydrostatic      (int ndep, double *tau, double *t, double *Pg, double *rho,
			  double *nel, double *pel,  double pgas_bound, float tol = 1.0e-5);

  void  hydrostatic      (int ndep, double *tau, double *t, double *Pg, double *rho,
			  double *nel, double *pel,  double *z, double *cmass, double pgas_bound, float tol = 1.0e-5);

  void hydrostatic_cmass(int ndep, double *tau, double *t, double *Pg, double *rho, double *nel,
			 double *z, double *cmass, double *ltau, double &pgas_bound);
  
  void  hydrostatic      (int ndep, float *tau, float *t, float *Pg, float *rho,
			  float *nel, float *pel, float pgas_bound, float tol = 1.0e-5);

  
  void store_partial_pressures(int ndep, int k, float na, float ne);
  void read_partial_pressures(int k, std::vector<float> &frac, std::vector<float> &part, float &xa, float &xe);
  void unique(void);
  void fill_densities(int ndep, double *t, double *pgas, double *rho, double *pel,
		      double *nne, int touse, int keep_nne = 0, float tol = 1.0e-5);

  // float test_pgas_from_rho(float T, float &Pg, float rho,  float &nne);
  float init_pe_from_T_pg(float t, float pg);

  void readAbund(std::string file);
  
};


#endif
