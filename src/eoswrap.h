#ifndef EOSWRAP_H
#define EOSWRAP_H

#include "physical_consts.h"
#include <vector>
#include <input.h>
#include <string>

struct iabund{
  char elem[3];
  float abund;
};

class eoswrap{
 public:
  static const int   MAX_ELEM = 99;
  static constexpr double bk = phyc::BK;

  size_t nspec, nz;
  float ABUND[MAX_ELEM];
  double tABUND;
  float xna, xne;
  int IXH1, IXHE1;
  
  eoswrap(){};
  ~eoswrap(){};
  //eoswrap(std::vector<line_t> &lines, std::string &abfile, double grav = 4.44);
  
  
  virtual double nne_from_T_Pg    (double T,  double Pg,  double &rho, double iPe = -1.0)= 0;
  virtual double nne_from_T_rho   (double T, double &iPg,  double rho) =0 ;
  virtual double rho_from_T_pel   (double T, double &iPg,   double Pe, float tol = 1.0e-5)=0;
  virtual double rho_from_T_nne   (double T, double &iPg,  double nne, float tol = 1.0e-5)=0;

  virtual double nne_from_T_Pg_nne (double T,  double Pg,  double &rho, double nne)=0;
  virtual double nne_from_T_rho_nne(double T, double &iPg,  double rho, double nne)=0;

  virtual void hydrostatic_cmass(int ndep, double *tau, double *t, double *Pg, double *rho, double *nel,
				 double *z, double *cmass, double *ltau, double &pgas_bound)=0;
  virtual void  hydrostatic      (int ndep, double *tau, double *t, double *Pg, double *rho,
				  double *nel, double *pel,  double *z, double *cmass, double pgas_bound, float tol = 1.0e-5) = 0;
  virtual void unique(void)=0;
  virtual void store_partial_pressures(int ndep, int k, float na, float ne)=0;
  virtual void read_partial_pressures(int k, std::vector<float> &frac, std::vector<float> &part, float &xa, float &xe)=0;

  virtual void contOpacity_TPg  (double T, double Pg, int nw, double *w, double *opac, double *scattering, double Pe=-1.0) = 0;
  virtual void  contOpacity_TRho  (double T, double rho, int nw, double *w, double *opac,
				   double *scattering, double Pe=-1.0) = 0;

  virtual void fill_densities(int ndep, double *t, double *pgas, double *rho, double *pel,
			double *nne, int touse,  int keep_nne, float tol) = 0;

  virtual   void contOpacity        (double T, int nw, double *w,
				     double *opac, double *scattering, std::vector<float> &frac, float na, float ne) = 0;
};



#endif
