/* ----
   WITTMANN EOS routines, adapted from the SIR code implementation to C++
   J. de la Cruz Rodriguez (ISP-SU)

   Added interface routines to the background opacity code.
   
   

   ---- */

#ifndef WITT_H
#define WITT_H

#include <vector>
#include <cstddef>
#include <string>
#include <cmath>
#include "cop.h"
#include "eoswrap.h"
#include "input.h"
#include "physical_consts.h"

namespace eos{

  static const int Nelements = 99;
  static const float ABUND_default[Nelements] = {
    -0.04048,-1.07,-10.95,-10.89, -9.44, -3.48, -3.99, 
    -3.11, -7.48, -3.95,  -5.71, -4.46, -5.57, -4.49, -6.59, -4.83, 
    -6.54, -5.48, -6.82,  -5.68, -8.94, -7.05, -8.04, -6.37, -6.65, 
    -4.50, -7.12, -5.79,  -7.83, -7.44, -9.16, -8.63, -9.67, -8.69, 
    -9.41, -8.81, -9.44,  -9.14, -9.80, -9.54,-10.62,-10.12,-20.00, 
    -10.20,-10.92,-10.35, -11.10,-10.18,-10.58,-10.04,-11.04, -9.80, 
    -10.53, -9.81,-10.92,  -9.91,-10.82,-10.49,-11.33,-10.54,-20.00, 
    -11.04,-11.53,-10.92, -11.94,-10.94,-11.78,-11.11,-12.04,-10.96, 
    -11.28,-11.16,-11.91, -10.93,-11.77,-10.59,-10.69,-10.24,-11.03, 
    -10.95,-11.14,-10.19, -11.33,-20.00,-20.00,-20.00,-20.00,-20.00, 
    -20.00,-11.92,-20.00, -12.51,-20.00,-20.00,-20.00,-20.00,-20.00, 
   -20.00,-20.00};
  
  
  
  class witt: public eoswrap{
  private:
    static const float prec;
    static const int ncontr;
    double Ab_others;

  public:
    double avweight, muH, rho_from_H, gravity;
    static const double saha_fac;
    
    witt(){};
    witt(std::vector<line_t> &lines, std::string &abFile, double grav = 4.44);
    ~witt(){};


    
    /* --- Prototypes of instantiated templates (in .cc file) --- */
    
    template <class T> void acota(T &x, T x0, T x1);
    template <class T> void acotasig(T &x, T x0, T x1);
    
    template <class T> T pe_from_pg(T temp, T Pgas, T *fe_out = NULL);
    template <class T> T pg_from_pe(T t, T Pe, T *fe_out = NULL);
    template <class T> T rho_from_pg(T temp, T Pgas, T &Pe);
    template <class T> T rho_from_pe(T temp, T Pe, T &Pg);
    template <class T> T pg_from_rho(T temp, T rho, T &Pe);

    template <class T> T pe_pg(T temp, T Pgas, T Pe, T *fe_out = NULL);
    template <class T> T saha(T theta, T chi, T u1, T u2, T Pe);
    template <class T> T dsaha(T theta, T chi, T du1, T du2);
    template <class T> void partition_f(int nel, T t, T &u1, T &u2, T &u3,
					T &du1, T &du2, T &du3);
    template <class T> void molecb(T X, T *Y, T *dy);
    template <class T> T init_pe_from_T_Pg(T t, T pg);
    //    template <class T> void gasb(T theta, T pe, T *p, T *dp, T *dpp);
    template <class T> void gasc(T t, T pe, T &pg, T *pp, T *fe_out = NULL);
    template <class T> T get_Nlow_over_U(int anum, int ion);
    template <class T> T getN_and_U(int iatom, int istage, T t, T Pgas, T Pe, T &u, bool divide_by_u = false);
    template <class T> T nsaha(T t, T xne, T u0, T u1, T eion);
    template <class T> void getBackgroundPartials(T t, T Pgas, T Pe, double *n_partial);
    template <class T> void contOpacity(T temp, T iPgas, T iPe, int nw, double *w, double *opac);
    template <class T> int getXpart(int iatom, T t, T Pgas, T Pe, T *xpa, T *u, T *ein,
				    bool divide_by_u = false);
    

    /* --- Prototypes for eoswrap --- */
    
    double Pe_from_Pg(double temp, double Pgas, double *fe_out = NULL){return pe_from_pg(temp, Pgas, fe_out);};
    double Pg_from_Pe(double temp, double Pe, double *fe_out = NULL){return pg_from_pe(temp, Pe, fe_out);};
    double Pg_from_Rho(double temp, double rho, double &Pe){return pg_from_rho(temp, rho, Pe);};
    void unique(void){};
    void store_partial_pressures(int ndep, int k, float na, float ne){};
    void read_partial_pressures(int k, std::vector<float> &frac, std::vector<float> &part, float &xa, float &xe){};
    
    /* --- header templates --- */
    
    template <class T> T sign(T a, T b){
      return fabs(a) * b / fabs(b);
    }


    // --- Prototypes and wrappers for STiC --- //
    
    double nne_from_T_Pg    (double T,  double Pg,  double &rho, double iPe = -1.0)
    {
      double Pe = 0.0;
      rho = rho_from_pg<double>(T, Pg, Pe);
      return Pe/(T*phyc::BK);
    }
    double nne_from_T_rho   (double T, double &iPg,  double rho)
    {
      double Pe;
      iPg = pg_from_rho<double>(T, rho, Pe);
      return Pe/(T*phyc::BK);
    }
    double rho_from_T_pel   (double T, double &iPg,   double Pe, float tol = 1.0e-5)
    {
      double rho = 0.0;
      rho_from_pe<double>(T, Pe, iPg);
      return rho;
    }
    double rho_from_T_nne   (double T, double &iPg,  double nne, float tol = 1.0e-5)
    {
      double Pe = nne * phyc::BK * T;
      double rho = rho_from_pe<double>(T, Pe, iPg);
      return rho;
    }
    double nne_from_T_Pg_nne (double T,  double Pg,  double &rho, double nne)
    {
      double Pe = 0.0;
      rho = rho_from_pg<double>(T, Pg, Pe);
      return nne;
    }
    double nne_from_T_rho_nne(double T, double &iPg,  double rho, double nne)
    {
      double Pe;
      iPg = pg_from_rho<double>(T, rho, Pe);
      return nne;
    }
    void contOpacity_TPg  (double T, double Pg, int nw, double *w, double *opac,
			   double *scattering, double Pe=-1.0)
    {
      Pe = pe_from_pg<double>(T,Pg, NULL);
      contOpacity<double>(T, Pg, Pe, nw, w, opac);
    }
    void  contOpacity_TRho  (double T, double rho, int nw, double *w, double *opac,
			     double *scattering, double Pe=-1.0)
    {
      double Pg = pg_from_rho<double>(T, rho, Pe);
      contOpacity<double>(T, Pg, Pe, nw, w, opac);
    }

    void  hydrostatic      (int ndep, double *tau, double *t, double *Pg, double *rho,
			    double *nel, double *pel,  double *z, double *cmass, double pgas_bound, float tol = 1.0e-5);
    void fill_densities(int ndep, double *t, double *pgas, double *rho, double *pel,
			double *nne, int touse,  int keep_nne, float tol);
    void readAbund(std::string file);
    void initAbundances(std::vector<iabund> &ab, bool verbose = false);
    void hydrostatic_cmass(int ndep, double *tau, double *t, double *Pg, double *rho, double *nel,
			   double *z, double *cmass, double *ltau, double &pgas_bound);

    void contOpacity        (double T, int nw, double *w,
			     double *opac, double *scattering, std::vector<float> &frac, float na, float ne){}



    
  };



  
}




#endif
