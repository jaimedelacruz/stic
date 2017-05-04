#ifndef RHF1DF_H
#define RHF1DF_H
#ifdef __cplusplus
extern "C" {
#endif
#include <rpc/types.h>
#include <stdbool.h>
#include <stdio.h>

  typedef struct{
    long *atom_file_pos;
  } rhinfo;
  
  typedef struct{
    bool_t allocated;
    double **chi_b, **eta_b, **sca_b, **chip_b;
  } rhbgmem;
  
  typedef struct{
    int nlambda, nrays;
    double *lambda, *I, *Q, *U, *V;
  } ospec;

  typedef struct{
    int nlambda, idx;
    double *rho;
  } crhprd;
  
  typedef struct{
    int nlevel, converged, nprd;
    double *n, *ntotal;
    crhprd *line;
  } crhatom;
  
  typedef struct{
    int nactive, ndep, nw;
    double *J, *J20, *lambda;
    crhatom *pop;
  } crhpop;
  
  typedef struct{
    int rank, verb, iter;
    bool_t stop;
    FILE *logfile;
    char filename[300];
  } MPI_t;
  
  void save_populations(crhpop *save_pop);
  void read_populations(crhpop *save_pop);
  void clean_saved_populations(crhpop *save_pop_ref);
  void UpdateAtmosDep(void);
  void Initvarious();
  void calculateRay(void);

  bool_t rhf1d(float muz, int rhs_ndep, double *rhs_T, double *rhs_rho, 
	       double *rhs_nne, double *rhs_vturb, double *rhs_v, 
	       double *rhs_B, double *rhs_inc, double *rhs_azi,
	       double *rhs_z, double *rhs_nhtot, double *rhs_ltau,
	       double *rhs_cmass, double gravity, bool_t stokes, ospec *sp,
	       crhpop *save_pop, int mynw, double *mylambda, int myrank, int savpop,
	       int iverbose, int *hydrostat);
#ifdef __cplusplus
}
# endif

#endif

