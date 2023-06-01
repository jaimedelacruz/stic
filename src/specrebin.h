#ifndef SPECREB
#define SPECREB

#include "spectral.h"

class specrebin: public instrument{
 public:
  region_t reg;
  int nt;
  bool firsttime;
  int reb;
  std::vector<spec_ft> ft;
  

    
  /* --- constructor/Destructor --- */

 specrebin():firsttime(true){};
  specrebin(region_t &in, int nthreads = 1);
  ~specrebin();
  
  /* --- Prototypes --- */
  void degrade_one(spec_ft &ft, double *dat, int ns);
  void degrade(mat<double> &syn, bool spectral = true, bool spatial = false, int ntt = 1);
  void degrade(double *syn, int ns);
  void init(int npsf, double const *psf);
  void update(size_t npsf, double *psf){init(int(npsf), psf);}
    
};



#endif
