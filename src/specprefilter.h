#ifndef SPECPREF
#define SPECPREF

class specprefilter: public instrument{
 public:
  region_t reg;
  bool firsttime;
  std::vector<double> iprof;
  

    
  /* --- constructor/Destructor --- */

 specprefilter(): reg(), firsttime(true), iprof(){};
  specprefilter(region_t &in, int nthreads = 1);
  ~specprefilter();
  
  /* --- Prototypes --- */
  //void degrade_one(spec_ft &ft, double *dat, int ns);
  //void degrade(mat<double> &syn, bool spectral = true, bool spatial = false, int ntt = 1);
  void degrade(double *syn, int ns);
  void init(int npsf, double const *psf);
  void update(size_t npsf, double *psf){init(int(npsf), psf);}
    
};

#endif
