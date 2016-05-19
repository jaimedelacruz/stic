ARCH     = $(shell uname -m)
ARCHF    = $(shell uname -a)

include makefiles/makefile.$(CPU).$(OS)

RHFOLD = rh/

INCLUDE = $(shell ncxx4-config --cflags) -I$(GSLPATH)/include -I./ -I$(RHFOLD) -I$(RHFOLD)/rh_1d
LIBS_FDENS = $(shell ncxx4-config --libs) -lstdc++ -L$(GSLPATH)/lib/
LIBS = $(LIBS_FDENS) -lgfortran -L./ -L$(RHFOLD)/rh_1d/ -lrhf1d -L$(GSLPATH) -lblas -llapack  -lgsl -lgslcblas -lfftw3 # -lifcore -static-intel #-fuse-linker-plugin 



FFILES = eos_math_special.o eos_eqns.o eos.o contop.o contop_3d.o
OFILES_SPARSE = input.o clm.o ceos.o comm.o depthmodel.o wavelet.o  hinode.o spectral.o atmosphere.o clte.o crh.o io.o slave.o  sparse.o master_sparse.o main_sparse.o 

FDENS = ceos.o io.o depthmodel.o fillDensities.o


.SUFFIXES: .o .f90 .cc

%.o: %.cc *.h
	$(MPCXX) $(CXXFLAGS) $(OPENMP) $(INCLUDE) -c  $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c  $< -o $@

RHinv: $(FFILES) $(OFILES_SPARSE)
#	$(MPFC)  -o RHinv.x $(FFILES) $(OFILES_SPARSE) $(INCLUDE) $(LIBS) -lmpi_cxx -fopenmp -nofor-main # Use with Intel Fortran
	$(MPCXX) -o RHinv.x  $(CXXFLAGS) $(OPENMP) $(FFILES) $(OFILES_SPARSE) $(INCLUDE) $(LIBS) #-nofor-main # USE with gfortran

clean:
	rm -f *.o *.mod

fillDensities: $(FFILES) $(FDENS)
	$(MPCXX) -o fillDensities.x $(CXXFLAGS) $(FFILES) $(OPENMP) $(FDENS) $(LIBS_FDENS)  $(INCLUDE) #-nofor-main

