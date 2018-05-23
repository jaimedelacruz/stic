CPU     = $(shell uname -m)
OS       = $(shell uname -s)

include makefiles/makefile.$(CPU).$(OS)

RHFOLD = rh/

INCLUDE = $(shell ncxx4-config --cflags) -I$(GSLPATH)/include -I./ -I$(RHFOLD) -I$(RHFOLD)/rh_1d
LIBS_FDENS = $(shell ncxx4-config --libs) -lstdc++ -L$(GSLPATH)/lib/
LIBS = $(LIBS_FDENS) -L./ -L$(RHFOLD)/rh_1d/ -lrhf1d -L$(GSLPATH) -lfftw3 -llapack  -lblas  -lgslcblas

MNAME = $(shell uname -n)
STMAC = STiC_$(MNAME).x

FFILES = eos_math_special.o eos_eqns.o eos.o 
OFILES_SPARSE = input.o clm.o cop.o ceos.o comm.o depthmodel.o  spectral.o fpigen.o fpi.o atmosphere.o clte.o crh.o io.o slave.o master_sparse.o main_sparse.o 

FDENS = cop.o ceos.o io.o depthmodel.o fillDensities.o


.SUFFIXES: .o .f90 .cc

%.o: %.cc *.h
	$(MPCXX) $(CXXFLAGS) $(OPENMP) $(INCLUDE) -c  $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c  $< -o $@

STiC: $(FFILES) $(OFILES_SPARSE)
	$(LINKER) $(LINKEROPTS) -o STiC.x  $(OPTS) $(OPENMP) $(FFILES) $(OFILES_SPARSE) $(INCLUDE) $(LIBS) $(LINKEROPTS)

machine: $(FFILES) $(OFILES_SPARSE)
	$(LINKER) $(LINKEROPTS) -o $(STMAC)  $(OPTS) $(OPENMP) $(FFILES) $(OFILES_SPARSE) $(INCLUDE) $(LIBS) $(LINKEROPTS)

clean:
	rm -f *.o *.mod

fillDensities: $(FFILES) $(FDENS)
	$(LINKER) -o fillDensities.x $(CXXFLAGS) $(FFILES) $(OPENMP) $(FDENS) $(LIBS_FDENS)  $(INCLUDE)  $(LINKEROPTS)

