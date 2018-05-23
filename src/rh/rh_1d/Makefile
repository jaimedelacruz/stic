
GEOMETRY = ONE_D_PLANE

include ../../makefiles/makefile.$(CPU).$(OS)

RHINC = -I../

## --- Define groups of object files --                -------------- ##

RH_OBS = ../hui_.o ../humlicek_.o ../abundance.o ../accelerate.o ../background.o ../backgropac_xdr.o ../barklem.o ../broad.o ../brs_xdr.o ../chemequil.o ../cocollisions.o ../collision.o ../complex.o ../cubeconvol.o ../duplicate.o ../error.o ../expint.o ../expspline.o ../fillgamma.o ../fixedrate.o ../fpehandler.o ../gammafunc.o ../gaussleg.o ../getcpu.o ../getlambda.o ../getline.o ../giigen.o ../h2collisions.o ../hunt.o ../humlicek.o ../hydrogen.o ../initial_xdr.o ../initscatter.o ../kurucz.o ../linear.o ../ltepops.o ../ludcmp.o ../matrix.o ../maxchange.o ../metal.o ../molzeeman.o ../nemetals.o ../ohchbf.o ../opacity.o ../options.o ../order.o ../parse.o ../paschen.o ../planck.o ../pops_xdr.o ../profile.o ../radrate_xdr.o ../readatom.o ../readb_xdr.o ../readj.o ../readvalue.o ../rayleigh.o ../readinput.o ../readmolecule.o  ../solvelinear.o ../solvene.o  ../spline.o ../statequil.o ../stokesopac.o ../stopreq.o ../thomson.o ../vacuumtoair.o ../voigt.o ../w3.o ../wigner.o ../writeatmos_xdr.o ../writeatom_xdr.o ../writecoll_xdr.o ../writedamp_xdr.o ../writeinput_xdr.o ../writemetal_xdr.o ../writemolec_xdr.o ../writeopac_xdr.o ../writespect_xdr.o ../zeeman.o  

ONE_D_OBJS = initial_j.o redistribute_j.o scatter_j.o sortlambda_j.o anglequad.o feautrier.o  formal.o  hydrostat.o  \
             piecestokes.o  piecewise.o bezier.o project.o  riiplane.o \
             writeflux_xdr.o  writegeom_xdr.o dummyatmos.o  \
	     background_j.o iterate_j.o calculateRay.o updateatmos.o interpol.o  rhf1d.o


.SUFFIXES: .o .f90 .cc

%.o: %.c *.h
	$(CC) $(RHINC) $(CFLAGS) -c  $< -o $@

%.o: %.f90
	$(FC) $(RHINC) $(FFLAGS) -c  $< -o $@

%.o: %.cc *.h
	$(CXX) $(RHINC) $(CXXFLAGS) -c  $< -o $@

librhf1d.a:  $(ONE_D_OBJS)
	$(AR) rcs $@  $(RH_OBS) $(ONE_D_OBJS) 

## --- Clean up --                                     -------------- ##

clean:
	rm -f *.o  librhf1d.a


## --- Explicit dependencies on include files --       -------------- ##


anglequad.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h

backgrcontr.o:          ../rh.h  ../atom.h  ../atmos.h  \
                        ../spectrum.h  geometry.h  ../background.h  \
                        ../inputs.h  ../error.h  ../statistics.h

background_j.o:         ../rh.h          ../atom.h        ../atmos.h       \
                         ../spectrum.h   ../constant.h    ../background.h  \
                        ../error.h       ../statistics.h  ../inputs.h      \
                        geometry.h           \
                        background_j.c

feautrier.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../error.h

formal.o:               ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../constant.h  \
                        ../background.h  ../inputs.h  ../error.h  ../xdr.h

hydrostat.o:            ../rh.h  ../atom.h  ../atmos.h  geometry.h \
                        ../background.h  ../constant.h  ../accelerate.h \
                        ../error.h

iter_1d.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../accelerate.h  ../error.h  \
                        ../statistics.h

multiatmos.o:           ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../constant.h  ../error.h  ../inputs.h  \
                        ../statistics.h  ../xdr.h

piecestokes.o:          ../rh.h ../error.h ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h

piecewise.o:            ../rh.h ../error.h ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h

project.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h

rhf1d.o:                ../rh.h  ../atom.h  ../atmos.h   \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../error.h  ../inputs.h  ../xdr.h

riiplane.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../constant.h  ../error.h

solve1d.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../accelerate.h  \
                        ../constant.h  ../statistics.h  ../inputs.h

solveray.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h  ../error.h  ../xdr.h

updateatmos.o:		../rh.h ../atmos.h ../atom.h geometry.h ../accelerate.h \
			../constant.h ../error.h ../inputs.h ../spectrum.h \
			rhf1d.h

writeflux_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../constant.h  \
                        ../error.h  ../xdr.h

writegeom_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../error.h  ../inputs.h  ../xdr.h

main.o:			io.h ceos.h rhf1d.h cmemt.h physical_consts.h
## ------- end ---------------------------- Makefile ---------------- ##
