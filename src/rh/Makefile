# ---- Makes objects for the RH library ---

include ../makefiles/makefile.$(CPU).$(OS)


OBS = abundance.o accelerate.o background.o backgropac_xdr.o barklem.o broad.o brs_xdr.o chemequil.o cocollisions.o collision.o complex.o cubeconvol.o duplicate.o error.o expint.o expspline.o fillgamma.o fixedrate.o fpehandler.o gammafunc.o gaussleg.o getcpu.o getlambda.o getline.o giigen.o h2collisions.o hunt.o humlicek.o hydrogen.o initial_xdr.o initscatter.o iterate.o kurucz.o linear.o ltepops.o ludcmp.o matrix.o maxchange.o metal.o molzeeman.o nemetals.o ohchbf.o opacity.o options.o order.o parse.o paschen.o planck.o pops_xdr.o profile.o radrate_xdr.o readatom.o readb_xdr.o readj.o readvalue.o rayleigh.o readinput.o readmolecule.o redistribute.o scatter.o solvelinear.o solvene.o sortlambda.o spline.o statequil.o stokesopac.o stopreq.o thomson.o vacuumtoair.o voigt.o w3.o wigner.o writeatmos_xdr.o writeatom_xdr.o writecoll_xdr.o writedamp_xdr.o writeinput_xdr.o writemetal_xdr.o writemolec_xdr.o writeopac_xdr.o writespect_xdr.o zeeman.o getcpu.o fpehandler.o hui_.o humlicek_.o


.SUFFIXES: .o .f90 .c

%.o: %.c *.h
	$(CC) $(CFLAGS) -c  $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c  $< -o $@

objects: $(OBS)

clean:
	rm -f *.o
