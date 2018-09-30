# STiC: the Stockholm inversion code
Developed by J. de la Cruz Rodriguez, J. Leenaarts, S. Danilovic & H. Uitenbroek.

STiC is a MPI-parallel non-LTE inversion code for observed full-Stokes observations.
It allows inverting lines from multiple active atoms, including partial redistribution effects of scattered photons.

## Support
The code is distributed "as it is".

Although we encourage people to submit bug-reports, dedicated support will only be provided within a scientific collaboration with us, as we have very limited man power.


## Citing STiC
If you have used STiC in your research, please add the following description in the manuscript (or a very similar one,
including all references) to acknowledge all the work that has been done by different scientists over the past years:

*STiC (de la Cruz Rodriguez et al. 2018; de la Cruz et al. 2016) is a MPI-parallel non-LTE inversion code*
*that utilises a modified version of RH (Uitenbroek 2001) to solve the atom population*
*densities assuming statistical equililibrium and plane-parallel geometry and it allows including partial*
*redistribution effects of scattered photons (Leenaarts et al. 2012). The radiative transport equation is*
*solved using cubic Bezier solvers (de la Cruz Rodriguez et al. 2013).*

*The inversion engine of STiC includes an equation of state extracted from the SME code (Valenti & Piskunov 2016).*

The references can be found here:
* [de la Cruz Rodriguez et al. (in prep.)](https://dubshen.astro.su.se/~jaime/stic.pdf)
* [de la Cruz Rodriguez et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...830L..30D)
* [de la Cruz Rodriguez et al. (2013)](http://adsabs.harvard.edu/abs/2013ApJ...764...33D)
* [Uitenbroek (2001)](http://adsabs.harvard.edu/abs/2001ApJ...557..389U)
* [Leenaarts et al. (2012)](http://adsabs.harvard.edu/abs/2012A%26A...543A.109L)
* [Valenti & Piskunov (2016)](http://adsabs.harvard.edu/abs/2017A%26A...597A..16P)


## Dependencies
We have tested STiC with the GCC compilers and the Intel compilers.

The code includes source written in C, C++-11 and Fortran.

It makes use of the following libraries: Eigen-3, FFTW-3, netCDF4-cxx4, openmpi-2 (or any other MPI-2 implementation).

The code has been extensively run in Linux and OSX.

## Compilation instructions

Very detailed compilation instructions for OSX and Linux can be found [here](http://dubshen.astro.su.se/~jaime/stic_compile.txt). However, we include a brief explanation,
assuming that all dependencies have already been installed:

We have prepared different makefiles for different platforms and operating systems.
in your $HOME/.bashrc file and start a new terminal to make the changes effective).
These variables are just labels for the makefiles:

In linux:
```
export OS=Linux
export CPU=x86_64
```

In OSX:
```
export OS=Darwin
export CPU=i386
```

Additionally, make openmpi aware of your compilers:
```
export OMPI_CC=gcc
export OMPI_CXX=g++
export OMPI_FC=gfortran

source ~/.bashrc
```

These variables are just labels for the makefile, the code will be compiled in 64 bit
mode anyway.

STiC is based on a modified version the excellent RH code (Uitenbroek 2001).
We have encapsulated RH in a module that needs to be compiled first:
```
cd stic/src/rh
make clean
make

cd rh_1d/
make clean
make

cd ../../
make clean
make
```
If everything went fine you will find the binary of STiC in the main src folder (STiC.x).

You can try to execute it and see if it starts. You should get something like this:
```
   SSSSSSSSSSSSSSS TTTTTTTTTTTTTTTTTTTTTTT  iiii         CCCCCCCCCCCCC
 SS:::::::::::::::ST:::::::::::::::::::::T i::::i     CCC::::::::::::C
S:::::SSSSSS::::::ST:::::::::::::::::::::T  iiii    CC:::::::::::::::C
S:::::S     SSSSSSST:::::TT:::::::TT:::::T         C:::::CCCCCCCC::::C
S:::::S            TTTTTT  T:::::T  TTTTTTiiiiiii C:::::C       CCCCCC
S:::::S                    T:::::T        i:::::iC:::::C
 S::::SSSS                 T:::::T         i::::iC:::::C
  SS::::::SSSSS            T:::::T         i::::iC:::::C
    SSS::::::::SS          T:::::T         i::::iC:::::C
       SSSSSS::::S         T:::::T         i::::iC:::::C
            S:::::S        T:::::T         i::::iC:::::C
            S:::::S        T:::::T         i::::i C:::::C       CCCCCC
SSSSSSS     S:::::S      TT:::::::TT      i::::::i C:::::CCCCCCCC::::C
S::::::SSSSSS:::::S      T:::::::::T      i::::::i  CC:::::::::::::::C
S:::::::::::::::SS       T:::::::::T      i::::::i    CCC::::::::::::C
 SSSSSSSSSSSSSSS         TTTTTTTTTTT      iiiiiiii       CCCCCCCCCCCCC

STIC: Initialized with 1 process(es)
file_check: ERROR, file  does not exist!
```

## Example inversion test
We have included an example to test the code (stic/example/).
The example is basically one spectrum from an observation of an active region on the Sun in the
Ca II K line, Ca II 8542 Å line, Fe I 6301 & 6302 Å lines. The observations were
acquired with the CRISP and CHROMIS instruments at the Swedish 1-m Solar Telescope.

To run the example we first need link all the python tools that are included in the repository:

```
cd stic/example/
ln -s ../pythontools/py2/* .
```

Now we can execute the script that prepares the data. This script takes the line profile in each of the
spectral lines, and it stores them in a netCDF4 container so the code can read them. It also generates
an ad-hoc initial model and the instrumental profiles of our instruments. The script is extensively commented.

```
python prepare_data.py
```

The main files to change the behavior of STiC are input.cfg, keyword.input, atoms.inputm kurucz.input.

Since our example is only one pixel, there is no need to run the parallel version.

You can run the default example by simply typing:
```
mpiexec -n 1 ../src/STiC.x
```
...and wait. The code will provide information of each Levenberg-Marquardt iteration as they are completed.

Once the inversion is completed, you can visualize the result with the plot.py script:
```
python plot.py
```

## Acknowledgements

We are gratefull to N. Piskunov for allowing us to use his excellent EOS in our code.

JdlCR is supported by grants from the Swedish Research Council (2015-03994), the Swedish National Space Board (128/15). This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (SUNMAG, grant agreement 759548).

SD and JdlCR are supported by a grant from the Swedish Civil Contingencies Agency (MSB).

This research was supported by the CHROMATIC grant of the Knut and Alice Wallenberg foundation.
