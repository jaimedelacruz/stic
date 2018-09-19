# STiC: the Stockholm inversion code
Developed by J. de la Cruz Rodriguez with contributions from J. Leenaarts, S. Danilovic & H. Uitenbroek.

STiC is a MPI-parallel non-LTE inversion code for observed full-Stokes observations.
It allows inverting lines from multiple active atoms, including partial redistribution effects of scattered photons.


## Support
The code is distributed "as it is".
Although we encourage people to submit bug-reports, dedicated support will only be provided within a scientific collaboration, as we have very limited man power.


## Citing STiC
If you have used STiC in your research, please add the following description in the manuscript (or a very similar one,
including all references) to acknowledge all the work that has been done by different scientists over the past years:
```
STiC (de la Cruz Rodriguez et al. 2018; de la Cruz et al. 2016) is a MPI-parallel non-LTE inversion code
that utilises a modified version of RH (Uitenbroek 2001) to solve the atom population
densities assuming statistical equililibrium and plane-parallel geometry and it allows including partial
redistribution effects of scattered photons (Leenaarts et al. 2012). The radiative transport equation is
solved using cubic Bezier solvers (de la Cruz Rodriguez et al. 2013).

The inversion engine of STiC includes an equation of state extracted from the SME code (Valenti & Piskunov 2016).
```
where:
* de la Cruz Rodriguez et al. (in prep.): no link yet
* [de la Cruz Rodriguez et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...830L..30D)
* [de la Cruz Rodriguez et al. (2013)](http://adsabs.harvard.edu/abs/2013ApJ...764...33D)
* [Uitenbroek (2001)](http://adsabs.harvard.edu/abs/2001ApJ...557..389U)
* [Leenaarts et al. (2012)](http://adsabs.harvard.edu/abs/2012A%26A...543A.109L)
* [Valenti & Piskunov (2016)](http://adsabs.harvard.edu/abs/2017A%26A...597A..16P)


## Dependencies
We have tested STiC with the GCC compilers and the Intel compilers.
The code includes source written in C, C++-11 and Fortran.

It makes use of the following libraries: Eigen-3, FFTW-3, netCDF4-cxx4.

