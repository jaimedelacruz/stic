#ifndef CONTOP3D_H
#define CONTOP3D_H

/* ---
   Background opacities, adapted from N. Piskunov's fortran and C routines.
   
   ---
   Modifications:
        2016-09-10, JdlCR: Created, as float64.
	2016-09-26, JdlCR: Minor bug fixed.
	
   --- */

double SEATON(double FREQ0, double XSECT, double POWER, double A, double FREQ);
double COULFF(double TLOG, double FREQLG, int NZ);
double COULX(int N, double freq, double Z);
double C1OP(double FREQ, double TKEV);
double Mg1OP(double FREQ, double FREQLG, double T, double TLOG);
double Al1OP(double freq);
double Si1OP(double FREQ, double FREQLG, double T, double TLOG);
double Fe1OP(double FREQ, double HKT);
double N1OP(double FREQ, double TKEV);
double O1OP(double FREQ);
double Mg2OP(double FREQ, double TKEV);
double Si2OP(double FREQ, double FREQLG, double T, double TLOG);
double Ca2OP(double FREQ, double TKEV);

void HOP(double &AHYD, double XNE, double XH1, double XH2, double FREQ, double FREQLG,
	 double T, double TLOG, double TKEV, double STIM, double EHVKT);
void HRAYOP(double &sigh, double XH1, double FREQ);
void H2PLOP(double &AH2P, double XH1, double XH2, double FREQ, double FREQLG,
	     double FREQ15, double TKEV, double STIM);
void HMINOP(double &AHMIN, double XH1, double XHMIN, double FREQ, double T,
	    double TKEV, double XNE, double EHVKT);
void HE1OP(double &AHE1, double XHE1, double XHE2, double XNE, double FREQ, double FREQLG,
	   double T, double TKEV, double TLOG, double EHVKT, double STIM);
void HE2OP(double &AHE2, double XHE2, double XHE3, double XNE, double FREQ, double FREQLG,
	   double T, double TKEV, double TLOG, double EHVKT, double STIM);
void HEMIOP(double &AHEMIN, double XHE1, double FREQ, double T, double XNE);
void HERAOP(double &SIGHE, double XHE1, double FREQ);
void COOLOP(double &ACOOL, double XC1, double XMg1, double XAl1, double XSi1, double XFe1, double STIM,
	    double FREQ, double FREQLG, double T, double TLOG, double TKEV, double HKT);
void LUKEOP(double &ALUKE,double XN1, double XO1, double XMg2, double XSi2, double XCa2, double STIM,
	    double FREQ, double FREQLG, double T, double TLOG, double TKEV);
void HOTOP(double &AHOT);
void ELECOP(double &SIGEL, double XNE);
void H2RAOP(double &SIGH2, double XH1, double FREQ, double T, double TKEV, double TLOG);


void cop(double T, double TKEV, double TK, double HKT, double TLOG,
	 double XNA, double XNE, double *WLGRID, double *OPACITY,
	 double *SCATTER,  
	 double H1, double H2, double HMIN, double HE1, double HE2,
	 double HE3, double C1, double AL1, double SI1, double SI2,
	 double CA1, double CA2, double MG1, double MG2, double FE1,
	 double N1, double O1, int nWLGRID, int NLINES, int NTOTALLIST);


#endif
