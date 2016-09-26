/* ---
   Background opacities, adapted from N. Piskunov's fortran and C routines.
   
   ---
   Modifications:
        2016-09-10, JdlCR: Created, as float64.
	2016-09-26, JdlCR: Minor bug fixed.

   --- */

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "cop.h"

/* ------------------------------------------------------------------------------ */

#define pow10(x) exp(2.302585093*(x))
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define sqr(a) (a*a)

/* ------------------------------------------------------------------------------ */

double SEATON(double FREQ0, double XSECT, double POWER, double A, double FREQ)
{
  return XSECT*(A+(1.-A)*(FREQ0/FREQ))*
         pow(FREQ0/FREQ, floor(2.*POWER+0.01)*0.5);
}

/* ------------------------------------------------------------------------------ */

double COULFF(double TLOG, double FREQLG, int NZ)
{
  static double Z4LOG[6]={0.,1.20412,1.90849,2.40824,2.79588,3.11261},
                A[12][11]={
     {5.53,5.49,5.46,5.43,5.40,5.25,5.00,4.69,4.48,4.16,3.85},
     {4.91,4.87,4.84,4.80,4.77,4.63,4.40,4.13,3.87,3.52,3.27},
     {4.29,4.25,4.22,4.18,4.15,4.02,3.80,3.57,3.27,2.98,2.70},
     {3.64,3.61,3.59,3.56,3.54,3.41,3.22,2.97,2.70,2.45,2.20},
     {3.00,2.98,2.97,2.95,2.94,2.81,2.65,2.44,2.21,2.01,1.81},
     {2.41,2.41,2.41,2.41,2.41,2.32,2.19,2.02,1.84,1.67,1.50},
     {1.87,1.89,1.91,1.93,1.95,1.90,1.80,1.68,1.52,1.41,1.30},
     {1.33,1.39,1.44,1.49,1.55,1.56,1.51,1.42,1.33,1.25,1.17},
     {0.90,0.95,1.00,1.08,1.17,1.30,1.32,1.30,1.20,1.15,1.11},
     {0.55,0.58,0.62,0.70,0.85,1.01,1.15,1.18,1.15,1.11,1.08},
     {0.33,0.36,0.39,0.46,0.59,0.76,0.97,1.09,1.13,1.10,1.08},
     {0.19,0.21,0.24,0.28,0.38,0.53,0.76,0.96,1.08,1.09,1.09}};
  double GAMLOG, HVKTLG, P, Q, CLFF;
  int IGAM, IHVKT;

/*  GAMLOG=log10(158000*Z*Z/T)*2 */

  GAMLOG=10.39638-TLOG/1.15129+Z4LOG[NZ-1];
  IGAM=min((int)(GAMLOG+7.),10); if(IGAM<1) IGAM=1;

/*  HVKTLG=2*log10(HVKT) */

  HVKTLG=(FREQLG-TLOG)/1.15129-20.63764;
  IHVKT=min((int)(HVKTLG+9.),11); if(IHVKT<1) IHVKT=1;
  P=GAMLOG-(IGAM-7);
  Q=HVKTLG-(IHVKT-9);
  CLFF=(1.-P)*((1.-Q)*A[IHVKT-1][IGAM-1]+Q*A[IHVKT][IGAM-1])+
       P*((1.-Q)*A[IHVKT-1][IGAM]+Q*A[IHVKT][IGAM]);
  return CLFF;
}

/* ------------------------------------------------------------------------------ */

double COULX(int N, double freq, double Z)
{
  static double A[6]={0.9916,1.105,1.101,1.101,1.102,1.0986},
                B[6]={2.719e3,-2.375e4,-9.863e3,-5.765e3,-3.909e3,-2.704e3},
                C[6]={-2.268e10,4.077e8,1.035e8,4.593e7,2.371e7,1.229e7};
  double CLX, FREQ1;
  int n;

  n=(N+1)*(N+1);
  if(freq>=Z*Z*3.28805e15/n)
  {
    FREQ1=freq*1.e-10;
    CLX=0.2815/FREQ1/FREQ1/FREQ1/n/n/(N+1)*Z*Z*Z*Z;
    if(N>=6) return CLX;
    CLX*=(A[N]+(B[N]+C[N]*(Z*Z/FREQ1))*(Z*Z/FREQ1));
    return CLX;
  }
  return 0.;
}

/* ------------------------------------------------------------------------------ */

void HOP(double &AHYD, double XNE, double XH1, double XH2, double FREQ, double FREQLG,
	 double T, double TLOG, double TKEV, double STIM, double EHVKT) /* REQUIRES FUNCTIONS COULX AND COULFF */
{
  int n1;
  double C, H, CFREE, FREQ3, FREET;
  double CONT[8], BOLT[8], EXLIM, BOLTEX, XR;

  FREQ3=FREQ*1.E-10; FREQ3 *= FREQ3*FREQ3;
  CFREE=3.6919E-22/FREQ3;

  for(int N = 0; N<8; N++){
    n1 = N+1; n1 *= n1;
    BOLT[N]=exp(-13.595*(1.-1./n1)/TKEV)*2.*n1*XH1;
  }
  
  FREET=XNE*CFREE*XH2/sqrt(T);
  XR=XH1/13.595*TKEV;
  BOLTEX=exp(-13.427/TKEV)*XR;
  EXLIM=exp(-13.595/TKEV)*XR;
  
  for(int N = 0; N<8; N++) CONT[N]=COULX(N,FREQ,1.0);
      
  C=0.2815/FREQ3;
  if(FREQ < 4.05933E13) BOLTEX=EXLIM/EHVKT;
  H=(CONT[6]*BOLT[6]+CONT[7]*BOLT[7]+(BOLTEX-EXLIM)*C+
     COULFF(TLOG,FREQLG,1)*FREET)*STIM;

  
  for(int N = 0; N<6; N++) H += CONT[N]*BOLT[N]*(1.-EHVKT);
  AHYD=H;

  return;
}

/* ------------------------------------------------------------------------------ */

void HRAYOP(double &sigh, double XH1, double FREQ)
{
  double WAVE, WW, SIG;

  WAVE=min(FREQ,2.463e15);
  WAVE=2.997925e18/WAVE;
  WW=WAVE*WAVE;
  SIG=(5.799e-13+1.422e-6/WW+2.784/(WW*WW))/(WW*WW);
  sigh=SIG*XH1*2.0;
  
  return;
}

/* ------------------------------------------------------------------------------ */

void  H2PLOP(double &AH2P, double XH1, double XH2, double FREQ, double FREQLG,
	     double FREQ15, double TKEV, double STIM)
{
  double ES,FR;
  
  AH2P=0.0;
  if(FREQ > 3.28805E15) return;
  FR=-3.0233E3+(3.7797E2+(-1.82496E1+(3.9207E-1-3.1672E-3*FREQLG)*
			  FREQLG)*FREQLG)*FREQLG;
  ES=-7.342E-3+(-2.409+(1.028+(-0.4230+(0.1224-0.01351*FREQ15)*
			       FREQ15)*FREQ15)*FREQ15)*FREQ15;
  AH2P=exp(-ES/TKEV+FR)*2.*XH1*XH2*STIM;
  
  return;
 }

/* ------------------------------------------------------------------------------ */

 void HMINOP(double &AHMIN, double XH1, double XHMIN, double FREQ, double T,
		   double TKEV, double XNE, double EHVKT)
 {

   double FREQ1,B,C,H,HMINBF,HMINFF,HMIN;

   FREQ1=FREQ*1.E-10;
   B=(1.3727E-15+4.3748/FREQ)/FREQ1;
   C=-2.5993E-7/pow(FREQ1,2);

   if(FREQ <= 1.8259E14) HMINBF=0.;
   else if(FREQ >= 2.111E14) HMINBF=6.801E-10+(5.358E-3+(1.481E3+(-5.519E7+4.808E11/FREQ1)/FREQ1)/FREQ1)/FREQ1;
   else HMINBF=3.695E-6+(-1.251E-1+1.052E3/FREQ1)/FREQ1;

   HMINFF=(B+C/T)*XH1*XNE*2.E-20;

   
  /*
    We use the number density / partition function for H-.
    The partition function for H- is 1 and not 2 as was mistakenly
    used before (fixed 2007-12-15, NP)! The EOS calculations are
    good up to 7730K, we use Kurucz approximation for higher T.
  */
   if(T < 7730.) HMIN=XHMIN;
   else HMIN=exp(0.7552/TKEV)/(2.*2.4148E15*T*sqrt(T))*XH1*XNE;
   
   H=HMINBF*(1-EHVKT)*HMIN*1.E-10;
   AHMIN=H+HMINFF;

   return;
 }

/* ------------------------------------------------------------------------------ */

void  HE1OP(double &AHE1, double XHE1, double XHE2, double XNE, double FREQ, double FREQLG,
	    double T, double TKEV, double TLOG, double EHVKT, double STIM)
{
  //     REQUIRES FUNCTION COULFF
  
  int NMIN;
  double BOLT[10],BOLTEX,XRLOG,HE1,EX,EXLIM;
  double FREQ3,FREET,CFREE,TRANS[10],C;
  
  static double G[10]={1.,3.,1.,9.,3.,3.,1.,9.,20.,3.};
  static double HEFREQ[10]={5.9452090e15,1.1528440e15,0.9803331e15,.8761076e15,
			    0.8147100e15,0.4519048e15,0.4030971e15,.8321191e15,
			    0.3660215e15,0.3627891e15};
  static double CHI[10]={0.,19.819,20.615,20.964,21.217,22.718,22.920,23.006,
			 23.073,23.086};

      
  for(int N=0; N<10; N++) BOLT[N]=exp(-CHI[N]/TKEV)*G[N]*XHE1;
  
  FREET=XNE*1.E-10*XHE2*1.E-10/sqrt(T)*1.E-10;
  XRLOG=log(XHE1*(2./13.595)*TKEV);
  BOLTEX=exp(-23.730/TKEV+XRLOG);
  EXLIM=exp(-24.587/TKEV+XRLOG);
  FREQ3=pow(FREQ*1.E-10,3);
  CFREE=3.6919E8/FREQ3;
  C=2.815E-1/FREQ3;
  
  for(NMIN=0; NMIN<10; NMIN++){
    TRANS[NMIN]=0.0;
    if(HEFREQ[NMIN] <= FREQ) break;
  }
  
  switch(NMIN)
    {
    case 0: TRANS[0]=exp(33.32-2.*FREQLG);
    case 1: TRANS[1]=exp(-390.026+(21.035-0.318*FREQLG)*FREQLG);
    case 2: TRANS[2]=exp(26.83-1.91*FREQLG);
    case 3: TRANS[3]=exp(61.21-2.9*FREQLG);
    case 4: TRANS[4]=exp(81.35-3.5*FREQLG);
    case 5: TRANS[5]=exp(12.69-1.54*FREQLG);
    case 6: TRANS[6]=exp(23.85-1.86*FREQLG);
    case 7: TRANS[7]=exp(49.30-2.60*FREQLG);
    case 8: TRANS[8]=exp(85.20-3.69*FREQLG);
    case 9: TRANS[9]=exp(58.81-2.89*FREQLG);
    default: break;
    }
  
  EX = BOLTEX;
  if(FREQ < 2.055E14) EX=EXLIM/EHVKT;
  HE1=(EX-EXLIM)*C;
  for(int N=0; N<10; N++) HE1 += TRANS[N]*BOLT[N];
  AHE1=(HE1+COULFF(TLOG,FREQLG,1)*FREET*CFREE)*STIM;

  return;
}

/* ------------------------------------------------------------------------------ */

void HE2OP(double &AHE2, double XHE2, double XHE3, double XNE, double FREQ, double FREQLG,
	   double T, double TKEV, double TLOG, double EHVKT, double STIM)
{
  /*
     REQUIRES FUNCTIONS COULX AND COULFF
     FREQUENCIES ARE 4X HYDROGEN,CHI ARE FOR ION POT=54.403
  */
  double HE2,C,FREQ3,EX,XR,CFREE,BOLT[9],CONT[9],
    EXLIM,FREET,BOLTEX;
  
  for(int N=0;N<9;N++) BOLT[N]=exp(-(54.403-54.403/sqr(N+1))/TKEV)*2.*sqr(N+1)*XHE2;
  
  FREET=XNE*XHE3/sqrt(T);
  XR=XHE2/13.595*TKEV;
  BOLTEX=exp(-53.859/TKEV)*XR;
  EXLIM=exp(-54.403/TKEV)*XR;
  
  for(int N=0;N<9;N++) CONT[N]=COULX(N,FREQ,2.0);
  
  FREQ3=pow(FREQ*1.E-5, 3);
  CFREE=3.6919E-07/FREQ3*4.0;
  C=2.815E14*2.0*2.0/FREQ3;
  EX=BOLTEX;
  
  if(FREQ < 1.31522E14) EX=EXLIM/EHVKT;
  HE2=(EX-EXLIM)*C;
  
  for(int N=0;N<9;N++) HE2+= CONT[N]*BOLT[N];
  
  HE2=(HE2+COULFF(TLOG,FREQLG,2)*CFREE*FREET)*STIM;
  
  if(HE2 >= 1.E-20) AHE2=HE2;
  else AHE2=0.0;
  
  return;
}

/* ------------------------------------------------------------------------------ */

void HEMIOP(double &AHEMIN, double XHE1, double FREQ, double T, double XNE)
{
  double A= 3.397E-26+(-5.216E-11+7.039E05/FREQ)/FREQ;
  double B=-4.116E-22+( 1.067E-06+8.135E09/FREQ)/FREQ;
  double C= 5.081E-17+(-8.724E-03-5.659E12/FREQ)/FREQ;
  AHEMIN=(A*T+B+C/T)*XNE*XHE1*1.E-20;
  return;
}

/* ------------------------------------------------------------------------------ */

void HERAOP(double &SIGHE, double XHE1, double FREQ)
{
  double WW=2.997925E+03/min(FREQ*1.E-15,5.15);
  WW *= WW;

  double arg = 1.+(2.44E5+5.94E10/(WW-2.90E5))/WW;
  double SIG=5.484E-14/WW/WW*arg*arg;
  SIGHE=SIG*XHE1;
  return;
}

/* ------------------------------------------------------------------------------ */

double Mg1OP(double FREQ, double FREQLG, double T, double TLOG)
{
  //     CROSS-SECTION TIMES THE PARTITION FUNCTION
  
  double DT,D,D1,XWL1,XWL2;
  int N, NT;
  static double PEACH[15][7]=
    {
      /* TEMP: 4000     5000     6000     7000     8000     9000    10000     WAVE(A) */
      {-42.474, -42.350, -42.109, -41.795, -41.467, -41.159, -40.883},/*  1500 */
      {-41.808, -41.735, -41.582, -41.363, -41.115, -40.866, -40.631},/*  1550 */
      {-41.273, -41.223, -41.114, -40.951, -40.755, -40.549, -40.347},/*  1621 */
      {-45.583, -44.008, -42.957, -42.205, -41.639, -41.198, -40.841},/*  1622 */
      {-44.324, -42.747, -41.694, -40.939, -40.370, -39.925, -39.566},/*  2513 */
      {-50.969, -48.388, -46.630, -45.344, -44.355, -43.568, -42.924},/*  2514 */
      {-50.633, -48.026, -46.220, -44.859, -43.803, -42.957, -42.264},/*  3756 */
      {-53.028, -49.643, -47.367, -45.729, -44.491, -43.520, -42.736},/*  3757 */
      {-51.785, -48.352, -46.050, -44.393, -43.140, -42.157, -41.363},/*  6549 */
      {-52.285, -48.797, -46.453, -44.765, -43.486, -42.480, -41.668},/*  6550 */
      {-52.028, -48.540, -46.196, -44.507, -43.227, -42.222, -41.408},/*  7234 */
      {-52.384, -48.876, -46.513, -44.806, -43.509, -42.488, -41.660},/*  7235 */
      {-52.363, -48.856, -46.493, -44.786, -43.489, -42.467, -41.639},/*  7291 */
      {-54.704, -50.772, -48.107, -46.176, -44.707, -43.549, -42.611},/*  7292 */
      {-54.359, -50.349, -47.643, -45.685, -44.198, -43.027, -42.418}};/* 9000 */
  
  static double FREQMG[7]={1.9341452e15,1.8488510e15,1.1925797e15,
                           7.9804046e14,4.5772110e14,4.1440977e14,
                           4.1113514e14};
  static double FLOG[9]={35.32123,35.19844,35.15334,34.71490,34.31318,
                         33.75728,33.65788,33.64994,33.43947};
  static double TLG[7]={8.29405,8.51719,8.69951,8.85367,
                        8.98720,9.10498,9.21034};
  
  NT=min(6,(int)floor(T/1000.)-3); if(NT<1) NT=1;
  
  DT=(TLOG-TLG[NT-1])/(TLG[NT]-TLG[NT-1]);
  for(N=0;N<7;N++) if(FREQ > FREQMG[N]) break;
  
  D=(FREQLG-FLOG[N])/(FLOG[N+1]-FLOG[N]);
  if(N > 1) N=2*N-1;
  D1=1.0-D;
  XWL1=PEACH[N+1][NT-1]*D + PEACH[N][NT-1]*D1;
  XWL2=PEACH[N+1][NT  ]*D + PEACH[N][NT  ]*D1;

  return exp(XWL1*(1.-DT)+XWL2*DT);
}

/* ------------------------------------------------------------------------------ */

double C1OP(double FREQ, double TKEV)
{
  // CROSS-SECTION TIMES THE PARTITION FUNCTION
  
  double C1240,C1444,X1444,X1240,X1100;
  
  C1240=5.*exp(-1.264/TKEV);
  C1444=exp(-2.683/TKEV);
  X1444=0.0;
  X1240=0.0;
  X1100=0.0;
  
  if(FREQ >= 2.7254E15) X1100=SEATON(2.7254E15,1.219E-17,2.0E0,3.317E0,FREQ);
  if(FREQ >= 2.4196E15) X1240=SEATON(2.4196E15,1.030E-17,1.5E0,2.789E0,FREQ);
  if(FREQ >= 2.0761E15) X1444=SEATON(2.0761E15,9.590E-18,1.5E0,3.501E0,FREQ);
  
  return X1100*9.+X1240*C1240+X1444*C1444;
}

/* ------------------------------------------------------------------------------ */

double Al1OP(double FREQ)
{
  if(FREQ > 1.443E15) return 2.1E-17*pow(1.443E15/FREQ,3)*6.0;
  else return 0.0;
}

/* ------------------------------------------------------------------------------ */

double Si1OP(double FREQ, double FREQLG, double T, double TLOG)
{
  //     CROSS-SECTION TIMES THE PARTITION FUNCTION
  
  int N,NT;
  double DT,DD,D,XWL1,XWL2;
  static double PEACH[19][9]=
  /* TEMP:4000   5000   6000   7000   8000   9000  10000  11000  12000   WAVE(A) */
    {{ 38.136,38.138,38.140,38.141,38.143,38.144,38.144,38.145,38.145},/* 1200 */
     { 37.834,37.839,37.843,37.847,37.850,37.853,37.855,37.857,37.858},/* 1400 */
     { 37.898,37.898,37.897,37.897,37.897,37.896,37.895,37.895,37.894},/* 1519 */
     { 40.737,40.319,40.047,39.855,39.714,39.604,39.517,39.445,39.385},/* 1520 */
     { 40.581,40.164,39.893,39.702,39.561,39.452,39.366,39.295,39.235},/* 1676 */
     { 45.521,44.456,43.753,43.254,42.878,42.580,42.332,42.119,41.930},/* 1677 */
     { 45.520,44.455,43.752,43.251,42.871,42.569,42.315,42.094,41.896},/* 1978 */
     { 55.068,51.783,49.553,47.942,46.723,45.768,44.997,44.360,43.823},/* 1979 */
     { 53.868,50.369,48.031,46.355,45.092,44.104,43.308,42.652,42.100},/* 5379 */
     { 54.133,50.597,48.233,46.539,45.261,44.262,43.456,42.790,42.230},/* 5380 */
     { 54.051,50.514,48.150,46.454,45.176,44.175,43.368,42.702,42.141},/* 5624 */
     { 54.442,50.854,48.455,46.733,45.433,44.415,43.592,42.912,42.340},/* 5625 */
     { 54.320,50.722,48.313,46.583,45.277,44.251,43.423,42.738,42.160},/* 6260 */
     { 55.691,51.965,49.444,47.615,46.221,45.119,44.223,43.478,42.848},/* 6261 */
     { 55.661,51.933,49.412,47.582,46.188,45.085,44.189,43.445,42.813},/* 6349 */
     { 55.973,52.193,49.630,47.769,46.349,45.226,44.314,43.555,42.913},/* 6350 */
     { 55.922,52.141,49.577,47.715,46.295,45.172,44.259,43.500,42.858},/* 6491 */
     { 56.828,52.821,50.110,48.146,46.654,45.477,44.522,43.730,43.061},/* 6492 */
     { 56.657,52.653,49.944,47.983,46.491,45.315,44.360,43.569,42.901}};/*6900 */
  /*     3P,1D,1S,1D,3D,3F,1D,3P */
  static double FREQSI[9]={2.1413750e15,1.97231650e15,1.7879689e15,
                           1.5152920e15,0.55723927e15,5.3295914e14,
                           4.7886458e14,4.72164220e14,4.6185133e14};
  static double FLOG[11]={35.45438,35.30022,35.21799,35.11986,34.95438,
                          33.95402,33.90947,33.80244,33.78835,33.76626,
                          33.70518};
  static double TLG[9]={8.29405,8.51719,8.69951,8.85367,8.98720,
                        9.10498,9.21034,9.30565,9.39266};
  
  NT=min(8,(int)floor(T/1000.)-3); if(NT<1) NT=1;
  DT=(TLOG-TLG[NT-1])/(TLG[NT]-TLG[NT-1]);
  
  for(N=0;N<9;N++) if(FREQ > FREQSI[N]) break;
      
  D=(FREQLG-FLOG[N])/(FLOG[N+1]-FLOG[N]);
    
  if(N>1) N=2*N-1;
  DD=1.-D;
  XWL1=PEACH[N+1][NT-1]*D+PEACH[N][NT-1]*DD;
  XWL2=PEACH[N+1][NT  ]*D+PEACH[N][NT  ]*DD;
  
  return exp(-(XWL1*(1.-DT)+XWL2*DT))*9.;
}

/* ------------------------------------------------------------------------------ */

double Fe1OP(double FREQ, double HKT)
{
  //     CROSS-SECTION TIMES PARTITION FUNCTION
  
  int I;
  double BOLT[48],XSECT[48],WAVENO, FE10P, XXX;
  
  static double G[48]={25.,35.,21.,15., 9.,35.,33.,21.,27.,49., 9.,21.,
                       27., 9., 9.,25.,33.,15.,35., 3., 5.,11.,15.,13.,
                       15., 9.,21.,15.,21.,25.,35., 9., 5.,45.,27.,21.,
                       15.,21.,15.,25.,21.,35., 5.,15.,45.,35.,55.,25.};
  static double E[48]={  500., 7500.,12500.,17500.,19000.,19500.,19500.,
			 21000.,22000.,23000.,23000.,24000.,24000.,24500.,
			 24500.,26000.,26500.,26500.,27000.,27500.,28500.,
			 29000.,29500.,29500.,29500.,30000.,31500.,31500.,
			 33500.,33500.,34000.,34500.,34500.,35000.,35500.,
			 37000.,37000.,37000.,38500.,40000.,40000.,41000.,
			 41000.,43000.,43000.,43000.,43000.,44000.};
  static double WNO[48]={63500.,58500.,53500.,59500.,45000.,44500.,44500.,
                         43000.,58000.,41000.,54000.,40000.,40000.,57500.,
                         55500.,38000.,57500.,57500.,37000.,54500.,53500.,
                         55000.,34500.,34500.,34500.,34000.,32500.,32500.,
                         32500.,32500.,32000.,29500.,29500.,31000.,30500.,
                         29000.,27000.,54000.,27500.,24000.,47000.,23000.,
                         44000.,42000.,42000.,21000.,42000.,42000.};
  
  WAVENO=FREQ/2.99792458E10;
  FE10P=0.0;
  if(WAVENO < 21000.) return FE10P;
  for(I=0; I<48; I++) BOLT[I]=G[I]*exp(-E[I]*2.99792458e10*HKT);
  for(I=0; I<48; I++)
  {
    XXX=((WNO[I]+3000.-WAVENO)/WNO[I]/.1);
    XSECT[I]=(WNO[I]<WAVENO)?3.e-18/(1.+XXX*XXX*XXX*XXX):0.;
  }
  for(I=0; I<48; I++) FE10P+=XSECT[I]*BOLT[I];
  
  return FE10P;
}

/* ------------------------------------------------------------------------------ */

void COOLOP(double &ACOOL, double XC1, double XMg1, double XAl1, double XSi1, double XFe1, double STIM,
	    double FREQ, double FREQLG, double T, double TLOG, double TKEV, double HKT)
{
  //     Si I, Mg I, Al I, C I, Fe I
    
  ACOOL=( C1OP(FREQ,TKEV          )*XC1 +
	  Mg1OP(FREQ,FREQLG,T,TLOG)*XMg1+
	  Al1OP(FREQ              )*XAl1+
	  Si1OP(FREQ,FREQLG,T,TLOG)*XSi1+
	  Fe1OP(FREQ,HKT          )*XFe1)*STIM;
  
  return;
}

/* ------------------------------------------------------------------------------ */

double N1OP(double FREQ, double TKEV)
{
  // CROSS-SECTION TIMES PARTITION FUNCTION

  double C1130,C1020,X1130,X1020,X853;
  
  C1130=6.*exp(-3.575/TKEV);
  C1020=10.*exp(-2.384/TKEV);
  X1130=0.;
  X1020=0.;
  X853=0.;
  
  if(FREQ >= 3.517915E15) X853 =SEATON(3.517915E15,1.142E-17,2.0E0,4.29E0,FREQ);
  if(FREQ >= 2.941534E15) X1020=SEATON(2.941534E15,4.410E-18,1.5E0,3.85E0,FREQ);
  if(FREQ >= 2.653317E15) X1130=SEATON(2.653317E15,4.200E-18,1.5E0,4.34E0,FREQ);
  
  return X853*4.+X1020*C1020+X1130*C1130;
}

/* ------------------------------------------------------------------------------ */

double O1OP(double FREQ)
{
  //     CROSS-SECTION TIMES PARTITION FUNCTION
  
  if(FREQ >= 3.28805E15) return 9.*SEATON(3.28805E15,2.94E-18,1.E0,2.66E0,FREQ);
  else return 0.0;
}

/* ------------------------------------------------------------------------------ */

double Mg2OP(double FREQ, double TKEV)
{
  //     CROSS-SECTION TIMES PARTITION FUNCTION
  
  double  C1169,X1169,X824;
  
  C1169=6.*exp(-4.43/TKEV);
  X1169=0.0;
  X824=0.0;
  
  if(FREQ >= 3.635492E15) X824 =SEATON(3.635492E15,1.40E-19,4.E0,6.7E0,FREQ);
  if(FREQ >= 2.564306E15) X1169=5.11E-19*pow(2.564306E15/FREQ,3);

  return X824*2.+X1169*C1169;
}

/* ------------------------------------------------------------------------------ */

double Si2OP(double FREQ, double FREQLG, double T, double TLOG)
{
  //     CROSS-SECTION TIMES THE PARTITION FUNCTION
  
  int NT,N;
  double DT,D,D1,XWL1,XWL2;
  static double PEACH[14][6]=
  /*    10000     12000     14000     16000     18000     20000       WAVE(A) */
    {{-43.8941, -43.8941, -43.8941, -43.8941, -43.8941, -43.8941},/*    500 */
     {-42.2444, -42.2444, -42.2444, -42.2444, -42.2444, -42.2444},/*    600 */
     {-40.6054, -40.6054, -40.6054, -40.6054, -40.6054, -40.6054},/*    759 */
     {-54.2389, -52.2906, -50.8799, -49.8033, -48.9485, -48.2490},/*    760 */
     {-50.4108, -48.4892, -47.1090, -46.0672, -45.2510, -44.5933},/*   1905 */
     {-52.0936, -50.0741, -48.5999, -47.4676, -46.5649, -45.8246},/*   1906 */
     {-51.9548, -49.9371, -48.4647, -47.3340, -46.4333, -45.6947},/*   1975 */
     {-54.2407, -51.7319, -49.9178, -48.5395, -47.4529, -46.5709},/*   1976 */
     {-52.7355, -50.2218, -48.4059, -47.0267, -45.9402, -45.0592},/*   3245 */
     {-53.5387, -50.9189, -49.0200, -47.5750, -46.4341, -45.5082},/*   3246 */
     {-53.2417, -50.6234, -48.7252, -47.2810, -46.1410, -45.2153},/*   3576 */
     {-53.5097, -50.8535, -48.9263, -47.4586, -46.2994, -45.3581},/*   3577 */
     {-54.0561, -51.2365, -49.1980, -47.6497, -46.4302, -45.4414},/*   3900 */
     {-53.8469, -51.0256, -48.9860, -47.4368, -46.2162, -45.2266}};/*  4200 */
  static double FREQSI[7]={4.9965417e15,3.9466738e15,1.5736321e15,
                           1.5171539e15,9.2378947e14,8.3825004e14,
                           7.6869872e14};
  /*     2P,2D,2P,2D,2P */
  static   double FLOG[9]={36.32984,36.14752,35.91165,34.99216,34.95561,
			   34.45941,34.36234,34.27572,34.20161};
  static double TLG[6]={9.21034,9.39266,9.54681,9.68034,9.79813,9.90349};

  
  NT=min(5,(int)floor(T/2000.)-4); if(NT<1) NT=1;
  DT=(TLOG-TLG[NT-1])/(TLG[NT]-TLG[NT-1]);

  for(N=0; N<7; N++) if(FREQ>FREQSI[N]) break;
  D=(FREQLG-FLOG[N])/(FLOG[N+1]-FLOG[N]);
  
  if(N>1) N=2*N-2;
  if(N==13) N=12;
  
  D1=1.-D;
  XWL1=PEACH[N+1][NT-1]*D+PEACH[N][NT-1]*D1;
  XWL2=PEACH[N+1][NT  ]*D+PEACH[N][NT  ]*D1;
  
  return exp(XWL1*(1.-DT)+XWL2*DT)*6.;
}

/* ------------------------------------------------------------------------------ */

double Ca2OP(double FREQ, double TKEV)
{
  //     CROSS-SECTION TIMES THE PARTITION FUNCTION
  double C1218,C1420,X1420,X1218,X1044,XXX;

  C1218=10.*exp(-1.697/TKEV);
  C1420=6.*exp(-3.142/TKEV);
  X1044=0.; X1218=0.; X1420=0.;
  
  if(FREQ>=2.870454e15)
  {
    XXX=(2.870454e15/FREQ); XXX=XXX*XXX*XXX; X1044=1.08e-19*XXX;
  }
  if(FREQ>=2.460127e15) X1218=1.64e-17*sqrt(2.460127e15/FREQ);
  if(FREQ>=2.110779e15) X1420=SEATON(2.110779e15,4.13e-18,3.,0.69, FREQ);
  
  return X1044+X1218*C1218+X1420*C1420;
}
  
/* ------------------------------------------------------------------------------ */

void LUKEOP(double &ALUKE,double XN1, double XO1, double XMg2, double XSi2, double XCa2,
	    double STIM, double FREQ, double FREQLG, double T, double TLOG, double TKEV)
{
  //     N I, O I, Si II, Mg II, Ca II

  ALUKE = ( N1OP(FREQ,TKEV          )*XN1 +
	    O1OP(FREQ               )*XO1 +
	    Mg2OP(FREQ,TKEV         )*XMg2+
	    Si2OP(FREQ,FREQLG,T,TLOG)*XSi2+
	    Ca2OP(FREQ,TKEV         )*XCa2)*STIM;
  
  return;
}

/* ------------------------------------------------------------------------------ */

void HOTOP(double &AHOT){AHOT=0.0;}; // dummy;

/* ------------------------------------------------------------------------------ */

void ELECOP(double &SIGEL, double XNE){ SIGEL = 0.6653E-24*XNE;};

/* ------------------------------------------------------------------------------ */

void H2RAOP(double &SIGH2, double XH1, double FREQ, double T, double TKEV, double TLOG)
{
  double ARG,H1,SIG,WW;
  
  WW=pow(2.997925E18/min(FREQ,2.922E15),2);
  
  SIG=(8.14E-13+1.28e-6/WW+1.61e0/(WW*WW))/(WW*WW);
  ARG=4.477/TKEV-4.6628E1+(1.8031E-3+(-5.023E-7+(8.1424E-11-5.0501E-15*T)*T)*T)*T-1.5*TLOG;
  H1=XH1*2.0;
  
  if(ARG > -80.0) SIGH2=exp(ARG)*H1*H1*SIG;
  else SIGH2=0.0;
  
  return;
}

/* ------------------------------------------------------------------------------ */

void cop(double T, double TKEV, double TK, double HKT, double TLOG,
	 double XNA, double XNE, double *WLGRID, double *OPACITY,
	 double *SCATTER,  
	 double H1, double H2, double HMIN, double HE1, double HE2,
	 double HE3, double C1, double AL1, double SI1, double SI2,
	 double CA1, double CA2, double MG1, double MG2, double FE1,
	 double N1, double O1, int nWLGRID, int NLINES, int NTOTALLIST)
{

  double AHYD,AH2P,AHMIN,AHE1,AHE2,AHEMIN,ACOOL,ALUKE,AHOT,A,B;
  double SIGH,SIGHE,SIGEL,SIGH2;

  for(int iWL = 0; iWL<nWLGRID; iWL++){

    double FREQ=2.997925E18/WLGRID[iWL];
    double FREQLG=log(FREQ);
    double FREQ15=FREQ*1.E-15;
    double EHVKT=exp(-FREQ*HKT);
    double STIM=1.0-EHVKT;
    
    AH2P=0.0, AHE1=0.0, AHE2=0.0, AHEMIN=0.0, ACOOL=0.0, ALUKE=0.0, AHOT=0.0;
    AHYD=0.0, AHMIN=0.0, SIGH=0.0, SIGHE=0.0, SIGEL=0.0, SIGH2=0.0;
    
    HOP(AHYD,XNE,H1,H2,FREQ,FREQLG,T,TLOG,TKEV,STIM,EHVKT);
    H2PLOP(AH2P,H1,H2,FREQ,FREQLG,FREQ15,TKEV,STIM);
    HMINOP(AHMIN,H1,HMIN,FREQ,T,TKEV,XNE,EHVKT);
    HRAYOP(SIGH,H1,FREQ);
    HE1OP(AHE1,HE1,HE2,XNE,FREQ,FREQLG,T,TKEV,TLOG,EHVKT,STIM);
    HE2OP(AHE2,HE2, HE3,XNE,FREQ,FREQLG,T, TKEV, TLOG, EHVKT,STIM);
    HEMIOP(AHEMIN,HE1,FREQ,T,XNE);
    HERAOP(SIGHE,HE1,FREQ);
    
    if(T < 12000.) COOLOP(ACOOL,C1 ,MG1, AL1,SI1, FE1,STIM,FREQ,FREQLG,T,TLOG,TKEV,HKT);
    if(T < 30000.) LUKEOP(ALUKE, N1, O1, MG2, SI2, CA2,STIM,FREQ,FREQLG,T,TLOG,TKEV);
    
    HOTOP(AHOT);
    ELECOP(SIGEL,XNE);
    H2RAOP(SIGH2,H1,FREQ,T,TKEV,TLOG);

    A=AHYD+AHMIN+AH2P+AHE1+AHE2+AHEMIN+ACOOL+ALUKE+AHOT;
    B=SIGH+SIGHE+SIGEL+SIGH2;
    
    OPACITY[iWL]=A+B;
    SCATTER[iWL]=B; 
  }

}

/* ------------------------------------------------------------------------------ */
