/* ------- file: -------------------------- matrix.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Dec 29 14:15:27 1999 --

       --------------------------                      ----------RH-- */

/* --- Routines to create 2-dimensional arrays that are contiguous in
       memory. Therefore the whole array pointed to by pointer **matrix_??
       can be addressed as either 1-dimensional array with pointer
       matrix_??[0], or as 2-dimensional array as matrix_??[i][j].

 Note: space is initialized to zeros by using calloc rather than malloc.

       --                                              -------------- */
 
#include <stdlib.h>

#include "rh.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern char   messageStr[];


/* ------- begin -------------------------- matrix_char.c ----------- */

char **matrix_char(int Nrow, int Ncol)
{
  register int i;

  char *theMatrix, **Matrix;
  int   typeSize = sizeof(char), pointerSize = sizeof(char *);

  theMatrix = (char *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (char **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_char.c ----------- */

/* ------- begin -------------------------- matrix_int.c ------------ */

int **matrix_int(int Nrow, int Ncol)
{
  register int i;
  static const int     typeSize = sizeof(double), pointerSize = sizeof(double *);

  int *theMatrix, **Matrix;

  theMatrix = (int *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (int **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_int.c ------------ */

/* ------- begin -------------------------- matrix_double.c --------- */

double **matrix_double(int Nrow, int Ncol)
{
  register int i;

  static const int     typeSize = sizeof(double), pointerSize = sizeof(double *);
  double *theMatrix, **Matrix;

  theMatrix = (double *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (double **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_double.c --------- */
float **matrix_float(int Nrow, int Ncol)
{
  register int i;

  static const int     typeSize = sizeof(float), pointerSize = sizeof(float *);
  float *theMatrix, **Matrix;

  theMatrix = (float *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (float **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- begin -------------------------- freeMatrix.c ------------ */

void freeMatrix(void **matrix)
{
  const char routineName[] = "freeMatrix";

  if (matrix == NULL  ||  matrix[0] == NULL) {
    sprintf(messageStr, "Trying to free NULL pointer");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    /* --- Free the memory allocated for matrix matrix -- ----------- */

    free(matrix[0]);
    free(matrix);
  }
}
/* ------- end ---------------------------- freeMatrix.c ------------ */

short **matrix_short(int Nrow, int Ncol)
{
  register int i;
  static const unsigned  typeSize = sizeof(short),  pointerSize = sizeof(short *);
 
  short *theMatrix, **Matrix;

  theMatrix = (short *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (short **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- freeMatrix.c ------------ */

unsigned char **matrix_uchar(int Nrow, int Ncol)
{
  register int i;
  static const unsigned  typeSize = sizeof(short),  pointerSize = sizeof(short *);
  
  unsigned char *theMatrix, **Matrix;

  theMatrix = (unsigned char *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (unsigned char **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
