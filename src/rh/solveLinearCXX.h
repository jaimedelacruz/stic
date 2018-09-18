#ifndef SOLL_H
#define SOLL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <rpc/types.h>
#include <stdbool.h>
  
  void solveLinearCXX(int N, double **A, double *b, bool_t improve);
 
    
#ifdef __cplusplus
}
# endif

#endif
