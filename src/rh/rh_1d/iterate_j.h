#ifndef ITER_J_H
#define ITER_J_H

void Iterate_j(int NmaxIter, double iterLimit, double *dpopmax_out);
double solveSpectrum(bool_t eval_operator, bool_t redistribute, int iter);
void *Formal_pthread(void *argument);


#endif
