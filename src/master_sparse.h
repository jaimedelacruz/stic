#ifndef MASTER_SPARSE_H
#define MASTER_SPARSE_H
//
#include "input.h"
#include "depthmodel.h"
#include "cmemt.h"
//
void do_master_sparse(int myrank, int nprocs,  char hostname[]);
void slaveInversion(iput_t &input, mdepthall_t &m, mat<double> &obs, mat<double> &pars, mat<double> &chi2);

#endif
