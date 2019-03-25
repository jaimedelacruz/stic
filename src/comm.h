/*
  Interface for MPI routines
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
*/
#ifndef COMM_H
#define COMM_H
//
#include <mpi.h>
#include <vector>
//#include <cmath>
#include "cmemt.h"
#include "input.h"
#include "depthmodel.h"
//
inline void comm_get_xy(const int t, const int nx, int &y, int &x){
  //double dum = double(t) / double(nx);
  //y = int(std::floor(dum));
  y = t / nx;
  x = t - (y*nx); 
}
//
void comm_get_buffer_size(iput_t &input);
void comm_send_parameters(iput_t &input);
void comm_recv_parameters(iput_t &input);
void comm_master_pack_data(iput_t &input, mat<double> &obs, mat<double> &model, 
			   unsigned long &ipix, int proc, mdepthall_t &m, int cgrad, int action = 1);
//void comm_master_unpack_data(int &iproc, iput_t input, mat<double> &obs, 
//			     mat<double> &pars, mat<double> &chi2);
void comm_master_unpack_data(int &iproc, iput_t input, mat<double> &obs, 
			     mat<double> &pars, mat<double> &chi2, unsigned long &irec,
			     mat<double> &dobs, int cgrad, mdepthall_t &m);

void comm_slave_unpack_data(iput_t &input, int &action, mat<double> &obs, mat<double> &pars, std::vector<mdepth_t> &m, int &cgrad);
void comm_kill_slaves(iput_t &input, int nprocs);
void comm_slave_pack_data(iput_t &input, mat<double> &obs, mat<double> &pars, mat<double> &dobs, int cgrad, std::vector<mdepth_t> &m); 
void comm_send_weights(iput_t &input, mat<double> &w);
int getNinstrumentData(std::vector<region_t> const &reg);

#endif /* COMM_H */ 
