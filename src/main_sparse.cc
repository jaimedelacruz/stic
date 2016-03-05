#include <mpi.h>
#include <iostream>
#include <cstdio>

#include "master_sparse.h"
#include "slave.h"
using namespace std;
//
int main(int narg, char* argv[]){
  
                                                      


  //
  // Init MPI
  //
  int nprocs = 1, myrank = 0, hlen = 0;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int status = 0;
  MPI_Init(&narg, &argv);


  //
  // Fill-in vars
  //
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);// Job number
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);// number of processors
  MPI_Get_processor_name(&hostname[0], &hlen);// Hostname


  /* --- Printout logo :-) --- */
  
  if(myrank == 0){
    const char title[6][100] = {"                         ______  _   _  _                                        ",                                 
				"                         | ___ \\| | | |(_)                                       ",
				" ______  ______  ______  | |_/ /| |_| | _  _ __  __   __  ______  ______  ______ ",
				"|______||______||______| |    / |  _  || || '_ \\ \\ \\ / / |______||______||______|",
				"                         | |\\ \\ | | | || || | | | \\ V /                          ",
				"                         \\_| \\_|\\_| |_/|_||_| |_|  \\_/                           "};
    
    
    
    for(int ii = 0; ii<6; ii++) fprintf(stderr,"%s\n", title[ii]);
    cerr<<endl;
  }               
  

  
  //
  // Split master/slave work
  //
  if(myrank == 0) do_master_sparse(myrank, nprocs, hostname);
  else do_slave(myrank, nprocs, hostname);

  
  //
  // Finalize work
  //
  MPI_Barrier(MPI_COMM_WORLD); // Wait until all processors reach this point
  MPI_Finalize();
  
}
