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
                                                                           
                                                                       
  if(myrank == 0){
    const char title[16][72] =    {
      "   SSSSSSSSSSSSSSS TTTTTTTTTTTTTTTTTTTTTTTIIIIIIIIII      CCCCCCCCCCCCC",
      " SS:::::::::::::::ST:::::::::::::::::::::TI::::::::I   CCC::::::::::::C",
      "S:::::SSSSSS::::::ST:::::::::::::::::::::TI::::::::I CC:::::::::::::::C",
      "S:::::S     SSSSSSST:::::TT:::::::TT:::::TII::::::IIC:::::CCCCCCCC::::C",
      "S:::::S            TTTTTT  T:::::T  TTTTTT  I::::I C:::::C       CCCCCC",
      "S:::::S                    T:::::T          I::::IC:::::C              ",
      " S::::SSSS                 T:::::T          I::::IC:::::C              ",
      "  SS::::::SSSSS            T:::::T          I::::IC:::::C              ",
      "    SSS::::::::SS          T:::::T          I::::IC:::::C              ",
      "       SSSSSS::::S         T:::::T          I::::IC:::::C              ",
      "            S:::::S        T:::::T          I::::IC:::::C              ",
      "            S:::::S        T:::::T          I::::I C:::::C       CCCCCC",
      "SSSSSSS     S:::::S      TT:::::::TT      II::::::IIC:::::CCCCCCCC::::C",
      "S::::::SSSSSS:::::S      T:::::::::T      I::::::::I CC:::::::::::::::C",
      "S:::::::::::::::SS       T:::::::::T      I::::::::I   CCC::::::::::::C",
      " SSSSSSSSSSSSSSS         TTTTTTTTTTT      IIIIIIIIII      CCCCCCCCCCCCC"};
    
                                                                  
    
    for(int ii = 0; ii<16; ii++) fprintf(stderr,"%s\n", title[ii]);
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
