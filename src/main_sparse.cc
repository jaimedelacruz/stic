#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <string>
#include <sys/time.h>

#include "master_sparse.h"
#include "slave.h"
using namespace std;
//
double master_getTime(double t0 = -1.0){
  struct timeval dum;
  gettimeofday(&dum, NULL);
  if(t0 < 0.0) return dum.tv_sec + dum.tv_usec * 1.0E-6;
  else return (dum.tv_sec + dum.tv_usec * 1.0E-6) - t0;
}
//
std::string double2time(double dt)
{
  int hh = int(dt/3600.);
  dt -= hh*3600.;
  int mm = int(dt/60.);
  dt -= mm*60.;
  
  char buff[13];
  sprintf(buff, "%02d:%02d:%06.3f", hh, mm, dt);

  return std::string(buff);
}

//
int main(int narg, char* argv[]){
  
                                                      
  double dt=master_getTime();

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
    const char title[16][71] =
      {                                                  
	"   SSSSSSSSSSSSSSS TTTTTTTTTTTTTTTTTTTTTTT  iiii         CCCCCCCCCCCCC",
	" SS:::::::::::::::ST:::::::::::::::::::::T i::::i     CCC::::::::::::C",
	"S:::::SSSSSS::::::ST:::::::::::::::::::::T  iiii    CC:::::::::::::::C",
	"S:::::S     SSSSSSST:::::TT:::::::TT:::::T         C:::::CCCCCCCC::::C",
	"S:::::S            TTTTTT  T:::::T  TTTTTTiiiiiii C:::::C       CCCCCC",
	"S:::::S                    T:::::T        i:::::iC:::::C              ",
	" S::::SSSS                 T:::::T         i::::iC:::::C              ",
	"  SS::::::SSSSS            T:::::T         i::::iC:::::C              ",
	"    SSS::::::::SS          T:::::T         i::::iC:::::C              ",
	"       SSSSSS::::S         T:::::T         i::::iC:::::C              ",
	"            S:::::S        T:::::T         i::::iC:::::C              ",
	"            S:::::S        T:::::T         i::::i C:::::C       CCCCCC",
	"SSSSSSS     S:::::S      TT:::::::TT      i::::::i C:::::CCCCCCCC::::C",
	"S::::::SSSSSS:::::S      T:::::::::T      i::::::i  CC:::::::::::::::C",
	"S:::::::::::::::SS       T:::::::::T      i::::::i    CCC::::::::::::C",
	" SSSSSSSSSSSSSSS         TTTTTTTTTTT      iiiiiiii       CCCCCCCCCCCCC"
      };
    
    cerr<<endl;
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
  if(myrank == 0){
    string ti = double2time(master_getTime(dt));
    printf("main: total elapsed time [%s]\n", ti.c_str());
  }
  
  MPI_Finalize();
  
}
