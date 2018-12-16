/*
###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

#include "cz.h"


int main(int argc, char *argv[]) {

  int myRank=0;

  if (argc != 6 && argc != 7 && argc != 9 && argc != 10) {
    if ( myRank == 0) {
      printf("\tUsage : ./cz-mpi gsz_x, gsz_y, gsz_z, linear_solver, IterationMax, [precond] [gdv_x, gdv_y, gdv_z]\n");
      printf("\t\tlinear_solver = {jacobi | psor | sor2sma | pbicgstab | lsor | ljacobi}\n");
      printf("\t\tprecond = {none | jacobi | psor | sor2sma}\n\n");
      printf("\t$ ./cz-mpi 64 64 64 jacobi 4000 2 2 1\n");
      printf("\t$ ./cz-mpi 64 64 64 psor 4000\n");
      printf("\t$ ./cz-mpi 64 64 64 pbicgstab 4000 sor2sma\n");
      printf("\t$ ./cz-mpi 64 64 64 pbicgstab 4000 sor2sma 2 1 3\n");
    }
    return 0;
  }


#ifndef DISABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  int mode = 0;
  mode = 1; // デバッグ時に定義

  CZ cz;
  cz.debug(mode);


  if( 0==cz.Evaluate(argc, argv) )
  {
    if ( myRank == 0) printf("\n\tSolver error.\n\n");
    #ifndef DISABLE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }


  #ifndef DISABLE_MPI
  MPI_Finalize();
  #endif

  return 0;
}
