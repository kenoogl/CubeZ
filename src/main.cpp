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

// デバッグ時に定義
#define __DEBUG

int main(int argc, char *argv[]) {

  int myRank=0;

  if (argc != 6 && argc != 9) {
    if ( myRank == 0) {
      printf("\tUsage : ./cz gsz_x, gsz_y, gsz_z, linear_solver, IterationMax [gdv_x, gdv_y, gdv_z]\n");
      printf("\t$ ./cz 64 64 64 pbicgstab 1000 2 2 2\n");
      printf("\t$ ./cz 64 64 64 pbicgstab 1000\n");
      printf("\t\tlinear_solver = {jacobi | sor | sor2sma | pbicgstab | lsor}\n");
    }
    return 0;
  }


#ifndef DISABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  int mode = 0;

#ifdef __DEBUG
  mode = 1;
#endif

  CZ cz;
  cz.debug(mode);


  if( 0==cz.Init(argc, argv) )
  {
    if ( myRank == 0) printf("\n\tSolver initialize error.\n\n");
    #ifndef DISABLE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }


  if ( 0==cz.Loop() )
  {
    if ( myRank == 0) printf("\n\tSolver MainLoop error.\n");
    #ifndef DISABLE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }


  if( 0==cz.Post() )
  {
    if ( myRank == 0) printf("\n\tSolver post error.\n");
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
