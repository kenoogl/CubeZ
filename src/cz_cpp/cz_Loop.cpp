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

// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 */
int CZ::Loop()
{
  double  res=0.0;
  int itr=0;

  switch (ls_type)
  {
    case LS_JACOBI:
      if ( 0 == (itr=JACOBI(res, P, RHS)) ) return 0;
      break;

    case LS_PSOR:
      if ( 0 == (itr=PSOR(res, P, RHS)) ) return 0;
      break;

    case LS_SOR2SMA:
      if ( 0 == (itr=RBSOR(res,P, RHS)) ) return 0;
      break;

    case LS_BICGSTAB:
      //ILAP = PBiCGstab(res_p, P, RHS);
      break;

    default:
      break;
  }

  Hostonly_ {
    printf("\n=================================\n");
    printf("Iter = %d  Res = %e\n", itr, res);
    printf("=================================\n");
  }


  return 1;
}
