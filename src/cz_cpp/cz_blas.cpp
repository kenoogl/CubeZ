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

double CZ::dot1(const int* sz, const int g, const int * idx,
                const REAL_TYPE* var, double& flop)
{
  __assume_aligned(d, ALIGN);
  double r = 0.0;

  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;
  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  double tmp;

  #pragma loop count (SdW-GUIDE-2)
  #pragma ivdep
  for (int k=1; k<nx; k++)
  {
    tmp = (double)(var[i] * msk[i]);
    r += tmp * tmp;
  }

  #pragma ivdep
  for (int i=bst; i<bed; i++)
  {
    tmp = (double)(var[i] * msk[i]);
    r += tmp * tmp;
  }

  #pragma loop count (SdW-GUIDE-2)
  #pragma ivdep
  for (int i=bed; i<nx; i++)
  {
    tmp = (double)(var[i] * msk[i]);
    r += tmp * tmp;
  }

  return r;
}

double CZ::dot1(const int* sz, const int g, const int * idx,
                const REAL_TYPE* var, double& flop)
{
  __assume_aligned(d, ALIGN);
  double r = 0.0;

  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;
  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  double tmp;

  #pragma loop count (SdW-GUIDE-2)
  #pragma ivdep
  for (int k=1; k<bst; k++)
  {
    tmp = (double)(var[i] * msk[i]);
    r += tmp * tmp;
  }

  #pragma ivdep
  for (int i=bst; i<bed; i++)
  {
    tmp = (double)(var[i] * msk[i]);
    r += tmp * tmp;
  }

  #pragma loop count (SdW-GUIDE-2)
  #pragma ivdep
  for (int i=bed; i<nx; i++)
  {
    tmp = (double)(var[i] * msk[i]);
    r += tmp * tmp;
  }

  return r;
}
