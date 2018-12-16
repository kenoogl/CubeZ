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

// @file cm_blas.cpp

/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     b    D vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 *       cz_dsolver tdma_0 と同等
 */
void CZ::tdma(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* b, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  d[0] = d[0]/b[0];
  w[0] = c[0]/b[0];

  for (int i=1; i<nx; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }

  for (int i=nx-2; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

}


/*
 * @brief Parallel Cyclic Reduction
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     b    D vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 *       cz_dsolver tdma_0 と同等
 */
void CZ::pcr(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* b, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  d[0] = d[0]/b[0];
  w[0] = c[0]/b[0];

  for (int i=1; i<nx; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }

  for (int i=nx-2; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

  int pn;
  if ( -1 == (pn=getNumStage(nx))) {
    printf("error : number of stage\n");
    exit(0);
  }
  printf("%d : pn=%d\n", nx, pn);

  int s;
  for (int p=1; p<=pn; p++)
  {
    s = 0x1 << (p-1);
    printf("%d %d\n", p, s);
  }

}
