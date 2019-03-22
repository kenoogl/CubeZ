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

// @file cz_pcr.cpp


/*
 * @brief Parallel Cyclic Reduction
 * @param [in      nx   配列長
 * @param [in      pn   nxを超える最小の2べき数の指数
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @note i方向に領域分割なしを想定
 *       tdma()とは異なり、解くべき範囲（内点）は[1,nx]

  nx = 23        // 解くべき方程式の次元数
  2^5=32 > nx    // nxを超える最小の2べき数
  pn=5           // その指数
  ss=2^{pn-1}=16 // PCRで参照するストライド s の最大値


 D1                         D2
 |                          |
 |                          |
 |--+--+--+--......+--+--+--|
 0  1  2  3  ......  22 23 24
 　　<-------------------->
　　　　 Inner points

 */

void CZ::pcr(const int nx, const int pn,
             REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c,
             REAL_TYPE* d1, REAL_TYPE* a1, REAL_TYPE* c1, double& flop)
{
  const int ss = 0x1 << (pn-1);
  
  //printArray(nx, a, "a");
  //printArray(nx, c, "c");
  //printArray(nx, d, "d");

  int s=0;
  for (int p=1; p<=pn; p++)
  //for (int p=1; p<pn; p++) // 一つ少なくすると2x2  or 3x3 の直接反転
  {
    s = 0x1 << (p-1); // s=2^{p-1}
    
    pcr_kernel(nx, s, d, a, c, d1, a1, c1, flop);
    
    REAL_TYPE* tmp;
    tmp = a; a = a1; a1 = tmp;
    tmp = c; c = c1; c1 = tmp;
    tmp = d; d = d1; d1 = tmp;
    
    //printArray(nx, a, "a");
    //printArray(nx, c, "c");
    //printArray(nx, d, "d");
  }
}


void CZ::pcr_kernel(const int nx, const int s,
                    REAL_TYPE* d,  REAL_TYPE* a,  REAL_TYPE* c,
                    REAL_TYPE* dn, REAL_TYPE* an, REAL_TYPE* cn,
                    double& flop)
{
  REAL_TYPE r, ap, cp;
  int iL, iR;

  flop += (double)(nx)*21.0;

  // #pragma omp parallel for simd private(iL, iR, ap, cp, r)
 #pragma omp parallel for private(iL, iR, ap, cp, r)
  for (int i=1; i<=nx; i++)
  {
    iL = std::max(i-s,0);
    iR = std::min(i+s,nx+1);

    ap = a[i];
    cp = c[i];
    r = 1.0 / ( 1.0 - ap * c[iL] - cp * a[iR] );
    an[i] = - r * ap * a[iL];
    cn[i] = - r * cp * c[iR];
    dn[i] =   r * ( d[i] - ap * d[iL] - cp * d[iR] );
  }
}

void CZ::printArray(int nx, REAL_TYPE* a, char* s)
{
  printf("%s : ", s);
  for (int i=0; i<=nx+1; i++)
  {
    printf("%6.3f ", a[i] );
  }
  printf("\n");
}
