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

// @file test2.cpp
// @note TDMA test with N/D bc and mask

#include "cz.h"

/*
    D1                N    D2          D3
    |                 |     |           |
    |                 |solid|           |
    |--+--+--+--+--+--+--+--+--+--+--+--|
    0  1  2  3  4  5  6  7  8  9 10 11 12

   -1  ..............-1  1  .......    2

  BC
    <Dirichlet> ==> 係数は境界点方向はゼロで、ソース項にマイナス値
    D1 = -1.0
    D2 =  1.0
    D3 =  2.0

    <Neumann> ==> 係数がゼロ
    N
 */

//  内点の数
#define N 11

int main(int argc, char *argv[])
{
  //                   0  1  2  3  4  5  6  7  8  9 10 11 12
  REAL_TYPE a[N+2] = { 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0};
  REAL_TYPE b[N+2] = { 0,-2,-2,-2,-2,-2,-1,-2,-2,-2,-2,-2, 0};
  REAL_TYPE c[N+2] = { 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0};
  REAL_TYPE d[N+2] = {-1, 1, 0, 0, 0, 0, 0, 0,-1, 0, 0,-2, 2};
  // d[]はRHSと解ベクトルに利用、置換は端点以外
  // 端点には境界条件、それ以外はソース項

  // マスク ...> 不要
  REAL_TYPE m[N+2] = { 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0};


  REAL_TYPE w[N+2], g[N+2];
  int n = N;

  // 内点の個数と、その先頭アドレス
  tdma_0_(&n, &d[1], &a[1], &b[1], &c[1], &w[1]);

  for (int i=0; i<n+2; i++) {
    printf("[%2d] %6.3f\n", i, d[i]);
  }

  return 0;
}
