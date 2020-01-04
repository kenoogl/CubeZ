/*
###################################################################################
#
# CubeZ
#
# Copyright (C) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

// @file test1.cpp
// @note TDMA test

#include "cz.h"

#define N 3

int main(int argc, char *argv[])
{
  REAL_TYPE a[N], b[N], c[N], d[N];
  REAL_TYPE w[N];
  int n = N;

 b[0] = 2; c[0] = 3;
 a[1] = 4; b[1] = 4; c[1] = -3;
           a[2] = 3; b[2] = -1;

 d[0] = 8.0;
 d[1] = 3.0;
 d[2] = 3.0;

 //tdma_0_(&n, d, a, b, c, w);
 CZ z;
 z.tdma(n, d, a, b, c, w);

 printf("%f %f %f\n", b[0], c[0], 0.0);
 printf("%f %f %f\n", a[1], b[1], c[1]);
 printf("%f %f %f\n", 0.0,  a[2], b[2]);
 printf("\n");
 printf("%f %f %f\n", d[0], d[1], d[2]);

 return 0;
}
