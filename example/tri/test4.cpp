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

// @file test4.cpp
// @note PCR test with Dirichlet BC

#include "cz.h"

/*


                     D1                         D2
                     |                          |
                     |                          |
    |--+--+-...-+--+-|--+--+--+--......+--+--+--|--+-...-+--+--|
  -15                0  1  2  3  ......  22 23 24              39
                        <-------------------->
                             Inner points
  BC
    <Dirichlet> ==> 係数は境界点方向はゼロで、ソース項にマイナス値
    D1 = -3.0
    D2 =  9.0
 */

//  内点の数
#define N 23

/*
  N              // 解くべき方程式の次元数
  2^5=32 > nx    // nxを超える最小の2べき数
  pn=5           // その指数
  ss=2^{pn-1}=16 // PCRで参照するストライド s の最大値
  1-ss=-15       // 配列参照の下限
  nx+ss=39       // 配列参照の上限
  nx+2*ss=55     // 配列長

 */

 REAL_TYPE* allocReal(const int sz, const int pn)
 {
   if ( !sz ) return NULL;

   int ss = 0x1 << (pn-1); // sz を超えない最大の2べき数
   size_t nx = sz + 2*ss;

   //printf("%d %d %d\n", sz, pn, ss);

   REAL_TYPE* var = new REAL_TYPE[nx];

 // #pragma omp parallel for
   for (int i=0; i<nx; i++) var[i]=0.0;

   return var;
 }

int main(int argc, char *argv[])
{
  //                     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
  //REAL_TYPE a[N+2] = { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};
  //REAL_TYPE b[N+2] = { 0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2, 0};
  //REAL_TYPE c[N+2] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0};
  //REAL_TYPE d[N+2] = {-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-9, 9};
  // d[]はRHSと解ベクトルに利用、置換は端点以外
  // 端点には境界条件、それ以外はソース項

  REAL_TYPE *a=NULL, *b=NULL, *c=NULL, *d=NULL, *p=NULL, *q=NULL, *r=NULL;
  int n = N;
  int pn;

  CZ z;

  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=z.getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  const int ss = 0x1 << (pn-1);

  //printf("n=%d : pn=%d , ss=%d \n", n, pn, ss);



  // array
  a = allocReal(n, pn);
  b = allocReal(n, pn);
  c = allocReal(n, pn);
  d = allocReal(n, pn);
  p = allocReal(n, pn);
  q = allocReal(n, pn);
  r = allocReal(n, pn);

  for (int i=2; i<=n; i++) a[_IDX_(i,ss)]=1.0;
  for (int i=1; i<=n; i++) b[_IDX_(i,ss)]=-2.0;
  for (int i=1; i<=n-1; i++) c[_IDX_(i,ss)]=1.0;
  d[_IDX_(0,ss)] = -3.0;
  d[_IDX_(1,ss)] =  3.0;
  d[_IDX_(n,ss)] = -9.0;
  d[_IDX_(n+1,ss)] =  9.0;

  // b[]=1.0となるように正規化
  for (int i=1; i<=n; i++)
  {
    a[_IDX_(i,ss)] /= -2.0;
    b[_IDX_(i,ss)] /= -2.0;
    c[_IDX_(i,ss)] /= -2.0;
    d[_IDX_(i,ss)] /= -2.0;
  }

  // copy
  for (int i=0; i<n+2*ss-1; i++) {
    p[i] = d[i];
    q[i] = a[i];
    r[i] = c[i];
  }

  z.printA(n, ss, a, "A");
  z.printA(n, ss, c, "C");
  printf("\n");

  // 内点の個数と、その先頭アドレス、外点の両端は境界値
  z.pcr(n, pn, d, a, c, p, q, r);

  // ok
  //z.tdma(n, &d[_IDX_(1,ss)], &a[_IDX_(1,ss)], &b[_IDX_(1,ss)], &c[_IDX_(1,ss)], &p[_IDX_(1,ss)]);

  FILE* fp;
  fp=fopen("result.txt", "w");

  for (int i=1-ss; i<=n+ss; i++) {
    fprintf(fp,"%2d, %6.3f\n", i, d[_IDX_(i,ss)]);
  }

  fclose(fp);

  return 0;
}
