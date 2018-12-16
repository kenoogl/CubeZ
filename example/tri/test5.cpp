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
    |--+--+--+--......+--+--+--|
    0  1  2  3  ......  22 23 24
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
 */

 REAL_TYPE* allocReal(const int sz)
 {
   if ( !sz ) return NULL;
   size_t nx = sz+2;

   REAL_TYPE* var = new REAL_TYPE[nx];

   #pragma omp parallel for
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

  REAL_TYPE *a=NULL, *b=NULL, *c=NULL, *d=NULL, *p=NULL;
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
/*
  int s=8;
  for (int i=0; i<=n; i++)
  {
    int iL = i - s;
    iL = std::max(iL,0);
    printf("%2d : %d\n", i, iL);
  }
*/

  // array
  a = allocReal(n);
  b = allocReal(n);
  c = allocReal(n);
  d = allocReal(n);
  p = allocReal(n);

  for (int i=2; i<=n; i++) a[i]=1.0;
  for (int i=1; i<=n; i++) b[i]=-2.0;
  for (int i=1; i<=n-1; i++) c[i]=1.0;
  d[0] = -3.0;
  d[1] =  3.0;
  d[n] = -9.0;
  d[n+1] =  9.0;

  // b[]=1.0となるように正規化

  for (int i=1; i<=n; i++)
  {
    a[i] /= -2.0;
    b[i] /= -2.0;
    c[i] /= -2.0;
    d[i] /= -2.0;
  }

  z.printB(n, a, "A");
  z.printB(n, c, "C");
  printf("\n");

  // 内点の個数と、その先頭アドレス、外点の両端は境界値
  z.pcr2(n, pn, d, a, c);

  // ok
  //z.tdma(n, &d[1], &a[1], &b[1], &c[1], &p[1]);

  FILE* fp;
  fp=fopen("result.txt", "w");

  for (int i=0; i<=n+1; i++) {
    fprintf(fp,"%2d, %6.3f\n", i, d[i]);
  }

  fclose(fp);

  return 0;
}
