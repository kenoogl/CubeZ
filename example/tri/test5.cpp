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

// @file test5.cpp
// @note PCR test with Dirichlet BC

#include "cz.h"

int order_of_PM_key=0;

pm_lib::PerfMonitor PM;

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

 void set_label(const std::string label, pm_lib::PerfMonitor::Type type)
 {
   bool exclusive=true;

   // 登録個数のチェック
   order_of_PM_key++;

   if ( order_of_PM_key > PM_NUM_MAX )
   {
     printf("\tThe number of labels for Performance monitor goes over limit.\n");
     exit(0);
   }

   // 文字数がTM_LABEL_MAX-1を超えるものはカット
   if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
   {
     printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
   }

   // Performance Monitorへの登録
   PM.setProperties(label, type, exclusive);
 }

 inline void TIMING_start(const std::string key) {
   PM.start(key);
 }

 inline void TIMING_stop(const std::string key, double flopPerTask=0.0, int iterationCount=1) {
   PM.stop(key, flopPerTask, (unsigned)iterationCount);
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

  REAL_TYPE *a=NULL, *c=NULL, *d=NULL;
  REAL_TYPE *a1=NULL, *c1=NULL, *d1=NULL;
  int n = N;
  int pn;
  int numThreads=0;
  
  CZ z;
  double flop = 0.0;

  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=z.getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  const int ss = 0x1 << (pn-1);
  printf("n=%d, pn=%d, ss=%d\n", n, pn, ss);


  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( 1 );
  std::string Parallel_str;
  if ( numThreads > 1 ) {
    Parallel_str = "OpenMP";
  }
  else {
    Parallel_str = "Serial";
  }
  PM.setParallelMode(Parallel_str, numThreads, 1);

  {
    using namespace pm_lib;

    set_label("PCR",  PerfMonitor::CALC);
  }


  // array
  a = allocReal(n);
  c = allocReal(n);
  d = allocReal(n);
  a1= allocReal(n);
  c1= allocReal(n);
  d1= allocReal(n);

  for (int i=2; i<=n; i++) a[i]=1.0;
  for (int i=1; i<=n-1; i++) c[i]=1.0;
  d[0] = -3.0;
  d[1] =  3.0;
  d[n] = -9.0;
  d[n+1] = 9.0;

  for (int i=1; i<=n; i++)
  {
    a[i] /= -2.0;
    c[i] /= -2.0;
    d[i] /= -2.0;
  }

  // 内点の個数と、その先頭アドレス、外点の両端は境界値
  TIMING_start("PCR");
  z.pcr(n, pn, d, a, c, d1, a1, c1, flop);
  TIMING_stop("PCR", flop);

  // ok
  //z.tdma(n, &d[1], &a[1], &b[1], &c[1], &p[1]);

  char str[100];
  sprintf(str, "PCR test");
  PM.print(stdout, "hoge", str);


  FILE* fp;
  fp=fopen("result.txt", "w");

  for (int i=0; i<=n+1; i++) {
    fprintf(fp,"%2d, %e %e\n", i, d[i], d[i]-(-3.0+12/(float)(n+1)*(float)i) );
  }

  fclose(fp);

  return 0;
}
