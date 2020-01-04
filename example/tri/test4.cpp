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

// @file test4.cpp
// @note multi-sysytem LU

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

// システムの数
#define Msystem 32

/*
 N              // 解くべき方程式の次元数
 */


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
  
  REAL_TYPE *a=NULL, *c=NULL, *d=NULL, *w=NULL, *b=NULL;
  int n = N;
  int numThreads=0;
  
  CZ z;
  double flop = 0.0;
  
  
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
    
    set_label("TDMAcoupled",  PerfMonitor::CALC);
  }
  
  
  // array
  int ms = Msystem;
  
  REAL_TYPE var_type=0;
  int nx = (n+2) * ms;
  
  a  = z.czAllocR(nx, var_type);
  b  = z.czAllocR(nx, var_type);
  c  = z.czAllocR(nx, var_type);
  d  = z.czAllocR(nx, var_type);
  w  = z.czAllocR(nx, var_type);
  
  
  for (int i=2; i<=n; i++) {
    for (int m=0; m<ms; m++) {
      a[m+ms*i]=1.0;
    }
  }
  for (int i=1; i<=n; i++) {
    for (int m=0; m<ms; m++) {
      b[m+ms*i]=-2.0;
    }
  }
  for (int i=1; i<=n-1; i++) {
    for (int m=0; m<ms; m++) {
      c[m+ms*i]=1.0;
    }
  }
  
  for (int m=0; m<ms; m++) {
    d[m+ms*0] = -3.0;
    d[m+ms*1] =  3.0;
    d[m+ms*n] = -9.0;
    d[m+ms*(n+1)] = 9.0;
  }
  
  for (int i=1; i<=n; i++)
  {
    for (int m=0; m<ms; m++) {
      a[m+ms*i] /= -2.0;
      b[m+ms*i] /= -2.0;
      c[m+ms*i] /= -2.0;
      d[m+ms*i] /= -2.0;
    }
  }
  
  flop = (double)(16.0*(double)(n-1)*ms);
  
  
  // 内点の個数と、その先頭アドレス、外点の両端は境界値
  TIMING_start("TDMAcoupled");
  //z.tdma1(n, &d[1], &a[1], &c[1], &w[1]);
  //tdma_p_(&n, &d[1], &a[1], &c[1], &w[1]);
  tdma_mp_(&n, &ms, d, a, c, w);
  TIMING_stop("TDMAcoupled", flop);
  
  // ok
  //z.tdma(n, &d[1], &a[1], &b[1], &c[1], &w[1]);
  
  char str[100];
  sprintf(str, "TDMA coupled");
  PM.print(stdout, "hoge", str);
  
  
  FILE* fp;
  fp=fopen("result.txt", "w");
  
  for (int i=0; i<=n+1; i++) {
    fprintf(fp,"%2d, ", i);
    for (int m=0; m<Msystem; m++) {
      fprintf(fp,"%6.2f, ", d[m+Msystem*i]);
    }
    fprintf(fp,"\n");
  }
  
  fclose(fp);
  
  z.czDelete(w);
  z.czDelete(a);
  z.czDelete(b);
  z.czDelete(c);
  z.czDelete(d);
  
  return 0;
}
