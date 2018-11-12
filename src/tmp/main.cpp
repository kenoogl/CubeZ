/*
 *  main.cpp
 *  Cube kernel
 *
 * Copyright (C) 2016-2017  Research Institute for Information Technology(RIIT), Kyushu University.
 * All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "kernel_def.h"
#include "FortFunc.h"
#include "LS.h"

using namespace std;
using namespace pm_lib;

int order_of_PM_key;      ///< PMlib用の登録番号カウンタ < PM_NUM_MAX

// #################################################################
REAL_TYPE* allocate_REAL(size_t size, REAL_TYPE init)
{
  REAL_TYPE *var;

  if ( !(var=new REAL_TYPE[size]) ) {
    printf("Failed to allocate memory %d (MB)", (unsigned)size/(1024*1024));
    assert(0);
  }
  for (int i=0; i<size; i++) *(var+i) = init;
  return var;
}

// #################################################################
unsigned* allocate_unsigned(size_t size, unsigned init)
{
  unsigned *var;

  if ( !(var=new unsigned[size]) ) {
    printf("Failed to allocate memory %d (MB)", (unsigned)size/(1024*1024));
    assert(0);
  }
  for (int i=0; i<size; i++) *(var+i) = init;
  return var;
}

// #################################################################
float* allocate_float(size_t size, unsigned init)
{
  float *var;

  if ( !(var=new float[size]) ) {
    printf("Failed to allocate memory %d (MB)", (unsigned)size/(1024*1024));
    assert(0);
  }
  for (int i=0; i<size; i++) *(var+i) = init;
  return var;
}



// #################################################################
int main (int argc, char * const argv[]) {
  if (argc != 6)
    {
      printf("\tUsage : ./cube sz_x, sz_y, sz_z, linear_solver, IterationMax\n");
      printf("\t$ ./cube 128 50 25 pbicgstab 1000\n");
      printf("\t\tlinear_solver = {jacobi | sor | sor2sma | pbicgstab | bicgstab | slor}\n");
      //printf("\t\tpadding_sw = {0-no padding, 1-add padding}\n");
      return 0;
    }

  float mem_size = 0.0;

  int sz_x = atoi(argv[1]);
  int sz_y = atoi(argv[2]);
  int sz_z = atoi(argv[3]);
  printf("\n===============================\n");
  printf("\n\tarray size = (%d %d %d)\n", sz_x, sz_y, sz_z);

  int gc = 2; // guide cell
  //int mode_padding = 0;
  //if ( strcasecmp("np",argv[6]) && strcasecmp("ap",argv[6]) ) exit(0);
  //if ( !strcasecmp("ap",argv[6]) ) mode_padding=1;

	int v_sz[4];
  v_sz[0] = sz_x;
  v_sz[1] = sz_y;
  v_sz[2] = sz_z;
  //int base = v_sz[0] + 2 * gc;
  //int padding = 0;
  //if (mode_padding==1) padding = (base % 2 == 0) ? 1 : 0;

  //v_sz[3] = v_sz[0] + padding;
  //printf("\n\tsz_x=%d gc=%d base=%d padding=%d sz_inner=%d\n",
  //        sz_x, gc, base, padding, v_sz[3]+2*gc);

  size_t m_size = (v_sz[0]+2*gc) * (v_sz[1]+2*gc) * (v_sz[2]+2*gc);
  //size_t m16_size= (v_sz[0]+2*gc) * (v_sz[1]+2*gc) * (v_sz[2]+2*gc) * 16;

  int KindOfPrecondition = 0; // 0-?, 1-??

  REAL_TYPE* p   = NULL;    // pressure
  REAL_TYPE* src0= NULL;    // source term 0
  REAL_TYPE* wrk = NULL;    // work array
  REAL_TYPE* exs = NULL;    // exact solution
  REAL_TYPE* gosa = NULL;
  REAL_TYPE* mask = NULL;

  REAL_TYPE* pcg_q  = NULL; // work for bicgstab
  REAL_TYPE* pcg_r  = NULL;
  REAL_TYPE* pcg_r0 = NULL;
  REAL_TYPE* pcg_p  = NULL;
  REAL_TYPE* pcg_p_ = NULL;
  REAL_TYPE* pcg_s  = NULL;
  REAL_TYPE* pcg_s_ = NULL;
  REAL_TYPE* pcg_t_ = NULL;

  // allocate array
  p   = allocate_REAL(m_size, 0.0);
  src0= allocate_REAL(m_size, 0.0);
  mask= allocate_REAL(m_size, 0.0);
  mem_size += (float)m_size * sizeof(REAL_TYPE) * 3.0;

#ifdef DEBUG
  exs = allocate_REAL(m_size, 0.0);
  gosa= allocate_REAL(m_size, 0.0);
  mem_size += (float)m_size * sizeof(REAL_TYPE) * 2.0;
#endif


  // constant
  double eps = 1.0e-5;    // convergence criteria
  REAL_TYPE ac1 = 1.1;    // acceleration coef. for SOR
  REAL_TYPE ac2 = 0.8;    // acceleration coef. for jacobi relaxation
  REAL_TYPE dh;           // grid width
  REAL_TYPE org[3];       // origin
  REAL_TYPE coef;

  dh  = 1.0/(REAL_TYPE)(sz_x+1);
  for (int i=0; i<3; i++) org[i] = 0.0;

  // type
  int ls_type=0;
  char* q = argv[4];
  printf("\tLinear solver = %s\n", q);

  char fname[20];
  memset(fname, 0, sizeof(char)*20);

  if      ( !strcasecmp(q, "jacobi") ) {
    ls_type = LS_JACOBI;
    strcpy(fname, "jacobi.txt");
    wrk = allocate_REAL(m_size, 0.0);
    coef = ac2;
  }

  else if ( !strcasecmp(q, "psor") ) {
    ls_type = LS_PSOR;
    strcpy(fname, "psor.txt");
    coef = ac1;
  }

  else if ( !strcasecmp(q, "lsor") ) {
    ls_type = LS_SLOR;
    strcpy(fname, "lsor.txt");
    coef = ac1;
  }

  else if ( !strcasecmp(q, "sor2sma") ) {
    ls_type = LS_SOR2SMA;
    strcpy(fname, "sor2sma.txt");
    coef = ac1;
  }

  else if ( !strcasecmp(q, "pbicgstab") ) {
    ls_type = LS_BICGSTAB;
    strcpy(fname, "pbicgstab.txt");
    coef = ac1;

    pcg_q  = allocate_REAL(m_size, 0.0);
    pcg_r  = allocate_REAL(m_size, 0.0);
    pcg_r0 = allocate_REAL(m_size, 0.0);
    pcg_p  = allocate_REAL(m_size, 0.0);
    pcg_p_ = allocate_REAL(m_size, 0.0);
    pcg_s  = allocate_REAL(m_size, 0.0);
    pcg_s_ = allocate_REAL(m_size, 0.0);
    pcg_t_ = allocate_REAL(m_size, 0.0);
  }

  else{
    printf("Invalid solver\n");
    exit(0);
  }

  mem_size /= (1024.0*1024.0*1024.0);
  printf("\nRequired memory size = %12.4f GB\n", mem_size);

  int ItrMax = atoi(argv[5]);

  // Apply BC
  bc_(v_sz, &gc, p, &dh);

#ifdef DEBUG
  // exact solution
  exact_(v_sz, &gc, exs, &dh);

  char fname2[20];
  memset(fname2, 0, sizeof(char)*20);
  strcpy(fname2, "exact.sph");
#endif

  // history title
  FILE* fp;
  if ( !(fp=fopen(fname, "w")) )
  {
    printf("\tSorry, can't open file.\n");
    assert(0);
  }

  //fprintf(fp,"Column_Data_00\n");
  fprintf(fp, "Itration      Residual\n");

  // timing
  int num_thread  = omp_get_max_threads();
  omp_set_num_threads(num_thread);

  PerfMonitor PM;
  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( 0 );
  PM.setParallelMode("OpenMP", num_thread, 1);
  set_timing_label(&PM);

  // mask
  init_mask_(mask, v_sz, &gc);

  // source term
  src_dirichlet_(src0, v_sz, &gc, &dh);

  // LS class
  LinearSolver LS(v_sz,
                  gc,
                  dh,
                  ItrMax,
                  coef,
                  eps,
                  &PM,
                  pcg_p,
                  pcg_p_,
                  pcg_r,
                  pcg_r0,
                  pcg_q,
                  pcg_s,
                  pcg_s_,
                  pcg_t_,
                  wrk,
                  exs,
                  gosa,
                  fp);

  int loop;
  double flop = 0.0;

  // scheme branch
  switch (ls_type) {
    case JACOBI:
      TIMING_start(&PM, "Jacobi");
      loop = LS.Jacobi(p, src0, ItrMax, flop);
      TIMING_stop(&PM, "Jacobi", flop);
      break;

    case SOR:
          TIMING_start(&PM, "PointSOR");
          loop = LS.PointSOR(p, src0, ItrMax, flop, false);
          TIMING_stop(&PM, "PointSOR", flop);
      break;

    case SLOR:
          TIMING_start(&PM, "SLOR");
          loop = LS.lineSOR(p, src0, ItrMax, flop, false);
          TIMING_stop(&PM, "SLOR", flop);
      break;

    case SOR2SMA:
          TIMING_start(&PM, "SOR2_SMA");
          loop = LS.SOR2_SMA(p, src0, ItrMax, flop, false);
          TIMING_stop(&PM, "SOR2_SMA", flop);
      break;

    case PBICGSTAB:
          TIMING_start(&PM, "PBiCGstab");
          loop = LS.PBiCGstab(p, src0, ItrMax, flop, true, KindOfPrecondition);
          TIMING_stop(&PM, "PBiCGstab", flop);
      break;

    case RES_SP:
          for (int i=1; i<=ItrMax; i++) {
            TIMING_start(&PM, "RESsp");
            flop = 0.0;
            loop = LS.residual_sp(p, mask, flop);
            TIMING_stop(&PM, "RESsp", flop);
          }
      break;

    case BICGSTAB:
          TIMING_start(&PM, "BiCGstab");
          loop = LS.PBiCGstab(p, src0, ItrMax, flop, false, KindOfPrecondition);
          TIMING_stop(&PM, "BiCGstab", flop);
      break;

    default:
      assert(0);
      break;
  }


  // close
  if ( !fp ) fclose(fp);

  // file out
#ifdef DEBUG
  strcpy(fname2, "p.sph"); //解ベクトルの出力
  fileout_(v_sz, &gc, p, &dh, org, fname2);

  strcpy(fname2, "e.sph");
  fileout_(v_sz, &gc, gosa, &dh, org, fname2);
#endif

  // profiling
  if ( !(fp=fopen("data-profiling.txt", "w")) )
  {
    printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
    assert(0);
  }

  // 測定結果の集計(gatherメソッドは全ノードで呼ぶこと)
  PM.gather();

  // 結果出力(排他測定のみ)
  printf("\n===============================\n");
  PM.print(stdout, "hoge", "foo");
  PM.print(fp, "hoge", "foo");

  if ( !fp ) fclose(fp);


  return 0;
}
