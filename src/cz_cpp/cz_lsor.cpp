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

/////////////////////////////////////////////////////////
void CZ::tdma_pre(REAL_TYPE* a,
                  REAL_TYPE* c,
                  REAL_TYPE* e,
                  REAL_TYPE* w,
                  double& flop)
{
  __assume_aligned(a,  ALIGN_SIZE);
  __assume_aligned(c,  ALIGN_SIZE);
  __assume_aligned(e,  ALIGN_SIZE);
  __assume_aligned(w,  ALIGN_SIZE);

  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  int nx = ked - kst + 1;

  flop += 14.0*(double)(nx-1);
  REAL_TYPE f;

  w[0] = c[0];

  // Forward:Peel
  for (int k=1; k<nx; k++)
  {
    f = 1.0 / (1.0 - a[k] * w[k-1]);
    w[k] = f * c[k];
    e[k] = f;
  }

}


/////////////////////////////////////////////////////////
void CZ::lsor_j(REAL_TYPE* d,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* e,
                REAL_TYPE* rhs,
                REAL_TYPE* msk,
                REAL_TYPE* d2,
                double &res,
                double &flop)
{
  int gc = GUIDE;
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  double flop_count = 0.0;
  res = 0.0;

  TIMING_start("LSOR_TDMA_Ex");
  #pragma omp parallel for reduction(+:res) schedule(dynamic,1) reduction(+:flop_count)
  for (int j=jst; j<=jed; j++) {

    lsor_lu_rhs_j_(d2, size, innerFidx, &gc, &j, x, rhs, &flop_count);

    lsor_lu_fwd_(d, size, innerFidx, &gc, &j, x, a, e, d2, msk, rhs, &flop_count);

    lsor_tdma_b_(d, size, innerFidx, &gc, &j, w, &flop_count);

    lsor_relax_(d, size, innerFidx, &gc, &j, x, msk, &ac1, &res, &flop_count);
  }
  TIMING_stop("LSOR_TDMA_Ex", flop_count);
}


/////////////////////////////////////////////////////////
void CZ::lsor_j4(REAL_TYPE* d,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* e,
                REAL_TYPE* rhs,
                REAL_TYPE* msk,
                REAL_TYPE* d2,
                double &res,
                double &flop)
{
  int gc = GUIDE;
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  double flop_count = 0.0;
  res = 0.0;

  TIMING_start("LSOR_TDMA_Ex");
  #pragma omp parallel for reduction(+:res) schedule(dynamic,1) reduction(+:flop_count)
  for (int j=jst; j<=jed; j++) {

    lsor_lu_rhs_j_(d2, size, innerFidx, &gc, &j, x, rhs, &flop_count);

    for (int l=0; l<2; l++) {
      lsor_tdma_fwd_(d, size, innerFidx, &gc, &j, d2, x, a, e, msk, rhs, &flop_count);
      lsor_tdma_b4_(d, size, innerFidx, &gc, &j, w, &flop_count);
      lsor_relax4_(d, size, innerFidx, &gc, &j, x, msk, &ac1, &res, &flop_count);
    }
  }
  TIMING_stop("LSOR_TDMA_Ex", flop_count);
}

/////////////////////////////////////////////////////////
void CZ::lsor_k(REAL_TYPE* d,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* e,
                REAL_TYPE* rhs,
                REAL_TYPE* msk,
                REAL_TYPE* d2,
                double &res,
                double &flop)
{
  constexpr int ITR_INNER = 2;

  int gc = GUIDE;
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  double residual[ITR_INNER];
  memset(residual, 0, sizeof(double)*ITR_INNER);

  double flop_count = 0.0;

  TIMING_start("LSOR_TDMA_Ex");
  #pragma omp parallel for reduction(+:res) schedule(dynamic,1) reduction(+:flop_count)
  for (int j=jst; j<=jed; j++) {
    lsor_lu_rhs_j_(d2, size, innerFidx, &gc, &j, x, rhs, &flop_count);

    for (int l=0; l<ITR_INNER; l++) {
      lsor_inner_(d, size, innerFidx, &gc, &j, d2, x, a, e, w, msk, rhs,
                     &ac1, &residual[l], &flop_count);
    }
  }
  TIMING_stop("LSOR_TDMA_Ex", flop_count);
/*
  for (int m=0; m<ITR_INNER; m++) {
    printf("%e ", residual[m]);
  }
  printf("\n");
*/
  res = residual[0];
}


/////////////////////////////////////////////////////////
void CZ::lsor_k2(REAL_TYPE* d,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* e,
                REAL_TYPE* rhs,
                REAL_TYPE* msk,
                REAL_TYPE* d2,
                double &res,
                double &flop)
{
  constexpr int ITR_INNER = 2;
  int ItrInner = ITR_INNER;
  int gc = GUIDE;
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  double resD[ITR_INNER];
  memset(resD, 0, sizeof(double)*ITR_INNER);

  double flop_count = 0.0;

  TIMING_start("LSOR_TDMA_Ex");
  //lsor_inner_b_(d, size, innerFidx, &gc, d2, x, a, e, w, msk, rhs,
  //                &ac1, &ItrInner, resD, &flop_count);
  //lsor_inner_rb_(d, size, innerFidx, &gc, d2, x, a, e, w, msk, rhs,
  //                &ac1, &ItrInner, resD, &flop_count);
  //lsor_inner_cb_(d, size, innerFidx, &gc, d2, x, a, e, w, msk, rhs,
  //                &ac1, &ItrInner, resD, &flop_count);
  lsor_inner_d_(d, size, innerFidx, &gc, d2, x, a, e, w, msk, rhs,
                  &ac1, &ItrInner, resD, &flop_count);
  TIMING_stop("LSOR_TDMA_Ex", flop_count);

  res = resD[0];
}

/////////////////////////////////////////////////////////
double CZ::lsor_k3(REAL_TYPE* d,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* e,
                REAL_TYPE* rhs,
                REAL_TYPE* msk,
                REAL_TYPE* d2,
                double* resD,
                double* rd,
                double &flop)
{
  int gc = GUIDE;
  int cl_sz = CL_SZ;
  double flop_count = 0.0;
  int ItrInner = 2;
  memset(resD, 0, sizeof(double)*CL_SZ*numThreads);
  memset(rd, 0, sizeof(double)*CL_SZ);

  TIMING_start("LSOR_TDMA_Ex");
  lsor_inner_c_(d, size, innerFidx, &gc, d2, x, a, e, w, msk, rhs,
                  &ac1, &cl_sz, &numThreads, &ItrInner, resD, &flop_count);

  for (int i=0; i<ItrInner; i++) {
    //printf("Itr=%d : ", i+1);
    for (int j=0; j<numThreads; j++) {
      int m = j*CL_SZ + i;
      rd[i] += resD[m];
      //printf("%e ", resD[m]);
    }
    //printf(": %e\n", rd[i]);
  }
  /*
  for (int i=0; i<ItrInner; i++) {
    rd[i] += resD[i];
    printf("%e ", rd[i]);
  }
  printf("\n");
  */
  TIMING_stop("LSOR_TDMA_Ex", flop_count);

  return rd[0];
}
