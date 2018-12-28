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

/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note simd intrinsic
 */
void CZ::tdma_pre(REAL_TYPE* a,
                  REAL_TYPE* c,
                  REAL_TYPE* e,
                  REAL_TYPE* w,
                  double& flop)
{
  __assume_aligned(a,  ALIGN);
  __assume_aligned(c,  ALIGN);
  __assume_aligned(e,  ALIGN);
  __assume_aligned(w,  ALIGN);

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


/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 */
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

  #pragma omp parallel for reduction(+:res) schedule(dynamic,1) private(flop_count)
  for (int j=jst; j<=jed; j++) {

    TIMING_start("LSOR_RHS_J");
    flop_count = 0.0;
    lsor_lu_rhs_j_(d2, size, innerFidx, &gc, &j, x, rhs, &flop_count);
    TIMING_stop("LSOR_RHS_J", flop_count);

    for (int k=kst; k<=ked; k++) {

      TIMING_start("LSOR_RHS_K");
      flop_count = 0.0;
      lsor_lu_rhs_k_(d, size, innerFidx, &gc, &j, &k, x, d2, msk, &flop_count);
      TIMING_stop("LSOR_RHS_K", flop_count);

      TIMING_start("LSOR_TDMA_BC");
      flop_count = 0.0;
      if (k == kst)
        lsor_lu_bc_kst_(d, size, innerFidx, &gc, &j, rhs, msk, &flop_count);
      if (k == ked)
        lsor_lu_bc_ked_(d, size, innerFidx, &gc, &j, rhs, msk, &flop_count);
      TIMING_stop("LSOR_TDMA_BC", flop_count);

      if (k>2)
      {
        TIMING_start("LSOR_TDMA_FWD");
        flop_count = 0.0;
        lsor_tdma_f_(d, size, innerFidx, &gc, &j, &k, a, e, &flop_count);
        TIMING_stop("LSOR_TDMA_FWD", flop_count);
      }
    }

    TIMING_start("LSOR_TDMA_BWD");
    flop_count = 0.0;
    lsor_tdma_b_(d, size, innerFidx, &gc, &j, w, &flop_count);
    TIMING_stop("LSOR_TDMA_BWD", flop_count);

    TIMING_start("LSOR_Relax_Ex");
    flop_count = 0.0;
    lsor_relax_(d, size, innerFidx, &gc, &j, x, msk, &ac1, &res, &flop_count);
    TIMING_stop("LSOR_Relax_Ex", flop_count);
  }
}



/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 * @note
 */
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

  #pragma omp parallel for reduction(+:res) schedule(dynamic,1) private(flop_count)
  for (int j=jst; j<=jed; j++) {

    TIMING_start("LSOR_RHS_J");
    flop_count = 0.0;
    lsor_lu_rhs_j4_(d2, size, innerFidx, &gc, &j, x, rhs, &flop_count);
    TIMING_stop("LSOR_RHS_J", flop_count);

    for (int l=0; l<2; l++) {

      TIMING_start("LSOR_TDMA_FWD");
      flop_count = 0.0;
      lsor_tdma_fwd_(d, size, innerFidx, &gc, &j, d2, x, a, e, msk, rhs, &flop_count);
      TIMING_stop("LSOR_TDMA_FWD", flop_count);

      TIMING_start("LSOR_TDMA_BWD");
      flop_count = 0.0;
      lsor_tdma_b4_(d, size, innerFidx, &gc, &j, w, &flop_count);
      TIMING_stop("LSOR_TDMA_BWD", flop_count);

      TIMING_start("LSOR_Relax_Ex");
      flop_count = 0.0;
      lsor_relax4_(d, size, innerFidx, &gc, &j, x, msk, &ac1, &res, &flop_count);
      TIMING_stop("LSOR_Relax_Ex", flop_count);
    }
  }
}


/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 * @note
 */
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

  #pragma omp parallel for reduction(+:res) schedule(dynamic,1) private(flop_count)
  for (int j=jst; j<=jed; j++) {

    TIMING_start("LSOR_RHS_J");
    flop_count = 0.0;
    lsor_lu_rhs_j4_(d2, size, innerFidx, &gc, &j, x, rhs, &flop_count);
    TIMING_stop("LSOR_RHS_J", flop_count);

    for (int l=0; l<ITR_INNER; l++) {
      TIMING_start("LSOR_TDMA_Ex");
      flop_count = 0.0;
      lsor_inner_(d, size, innerFidx, &gc, &j, d2, x, a, e, w, msk, rhs,
                     &ac1, &residual[l], &flop_count);
      TIMING_stop("LSOR_TDMA_Ex", flop_count);
    }
  }
/*
  for (int m=0; m<ITR_INNER; m++) {
    printf("%e ", residual[m]);
  }
  printf("\n");
*/
  res = residual[0];
}
