//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_LS.C
 * @brief  LS Class
 * @author aics
 */


#include "LS.h"


// #################################################################
double LinearSolver::Fdot1(REAL_TYPE* x, double& flop)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy=0.0;

  TIMING_start(PM, "Blas_dot1");
  blas_dot1_(&xy, x, size, &guide, &flop_count);
  TIMING_stop(PM, "Blas_dot1", flop_count);
  flop += flop_count;

  return xy;
}


// #################################################################
double LinearSolver::Fdot2(REAL_TYPE* x, REAL_TYPE* y, double& flop)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy=0.0;

  TIMING_start(PM, "Blas_dot2");
  blas_dot2_(&xy, x, y, size, &guide, &flop_count);
  TIMING_stop(PM, "Blas_dot2", flop_count);
  flop += flop_count;

  return xy;
}

// #################################################################
int LinearSolver::Jacobi(REAL_TYPE* x,
                         REAL_TYPE* b,
                         const int itrMax,
                         double& flop,
                         bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = coef_ac;
  int lc=0;                       /// ループカウント

  // x     圧力 p^{n+1}
  // b     RHS vector

  for (lc=1; lc<=itrMax; lc++)
  {
    double res=0.0;

    // 反復処理
    TIMING_start(PM, "JACOBI_kernel");
    flop_count = 0.0;
    jacobi_(x, size, &guide, cf, &omg, b, &res, wk, &flop_count);
    TIMING_stop(PM, "JACOBI_kernel", flop_count);


    // 境界条件
    if ( converge_check ) {
      TIMING_start(PM, "BoundaryCondition");
      bc_(size, &guide, x, &dh);
      TIMING_stop(PM, "BoundaryCondition");
    }

#ifdef DEBUG
    if ( converge_check ) fprintf(fp, "%8d %12.6e\n", lc, res);
#endif

    flop += flop_count;

    // 収束判定
    if ( converge_check ) {
      if ( res < eps ) break; //epsで比較
    }
  }

  return lc;
}

// #################################################################
int LinearSolver::lineSOR(REAL_TYPE* x,
                           REAL_TYPE* b,
                           const int itrMax,
                           double& flop,
                           bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = coef_ac;
  int lc=0;                       /// ループカウント

  // x     圧力 p^{n+1}
  // b     RHS vector

  for (lc=1; lc<=itrMax; lc++) {
    double res=0.0;

    // 反復処理
    TIMING_start(PM, "SLOR_kernel");
    flop_count = 0.0;
    slor_(x, size, &guide, cf, &omg, b, &res, &flop_count);
    TIMING_stop(PM, "SLOR_kernel", flop_count);


    // 境界条件
    if ( converge_check ) {
      TIMING_start(PM, "BoundaryCondition");
      bc_(size, &guide, x, &dh);
      TIMING_stop(PM, "BoundaryCondition");
    }

#ifdef DEBUG
    if ( converge_check ) fprintf(fp, "%8d %12.6e\n", lc, res);
#endif

    flop += flop_count;

    // 収束判定
    if ( converge_check ) {
      if ( res < eps ) break; //epsで比較
    }
  }

  return lc;
}

// #################################################################
int LinearSolver::PointSOR(REAL_TYPE* x,
                           REAL_TYPE* b,
                           const int itrMax,
                           double& flop,
                           bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = coef_ac;
  int lc=0;                       /// ループカウント

  // x     圧力 p^{n+1}
  // b     RHS vector

  for (lc=1; lc<=itrMax; lc++) {
    double res=0.0;

    // 反復処理
    TIMING_start(PM, "SOR_kernel");
    flop_count = 0.0;
    psor_(x, size, &guide, cf, &omg, b, &res, &flop_count);
    TIMING_stop(PM, "SOR_kernel", flop_count);


    // 境界条件
    if ( converge_check ) {
      TIMING_start(PM, "BoundaryCondition");
      bc_(size, &guide, x, &dh);
      TIMING_stop(PM, "BoundaryCondition");
    }

#ifdef DEBUG
    if ( converge_check ) fprintf(fp, "%8d %12.6e\n", lc, res);
#endif

    flop += flop_count;

    // 収束判定
    if ( converge_check ) {
      if ( res < eps ) break; //epsで比較
    }
  }

  return lc;
}


// #################################################################
int LinearSolver::SOR2_SMA(REAL_TYPE* x,
                           REAL_TYPE* b,
                           const int itrMax,
                           double& flop,
                           bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  int lc=0;                       /// ループカウント
  REAL_TYPE omg = coef_ac;

  // x     圧力 p^{n+1}
  // b     RHS vector

  for (lc=1; lc<=itrMax; lc++) {
    // 2色のマルチカラー(Red&Black)のセットアップ

    // ip = 0 基点(1,1,1)がRからスタート
    //    = 1 基点(1,1,1)がBからスタート
    int ip = 0;
    double res=0.0;

    // 各カラー毎の間に同期, 残差は色間で積算する
    // R - color=0 / B - color=1
    for (int color=0; color<2; color++) {

      TIMING_start(PM, "SOR2SMA_kernel");
      flop_count = 0.0; // 色間で積算しない
      psor2sma_core_(x, size, &guide, cf, &ip, &color, &omg, b, &res, &flop_count);
      TIMING_stop(PM, "SOR2SMA_kernel", flop_count);

      // 境界条件 >> PBiCGstabで前処理のときには境界条件を呼ばない。呼んではならぬ？ bのBCではないから？
      if ( converge_check ) {
        TIMING_start(PM, "BoundaryCondition");
        bc_(size, &guide, x, &dh);
        TIMING_stop(PM, "BoundaryCondition");
      }

      flop += flop_count;
    }

#ifdef DEBUG
    if ( converge_check ) fprintf(fp, "%8d %12.6e\n", lc, res);
#endif

    if ( converge_check ) {
      if ( res < eps ) break; //epsで比較
    }
  }

  return lc;
}



// #################################################################
void LinearSolver::Preconditioner(REAL_TYPE* x,
                                  REAL_TYPE* b,
                                  double& flop,
                                  const bool isPrecond,
                                  int KindOfPrecondition)
{
  // 前処理なし(コピー)
  if ( !isPrecond )
    {
      TIMING_start(PM, "Blas_copy");
      blas_copy_(x, b, size, &guide);
      TIMING_stop(PM, "Blas_copy");
      return;
    }

  int lc_max = 4;

  // 前処理
  if(KindOfPrecondition==0)
  {
    SOR2_SMA(x, b, lc_max, flop, false);
  }
  else if(KindOfPrecondition==1)
  {
    PointSOR(x, b, lc_max, flop, false);
  }
}

// #################################################################
// PBiCBSTAB 収束判定は残差
// @note 反復回数がマシンによって異なる現象がある．
// Xeon E5では同じ反復回数になるのに対して，Core i7では反復回数が試行毎に異なる．
// 内積のOpenMP並列のためか？
int LinearSolver::PBiCGstab(REAL_TYPE* x,
                            REAL_TYPE* b,
                            const int ItrMax,
                            double& flop,
                            const bool isPrecond,
                            int KindOfPrecondition)
{
  double flop_count = 0.0;

  TIMING_start(PM, "Blas_clear");
  blas_clear_(pcg_q, size, &guide);
  TIMING_stop(PM, "Blas_clear");

  TIMING_start(PM, "Blas_residual");
  flop_count = 0.0;
  blas_calc_rk_(pcg_r, x, b, size, &guide, cf, &flop_count);
  TIMING_stop(PM, "Blas_residual", flop_count);
  flop += flop_count;

  TIMING_start(PM, "Blas_copy");
  blas_copy_(pcg_r0, pcg_r, size, &guide);
  TIMING_stop(PM, "Blas_copy");

  double rho_old = 1.0;
  double alpha   = 0.0;
  double omega   = 1.0;
  double r_omega = -omega;
  int lc=0;                      /// ループカウント


  for (lc=1; lc<=ItrMax; lc++)
    {
      flop_count = 0.0;
      double rho = Fdot2(pcg_r, pcg_r0, flop_count);
      flop += flop_count;

      if( fabs(rho) < FLT_MIN )
      {
        lc = 0;
        break;
      }

      if( lc == 1 )
      {
        TIMING_start(PM, "Blas_copy");
        blas_copy_(pcg_p, pcg_r, size, &guide);
        TIMING_stop(PM, "Blas_copy");
      }
      else
      {
        double beta = rho / rho_old * alpha / omega;

        TIMING_start(PM, "Blas_bicg1");
        flop_count = 0.0;
        blas_bicg_1_(pcg_p, pcg_r, pcg_q, &beta, &omega, size, &guide, &flop_count);
        TIMING_stop(PM, "Blas_bicg1", flop_count);
        flop += flop_count;
      }

      TIMING_start(PM, "Blas_clear");
      blas_clear_(pcg_p_, size, &guide);
      TIMING_stop(PM, "Blas_clear");

      flop_count = 0.0;
      Preconditioner(pcg_p_, pcg_p, flop_count, isPrecond, KindOfPrecondition);
      flop += flop_count;

      TIMING_start(PM, "Blas_ax");
      flop_count = 0.0;
      blas_calc_ax_(pcg_q, pcg_p_, size, &guide, cf, &flop_count);
      TIMING_stop(PM, "Blas_ax", flop_count);
      flop += flop_count;

      flop_count = 0.0;
      alpha = rho / Fdot2(pcg_q, pcg_r0, flop_count);
      flop += flop_count;
      double r_alpha = -alpha;

      TIMING_start(PM, "Blas_triad");
      flop_count = 0.0;
      blas_triad_(pcg_s, pcg_q, pcg_r, &r_alpha, size, &guide, &flop_count);
      TIMING_stop(PM, "Blas_triad", flop_count);
      flop += flop_count;

      TIMING_start(PM, "Blas_clear");
      blas_clear_(pcg_s_, size, &guide);
      TIMING_stop(PM, "Blas_clear");

      flop_count = 0.0;
      Preconditioner(pcg_s_, pcg_s, flop_count, isPrecond,KindOfPrecondition);
      flop += flop_count;

      TIMING_start(PM, "Blas_ax");
      flop_count = 0.0;
      blas_calc_ax_(pcg_t_, pcg_s_, size, &guide, cf, &flop_count);
      TIMING_stop(PM, "Blas_ax", flop_count);
      flop += flop_count;

      flop_count = 0.0;
      omega = Fdot2(pcg_t_, pcg_s, flop_count) / Fdot1(pcg_t_, flop_count);
      r_omega = -omega;
      flop += flop_count;

      TIMING_start(PM, "Blas_bicg2");
      flop_count = 0.0;
      blas_bicg_2_(x, pcg_p_, pcg_s_, &alpha , &omega, size, &guide, &flop_count);
      TIMING_stop(PM, "Blas_bicg2", flop_count);
      flop += flop_count;

      TIMING_start(PM, "Blas_triad");
      flop_count = 0.0;
      blas_triad_(pcg_r, pcg_t_, pcg_s, &r_omega, size, &guide, &flop_count);
      TIMING_stop(PM, "Blas_triad", flop_count);
      flop += flop_count;

      flop_count = 0.0;
      double res = Fdot1(pcg_r, flop_count);
      flop += flop_count;

#ifdef DEBUG
      fprintf(fp, "%8d %12.6e\n", lc, res);
#endif

      //if ( res < eps ) break; //epsで比較

      rho_old = rho;
    }

  // 境界条件
  TIMING_start(PM, "BoundaryCondition");
  bc_(size, &guide, x, &dh);
  TIMING_stop(PM, "BoundaryCondition");

  return lc;
}

// #################################################################
// 挙動のおかしなresidualのみ
int LinearSolver::residual_sp(REAL_TYPE* x,
                              REAL_TYPE* b,
                              double& flop)
{
  double flop_count = 0.0;

  TIMING_start(PM, "Blas_residual");
  flop_count = 0.0;
  //blas_calc_rk_(pcg_r, x, b, size, &guide, cf, &flop_count);
  //blas_calc_rk_wob_(pcg_r, x, size, &guide, cf, &flop_count);
  blas_calc_rk_1d_(pcg_r, x, b, size, &guide, cf, &flop_count);
  TIMING_stop(PM, "Blas_residual", flop_count);
  flop += flop_count;

}
