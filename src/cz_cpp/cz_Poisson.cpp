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

/**
 * @file   cz_Poisson.cpp
 * @brief  RIAM class
 * @author RIIT
 */

 /* #################################################################
 * @brief JACOBI反復
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 */
 int CZ::JACOBI(double& res, REAL_TYPE* X, REAL_TYPE* B)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;

    for (itr=1; itr<=ItrMax; itr++)
    {
      res = 0.0;

      TIMING_start("JACOBI_kernel");
      flop_count = 0.0;
      jacobi_(X, size, innerFidx, &gc, cf, &ac2, B, &res, WRK, &flop_count);
      TIMING_stop("JACOBI_kernel", flop_count);


      if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

      if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;


      TIMING_start("BoundaryCondition");
      bc_(size, &gc, X, pitch, origin, nID);
      TIMING_stop("BoundaryCondition");

      res *= res_normal;
      res = sqrt(res);
      Hostonly_ fprintf(fph, "%6d  %13.6e\n", itr, res);

      if ( res < eps ) break;
    }

    return itr;
  }


 /* #################################################################
 * @brief SOR反復
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 */
 int CZ::PSOR(double& res, REAL_TYPE* X, REAL_TYPE* B)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;

    for (itr=1; itr<=ItrMax; itr++)
    {
      res = 0.0;

      TIMING_start("SOR_kernel");
      flop_count = 0.0;
      psor_(X, size, innerFidx, &gc, cf, &ac1, B, &res, &flop_count);
      TIMING_stop("SOR_kernel", flop_count);


      if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;
      if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;


      TIMING_start("BoundaryCondition");
      bc_(size, &gc, X, pitch, origin, nID);
      TIMING_stop("BoundaryCondition");

      res *= res_normal;
      res = sqrt(res);
      fprintf(fph, "%6d  %13.6e\n", itr, res);

      if ( res < eps ) break;

    } // Iteration

    return itr;
  }


/* #################################################################
* @brief SOR2SMA反復
* @param [in,out] X      解ベクトル
* @param [in]     B      RHSベクトル
*/
int CZ::RBSOR(double& res, REAL_TYPE* X, REAL_TYPE* B)
 {
   int itr;
   double flop_count = 0.0;
   int gc = GUIDE;

   for (itr=1; itr<=ItrMax; itr++)
   {
     res = 0.0;

     // 2色のマルチカラー(Red&Black)のセットアップ
     // ip = 0 基点(1,1,1)が Rからスタート
     //    = 1 基点(1,1,1)が Bからスタート
     int ip=0;
     if ( numProc > 1 )
     {
       ip = (head[0]+head[1]+head[2]+1) % 2;
     }
     else
     {
       ip = 0;
     }

     // 各カラー毎の間に同期, 残差は色間で積算する
     // R - color=0 / B - color=1
     TIMING_start("SOR2SMA_kernel");
     flop_count = 0.0;
     for (int color=0; color<2; color++)
     {
       // res_p >> 反復残差の二乗和
       psor2sma_core_(X, size, innerFidx, &gc, cf, &ip, &color, &ac1, B, &res, &flop_count);
     }
     TIMING_stop("SOR2SMA_kernel", flop_count);


     if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;
     if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;


     TIMING_start("BoundaryCondition");
     bc_(size, &gc, X, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     res *= res_normal;
     res = sqrt(res);
     fprintf(fph, "%6d  %13.6e\n", itr, res);

     if ( res < eps ) break;

   } // Iteration

   return itr;
 }

/*
 // #################################################################
 double CZ::Fdot1(REAL_TYPE* x)
 {
   double flop_count=0.0;          /// 浮動小数点演算数
   double xy = 0.0;

   TIMING_start("Dot1");
   blas_dot1_(&xy, x, size, &flop_count);
   TIMING_stop("Dot1", flop_count);

   if ( !Comm_SUM_1(&xy, "A_R_Dot") ) Exit(0);

   return xy;
 }

 // #################################################################
 double RIAMC::Fdot2(REAL_TYPE* x, REAL_TYPE* y)
 {
   double flop_count=0.0;          /// 浮動小数点演算数
   double xy = 0.0;

   TIMING_start("Dot2");
   blas_dot2_(&xy, x, y, size, &flop_count);
   TIMING_stop("Dot2", flop_count);

   if ( !Comm_SUM_1(&xy, "A_R_Dot") ) Exit(0);

   return xy;
 }

 // #################################################################
 void RIAMC::Preconditioner(REAL_TYPE* X, REAL_TYPE* B)
 {
   double res_p = 0.0; // res_p >> 反復残差の二乗和
   double rms_p = 0.0; // rms_p >> 圧力修正量の二乗和

   // 前処理なし(コピー)
   if ( !ItrCtl.precondition )
   {
     TIMING_start("Blas_Copy");
     blas_copy_(X, B, size);
     TIMING_stop("Blas_Copy");
     return;
   }

   int l = SOR(rms_p, res_p, ItrCtl.InnerItr, X, B);
 }


// #################################################################
// @brief PBiCGSTAB反復
// @param [in,out] res    残差
// @param [in,out] rms    修正量のRMS値
// @param [in]     ItrMax 最大反復数
// @param [in,out] X      解ベクトル
// @param [in]     B      RHSベクトル
 int RIAMC::PBiCGstab(double& res,
                      double& rms,
                      const int ItrMax,
                      REAL_TYPE* X,
                      REAL_TYPE* B)
 {
   REAL_TYPE Omega = ItrCtl.omg;
   REAL_TYPE dh = (REAL_TYPE)pitch[0];
   double flop = 0.0;
   res = 0.0;
   rms = 0.0;

   TIMING_start("Blas_Clear");
   blas_clear_(pcg_q , size);
   TIMING_stop("Blas_Clear");

   TIMING_start("Blas_Residual");
   flop = 0.0;
   calc_rk_maf_(size, innerFidx, pcg_r, B, X, Z, &dh, &flop);
   TIMING_stop("Blas_Residual", flop);

   if ( !Comm_S(pcg_r, 1, "Comm_Poisson") ) return 0;

   TIMING_start("Blas_Copy");
   blas_copy_(pcg_r0, pcg_r, size);
   TIMING_stop("Blas_Copy");

   double rho_old = 1.0;
   double alpha = 0.0;
   double omega  = 1.0;
   double r_omega = -omega;
   int lc=0;                      /// ループカウント

   for (lc=1; lc<ItrMax; lc++)
   {
     double rho = Fdot2(pcg_r, pcg_r0);

     if( fabs(rho) < FLT_MIN )
     {
       lc = 0;
       break;
     }

     if( lc == 1 )
     {
       TIMING_start("Blas_Copy");
       blas_copy_(pcg_p, pcg_r, size);
       TIMING_stop("Blas_Copy");
     }
     else
     {
       double beta = rho / rho_old * alpha / omega;

       TIMING_start("Blas_BiCG_1");
       flop = 0.0;
       blas_bicg_1_(pcg_p, pcg_r, pcg_q, &beta, &omega, size, &flop);
       TIMING_stop("Blas_BiCG_1", flop);
     }

     if ( !Comm_S(pcg_p, 1, "Comm_Poisson") ) return 0;

     TIMING_start("Blas_Clear");
     blas_clear_(pcg_p_ , size);
     TIMING_stop("Blas_Clear");

     Preconditioner(pcg_p_, pcg_p);

     TIMING_start("Blas_AX");
     flop = 0.0;
     calc_ax_maf_(size, innerFidx, pcg_q, pcg_p_, Z, &dh, &flop);
     TIMING_stop("Blas_AX", flop);

     alpha = rho / Fdot2(pcg_q, pcg_r0);

     double r_alpha = -alpha;
     TIMING_start("Blas_TRIAD");
     flop = 0.0;
     blas_triad_(pcg_s, pcg_q, pcg_r, &r_alpha, size, &flop);
     TIMING_stop("Blas_TRIAD", flop);

     if ( !Comm_S(pcg_s, 1, "Comm_Res_Poisson") ) return 0;

     TIMING_start("Blas_Clear");
     blas_clear_(pcg_s_ , size);
     TIMING_stop("Blas_Clear");

     Preconditioner(pcg_s_, pcg_s);

     TIMING_start("Blas_AX");
     flop = 0.0;
     calc_ax_maf_(size, innerFidx, pcg_t_, pcg_s_, Z, &dh, &flop);
     TIMING_stop("Blas_AX", flop);

     omega = Fdot2(pcg_t_, pcg_s) / Fdot1(pcg_t_);
     r_omega = -omega;

     TIMING_start("Blas_BiCG_2");
     flop = 0.0;
     blas_bicg_2_(P, pcg_p_, pcg_s_, &alpha , &omega, size, &rms, &flop);
     TIMING_stop("Blas_BiCG_2", flop);

     TIMING_start("Blas_TRIAD");
     flop = 0.0;
     blas_triad_  (pcg_r, pcg_t_, pcg_s, &r_omega, size, &flop);
     TIMING_stop("Blas_TRIAD", flop);

     res = Fdot1(pcg_r);


     if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;
     if ( !Comm_SUM_2(&rms, &res, "Comm_Res_Poisson") ) return 0;


     // RMS \sqrt{ \frac{1}{N} \sum_{res^2} }
     rms *= rms_normal;
     res *= rms_normal;

     rms = sqrt(rms);
     res = sqrt(res);

     // 初期残差
     if (lc==1) res_r0 = res;


     // convergence check
     if ( ItrCtl.ResNorm == nrm_r_rms )
     {
       if (rms <= ItrCtl.eps) break;
     }
     else if ( ItrCtl.ResNorm == nrm_r_res )
     {
       if (res <= ItrCtl.eps) break;
     }
     else if ( ItrCtl.ResNorm == nrm_r_r0 )
     {
       res = (res_r0>ItrCtl.eps) ? res/res_r0 : res/ItrCtl.eps;
       if (res <= ItrCtl.eps) break;
     }

     rho_old = rho;
   } // lc

   //TIMING_start("Poisson_BC");
   //bc_new_prs_(size, innerFidx, &RE, x, Z, &dh, V, nID, &flop_count);
   //TIMING_stop("Poisson_BC");

   if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

   return lc;
 }
*/
