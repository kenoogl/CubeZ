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
 * @brief  CZ class
 * @author RIIT
 */

 /* #################################################################
 * @brief JACOBI反復
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 * @param [in]     s_type ソルバーの指定
 * @param [in]     converge_check 0のとき、収束判定しない
 */
 int CZ::JACOBI(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop,
                int s_type,
                bool converge_check)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;

    for (itr=1; itr<=itr_max; itr++)
    {
      res = 0.0;

      if (s_type==LS_JACOBI_MAF)
      {
        TIMING_start("JACOBI_MAF_kernel");
        flop_count = 0.0;
        jacobi_maf_(X, size, innerFidx, &gc, xc, yc, zc, &ac1, B, &res, WRK, &flop_count);
        TIMING_stop("JACOBI_MAF_kernel", flop_count);
      }
      else
      {
        TIMING_start("JACOBI_kernel");
        flop_count = 0.0;
        jacobi_(X, size, innerFidx, &gc, cf, &ac1, B, &res, WRK, &flop_count);
        TIMING_stop("JACOBI_kernel", flop_count);
      }
      flop += flop_count;

      if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;


      if ( converge_check ) {
        if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

        res *= res_normal;
        res = sqrt(res);
        Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

        TIMING_start("BoundaryCondition");
        bc_k_(size, &gc, X, pitch, origin, nID);
        TIMING_stop("BoundaryCondition");

        if ( res < eps ) break;
      }
    }

    return itr;
  }


 /* #################################################################
 * @brief SOR反復
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 * @param [in]     s_type ソルバーの指定
 * @param [in]     converge_check 0のとき、収束判定しない
 */
 int CZ::PSOR(double& res, REAL_TYPE* X, REAL_TYPE* B,
              const int itr_max, double& flop,
              int s_type,
              bool converge_check)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;

    for (itr=1; itr<=itr_max; itr++)
    {
      res = 0.0;

      if (s_type==LS_PSOR_MAF)
      {
        TIMING_start("SOR_MAF_kernel");
        flop_count = 0.0;
        psor_maf_(X, size, innerFidx, &gc, xc, yc, zc, &ac1, B, &res, &flop_count);
        TIMING_stop("SOR_MAF_kernel", flop_count);
      }
      else
      {
        TIMING_start("SOR_kernel");
        flop_count = 0.0;
        psor_(X, size, innerFidx, &gc, cf, &ac1, B, &res, &flop_count);
        TIMING_stop("SOR_kernel", flop_count);
      }
      flop += flop_count;

      if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;



      if ( converge_check ) {
        if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

        res *= res_normal;
        res = sqrt(res);

        Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

        TIMING_start("BoundaryCondition");
        bc_k_(size, &gc, X, pitch, origin, nID);
        TIMING_stop("BoundaryCondition");

        if ( res < eps ) break;
      }

    } // Iteration

    return itr;
  }


/* #################################################################
* @brief SOR2SMA反復
* @param [in,out] res    残差
* @param [in,out] X      解ベクトル
* @param [in]     B      RHSベクトル
* @param [in]     itr_max 最大反復数
* @param [in]     flop   浮動小数点演算数
* @param [in]     s_type ソルバーの指定
* @param [in]     converge_check 0のとき、収束判定しない
*/
int CZ::RBSOR(double& res, REAL_TYPE* X, REAL_TYPE* B,
              const int itr_max, double& flop,
              int s_type,
              bool converge_check)
 {
   int itr;
   double flop_count = 0.0;
   int gc = GUIDE;

   for (itr=1; itr<=itr_max; itr++)
   {
     flop_count = 0.0;
     res = 0.0;

     // 2色のマルチカラー(Red&Black)のセットアップ
     // ip = 0 基点(2,2,2)が Rからスタート
     //    = 1 基点(2,2,2)が Bからスタート
     int ip=0;
     if ( numProc > 1 )
     {
       ip = (head[0] + head[1] + head[2]+1) % 2;
     }
     else
     {
       ip = 0;
     }

     // 各カラー毎の間に同期, 残差は色間で積算する
     // R - color=0 / B - color=1
     if (s_type==LS_SOR2SMA_MAF)
     {
       TIMING_start("SOR2SMA_MAF_kernel");
       flop_count = 0.0;
       for (int color=0; color<2; color++)
       {
         // res_p >> 反復残差の二乗和
         psor2sma_core_maf_(X, size, innerFidx, &gc, xc, yc, zc, &ip, &color, &ac1, B, &res, &flop_count);
       }
       TIMING_stop("SOR2SMA_MAF_kernel", flop_count);
     }
     else
     {
       TIMING_start("SOR2SMA_kernel");
       flop_count = 0.0;
       for (int color=0; color<2; color++)
       {
         // res_p >> 反復残差の二乗和
         psor2sma_core_(X, size, innerFidx, &gc, cf, &ip, &color, &ac1, B, &res, &flop_count);
       }
       TIMING_stop("SOR2SMA_kernel", flop_count);
     }
     flop += flop_count;
     

     if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;


     if ( converge_check ) {
       if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

       res *= res_normal;
       res = sqrt(res);
       Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

       TIMING_start("BoundaryCondition");
       bc_k_(size, &gc, X, pitch, origin, nID);
       TIMING_stop("BoundaryCondition");

       if ( res < eps ) break;
     }

   } // Iteration

   return itr;
 }


 // #################################################################
 REAL_TYPE CZ::Fdot1(REAL_TYPE* x, double& flop)
 {
   double flop_count=0.0;          /// 浮動小数点演算数
   REAL_TYPE xy = 0.0;
   int gc = GUIDE;

   TIMING_start("Dot1");
   blas_dot1_(&xy, x, size, innerFidx, &gc, &flop_count);
   TIMING_stop("Dot1", flop_count);
   flop += flop_count;

   if ( !Comm_SUM_1(&xy, "A_R_Dot") ) Exit(0);

   return xy;
 }

 // #################################################################
 REAL_TYPE CZ::Fdot2(REAL_TYPE* x, REAL_TYPE* y, double& flop)
 {
   double flop_count=0.0;          /// 浮動小数点演算数
   REAL_TYPE xy = 0.0;
   int gc = GUIDE;

   TIMING_start("Dot2");
   blas_dot2_(&xy, x, y, size, innerFidx, &gc, &flop_count);
   TIMING_stop("Dot2", flop_count);
   flop += flop_count;

   if ( !Comm_SUM_1(&xy, "A_R_Dot") ) Exit(0);

   return xy;
 }

 // #################################################################
 void CZ::Preconditioner(REAL_TYPE* xx,
                         REAL_TYPE* bb,
                         double& flop,
                         int s_type)
 {
   int gc = GUIDE;
   double res = 0.0;
   int lc_max = 4;

   switch (s_type)
   {
     case LS_JACOBI:
     case LS_JACOBI_MAF:
       JACOBI(res, xx, bb, lc_max, flop, s_type, false);
       break;

     case LS_PSOR:
     case LS_PSOR_MAF:
       PSOR(res, xx, bb, lc_max, flop, s_type, false);
       break;

     case LS_SOR2SMA:
     case LS_SOR2SMA_MAF:
       RBSOR(res, xx, bb, lc_max, flop, s_type, false);
       break;
     
     case LS_PCR:
     case LS_PCR_MAF:
       LSOR_PCR(res, xx, bb, lc_max, flop, s_type, false);
       break;
       
     case LS_PCR_RB:
     case LS_PCR_RB_MAF:
       LSOR_PCR_RB(res, xx, bb, lc_max, flop, s_type, false);
       break;

     default:
       blas_copy_(xx, bb, size, &gc);
   }
 }


// #################################################################
// @brief PBiCGSTAB反復
// @param [in,out] res    残差
// @param [in,out] X      解ベクトル
// @param [in]     B      RHSベクトル
// @param [in]     flop   浮動小数点演算数
// @param [in]     s_type ソルバーの指定
 int CZ::PBiCGSTAB(double& res,
                   REAL_TYPE* X,
                   REAL_TYPE* B,
                   double& flop,
                   int s_type)
 {
   int itr;
   double flop_count = 0.0;
   int gc = GUIDE;
   res = 0.0;

   TIMING_start("Blas_Clear");
   blas_clear_(pcg_q , size, &gc);
   TIMING_stop("Blas_Clear");

   
   TIMING_start("Blas_Residual");
   flop_count = 0.0;
   if (s_type==LS_BICGSTAB_MAF)
   {
     calc_rk_maf_(pcg_r, X, B, size, innerFidx, &gc, xc, yc, zc, pvt, &flop_count);
   }
   else
   {
     blas_calc_rk_(pcg_r, X, B, size, innerFidx, &gc, cf, &flop_count);
   }
   TIMING_stop("Blas_Residual", flop_count);
   flop += flop_count;

   
   if ( !Comm_S(pcg_r, 1, "Comm_Poisson") ) return 0;

   TIMING_start("Blas_Copy");
   blas_copy_(pcg_r0, pcg_r, size, &gc);
   TIMING_stop("Blas_Copy");

   REAL_TYPE rho_old = 1.0;
   REAL_TYPE alpha = 0.0;
   REAL_TYPE omega  = 1.0;
   REAL_TYPE r_omega = -omega;

   for (itr=1; itr<ItrMax; itr++)
   {
     flop_count = 0.0;
     REAL_TYPE rho = Fdot2(pcg_r, pcg_r0, flop_count);
     flop += flop_count;

     if( fabs(rho) < FLT_MIN )
     {
       itr = 0;
       break;
     }

     if( itr == 1 )
     {
       TIMING_start("Blas_Copy");
       blas_copy_(pcg_p, pcg_r, size, &gc);
       TIMING_stop("Blas_Copy");
     }
     else
     {
       REAL_TYPE beta = rho / rho_old * alpha / omega;

       TIMING_start("Blas_BiCG_1");
       flop_count = 0.0;
       blas_bicg_1_(pcg_p, pcg_r, pcg_q, &beta, &omega, size, innerFidx, &gc, &flop_count);
       TIMING_stop("Blas_BiCG_1", flop_count);
       flop += flop_count;
     }

     if ( !Comm_S(pcg_p, 1, "Comm_Poisson") ) return 0;

     TIMING_start("Blas_Clear");
     blas_clear_(pcg_p_ , size, &gc);
     TIMING_stop("Blas_Clear");

     flop_count = 0.0;
     Preconditioner(pcg_p_, pcg_p, flop_count, pc_type);
     flop += flop_count;

     
     TIMING_start("Blas_AX");
     flop_count = 0.0;
     if (s_type==LS_BICGSTAB_MAF)
     {
       calc_ax_maf_(pcg_q, pcg_p_, size, innerFidx, &gc, xc, yc, zc, pvt, &flop_count);
     }
     else
     {
       blas_calc_ax_(pcg_q, pcg_p_, size, innerFidx, &gc, cf, &flop_count);
     }
     TIMING_stop("Blas_AX", flop_count);
     flop += flop_count;

     flop_count = 0.0;
     alpha = rho / Fdot2(pcg_q, pcg_r0, flop_count);
     flop += flop_count;

     
     REAL_TYPE r_alpha = -alpha;
     TIMING_start("Blas_TRIAD");
     flop_count = 0.0;
     blas_triad_(pcg_s, pcg_q, pcg_r, &r_alpha, size, innerFidx, &gc, &flop_count);
     TIMING_stop("Blas_TRIAD", flop_count);
     flop += flop_count;

     if ( !Comm_S(pcg_s, 1, "Comm_Res_Poisson") ) return 0;

     TIMING_start("Blas_Clear");
     blas_clear_(pcg_s_ , size, &gc);
     TIMING_stop("Blas_Clear");

     flop_count = 0.0;
     Preconditioner(pcg_s_, pcg_s, flop_count, pc_type);
     flop += flop_count;

     
     TIMING_start("Blas_AX");
     flop_count = 0.0;
     if (s_type==LS_BICGSTAB_MAF)
     {
       calc_ax_maf_(pcg_t_, pcg_s_, size, innerFidx, &gc, xc, yc, zc, pvt, &flop_count);
     }
     else
     {
       blas_calc_ax_(pcg_t_, pcg_s_, size, innerFidx, &gc, cf, &flop_count);
     }
     TIMING_stop("Blas_AX", flop_count);
     flop += flop_count;

     
     flop_count = 0.0;
     omega = Fdot2(pcg_t_, pcg_s, flop_count) / Fdot1(pcg_t_, flop_count);
     r_omega = -omega;
     flop += flop_count;

     TIMING_start("Blas_BiCG_2");
     flop_count = 0.0;
     blas_bicg_2_(X, pcg_p_, pcg_s_, &alpha , &omega, size, innerFidx, &gc, &flop_count);
     TIMING_stop("Blas_BiCG_2", flop_count);
     flop += flop_count;

     TIMING_start("Blas_TRIAD");
     flop_count = 0.0;
     blas_triad_(pcg_r, pcg_t_, pcg_s, &r_omega, size, innerFidx, &gc, &flop_count);
     TIMING_stop("Blas_TRIAD", flop_count);
     flop += flop_count;

     flop_count = 0.0;
     res = Fdot1(pcg_r, flop_count);
     flop += flop_count;


     if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

     if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

     res *= res_normal;
     res = sqrt(res);
     Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

     if ( res < eps ) break;

     rho_old = rho;
   } // itr

   TIMING_start("BoundaryCondition");
   bc_k_(size, &gc, X, pitch, origin, nID);
   TIMING_stop("BoundaryCondition");

   return itr;
 }



/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 * @param [in]     s_type ソルバーの指定
 * @note LSOR_P6から、マルチカラー化
 */
int CZ::LSOR_PCR_RB(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop,
                int s_type,
                bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int gc = GUIDE;
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int n = ked - kst + 1;
  int pn;
  
  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  
  
  for (itr=1; itr<=itr_max; itr++)
  {
    flop_count = 0.0;
    res = 0.0;
    
    // 2色のマルチカラー(Red&Black)のセットアップ
    int ip=0;
    if ( numProc > 1 )
    {
      ip = (head[0] + head[1] + head[2]+1) % 2;
    }
    else
    {
      ip = 0;
    }
    
    if (s_type==LS_PCR_RB_MAF)
    {
      TIMING_start("PCR_RB_MAF");
      for (int color=0; color<2; color++)
      {
        pcr_rb_maf_(size, innerFidx, &gc, &pn, &ip, &color, X, MSK, B, xc, yc, zc,
                    WA, WC, WD, WAA, WCC, WDD,
                    &ac1, &res, &flop_count);
      }
      TIMING_stop("PCR_RB_MAF", flop_count);
    }
    else
    {
      TIMING_start("PCR_RB");
      for (int color=0; color<2; color++)
      {
        pcr_rb_(size, innerFidx, &gc, &pn, &ip, &color, X, MSK, B, 
                WA, WC, WD, WAA, WCC, WDD,
                &ac1, &res, &flop_count);
      }
      TIMING_stop("PCR_RB", flop_count);
    }
    flop += flop_count;
    
    
    if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;
    
    if ( converge_check ) {
      if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;
      
      res *= res_normal;
      res = sqrt(res);
      Hostonly_ {
        fprintf(fph, "%6d, %13.6e\n", itr, res);
        fflush(fph);
      }
      
      if ( res < eps ) break;
    }
    
  } // Iteration
  
  return itr;
}


/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 * @param [in]     s_type ソルバーの指定
 */
int CZ::LSOR_PCR(double& res, REAL_TYPE* X, REAL_TYPE* B,
                    const int itr_max, double& flop,
                    int s_type,
                    bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int gc = GUIDE;
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int n = ked - kst + 1;
  int pn;
  
  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  
  
  for (itr=1; itr<=itr_max; itr++)
  {
    flop_count = 0.0;
    res = 0.0;
    
    if (s_type==LS_PCR_MAF)
    {
      TIMING_start("PCR_MAF");
      pcr_maf_(size, innerFidx, &gc, &pn, X, MSK, B, xc, yc, zc,
                 WA, WC, WD, WAA, WCC, WDD,
                 &ac1, &res, &flop_count);
      TIMING_stop("PCR_MAF", flop_count);
    }
    else
    {
      TIMING_start("PCR");
      pcr_(size, innerFidx, &gc, &pn, X, MSK, B,
                WA, WC, WD, WAA, WCC, WDD,
                &ac1, &res, &flop_count);
      TIMING_stop("PCR", flop_count);
    }
    flop += flop_count;
    
    
    if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;
    
    if ( converge_check ) {
      if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;
      
      res *= res_normal;
      res = sqrt(res);
      Hostonly_ {
        fprintf(fph, "%6d, %13.6e\n", itr, res);
        fflush(fph);
      }
      
      if ( res < eps ) break;
    }
    
  } // Iteration
  
  return itr;
}

