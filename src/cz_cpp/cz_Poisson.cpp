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
 * @param [in]     converge_check 0のとき、収束判定しない
 */
 int CZ::JACOBI(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;

    for (itr=1; itr<=itr_max; itr++)
    {
      res = 0.0;

      TIMING_start("JACOBI_kernel");
      flop_count = 0.0;
      jacobi_(X, size, innerFidx, &gc, cf, &ac1, B, &res, WRK, &flop_count);
      TIMING_stop("JACOBI_kernel", flop_count);
      flop += flop_count;

      if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;


      if ( converge_check ) {
        if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

        res *= res_normal;
        res = sqrt(res);
        Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

        TIMING_start("BoundaryCondition");
        bc_(size, &gc, X, pitch, origin, nID);
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
 * @param [in]     converge_check 0のとき、収束判定しない
 */
 int CZ::PSOR(double& res, REAL_TYPE* X, REAL_TYPE* B,
              const int itr_max, double& flop, bool converge_check)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;

    for (itr=1; itr<=itr_max; itr++)
    {
      res = 0.0;

      TIMING_start("SOR_kernel");
      flop_count = 0.0;
      psor_(X, size, innerFidx, &gc, cf, &ac1, B, &res, &flop_count);
      TIMING_stop("SOR_kernel", flop_count);

      if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;



      if ( converge_check ) {
        if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

        res *= res_normal;
        res = sqrt(res);

        Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

        TIMING_start("BoundaryCondition");
        bc_(size, &gc, X, pitch, origin, nID);
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
* @param [in]     converge_check 0のとき、収束判定しない
*/
int CZ::RBSOR(double& res, REAL_TYPE* X, REAL_TYPE* B,
              const int itr_max, double& flop, bool converge_check)
 {
   int itr;
   double flop_count = 0.0;
   int gc = GUIDE;

   for (itr=1; itr<=itr_max; itr++)
   {
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
     TIMING_start("SOR2SMA_kernel");
     flop_count = 0.0;
     for (int color=0; color<2; color++)
     {
       // res_p >> 反復残差の二乗和
       psor2sma_core_(X, size, innerFidx, &gc, cf, &ip, &color, &ac1, B, &res, &flop_count);
     }
     TIMING_stop("SOR2SMA_kernel", flop_count);

     if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;


     if ( converge_check ) {
       if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

       res *= res_normal;
       res = sqrt(res);
       Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

       TIMING_start("BoundaryCondition");
       bc_(size, &gc, X, pitch, origin, nID);
       TIMING_stop("BoundaryCondition");

       if ( res < eps ) break;
     }

   } // Iteration

   return itr;
 }


 // #################################################################
 double CZ::Fdot1(REAL_TYPE* x, double& flop)
 {
   double flop_count=0.0;          /// 浮動小数点演算数
   double xy = 0.0;
   int gc = GUIDE;

   TIMING_start("Dot1");
   blas_dot1_(&xy, x, size, innerFidx, &gc, &flop_count);
   TIMING_stop("Dot1", flop_count);
   flop += flop_count;

   if ( !Comm_SUM_1(&xy, "A_R_Dot") ) Exit(0);

   return xy;
 }

 // #################################################################
 double CZ::Fdot2(REAL_TYPE* x, REAL_TYPE* y, double& flop)
 {
   double flop_count=0.0;          /// 浮動小数点演算数
   double xy = 0.0;
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
                         double& flop)
 {
   int gc = GUIDE;
   double res = 0.0;
   int lc_max = 4;

   switch (pc_type)
   {
     case LS_JACOBI:
       JACOBI(res, xx, bb, lc_max, flop, false);
       break;

     case LS_PSOR:
       PSOR(res, xx, bb, lc_max, flop, false);
       break;

     case LS_SOR2SMA:
       RBSOR(res, xx, bb, lc_max, flop, false);
       break;

     case LS_LSOR_A:
       LSOR_A(res, xx, bb, lc_max, flop, false);
       break;

     case LS_LSOR_E:
       LSOR_E(res, xx, bb, lc_max, flop, false);
       break;

     case LS_LSOR_F:
       LSOR_F(res, xx, bb, lc_max, flop, false);
       break;

     case LS_LJCB_E:
       LJCB_E(res, xx, bb, lc_max, flop, false);
       break;

     case LS_LSOR_J:
     case LS_LSOR_J4:
     case LS_LSOR_K:
     case LS_LSOR_K2:
       LSOR_J(res, xx, bb, lc_max, flop, false);
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
 int CZ::PBiCGSTAB(double& res,
                   REAL_TYPE* X,
                   REAL_TYPE* B,
                   double& flop)
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
   blas_calc_rk_(pcg_r, X, B, size, innerFidx, &gc, cf, &flop_count);
   TIMING_stop("Blas_Residual", flop_count);
   flop += flop_count;

   if ( !Comm_S(pcg_r, 1, "Comm_Poisson") ) return 0;

   TIMING_start("Blas_Copy");
   blas_copy_(pcg_r0, pcg_r, size, &gc);
   TIMING_stop("Blas_Copy");

   double rho_old = 1.0;
   double alpha = 0.0;
   double omega  = 1.0;
   double r_omega = -omega;

   for (itr=1; itr<ItrMax; itr++)
   {
     flop_count = 0.0;
     double rho = Fdot2(pcg_r, pcg_r0, flop_count);
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
       double beta = rho / rho_old * alpha / omega;

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
     Preconditioner(pcg_p_, pcg_p, flop_count);
     flop += flop_count;

     TIMING_start("Blas_AX");
     flop_count = 0.0;
     blas_calc_ax_(pcg_q, pcg_p_, size, innerFidx, &gc, cf, &flop_count);
     TIMING_stop("Blas_AX", flop_count);
     flop += flop_count;

     flop_count = 0.0;
     alpha = rho / Fdot2(pcg_q, pcg_r0, flop_count);
     flop += flop_count;

     double r_alpha = -alpha;
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
     Preconditioner(pcg_s_, pcg_s, flop_count);
     flop += flop_count;

     TIMING_start("Blas_AX");
     flop_count = 0.0;
     blas_calc_ax_(pcg_t_, pcg_s_, size, innerFidx, &gc, cf, &flop_count);
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
   bc_(size, &gc, X, pitch, origin, nID);
   TIMING_stop("BoundaryCondition");

   return itr;
 }



   /* #################################################################
   * @brief Line SOR Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   * @note tdma_sor_bから、TDMAをkループに変更、kは最外側のまま
   */
   int CZ::LSOR_E(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* c;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     memcpy(q, B, sizeof(REAL_TYPE)*(
       (size[0]+2*GUIDE)*(size[1]+2*GUIDE)*(size[2]+2*GUIDE)
     ));

     a = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=0; i<size[2]+2*GUIDE; i++) {
       a[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }



     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_MS_kernel");
       flop_count = 0.0;
       tdma_lsor_e_(q, size, innerFidx, &gc, X, w, a, c, B, &ac1, &res, &flop_count);
       TIMING_stop("LSOR_MS_kernel", flop_count);

       if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

       if ( converge_check ) {
         if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

         res *= res_normal;
         res = sqrt(res);
         Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

         // BCはq[]に与えられているので、不要
         if ( res < eps ) break;
       }

     } // Iteration

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(c);

     return itr;
   }

   /* #################################################################
   * @brief Line SOR Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   * @note tdma_sor_bから、TDMAをkループに変更、kは最外側のまま
   */
   int CZ::LSOR_F(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;
     int kst = innerFidx[K_minus];

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* c;
     REAL_TYPE* e;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     memcpy(q, B, sizeof(REAL_TYPE)*(
       (size[0]+2*GUIDE)*(size[1]+2*GUIDE)*(size[2]+2*GUIDE)
     ));

     a = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);
     e = czAllocR(size[2]+2*GUIDE, var_type);
     w = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }

     TIMING_start("LSOR_LU_decomp");
     flop_count = 0.0;
     tdma_pre(&a[kst+GUIDE-1],
              &c[kst+GUIDE-1],
              &e[kst+GUIDE-1],
              &w[kst+GUIDE-1],
              flop_count);
     TIMING_stop("LSOR_LU_decomp", flop_count);


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_MS_kernel");
       flop_count = 0.0;
       tdma_lsor_f_(q, size, innerFidx, &gc, X, w, e, a, B, &ac1, &res, &flop_count);
       TIMING_stop("LSOR_MS_kernel", flop_count);

       if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

       if ( converge_check ) {
         if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

         res *= res_normal;
         res = sqrt(res);
         Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

         // BCはq[]に与えられているので、不要
         if ( res < eps ) break;
       }

     } // Iteration

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(c);

     return itr;
   }


   /* #################################################################
   * @brief Line JACOBI Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   * @note LUの係数の再利用
   */
   int CZ::LJCB_E(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;
     int kst = innerFidx[K_minus];

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* c;
     REAL_TYPE* e;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     blas_copy_(q, B, size, &gc);

     a = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);
     e = czAllocR(size[2]+2*GUIDE, var_type);
     w = czAllocR(size[2]+2*GUIDE, var_type);


     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }

     TIMING_start("LSOR_LU_decomp");
     flop_count = 0.0;
     tdma_pre(&a[kst+GUIDE-1],
              &c[kst+GUIDE-1],
              &e[kst+GUIDE-1],
              &w[kst+GUIDE-1],
              flop_count);
     TIMING_stop("LSOR_LU_decomp", flop_count);



     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_f0_kernel");
       flop_count = 0.0;
       ljcb_f0_(q, size, innerFidx, &gc, X, B, &flop_count);
       TIMING_stop("LJCB_f0_kernel", flop_count);

       TIMING_start("LJCB_f1_kernel");
       flop_count = 0.0;
       ljcb_g1_(q, size, innerFidx, &gc, B, &flop_count);
       TIMING_stop("LJCB_f1_kernel", flop_count);

       TIMING_start("LJCB_f2_kernel");
       flop_count = 0.0;
       ljcb_g2_(q, size, innerFidx, &gc, e, a, &flop_count);
       TIMING_stop("LJCB_f2_kernel", flop_count);

       TIMING_start("LJCB_f3_kernel");
       flop_count = 0.0;
       ljcb_g3_(q, size, innerFidx, &gc, w, &flop_count);
       TIMING_stop("LJCB_f3_kernel", flop_count);

       TIMING_start("LJCB_f4_kernel");
       flop_count = 0.0;
       ljcb_f4_(q, size, innerFidx, &gc, X, &ac1, &res, &flop_count);
       TIMING_stop("LJCB_f4_kernel", flop_count);


       if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

       if ( converge_check ) {
         if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

         res *= res_normal;
         res = sqrt(res);
         Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

         // BCはq[]に与えられているので、不要
         if ( res < eps ) break;
       }

     } // Iteration

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(c);
     czDelete(e);

     return itr;
   }



   /* #################################################################
   * @brief Line SOR Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   */
   int CZ::LSOR_J(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* c;
     REAL_TYPE* e;
     REAL_TYPE* q2;

     int kst = innerFidx[K_minus];

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (q2= czAllocR_S3D(size,var_type)) == NULL ) return 0;

     memcpy(q, B, sizeof(REAL_TYPE)*(
       (size[0]+2*GUIDE)*(size[1]+2*GUIDE)*(size[2]+2*GUIDE)
     ));

     a = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);
     e = czAllocR(size[2]+2*GUIDE, var_type);
     w = czAllocR(size[2]+2*GUIDE, var_type);

     double* resD = czAllocR(CL_SZ*numThreads, flop_count);
     double* rd = czAllocR(CL_SZ, flop_count);


     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }

     TIMING_start("LSOR_LU_decomp");
     flop_count = 0.0;
     tdma_pre(&a[kst+GUIDE-1],
              &c[kst+GUIDE-1],
              &e[kst+GUIDE-1],
              &w[kst+GUIDE-1],
              flop_count);
     TIMING_stop("LSOR_LU_decomp", flop_count);


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       flop_count = 0.0;
       if (ls_type==LS_LSOR_J) {
         lsor_j(q, X, w, a, e, B, MSK, q2, res, flop_count);
       }
       else if (ls_type==LS_LSOR_J4) {
         lsor_j4(q, X, w, a, e, B, MSK, q2, res, flop_count);
       }
       else if (ls_type==LS_LSOR_K) {
         lsor_k(q, X, w, a, e, B, MSK, q2, res, flop_count);
       }
       else if (ls_type==LS_LSOR_K2) {
         lsor_k2(q, X, w, a, e, B, MSK, q2, res, flop_count);
       }
       else if (ls_type==LS_LSOR_K3) {
         res = lsor_k3(q, X, w, a, e, B, MSK, q2, resD, rd, flop_count);
       }


       if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

       if ( converge_check ) {
         if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

         res *= res_normal;
         res = sqrt(res);
         Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

         if ( res < eps ) break;
       }

     } // Iteration

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(c);
     czDelete(e);
     czDelete(q2);
     czDelete(resD);
     czDelete(rd);

     return itr;
   }


/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 */
int CZ::LSOR_P1(double& res, REAL_TYPE* X, REAL_TYPE* B,
               const int itr_max, double& flop, bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int gc = GUIDE;
  REAL_TYPE var_type=0;
  
  REAL_TYPE* a;
  REAL_TYPE* c;
  REAL_TYPE* d;
  REAL_TYPE* a1;
  REAL_TYPE* c1;
  REAL_TYPE* d1;
  
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int n = ked - kst + 1;
  int pn;
  
  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }

  a  = czAllocR_S3D(size, var_type);
  c  = czAllocR_S3D(size, var_type);
  d  = czAllocR_S3D(size, var_type);
  a1 = czAllocR_S3D(size, var_type);
  c1 = czAllocR_S3D(size, var_type);
  d1 = czAllocR_S3D(size, var_type);
  
  
  for (itr=1; itr<=itr_max; itr++)
  {
    flop_count = 0.0;
    res = 0.0;
    TIMING_start("LSOR_PCR");
    
    lsor_pcr_kij_(size, innerFidx, &gc, &pn, X, a, c, d, a1, c1, d1, MSK, B, &ac1, &res, &flop_count);

    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
  
  czDelete(a);
  czDelete(c);
  czDelete(d);
  czDelete(a1);
  czDelete(c1);
  czDelete(d1);
  
  return itr;
}


/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 */
int CZ::LSOR_P2(double& res, REAL_TYPE* X, REAL_TYPE* B,
               const int itr_max, double& flop, bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int gc = GUIDE;
  REAL_TYPE var_type=0;
  
  REAL_TYPE* a;
  REAL_TYPE* c;
  REAL_TYPE* d;
  REAL_TYPE* a1;
  REAL_TYPE* c1;
  REAL_TYPE* d1;
  
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int n = ked - kst + 1;
  int pn;
  
  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  
  a  = czAllocR_S3D(size, var_type);
  c  = czAllocR_S3D(size, var_type);
  d  = czAllocR_S3D(size, var_type);
  a1 = czAllocR_S3D(size, var_type);
  c1 = czAllocR_S3D(size, var_type);
  d1 = czAllocR_S3D(size, var_type);
  
  
  for (itr=1; itr<=itr_max; itr++)
  {
    flop_count = 0.0;
    res = 0.0;
    TIMING_start("LSOR_PCR");
    
    lsor_pcr_kij2_(size, innerFidx, &gc, &pn, X, a, c, d, a1, c1, d1, MSK, B, &ac1, &res, &flop_count);
    
    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
  
  czDelete(a);
  czDelete(c);
  czDelete(d);
  czDelete(a1);
  czDelete(c1);
  czDelete(d1);
  
  return itr;
}

/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 */
int CZ::LSOR_P3(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int gc = GUIDE;
  REAL_TYPE var_type=0;
  
  REAL_TYPE* a;
  REAL_TYPE* c;
  REAL_TYPE* d;
  REAL_TYPE* a1;
  REAL_TYPE* c1;
  REAL_TYPE* d1;
  
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int n = ked - kst + 1;
  int pn;
  
  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  
  a  = czAllocR_S3D(size, var_type);
  c  = czAllocR_S3D(size, var_type);
  d  = czAllocR_S3D(size, var_type);
  a1 = czAllocR_S3D(size, var_type);
  c1 = czAllocR_S3D(size, var_type);
  d1 = czAllocR_S3D(size, var_type);
  
  
  for (itr=1; itr<=itr_max; itr++)
  {
    flop_count = 0.0;
    res = 0.0;
    TIMING_start("LSOR_PCR");
    
    lsor_pcr_kij3_(size, innerFidx, &gc, &pn, X, a, c, d, a1, c1, d1, MSK, B, &ac1, &res, &flop_count);
    
    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
  
  czDelete(a);
  czDelete(c);
  czDelete(d);
  czDelete(a1);
  czDelete(c1);
  czDelete(d1);
  
  return itr;
}

/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 */
int CZ::LSOR_P4(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int gc = GUIDE;
  REAL_TYPE var_type=0;
  
  REAL_TYPE* a;
  REAL_TYPE* c;
  REAL_TYPE* d;
  REAL_TYPE* a1;
  REAL_TYPE* c1;
  REAL_TYPE* d1;
  
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int n = ked - kst + 1;
  int pn;
  
  // Nを超える最小の2べき数の乗数 pn
  if ( -1 == (pn=getNumStage(n))) {
    printf("error : number of stage\n");
    exit(0);
  }
  
  a  = czAllocR_S3D(size, var_type);
  c  = czAllocR_S3D(size, var_type);
  d  = czAllocR_S3D(size, var_type);
  a1 = czAllocR_S3D(size, var_type);
  c1 = czAllocR_S3D(size, var_type);
  d1 = czAllocR_S3D(size, var_type);
  
  
  for (itr=1; itr<=itr_max; itr++)
  {
    flop_count = 0.0;
    res = 0.0;
    TIMING_start("LSOR_PCR");
    
    lsor_pcr_kij4_(size, innerFidx, &gc, &pn, X, a, c, d, a1, c1, d1, MSK, B, &ac1, &res, &flop_count);
    
    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
  
  czDelete(a);
  czDelete(c);
  czDelete(d);
  czDelete(a1);
  czDelete(c1);
  czDelete(d1);
  
  return itr;
}

/* #################################################################
 * @brief Line SOR PCR
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数
 */
int CZ::LSOR_P5(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
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
    TIMING_start("LSOR_PCR");
    
    lsor_pcr_kij5_(size, innerFidx, &gc, &pn, X, MSK, B, &ac1, &res, &flop_count);
    
    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
 */
int CZ::LSOR_P6(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
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
    TIMING_start("LSOR_PCR");
    
    lsor_pcr_kij6_(size, innerFidx, &gc, &pn, X, MSK, B, &ac1, &res, &flop_count);
    
    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
 */
int CZ::LSOR_P7(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
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
    
    TIMING_start("LSOR_PCR");
    for (int color=0; color<2; color++)
    {
      lsor_pcr_kij7_(size, innerFidx, &gc, &pn, &ip, &color, X, MSK, B, &ac1, &res, &flop_count);
    }
    TIMING_stop("LSOR_PCR", flop_count);
    
    
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
