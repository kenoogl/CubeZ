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

   switch (KindOfPrecondition)
   {
     case 0:
       TIMING_start("Blas_Copy");
       blas_copy_(xx, bb, size, &gc);
       TIMING_stop("Blas_Copy");
       break;

     case 1:
       JACOBI(res, xx, bb, lc_max, flop, false);
       break;

     case 2:
       PSOR(res, xx, bb, lc_max, flop, false);
       break;

     case 3:
       RBSOR(res, xx, bb, lc_max, flop, false);
       break;

     case 4:
       LJCB_MSD(res, xx, bb, lc_max, flop, false);
       break;
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
 * @brief Line Jacobi反復
 * @param [in,out] res    残差
 * @param [in,out] X      解ベクトル
 * @param [in]     B      RHSベクトル
 * @param [in]     itr_max 最大反復数
 * @param [in]     flop   浮動小数点演算数

 int CZ::LJacobi(double& res, REAL_TYPE* X, REAL_TYPE* B,
              const int itr_max, double& flop, bool converge_check)
  {
    int itr;
    double flop_count = 0.0;
    int gc = GUIDE;
    REAL_TYPE mat[3]={-1.0/6.0, 1.0, -1.0/6.0};

    REAL_TYPE* q;  // RHS
    REAL_TYPE* w;  // work

    if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;

    TIMING_start("BoundaryCondition");
    bc_(size, &gc, q, pitch, origin, nID);
    TIMING_stop("BoundaryCondition");


    for (itr=1; itr<=itr_max; itr++)
    {
      res = 0.0;

      /*
      TIMING_start("TDMA_rhs");
      flop_count = 0.0;
      tdma_rhs_(q, size, innerFidx, &gc, X, &flop_count);
      TIMING_stop("TDMA_rhs", flop_count);

      TIMING_start("TDMA_kernel");
      flop_count = 0.0;
      tdma_wrap_(q, size, innerFidx, &gc, w, mat, &flop_count);
      TIMING_stop("TDMA_kernel", flop_count);

      TIMING_start("LJacobi_kernel");
      flop_count = 0.0;
      tdma_sor_(X, size, innerFidx, &gc, &ac1, q, &res, &flop_count);
      TIMING_stop("LJacobi_kernel", flop_count);


      TIMING_start("LJacobi_kernel");
      flop_count = 0.0;
      tdma_ljcb_(q, size, innerFidx, &gc, X, w, mat, &ac1, &res, &flop_count);
      TIMING_stop("LJacobi_kernel", flop_count);

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

    if (q) delete [] q;
    if (w) delete [] w;

    return itr;
  }*/

  /* #################################################################
  * @brief Line SOR反復
  * @param [in,out] res    残差
  * @param [in,out] X      解ベクトル
  * @param [in]     B      RHSベクトル
  * @param [in]     itr_max 最大反復数
  * @param [in]     flop   浮動小数点演算数
  */
  int CZ::LSOR(double& res, REAL_TYPE* X, REAL_TYPE* B,
               const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     //REAL_TYPE mat[3]={-1.0/6.0, 1.0, -1.0/6.0};

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;
     a = new REAL_TYPE [size[0]+2*gc];
     b = new REAL_TYPE [size[0]+2*gc];
     c = new REAL_TYPE [size[0]+2*gc];

     for (int i=0; i<size[0]+2*gc; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=1; i<size[0]; i++) {
       a[i+gc] = -1.0/6.0;
     }
     for (int i=0; i<size[0]; i++) {
       b[i+gc] = 1.0;
     }
     for (int i=0; i<size[0]-1; i++) {
       c[i+gc] = -1.0/6.0;
     }

     TIMING_start("BoundaryCondition");
     bc_src_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_kernel");
       flop_count = 0.0;
       tdma_lsor_(q, size, innerFidx, &gc, X, w, a, b, c, &ac1, &res, &flop_count);
       TIMING_stop("LSOR_kernel", flop_count);

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

     if (q) delete [] q;
     if (w) delete [] w;
     if (a) delete [] a;
     if (b) delete [] b;
     if (c) delete [] c;

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
   int CZ::LSOR_MS(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE mat[3]={-1.0/6.0, 1.0, -1.0/6.0};

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;

     TIMING_start("BoundaryCondition");
     bc_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_MS_kernel");
       flop_count = 0.0;
       lsor_ms(q, X, w, res, flop_count);
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

     if (q) delete [] q;
     if (w) delete [] w;

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
   int CZ::LSOR_MSB(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE mat[3]={-1.0/6.0, 1.0, -1.0/6.0};

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;
     a = new REAL_TYPE [size[0]+2*gc];
     b = new REAL_TYPE [size[0]+2*gc];
     c = new REAL_TYPE [size[0]+2*gc];

     for (int i=0; i<size[0]+2*gc; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=2; i<size[0]; i++) {
       a[i+gc] = -1.0/6.0;
     }
     for (int i=0; i<size[0]; i++) {
       b[i+gc] = 1.0;
     }
     for (int i=0; i<size[0]-2; i++) {
       c[i+gc] = -1.0/6.0;
     }

     TIMING_start("BoundaryCondition");
     bc_src_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_MS_kernel");
       flop_count = 0.0;
       tdma_lsor_b_(q, size, innerFidx, &gc, X, w, a, b, c, &ac1, &res, &flop_count);
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

     if (q) delete [] q;
     if (w) delete [] w;
     if (a) delete [] a;
     if (b) delete [] b;
     if (c) delete [] c;

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
   int CZ::LSOR_MSC(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     //REAL_TYPE mat[3]={-1.0/6.0, 1.0, -1.0/6.0};

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;
     a = new REAL_TYPE [size[2]+2*gc];
     b = new REAL_TYPE [size[2]+2*gc];
     c = new REAL_TYPE [size[2]+2*gc];

     for (int i=0; i<size[2]+2*gc; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+gc-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+gc-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+gc-1] = -1.0/6.0;
     }

     TIMING_start("BoundaryCondition");
     bc_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_MS_kernel");
       flop_count = 0.0;
       tdma_lsor_c_(q, size, innerFidx, &gc, X, w, a, b, c, &ac1, &res, &flop_count);
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

     if (q) delete [] q;
     if (w) delete [] w;
     if (a) delete [] a;
     if (b) delete [] b;
     if (c) delete [] c;

     return itr;
   }

   /* #################################################################
   * @brief Line JACOBI Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   */
   int CZ::LJCB_MSD(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;
     a = new REAL_TYPE [size[2]+2*gc];
     b = new REAL_TYPE [size[2]+2*gc];
     c = new REAL_TYPE [size[2]+2*gc];

     for (int i=0; i<size[2]+2*gc; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+gc-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+gc-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+gc-1] = -1.0/6.0;
     }

     TIMING_start("BoundaryCondition");
     bc_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_MS_kernel");
       flop_count = 0.0;
       tdma_ljcb_d_(q, size, innerFidx, &gc, X, w, a, b, c, &ac1, &res, &flop_count);
       TIMING_stop("LJCB_MS_kernel", flop_count);

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

     if (q) delete [] q;
     if (w) delete [] w;
     if (a) delete [] a;
     if (b) delete [] b;
     if (c) delete [] c;

     return itr;
   }

   /* #################################################################
   * @brief Line JACOBI Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   */
   int CZ::LJCB_MSE(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;
     a = new REAL_TYPE [size[2]+2*gc];
     b = new REAL_TYPE [size[2]+2*gc];
     c = new REAL_TYPE [size[2]+2*gc];

     for (int i=0; i<size[2]+2*gc; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+gc-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+gc-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+gc-1] = -1.0/6.0;
     }

     TIMING_start("BoundaryCondition");
     bc_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_MS_kernel");
       flop_count = 0.0;
       tdma_ljcb_e_(q, size, innerFidx, &gc, X, w, a, b, c, MSK, &ac1, &res, &flop_count);
       TIMING_stop("LJCB_MS_kernel", flop_count);

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

     if (q) delete [] q;
     if (w) delete [] w;
     if (a) delete [] a;
     if (b) delete [] b;
     if (c) delete [] c;

     return itr;
   }


   /* #################################################################
   * @brief Line JACOBI Multi System
   * @param [in,out] res    残差
   * @param [in,out] X      解ベクトル
   * @param [in]     B      RHSベクトル
   * @param [in]     itr_max 最大反復数
   * @param [in]     flop   浮動小数点演算数
   */
   int CZ::LJCB_MSF(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = Alloc_Real_S3D(size)) == NULL ) return 0;
     if( (w = Alloc_Real_S3D(size)) == NULL ) return 0;
     a = new REAL_TYPE [size[2]+2*gc];
     b = new REAL_TYPE [size[2]+2*gc];
     c = new REAL_TYPE [size[2]+2*gc];

     for (int i=0; i<size[2]+2*gc; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+gc-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+gc-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+gc-1] = -1.0/6.0;
     }

     TIMING_start("BoundaryCondition");
     bc_(size, &gc, q, pitch, origin, nID);
     TIMING_stop("BoundaryCondition");


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_f0_kernel");
       flop_count = 0.0;
       ljcb_f0_(q, size, innerFidx, &gc, X, &flop_count);
       TIMING_stop("LJCB_f0_kernel", flop_count);

       TIMING_start("LJCB_f1_kernel");
       flop_count = 0.0;
       ljcb_f1_(q, size, innerFidx, &gc, w, b, c, &flop_count);
       TIMING_stop("LJCB_f1_kernel", flop_count);

       TIMING_start("LJCB_f2_kernel");
       flop_count = 0.0;
       ljcb_f2_(q, size, innerFidx, &gc, w, a, b, c, &flop_count);
       TIMING_stop("LJCB_f2_kernel", flop_count);

       TIMING_start("LJCB_f3_kernel");
       flop_count = 0.0;
       ljcb_f3_(q, size, innerFidx, &gc, w, &flop_count);
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

     if (q) delete [] q;
     if (w) delete [] w;
     if (a) delete [] a;
     if (b) delete [] b;
     if (c) delete [] c;

     return itr;
   }
