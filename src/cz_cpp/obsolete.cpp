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

/* #################################################################
* @brief Line SOR反復
* @param [in,out] res    残差
* @param [in,out] X      解ベクトル
* @param [in]     B      RHSベクトル
* @param [in]     itr_max 最大反復数
* @param [in]     flop   浮動小数点演算数
* @note LSORのリファレンス実装
        カーネルは全てFortran
        内側ループはi
*/
int CZ::LSOR_A(double& res, REAL_TYPE* X, REAL_TYPE* B,
             const int itr_max, double& flop, bool converge_check)
 {
   int itr;
   double flop_count = 0.0;
   int gc = GUIDE;
   REAL_TYPE var_type=0;

   REAL_TYPE* q;  // RHS
   REAL_TYPE* w;  // work
   REAL_TYPE* a;
   REAL_TYPE* b;
   REAL_TYPE* c;

   if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
   if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

   blas_copy_(q, B, size, &gc);

   a = czAllocR(size[0]+2*GUIDE, var_type);
   b = czAllocR(size[0]+2*GUIDE, var_type);
   c = czAllocR(size[0]+2*GUIDE, var_type);

   for (int i=0; i<size[0]+2*GUIDE; i++) {
     a[i] = 0.0;
     b[i] = 0.0;
     c[i] = 0.0;
   }
   for (int i=1; i<size[0]; i++) {
     a[i+GUIDE] = -1.0/6.0;
   }
   for (int i=0; i<size[0]; i++) {
     b[i+GUIDE] = 1.0;
   }
   for (int i=0; i<size[0]-1; i++) {
     c[i+GUIDE] = -1.0/6.0;
   }



   for (itr=1; itr<=itr_max; itr++)
   {
     res = 0.0;

     TIMING_start("LSOR_kernel");
     flop_count = 0.0;
     tdma_lsor_a_(q, size, innerFidx, &gc, X, w, a, b, c, B, &ac1, &res, &flop_count);
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

   czDelete(q);
   czDelete(w);
   czDelete(a);
   czDelete(b);
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
    * @note カーネルをCで実装
            TDMAの係数がスカラー固定なのでB/F異なる
            内側ループはi
    */
    int CZ::LSOR_C(double& res, REAL_TYPE* X, REAL_TYPE* B,
                 const int itr_max, double& flop, bool converge_check)
    {
      int itr;
      double flop_count = 0.0;
      int gc = GUIDE;
      REAL_TYPE mat[3]={-1.0/6.0, 1.0, -1.0/6.0};
      REAL_TYPE var_type=0;


      REAL_TYPE* q;  // RHS
      REAL_TYPE* w;  // work

      if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
      if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

      blas_copy_(q, B, size, &gc);


      for (itr=1; itr<=itr_max; itr++)
      {
        res = 0.0;

        TIMING_start("LSOR_MS_kernel");
        flop_count = 0.0;
        lsor_ms(q, X, w, B, res, flop_count);
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
   int CZ::LSOR_D(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     blas_copy_(q, B, size, &gc);

     a = czAllocR(size[2]+2*GUIDE, var_type);
     b = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=0; i<size[2]+2*GUIDE; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+GUIDE-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }



     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LSOR_MS_kernel");
       flop_count = 0.0;
       tdma_lsor_d_(q, size, innerFidx, &gc, X, w, a, b, c, B, &ac1, &res, &flop_count);
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
     czDelete(b);
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
   * @note TDMAはk-loop(外側)で、内側のiループでベクトル化
   */
   int CZ::LJCB_A(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     blas_copy_(q, B, size, &gc);

     a = czAllocR(size[2]+2*GUIDE, var_type);
     b = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=0; i<size[2]+2*GUIDE; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+GUIDE-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }



     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_MS_kernel");
       flop_count = 0.0;
       tdma_ljcb_a_(q, size, innerFidx, &gc, X, w, a, b, c, B, &ac1, &res, &flop_count);
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

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(b);
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
   * @note tdma_ljcb_a_のijループを一重化
   */
   int CZ::LJCB_B(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     blas_copy_(q, B, size, &gc);

     a = czAllocR(size[2]+2*GUIDE, var_type);
     b = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=0; i<size[2]+2*GUIDE; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+GUIDE-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }


     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_MS_kernel");
       flop_count = 0.0;
       tdma_ljcb_b_(q, size, innerFidx, &gc, X, w, a, b, c, B, MSK, &ac1, &res, &flop_count);
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

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(b);
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
   * @note tdma_ljcb_a_の各カーネルを独立サブルーチンにして計時
   */
   int CZ::LJCB_C(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     blas_copy_(q, B, size, &gc);

     a = czAllocR(size[2]+2*GUIDE, var_type);
     b = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=0; i<size[2]+2*GUIDE; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+GUIDE-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }



     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_f0_kernel");
       flop_count = 0.0;
       ljcb_f0_(q, size, innerFidx, &gc, X, B, &flop_count);
       TIMING_stop("LJCB_f0_kernel", flop_count);

       TIMING_start("LJCB_f1_kernel");
       flop_count = 0.0;
       ljcb_f1_(q, size, innerFidx, &gc, w, b, c, B, &flop_count);
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

     czDelete(q);
     czDelete(w);
     czDelete(a);
     czDelete(b);
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
   * @note tdma_ljcb_a_の各k-loopを最内側に
   */
   int CZ::LJCB_D(double& res, REAL_TYPE* X, REAL_TYPE* B,
                const int itr_max, double& flop, bool converge_check)
   {
     int itr;
     double flop_count = 0.0;
     int gc = GUIDE;
     REAL_TYPE var_type=0;

     REAL_TYPE* q;  // RHS
     REAL_TYPE* w;  // work
     REAL_TYPE* a;
     REAL_TYPE* b;
     REAL_TYPE* c;

     if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
     if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

     //blas_copy_(q, B, size, &gc);
     memcpy(q, B, sizeof(REAL_TYPE)*(
       (size[0]+2*GUIDE)*(size[1]+2*GUIDE)*(size[2]+2*GUIDE)
     ));

     a = czAllocR(size[2]+2*GUIDE, var_type);
     b = czAllocR(size[2]+2*GUIDE, var_type);
     c = czAllocR(size[2]+2*GUIDE, var_type);

     for (int i=0; i<size[2]+2*GUIDE; i++) {
       a[i] = 0.0;
       b[i] = 0.0;
       c[i] = 0.0;
     }
     for (int i=3; i<=size[2]-1; i++) {
       a[i+GUIDE-1] = -1.0/6.0;
     }
     for (int i=2; i<=size[2]-1; i++) {
       b[i+GUIDE-1] = 1.0;
     }
     for (int i=2; i<=size[2]-2; i++) {
       c[i+GUIDE-1] = -1.0/6.0;
     }



     for (itr=1; itr<=itr_max; itr++)
     {
       res = 0.0;

       TIMING_start("LJCB_f0_kernel");
       flop_count = 0.0;
       ljcb_f0t_(q, size, innerFidx, &gc, X, B, &flop_count);
       TIMING_stop("LJCB_f0_kernel", flop_count);

       TIMING_start("LJCB_f1_kernel");
       flop_count = 0.0;
       ljcb_f1t_(q, size, innerFidx, &gc, w, b, c, B, &flop_count);
       TIMING_stop("LJCB_f1_kernel", flop_count);

       TIMING_start("LJCB_f2_kernel");
       flop_count = 0.0;
       ljcb_f2t_(q, size, innerFidx, &gc, w, a, b, c, &flop_count);
       TIMING_stop("LJCB_f2_kernel", flop_count);

       TIMING_start("LJCB_f3_kernel");
       flop_count = 0.0;
       ljcb_f3t_(q, size, innerFidx, &gc, w, &flop_count);
       TIMING_stop("LJCB_f3_kernel", flop_count);

       TIMING_start("LJCB_f4_kernel");
       flop_count = 0.0;
       ljcb_f4t_(q, size, innerFidx, &gc, X, &ac1, &res, &flop_count);
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
     czDelete(b);
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
*/
int CZ::LSOR_SIMD(double& res, REAL_TYPE* X, REAL_TYPE* B,
             const int itr_max, double& flop, bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int gc = GUIDE;
  REAL_TYPE var_type=0;

  REAL_TYPE* q;  // RHS
  REAL_TYPE* w;  // work
  REAL_TYPE* a;
  REAL_TYPE* b;
  REAL_TYPE* c;

  if( (q = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  if( (w = czAllocR_S3D(size,var_type)) == NULL ) return 0;

  memcpy(q, B, sizeof(REAL_TYPE)*(
    (size[0]+2*GUIDE)*(size[1]+2*GUIDE)*(size[2]+2*GUIDE)
  ));

  a = czAllocR(size[2]+2*GUIDE, var_type);
  b = czAllocR(size[2]+2*GUIDE, var_type);
  c = czAllocR(size[2]+2*GUIDE, var_type);

  for (int i=0; i<size[2]+2*GUIDE; i++) {
    a[i] = 0.0;
    b[i] = 0.0;
    c[i] = 0.0;
  }
  for (int i=3; i<=size[2]-1; i++) {
    a[i+GUIDE-1] = -1.0/6.0;
  }
  for (int i=2; i<=size[2]-1; i++) {
    b[i+GUIDE-1] = 1.0;
  }
  for (int i=2; i<=size[2]-2; i++) {
    c[i+GUIDE-1] = -1.0/6.0;
  }


  TIMING_start("LSOR_simd_Itr");
  for (itr=1; itr<=itr_max; itr++)
  {
    res = 0.0;

    TIMING_start("LSOR_simd_kernel");
    flop_count = 0.0;
    lsor_simd(q, X, w, a, b, c, B, res, flop_count);
    //lsor_simd2(q, X, w, a, c, B, MSK, res, flop_count);
    //lsor_simd3(q, X, w, a, c, B, MSK, res, flop_count);
    TIMING_stop("LSOR_simd_kernel", flop_count);

    if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

    if ( converge_check ) {
      if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

      res *= res_normal;
      res = sqrt(res);
      Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

      if ( res < eps ) break;
    }

  } // Iteration
  TIMING_stop("LSOR_simd_Itr");

  czDelete(q);
  czDelete(w);
  czDelete(a);
  czDelete(b);
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
*/
int CZ::LSOR_SIMD2(double& res, REAL_TYPE* X, REAL_TYPE* B,
             const int itr_max, double& flop, bool converge_check)
{
  int itr;
  double flop_count = 0.0;
  int gc = GUIDE;
  REAL_TYPE var_type=0;

  REAL_TYPE* q;  // RHS
  REAL_TYPE* q2;
  REAL_TYPE* w;  // work
  REAL_TYPE* a;
  REAL_TYPE* c;
  REAL_TYPE* e;

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

  for (int i=0; i<size[2]+2*GUIDE; i++) {
    a[i] = 0.0;
    c[i] = 0.0;
    e[i] = 0.0;
    w[i] = 0.0;
  }
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


  TIMING_start("LSOR_simd_Itr");
  for (itr=1; itr<=itr_max; itr++)
  {
    res = 0.0;

    TIMING_start("LSOR_simd_kernel");
    flop_count = 0.0;
    lsor_simd4(q, X, w, a, e, B, MSK, res, flop_count);
    //lsor_simd5(q, X, w, a, e, B, MSK, q2, res, flop_count);
    //lsor_simd6(q, X, w, a, e, B, MSK, q2, res, flop_count);
    TIMING_stop("LSOR_simd_kernel", flop_count);

    if ( !Comm_S(X, 1, "Comm_Poisson") ) return 0;

    if ( converge_check ) {
      if ( !Comm_SUM_1(&res, "Comm_Res_Poisson") ) return 0;

      res *= res_normal;
      res = sqrt(res);
      Hostonly_ fprintf(fph, "%6d, %13.6e\n", itr, res);

      if ( res < eps ) break;
    }

  } // Iteration
  TIMING_stop("LSOR_simd_Itr");

  czDelete(q);
  czDelete(w);
  czDelete(a);
  czDelete(c);
  czDelete(e);
  czDelete(q2);

  return itr;
}



void CZ::ms_rhs8v(const int* ia,
                  const int* ja,
                  const int kst,
                  const int ked,
                  REAL_TYPE* d,
                  REAL_TYPE* x,
                  REAL_TYPE* rhs,
                  REAL_TYPE* msk,
                  double& flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(rhs, ALIGN);
  __assume_aligned(msk, ALIGN);

  flop += 48.0*(double)(ked-kst+2);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  REAL_TYPE r = 1.0/6.0;

  TIMING_start("LSOR_RHS_Body");
  #pragma vector always
  #pragma ivdep
  #pragma unroll(2)
  for (int k=kst-1; k<ked; k++) {
    size_t m0 = _IDX_S3D(k,ia[0],  ja[0]  ,NK, NI, GUIDE);
    d[m0] =(( x[_IDX_S3D(k,ia[0]-1,ja[0]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[0]+1,ja[0]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[0]  ,ja[0]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[0]  ,ja[0]+1,NK, NI, GUIDE)]
          ) * r + rhs[m0]
          ) *     msk[m0];

    size_t m1 = _IDX_S3D(k,ia[1]  ,ja[1],  NK, NI, GUIDE);
    d[m1] =(( x[_IDX_S3D(k,ia[1]-1,ja[1]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[1]+1,ja[1]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[1]  ,ja[1]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[1]  ,ja[1]+1,NK, NI, GUIDE)]
          ) * r + rhs[m1]
          ) *     msk[m1];

    size_t m2 = _IDX_S3D(k,ia[2]  ,ja[2]  ,NK, NI, GUIDE);
    d[m2] =(( x[_IDX_S3D(k,ia[2]-1,ja[2]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[2]+1,ja[2]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[2]  ,ja[2]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[2]  ,ja[2]+1,NK, NI, GUIDE)]
          ) * r + rhs[m2]
          ) *     msk[m2];

    size_t m3 = _IDX_S3D(k,ia[3]  ,ja[3]  ,NK, NI, GUIDE);
    d[m3] =(( x[_IDX_S3D(k,ia[3]-1,ja[3]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[3]+1,ja[3]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[3]  ,ja[3]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[3]  ,ja[3]+1,NK, NI, GUIDE)]
          ) * r + rhs[m3]
          ) *     msk[m3];

    size_t m4 = _IDX_S3D(k,ia[4]  ,ja[4]  ,NK, NI, GUIDE);
    d[m4] =(( x[_IDX_S3D(k,ia[4]-1,ja[4]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[4]+1,ja[4]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[4]  ,ja[4]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[4]  ,ja[4]+1,NK, NI, GUIDE)]
          ) * r + rhs[m4]
          ) *     msk[m4];

    size_t m5 = _IDX_S3D(k,ia[5]  ,ja[5]  ,NK, NI, GUIDE);
    d[m5] =(( x[_IDX_S3D(k,ia[5]-1,ja[5]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[5]+1,ja[5]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[5]  ,ja[5]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[5]  ,ja[5]+1,NK, NI, GUIDE)]
          ) * r + rhs[m5]
          ) *     msk[m5];

    size_t m6 = _IDX_S3D(k,ia[6]  ,ja[6]  ,NK, NI, GUIDE);
    d[m6] =(( x[_IDX_S3D(k,ia[6]-1,ja[6]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[6]+1,ja[6]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[6]  ,ja[6]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[6]  ,ja[6]+1,NK, NI, GUIDE)]
          ) * r + rhs[m6]
          ) *     msk[m6];

    size_t m7 = _IDX_S3D(k,ia[7]  ,ja[7]  ,NK, NI, GUIDE);
    d[m7] =(( x[_IDX_S3D(k,ia[7]-1,ja[7]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[7]+1,ja[7]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[7]  ,ja[7]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[7]  ,ja[7]+1,NK, NI, GUIDE)]
          ) * r + rhs[m7]
          ) *     msk[m7];
  }
  TIMING_stop("LSOR_RHS_Body", 48.0*(double)(ked-kst+2));
}


void CZ::ms_rhs4v(const int* ia,
                  const int* ja,
                  const int kst,
                  const int ked,
                  REAL_TYPE* d,
                  REAL_TYPE* x,
                  REAL_TYPE* rhs,
                  REAL_TYPE* msk,
                  double& flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(rhs, ALIGN);
  __assume_aligned(msk, ALIGN);

  flop += 24.0*(double)(ked-kst+2);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  REAL_TYPE r = 1.0/6.0;

  TIMING_start("LSOR_RHS_Body");
  #pragma vector always
  #pragma ivdep
  for (int k=kst-1; k<ked; k++) {
    size_t m0 = _IDX_S3D(k,ia[0],  ja[0]  ,NK, NI, GUIDE);
    d[m0] =(( x[_IDX_S3D(k,ia[0]-1,ja[0]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[0]+1,ja[0]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[0]  ,ja[0]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[0]  ,ja[0]+1,NK, NI, GUIDE)]
          ) * r + rhs[m0]
          ) *     msk[m0];

    size_t m1 = _IDX_S3D(k,ia[1]  ,ja[1],  NK, NI, GUIDE);
    d[m1] =(( x[_IDX_S3D(k,ia[1]-1,ja[1]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[1]+1,ja[1]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[1]  ,ja[1]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[1]  ,ja[1]+1,NK, NI, GUIDE)]
          ) * r + rhs[m1]
          ) *     msk[m1];

    size_t m2 = _IDX_S3D(k,ia[2]  ,ja[2]  ,NK, NI, GUIDE);
    d[m2] =(( x[_IDX_S3D(k,ia[2]-1,ja[2]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[2]+1,ja[2]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[2]  ,ja[2]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[2]  ,ja[2]+1,NK, NI, GUIDE)]
          ) * r + rhs[m2]
          ) *     msk[m2];

    size_t m3 = _IDX_S3D(k,ia[3]  ,ja[3]  ,NK, NI, GUIDE);
    d[m3] =(( x[_IDX_S3D(k,ia[3]-1,ja[3]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[3]+1,ja[3]  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[3]  ,ja[3]-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,ia[3]  ,ja[3]+1,NK, NI, GUIDE)]
          ) * r + rhs[m3]
          ) *     msk[m3];
  }
  TIMING_stop("LSOR_RHS_Body", 24.0*(double)(ked-kst+2));
}


void CZ::ms_rhs4(const int i,
                 const int j,
                 const int kst,
                 const int ked,
                 REAL_TYPE* d,
                 REAL_TYPE* x,
                 REAL_TYPE* rhs,
                 REAL_TYPE* msk,
                 double& flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(rhs, ALIGN);
  __assume_aligned(msk, ALIGN);

  flop += 24.0*(double)(ked-kst+2);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  REAL_TYPE r = 1.0/6.0;
  size_t    m0, m1, m2, m3;

  TIMING_start("LSOR_RHS_Body");
  #pragma vector always
  #pragma ivdep
  for (int k=kst-1; k<ked; k++) {
    m0 = _IDX_S3D(k,i,j,NK,NI,GUIDE);
    d[m0] =(( x[_IDX_S3D(k,i-1,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+1,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i  ,j-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i  ,j+1,NK, NI, GUIDE)]
          ) * r + rhs[m0]
          ) *     msk[m0];
    m1 = _IDX_S3D(k,i+1,j,NK,NI,GUIDE);
    d[m1] =(( x[_IDX_S3D(k,i  ,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+2,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+1,j-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+1,j+1,NK, NI, GUIDE)]
          ) * r + rhs[m1]
          ) *     msk[m1];

    m2 = _IDX_S3D(k,i+2,j,NK,NI,GUIDE);
    d[m2] =(( x[_IDX_S3D(k,i+1,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+3,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+2,j-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+2,j+1,NK, NI, GUIDE)]
          ) * r + rhs[m2]
          ) *     msk[m2];

    m3 = _IDX_S3D(k,i+3,j,NK,NI,GUIDE);
    d[m3] =(( x[_IDX_S3D(k,i+2,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+4,j  ,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+3,j-1,NK, NI, GUIDE)]
            + x[_IDX_S3D(k,i+3,j+1,NK, NI, GUIDE)]
          ) * r + rhs[m3]
          ) *     msk[m3];
  }
  TIMING_stop("LSOR_RHS_Body", 24.0*(double)(ked-kst+2));
}
