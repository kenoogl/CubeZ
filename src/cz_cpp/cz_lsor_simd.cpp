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
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 */
void CZ::lsor_simd(REAL_TYPE* d,
                   REAL_TYPE* x,
                   REAL_TYPE* w,
                   REAL_TYPE* a,
                   REAL_TYPE* b,
                   REAL_TYPE* c,
                   REAL_TYPE* rhs,
                   double &res,
                   double &flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(w, ALIGN);
  __assume_aligned(a, ALIGN);
  __assume_aligned(b, ALIGN);
  __assume_aligned(c, ALIGN);
  __assume_aligned(rhs, ALIGN);


  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int ist = innerFidx[I_minus];
  int ied = innerFidx[I_plus];
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  int nn = ked - kst + 1;

  REAL_TYPE r = 1.0/6.0;
  REAL_TYPE omg = ac1;
  REAL_TYPE pp, dp, pn;

  double f1 = (double)(nn)*5.0;
  double f2 = (double)(nn)*16.0 + 16.0;
  double f3 = (double)(nn)*5.0;

  flop += (double)( (ied-ist+1)*(jed-jst+1)*(f1+f2+f3+4.0) );
  res = 0.0;


  #pragma omp parallel for private(pp,dp,pn) reduction(+:res)
  for (int j=jst-1; j<jed; j++) {
  for (int i=ist-1; i<ied; i++) {

    TIMING_start("LSOR_RHS_Body");
    #pragma vector always
    #pragma ivdep
    for (int k=kst-1; k<ked; k++) {
      size_t m = _IDX_S3D(k,i,j,NK,NI,GUIDE);
      d[m] = ( x[_IDX_S3D(k,i-1,j  ,NK, NI, GUIDE)]
            +  x[_IDX_S3D(k,i+1,j  ,NK, NI, GUIDE)]
            +  x[_IDX_S3D(k,i  ,j-1,NK, NI, GUIDE)]
            +  x[_IDX_S3D(k,i  ,j+1,NK, NI, GUIDE)]
          ) * r + rhs[m];
    }
    TIMING_stop("LSOR_RHS_Body", f1);

    TIMING_start("LSOR_TDMA_BC");
    d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)] +=  rhs[_IDX_S3D(kst-2,i,j,NK,NI,GUIDE)]*r;
    d[_IDX_S3D(ked-1,i,j,NK,NI,GUIDE)] +=  rhs[_IDX_S3D(ked  ,i,j,NK,NI,GUIDE)]*r;
    TIMING_stop("LSOR_TDMA_BC", 4.0);

    TIMING_start("LSOR_TDMA");
    /*
    tdma(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)],
           &a[kst+GUIDE-1],
           &b[kst+GUIDE-1],
           &c[kst+GUIDE-1],
           &w[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)]);

    tdma1(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)],
           &a[kst+GUIDE-1],
           &c[kst+GUIDE-1],
           &w[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)]);

    tdma2(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)],
           &a[kst+GUIDE-1],
           &c[kst+GUIDE-1],
           &w[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)]);
    */
    tdma3(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)],
           &a[kst+GUIDE-1],
           &c[kst+GUIDE-1],
           &w[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)]);
    TIMING_stop("LSOR_TDMA", f2);


    TIMING_start("LSOR_Relax_Body");
    #pragma vector always
    #pragma ivdep
    for (int k=kst-1; k<ked; k++) {
      size_t m = _IDX_S3D(k,i,j,NK,NI,GUIDE);
      pp = x[m];
      dp = ( d[m] - pp ) * omg;
      pn = pp + dp;
      x[m] = pn;
      res += dp * dp;
    }
    TIMING_stop("LSOR_Relax_Body", f3);

  }}

}


/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 */
void CZ::lsor_simd2(REAL_TYPE* d,
                   REAL_TYPE* x,
                   REAL_TYPE* w,
                   REAL_TYPE* a,
                   REAL_TYPE* c,
                   REAL_TYPE* rhs,
                   REAL_TYPE* msk,
                   double &res,
                   double &flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(w, ALIGN);
  __assume_aligned(a, ALIGN);
  __assume_aligned(c, ALIGN);
  __assume_aligned(msk, ALIGN);
  __assume_aligned(rhs, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int ist = innerFidx[I_minus];
  int ied = innerFidx[I_plus];
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  int nn  = ked - kst + 1;
  double flop_count;
  REAL_TYPE r = 1.0/6.0;

  res = 0.0;

  #pragma omp parallel for schedule(dynamic,1) \
              reduction(+:res) private(flop_count)
  for (int j=jst-1; j<jed; j++) {
  for (int i=ist-1; i<ied; i+=4) { // AVX512 > ストライド４ではマスクできない

    TIMING_start("LSOR_RHS");
    flop_count = 0.0;
    ms_rhs4(i, j, kst, ked, d, x, rhs, msk, flop_count);
    TIMING_stop("LSOR_RHS", flop_count);



    TIMING_start("LSOR_TDMA_BC");
    d[_IDX_S3D(kst-1,i  ,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,i  ,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(kst-2,i  ,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(kst-1,i  ,j,NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,i  ,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,i  ,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(ked  ,i  ,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(ked-1,i  ,j,NK,NI,GUIDE)];
    d[_IDX_S3D(kst-1,i+1,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,i+1,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(kst-2,i+1,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(kst-1,i+1,j,NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,i+1,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,i+1,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(ked  ,i+1,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(ked-1,i+1,j,NK,NI,GUIDE)];

    d[_IDX_S3D(kst-1,i+2,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,i+2,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(kst-2,i+2,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(kst-1,i+2,j,NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,i+2,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,i+2,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(ked  ,i+2,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(ked-1,i+2,j,NK,NI,GUIDE)];

    d[_IDX_S3D(kst-1,i+3,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,i+3,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(kst-2,i+3,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(kst-1,i+3,j,NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,i+3,j,NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,i+3,j,NK,NI,GUIDE)]
                                         + rhs[_IDX_S3D(ked  ,i+3,j,NK,NI,GUIDE)] * r )
                                         * msk[_IDX_S3D(ked-1,i+3,j,NK,NI,GUIDE)];
    TIMING_stop("LSOR_TDMA_BC", 24.0);


    TIMING_start("LSOR_TDMA");
    flop_count = 0.0;
    tdma4(nn,
         &d[_IDX_S3D(kst-1,i  ,j,NK,NI,GUIDE)],
         &d[_IDX_S3D(kst-1,i+1,j,NK,NI,GUIDE)],
         &d[_IDX_S3D(kst-1,i+2,j,NK,NI,GUIDE)],
         &d[_IDX_S3D(kst-1,i+3,j,NK,NI,GUIDE)],
         &a[kst+GUIDE-1],
         &c[kst+GUIDE-1],
         &w[_IDX_S3D(kst-1,i  ,j,NK,NI,GUIDE)],
         &w[_IDX_S3D(kst-1,i+1,j,NK,NI,GUIDE)],
         &w[_IDX_S3D(kst-1,i+2,j,NK,NI,GUIDE)],
         &w[_IDX_S3D(kst-1,i+3,j,NK,NI,GUIDE)],
         flop_count
       );
    TIMING_stop("LSOR_TDMA", flop_count);


    TIMING_start("LSOR_Relax");
    flop_count = 0.0;
    res += relax4(i, j, kst, ked, d, x, msk, flop_count);
    //res += relax4s(i, j, kst, ked, d, x, msk, flop_count);
    TIMING_stop("LSOR_Relax", flop_count);
  }}

}


/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 */
void CZ::lsor_simd3(REAL_TYPE* d,
                   REAL_TYPE* x,
                   REAL_TYPE* w,
                   REAL_TYPE* a,
                   REAL_TYPE* c,
                   REAL_TYPE* rhs,
                   REAL_TYPE* msk,
                   double &res,
                   double &flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(w, ALIGN);
  __assume_aligned(a, ALIGN);
  __assume_aligned(c, ALIGN);
  __assume_aligned(msk, ALIGN);
  __assume_aligned(rhs, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int ist = innerFidx[I_minus];
  int ied = innerFidx[I_plus];
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  double flop_count = 0.0;
  REAL_TYPE r = 1.0/6.0;

  // (i,j)点をサンプルするためのインデクス計算
  int QI = ied - ist + 1;
  int QJ = jed - jst + 1;
  int QK = ked - kst + 1;

  res = 0.0;

  #pragma omp parallel for schedule(dynamic,1) reduction(+:res) private(flop_count)
  for (int l=0; l<QI*QJ; l+=4)
  {
    int ia[4], ja[4];
    sIndex(ia[0], ja[0], l  , QI, ist, jst);
    sIndex(ia[1], ja[1], l+1, QI, ist, jst);
    sIndex(ia[2], ja[2], l+2, QI, ist, jst);
    sIndex(ia[3], ja[3], l+3, QI, ist, jst);

    TIMING_start("LSOR_RHS");
    flop_count = 0.0;
    ms_rhs4v(ia, ja, kst, ked, d, x, rhs, msk, flop_count);
    TIMING_stop("LSOR_RHS", flop_count);


    TIMING_start("LSOR_TDMA_BC");
    d[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(kst-2,ia[0],ja[0],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,ia[0],ja[0],NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,ia[0],ja[0],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(ked  ,ia[0],ja[0],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(ked-1,ia[0],ja[0],NK,NI,GUIDE)];
    d[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(kst-2,ia[1],ja[1],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,ia[1],ja[1],NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,ia[1],ja[1],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(ked  ,ia[1],ja[1],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(ked-1,ia[1],ja[1],NK,NI,GUIDE)];

    d[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(kst-2,ia[2],ja[2],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,ia[2],ja[2],NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,ia[2],ja[2],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(ked  ,ia[2],ja[2],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(ked-1,ia[2],ja[2],NK,NI,GUIDE)];

    d[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(kst-2,ia[3],ja[3],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)];
    d[_IDX_S3D(ked-1,ia[3],ja[3],NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,ia[3],ja[3],NK,NI,GUIDE)]
                                               + rhs[_IDX_S3D(ked  ,ia[3],ja[3],NK,NI,GUIDE)] * r )
                                               * msk[_IDX_S3D(ked-1,ia[3],ja[3],NK,NI,GUIDE)];
    TIMING_stop("LSOR_TDMA_BC", 24.0);


    TIMING_start("LSOR_TDMA");
    flop_count = 0.0;
    tdma5(QK,
         &d[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)],
         &d[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)],
         &d[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)],
         &d[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)],
         &a[kst+GUIDE-1],
         &c[kst+GUIDE-1],
         &w[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)],
         &w[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)],
         &w[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)],
         &w[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)],
         flop_count
       );
    TIMING_stop("LSOR_TDMA", flop_count);


    TIMING_start("LSOR_Relax");
    flop_count = 0.0;
    res += relax4c(ia, ja, kst, ked, d, x, msk, flop_count);
    //res += relax_256(ia, ja, kst, ked, d, x, msk, flop_count);
    //res += relax_256s(ia, ja, kst, ked, d, x, msk, flop_count);
    TIMING_stop("LSOR_Relax", flop_count);
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
void CZ::lsor_simd4(REAL_TYPE* d,
                    REAL_TYPE* x,
                    REAL_TYPE* w,
                    REAL_TYPE* a,
                    REAL_TYPE* e,
                    REAL_TYPE* rhs,
                    REAL_TYPE* msk,
                    double &res,
                    double &flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(w, ALIGN);
  __assume_aligned(a, ALIGN);
  __assume_aligned(e, ALIGN);
  __assume_aligned(msk, ALIGN);
  __assume_aligned(rhs, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int ist = innerFidx[I_minus];
  int ied = innerFidx[I_plus];
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  double flop_count = 0.0;
  REAL_TYPE r = 1.0/6.0;

  // (i,j)点をサンプルするためのインデクス計算
  int QI = ied - ist + 1;
  int QJ = jed - jst + 1;
  int QK = ked - kst + 1;

  res = 0.0;

  for (int col=0; col<2; col++)
  {
    #pragma omp parallel for schedule(dynamic,1) reduction(+:res) private(flop_count)
    for (int l=col; l<QI*QJ; l+=8)
    {
      int ia[4], ja[4];
      sIndex(ia[0], ja[0], l  ,  QI, ist, jst);
      sIndex(ia[1], ja[1], l+2,  QI, ist, jst);
      sIndex(ia[2], ja[2], l+4,  QI, ist, jst);
      sIndex(ia[3], ja[3], l+6,  QI, ist, jst);


      TIMING_start("LSOR_RHS");
      flop_count = 0.0;
      ms_rhs4v(ia, ja, kst, ked, d, x, rhs, msk, flop_count);
      TIMING_stop("LSOR_RHS", flop_count);


      TIMING_start("LSOR_TDMA_BC");
      for (int t=0; t<4; t++)
      {
        d[_IDX_S3D(kst-1,ia[t],ja[t],NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,ia[t],ja[t],NK,NI,GUIDE)]
                                                   + rhs[_IDX_S3D(kst-2,ia[t],ja[t],NK,NI,GUIDE)] * r )
                                                   * msk[_IDX_S3D(kst-1,ia[t],ja[t],NK,NI,GUIDE)];
        d[_IDX_S3D(ked-1,ia[t],ja[t],NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,ia[t],ja[t],NK,NI,GUIDE)]
                                                   + rhs[_IDX_S3D(ked  ,ia[t],ja[t],NK,NI,GUIDE)] * r )
                                                   * msk[_IDX_S3D(ked-1,ia[t],ja[t],NK,NI,GUIDE)];
      }
      TIMING_stop("LSOR_TDMA_BC", 24.0);


      TIMING_start("LSOR_TDMA");
      flop_count = 0.0;
      tdma6(ia, ja, a, e, w, d, flop_count);
      TIMING_stop("LSOR_TDMA", flop_count);


      TIMING_start("LSOR_Relax");
      flop_count = 0.0;
      res += relax4c(ia, ja, kst, ked, d, x, msk, flop_count);
      TIMING_stop("LSOR_Relax", flop_count);
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
 */
void CZ::lsor_simd5(REAL_TYPE* d,
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
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(w, ALIGN);
  __assume_aligned(a, ALIGN);
  __assume_aligned(e, ALIGN);
  __assume_aligned(msk, ALIGN);
  __assume_aligned(rhs, ALIGN);
  __assume_aligned(d2, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int ist = innerFidx[I_minus];
  int ied = innerFidx[I_plus];
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  double flop_count = 0.0;
  REAL_TYPE r = 1.0/6.0;

  // (i,j)点をサンプルするためのインデクス計算
  int QI = ied - ist + 1;
  int QJ = jed - jst + 1;


  res = 0.0;

  for (int col=0; col<2; col++)
  {
    #pragma omp parallel for schedule(dynamic,1) reduction(+:res) private(flop_count)
    for (int l=col; l<QI*QJ; l+=16)
    {
      int ia[8], ja[8];
      sIndex(ia[0], ja[0], l  ,  QI, ist, jst);
      sIndex(ia[1], ja[1], l+2,  QI, ist, jst);
      sIndex(ia[2], ja[2], l+4,  QI, ist, jst);
      sIndex(ia[3], ja[3], l+6,  QI, ist, jst);
      sIndex(ia[4], ja[4], l+8,  QI, ist, jst);
      sIndex(ia[5], ja[5], l+10, QI, ist, jst);
      sIndex(ia[6], ja[6], l+12, QI, ist, jst);
      sIndex(ia[7], ja[7], l+14, QI, ist, jst);


      TIMING_start("LSOR_RHS");
      flop_count = 0.0;
      ms_rhs8v(ia, ja, kst, ked, d, x, rhs, msk, flop_count);
      TIMING_stop("LSOR_RHS", flop_count);


      TIMING_start("LSOR_TDMA_BC");
      for (int t=0; t<8; t++)
      {
        d[_IDX_S3D(kst-1,ia[t],ja[t],NK,NI,GUIDE)] = ( d[_IDX_S3D(kst-1,ia[t],ja[t],NK,NI,GUIDE)]
                                                   + rhs[_IDX_S3D(kst-2,ia[t],ja[t],NK,NI,GUIDE)] * r )
                                                   * msk[_IDX_S3D(kst-1,ia[t],ja[t],NK,NI,GUIDE)];
        d[_IDX_S3D(ked-1,ia[t],ja[t],NK,NI,GUIDE)] = ( d[_IDX_S3D(ked-1,ia[t],ja[t],NK,NI,GUIDE)]
                                                   + rhs[_IDX_S3D(ked  ,ia[t],ja[t],NK,NI,GUIDE)] * r )
                                                   * msk[_IDX_S3D(ked-1,ia[t],ja[t],NK,NI,GUIDE)];
      }
      TIMING_stop("LSOR_TDMA_BC", 48.0);


      TIMING_start("LSOR_TDMA");
      flop_count = 0.0;
      tdma7(ia, ja, a, e, w, d, d2, flop_count);
/*
      tdma8(QK,
           &a[kst+GUIDE-1],
           &e[kst+GUIDE-1],
           &w[kst+GUIDE-1],
           &d[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[4],ja[4],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[5],ja[5],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[6],ja[6],NK,NI,GUIDE)],
           &d[_IDX_S3D(kst-1,ia[7],ja[7],NK,NI,GUIDE)],
           flop_count
         );
         */
      TIMING_stop("LSOR_TDMA", flop_count);

      REAL_TYPE* swp = NULL;
      swp = d;
      d = d2;
      d2 = swp;

      TIMING_start("LSOR_Relax");
      flop_count = 0.0;
      res += relax8c(ia, ja, kst, ked, d, x, msk, flop_count);
      TIMING_stop("LSOR_Relax", flop_count);
    }
  }
}
