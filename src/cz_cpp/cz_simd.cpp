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

    TIMING_start("LSOR_simd_f1");
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
    TIMING_stop("LSOR_simd_f1", f1);

    d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)] +=  rhs[_IDX_S3D(kst-2,i,j,NK,NI,GUIDE)]*r;
    d[_IDX_S3D(ked-1,i,j,NK,NI,GUIDE)] +=  rhs[_IDX_S3D(ked  ,i,j,NK,NI,GUIDE)]*r;


    TIMING_start("LSOR_simd_f2_tdma");
    tdma2(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)],
           &a[kst+GUIDE-1],
           &b[kst+GUIDE-1],
           &c[kst+GUIDE-1],
           &w[_IDX_S3D(kst-1,i,j,NK,NI,GUIDE)]);
    TIMING_stop("LSOR_simd_f2_tdma", f2);


    TIMING_start("LSOR_simd_f3");
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
    TIMING_stop("LSOR_simd_f3", f3);

  }}

}

/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     b    D vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 *       cz_dsolver tdma_0 と同等
 */
void CZ::tdma2(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* b, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 14.0;
  const double f2 = 2.0;

  d[0] = d[0]/b[0];
  w[0] = c[0]/b[0];


  // Forward:Peel
  TIMING_start("TDMA_simd_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=1; i<bst; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_simd_F_peel", f1*(double)(bst));

  // Forward:SIMD body
  TIMING_start("TDMA_simd_F_body");
  for (int i=bst; i<bed; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_simd_F_body", f1*(double)(bed-bst+1));

  // Forward:Reminder
  TIMING_start("TDMA_simd_F_remainder");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=bed; i<nx; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_simd_F_remainder", f1*(double)(nx-bed+1));


  // Backward:Reminder
  TIMING_start("TDMA_simd_R_remainder");
  #pragma loop count (SdW-GUIDE-3)
  for (int i=nx-2; i>=bed; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_simd_R_remainder", f2*(double)(nx-1-bed));

  // Backward:SIMD body
  TIMING_start("TDMA_simd_R_body");
  for (int i=bed-1; i>=bst; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_simd_R_body", f2*(double)(bed-bst));

  // Backward:Peel
  TIMING_start("TDMA_simd_R_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=bst-1; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_simd_R_peel", f2*(double)(bst));

}

/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 *       cz_dsolver tdma_0 と同等
 */
void CZ::tdma3(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 14.0;
  const double f2 = 2.0;

  //d[0] = d[0];
  w[0] = c[0];


  // Forward:Peel
  TIMING_start("TDMA_simd_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=1; i<bst; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_simd_F_peel", f1*(double)(bst));

  // Forward:SIMD body
  TIMING_start("TDMA_simd_F_body");
  for (int i=bst; i<bed; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_simd_F_body", f1*(double)(bed-bst+1));

  // Forward:Reminder
  TIMING_start("TDMA_simd_F_remainder");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=bed; i<nx; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_simd_F_remainder", f1*(double)(nx-bed+1));


  // Backward:Reminder
  TIMING_start("TDMA_simd_R_remainder");
  #pragma loop count (SdW-GUIDE-3)
  for (int i=nx-2; i>=bed; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_simd_R_remainder", f2*(double)(nx-1-bed));

  // Backward:SIMD body
  TIMING_start("TDMA_simd_R_body");
  for (int i=bed-1; i>=bst; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_simd_R_body", f2*(double)(bed-bst));

  // Backward:Peel
  TIMING_start("TDMA_simd_R_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=bst-1; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_simd_R_peel", f2*(double)(bst));

}


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 */
void CZ::tdma4(const int nx,
               const int i,
               const int j,
               REAL_TYPE* d0,
               REAL_TYPE* d1,
               REAL_TYPE* d2,
               REAL_TYPE* d3,
               REAL_TYPE* a,
               REAL_TYPE* c,
               REAL_TYPE* w0,
               REAL_TYPE* w1,
               REAL_TYPE* w2,
               REAL_TYPE* w3,
               double& flop
             )
{
  __assume_aligned(d0, ALIGN);
  __assume_aligned(d1, ALIGN);
  __assume_aligned(d2, ALIGN);
  __assume_aligned(d3, ALIGN);
  __assume_aligned(w0, ALIGN);
  __assume_aligned(w1, ALIGN);
  __assume_aligned(w2, ALIGN);
  __assume_aligned(w3, ALIGN);
  __assume_aligned(a,  ALIGN);
  __assume_aligned(c,  ALIGN);

  REAL_TYPE e0, e1, e2, e3;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 56.0*(double)bst;
  const double f2 = 56.0*(double)(nx-bst+1);
  const double f3 = 8.0*(double)(nx-1-bed);
  const double f4 = 8.0*(double)bed;


  //d[0] = d[0];
  w0[0] = c[0];
  w1[0] = c[0];
  w2[0] = c[0];
  w3[0] = c[0];


  // Forward:Peel
  TIMING_start("TDMA_simd_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=1; i<bst; i++)
  {
    e0 = 1.0 / (1.0 - a[i] * w0[i-1]);
    w0[i] = e0 * c[i];
    d0[i] = (d0[i] - a[i] * d0[i-1]) * e0;

    e1 = 1.0 / (1.0 - a[i] * w1[i-1]);
    w1[i] = e1 * c[i];
    d1[i] = (d1[i] - a[i] * d1[i-1]) * e1;

    e2 = 1.0 / (1.0 - a[i] * w2[i-1]);
    w2[i] = e2 * c[i];
    d2[i] = (d2[i] - a[i] * d2[i-1]) * e2;

    e3 = 1.0 / (1.0 - a[i] * w3[i-1]);
    w3[i] = e3 * c[i];
    d3[i] = (d3[i] - a[i] * d3[i-1]) * e3;
  }
  TIMING_stop("TDMA_simd_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_simd_F_body");
  for (int i=bst; i<nx; i++)
  {
    e0 = 1.0 / (1.0 - a[i] * w0[i-1]);
    w0[i] = e0 * c[i];
    d0[i] = (d0[i] - a[i] * d0[i-1]) * e0;

    e1 = 1.0 / (1.0 - a[i] * w1[i-1]);
    w1[i] = e1 * c[i];
    d1[i] = (d1[i] - a[i] * d1[i-1]) * e1;

    e2 = 1.0 / (1.0 - a[i] * w2[i-1]);
    w2[i] = e2 * c[i];
    d2[i] = (d2[i] - a[i] * d2[i-1]) * e2;

    e3 = 1.0 / (1.0 - a[i] * w3[i-1]);
    w3[i] = e3 * c[i];
    d3[i] = (d3[i] - a[i] * d3[i-1]) * e3;
  }
  TIMING_stop("TDMA_simd_F_body", f2);
  flop += f2;


  // Backward:Reminder
  TIMING_start("TDMA_simd_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int i=nx-2; i>=bed; i--)
  {
    d0[i] = d0[i] - w0[i] * d0[i+1];
    d1[i] = d1[i] - w1[i] * d1[i+1];
    d2[i] = d2[i] - w2[i] * d2[i+1];
    d3[i] = d3[i] - w3[i] * d3[i+1];
  }
  TIMING_stop("TDMA_simd_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_simd_R_body");
  for (int i=bed-1; i>=0; i--)
  {
    d0[i] = d0[i] - w0[i] * d0[i+1];
    d1[i] = d1[i] - w1[i] * d1[i+1];
    d2[i] = d2[i] - w2[i] * d2[i+1];
    d3[i] = d3[i] - w3[i] * d3[i+1];
  }
  TIMING_stop("TDMA_simd_R_body", f4);
  flop += f4;

}


double CZ::relax(const int i,
                 const int j,
                 const int kst,
                 const int ked,
                 REAL_TYPE* d,
                 REAL_TYPE* x,
                 REAL_TYPE* m,
                 double& flop)
{
  __assume_aligned(d, ALIGN);
  __assume_aligned(x, ALIGN);
  __assume_aligned(m, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  REAL_TYPE res=0.0;
  REAL_TYPE omg = ac1;
  REAL_TYPE pp0, dp0, pn0;
  REAL_TYPE pp1, dp1, pn1;
  REAL_TYPE pp2, dp2, pn2;
  REAL_TYPE pp3, dp3, pn3;
  size_t    m0, m1, m2, m3;

  flop += 24.0*(double)(ked-kst+2);

  #pragma vector always
  #pragma ivdep
  for (int k=kst-1; k<ked; k++) {
    m0 = _IDX_S3D(k,i,j,NK,NI,GUIDE);
    pp0 = x[m0];
    dp0 = ( d[m0] - pp0 ) * omg * m[m0];
    pn0 = pp0 + dp0;
    x[m0] = pn0;

    m1 = _IDX_S3D(k,i+1,j,NK,NI,GUIDE);
    pp1 = x[m1];
    dp1 = ( d[m1] - pp1 ) * omg * m[m1];
    pn1 = pp1 + dp1;
    x[m1] = pn1;

    m2 = _IDX_S3D(k,i+2,j,NK,NI,GUIDE);
    pp2 = x[m2];
    dp2 = ( d[m2] - pp2 ) * omg * m[m2];
    pn2 = pp2 + dp2;
    x[m2] = pn2;

    m3 = _IDX_S3D(k,i+3,j,NK,NI,GUIDE);
    pp3 = x[m3];
    dp3 = ( d[m3] - pp3 ) * omg * m[m3];
    pn3 = pp3 + dp3;
    x[m3] = pn3;

    res += dp0 * dp0 + dp1 * dp1 + dp2 * dp2 + dp3 * dp3;
  }

  return (double)res;
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

  REAL_TYPE r = 1.0/6.0;
  size_t    m0, m1, m2, m3;

  int nn  = ked - kst + 1;
  double flop_count;

  double f1 = 24.0*(double)(ked-kst+2);

  res = 0.0;

  #pragma omp parallel for schedule(dynamic,1) \
              reduction(+:res) \
              private(m0, m1, m2, m3)
  for (int j=jst-1; j<jed; j++) {
  for (int i=ist-1; i<ied; i+=4) { // AVX512 > ストライド４ではマスクできない

    TIMING_start("LSOR_simd_f1");
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
    TIMING_stop("LSOR_simd_f1", f1);
    flop += f1;



    TIMING_start("LSOR_simd_bc");
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
    TIMING_stop("LSOR_simd_bc", 24.0);
    flop += 24.0;


    TIMING_start("LSOR_simd_f2");
    flop_count = 0.0;
    tdma4(nn, i, j,
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
    TIMING_stop("LSOR_simd_f2", flop_count);
    flop += flop_count;


    TIMING_start("LSOR_simd_f3");
    flop_count = 0.0;
    res += relax(i, j, kst, ked, d, x, msk, flop_count);
    TIMING_stop("LSOR_simd_f3", flop_count);
    flop += flop_count;

  }}

}


/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 * @note remainder loopをbodyに含める
 *                 マスクで調整
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

  REAL_TYPE r = 1.0/6.0;
  REAL_TYPE omg = ac1;
  REAL_TYPE pp0, dp0, pn0;
  REAL_TYPE pp1, dp1, pn1;
  REAL_TYPE pp2, dp2, pn2;
  REAL_TYPE pp3, dp3, pn3;
  REAL_TYPE e0, e1, e2, e3;
  size_t    m0, m1, m2, m3;
  int nn = ked - kst + 1;
  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;

  double f1 = 24.0*(double)(ked-kst+2);
  double f2 = 56.0*(double)(SdW-GUIDE-2);
  double f3 = 24.0*(double)(ked-kst+2);
  double f4 = 56.0*(double)(SdW*SdB+SdW-GUIDE-2);
  double f5 = 8.0*(double)(SdW-GUIDE-2);
  double f6 = 8.0*(double)(SdW*SdB+SdW-GUIDE-2);

  flop += (double)( (ied-ist+1)*(jed-jst+1)/4.0
             *(f1 + f2 + f3 + f4 + f5 + f6 + 24.0) );
  res = 0.0;

  REAL_TYPE c0 = c[kst+GUIDE-1];


  #pragma omp parallel for schedule(dynamic,1) \
              reduction(+:res) \
              private(pp0, dp0, pn0) \
              private(pp1, dp1, pn1) \
              private(pp2, dp2, pn2) \
              private(pp3, dp3, pn3) \
              private(m0, m1, m2, m3) \
              private(e0, e1, e2, e3)
  for (int j=jst-1; j<jed; j++) {
  for (int i=ist-1; i<ied; i+=4) { // AVX512 > ストライド４ではマスクできない

    TIMING_start("LSOR_simd_f1");
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
    TIMING_stop("LSOR_simd_f1", f1);



    TIMING_start("LSOR_simd_bc");
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
    TIMING_stop("LSOR_simd_bc", 24.0);


    w[_IDX_S3D(kst-1,i  ,j,NK,NI,GUIDE)] = c0;
    w[_IDX_S3D(kst-1,i+1,j,NK,NI,GUIDE)] = c0;
    w[_IDX_S3D(kst-1,i+2,j,NK,NI,GUIDE)] = c0;
    w[_IDX_S3D(kst-1,i+3,j,NK,NI,GUIDE)] = c0;


    // Forward:Peel
    TIMING_start("TDMA_simd_F_peel");
    #pragma loop count (SdW-GUIDE-2)
    for (int k=1; k<bst; k++)
    {
      m0 = _IDX_S3D(k+1,i  ,j,NK,NI,GUIDE);
      e0 = 1.0 / (1.0 - a[k] * w[m0-1]);
      w[m0] = e0 * c[k];
      d[m0] = (d[m0] - a[k] * d[m0-1]) * e0;

      m1 = _IDX_S3D(k+1,i+1,j,NK,NI,GUIDE);
      e1 = 1.0 / (1.0 - a[k] * w[m1-1]);
      w[m1] = e1 * c[k];
      d[m1] = (d[m1] - a[k] * d[m1-1]) * e1;

      m2 = _IDX_S3D(k+1,i+2,j,NK,NI,GUIDE);
      e2 = 1.0 / (1.0 - a[k] * w[m2-1]);
      w[m2] = e2 * c[k];
      d[m2] = (d[m2] - a[k] * d[m2-1]) * e2;

      m3 = _IDX_S3D(k+1,i+3,j,NK,NI,GUIDE);
      e3 = 1.0 / (1.0 - a[k] * w[m3-1]);
      w[m3] = e3 * c[k];
      d[m3] = (d[m3] - a[k] * d[m3-1]) * e3;
    }
    TIMING_stop("TDMA_simd_F_peel", f2);


    // Forward:SIMD body
    TIMING_start("TDMA_simd_F_body");
    for (int k=bst; k<nn; k++)
    {
      m0 = _IDX_S3D(k+1,i  ,j,NK,NI,GUIDE);
      e0 = 1.0 / (1.0 - a[k] * w[m0-1]);
      w[m0] = e0 * c[k];
      d[m0] = (d[m0] - a[k] * d[m0-1]) * e0;

      m1 = _IDX_S3D(k+1,i+1,j,NK,NI,GUIDE);
      e1 = 1.0 / (1.0 - a[k] * w[m1-1]);
      w[m1] = e1 * c[k];
      d[m1] = (d[m1] - a[k] * d[m1-1]) * e1;

      m2 = _IDX_S3D(k+1,i+2,j,NK,NI,GUIDE);
      e2 = 1.0 / (1.0 - a[k] * w[m2-1]);
      w[m2] = e2 * c[k];
      d[m2] = (d[m2] - a[k] * d[m2-1]) * e2;

      m3 = _IDX_S3D(k+1,i+3,j,NK,NI,GUIDE);
      e3 = 1.0 / (1.0 - a[k] * w[m3-1]);
      w[m3] = e3 * c[k];
      d[m3] = (d[m3] - a[k] * d[m3-1]) * e3;
    }
    TIMING_stop("TDMA_simd_F_body", f4);


    // Backward:Peel
    TIMING_start("TDMA_simd_R_peel");
    #pragma loop count (SdW-GUIDE-2)
    for (int k=nn-2; k>=bed; k--)
    {
      m0 = _IDX_S3D(k+1,i  ,j,NK,NI,GUIDE);
      d[m0] = d[m0] - w[m0] * d[m0+1];

      m1 = _IDX_S3D(k+1,i+1,j,NK,NI,GUIDE);
      d[m1] = d[m1] - w[m1] * d[m1+1];

      m2 = _IDX_S3D(k+1,i+2,j,NK,NI,GUIDE);
      d[m2] = d[m2] - w[m2] * d[m2+1];

      m3 = _IDX_S3D(k+1,i+3,j,NK,NI,GUIDE);
      d[m3] = d[m3] - w[m3] * d[m3+1];
    }
    TIMING_stop("TDMA_simd_R_peel", f5);


    // Backward:SIMD body
    TIMING_start("TDMA_simd_R_body");
    for (int k=bed-1; k>=0; k--)
    {
      m0 = _IDX_S3D(k+1,i  ,j,NK,NI,GUIDE);
      d[m0] = d[m0] - w[m0] * d[m0+1];

      m1 = _IDX_S3D(k+1,i+1,j,NK,NI,GUIDE);
      d[m1] = d[m1] - w[m1] * d[m1+1];

      m2 = _IDX_S3D(k+1,i+2,j,NK,NI,GUIDE);
      d[m2] = d[m2] - w[m2] * d[m2+1];

      m3 = _IDX_S3D(k+1,i+3,j,NK,NI,GUIDE);
      d[m3] = d[m3] - w[m3] * d[m3+1];
    }
    TIMING_stop("TDMA_simd_R_body", f6);


    TIMING_start("LSOR_simd_f3");
    #pragma vector always
    #pragma ivdep
    for (int k=kst-1; k<ked; k++) {
      m0 = _IDX_S3D(k,i,j,NK,NI,GUIDE);
      pp0 = x[m0];
      dp0 = ( d[m0] - pp0 ) * omg * msk[m0];
      pn0 = pp0 + dp0;
      x[m0] = pn0;

      m1 = _IDX_S3D(k,i+1,j,NK,NI,GUIDE);
      pp1 = x[m1];
      dp1 = ( d[m1] - pp1 ) * omg * msk[m1];
      pn1 = pp1 + dp1;
      x[m1] = pn1;

      m2 = _IDX_S3D(k,i+2,j,NK,NI,GUIDE);
      pp2 = x[m2];
      dp2 = ( d[m2] - pp2 ) * omg * msk[m2];
      pn2 = pp2 + dp2;
      x[m2] = pn2;

      m3 = _IDX_S3D(k,i+3,j,NK,NI,GUIDE);
      pp3 = x[m3];
      dp3 = ( d[m3] - pp3 ) * omg * msk[m3];
      pn3 = pp3 + dp3;
      x[m3] = pn3;

      res += dp0 * dp0 + dp1 * dp1 + dp2 * dp2 + dp3 * dp3;
    }
    TIMING_stop("LSOR_simd_f3", f3);

  }}

}
