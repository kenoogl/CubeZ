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
  int gc = GUIDE;
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
      size_t m = _IDX_S3D(k,i,j,NK,NI,gc);
      d[m] = ( x[_IDX_S3D(k,i-1,j  ,NK, NI, gc)]
            +  x[_IDX_S3D(k,i+1,j  ,NK, NI, gc)]
            +  x[_IDX_S3D(k,i  ,j-1,NK, NI, gc)]
            +  x[_IDX_S3D(k,i  ,j+1,NK, NI, gc)]
          ) * r + rhs[m];
    }
    TIMING_stop("LSOR_simd_f1", f1);

    d[_IDX_S3D(kst-1,i,j,NK,NI,gc)] +=  rhs[_IDX_S3D(kst-2,i,j,NK,NI,gc)]*r;
    d[_IDX_S3D(ked-1,i,j,NK,NI,gc)] +=  rhs[_IDX_S3D(ked,i,j,NK,NI,gc)]*r;

    TIMING_start("LSOR_simd_f2");

    tdma2(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,gc)],
           &a[kst+gc-1],
           &b[kst+gc-1],
           &c[kst+gc-1],
           &w[_IDX_S3D(kst-1,i,j,NK,NI,gc)]);

    TIMING_stop("LSOR_simd_f2", f2);

    TIMING_start("LSOR_simd_f3");
    #pragma vector always
    #pragma ivdep
    for (int k=kst-1; k<ked; k++) {
      size_t m = _IDX_S3D(k,i,j,NK,NI,gc);
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
  int gc = GUIDE;

  int bst = NumSW-gc-1;
  int bed = NumSW*(NumSB+1)-gc-1;
  //printf("st=%d ed =%d\n", bst, bed);

  d[0] = d[0]/b[0];
  w[0] = c[0]/b[0];

  // Forward:Peel
  for (int i=1; i<bst; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }

  // Forward:SIMD body
  for (int i=bst; i<bed; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }

  // Forward:Reminder
  for (int i=bed; i<nx; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }


  // Backward:Reminder
  for (int i=nx-2; i>=bed; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

  // Backward:SIMD body
  for (int i=bed-1; i>=bst; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

  // Backward:Peel
  for (int i=bst-1; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

}
