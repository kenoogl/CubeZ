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
 * @param [in]     a    係数 L_1
 * @param [in]     b    係数 D
 * @param [in]     c    係数 U_1
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 *       cz_dsolver tdma_1 と同等

void CZ::tdma_m(int nx,
                REAL_TYPE* d,
                const REAL_TYPE a,
                const REAL_TYPE b,
                const REAL_TYPE c,
                REAL_TYPE* w)
{
  int ist = innerFidx[I_minus];
  int ied = innerFidx[I_plus];
  int jst = innerFidx[J_minus];
  int jed = innerFidx[J_plus];
  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];

  REAL_TYPE e;
  size_t m;

  #pragma omp parallel for collapse(2) private(m)
  for (int k=kst-1; k<ked; k++) {
  for (int j=jst-1; j<jed; j++) {
    m = _IDX_S3D(0,j,k,NI,NJ,GUIDE);
    d[m] = d[m]/b;
    w[m] = c/b;
  }}


  for (int i=1; i<nx; i++)
  {
    e = 1.0 / (b - a * w[i-1]);
    w[i] = e * c;
    d[i] = (d[i] - a * d[i-1]) * e;
  }

  for (int i=nx-2; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

}*/


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    係数 L_1
 * @param [in]     b    係数 D
 * @param [in]     c    係数 U_1
 * @param [in]     w    work vector (U_1)
 * @note i方向に領域分割なしを想定
 *       cz_dsolver tdma_1 と同等
 */
void CZ::tdma_s(int nx,
                REAL_TYPE* d,
                const REAL_TYPE a,
                const REAL_TYPE b,
                const REAL_TYPE c,
                REAL_TYPE* w)
{
  REAL_TYPE e;

  d[0] = d[0]/b;
  w[0] = c/b;

  for (int i=1; i<nx; i++)
  {
    e = 1.0 / (b - a * w[i-1]);
    w[i] = e * c;
    d[i] = (d[i] - a * d[i-1]) * e;
  }

  for (int i=nx-2; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }

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
void CZ::tdma(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* b, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  TIMING_start("TDMA_F_body");
  d[0] = d[0]/b[0];
  w[0] = c[0]/b[0];

  for (int i=1; i<nx; i++)
  {
    e = 1.0 / (b[i] - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }

  for (int i=nx-2; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_F_body", 16.0*(double)(nx-1));
}

// tdma()のb[]=1.0でスケーリング、省略
void CZ::tdma1(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  TIMING_start("TDMA_F_body");
  w[0] = c[0];

  for (int i=1; i<nx; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }

  for (int i=nx-2; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_F_body", 16.0*(double)(nx-1));
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
void CZ::tdma2(int nx, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c, REAL_TYPE* w)
{
  REAL_TYPE e;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 14.0;
  const double f2 = 2.0;

  w[0] = c[0];


  // Forward:Peel
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=1; i<bst; i++)
  {
    e = 1.0 / (1.0- a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_F_peel", f1*(double)(bst));

  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int i=bst; i<bed; i++)
  {
    e = 1.0 / (1.0- a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_F_body", f1*(double)(bed-bst+1));

  // Forward:Reminder
  TIMING_start("TDMA_F_remainder");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=bed; i<nx; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_F_remainder", f1*(double)(nx-bed+1));


  // Backward:Reminder
  TIMING_start("TDMA_R_remainder");
  #pragma loop count (SdW-GUIDE-3)
  for (int i=nx-2; i>=bed; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_R_remainder", f2*(double)(nx-1-bed));

  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  for (int i=bed-1; i>=bst; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_R_body", f2*(double)(bed-bst));

  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=bst-1; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_R_peel", f2*(double)(bst));

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
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int i=1; i<bst; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_F_peel", f1*(double)(bst));

  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int i=bst; i<nx; i++)
  {
    e = 1.0 / (1.0 - a[i] * w[i-1]);
    w[i] = e * c[i];
    d[i] = (d[i] - a[i] * d[i-1]) * e;
  }
  TIMING_stop("TDMA_F_body", f1*(double)(nx-bst+1));


  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int i=nx-2; i>=bed; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_R_peel", f2*(double)(nx-1-bed));

  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  for (int i=bed-1; i>=0; i--)
  {
    d[i] = d[i] - w[i] * d[i+1];
  }
  TIMING_stop("TDMA_R_body", f2*(double)(bed));

}


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note
 */
void CZ::tdma4(const int nx,
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
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    e0 = 1.0 / (1.0 - a[k] * w0[k-1]);
    w0[k] = e0 * c[k];
    d0[k] = (d0[k] - a[k] * d0[k-1]) * e0;

    e1 = 1.0 / (1.0 - a[k] * w1[k-1]);
    w1[k] = e1 * c[k];
    d1[k] = (d1[k] - a[k] * d1[k-1]) * e1;

    e2 = 1.0 / (1.0 - a[k] * w2[k-1]);
    w2[k] = e2 * c[k];
    d2[k] = (d2[k] - a[k] * d2[k-1]) * e2;

    e3 = 1.0 / (1.0 - a[k] * w3[k-1]);
    w3[k] = e3 * c[k];
    d3[k] = (d3[k] - a[k] * d3[k-1]) * e3;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int k=bst; k<nx; k++)
  {
    e0 = 1.0 / (1.0 - a[k] * w0[k-1]);
    w0[k] = e0 * c[k];
    d0[k] = (d0[k] - a[k] * d0[k-1]) * e0;

    e1 = 1.0 / (1.0 - a[k] * w1[k-1]);
    w1[k] = e1 * c[k];
    d1[k] = (d1[k] - a[k] * d1[k-1]) * e1;

    e2 = 1.0 / (1.0 - a[k] * w2[k-1]);
    w2[k] = e2 * c[k];
    d2[k] = (d2[k] - a[k] * d2[k-1]) * e2;

    e3 = 1.0 / (1.0 - a[k] * w3[k-1]);
    w3[k] = e3 * c[k];
    d3[k] = (d3[k] - a[k] * d3[k-1]) * e3;
  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;


  // Backward:Reminder
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int k=nx-2; k>=bed; k--)
  {
    d0[k] = d0[k] - w0[k] * d0[k+1];
    d1[k] = d1[k] - w1[k] * d1[k+1];
    d2[k] = d2[k] - w2[k] * d2[k+1];
    d3[k] = d3[k] - w3[k] * d3[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  for (int k=bed-1; k>=0; k--)
  {
    d0[k] = d0[k] - w0[k] * d0[k+1];
    d1[k] = d1[k] - w1[k] * d1[k+1];
    d2[k] = d2[k] - w2[k] * d2[k+1];
    d3[k] = d3[k] - w3[k] * d3[k+1];
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;

}


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note
 */
void CZ::tdma5(const int nx,
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
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    e0 = 1.0 / (1.0 - a[k] * w0[k-1]);
    w0[k] = e0 * c[k];
    d0[k] = (d0[k] - a[k] * d0[k-1]) * e0;

    e1 = 1.0 / (1.0 - a[k] * w1[k-1]);
    w1[k] = e1 * c[k];
    d1[k] = (d1[k] - a[k] * d1[k-1]) * e1;

    e2 = 1.0 / (1.0 - a[k] * w2[k-1]);
    w2[k] = e2 * c[k];
    d2[k] = (d2[k] - a[k] * d2[k-1]) * e2;

    e3 = 1.0 / (1.0 - a[k] * w3[k-1]);
    w3[k] = e3 * c[k];
    d3[k] = (d3[k] - a[k] * d3[k-1]) * e3;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int k=bst; k<nx; k++)
  {
    e0 = 1.0 / (1.0 - a[k] * w0[k-1]);
    w0[k] = e0 * c[k];
    d0[k] = (d0[k] - a[k] * d0[k-1]) * e0;

    e1 = 1.0 / (1.0 - a[k] * w1[k-1]);
    w1[k] = e1 * c[k];
    d1[k] = (d1[k] - a[k] * d1[k-1]) * e1;

    e2 = 1.0 / (1.0 - a[k] * w2[k-1]);
    w2[k] = e2 * c[k];
    d2[k] = (d2[k] - a[k] * d2[k-1]) * e2;

    e3 = 1.0 / (1.0 - a[k] * w3[k-1]);
    w3[k] = e3 * c[k];
    d3[k] = (d3[k] - a[k] * d3[k-1]) * e3;
  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;


  // Backward:Reminder
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int k=nx-2; k>=bed; k--)
  {
    d0[k] = d0[k] - w0[k] * d0[k+1];
    d1[k] = d1[k] - w1[k] * d1[k+1];
    d2[k] = d2[k] - w2[k] * d2[k+1];
    d3[k] = d3[k] - w3[k] * d3[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  for (int k=bed-1; k>=0; k--)
  {
    d0[k] = d0[k] - w0[k] * d0[k+1];
    d1[k] = d1[k] - w1[k] * d1[k+1];
    d2[k] = d2[k] - w2[k] * d2[k+1];
    d3[k] = d3[k] - w3[k] * d3[k+1];
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;
}

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
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note LUの共通化
 */
void CZ::tdma6(const int* ia,
               const int* ja,
               REAL_TYPE* ap,
               REAL_TYPE* ep,
               REAL_TYPE* wp,
               REAL_TYPE* dp,
               double& flop)
{
  __assume_aligned(ap, ALIGN);
  __assume_aligned(ep, ALIGN);
  __assume_aligned(wp, ALIGN);
  __assume_aligned(dp, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 12.0*(double)bst;
  const double f2 = 12.0*(double)(nx-bst+1);
  const double f3 = 8.0*(double)(nx-1-bed);
  const double f4 = 8.0*(double)bed;

  REAL_TYPE ee, aa, ww;

  float *a, *e, *w;
  float *d0, *d1, *d2, *d3;

  a = &ap[kst+GUIDE-1];
  e = &ep[kst+GUIDE-1];
  w = &wp[kst+GUIDE-1];

  d0 = &dp[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
  d1 = &dp[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)],
  d2 = &dp[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)],
  d3 = &dp[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)],

  // Forward:Peel
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int k=bst; k<nx; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;


  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int k=nx-2; k>=bed; k--)
  {
    ww = w[k];
    d0[k] -= ww * d0[k+1];
    d1[k] -= ww * d1[k+1];
    d2[k] -= ww * d2[k+1];
    d3[k] -= ww * d3[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  for (int k=bed-1; k>=0; k--)
  {
    ww = w[k];
    d0[k] -= ww * d0[k+1];
    d1[k] -= ww * d1[k+1];
    d2[k] -= ww * d2[k+1];
    d3[k] -= ww * d3[k+1];
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;
}


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note tdma6() の8段
 */
void CZ::tdma6_8(const int* ia,
                 const int* ja,
                 REAL_TYPE* restrict ap,
                 REAL_TYPE* restrict ep,
                 REAL_TYPE* restrict wp,
                 REAL_TYPE* restrict dp,
                 double& flop)
{
  __assume_aligned(ap, ALIGN);
  __assume_aligned(ep, ALIGN);
  __assume_aligned(wp, ALIGN);
  __assume_aligned(dp, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 24.0*(double)bst;
  const double f2 = 24.0*(double)(nx-bst+1);
  const double f3 = 16.0*(double)(nx-1-bed);
  const double f4 = 16.0*(double)bed;

  REAL_TYPE ee, aa, ww;

  float *a, *e, *w;
  float *d0, *d1, *d2, *d3, *d4, *d5, *d6, *d7;

  a = &ap[kst+GUIDE-1];
  e = &ep[kst+GUIDE-1];
  w = &wp[kst+GUIDE-1];

  d0 = &dp[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
  d1 = &dp[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)];
  d2 = &dp[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)];
  d3 = &dp[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)];
  d4 = &dp[_IDX_S3D(kst-1,ia[4],ja[4],NK,NI,GUIDE)];
  d5 = &dp[_IDX_S3D(kst-1,ia[5],ja[5],NK,NI,GUIDE)];
  d6 = &dp[_IDX_S3D(kst-1,ia[6],ja[6],NK,NI,GUIDE)];
  d7 = &dp[_IDX_S3D(kst-1,ia[7],ja[7],NK,NI,GUIDE)];

  // Forward:Peel
  TIMING_start("TDMA_F_peel");
  #pragma vector always
  #pragma ivdep
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
    d4[k] = (d4[k] - aa * d4[k-1]) * ee;
    d5[k] = (d5[k] - aa * d5[k-1]) * ee;
    d6[k] = (d6[k] - aa * d6[k-1]) * ee;
    d7[k] = (d7[k] - aa * d7[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  #pragma vector always
  #pragma ivdep
  #pragma unroll(4)
  for (int k=bst; k<nx; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
    d4[k] = (d4[k] - aa * d4[k-1]) * ee;
    d5[k] = (d5[k] - aa * d5[k-1]) * ee;
    d6[k] = (d6[k] - aa * d6[k-1]) * ee;
    d7[k] = (d7[k] - aa * d7[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;


  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  #pragma vector always
  #pragma ivdep
  for (int k=nx-2; k>=bed; k--)
  {
    ww = w[k];
    d0[k] -= ww * d0[k+1];
    d1[k] -= ww * d1[k+1];
    d2[k] -= ww * d2[k+1];
    d3[k] -= ww * d3[k+1];
    d4[k] -= ww * d4[k+1];
    d5[k] -= ww * d5[k+1];
    d6[k] -= ww * d6[k+1];
    d7[k] -= ww * d7[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  #pragma vector always
  #pragma ivdep
  #pragma unroll(4)
  for (int k=bed-1; k>=0; k--)
  {
    ww = w[k];
    d0[k] -= ww * d0[k+1];
    d1[k] -= ww * d1[k+1];
    d2[k] -= ww * d2[k+1];
    d3[k] -= ww * d3[k+1];
    d4[k] -= ww * d4[k+1];
    d5[k] -= ww * d5[k+1];
    d6[k] -= ww * d6[k+1];
    d7[k] -= ww * d7[k+1];
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;
}


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note tdma6() の8段
 */
void CZ::tdma6_8_4(const int* ia,
                   const int* ja,
                   REAL_TYPE* restrict ap,
                   REAL_TYPE* restrict ep,
                   REAL_TYPE* restrict wp,
                   REAL_TYPE* restrict dp,
                   double& flop)
{
  __assume_aligned(ap, ALIGN);
  __assume_aligned(ep, ALIGN);
  __assume_aligned(wp, ALIGN);
  __assume_aligned(dp, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 24.0*(double)bst;
  const double f2 = 24.0*(double)(nx-bst+1);
  const double f3 = 16.0*(double)(nx-1-bed);
  const double f4 = 16.0*(double)bed;

  REAL_TYPE e0, a0, w0;
  REAL_TYPE e1, a1, w1;
  REAL_TYPE e2, a2, w2;
  REAL_TYPE e3, a3, w3;

  float *a, *e, *w;
  float *d0, *d1, *d2, *d3, *d4, *d5, *d6, *d7;

  a = &ap[kst+GUIDE-1];
  e = &ep[kst+GUIDE-1];
  w = &wp[kst+GUIDE-1];

  d0 = &dp[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
  d1 = &dp[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)];
  d2 = &dp[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)];
  d3 = &dp[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)];
  d4 = &dp[_IDX_S3D(kst-1,ia[4],ja[4],NK,NI,GUIDE)];
  d5 = &dp[_IDX_S3D(kst-1,ia[5],ja[5],NK,NI,GUIDE)];
  d6 = &dp[_IDX_S3D(kst-1,ia[6],ja[6],NK,NI,GUIDE)];
  d7 = &dp[_IDX_S3D(kst-1,ia[7],ja[7],NK,NI,GUIDE)];

  // Forward:Peel
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    a0 = a[k];
    e0 = e[k];
    d0[k] = (d0[k] - a0 * d0[k-1]) * e0;
    d1[k] = (d1[k] - a0 * d1[k-1]) * e0;
    d2[k] = (d2[k] - a0 * d2[k-1]) * e0;
    d3[k] = (d3[k] - a0 * d3[k-1]) * e0;
    d4[k] = (d4[k] - a0 * d4[k-1]) * e0;
    d5[k] = (d5[k] - a0 * d5[k-1]) * e0;
    d6[k] = (d6[k] - a0 * d6[k-1]) * e0;
    d7[k] = (d7[k] - a0 * d7[k-1]) * e0;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int k=bst; k<nx; k+=2)
  {
    a0 = a[k];
    e0 = e[k];
    d0[k] = (d0[k] - a0 * d0[k-1]) * e0;
    d1[k] = (d1[k] - a0 * d1[k-1]) * e0;
    d2[k] = (d2[k] - a0 * d2[k-1]) * e0;
    d3[k] = (d3[k] - a0 * d3[k-1]) * e0;
    d4[k] = (d4[k] - a0 * d4[k-1]) * e0;
    d5[k] = (d5[k] - a0 * d5[k-1]) * e0;
    d6[k] = (d6[k] - a0 * d6[k-1]) * e0;
    d7[k] = (d7[k] - a0 * d7[k-1]) * e0;
    a1 = a[k+1];
    e1 = e[k+1];
    d0[k+1] = (d0[k+1] - a1 * d0[k]) * e1;
    d1[k+1] = (d1[k+1] - a1 * d1[k]) * e1;
    d2[k+1] = (d2[k+1] - a1 * d2[k]) * e1;
    d3[k+1] = (d3[k+1] - a1 * d3[k]) * e1;
    d4[k+1] = (d4[k+1] - a1 * d4[k]) * e1;
    d5[k+1] = (d5[k+1] - a1 * d5[k]) * e1;
    d6[k+1] = (d6[k+1] - a1 * d6[k]) * e1;
    d7[k+1] = (d7[k+1] - a1 * d7[k]) * e1;
    /*
    a2 = a[k+2];
    e2 = e[k+2];
    d0[k+2] = (d0[k+2] - a2 * d0[k+1]) * e2;
    d1[k+2] = (d1[k+2] - a2 * d1[k+1]) * e2;
    d2[k+2] = (d2[k+2] - a2 * d2[k+1]) * e2;
    d3[k+2] = (d3[k+2] - a2 * d3[k+1]) * e2;
    d4[k+2] = (d4[k+2] - a2 * d4[k+1]) * e2;
    d5[k+2] = (d5[k+2] - a2 * d5[k+1]) * e2;
    d6[k+2] = (d6[k+2] - a2 * d6[k+1]) * e2;
    d7[k+2] = (d7[k+2] - a2 * d7[k+1]) * e2;
    a3 = a[k+3];
    e3 = e[k+3];
    d0[k+3] = (d0[k+3] - a3 * d0[k+2]) * e3;
    d1[k+3] = (d1[k+3] - a3 * d1[k+2]) * e3;
    d2[k+3] = (d2[k+3] - a3 * d2[k+2]) * e3;
    d3[k+3] = (d3[k+3] - a3 * d3[k+2]) * e3;
    d4[k+3] = (d4[k+3] - a3 * d4[k+2]) * e3;
    d5[k+3] = (d5[k+3] - a3 * d5[k+2]) * e3;
    d6[k+3] = (d6[k+3] - a3 * d6[k+2]) * e3;
    d7[k+3] = (d7[k+3] - a3 * d7[k+2]) * e3;
    */
  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;


  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int k=nx-2; k>=bed; k--)
  {
    w0 = w[k];
    d0[k] -= w0 * d0[k+1];
    d1[k] -= w0 * d1[k+1];
    d2[k] -= w0 * d2[k+1];
    d3[k] -= w0 * d3[k+1];
    d4[k] -= w0 * d4[k+1];
    d5[k] -= w0 * d5[k+1];
    d6[k] -= w0 * d6[k+1];
    d7[k] -= w0 * d7[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  #pragma unroll(8)
  for (int k=bed-1; k>=0; k-=2)
  {
    w0 = w[k];
    d0[k] -= w0 * d0[k+1];
    d1[k] -= w0 * d1[k+1];
    d2[k] -= w0 * d2[k+1];
    d3[k] -= w0 * d3[k+1];
    d4[k] -= w0 * d4[k+1];
    d5[k] -= w0 * d5[k+1];
    d6[k] -= w0 * d6[k+1];
    d7[k] -= w0 * d7[k+1];
    w1 = w[k+1];
    d0[k+1] -= w1 * d0[k+2];
    d1[k+1] -= w1 * d1[k+2];
    d2[k+1] -= w1 * d2[k+2];
    d3[k+1] -= w1 * d3[k+2];
    d4[k+1] -= w1 * d4[k+2];
    d5[k+1] -= w1 * d5[k+2];
    d6[k+1] -= w1 * d6[k+2];
    d7[k+1] -= w1 * d7[k+2];
    /*
    w2 = w[k+2];
    d0[k+2] -= w2 * d0[k+3];
    d1[k+2] -= w2 * d1[k+3];
    d2[k+2] -= w2 * d2[k+3];
    d3[k+2] -= w2 * d3[k+3];
    d4[k+2] -= w2 * d4[k+3];
    d5[k+2] -= w2 * d5[k+3];
    d6[k+2] -= w2 * d6[k+3];
    d7[k+2] -= w2 * d7[k+3];
    w3 = w[k+3];
    d0[k+3] -= w3 * d0[k+4];
    d1[k+3] -= w3 * d1[k+4];
    d2[k+3] -= w3 * d2[k+4];
    d3[k+3] -= w3 * d3[k+4];
    d4[k+3] -= w3 * d4[k+4];
    d5[k+3] -= w3 * d5[k+4];
    d6[k+3] -= w3 * d6[k+4];
    d7[k+3] -= w3 * d7[k+4];
    */
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;
}


/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note tdma6()と値の受け渡しが異なるもの
 */
void CZ::tdma7(const int* restrict ia,
               const int* restrict ja,
               REAL_TYPE* restrict ap,
               REAL_TYPE* restrict ep,
               REAL_TYPE* restrict wp,
               REAL_TYPE* restrict dp,
               REAL_TYPE* restrict dw,
               double& flop)
{
  __assume_aligned(ap, ALIGN);
  __assume_aligned(ep, ALIGN);
  __assume_aligned(wp, ALIGN);
  __assume_aligned(dp, ALIGN);
  __assume_aligned(dw, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;

  const double f1 = 24.0*(double)bst;
  const double f2 = 24.0*(double)(nx-bst+1);
  const double f3 = 16.0*(double)(nx-1-bed);
  const double f4 = 16.0*(double)bed;

  REAL_TYPE ee, aa, ww;
  float *a, *e, *w;
  float *d0, *d1, *d2, *d3, *d4, *d5, *d6, *d7;
  float *x0, *x1, *x2, *x3, *x4, *x5, *x6, *x7;

  a = &ap[kst+GUIDE-1];
  e = &ep[kst+GUIDE-1];
  w = &wp[kst+GUIDE-1];

  d0 = &dp[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
  d1 = &dp[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)];
  d2 = &dp[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)];
  d3 = &dp[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)];
  d4 = &dp[_IDX_S3D(kst-1,ia[4],ja[4],NK,NI,GUIDE)];
  d5 = &dp[_IDX_S3D(kst-1,ia[5],ja[5],NK,NI,GUIDE)];
  d6 = &dp[_IDX_S3D(kst-1,ia[6],ja[6],NK,NI,GUIDE)];
  d7 = &dp[_IDX_S3D(kst-1,ia[7],ja[7],NK,NI,GUIDE)];

  x0 = &dw[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
  x1 = &dw[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)];
  x2 = &dw[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)];
  x3 = &dw[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)];
  x4 = &dw[_IDX_S3D(kst-1,ia[4],ja[4],NK,NI,GUIDE)];
  x5 = &dw[_IDX_S3D(kst-1,ia[5],ja[5],NK,NI,GUIDE)];
  x6 = &dw[_IDX_S3D(kst-1,ia[6],ja[6],NK,NI,GUIDE)];
  x7 = &dw[_IDX_S3D(kst-1,ia[7],ja[7],NK,NI,GUIDE)];

  // Forward:Peel
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
    d4[k] = (d4[k] - aa * d4[k-1]) * ee;
    d5[k] = (d5[k] - aa * d5[k-1]) * ee;
    d6[k] = (d6[k] - aa * d6[k-1]) * ee;
    d7[k] = (d7[k] - aa * d7[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  for (int k=bst; k<nx; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
    d4[k] = (d4[k] - aa * d4[k-1]) * ee;
    d5[k] = (d5[k] - aa * d5[k-1]) * ee;
    d6[k] = (d6[k] - aa * d6[k-1]) * ee;
    d7[k] = (d7[k] - aa * d7[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;

  x0[nx-1] = d0[nx-1];
  x1[nx-1] = d1[nx-1];
  x2[nx-1] = d2[nx-1];
  x3[nx-1] = d3[nx-1];
  x4[nx-1] = d4[nx-1];
  x5[nx-1] = d5[nx-1];
  x6[nx-1] = d6[nx-1];
  x7[nx-1] = d7[nx-1];

  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int k=nx-2; k>=bed; k--)
  {
    ww = w[k];
    x0[k] = d0[k] - ww * x0[k+1];
    x1[k] = d1[k] - ww * x1[k+1];
    x2[k] = d2[k] - ww * x2[k+1];
    x3[k] = d3[k] - ww * x3[k+1];
    x4[k] = d4[k] - ww * x4[k+1];
    x5[k] = d5[k] - ww * x5[k+1];
    x6[k] = d6[k] - ww * x6[k+1];
    x7[k] = d7[k] - ww * x7[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // Backward:SIMD body
  TIMING_start("TDMA_R_body");
  for (int k=bed-1; k>=0; k--)
  {
    ww = w[k];
    x0[k] = d0[k] - ww * x0[k+1];
    x1[k] = d1[k] - ww * x1[k+1];
    x2[k] = d2[k] - ww * x2[k+1];
    x3[k] = d3[k] - ww * x3[k+1];
    x4[k] = d4[k] - ww * x4[k+1];
    x5[k] = d5[k] - ww * x5[k+1];
    x6[k] = d6[k] - ww * x6[k+1];
    x7[k] = d7[k] - ww * x7[k+1];
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;
}



/*
 * @brief Thomas Algorithm
 * @param [in      nx   配列長
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @param [in]     w    work vector (U_1)
 * @note Xeon Gold 6140 でうまく動いている、COre i7では結果不正
*/
void CZ::tdma8(const int* restrict ia,
               const int* restrict ja,
               REAL_TYPE* restrict ap,
               REAL_TYPE* restrict ep,
               REAL_TYPE* restrict wp,
               REAL_TYPE* restrict dp,
               double& flop)
{
  __assume_aligned(ap, ALIGN);
  __assume_aligned(ep, ALIGN);
  __assume_aligned(wp, ALIGN);
  __assume_aligned(dp, ALIGN);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  int kst = innerFidx[K_minus];
  int ked = innerFidx[K_plus];
  int nx = ked - kst + 1;

  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;

  const double f1 = 24.0*(double)bst;
  const double f2 = 24.0*(double)(nx-bst+1)/8;
  const double f3 = 16.0*(double)(nx-1-bed);
  const double f4 = 16.0*(double)bed/8;

  float aa, ee, ww;
  float* restrict a;
  float* restrict e;
  float* restrict w;

  float* restrict d0;
  float* restrict d1;
  float* restrict d2;
  float* restrict d3;
  float* restrict d4;
  float* restrict d5;
  float* restrict d6;
  float* restrict d7;

  a = &ap[kst+GUIDE-1];
  e = &ep[kst+GUIDE-1];
  w = &wp[kst+GUIDE-1];

  d0 = &dp[_IDX_S3D(kst-1,ia[0],ja[0],NK,NI,GUIDE)];
  d1 = &dp[_IDX_S3D(kst-1,ia[1],ja[1],NK,NI,GUIDE)];
  d2 = &dp[_IDX_S3D(kst-1,ia[2],ja[2],NK,NI,GUIDE)];
  d3 = &dp[_IDX_S3D(kst-1,ia[3],ja[3],NK,NI,GUIDE)];
  d4 = &dp[_IDX_S3D(kst-1,ia[4],ja[4],NK,NI,GUIDE)];
  d5 = &dp[_IDX_S3D(kst-1,ia[5],ja[5],NK,NI,GUIDE)];
  d6 = &dp[_IDX_S3D(kst-1,ia[6],ja[6],NK,NI,GUIDE)];
  d7 = &dp[_IDX_S3D(kst-1,ia[7],ja[7],NK,NI,GUIDE)];

  //check_align(&d0[bst], "d0");
  //check_align(&d1[bst], "d1");


  // Forward:Peel
  TIMING_start("TDMA_F_peel");
  #pragma loop count (SdW-GUIDE-2)
  for (int k=1; k<bst; k++)
  {
    aa = a[k];
    ee = e[k];
    d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    d1[k] = (d1[k] - aa * d1[k-1]) * ee;
    d2[k] = (d2[k] - aa * d2[k-1]) * ee;
    d3[k] = (d3[k] - aa * d3[k-1]) * ee;
    d4[k] = (d4[k] - aa * d4[k-1]) * ee;
    d5[k] = (d5[k] - aa * d4[k-1]) * ee;
    d6[k] = (d6[k] - aa * d6[k-1]) * ee;
    d7[k] = (d7[k] - aa * d7[k-1]) * ee;
  }
  TIMING_stop("TDMA_F_peel", f1);
  flop += f1;


  // 最初のd[k-1]を作っておく(transposed) unalign
  // d0[], d1[], d2[], d3[], d4[], d5[], d6[], d7[]
  // gather命令可能
  __m256 t0 = _mm256_set_ps(
            d7[bst-1], d6[bst-1], d5[bst-1], d4[bst-1],
            d3[bst-1], d2[bst-1], d1[bst-1], d0[bst-1]
         );

  // Forward:SIMD body
  TIMING_start("TDMA_F_body");
  #pragma ivdep
  for (int k=bst; k<nx; k+=8)
  {
    __m256 aaa = _mm256_load_ps(&a[k]);
    __m256 eee = _mm256_load_ps(&e[k]);

    __m256 t1[8], dd[8];

    t1[0] = _mm256_load_ps(&d0[k]);
    t1[1] = _mm256_load_ps(&d1[k]);
    t1[2] = _mm256_load_ps(&d2[k]);
    t1[3] = _mm256_load_ps(&d3[k]);
    t1[4] = _mm256_load_ps(&d4[k]);
    t1[5] = _mm256_load_ps(&d5[k]);
    t1[6] = _mm256_load_ps(&d6[k]);
    t1[7] = _mm256_load_ps(&d7[k]);


    _mm256_transpose_8x8_ps(dd, t1);

    // aa = a[k];
    // ee = e[k];
    // d0[k] = (d0[k] - aa * d0[k-1]) * ee;
    // fnmadd_ps(a,b,c) => -a*b+c

    __m256 ac = _mm256_broadcastss_ps( _mm256_extractf128_ps(aaa, 0) );
    dd[0] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, t0,    dd[0]),
            _mm256_broadcastss_ps( _mm256_extractf128_ps(eee, 0) )
            );

    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permute_ps(aaa, _MM_SHUFFLE(3,2,1,1)), 0) );
    dd[1] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[0], dd[1]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permute_ps(eee, _MM_SHUFFLE(3,2,1,1)), 0) )
            );

    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permute_ps(aaa, _MM_SHUFFLE(3,2,1,2)), 0) );
    dd[2] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[1], dd[2]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permute_ps(eee, _MM_SHUFFLE(3,2,1,2)), 0) )
            );

    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permute_ps(aaa, _MM_SHUFFLE(3,2,1,3)), 0) );
    dd[3] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[2], dd[3]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permute_ps(eee, _MM_SHUFFLE(3,2,1,3)), 0) )
            );

    eee = _mm256_set1_ps(e[k+4]);
    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(aaa, _mm256_set_epi32(7,6,5,4,3,2,1,4)),
           0) );
    dd[4] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[3], dd[4]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permutevar8x32_ps(eee, _mm256_set_epi32(7,6,5,4,3,2,1,4)),
                   0) )
          );

    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(aaa, _mm256_set_epi32(7,6,5,4,3,2,1,5)),
           0) );
    dd[5] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[4], dd[5]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permutevar8x32_ps(eee, _mm256_set_epi32(7,6,5,4,3,2,1,5)),
                   0) )
            );

    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(aaa, _mm256_set_epi32(7,6,5,4,3,2,1,6)),
           0) );
    dd[6] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[5], dd[6]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permutevar8x32_ps(eee, _mm256_set_epi32(7,6,5,4,3,2,1,6)),
                   0) )
            );

    ac = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(aaa, _mm256_set_epi32(7,6,5,4,3,2,1,7)),
           0) );
    dd[7] = _mm256_mul_ps(
            _mm256_fnmadd_ps(ac, dd[6], dd[7]),
            _mm256_broadcastss_ps(
                   _mm256_extractf128_ps(
                     _mm256_permutevar8x32_ps(eee, _mm256_set_epi32(7,6,5,4,3,2,1,7)),
                   0) )
            );

    // copy dd[7] > t0, 以降の処理の前に
    t0 = dd[7];

    // 出力のため転置
    _mm256_transpose_8x8_ps(t1, dd);

    _mm256_store_ps( d0+k , t1[0] );
    _mm256_store_ps( d1+k , t1[1] );
    _mm256_store_ps( d2+k , t1[2] );
    _mm256_store_ps( d3+k , t1[3] );
    _mm256_store_ps( d4+k , t1[4] );
    _mm256_store_ps( d5+k , t1[5] );
    _mm256_store_ps( d6+k , t1[6] );
    _mm256_store_ps( d7+k , t1[7] );

  }
  TIMING_stop("TDMA_F_body", f2);
  flop += f2;


  // Backward:Peel
  TIMING_start("TDMA_R_peel");
  #pragma loop count (SdW-GUIDE-3)
  for (int k=nx-2; k>=bed; k--)
  {
    ww = w[k];
    d0[k] = d0[k] - ww * d0[k+1];
    d1[k] = d1[k] - ww * d1[k+1];
    d2[k] = d2[k] - ww * d2[k+1];
    d3[k] = d3[k] - ww * d3[k+1];
    d4[k] = d4[k] - ww * d4[k+1];
    d5[k] = d5[k] - ww * d5[k+1];
    d6[k] = d6[k] - ww * d6[k+1];
    d7[k] = d7[k] - ww * d7[k+1];
  }
  TIMING_stop("TDMA_R_peel", f3);
  flop += f3;


  // x[k+1]の内容を作っておく
  t0 = _mm256_set_ps(
            d7[bed], d6[bed], d5[bed], d4[bed],
            d3[bed], d2[bed], d1[bed], d0[bed]
      );

  // Backward:SIMD body
  // 降順にアクセス 開始はbed-1ではなく、64バイト境界のbed-8から-GUIDE-1(先頭)まで
  // 計算順序はt1[7]から
  TIMING_start("TDMA_R_body");
  for (int k=bed-8; k>=-GUIDE-1; k-=8)
  {
    // ww = w[k];
    __m256 www = _mm256_load_ps(&w[k]);
    __m256 t1[8], dd[8], wc;

    t1[0] = _mm256_load_ps(&d0[k]);
    t1[1] = _mm256_load_ps(&d1[k]);
    t1[2] = _mm256_load_ps(&d2[k]);
    t1[3] = _mm256_load_ps(&d3[k]);
    t1[4] = _mm256_load_ps(&d4[k]);
    t1[5] = _mm256_load_ps(&d5[k]);
    t1[6] = _mm256_load_ps(&d6[k]);
    t1[7] = _mm256_load_ps(&d7[k]);

    _mm256_transpose_8x8_ps(dd, t1);


    // d0[k] = d0[k] - ww * d0[k+1];
    // fnmadd_ps(a,b,c) => -a*b+c

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(www, _mm256_set_epi32(7,6,5,4,3,2,1,7)),
           0) );
    dd[7] = _mm256_fnmadd_ps(wc, t0   , dd[7]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(www, _mm256_set_epi32(7,6,5,4,3,2,1,6)),
           0) );
    dd[6] = _mm256_fnmadd_ps(wc, dd[7], dd[6]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(www, _mm256_set_epi32(7,6,5,4,3,2,1,5)),
           0) );
    dd[5] = _mm256_fnmadd_ps(wc, dd[6], dd[5]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permutevar8x32_ps(www, _mm256_set_epi32(7,6,5,4,3,2,1,4)),
           0) );
    dd[4] = _mm256_fnmadd_ps(wc, dd[5], dd[4]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permute_ps(www, _MM_SHUFFLE(3,2,1,3)), 0) );
    dd[3] = _mm256_fnmadd_ps(wc, dd[4], dd[3]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permute_ps(www, _MM_SHUFFLE(3,2,1,2)), 0) );
    dd[2] = _mm256_fnmadd_ps(wc, dd[3], dd[2]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(
             _mm256_permute_ps(www, _MM_SHUFFLE(3,2,1,1)), 0) );
    dd[1] = _mm256_fnmadd_ps(wc, dd[2], dd[1]);

    wc = _mm256_broadcastss_ps(
           _mm256_extractf128_ps(www, 0) );
    dd[0] = _mm256_fnmadd_ps(wc, dd[1], dd[0]);

    t0 = dd[0];

    _mm256_transpose_8x8_ps(t1, dd);

    _mm256_store_ps( d0+k , t1[0] );
    _mm256_store_ps( d1+k , t1[1] );
    _mm256_store_ps( d2+k , t1[2] );
    _mm256_store_ps( d3+k , t1[3] );
    _mm256_store_ps( d4+k , t1[4] );
    _mm256_store_ps( d5+k , t1[5] );
    _mm256_store_ps( d6+k , t1[6] );
    _mm256_store_ps( d7+k , t1[7] );
  }
  TIMING_stop("TDMA_R_body", f4);
  flop += f4;
}
