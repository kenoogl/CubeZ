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
 * @note simd intrinsic
 */
void CZ::tdma6(const int nx,
               REAL_TYPE* a,
               REAL_TYPE* e,
               REAL_TYPE* w,
               REAL_TYPE* d0,
               REAL_TYPE* d1,
               REAL_TYPE* d2,
               REAL_TYPE* d3,
               double& flop)
{
  __assume_aligned(a,  ALIGN);
  __assume_aligned(e,  ALIGN);
  __assume_aligned(w,  ALIGN);
  __assume_aligned(d0, ALIGN);
  __assume_aligned(d1, ALIGN);
  __assume_aligned(d2, ALIGN);
  __assume_aligned(d3, ALIGN);


  int bst = SdW-GUIDE-1;
  int bed = SdW*(SdB+1)-GUIDE-1;
  const double f1 = 12.0*(double)bst;
  const double f2 = 12.0*(double)(nx-bst+1);
  const double f3 = 8.0*(double)(nx-1-bed);
  const double f4 = 8.0*(double)bed;

  REAL_TYPE ee, aa, ww;



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
