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

// @file cm_blas.cpp

/*
 * @brief LSOR Multi-System
 * @param [in]     d    RHS vector
 * @param [in,out] x    solution vector
 * @param [in]     w    work vector (U_1)
 * @param [in]     rhs
 * @param [out]    res  residual
 */
void CZ::lsor_ms(REAL_TYPE* d, REAL_TYPE* x, REAL_TYPE* w, REAL_TYPE* rhs, double &res, double &flop)
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

  int nn = ied - ist + 1;

  REAL_TYPE dd = 6.0;
  REAL_TYPE r = 1.0/6.0;
  REAL_TYPE a = -1.0/6.0;
  REAL_TYPE b = 1.0;
  REAL_TYPE c = -1.0/6.0;
  REAL_TYPE omg = ac1;
  REAL_TYPE pp, dp, pn;

  flop += (4.0 + 16.0 + 5.0)
        * (double)( (ied-ist+1)*(jed-jst+1)*(ked-kst+1) );
  res = 0.0;

  #pragma omp parallel for private(pp,dp,pn) reduction(+:res) collapse(2)
  for (int k=kst-1; k<ked; k++) {
  for (int j=jst-1; j<jed; j++) {

    #pragma vector always
    #pragma ivdep
    for (int i=ist-1; i<ied; i++) {
      d[_IDX_S3D(i,j,k,NI, NJ, gc)] = ( x[_IDX_S3D(i,j-1,k  ,NI, NJ, gc)]
                                     +  x[_IDX_S3D(i,j+1,k  ,NI, NJ, gc)]
                                     +  x[_IDX_S3D(i,j  ,k-1,NI, NJ, gc)]
                                     +  x[_IDX_S3D(i,j  ,k+1,NI, NJ, gc)]
                                   ) * r + rhs[_IDX_S3D(i,j,k,NI, NJ, gc)];
    }

    d[_IDX_S3D(ist,j,k,NI,NJ,gc)] +=  rhs[_IDX_S3D(ist-1,j,k,NI,NJ,gc)]*r;
    d[_IDX_S3D(ied,j,k,NI,NJ,gc)] +=  rhs[_IDX_S3D(ied+1,j,k,NI,NJ,gc)]*r;

    tdma_s(nn,
           &d[_IDX_S3D(ist-1,j,k,NI,NJ,gc)],
           a, b, c,
           &w[_IDX_S3D(ist-1,j,k,NI,NJ,gc)]);


    //tdma_m(nn, d, a, b, c, w);

    #pragma vector always
    #pragma ivdep
    for (int i=ist-1; i<ied; i++) {
      pp = x[_IDX_S3D(i,j,k,NI,NJ,gc)];
      dp = ( d[_IDX_S3D(i,j,k,NI,NJ,gc)] - pp ) * omg;
      pn = pp + dp;
      x[_IDX_S3D(i,j,k,NI,NJ,gc)] = pn;
      res += dp * dp;
    }

  }}

}


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
    m = _IDX_S3D(0,j,k,NI,NJ,gc);
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

}


/*
 * @brief Parallel Cyclic Reduction
 * @param [in      nx   配列長
 * @param [in      pn   nxを超える最小の2べき数の指数
 * @param [in,out] d    RHS/解ベクトル X[nx]
 * @param [in]     a    L_1 vector
 * @param [in]     c    U_1 vector
 * @note i方向に領域分割なしを想定
 *       tdma()とは異なり、解くべき範囲（内点）は[1,nx]

  nx = 12        // 解くべき方程式の次元数
  2^4=16 > nx    // nxを超える最小の2べき数
  pn=4           // その指数
  ss=2^{pn-1}=8  // PCRで参照するストライド s の最大値
  1-ss=-7        // 配列参照の下限
  nx+ss=20       // 配列参照の上限
  nx+2*ss        // 配列長

                          <-------------------------------->
  |--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--|
 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20

  #define _IDX_(_I,_SS)  ( (_I) - (_SS) ) で配列にアクセス

 */
void CZ::pcr(const int nx, const int pn, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c,
                                         REAL_TYPE* d1, REAL_TYPE* a1, REAL_TYPE* c1)
{
  REAL_TYPE r, ap, cp;

  const int ss = 0x1 << (pn-1);

  int s=0;
  for (int p=1; p<=pn; p++)
  //for (int p=1; p<pn; p++) // 一つ少なくすると2x2  or 3x3 の直接反転
  {
    s = 0x1 << (p-1); // s=2^{p-1}

    if (p==1)
    {
      int ip=0;
      pcr_kernel_2(nx, s, ss, ip, d, a, c);

      ip=1;
      pcr_kernel_2(nx, s, ss, ip, d1, a1, c1);

      pcr_merge(nx, ss, d, a, c, d1, a1, c1);
    }
    else
    {
      pcr_kernel_1(nx, s, ss, d, a, c);
    }

//pcr_kernel_1(nx, s, ss, d, a, c);
  }



}

void CZ::pcr2(const int nx, const int pn, REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c)
{

  const int ss = 0x1 << (pn-1);

  int s=0;
  for (int p=1; p<=pn; p++)
  //for (int p=1; p<pn; p++) // 一つ少なくすると2x2  or 3x3 の直接反転
  {
    s = 0x1 << (p-1); // s=2^{p-1}
/*
    if (p==1)
    {
      int ip=0;
      pcr_kernel_2(nx, s, ss, ip, d, a, c);

      ip=1;
      pcr_kernel_2(nx, s, ss, ip, d1, a1, c1);

      pcr_merge(nx, ss, d, a, c, d1, a1, c1);
    }
    else
    {
      pcr_kernel_1(nx, s, ss, d, a, c);
    }
*/
    pcr_kernel_3(nx, s, ss, d, a, c);
  }

}

void CZ::pcr_merge(const int nx, const int ss,
                   REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c,
                   REAL_TYPE* d1, REAL_TYPE* a1, REAL_TYPE* c1)
{
  for (int i=2; i<=nx; i+=2)
  {
    a[_IDX_(i,ss)] = a1[_IDX_(i,ss)];
    c[_IDX_(i,ss)] = c1[_IDX_(i,ss)];
    d[_IDX_(i,ss)] = d1[_IDX_(i,ss)];
  }

  printA(nx, ss, a, "A");
  printA(nx, ss, c, "C");
  printf("\n");

}

void CZ::printA(int nx, int ss, REAL_TYPE* a, char* s)
{
  printf("%s : ", s);
  for (int i=0; i<=nx+1; i++)
  {
    printf("%5.3f ", a[_IDX_(i,ss)] );
  }
  printf("\n");
}

void CZ::printB(int nx, REAL_TYPE* a, char* s)
{
  printf("%s : ", s);
  for (int i=0; i<=nx+1; i++)
  {
    printf("%7.1e ", a[i] );
  }
  printf("\n");
}


void CZ::pcr_kernel_1(const int nx, const int s, const int ss,
                    REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c)
{
  REAL_TYPE r, ap, cp;

  for (int i=1; i<=nx; i++)
  {
    ap = a[_IDX_(i,ss)];
    cp = c[_IDX_(i,ss)];
    r = 1.0 / ( 1.0 - ap * c[_IDX_(i-s,ss)] - cp * a[_IDX_(i+s,ss)] );
    a[_IDX_(i,ss)] = - r * ap * a[_IDX_(i-s,ss)];
    c[_IDX_(i,ss)] = - r * cp * c[_IDX_(i+s,ss)];
    d[_IDX_(i,ss)] =   r * ( d[_IDX_(i,ss)] - ap * d[_IDX_(i-s,ss)] - cp * d[_IDX_(i+s,ss)] );
  }
  printA(nx, ss, a, "A");
  printA(nx, ss, c, "C");
  printf("\n");
}

void CZ::pcr_kernel_2(const int nx, const int s, const int ss, const int ip,
                    REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c)
{
  REAL_TYPE r, ap, cp;

  for (int i=1+ip; i<=nx; i+=2)
  {
    ap = a[_IDX_(i,ss)];
    cp = c[_IDX_(i,ss)];
    r = 1.0 / ( 1.0 - ap * c[_IDX_(i-s,ss)] - cp * a[_IDX_(i+s,ss)] );
    a[_IDX_(i,ss)] = - r * ap * a[_IDX_(i-s,ss)];
    c[_IDX_(i,ss)] = - r * cp * c[_IDX_(i+s,ss)];
    d[_IDX_(i,ss)] =   r * ( d[_IDX_(i,ss)] - ap * d[_IDX_(i-s,ss)] - cp * d[_IDX_(i+s,ss)] );
  }
}


void CZ::pcr_kernel_3(const int nx, const int s, const int ss,
                    REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c)
{
  REAL_TYPE r, ap, cp;
  int iL, iR;

  for (int i=1; i<=nx; i++)
  {
    iL = i-s;
    iL = std::max(iL,0);
    iR = i+s;
		iR = std::min(iR,nx+1);

    ap = a[i];
    cp = c[i];
    r = 1.0 / ( 1.0 - ap * c[iL] - cp * a[iR] );
    a[i] = - r * ap * a[iL];
    c[i] = - r * cp * c[iR];
    d[i] =   r * ( d[i] - ap * d[iL] - cp * d[iR] );
  }
  printB(nx, a, "A");
  printB(nx, c, "C");
  printf("\n");
}
