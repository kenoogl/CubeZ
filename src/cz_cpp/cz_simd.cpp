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
    /*
    tdma(nn,
           &d[_IDX_S3D(kst-1,i,j,NK,NI,gc)],
           a,
           b,
           c,
           &w[_IDX_S3D(kst-1,i,j,NK,NI,gc)]);
           */
    tdma_0_(&nn,
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
