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

double CZ::relax4(const int i,
                 const int j,
                 const int kst,
                 const int ked,
                 REAL_TYPE* d,
                 REAL_TYPE* x,
                 REAL_TYPE* m,
                 double& flop)
{
  __assume_aligned(d, ALIGN_SIZE);
  __assume_aligned(x, ALIGN_SIZE);
  __assume_aligned(m, ALIGN_SIZE);

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

  TIMING_start("LSOR_Relax_Body");
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
  TIMING_stop("LSOR_Relax_Body", 24.0*(double)(ked-kst+2));

  return (double)res;
}


double CZ::relax4s(const int i,
                   const int j,
                   const int kst,
                   const int ked,
                   REAL_TYPE* d,
                   REAL_TYPE* x,
                   REAL_TYPE* m,
                   double& flop)
{
  __assume_aligned(d, ALIGN_SIZE);
  __assume_aligned(x, ALIGN_SIZE);
  __assume_aligned(m, ALIGN_SIZE);

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

  int km = kst+SdW-GUIDE-2;

  double f1 = 24.0*(double)(km-kst+1);
  double f2 = 24.0*(double)(ked-km+1);

  flop += (f1+f2);

  TIMING_start("LSOR_Relax_Peel");
  #pragma loop count (SdW-GUIDE-2)
  #pragma ivdep
  for (int k=kst-1; k<km; k++) {
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
  TIMING_stop("LSOR_Relax_Peel", f1);


  TIMING_start("LSOR_Relax_Body");
  #pragma vector always
  #pragma ivdep
  for (int k=km; k<ked; k++) {
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
  TIMING_stop("LSOR_Relax_Body", f2);

  return (double)res;
}

double CZ::relax4c(const int* ia,
                   const int* ja,
                   const int kst,
                   const int ked,
                   REAL_TYPE* restrict d,
                   REAL_TYPE* restrict x,
                   REAL_TYPE* restrict m,
                   double& flop)
{
  __assume_aligned(d, ALIGN_SIZE);
  __assume_aligned(x, ALIGN_SIZE);
  __assume_aligned(m, ALIGN_SIZE);

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

  TIMING_start("LSOR_Relax_Body");
  #pragma vector always
  #pragma ivdep
  for (int k=kst-1; k<ked; k++) {
    m0 = _IDX_S3D(k,ia[0],ja[0],NK,NI,GUIDE);
    pp0 = x[m0];
    dp0 = ( d[m0] - pp0 ) * omg * m[m0];
    pn0 = pp0 + dp0;
    x[m0] = pn0;

    m1 = _IDX_S3D(k,ia[1],ja[1],NK,NI,GUIDE);
    pp1 = x[m1];
    dp1 = ( d[m1] - pp1 ) * omg * m[m1];
    pn1 = pp1 + dp1;
    x[m1] = pn1;

    m2 = _IDX_S3D(k,ia[2],ja[2],NK,NI,GUIDE);
    pp2 = x[m2];
    dp2 = ( d[m2] - pp2 ) * omg * m[m2];
    pn2 = pp2 + dp2;
    x[m2] = pn2;

    m3 = _IDX_S3D(k,ia[3],ja[3],NK,NI,GUIDE);
    pp3 = x[m3];
    dp3 = ( d[m3] - pp3 ) * omg * m[m3];
    pn3 = pp3 + dp3;
    x[m3] = pn3;

    res += dp0 * dp0 + dp1 * dp1 + dp2 * dp2 + dp3 * dp3;
  }
  TIMING_stop("LSOR_Relax_Body", 24.0*(double)(ked-kst+2));

  return (double)res;
}

double CZ::relax8c(const int* ia,
                   const int* ja,
                   const int kst,
                   const int ked,
                   REAL_TYPE* d,
                   REAL_TYPE* x,
                   REAL_TYPE* m,
                   double& flop)
{
  __assume_aligned(d, ALIGN_SIZE);
  __assume_aligned(x, ALIGN_SIZE);
  __assume_aligned(m, ALIGN_SIZE);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  REAL_TYPE res=0.0;
  REAL_TYPE omg = ac1;
  REAL_TYPE pp0, dp0, pn0;
  REAL_TYPE pp1, dp1, pn1;
  REAL_TYPE pp2, dp2, pn2;
  REAL_TYPE pp3, dp3, pn3;
  REAL_TYPE pp4, dp4, pn4;
  REAL_TYPE pp5, dp5, pn5;
  REAL_TYPE pp6, dp6, pn6;
  REAL_TYPE pp7, dp7, pn7;
  size_t    m0, m1, m2, m3, m4, m5, m6, m7;

  flop += 48.0*(double)(ked-kst+2);

  TIMING_start("LSOR_Relax_Body");
  #pragma vector always
  #pragma ivdep
  #pragma unroll(2)
  for (int k=kst-1; k<ked; k++) {
    m0 = _IDX_S3D(k,ia[0],ja[0],NK,NI,GUIDE);
    pp0 = x[m0];
    dp0 = ( d[m0] - pp0 ) * omg * m[m0];
    pn0 = pp0 + dp0;
    x[m0] = pn0;

    m1 = _IDX_S3D(k,ia[1],ja[1],NK,NI,GUIDE);
    pp1 = x[m1];
    dp1 = ( d[m1] - pp1 ) * omg * m[m1];
    pn1 = pp1 + dp1;
    x[m1] = pn1;

    m2 = _IDX_S3D(k,ia[2],ja[2],NK,NI,GUIDE);
    pp2 = x[m2];
    dp2 = ( d[m2] - pp2 ) * omg * m[m2];
    pn2 = pp2 + dp2;
    x[m2] = pn2;

    m3 = _IDX_S3D(k,ia[3],ja[3],NK,NI,GUIDE);
    pp3 = x[m3];
    dp3 = ( d[m3] - pp3 ) * omg * m[m3];
    pn3 = pp3 + dp3;
    x[m3] = pn3;

    m4 = _IDX_S3D(k,ia[4],ja[4],NK,NI,GUIDE);
    pp4 = x[m4];
    dp4 = ( d[m4] - pp4 ) * omg * m[m4];
    pn4 = pp4 + dp4;
    x[m4] = pn4;

    m5 = _IDX_S3D(k,ia[5],ja[5],NK,NI,GUIDE);
    pp5 = x[m5];
    dp5 = ( d[m5] - pp5 ) * omg * m[m5];
    pn5 = pp5 + dp5;
    x[m5] = pn5;

    m6 = _IDX_S3D(k,ia[6],ja[6],NK,NI,GUIDE);
    pp6 = x[m6];
    dp6 = ( d[m6] - pp6 ) * omg * m[m6];
    pn6 = pp6 + dp6;
    x[m6] = pn6;

    m7 = _IDX_S3D(k,ia[7],ja[7],NK,NI,GUIDE);
    pp7 = x[m7];
    dp7 = ( d[m7] - pp7 ) * omg * m[m7];
    pn7 = pp7 + dp7;
    x[m7] = pn7;

    res += dp0 * dp0
         + dp1 * dp1
         + dp2 * dp2
         + dp3 * dp3
         + dp4 * dp4
         + dp5 * dp5
         + dp6 * dp6
         + dp7 * dp7;
  }
  TIMING_stop("LSOR_Relax_Body", 48.0*(double)(ked-kst+2));

  return (double)res;
}


// @note relax4c()をSIMD化
double CZ::relax_256(const int* ia,
                     const int* ja,
                     const int kst,
                     const int ked,
                     REAL_TYPE* d,
                     REAL_TYPE* x,
                     REAL_TYPE* msk,
                     double& flop)
{
  __assume_aligned(d, ALIGN_SIZE);
  __assume_aligned(x, ALIGN_SIZE);
  __assume_aligned(msk, ALIGN_SIZE);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  float coef = ac1;
  float tmp0 = 0.0;
  //float tmp1 = 0.0;
  //float tmp2 = 0.0;
  //float tmp3 = 0.0;

  double f2 = (20.0+64.0)*(double)(ked-kst+2)/8.0;

  flop += f2;


  __attribute__((aligned(32))) float t0[8] = {0};
  __attribute__((aligned(32))) float t1[8] = {0};
  __attribute__((aligned(32))) float t2[8] = {0};
  __attribute__((aligned(32))) float t3[8] = {0};

  //_mm256_store_ps(t0, omg);
  //printf("omg= %f %f %f %f %f %f %f %f \n", t0[0],t0[1],t0[2],t0[3],t0[4],t0[5],t0[6],t0[7]);

  TIMING_start("LSOR_Relax_Body");
  #pragma vector always
  #pragma ivdep
  for (int k=kst-1; k<ked; k+=8) {
    size_t m0 = _IDX_S3D(k,ia[0],ja[0],NK,NI,GUIDE);
    size_t m1 = _IDX_S3D(k,ia[1],ja[1],NK,NI,GUIDE);
    size_t m2 = _IDX_S3D(k,ia[2],ja[2],NK,NI,GUIDE);
    size_t m3 = _IDX_S3D(k,ia[3],ja[3],NK,NI,GUIDE);

    // pp0 = x[m0];
    __m256 pp0 = _mm256_loadu_ps(&x[m0]);
    __m256 pp1 = _mm256_loadu_ps(&x[m1]);
    __m256 pp2 = _mm256_loadu_ps(&x[m2]);
    __m256 pp3 = _mm256_loadu_ps(&x[m3]);

    // d[m0]
    __m256 d0 = _mm256_loadu_ps(&d[m0]);
    __m256 d1 = _mm256_loadu_ps(&d[m1]);
    __m256 d2 = _mm256_loadu_ps(&d[m2]);
    __m256 d3 = _mm256_loadu_ps(&d[m3]);

    // msk, tmp
    __m256 y0 = _mm256_loadu_ps(&msk[m0]);
    __m256 y1 = _mm256_loadu_ps(&msk[m1]);
    __m256 y2 = _mm256_loadu_ps(&msk[m2]);
    __m256 y3 = _mm256_loadu_ps(&msk[m3]);

    // dp0 = ( d[m0] - pp0 ) * omg * m[m0];
    __m256 omg = _mm256_set1_ps(coef);

    d0 = _mm256_mul_ps(
         _mm256_fmsub_ps( d0, omg, _mm256_mul_ps( pp0, omg ) ), y0
        );

    d1 = _mm256_mul_ps(
         _mm256_fmsub_ps( d1, omg, _mm256_mul_ps( pp1, omg ) ), y1
        );

    d2 = _mm256_mul_ps(
         _mm256_fmsub_ps( d2, omg, _mm256_mul_ps( pp2, omg ) ), y2
        );

    d3 = _mm256_mul_ps(
         _mm256_fmsub_ps( d3, omg, _mm256_mul_ps( pp3, omg ) ), y3
        );

    // pn0 = pp0 + dp0;
    y0 = _mm256_add_ps( pp0, d0 );
    y1 = _mm256_add_ps( pp1, d1 );
    y2 = _mm256_add_ps( pp2, d2 );
    y3 = _mm256_add_ps( pp3, d3 );

    // x[m0] = pn0;
    _mm256_storeu_ps( x+m0 , y0 );
    _mm256_storeu_ps( x+m1 , y1 );
    _mm256_storeu_ps( x+m2 , y2 );
    _mm256_storeu_ps( x+m3 , y3 );

    // res += dp0 * dp0 + dp1 * dp1 + dp2 * dp2 + dp3 * dp3;
    _mm256_store_ps(t0, _mm256_dp_ps(d0, d0, 0xFF));
    _mm256_store_ps(t1, _mm256_dp_ps(d1, d1, 0xFF));
    _mm256_store_ps(t2, _mm256_dp_ps(d2, d2, 0xFF));
    _mm256_store_ps(t3, _mm256_dp_ps(d3, d3, 0xFF));
    tmp0 += t0[0] + t0[4]
          + t1[0] + t1[4]
          + t2[0] + t2[4]
          + t3[0] + t3[4];

    /*
    _mm256_store_ps(t0, d0);
    tmp0 += t0[0] * t0[0]
          + t0[1] * t0[1]
          + t0[2] * t0[2]
          + t0[3] * t0[3]
          + t0[4] * t0[4]
          + t0[5] * t0[5]
          + t0[6] * t0[6]
          + t0[7] * t0[7];

    _mm256_store_ps(t1, d1);
    tmp1 += t1[0] * t1[0]
          + t1[1] * t1[1]
          + t1[2] * t1[2]
          + t1[3] * t1[3]
          + t1[4] * t1[4]
          + t1[5] * t1[5]
          + t1[6] * t1[6]
          + t1[7] * t1[7];

    _mm256_store_ps(t2, d2);
    tmp2 += t2[0] * t2[0]
          + t2[1] * t2[1]
          + t2[2] * t2[2]
          + t2[3] * t2[3]
          + t2[4] * t2[4]
          + t2[5] * t2[5]
          + t2[6] * t2[6]
          + t2[7] * t2[7];

    _mm256_store_ps(t3, d3);
    tmp3 += t3[0] * t3[0]
          + t3[1] * t3[1]
          + t3[2] * t3[2]
          + t3[3] * t3[3]
          + t3[4] * t3[4]
          + t3[5] * t3[5]
          + t3[6] * t3[6]
          + t3[7] * t3[7];
    */
  }
  TIMING_stop("LSOR_Relax_Body", f2);

  //return (double)( tmp0 + tmp1 + tmp2 + tmp3 );
  return (double)( tmp0 );
}


// @note relax_256()をPeel, _mm256_load_ps()
double CZ::relax_256s(const int* ia,
                      const int* ja,
                      const int kst,
                      const int ked,
                      REAL_TYPE* d,
                      REAL_TYPE* x,
                      REAL_TYPE* msk,
                      double& flop)
{
  __assume_aligned(d, ALIGN_SIZE);
  __assume_aligned(x, ALIGN_SIZE);
  __assume_aligned(msk, ALIGN_SIZE);

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];

  float coef = ac1;
  float tmp0 = 0.0;
  float tmp1 = 0.0;
  float tmp2 = 0.0;
  float tmp3 = 0.0;
  REAL_TYPE pp0, dp0, pn0;
  REAL_TYPE pp1, dp1, pn1;
  REAL_TYPE pp2, dp2, pn2;
  REAL_TYPE pp3, dp3, pn3;
  size_t    m0, m1, m2, m3;


  //_mm256_store_ps(t0, omg);
  //printf("omg= %f %f %f %f %f %f %f %f \n", t0[0],t0[1],t0[2],t0[3],t0[4],t0[5],t0[6],t0[7]);

  int km = kst+SdW-GUIDE-2;

  double f1 = 24.0*(double)(km-kst+1);
  double f2 = (16.0+64.0)*(double)(ked-km+1)/8.0;

  flop += (f1+f2);

  TIMING_start("LSOR_Relax_Peel");
  #pragma loop count (SdW-GUIDE-2)
  #pragma ivdep
  for (int k=kst-1; k<km; k++) {
    m0 = _IDX_S3D(k,ia[0],ja[0],NK,NI,GUIDE);
    pp0 = x[m0];
    dp0 = ( d[m0] - pp0 ) * coef * msk[m0];
    pn0 = pp0 + dp0;
    x[m0] = pn0;

    m1 = _IDX_S3D(k,ia[1],ja[1],NK,NI,GUIDE);
    pp1 = x[m1];
    dp1 = ( d[m1] - pp1 ) * coef * msk[m1];
    pn1 = pp1 + dp1;
    x[m1] = pn1;

    m2 = _IDX_S3D(k,ia[2],ja[2],NK,NI,GUIDE);
    pp2 = x[m2];
    dp2 = ( d[m2] - pp2 ) * coef * msk[m2];
    pn2 = pp2 + dp2;
    x[m2] = pn2;

    m3 = _IDX_S3D(k,ia[3],ja[3],NK,NI,GUIDE);
    pp3 = x[m3];
    dp3 = ( d[m3] - pp3 ) * coef * msk[m3];
    pn3 = pp3 + dp3;
    x[m3] = pn3;

    tmp0 += dp0 * dp0 + dp1 * dp1 + dp2 * dp2 + dp3 * dp3;
  }
  TIMING_stop("LSOR_Relax_Peel", f1);


  __attribute__((aligned(32))) float t0[8] = {0};
  __attribute__((aligned(32))) float t1[8] = {0};
  __attribute__((aligned(32))) float t2[8] = {0};
  __attribute__((aligned(32))) float t3[8] = {0};

  __m256 omg = _mm256_set1_ps(coef);

  TIMING_start("LSOR_Relax_Body");
  for (int k=km; k<ked; k+=8) {
    size_t m0 = _IDX_S3D(k,ia[0],ja[0],NK,NI,GUIDE);
    size_t m1 = _IDX_S3D(k,ia[1],ja[1],NK,NI,GUIDE);
    size_t m2 = _IDX_S3D(k,ia[2],ja[2],NK,NI,GUIDE);
    size_t m3 = _IDX_S3D(k,ia[3],ja[3],NK,NI,GUIDE);

    // pp0 = x[m0];
    __m256 pp0 = _mm256_load_ps(&x[m0]);
    __m256 pp1 = _mm256_load_ps(&x[m1]);
    __m256 pp2 = _mm256_load_ps(&x[m2]);
    __m256 pp3 = _mm256_load_ps(&x[m3]);

    // d[m0]
    __m256 d0 = _mm256_load_ps(&d[m0]);
    __m256 d1 = _mm256_load_ps(&d[m1]);
    __m256 d2 = _mm256_load_ps(&d[m2]);
    __m256 d3 = _mm256_load_ps(&d[m3]);

    // msk, tmp
    __m256 y0 = _mm256_load_ps(&msk[m0]);
    __m256 y1 = _mm256_load_ps(&msk[m1]);
    __m256 y2 = _mm256_load_ps(&msk[m2]);
    __m256 y3 = _mm256_load_ps(&msk[m3]);

    // dp0 = ( d[m0] - pp0 ) * omg * m[m0];
    __m256 dp0 = _mm256_mul_ps(
                    _mm256_fmsub_ps( d0, omg, _mm256_mul_ps( pp0, omg ) ), y0
                  );

    __m256 dp1 = _mm256_mul_ps(
                    _mm256_fmsub_ps( d1, omg, _mm256_mul_ps( pp1, omg ) ), y1
                  );

    __m256 dp2 = _mm256_mul_ps(
                    _mm256_fmsub_ps( d2, omg, _mm256_mul_ps( pp2, omg ) ), y2
                  );

    __m256 dp3 = _mm256_mul_ps(
                    _mm256_fmsub_ps( d3, omg, _mm256_mul_ps( pp3, omg ) ), y3
                  );

    // pn0 = pp0 + dp0;
    y0 = _mm256_add_ps( pp0, dp0 );
    y1 = _mm256_add_ps( pp1, dp1 );
    y2 = _mm256_add_ps( pp2, dp2 );
    y3 = _mm256_add_ps( pp3, dp3 );

    // x[m0] = pn0;
    _mm256_store_ps( x+m0 , y0 );
    _mm256_store_ps( x+m1 , y1 );
    _mm256_store_ps( x+m2 , y2 );
    _mm256_store_ps( x+m3 , y3 );

    // res += dp0 * dp0 + dp1 * dp1 + dp2 * dp2 + dp3 * dp3;
    _mm256_store_ps(t0, dp0);
    tmp0 += t0[0] * t0[0]
          + t0[1] * t0[1]
          + t0[2] * t0[2]
          + t0[3] * t0[3]
          + t0[4] * t0[4]
          + t0[5] * t0[5]
          + t0[6] * t0[6]
          + t0[7] * t0[7];

    _mm256_store_ps(t1, dp1);
    tmp1 += t1[0] * t1[0]
          + t1[1] * t1[1]
          + t1[2] * t1[2]
          + t1[3] * t1[3]
          + t1[4] * t1[4]
          + t1[5] * t1[5]
          + t1[6] * t1[6]
          + t1[7] * t1[7];

    _mm256_store_ps(t2, dp2);
    tmp2 += t2[0] * t2[0]
          + t2[1] * t2[1]
          + t2[2] * t2[2]
          + t2[3] * t2[3]
          + t2[4] * t2[4]
          + t2[5] * t2[5]
          + t2[6] * t2[6]
          + t2[7] * t2[7];

    _mm256_store_ps(t3, dp3);
    tmp3 += t3[0] * t3[0]
          + t3[1] * t3[1]
          + t3[2] * t3[2]
          + t3[3] * t3[3]
          + t3[4] * t3[4]
          + t3[5] * t3[5]
          + t3[6] * t3[6]
          + t3[7] * t3[7];
        }
  TIMING_stop("LSOR_Relax_Body", f2);

  return (double)( tmp0 + tmp1 + tmp2 + tmp3 );
}
