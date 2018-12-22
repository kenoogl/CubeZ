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
