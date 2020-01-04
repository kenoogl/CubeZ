#ifndef _CZ_DEFINE_H_
#define _CZ_DEFINE_H_

/*
###################################################################################
#
# CubeZ
#
# Copyright (C) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

/**
 * @file   cz_Define.h
 * @brief  CubeZ Definition Header
 */

#ifndef DISABLE_MPI
#include <CB_Define.h>
#else
#include "CB_Define_stub.h"
#endif


// precision
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#else
// /** 実数型の指定
//  * - デフォルトでは、REAL_TYPE=float
//   * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
//    *   REAL_TYPE=doubleになる
//     */
#define REAL_TYPE float
#endif


#define GUIDE 2  ///< ガイドセル数

#define DETAIL      2

#ifdef _OPENMP
#include <omp.h>
#endif

#include <float.h>

#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE_EPSILON DBL_MIN
#else
#define REAL_TYPE_EPSILON FLT_MIN
#endif


// PCRのアクセス
#define _IDX_(_I,_SS)  ( _I+_SS ) // SSは最大ストライド幅


// PMlibの登録ラベル個数
#define PM_NUM_MAX 200

// ラベルの最大文字数
#define TM_LABEL_MAX 24


enum LinearSolver
{
  LS_PSOR=1,
  LS_SOR2SMA,
  LS_BICGSTAB,
  LS_JACOBI,
  LS_PCR, // 5
  LS_PCR_EDA,
  LS_PCR_ESA,
  LS_PCR_RB,
  LS_PCR_RB_ESA,
  LS_PCR_J_ESA,
  LS_PSOR_MAF,
  LS_SOR2SMA_MAF, // 11
  LS_BICGSTAB_MAF,
  LS_JACOBI_MAF,
  LS_PCR_MAF,
  LS_PCR_EDA_MAF,
  LS_PCR_ESA_MAF, // 16
  LS_PCR_RB_MAF,
  LS_PCR_RB_ESA_MAF
};


#endif // _CZ_DEFINE_H_
