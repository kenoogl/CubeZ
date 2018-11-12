#ifndef _FFV_LS_H_
#define _FFV_LS_H_
//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (C) 2016-2017  Research Institute for Information Technology(RIIT), Kyushu University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   LS.h
 * @brief  LS Class header
 * @author riit
 */

 #ifndef DISABLE_MPI
 #include <CB_SubDomain.h> // 先頭にmpi.h
 #include <CB_Comm.h>
 #endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cz_Define.h"
#include "DomainInfo.h"
#include "cz_Ffunc.h"


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace pm_lib;

class LinearSolver : public DomainInfo {

private:
  int size[4];
  int guide;
  int IterationMax;
  REAL_TYPE dh;
  REAL_TYPE coef_ac;
  double eps;

  PerfMonitor* PM;   ///< PerfMonitor class

  REAL_TYPE* pcg_p;  ///< work for BiCGstab
  REAL_TYPE* pcg_p_; ///< work for BiCGstab
  REAL_TYPE* pcg_r;  ///< work for BiCGstab
  REAL_TYPE* pcg_r0; ///< work for BiCGstab
  REAL_TYPE* pcg_q ; ///< work for BiCGstab
  REAL_TYPE* pcg_s;  ///< work for BiCGstab
  REAL_TYPE* pcg_s_; ///< work for BiCGstab
  REAL_TYPE* pcg_t_; ///< work for BiCGstab
  REAL_TYPE* wk;
  REAL_TYPE* exs;
  REAL_TYPE* gosa;

  REAL_TYPE cf[7];

  FILE* fp;

public:

  /** デフォルトコンストラクタ */
  LinearSolver() {
    size[0] = size[1] = size[2] = size[3] = 0;
    guide = 0;
    IterationMax = 0;
    dh = 0.0;
    coef_ac = 0.0;
    eps = 0.0;

    PM  = NULL;
    pcg_p  = NULL;
    pcg_p_ = NULL;
    pcg_r  = NULL;
    pcg_r0 = NULL;
    pcg_q  = NULL;
    pcg_s  = NULL;
    pcg_s_ = NULL;
    pcg_t_ = NULL;
    wk = NULL;
    fp = NULL;

    // dummy coef
    for (int i=0; i<6; i++) {
      cf[i] = 1.0;
    }
    cf[6] = 6.0;
  }

  /** コンストラクタ **/
  LinearSolver(const int sz[4],
               const int guide,
               const REAL_TYPE dh,
               const int IterationMax,
               const REAL_TYPE coef_ac,
               const double eps,
               PerfMonitor* PM,
               REAL_TYPE* pcg_p,
               REAL_TYPE* pcg_p_,
               REAL_TYPE* pcg_r,
               REAL_TYPE* pcg_r0,
               REAL_TYPE* pcg_q,
               REAL_TYPE* pcg_s,
               REAL_TYPE* pcg_s_,
               REAL_TYPE* pcg_t_,
               REAL_TYPE* wk,
               REAL_TYPE* exs,
               REAL_TYPE* gosa,
               FILE* fp
               ) {
    this->guide = guide;
    this->dh = dh;
    this->IterationMax = IterationMax;
    this->coef_ac= coef_ac;
    this->eps = eps;

    this->PM     = PM;
    this->pcg_p  = pcg_p;
    this->pcg_p_ = pcg_p_;
    this->pcg_r  = pcg_r;
    this->pcg_r0 = pcg_r0;
    this->pcg_q  = pcg_q;
    this->pcg_s  = pcg_s;
    this->pcg_s_ = pcg_s_;
    this->pcg_t_ = pcg_t_;
    this->wk     = wk;
    this->exs    = exs;
    this->gosa   = gosa;

    this->fp = fp;

    for (int i=0; i<4; i++)
    {
      this->size[i] = sz[i];
    }

    for (int i=0; i<6; i++) {
      cf[i] = 1.0;
    }
    cf[6] = 6.0;
  }

  /** デストラクタ */
  ~LinearSolver() {}



protected:


  /**
   * @brief Fdot for 1 array
   * @retval  内積値
   * @param [in]     x     vector1
   * @param [in,out] flop  浮動小数点演算数
   */
  double Fdot1(REAL_TYPE* x, double& flop);


  /**
   * @brief Fdot for 2 arrays
   * @retval  内積値
   * @param [in]     x    vector1
   * @param [in]     y    vector2
   * @param [in,out] flop 浮動小数点演算数
   */
  double Fdot2(REAL_TYPE* x, REAL_TYPE* y, double& flop);

  /**
   * @brief Preconditioner
   * @param [in,out] x         解ベクトル
   * @param [in]     b         RHS vector
   * @param [in,out] flop      浮動小数点演算数
   * @param [in]     isPrecond 前処理あり(true)
   */
  void Preconditioner(REAL_TYPE* x,
                      REAL_TYPE* b,
                      double& flop,
                      const bool isPrecond,
                      int KindOfPrecond);



public:

  /**
   * @brief Jacobi relaxation
   * @retval 反復数
   * @param [in,out] x      解ベクトル
   * @param [in]     b      RHS vector
   * @param [in]     itrMax 反復最大値
   * @param [in,out] flop   浮動小数点演算数
   * @param [in]     converge_check 収束判定あり(true)
   */
  int Jacobi(REAL_TYPE* x,
             REAL_TYPE* b,
             const int itrMax,
             double& flop,
             bool converge_check=true);

  /**
   * @brief SOR法
   * @retval 反復数
   * @param [in,out] x      解ベクトル
   * @param [in]     b      RHS vector
   * @param [in]     itrMax 反復最大値
   * @param [in,out] flop   浮動小数点演算数
   * @param [in]     converge_check 収束判定あり(true)
   */
  int PointSOR(REAL_TYPE* x,
               REAL_TYPE* b,
               const int itrMax,
               double& flop,
               bool converge_check=true);

/**
 * @brief SLOR法
 * @retval 反復数
 * @param [in,out] x      解ベクトル
 * @param [in]     b      RHS vector
 * @param [in]     itrMax 反復最大値
 * @param [in,out] flop   浮動小数点演算数
 * @param [in]     converge_check 収束判定あり(true)
 */
    int lineSOR(REAL_TYPE* x,
             REAL_TYPE* b,
             const int itrMax,
             double& flop,
             bool converge_check=true);

  /**
   * @brief 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     itrMax         反復最大値
   * @param [in,out] flop           浮動小数点演算数
   * @param [in]     converge_check 収束判定あり(true)
   */
  int SOR2_SMA(REAL_TYPE* x,
               REAL_TYPE* b,
               const int itrMax,
               double& flop,
               bool converge_check=true);


  /**
   * @brief 前処理つきBiCGstab
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     itrMax  反復最大値
   * @param [in,out] flop    浮動小数点演算数
   * @param [in]     isPrecond 前処理あり(true)
   */
  int PBiCGstab(REAL_TYPE* x,
                REAL_TYPE* b,
                const int itrMax,
                double& flop,
                const bool isPrecond,
                int KindOfPrecond);


  int residual_sp(REAL_TYPE* x,
                  REAL_TYPE* b,
                  double& flop);

};

#endif // _FFV_LS_H_
