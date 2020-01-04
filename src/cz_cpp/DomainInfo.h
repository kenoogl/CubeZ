#ifndef _CZ_DOMAIN_INFO_H_
#define _CZ_DOMAIN_INFO_H_

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
 * @file   DomainInfo.h
 * @brief  FlowBase DomainInfo class Header
 */

#include "cz_Define.h"
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

class DomainInfo {

public:
  int procGrp;            ///< プロセスグループ番号
  int myRank;             ///< 自ノードのランク番号
  int numProc;            ///< 全ランク数
  int numThreads;         ///< スレッド数

  int nID[NOFACE];        ///< 隣接ブロックのランク番号
  int head[3];            ///< 開始インデクス（グローバルインデクス, Fortran）
  int G_div[3];           ///< プロセス分割数

  REAL_TYPE pitch[3];     ///< 格子幅 (Non-dimensional)

  int size[3];            ///< 格子数 (Local, Non-dimensional)
  REAL_TYPE origin[3];    ///< 領域基点   (Local, Non-dimensional)
  REAL_TYPE region[3];    ///< 領域サイズ (Local, Non-dimensional)

  int G_size[3];          ///< 領域分割数 (Global, Non-dimensional)
  REAL_TYPE G_origin[3];  ///< 領域基点   (Global, Non-dimensional)
  REAL_TYPE G_region[3];  ///< 領域サイズ (Global, Non-dimensional)

  int innerFidx[6];       ///< 内部領域の開始終了インデクス(Fortran)

  std::string HostName;   ///< ホスト名



  /** コンストラクタ */
  DomainInfo() {
    procGrp = 0;
    myRank  = -1;
    numProc = 1;
    numThreads = 1;
    for (int i=0; i<NOFACE; i++) nID[i] = -1;

    for (int i=0; i<3; i++)
    {
      head[i]       = 0;
      size[i]       = 0;
      G_size[i]     = 0;
      G_div[i]      = 0;
      pitch[i]      = 0.0;
      origin[i]     = 0.0;
      region[i]     = 0.0;
      G_origin[i]   = 0.0;
      G_region[i]   = 0.0;
    }

    for (int i=0; i<6; i++) innerFidx[i]=0;
  }

  /** デストラクタ */
  virtual ~DomainInfo() {}


public:

  void setVar_Parallel(const int m_myRank,
                       const int m_numProc,
                       const int m_numThreads)
  {
    myRank = m_myRank;
    numProc= m_numProc;
    numThreads = m_numThreads;
  }

  void setVar_Domain(const int m_G_div[3],
                     const int m_head[3],
                     const int m_size[3],
                     const int m_G_size[3])
  {
    for (int i=0; i<3; i++) {
      G_div[i] = m_G_div[i];
      head[i]  = m_head[i];
      size[i]  = m_size[i];
      G_size[i]= m_G_size[i];
    }
  }

  void setPitch(const REAL_TYPE m_pch[3]) {
    pitch[0] = m_pch[0];
    pitch[1] = m_pch[1];
    pitch[2] = m_pch[2];
  }

  // #################################################################
  /**
   * @brief S3D配列のアロケート
   * @param [in] sz 配列サイズ
   * @ret pointer
   */
  REAL_TYPE* Alloc_Real_S3D(const int* sz)
  {
    if ( !sz ) return NULL;

    size_t dims[3], nx;

    dims[0] = (size_t)(sz[0] + 2*GUIDE);
    dims[1] = (size_t)(sz[1] + 2*GUIDE);
    dims[2] = (size_t)(sz[2] + 2*GUIDE);

    nx = dims[0] * dims[1] * dims[2];

    REAL_TYPE* var = new REAL_TYPE[nx];

  #pragma omp parallel for
    for (int i=0; i<nx; i++) var[i]=0.0;

    return var;
  }

};

#endif // _CZ_DOMAIN_INFO_H_
