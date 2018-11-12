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

#ifndef _CZ_H_
#define _CZ_H_

#ifndef DISABLE_MPI
#include <CB_SubDomain.h> // 先頭にmpi.h
#include <CB_Comm.h>
#endif

#include <stdio.h>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "cz_Define.h"
#include "cz_Ffunc.h"
#include "czVersion.h"
#include "DomainInfo.h"


// FX10 profiler
#if defined __K_FPCOLL
#include <fjcoll.h>
#elif defined __FX_FAPP
#include <fj_tool/fjcoll.h>
#endif

#include <PerfMonitor.h>


using namespace std;


class CZ : public DomainInfo {

private:

  double comm_size;
  int debug_mode;          ///< if 1, Debug mode
  int ItrMax;              ///< 最大反復回数
  int ls_type;             ///< 線形ソルバ種類
  double eps;              ///< convergence criteria
  REAL_TYPE ac1;           ///< acceleration coef. for SOR
  REAL_TYPE ac2;           ///< acceleration coef. for jacobi relaxation
  double res_normal;       ///< 全計算点数
  REAL_TYPE cf[7];         ///< 係数

  // PMlib
  pm_lib::PerfMonitor PM;  ///< 性能モニタクラス
  int order_of_PM_key;     ///< PMlib用の登録番号カウンタ < PM_NUM_MAX
  string Parallel_str;     ///< 並列モードの文字列

#ifndef DISABLE_MPI
  SubDomain D;               ///< 領域分割情報保持クラス
  BrickComm CM;              ///< 通信クラス
  MPI_Request req[NOFACE*2]; ///< Communication identifier for nonblocking
#endif

  FILE* fph;


public:
  REAL_TYPE* WRK;    ///< ワーク配列
  REAL_TYPE* P;      ///< 圧力
  REAL_TYPE* RHS;    ///< Poissonのソース項
  REAL_TYPE* MSK;    ///< マスク配列

  REAL_TYPE* EXS;    ///< 厳密解
  REAL_TYPE* ERR;    ///< 誤差

  REAL_TYPE* pcg_p;  ///< work for BiCGstab
  REAL_TYPE* pcg_p_; ///< work for BiCGstab
  REAL_TYPE* pcg_r;  ///< work for BiCGstab
  REAL_TYPE* pcg_r0; ///< work for BiCGstab
  REAL_TYPE* pcg_q ; ///< work for BiCGstab
  REAL_TYPE* pcg_s;  ///< work for BiCGstab
  REAL_TYPE* pcg_s_; ///< work for BiCGstab
  REAL_TYPE* pcg_t ; ///< work for BiCGstab
  REAL_TYPE* pcg_t_; ///< work for BiCGstab


public:
  // コンストラクタ
  CZ()
  {
    order_of_PM_key = 0;
    comm_size = 0.0;
    debug_mode = 0;
    ItrMax = 0;
    ls_type = 0;
    eps = 1.0e-5;
    ac1 = 1.1;
    ac2 = 0.8;
    res_normal = 0.0;

    for (int i=0; i<6; i++) {
      cf[i] = 1.0;
    }
    cf[6] = 6.0;

#ifndef DISABLE_MPI
    if (NOFACE != 6) {
      printf("Error : NOFACE !=6\n");
      exit(0);
    }
    for (int i=0; i<NOFACE*2; i++) req[i] = MPI_REQUEST_NULL;
#endif
  }


  // デストラクタ
  ~CZ() {}


public:
  int Init(int argc, char **argv);
  int Loop();
  int Post();
  void debug(int m_mode) {
      debug_mode = m_mode;
  }



private:
  bool Comm_S(REAL_TYPE* sa, const int gc, const string label="");
  bool Comm_V(REAL_TYPE* va, const int gc, const string label="");
  bool Comm_SUM_1(int* var, const string label="");
  bool Comm_SUM_1(double* var, const string label="");
  bool Comm_MIN_1(double* var, const string label="");
  bool Comm_MAX_1(double* var, const string label="");
  bool Comm_SUM_2(double* var1, double* var2, const string label="");

  bool displayMemoryInfo(FILE* fp, double& G_mem, double L_mem, const char* str);


  // 計算する内点のインデクス範囲と点数
  double range_inner_index();




  int JACOBI(double& res, REAL_TYPE* X, REAL_TYPE* B);
  int PSOR  (double& res, REAL_TYPE* X, REAL_TYPE* B);
  int RBSOR (double& res, REAL_TYPE* X, REAL_TYPE* B);
  int PBiCGstab(double& res, REAL_TYPE* X, REAL_TYPE* B);

  /**
   * @brief Fdot for 1 array
   * @retval  内積値
   * @param [in]   x   vector1
   *
  double Fdot1(REAL_TYPE* x);


  /**
   * @brief Fdot for 2 arrays
   * @retval  内積値
   * @param [in]   x   vector1
   * @param [in]   y   vector2
   *
  double Fdot2(REAL_TYPE* x, REAL_TYPE* y);


  /**
   * @brief Preconditioner
   * @param [in,out] X  解ベクトル
   * @param [in]     B  RHS vector
   *
  void Preconditioner(REAL_TYPE* X, REAL_TYPE* B);
*/

  // タイミング測定区間にラベルを与えるラッパー
  void set_label(const string label, pm_lib::PerfMonitor::Type type, bool exclusive=true);


  // タイミング測定区間にラベルを与える
  void set_timing_label();


  /**
   * @brief タイミング測定開始
   * @param [in] key ラベル
   */
  inline void TIMING_start(const string key)
  {
#ifndef DISABLE_PMLIB
    // PMlib Intrinsic profiler
    PM.start(key);

    const char* s_label = key.c_str();

    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( s_label );
#elif defined __FX_FAPP
    fapp_start( s_label, 0, 0);
#endif

#endif // DISABLE_PMLIB
  }


  /**
   * @brief タイミング測定終了
   * @param [in] key             ラベル
   * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
   * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
   */
  inline void TIMING_stop(const string key, double flopPerTask=0.0, int iterationCount=1)
  {
#ifndef DISABLE_PMLIB
    // Venus FX profiler
    const char* s_label = key.c_str();

#if defined __K_FPCOLL
    stop_collection( s_label );
#elif defined __FX_FAPP
    fapp_stop( s_label, 0, 0);
#endif

    // PMlib Intrinsic profiler
    PM.stop(key, flopPerTask, (unsigned)iterationCount);
#endif // DISABLE_PMLIB
  }




void setParallelism()
{
  if( numProc > 1 )
  {
    if ( numThreads > 1 )
    {
      Parallel_str = "Hybrid";
    }
    else
    {
      Parallel_str = "FlatMPI";
    }
  }
  else
  {
    if ( numThreads > 1 )
    {
      Parallel_str = "OpenMP";
    }
    else
    {
      Parallel_str = "Serial";
    }
  }
}

/**
 * @brief メモリ使用量を表示する
 * @param [in] mode     処理モード
 * @param [in] Memory   必要メモリ量
 * @param [in] l_memory local
 * @param [in] fp       ファイルポインタ
 */
void MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp);


std::string GetHostName()
{
  char name[512];
  memset(name, 0x00, sizeof(char)*512);
  if( gethostname(name, 512) != 0 ) return std::string("");
  return std::string(name);
}

};


#endif // _CZ_H_
