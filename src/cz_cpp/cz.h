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
#include <string.h>

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

#ifndef DISABLE_PMLIB
#include <PerfMonitor.h>
#endif

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> added by nVidia 2019-08-06 for profiling?
#ifdef USE_NVTX

#include "nvToolsExt.h"

const uint32_t colors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000, 0xffffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

#define PUSH_RANGE(name,cid) { \
int color_id = cid; \
color_id = color_id%num_colors;\
nvtxEventAttributes_t eventAttrib = {0}; \
eventAttrib.version = NVTX_VERSION; \
eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
eventAttrib.colorType = NVTX_COLOR_ARGB; \
eventAttrib.color = colors[color_id]; \
eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
eventAttrib.message.ascii = name; \
nvtxRangePushEx(&eventAttrib); \
}
#define POP_RANGE nvtxRangePop();

#else

#define PUSH_RANGE(name,cid)
#define POP_RANGE

#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<< nVidia


using namespace std;


class CZ : public DomainInfo {

private:

  double comm_size;
  int debug_mode;          ///< if 1, Debug mode
  int ItrMax;              ///< 最大反復回数
  int ls_type;             ///< 線形ソルバ種類
  int pc_type;             ///< 前処理ソルバ種類
  double eps;              ///< convergence criteria
  REAL_TYPE ac1;           ///< acceleration coef.
  double res_normal;       ///< 全計算点数
  REAL_TYPE cf[7];         ///< 係数
  std::string precon;      ///< 前処理文字列
  int SW_maf;
  int SW_esa;


  int order_of_PM_key;     ///< PMlib用の登録番号カウンタ < PM_NUM_MAX
  string Parallel_str;     ///< 並列モードの文字列

#ifndef DISABLE_PMLIB
  // PMlib
  pm_lib::PerfMonitor PM;  ///< 性能モニタクラス
#endif

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
  
  REAL_TYPE* xc; ///< 格子
  REAL_TYPE* yc;
  REAL_TYPE* zc;
  REAL_TYPE* vrtmp; ///< temporary for supressing vector reduction
  REAL_TYPE* pvt; ///< 行の最大要素

  REAL_TYPE* WA;
  REAL_TYPE* WC;
  REAL_TYPE* WD;
  REAL_TYPE* WAA;
  REAL_TYPE* WCC;
  REAL_TYPE* WDD;
  
  REAL_TYPE* SA;
  REAL_TYPE* SC;
  REAL_TYPE* SD;

public:
  // コンストラクタ
  CZ()
  {
    order_of_PM_key = 0;
    comm_size = 0.0;
    debug_mode = 0;
    ItrMax = 0;
    ls_type = 0;
    pc_type = 0;
    eps = 1.0e-5;
    ac1 = 0.0;
    res_normal = 0.0;
    SW_maf = 0;
    SW_esa = 0;
    
    
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
  int Evaluate(int argc, char **argv);

  void debug(int m_mode) {
      debug_mode = m_mode;
  }

/*
  void print_m256(__m256 x) {
    printf("%f %f %f %f %f %f ^%f %f\n",
    x[7], x[6], x[5], x[4],
    x[3], x[2], x[1], x[0]);
  }
*/

  // #################################################################
  /**
   * @brief S3D配列のアロケート
   * @param [in] sz 配列サイズ
   * @ret pointer
   */
  template <typename T>
  T* czAllocR_S3D(const int* sz, T type)
  {
    if ( !sz ) return NULL;
    
    size_t dims[3], nx;
    
    dims[0] = (size_t)(sz[0] + 2*GUIDE);
    dims[1] = (size_t)(sz[1] + 2*GUIDE);
    dims[2] = (size_t)(sz[2] + 2*GUIDE);
    
    nx = dims[0] * dims[1] * dims[2];
    T* var = new T[nx];
    
#ifndef __NEC__
#pragma omp parallel for schedule(static)
#endif
#ifdef _OPENACC
#pragma acc kernels
#endif
    for (int i=0; i<nx; i++) var[i]=0;
    
    return var;
  }
  
  // #################################################################
  template <typename T>
  T* czAllocR(const int sz, T type)
  {
    if ( !sz ) return NULL;
    
    size_t nx = sz;
    T* var = new T[nx];
    
#ifndef __NEC__
#pragma omp parallel for schedule(static)
#endif
#ifdef _OPENACC
#pragma acc kernels
#endif
    for (int i=0; i<nx; i++) var[i]=0;
    
    return var;
  }
  
  // #################################################################
  template <typename T>
  void czDelete(T* ptr)
  {
    if (ptr) {
      delete [] ptr;
      ptr = NULL;
    }
  }
  

  // #################################################################
  template <typename T>
  void check_align (T* var, std::string str)
  {
    long int li = (long int)var;
    printf("\t%10s : addrs= %08x align(4=%2u, 8=%2u, 16=%2u, 32=%2d, 64=%2d)\n\n",
           str.c_str(),li, li%4, li%8, li%16, li%32, li%64);
  }

  // #################################################################
  // cz_tdma.cpp

  void tdma(int nx,
            REAL_TYPE* d,
            REAL_TYPE* a,
            REAL_TYPE* b,
            REAL_TYPE* c,
            REAL_TYPE* w);

  void tdma1(int nx,
            REAL_TYPE* d,
            REAL_TYPE* a,
            REAL_TYPE* c,
            REAL_TYPE* w);
  
  
  // @param [in] n 方程式の次元数
  // @retval nを超える最小の2べき数の乗数
  int getNumStage(int n) {
    int b = 1;
    for (int i=1; i<20; i++) { // とりあえず20乗まで試しておけば十分
      b *= 2;
      if (n<b) return i;
    }
    return -1;
  }


  void pcr(const int nx,
           const int pn,
           REAL_TYPE* d,
           REAL_TYPE* a,
           REAL_TYPE* c,
           REAL_TYPE* d1,
           REAL_TYPE* a1,
           REAL_TYPE* c1,
           double& flop);

  void pcr_kernel(const int nx,
                  const int s,
                  REAL_TYPE* d,
                  REAL_TYPE* a,
                  REAL_TYPE* c,
                  REAL_TYPE* dn,
                  REAL_TYPE* an,
                  REAL_TYPE* cn,
                  double& flop);
  
  void printArray(int nx, REAL_TYPE* a, char* s);
  
  
private:
  inline static void cIndex(int& i, int& j,
                const int l, const int ni, const int is, const int js) {
    int jj = l / ni;
    int ii = l - jj*ni;
    j = jj + js - 1;
    i = ii + is - 1;
  };

  inline static void fIndex(int& i, int& j,
                const int l, const int ni, const int is, const int js) {
    int jj = l / ni;
    int ii = l - jj*ni;
    j = jj + js;
    i = ii + is;
  };


  inline void matx2(REAL_TYPE* d, REAL_TYPE a, REAL_TYPE c)
  {
    /*   Ax = d    8+8 fp

         A=|1 c|   x={x1, x2}, d={d1, d2}
           |a 1| ,
     */

     REAL_TYPE J = 1.0 / (1.0 - a * c);
     REAL_TYPE d1 = d[0];
     REAL_TYPE d2 = d[1];
     d[0] = (d1 - c * d2) * J;
     d[1] = (d2 - a * d1) * J;
  }

  inline void matx3(REAL_TYPE* d, REAL_TYPE* a, REAL_TYPE* c)
  {
    /*  |  1 c1  0 |  8+24 fp
        | a2  1 c2 |
        |  0 a3  1 |
    */

     REAL_TYPE a2 = a[1];
     REAL_TYPE a3 = a[2];
     REAL_TYPE c1 = c[0];
     REAL_TYPE c2 = c[1];
     REAL_TYPE d1 = d[0];
     REAL_TYPE d2 = d[1];
     REAL_TYPE d3 = d[2];
     REAL_TYPE J = 1.0 / (1.0 - c2 * a3 - c1 * a2);
     d[0] = ( d1 * (3.0-c2*a3) - c1*d2 ) * J;
     d[1] = (1.0 - d1*a2 + 2.0*d2 - c2*d3) * J;
     d[2] = (1.0 + 2.0*d3 - a3*d2 - a2*c1) * J;
  }



  bool Comm_S(REAL_TYPE* sa, const int gc, const string label="");
  bool Comm_V(REAL_TYPE* va, const int gc, const string label="");
  
  bool Comm_SUM_1(int* var, const string label="");
  bool Comm_SUM_1(double* var, const string label="");
  bool Comm_SUM_1(float*  var, const string label="");
  
  bool Comm_MIN_1(double* var, const string label="");
  bool Comm_MIN_1(float*  var, const string label="");
  bool Comm_MAX_1(double* var, const string label="");
  bool Comm_MAX_1(float*  var, const string label="");
  bool Comm_SUM_2(double* var1, double* var2, const string label="");
  bool Comm_SUM_2(float*  var1, float*  var2, const string label="");

  bool displayMemoryInfo(FILE* fp, double& G_mem, double L_mem, const char* str);


  // 計算する内点のインデクス範囲と点数
  double range_inner_index();


  int JACOBI(double& res,
             REAL_TYPE* X,
             REAL_TYPE* B,
             const int itrMax,
             double& flop,
             int s_type,
             bool converge_check=true);

  int PSOR  (double& res,
             REAL_TYPE* X,
             REAL_TYPE* B,
             const int itrMax,
             double& flop,
             int s_type,
             bool converge_check=true);

  int RBSOR (double& res,
             REAL_TYPE* X,
             REAL_TYPE* B,
             const int itrMax,
             double& flop,
             int s_type,
             bool converge_check=true);

  int PBiCGSTAB(double& res,
                REAL_TYPE* X,
                REAL_TYPE* B,
                double& flop,
                int s_type);

  int LSOR_PCR(double& res,
               REAL_TYPE* X,
               REAL_TYPE* B,
               const int itr_max,
               double& flop,
               int s_type,
               bool converge_check=true);
  
  int LSOR_PCR_EDA(double& res,
                REAL_TYPE* X,
                REAL_TYPE* B,
                const int itr_max,
                double& flop,
                int s_type,
                bool converge_check=true);
  
  int LSOR_PCR_ESA(double& res,
                   REAL_TYPE* X,
                   REAL_TYPE* B,
                   const int itr_max,
                   double& flop,
                   int s_type,
                   bool converge_check=true);
  
  int LSOR_PCR_RB(double& res,
              REAL_TYPE* X,
              REAL_TYPE* B,
              const int itr_max,
              double& flop,
              int s_type,
              bool converge_check=true);

  int LSOR_PCR_RB_ESA(double& res,
                  REAL_TYPE* X,
                  REAL_TYPE* B,
                  const int itr_max,
                  double& flop,
                  int s_type,
                  bool converge_check=true);
  
  REAL_TYPE Fdot1(REAL_TYPE* x, double& flop);

  REAL_TYPE Fdot2(REAL_TYPE* x, REAL_TYPE* y, double& flop);

  void Preconditioner(REAL_TYPE* xx,
                      REAL_TYPE* bb,
                      double& flop,
                      int s_type);

  void setStrPre();
  void setLS(char* q, char* fname);


#ifndef DISABLE_PMLIB
  // タイミング測定区間にラベルを与えるラッパー
  void set_label(const string label, pm_lib::PerfMonitor::Type type, bool exclusive=true);
#endif


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

  std::string printMethod(int type);

};


#endif // _CZ_H_
