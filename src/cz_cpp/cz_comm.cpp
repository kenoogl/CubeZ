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

// #################################################################
/*
 * @brief スカラー配列の同期
 * @param [in,out] sa     Scalar array
 * @param [in]     gc     通信するガイドセル幅
 * @param [in]     label  PMlibラベル
 * @retval true/false
 * @note 同期通信的な利用
 */
bool CZ::Comm_S(REAL_TYPE* sa, int gc, const string label)
{
  if ( numProc == 1 ) return true;

  bool flag = true;
  if (!label.empty()) TIMING_start(label);

#ifndef DISABLE_MPI
  if ( !CM.Comm_S_nonblocking(sa, gc, req) ) flag=false;
  if ( !CM.Comm_S_wait_nonblocking(sa, gc, req) ) flag=false;
#endif

  if (!label.empty()) TIMING_stop(label, comm_size);

  return (flag)?true:false;
}


// #################################################################
/*
 * @brief ベクトル配列の同期
 * @param [in,out] va     Vector array
 * @param [in]     gc     通信するガイドセル幅
 * @param [in]     label  PMlibラベル
 * @retval true/false
 * @note 同期通信的な利用
 */
bool CZ::Comm_V(REAL_TYPE* va, int gc, const string label)
{
  if ( numProc == 1 ) return true;

  bool flag = true;
  if (!label.empty()) TIMING_start(label);

#if 0
#ifndef DISABLE_MPI
  // 各成分の先頭アドレス
  size_t p_u = 0;
  size_t p_v = (size[0]+2*GUIDE) * (size[1]+2*GUIDE) * (size[2]+2*GUIDE);
  size_t p_w = (size[0]+2*GUIDE) * (size[1]+2*GUIDE) * (size[2]+2*GUIDE) * 2;

  if ( !CM.Comm_S_nonblocking(&va[p_u], gc, req) ) return false;
  if ( !CM.Comm_S_wait_nonblocking(&va[p_u], gc, req) ) return false;
  if ( MPI_SUCCESS != MPI_Barrier(MPI_COMM_WORLD) ) return false;

  if ( !CM.Comm_S_nonblocking(&va[p_v], gc, req) ) return false;
  if ( !CM.Comm_S_wait_nonblocking(&va[p_v], gc, req) ) return false;
  if ( MPI_SUCCESS != MPI_Barrier(MPI_COMM_WORLD) ) return false;

  if ( !CM.Comm_S_nonblocking(&va[p_w], gc, req) ) return false;
  if ( !CM.Comm_S_wait_nonblocking(&va[p_w], gc, req) ) return false;
  if ( MPI_SUCCESS != MPI_Barrier(MPI_COMM_WORLD) ) return false;
#endif
#endif


#ifndef DISABLE_MPI
  if ( !CM.Comm_V_nonblocking(va, gc, req) ) flag=false;
  if ( !CM.Comm_V_wait_nonblocking(va, gc, req) ) flag=false;
#endif

  if (!label.empty()) TIMING_stop(label, comm_size*3.0);
  return (flag)?true:false;
}


// #################################################################
/*
 * @brief int型1変数のAllreduce
 * @param [in,out] var     対象変数
 * @param [in]     label   PMlibラベル
 * @retval true/false
 */
bool CZ::Comm_SUM_1(int* var, const string label)
{
  if ( numProc == 1 ) return true;

#ifndef DISABLE_MPI
  int tmp = *var;
  bool flag = true;

  if (!label.empty()) TIMING_start(label);
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    var,
                                    1,
                                    MPI_INT,
                                    MPI_SUM,
                                    MPI_COMM_WORLD) ) flag=false;
  if (!label.empty()) TIMING_stop(label, 2.0*numProc*sizeof(int));
  return (flag)?true:false;
#endif
}


// #################################################################
/*
 * @brief double型1変数のAllreduce
 * @param [in,out] var     対象変数
 * @param [in]     label   PMlibラベル
 * @retval true/false
 */
bool CZ::Comm_SUM_1(double* var, const string label)
{
  if ( numProc == 1 ) return true;

#ifndef DISABLE_MPI
  double tmp = *var;
  bool flag = true;

  if (!label.empty()) TIMING_start(label);
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    var,
                                    1,
                                    MPI_DOUBLE,
                                    MPI_SUM,
                                    MPI_COMM_WORLD) ) flag=false;
  if (!label.empty()) TIMING_stop(label, 2.0*numProc*sizeof(double));
  return (flag)?true:false;
#endif
}

// #################################################################
/*
 * @brief double型1変数のAllreduce
 * @param [in,out] var     対象変数
 * @param [in]     label   PMlibラベル
 * @retval true/false
 */
bool CZ::Comm_MIN_1(double* var, const string label)
{
  if ( numProc == 1 ) return true;

#ifndef DISABLE_MPI
  double tmp = *var;
  bool flag = true;

  if (!label.empty()) TIMING_start(label);
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    var,
                                    1,
                                    MPI_DOUBLE,
                                    MPI_MIN,
                                    MPI_COMM_WORLD) ) flag=false;
  if (!label.empty()) TIMING_stop(label, 2.0*numProc*sizeof(double));
  return (flag)?true:false;
#endif
}

// #################################################################
/*
 * @brief double型1変数のAllreduce
 * @param [in,out] var     対象変数
 * @param [in]     label   PMlibラベル
 * @retval true/false
 */
bool CZ::Comm_MAX_1(double* var, const string label)
{
  if ( numProc == 1 ) return true;

#ifndef DISABLE_MPI
  double tmp = *var;
  bool flag = true;

  if (!label.empty()) TIMING_start(label);
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    var,
                                    1,
                                    MPI_DOUBLE,
                                    MPI_MAX,
                                    MPI_COMM_WORLD) ) flag=false;
  if (!label.empty()) TIMING_stop(label, 2.0*numProc*sizeof(double));
  return (flag)?true:false;
#endif
}

// #################################################################
/*
 * @brief double型2変数のAllreduce
 * @param [in,out] var1    対象変数
 * @param [in,out] var2    対象変数
 * @param [in]     label   PMlibラベル
 * @retval true/false
 */
bool CZ::Comm_SUM_2(double* var1, double* var2, const string label)
{
  if ( numProc == 1 ) return true;

#ifndef DISABLE_MPI
  double buf[2], tmp[2];
  tmp[0] = buf[0] = *var1;
  tmp[1] = buf[1] = *var2;
  bool flag = true;

  if (!label.empty()) TIMING_start(label);
  if ( MPI_SUCCESS != MPI_Allreduce(tmp,
                                    buf,
                                    2,
                                    MPI_DOUBLE,
                                    MPI_SUM,
                                    MPI_COMM_WORLD) ) flag=false;
  if (!label.empty()) TIMING_stop(label, 4.0*numProc*sizeof(double));
  *var1 = buf[0];
  *var2 = buf[1];

  return (flag)?true:false;
#endif
}
