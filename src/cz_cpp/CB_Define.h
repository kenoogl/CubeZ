#ifndef _CB_DEFINE_H_
#define _CB_DEFINE_H_

/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2018 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/**
 * @file   CB_Define.h
 * @brief  CBrick Definition Header
 */

#include <float.h>
#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// precision
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#else
/** 実数型の指定
 * - デフォルトでは、REAL_TYPE=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   REAL_TYPE=doubleになる
 */
#define REAL_TYPE float
#endif

#define _SIZE_DOUBLE_ 8

#ifndef _DIAGONAL_COMM
#define NOFACE 6
#else
#define NOFACE 26
#endif

#define ON 1
#define OFF 0


enum DIRection {
  I_minus=0,
  I_plus,
  J_minus,
  J_plus,
  K_minus,
  K_plus,
#ifdef _DIAGONAL_COMM
  // edge
  E_mYmZ,
  E_pYmZ,
  E_mYpZ,
  E_pYpZ,
  E_mXmZ,
  E_pXmZ,
  E_mXpZ,
  E_pXpZ,
  E_mXmY,
  E_pXmY,
  E_mXpY,
  E_pXpY,
  // corner
  C_mXmYmZ,
  C_pXmYmZ,
  C_mXpYmZ,
  C_pXpYmZ,
  C_mXmYpZ,
  C_pXmYpZ,
  C_mXpYpZ,
  C_pXpYpZ,
#endif
};





/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ [C version]
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
 #define _IDX_S3D(_I,_J,_K,_NI,_NJ,_VC) \
 ( (_K+(_VC)) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) \
 + (_J+(_VC)) * (_NI+2*(_VC)) \
 + (_I+(_VC)) \
 )


/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記
 *  @param [in] _I  i方向インデクス (0,NI-1)
 *  @param [in] _J  j方向インデクス (0, NJ-1)
 *  @param [in] _K  k方向インデクス (0, NK-1)
 *  @param [in] _NI i方向インデクスサイズ NI
 *  @param [in] _NJ j方向インデクスサイズ NJ
 *  @param [in] _NK k方向インデクスサイズ NK
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
 #define _F_IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
 ( (_K+(_VC)-1) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) \
 + (_J+(_VC)-1) * (_NI+2*(_VC)) \
 + (_I+(_VC)-1) \
 )


/** 3次元インデクス(i,j,k,l) -> 1次元インデクス変換マクロ [C version]
 *  @param [in] _I  i方向インデクス (1,NI)
 *  @param [in] _J  j方向インデクス (1,NJ)
 *  @param [in] _K  k方向インデクス (1,NK)
 *  @param [in] _L  ベクトル成分インデクス {1,2,3}
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
 #define _F_IDX_V3D(_I,_J,_K,_L,_NI,_NJ,_NK,_VC) \
 ( ((_L)-1) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) * (_NK+2*(_VC))  \
 + (_K+(_VC)-1) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) \
 + (_J+(_VC)-1) * (_NI+2*(_VC)) \
 + (_I+(_VC)-1) \
 )

 /** 3次元インデクス(i,j,k,l) -> 1次元インデクス変換マクロ [C version]
  *  @param [in] _I  i方向インデクス
  *  @param [in] _J  j方向インデクス
  *  @param [in] _K  k方向インデクス
  *  @param [in] _L  ベクトル成分インデクス {0,1,2}
  *  @param [in] _NI i方向インデクスサイズ
  *  @param [in] _NJ j方向インデクスサイズ
  *  @param [in] _NK k方向インデクスサイズ
  *  @param [in] _VC 仮想セル数
  *  @return 1次元インデクス
  */
  #define _IDX_V3D(_I,_J,_K,_L,_NI,_NJ,_NK,_VC) \
  ( (_L) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) * (_NK+2*(_VC))  \
  + (_K+(_VC)) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) \
  + (_J+(_VC)) * (_NI+2*(_VC)) \
  + (_I+(_VC)) \
  )


#define stamped_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#define stamped_fprintf fprintf(fp, "%s (%d):  ",__FILE__, __LINE__), fprintf
#define mark() printf("%s (%d) [%d]:\n",__FILE__, __LINE__, myRank)
#define Hostonly_ if(myRank==0)

#define Exit(x) \
((void)printf("exit at %s:%u\n", __FILE__, __LINE__), exit((x)))

#endif // _CB_DEFINE_H_
