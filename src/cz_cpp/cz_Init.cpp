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
#include "czVersion.h"
#include <climits>

int CZ::Init(int argc, char **argv)
{
  int div_type = 0;         ///< 分割指定 (0-自動、1-指定)
  double G_Memory = 0.0;    ///< 計算に必要なメモリ量（グローバル）
  double L_Memory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）

  procGrp = 0;

#ifdef _OPENMP
	char* c_env = std::getenv("OMP_NUM_THREADS");
	if (c_env == NULL) {
		omp_set_num_threads(1);	// OMP_NUM_THREADS was not defined. set as 1
	}
	numThreads  = omp_get_max_threads();
#else
	numThreads  = 1;
#endif


#ifndef DISABLE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#else
  myRank = 0;
  numProc= 1;
#endif


  double flop_count    = 0.0;  ///< flops計算用

  G_size[0] = atoi(argv[1]);
  G_size[1] = atoi(argv[2]);
  G_size[2] = atoi(argv[3]);

  // 分割数が指定されている場合
  if (argc == 9) {
    div_type=1;
    G_div[0] = atoi(argv[6]);
    G_div[1] = atoi(argv[7]);
    G_div[2] = atoi(argv[8]);
  }

  // 等方
  pitch[0] = pitch[1] = pitch[2] = 1.0/(REAL_TYPE)(G_size[0]-1);


  // 分割数のチェック
  if ( div_type == 1 && G_div[0]*G_div[1]*G_div[2] != numProc) {
    printf("\tThe number of proceees does not agree with the division size.\n");
    return 0;
  }



  // 領域分割

  if ( numProc > 1 )
  {
#ifndef DISABLE_MPI

    Hostonly_ printf("\n++++++++++++++++++++++++++++ >> CBrick\n\n");

    int gc = GUIDE;
    D.setSubDomain(G_size, gc, numProc, myRank, procGrp,
                    MPI_COMM_WORLD, "node", "Findex");

    // 分割数指定の場合
    if (div_type == 1)
    {
      if ( !D.setDivision(G_div) ) return 0;
    }

    // 引数 0 は3次元分割
    if ( !D.findOptimalDivision(0) ) return 0;

    if ( !D.createRankTable() ) return 0;

    // 領域分割数取得 （自動分割の場合）
    D.getGlobalDivision(G_div);

    //自ランクのHeadIndexの取得
    D.getLocalHead(head);

    //自ランクの格子数取得
    D.getLocalSize(size);

    // 自ランクの基点座標 ノード点座標
    // G_origin[]は 0.0
    origin[0] = G_origin[0] + (head[0]-1)*pitch[0];
    origin[1] = G_origin[1] + (head[1]-1)*pitch[1];
    origin[2] = G_origin[2] + (head[2]-1)*pitch[2];

    //自ランクの隣接ランク番号を取得
    D.getCommTable(nID);


    // 通信クラス設定
    if ( !CM.setBrickComm(size, gc, MPI_COMM_WORLD, nID, "node") ) {
      stamped_printf("\tBrickComm settng error.\n");
      return 0;
    }

    // 通信バッファ確保
    if  ( !CM.init(1) ) {
      stamped_printf("\tBrickComm initialize error.\n");
      return 0;
    }

    Hostonly_ printf("++++++++++++++++++++++++++++ << CBrick\n\n");

#endif // DISABLE_MPI
  }
  else // Serial
  {
    G_div[0] = 1;
    G_div[1] = 1;
    G_div[2] = 1;

    size[0] = G_size[0];
    size[1] = G_size[1];
    size[2] = G_size[2];

    head[0] = 1;
    head[1] = 1;
    head[2] = 1;

    origin[0] = G_origin[0];
    origin[1] = G_origin[1];
    origin[2] = G_origin[2];
  }


  // 通信量 双方向 x ２面
  comm_size = (double)( (size[0]+2*GUIDE) * (size[1]+2*GUIDE)
                      + (size[1]+2*GUIDE) * (size[2]+2*GUIDE)
                      + (size[0]+2*GUIDE) * (size[2]+2*GUIDE) )
                      * 2.0 * 2.0 * sizeof(REAL_TYPE);


  // 線形ソルバ
  char* q = argv[4];
  char fname[20];
  memset(fname, 0, sizeof(char)*20);

  if ( !strcasecmp(q, "jacobi") ) {
    ls_type = LS_JACOBI;
    strcpy(fname, "jacobi.txt");
  }

  else if ( !strcasecmp(q, "psor") ) {
    ls_type = LS_PSOR;
    strcpy(fname, "psor.txt");
  }

  else if ( !strcasecmp(q, "lsor") ) {
    ls_type = LS_LSOR;
    strcpy(fname, "lsor.txt");
  }

  else if ( !strcasecmp(q, "sor2sma") ) {
    ls_type = LS_SOR2SMA;
    strcpy(fname, "sor2sma.txt");
  }

  else if ( !strcasecmp(q, "pbicgstab") ) {
    ls_type = LS_BICGSTAB;
    strcpy(fname, "pbicgstab.txt");
  }

  else{
    printf("Invalid solver\n");
    exit(0);
  }


  // history title
  Hostonly_ {
    if ( !(fph=fopen(fname, "w")) )
    {
      printf("\tSorry, can't open file.\n");
      assert(0);
    }

    fprintf(fph, "Itration      Residual\n");
  }



  // 配列のアロケート
  double array_size = (size[0]+2*GUIDE) * (size[1]+2*GUIDE) * (size[2]+2*GUIDE);

  L_Memory += ( array_size * 3 ) * (double)sizeof(REAL_TYPE);

  if( (RHS = Alloc_Real_S3D(size)) == NULL ) return false;
  if( (P   = Alloc_Real_S3D(size)) == NULL ) return false;
  if( (WRK = Alloc_Real_S3D(size)) == NULL ) return false;
  //if( (MSK = Alloc_Real_S3D(size)) == NULL ) return false;

  if (debug_mode == 1) {
    L_Memory += ( array_size * 1 ) * (double)sizeof(REAL_TYPE);

    //if( (EXS = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (ERR = Alloc_Real_S3D(size)) == NULL ) return false;
  }


  if (ls_type == LS_BICGSTAB)
  {
    L_Memory += ( array_size * 9 ) * (double)sizeof(REAL_TYPE);

    if( (pcg_p  = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_p_ = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_r  = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_r0 = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_q  = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_s  = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_s_ = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_t  = Alloc_Real_S3D(size)) == NULL ) return false;
    if( (pcg_t_ = Alloc_Real_S3D(size)) == NULL ) return false;
  }


  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf(    "\n----------\n\n");
  }

  G_Memory = L_Memory;
  if ( !displayMemoryInfo(stdout, G_Memory, L_Memory, "Solver") ) return 0;



  // 計算するインデクス範囲の決定
  double sum_r = range_inner_index();
  if ( !Comm_SUM_1(&sum_r) ) return 0;
  res_normal = 1.0/(double)sum_r;
  Hostonly_ printf("Sum of inner = %e\n", sum_r);

  // 最大反復回数
  ItrMax = atoi(argv[5]);

  // Apply BC
  int gc = GUIDE;
  bc_(size, &gc, P, pitch, origin, nID);


  // mask
  //init_mask_(MSK, size, &gc);

  // source term
  src_dirichlet_(RHS, size, innerFidx, &gc, pitch);

  if ( !Comm_S(RHS, 1) ) return 0;



  // タイミング測定の初期化
  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( myRank );
  setParallelism();
  PM.setParallelMode(Parallel_str, numThreads, numProc);
  set_timing_label();

  return 1;
}



// #################################################################
/* @brief 計算する内点のインデクス範囲と点数
 * @retval 計算点数
 */
double CZ::range_inner_index()
{
  int ist, jst, kst, ied, jed, ked;

  ist = 2;
  jst = 2;
  kst = 2;
  ied = size[0];
  jed = size[1];
  ked = size[2];

  if (nID[I_plus] < 0)  ied = size[0] - 1;
  if (nID[J_plus] < 0)  jed = size[1] - 1;
  if (nID[K_plus] < 0)  ked = size[2] - 1;

  innerFidx[I_minus] = ist;
  innerFidx[I_plus]  = ied;
  innerFidx[J_minus] = jst;
  innerFidx[J_plus]  = jed;
  innerFidx[K_minus] = kst;
  innerFidx[K_plus]  = ked;

  double sum = 0.0;
  sum = (double)(ied-ist+1)
      * (double)(jed-jst+1)
      * (double)(ked-kst+1);

  printf("[%d] :%3d %3d %3d %3d %3d %3d : %3d %3d %3d\n",
    myRank, ist, ied, jst, jed, kst, ked,
    size[0], size[1], size[2]);

  return sum;
}

// #################################################################
/* @brief メモリ消費情報を表示
 * @param [in]     fp    ファイルポインタ
 * @param [in]     G_mem グローバルメモリサイズ
 * @param [in]     L_mem ローカルメモリサイズ
 * @param [in]     str   表示用文字列
 */
bool CZ::displayMemoryInfo(FILE* fp, double& G_mem, double L_mem, const char* str)
{
  if ( !Comm_SUM_1(&G_mem) ) return false;

  Hostonly_
  {
    MemoryRequirement(str, G_mem, L_mem, fp);
    fprintf(fp, "\n\n");
  }

  return true;
}


// #################################################################
// メモリ使用量を表示する
void CZ::MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp)
{
  const double mem = Memory;
  const double lmem= l_memory;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional

  fprintf (fp,"\t>> Memory required for %s : ", mode);

  // Global memory
  fprintf (fp," Global=");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }

  // Local memory
  fprintf (fp," : Local=");
  if ( lmem > PB ) {
    fprintf (fp,"%6.2f (PB)\n", lmem / PB *factor);
  }
  else if ( lmem > TB ) {
    fprintf (fp,"%6.2f (TB)\n", lmem / TB *factor);
  }
  else if ( lmem > GB ) {
    fprintf (fp,"%6.2f (GB)\n", lmem / GB *factor);
  }
  else if ( lmem > MB ) {
    fprintf (fp,"%6.2f (MB)\n", lmem / MB *factor);
  }
  else if ( lmem > KB ) {
    fprintf (fp,"%6.2f (KB)\n", lmem / KB *factor);
  }
  else if ( lmem <= KB ){
    fprintf (fp,"%6.2f (B)\n", lmem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)\n", (int)(lmem *factor) );
  }

  fflush(fp);
}


// #################################################################
/**
 * @brief タイミング測定区間にラベルを与えるラッパー
 * @param [in] label     ラベル
 * @param [in] type      測定対象タイプ(COMM or CALC)
 * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
 */
void CZ::set_label(const string label, pm_lib::PerfMonitor::Type type, bool exclusive)
{
  // 登録個数のチェック
  order_of_PM_key++;

  if ( order_of_PM_key > PM_NUM_MAX )
  {
    printf("\tThe number of labels for Performance monitor goes over limit.\n");
    exit(0);
  }

  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
  {
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }

  // Performance Monitorへの登録
  PM.setProperties(label, type, exclusive);
}


// #################################################################
/**
 * @brief タイミング測定区間にラベルを与える
 */
void CZ::set_timing_label()
{
using namespace pm_lib;

  set_label("Comm_RHS",            PerfMonitor::COMM);

  set_label("Dot1"   ,             PerfMonitor::CALC);
  set_label("Dot2"   ,             PerfMonitor::CALC);
  set_label("A_R_Dot",             PerfMonitor::COMM);
  set_label("Blas_Copy",           PerfMonitor::CALC);
  set_label("Blas_Clear",          PerfMonitor::CALC);
  set_label("Blas_Residual",       PerfMonitor::CALC);
  set_label("Blas_BiCG_1",         PerfMonitor::CALC);
  set_label("Blas_BiCG_2",         PerfMonitor::CALC);
  set_label("Blas_AX",             PerfMonitor::CALC);
  set_label("Blas_TRIAD",          PerfMonitor::CALC);


  // 非排他, 計算
  set_label("Jacobi",           PerfMonitor::CALC, false);
  set_label("PSOR",             PerfMonitor::CALC, false);
  set_label("SOR2_SMA",         PerfMonitor::CALC, false);
  set_label("PBiCGstab",        PerfMonitor::CALC, false);
  set_label("BiCGstab",         PerfMonitor::CALC, false);

  set_label("JACOBI_kernel",    PerfMonitor::CALC, true);
  set_label("SOR_kernel",       PerfMonitor::CALC, true);
  set_label("SOR2SMA_kernel",   PerfMonitor::CALC, true);
  set_label("LSOR_kernel",      PerfMonitor::CALC, true);

  set_label("Comm_Poisson",     PerfMonitor::COMM);
  set_label("Comm_Res_Poisson", PerfMonitor::COMM);

  set_label("BoundaryCondition",PerfMonitor::CALC);
}
