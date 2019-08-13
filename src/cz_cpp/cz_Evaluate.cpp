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
#include <climits>

/*
  \nabla^2 p = 0 を解く
   Z=0, 1に境界条件を与える
*/


int CZ::Evaluate(int argc, char **argv)
{
  int div_type = 0;         ///< 分割指定 (0-自動、1-指定)
  double G_Memory = 0.0;    ///< 計算に必要なメモリ量（グローバル）
  double L_Memory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  int gc = GUIDE;

  procGrp = 0;
  
#if defined(_OPENACC)
  printf("OpenACC is defined\n");
#endif

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

  // pbicgstab preconditoner
  char* q = argv[4];

  if ( !strcasecmp(q, "pbicgstab") || !strcasecmp(q, "pbicgstab_maf") )
    {
      if (argc!=8 && argc!=11) {
        Hostonly_ printf("command line error : pbicgstab\n");
        exit(0);
      }
      precon = argv[7];
    }


  if (argc == 10) {
    div_type=1;
    G_div[0] = atoi(argv[7]);
    G_div[1] = atoi(argv[8]);
    G_div[2] = atoi(argv[9]);
  }

  if (argc == 11) {
    div_type=1;
    G_div[0] = atoi(argv[8]);
    G_div[1] = atoi(argv[9]);
    G_div[2] = atoi(argv[10]);
  }

  // Z方向を基準に等方
  pitch[0] = pitch[1] = pitch[2] = 1.0/(REAL_TYPE)(G_size[2]-1);
  


  // 分割数のチェック
  if ( div_type == 1 && G_div[0]*G_div[1]*G_div[2] != numProc) {
    printf("\tThe number of proceees does not agree with the division size.\n");
    return 0;
  }

  // 係数
  ac1 = atof(argv[6]);

  // 領域分割

  if ( numProc > 1 )
  {
#ifndef DISABLE_MPI

    Hostonly_ printf("\n++++++++++++++++++++++++++++ >> CBrick\n\n");


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
  char fname[20];
  memset(fname, 0, sizeof(char)*20);

  setLS(q, fname);
  

  printf("Iterative Mehtod = %s\n", printMethod(ls_type).c_str() );
  if (ls_type==LS_BICGSTAB || ls_type==LS_BICGSTAB_MAF)
  {
    printf("Preconditioner = %s\n", printMethod(pc_type).c_str() );
  }



  /* 逐次のみ、k方向を内側にしているので通信面を変更
  else if ( !strcasecmp(q, "lsor_simd") ) {

  }
  */


  // history title
  Hostonly_ {
    if ( !(fph=fopen(fname, "w")) )
    {
      printf("\tSorry, can't open file.\n");
      exit(0);
    }

    fprintf(fph, "Itration      Residual\n");
  }


  // 計算するインデクス範囲の決定
  double sum_r = range_inner_index();
  if ( !Comm_SUM_1(&sum_r) ) return 0;
  res_normal = 1.0/(double)sum_r;
  //Hostonly_ printf("Sum of inner = %e\n", sum_r);
  

  // 配列のアロケート
  double array_size = (size[0]+2*GUIDE) * (size[1]+2*GUIDE) * (size[2]+2*GUIDE);

  L_Memory += ( array_size * 3 ) * (double)sizeof(REAL_TYPE);

  // アロケートのためのダミー型
  REAL_TYPE var_type=0;
  
  // nVidiaTools 定義されなければ何もしない
  PUSH_RANGE("memory allocation", 1);

  if( (RHS = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  if( (P   = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  if( (WRK = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  if( (MSK = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  if( (pvt = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  
  if( (xc = czAllocR(size[0]+2*GUIDE, var_type)) == NULL ) return 0;
  if( (yc = czAllocR(size[1]+2*GUIDE, var_type)) == NULL ) return 0;
  if( (zc = czAllocR(size[2]+2*GUIDE, var_type)) == NULL ) return 0;
  
  
  if (SW_maf == 1 )
  {
    if( (vrtmp = czAllocR(size[2]+2*GUIDE, var_type)) == NULL ) return 0;
  }
  
  int tmpz = (size[2]+2*GUIDE);
  
  if( (WAA= czAllocR(tmpz, var_type)) == NULL ) return 0;
  if( (WCC= czAllocR(tmpz, var_type)) == NULL ) return 0;
  if( (WDD= czAllocR(tmpz, var_type)) == NULL ) return 0;
  if( (WA = czAllocR(tmpz, var_type)) == NULL ) return 0;
  if( (WC = czAllocR(tmpz, var_type)) == NULL ) return 0;
  if( (WD = czAllocR(tmpz, var_type)) == NULL ) return 0;

  if (debug_mode == 1) {
    L_Memory += ( array_size * 1 ) * (double)sizeof(REAL_TYPE);

    //if( (EXS = czAllocR_S3D(size)) == NULL ) return 0;
    if( (ERR = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  }

  //check_align(RHS, "rhs");


  if (ls_type == LS_BICGSTAB || ls_type == LS_BICGSTAB_MAF)
  {
    L_Memory += ( array_size * 9 ) * (double)sizeof(REAL_TYPE);

    if( (pcg_p  = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_p_ = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_r  = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_r0 = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_q  = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_s  = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_s_ = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_t  = czAllocR_S3D(size,var_type)) == NULL ) return 0;
    if( (pcg_t_ = czAllocR_S3D(size,var_type)) == NULL ) return 0;
  }
  
  // PCR用の配列確保
  if ( SW_esa == 1 )
  {
    int kst = innerFidx[K_minus];
    int ked = innerFidx[K_plus];
    int n = ked - kst + 1;
    int pn;
    
    // Nを超える最小の2べき数の乗数 pn
    if ( -1 == (pn=getNumStage(n))) {
      printf("error : number of stage\n");
      exit(0);
    }
    
    int s = pow(2, pn-2);
    int kk = ked - kst+ 2*s + 1;
    printf("s = %d  kk= %d\n", s, kk);
    
    // スレッド用
    if( (SA = czAllocR(kk, var_type)) == NULL ) return 0;
    if( (SC = czAllocR(kk, var_type)) == NULL ) return 0;
    if( (SD = czAllocR(kk, var_type)) == NULL ) return 0;
    L_Memory += ( kk * 3 ) * (double)sizeof(REAL_TYPE);
  }
  
  POP_RANGE; // nVIdiaTools
  

  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf(    "\n----------\n\n");
  }

  G_Memory = L_Memory;
  if ( !displayMemoryInfo(stdout, G_Memory, L_Memory, "Solver") ) return 0;



  // 最大反復回数
  ItrMax = atoi(argv[5]);
  
  
  // nVidiaTools
  PUSH_RANGE("initialize grid", 2);
  
  
  // 一次元格子  何か値をいれておく
#ifdef _OPENACC
#pragma acc kernels
#pragma acc loop independent
#endif
  for (int i=0; i<size[0]+2*GUIDE; i++)
  {
    xc[i] = (REAL_TYPE)(i-1) * pitch[0];
  }
  
#ifdef _OPENACC
#pragma acc kernels
#pragma acc loop independent
#endif
  for (int i=0; i<size[1]+2*GUIDE; i++)
  {
    yc[i] = (REAL_TYPE)(i-1) * pitch[1];
  }
  
#ifdef _OPENACC
#pragma acc kernels
#pragma acc loop independent
#endif
  for (int i=0; i<size[2]+2*GUIDE; i++)
  {
    zc[i] = (REAL_TYPE)(i-1) * pitch[2];
  }
  
  POP_RANGE;  // nVidiaTools
  
  // 行の最大値の逆数
  PUSH_RANGE("search_pivot", 3);
  search_pivot_(pvt, size, innerFidx, &gc, xc, yc, zc);
  POP_RANGE;
  

  // Apply BC
  PUSH_RANGE("bc_k", 4);
  bc_k_(size, &gc, P, pitch, origin, nID);
  POP_RANGE;
  
  if ( !Comm_S(P, 1) ) return 0;
  
  
  // source term >> ソース項ゼロ
  PUSH_RANGE("bc_k", 4);
  bc_k_(size, &gc, RHS, pitch, origin, nID);
  POP_RANGE;
  
  if ( !Comm_S(RHS, 1) ) return 0;
  
  
  PUSH_RANGE("imask_k", 5);
  imask_k_(MSK, size, innerFidx, &gc);
  POP_RANGE;


#ifndef DISABLE_PMLIB
  // タイミング測定の初期化
  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( myRank );
  setParallelism();
  PM.setParallelMode(Parallel_str, numThreads, numProc);
  set_timing_label();
#endif



  /////////////////////////////////////////////////////////////
  // Loop

  double res=0.0;
  int itr=0;
  double flop=0.0; // dummy


  PUSH_RANGE("main loop",6);
  switch (ls_type)
  {
    case LS_JACOBI:
    case LS_JACOBI_MAF:
      TIMING_start("JACOBI");
      if ( 0 == (itr=JACOBI(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
      TIMING_stop("JACOBI", flop);
      break;

    case LS_PSOR:
    case LS_PSOR_MAF:
      TIMING_start("PSOR");
      if ( 0 == (itr=PSOR(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
      TIMING_stop("PSOR", flop);
      break;

    case LS_SOR2SMA:
    case LS_SOR2SMA_MAF:
      TIMING_start("SOR2SMA");
      if ( 0 == (itr=RBSOR(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
      TIMING_stop("SOR2SMA", flop);
      break;

    case LS_BICGSTAB:
    case LS_BICGSTAB_MAF:
      TIMING_start("PBiCGSTAB");
      if ( 0 == (itr=PBiCGSTAB(res, P, RHS, flop, ls_type)) ) return 0;
      TIMING_stop("PBiCGSTAB", flop);
      break;

    case LS_PCR:
    case LS_PCR_MAF:
      TIMING_start("LSOR");
      if ( 0 == (itr=LSOR_PCR(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
      TIMING_stop("LSOR", flop);
      break;
    
    case LS_PCR_EDA:
    case LS_PCR_EDA_MAF:
    TIMING_start("LSOR");
    if ( 0 == (itr=LSOR_PCR_EDA(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
    TIMING_stop("LSOR", flop);
    break;
    
    case LS_PCR_ESA:
    case LS_PCR_ESA_MAF:
    TIMING_start("LSOR");
    if ( 0 == (itr=LSOR_PCR_ESA(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
    TIMING_stop("LSOR", flop);
    break;
      
    case LS_PCR_RB:
    case LS_PCR_RB_MAF:
      TIMING_start("LSOR");
      if ( 0 == (itr=LSOR_PCR_RB(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
      TIMING_stop("LSOR", flop);
      break;
      
    case LS_PCR_RB_ESA:
    case LS_PCR_RB_ESA_MAF:
      TIMING_start("LSOR");
      if ( 0 == (itr=LSOR_PCR_RB_ESA(res, P, RHS, ItrMax, flop, ls_type)) ) return 0;
      TIMING_stop("LSOR", flop);
      break;
      
    default:
      break;
  }
  POP_RANGE;
  

  Hostonly_ {
    printf("\n=================================\n");
    printf("Iter = %d  Res = %e\n", itr, res);
    printf("=================================\n");
  }


  /////////////////////////////////////////////////////////////
  // post

  Hostonly_ {
    if (!fph) fclose(fph);
  }

#ifndef DISABLE_PMLIB

  FILE *fp = NULL;

  Hostonly_
  {
    if ( !(fp=fopen("profiling.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
      Exit(0);
    }
  }


  // 測定結果の集計(gatherメソッドは全ノードで呼ぶこと)
  PM.gather();

  char str[100];
  sprintf(str, "CubeZ %s", CZ_VERSION);

  string HostName = GetHostName();

  // 結果出力(排他測定のみ)
  PM.print(stdout, HostName, str);
  PM.print(fp, HostName, str);
  PM.printDetail(fp);

  for (int i=0; i<numProc; i++)
  {
    PM.printThreads(fp, i, 0); // time cost order
  }

	PM.printLegend(fp);
  Hostonly_ {
    fflush(fp);
    fclose(fp);
  }

#endif


  char tmp_fname[30];
  int loc[3];

  if (debug_mode==1) {

    double errmax = 0.0;

    sprintf( tmp_fname, "p_%05d.sph", myRank );
    fileout_t_(size, &gc, P, pitch, origin, tmp_fname);
    exact_t_(size, &gc, ERR, pitch, origin);
    err_t_  (size, innerFidx, &gc, &errmax, P, ERR, loc);
    if ( !Comm_MAX_1(&errmax, "Comm_Res_Poisson") ) return 0;
    Hostonly_ printf("\nError max = %e at (%d %d %d)\n\n", errmax, loc[0],loc[1],loc[2]);
    sprintf( tmp_fname, "e_%05d.sph", myRank );
    fileout_t_(size, &gc, ERR, pitch, origin, tmp_fname);

  } // debug


  return 1;
}



void CZ::setStrPre()
{
  if ( !strcasecmp(precon.c_str(), "jacobi") ) {
    pc_type = LS_JACOBI;
  }
  else if ( !strcasecmp(precon.c_str(), "psor") ) {
    pc_type = LS_PSOR;
  }
  else if ( !strcasecmp(precon.c_str(), "sor2sma") ) {
    pc_type = LS_SOR2SMA;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_rb") ) {
    pc_type = LS_PCR_RB;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_rb_esa") ) {
    pc_type = LS_PCR_RB_ESA;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr") ) {
    pc_type = LS_PCR;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_eda") ) {
    pc_type = LS_PCR_EDA;
  }
  else if ( !strcasecmp(precon.c_str(), "jacobi_maf") ) {
    pc_type = LS_JACOBI_MAF;
  }
  else if ( !strcasecmp(precon.c_str(), "psor_maf") ) {
    pc_type = LS_PSOR_MAF;
  }
  else if ( !strcasecmp(precon.c_str(), "sor2sma_maf") ) {
    pc_type = LS_SOR2SMA_MAF;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_rb_maf") ) {
    pc_type = LS_PCR_RB_MAF;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_rb_esa_maf") ) {
    pc_type = LS_PCR_RB_ESA_MAF;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_maf") ) {
    pc_type = LS_PCR_MAF;
  }
  else if ( !strcasecmp(precon.c_str(), "pcr_eda_maf") ) {
    pc_type = LS_PCR_EDA_MAF;
  }
  else printf("precon=%s\n", precon.c_str());
}



std::string CZ::printMethod(int type)
{
  std::string str;
  
  if ( type == LS_JACOBI ) {
    str = "jacobi";
  }
  else if ( type == LS_PSOR ) {
    str = "psor";
  }
  else if ( type == LS_SOR2SMA ) {
    str = "sor2sma";
  }
  else if ( type == LS_PCR ) {
    str = "pcr";
  }
  else if ( type == LS_PCR_EDA ) {
    str = "pcr_eda";
  }
  else if ( type == LS_PCR_ESA ) {
    str = "pcr_esa";
  }
  else if ( type == LS_PCR_RB ) {
    str = "pcr_rb";
  }
  else if ( type == LS_PCR_RB_ESA ) {
    str = "pcr_rb_esa";
  }
  else if ( type == LS_JACOBI_MAF ) {
    str = "jacobi_maf";
  }
  else if ( type == LS_PSOR_MAF ) {
    str = "psor_maf";
  }
  else if ( type == LS_SOR2SMA_MAF ) {
    str = "sor2sma_maf";
  }
  else if ( type == LS_PCR_RB_MAF ) {
    str = "pcr_rb_maf";
  }
  else if ( type == LS_PCR_RB_ESA_MAF ) {
    str = "pcr_rb_esa_maf";
  }
  else if ( type == LS_PCR_MAF ) {
    str = "pcr_maf";
  }
  else if ( type == LS_PCR_EDA_MAF ) {
    str = "pcr_eda_maf";
  }
  else if ( type == LS_PCR_ESA_MAF ) {
    str = "pcr_esa_maf";
  }
  
  return str;
}



void CZ::setLS(char* q, char* fname)
{
  if ( !strcasecmp(q, "jacobi") ) {
    ls_type = LS_JACOBI;
    strcpy(fname, "jacobi.txt");
  }
  
  else if ( !strcasecmp(q, "psor") ) {
    ls_type = LS_PSOR;
    strcpy(fname, "psor.txt");
  }
  
  else if ( !strcasecmp(q, "sor2sma") ) {
    ls_type = LS_SOR2SMA;
    strcpy(fname, "sor2sma.txt");
  }
  
  else if ( !strcasecmp(q, "pbicgstab") ) {
    ls_type = LS_BICGSTAB;
    strcpy(fname, "pbicgstab.txt");
    setStrPre();
  }
  
  else if ( !strcasecmp(q, "pcr_rb") ) {
    ls_type = LS_PCR_RB;
    strcpy(fname, "pcr_rb.txt");
  }
  
  else if ( !strcasecmp(q, "pcr_rb_esa") ) {
    ls_type = LS_PCR_RB_ESA;
    strcpy(fname, "pcr_rb_esa.txt");
    SW_esa = 1;
  }
  
  else if ( !strcasecmp(q, "pcr") ) {
    ls_type = LS_PCR;
    strcpy(fname, "pcr.txt");
  }
  
  else if ( !strcasecmp(q, "pcr_eda") ) {
    ls_type = LS_PCR_EDA;
    strcpy(fname, "pcr_eda.txt");
  }
  
  else if ( !strcasecmp(q, "pcr_esa") ) {
    ls_type = LS_PCR_ESA;
    strcpy(fname, "pcr_esa.txt");
    SW_esa = 1;
  }
  
  
  // MAF
  else if ( !strcasecmp(q, "jacobi_maf") ) {
    ls_type = LS_JACOBI_MAF;
    strcpy(fname, "jacobi_maf.txt");
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "psor_maf") ) {
    ls_type = LS_PSOR_MAF;
    strcpy(fname, "psor_maf.txt");
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "sor2sma_maf") ) {
    ls_type = LS_SOR2SMA_MAF;
    strcpy(fname, "sor2sma_maf.txt");
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "pbicgstab_maf") ) {
    ls_type = LS_BICGSTAB_MAF;
    strcpy(fname, "pbicgstab_maf.txt");
    setStrPre();
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "pcr_rb_maf") ) {
    ls_type = LS_PCR_RB_MAF;
    strcpy(fname, "pcr_rb_maf.txt");
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "pcr_rb_esa_maf") ) {
    ls_type = LS_PCR_RB_ESA_MAF;
    strcpy(fname, "pcr_rb_esa_maf.txt");
    SW_maf = 1;
    SW_esa = 1;
  }
  
  else if ( !strcasecmp(q, "pcr_maf") ) {
    ls_type = LS_PCR_MAF;
    strcpy(fname, "pcr_maf.txt");
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "pcr_eda_maf") ) {
    ls_type = LS_PCR_EDA_MAF;
    strcpy(fname, "pcr_eda_maf.txt");
    SW_maf = 1;
  }
  
  else if ( !strcasecmp(q, "pcr_esa_maf") ) {
    ls_type = LS_PCR_ESA_MAF;
    strcpy(fname, "pcr_esa_maf.txt");
    SW_maf = 1;
    SW_esa = 1;
  }
  
  else {
    printf("Invalid solver\n");
    exit(0);
  }
}
