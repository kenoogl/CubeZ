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

int CZ::Evaluate(int argc, char **argv)
{
  int div_type = 0;         ///< 分割指定 (0-自動、1-指定)
  double G_Memory = 0.0;    ///< 計算に必要なメモリ量（グローバル）
  double L_Memory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  int gc = GUIDE;

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
  char* q = argv[4];
  if ( !strcasecmp(q, "pbicgstab") )
  {
    if (argc!=7 && argc!=10) {
      Hostonly_ printf("command line error : pbicgstab\n");
      exit(0);
    }
    char* w = argv[6];
    setPrecond(w);
  }

  if (argc == 9) {
    div_type=1;
    G_div[0] = atoi(argv[6]);
    G_div[1] = atoi(argv[7]);
    G_div[2] = atoi(argv[8]);
  }

  if (argc == 10) {
    div_type=1;
    G_div[0] = atoi(argv[7]);
    G_div[1] = atoi(argv[8]);
    G_div[2] = atoi(argv[9]);
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

  else if ( !strcasecmp(q, "ljacobi") ) {
    ls_type = LS_LJACOBI;
    strcpy(fname, "ljacobi.txt");
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

  if( (RHS = Alloc_Real_S3D(size)) == NULL ) return 0;
  if( (P   = Alloc_Real_S3D(size)) == NULL ) return 0;
  if( (WRK = Alloc_Real_S3D(size)) == NULL ) return 0;
  //if( (MSK = Alloc_Real_S3D(size)) == NULL ) return 0;

  if (debug_mode == 1) {
    L_Memory += ( array_size * 1 ) * (double)sizeof(REAL_TYPE);

    //if( (EXS = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (ERR = Alloc_Real_S3D(size)) == NULL ) return 0;
  }


  if (ls_type == LS_BICGSTAB)
  {
    L_Memory += ( array_size * 9 ) * (double)sizeof(REAL_TYPE);

    if( (pcg_p  = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_p_ = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_r  = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_r0 = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_q  = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_s  = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_s_ = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_t  = Alloc_Real_S3D(size)) == NULL ) return 0;
    if( (pcg_t_ = Alloc_Real_S3D(size)) == NULL ) return 0;
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
  bc_(size, &gc, P, pitch, origin, nID);

  // source term
  //src_dirichlet_(RHS, size, &gc, pitch, nID);

  if ( !Comm_S(RHS, 1) ) return 0;



  // タイミング測定の初期化
  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( myRank );
  setParallelism();
  PM.setParallelMode(Parallel_str, numThreads, numProc);
  set_timing_label();



  /////////////////////////////////////////////////////////////
  // Loop

  double  res=0.0;
  int itr=0;
  double flop=0.0; // dummy

  TIMING_start("BoundaryCondition");
  bc_(size, &gc, P, pitch, origin, nID);
  TIMING_stop("BoundaryCondition");

  switch (ls_type)
  {
    case LS_JACOBI:
      TIMING_start("JACOBI");
      if ( 0 == (itr=JACOBI(res, P, RHS, ItrMax, flop)) ) return 0;
      TIMING_stop("JACOBI", flop);
      break;

    case LS_PSOR:
      TIMING_start("PSOR");
      if ( 0 == (itr=PSOR(res, P, RHS, ItrMax, flop)) ) return 0;
      TIMING_stop("PSOR", flop);
      break;

    case LS_SOR2SMA:
      TIMING_start("SOR2SMA");
      if ( 0 == (itr=RBSOR(res, P, RHS, ItrMax, flop)) ) return 0;
      TIMING_stop("SOR2SMA", flop);
      break;

    case LS_BICGSTAB:
      TIMING_start("PBiCGSTAB");
      if ( 0 == (itr=PBiCGSTAB(res, P, RHS, flop)) ) return 0;
      TIMING_stop("PBiCGSTAB", flop);
      break;

    case LS_LSOR:
      TIMING_start("LSOR");
      if ( 0 == (itr=LSOR(res, P, RHS, ItrMax, flop)) ) return 0;
      TIMING_stop("LSOR", flop);
      break;

    case LS_LJACOBI:
      TIMING_start("LJacobi");
      if ( 0 == (itr=LJacobi(res, P, RHS, ItrMax, flop)) ) return 0;
      TIMING_stop("LJacobi", flop);
      break;

    default:
      break;
  }

  Hostonly_ {
    printf("\n=================================\n");
    printf("Iter = %d  Res = %e\n", itr, res);
    printf("=================================\n");
  }


  /////////////////////////////////////////////////////////////
  // post

  if (!fph) fclose(fph);

  char tmp_fname[30];

  if (debug_mode==1) {
    sprintf( tmp_fname, "p_%05d.sph", myRank );
    fileout_(size, &gc, P, pitch, origin, tmp_fname);

    exact_(size, &gc, ERR, pitch, origin);
    double errmax = 0.0;
    err_  (size, innerFidx, &gc, &errmax, P, ERR);
    if ( !Comm_MAX_1(&errmax, "Comm_Res_Poisson") ) return 0;
    Hostonly_ printf("Error max = %e\n", errmax);
    sprintf( tmp_fname, "e_%05d.sph", myRank );
    fileout_(size, &gc, ERR, pitch, origin, tmp_fname);
  }


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


  return 1;
}
