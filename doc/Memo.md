# Memo for CubeZ

#### Version 1.2.1

#####  sor2smaの実装 >> インデクスでの実装でベクトル化する
1) ベクトル化されている
~~~
Y------>  do j=jst,jed
|+----->  do i=ist,ied
||V---->  do k=kst+mod(i+j+kp,2), ked, 2
~~~

2) ベクトル化できていない
~~~
Y------>  do j=jst,jed
|+----->  do i=ist,ied
||+---->  do k=kst,ked
|||       if(mod(i+j+k,2) /= color) cycle
~~~

#####  pcr_rbの実装　Auroraではif文の実装が性能高い

1) インデクスループで計算 >> 部分並列

~~~
Y------>  do j=jst, jed
|+----->   !do i=ist+mod(j+ip,2), ied, 2
~~~

2) if文 >> 並列化されている

~~~
P------>  do j=jst, jed
|+----->  do i=ist, ied
||        if(mod(i+j,2) /= color) cycle
~~~

- pcr_rb系の実装はcycleｂ文で実装


#### Version 1.2.0
- mretain-none復活
- SVRの初期化はFortran内で

####  Version 1.1.9
- CMakeLists.txt のACCからnofmaを削除
- nVidiaTools用のマクロ挿入 PUSH_RANGE, POP_RANGE
- メモリアロケート時点以降の全てにOpenACCディレクティブを入れる
- nvtx.f90
- cz_lsor.f90を外す
- openaccでOMPと同じくfirstprivateに変更

####  Version 1.1.8
- Autrora maf用のコンパイラオプション -static-nec

####  Version 1.1.7
- Autrora maf用のVectorReductionの配列の初期化コスト >> "VRtmp_Init"
- 対象コードをmafのみに修正（バグ）


####  Version 1.1.6
- pcr_eda, pcr_eda_maf のprivate(kl, kr) を外す
- pcr_esa / pcr_esa_maf, pcr_rb_esa/pcr_rb_esa_mafの確保する配列の大きさは1スレッド分として、サブルーチン内でprivate

####  Version 1.1.5
- 全てのmaf サブルーチンにvector reductionを実装
- OpenMPとOpenACCの排他制御、OpenACCが優先
- サブルーチンの整理
  - pcr + rb + [eda, esa]
  - pcrv > pcr_eda (extended dynamic array)
  - pcrv_sa > pcr_esa (extended static array)
- LS_PCR_RB_ESA, LS_PCR_RB_ESA_MAF追加
- pcr_rb_esa, pce_rb_esa_maf
- esa
  - firstprivate(id) > private(id)
  - private(a, c, d) ?
  - 変数sを削除 >> ss
  - pcr_rbでres　を　res1でreductionするとエラーになる件、原因不明

####  Version 1.1.4
- Auroraのコンパイルオプションの"-mretain-none"は不要
- printMethod(), setStrPre() 追加
- enable_VectorReduction >> jacobi_maf()のベクトルリダクションの制御 >> 　Aurorade
maf系のみ問題になっている
- czAllocR_S3D(), czAllocR(), czAllocR2()へのOpenACCディレクティブ追加
- thread_max >> numThreadsへ
- 1.1.3のjacobi_mafのaurora用の実装では、ベクトルリダクションを抑制するためにtmp配列をFortran側でアロケートしていたが、これを予め確保し、ゼロ初期化して使う実装に変更  553 > 742 GFLOPS @8 threads


####  Version 1.1.3
- NECからのフィードバック
- search_pivotでssはreduction対象ではないのでprivate(ss)に変更｀
- Auroraのコンパイルオプションに"-mretain-none"を追加
- src/CMakeLists.txtでAurora用のFortranコンパイラでリンクする指定
- bc_k()でスレッド間同期のオーバーヘッドを削減するため、各境界に関する処理を単一のparallel region内で行い、
更にスレッド間同期が不要な箇所にNOWAITを挿入．
- Auroraではjacobi_maf_()などにあるベクトル総和処理 "res1 = res1 + dp * dp" が非常に高コストとなっているので，k方向のサイズの配列tmpを定義し、主要ループ内ではこの配列への足しこみを行う実装にする


####  Version 1.1.1
- pcrの配列a,c, dはprivateにすると、フォーク時にアロケーションする？　スレッド分アロケートしておいて、アドレスを渡す

####  Version 1.1.0
- pcr系のflop count 修正
- pcrv_maf


####  Version 1.0.9
- Aurora用のmake >> Fortranでリンク
- cz_solver.f90のpcr()でインデクスのmin, maxを外すように配列領域を拡大 >> pcrv()
- pcrの最終段の計算は簡単になる. 最初の2x2だけでよい
- NECのFBからpcrvを実装 a,c,dはallocateを使う実装
- 反復回数、履歴ともpcrと同じになる。mac gccでの時間は 163.3  >> 54.3 secとなり 1/3


#### ver 1.0.8
- PGI compilerで最適化されない件（IntelはOKだが）
  - Loop not vectorized: mixed data typesのメッセージ
  - ループ内で単精度と倍精度を混在させている場合に最適化できないので、残差のreductionを単精度に変更

- 係数を単精度に変更
  - blas_triad, blas_dot1, blas_dot2, blas_bicg_1, blas_bicg_1, Fdot1, Fdot2
  
  - floatを追加: Comm_SUM_1, Comm_SUM_2, Comm_MIN_1, Comm_MAX_1
  
- メンテ中のサブルーチンは
  - cz_solver.f90, cz_blas.f90, cz_maf.f90


#### ver 1.0.7
- openACCでのビルド問題解決（とりあえず）
- cz_f90とsrc直下のCMakeLists.txtで
  ~~~
  set_target_properties(${cz_target} PROPERTIES LINKER_LANGUAGE Fortran)
  ~~~
  
- CMakeLists.txt
  
  ~~~
  if(with_ACC)
  if (with_ACC STREQUAL "Pascal")
  SET(ACC "-ta=tesla,cc60,cuda10.1,nofma,managed")
  elseif (with_ACC STREQUAL "Volta")
  SET(ACC "-ta=tesla,cc70,cuda10.1,nofma,managed")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc ${ACC}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc ${ACC}")
  SET(CMAKE_EXE_LINKER_FLAGS "-Mnomain")
  endif()
  ~~~

- `-tp=...managed`として，f90のソースはopenACC任せにする


#### ver 1.0.4
- AVXのオプション処理を変更 >> CompileOptionSelector.cmake
- gcc error 対応 >> collapse(2)節をいれない

~~~
cz_solver.f90:414:28:

do i=ist+mod(j+ip,2), ied, 2
1
Error: !$OMP DO collapsed loops don't form rectangular iteration space at (1)
~~~

#### ver 1.0.3
- OpenACCのオプションを追加 with_ACC={off|Pascal|Volta}
- with_AVX512 >> with_SIMD={OFF|256|512}
- INTEL_SKLを削除


#### ver 1.0.1
- ftraceオプションを追加（Auroraの場合のみ）`-fpp` を指定すると`_aurora_(=1)`が定義されることを利用
- IVDEPディレクティブを追加し、効果あり
- shortloop, -mretain-noneは効果なし


#### ver 1.0.0
- aurora でコンパイル
- ファイル出力は マクロ`_aurora_` を使ってコンパイルしないように制御


#### ver 0.9.9
- aurora ブランチを作成、auroraのビルド対応
- PMlibを外すときの処理　> `DISABLE_PMLIB`
- pcr / pcr_rbの変数 `a,c,d,a1,c1,d1`を渡すように変更 >> `NEC2003f_Alloc_var`エラーがでるため
- ALIGNMENT, SIMD_WIDTH削除



#### ver 0.9.8
- PCR, PCR_RBにフォーカスした構成

#### 2019-06-16 ver 0.9.6
- ピボット計算の !OMPreduction(max:ss) を追加
- `lsor_p7_maf`で係数の計算を修正。lsor_p7と一致を確認


#### 2019-06-10 ver 0.9.5
- pbicgstabとpbicgstab_mafでmafの方が収束が遅いので調査。前者は係数が1，後者は10^4のオーダー。スケーリングを導入し、`calc_ax_maf`, `calc_rk_maf`でテスト。同程度になることを確認。
- psorなどはrpをddで除くとpvt()が相殺して1になるので不要


#### 2019-06-09 ver 0.9.4
- cz_maf.f90 pbicgstab_mafの実装


#### 2019-06-09 ver 0.9.3
- cz_lsor.f90  >> 最終的にはテストしないソース
- cz_solver.f90 >> 比較対象コード　jacobi, sor2sma, pbicgstab, lsor_p7
- cz_maf.f90 MAFバージョン (jacobi, psor, sor2sma)


#### 2019-06-09 ver 0.9.2
-  lsor_p7を前処理に追加
- 比較ターゲットはjacobi, sor2sma, pbicgstab, lsor_p7



#### 2019-06-09 ver 0.9.1
-  SIMD_AVX512へ変更

- 各メソッドにコメント追加(LSOR_P1~P7)
  - LSOR_P1     PCRのベースライン、ワークは全て3D配列
  - LSOR_P2     LSOR_P1から、最終段を直接反転
  - LSOR_P3     LSOR_P2から、最終段の直接反転を分割
  - LSOR_P4     LSOR_P3から、最終段の直接反転を関数か呼び出しから直接展開へ
  - LSOR_P5     LSOR_P4から、一次元配列利用、テンポラリに倍精度
  - LSOR_P6     LSOR_P5から、テンポラリを単精度にして変換コストを下げる
  - LSOR_P7     LSOR_P6から、マルチカラー化
  
  - LSORはSOR系に比べて、加速係数を大きくとれるので、反復回数が少なくて済む


#### 2019-05-22 ver 0.9.0

- Alignment=64に固定。
- SIMD=256デフォルト
- 配列は KIJ のみにする。
- 不要コードを削除。

	~~~
	LS_LSOR_A
	LS_LSOR_B
	LS_LSOR_C
	LS_LSOR_D
	LS_LSOR_E
	LS_LSOR_F
	LS_LSOR_J
    LS_LSOR_J4
    LS_LSOR_K
    LS_LSOR_K2
    LS_LSOR_K3
    LS_LJCB_A
    LS_LJCB_B
    LS_LJCB_C
    LS_LJCB_D
    LS_LJCB_E
    LS_LJCB_F
    LS_LSOR_SIMD
	~~~


#### 2019-05-21 ver 0.9.0

- 逐次計算の実装を確認している。
- 逐次計算のコンパイルのために、CBrickからヘッダをコピーしている。
- コードは `main >> cz_Evaluate >> cz_Poisson`
- Alignment情報は、size[0] をみている。これは組み込みSIMD関数のためのトライアルで、現時点の実装ではなくてもよい。
- 選択スキームによって、IJK, KIJ, IKJの3種類ある。

