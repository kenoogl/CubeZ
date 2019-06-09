# Memo for CubeZ



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

