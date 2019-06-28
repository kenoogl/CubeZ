# Memo for CubeZ

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

