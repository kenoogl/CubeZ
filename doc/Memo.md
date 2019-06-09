# Memo for CubeZ


#### 2019-05-22 ver 0.9.1

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

