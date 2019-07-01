# CubeZ

## todo
BiCGSTABの境界条件は係数に織り込むFFVC方式でないと、導入できない
現状はディリクレ条件でのテストのみ



## REVISION HISTORY

---
- 2019-7-1 Version 1.0.7
- add -ta=,,,managed


---
- 2019-7-1 Version 1.0.6
- add _OPENACC macro
- add "-cpp" option for gfortran


---
- 2019-7-1 Version 1.0.5
- modify OpenACC directive and CMakefile.txt


---
- 2019-7-1 Version 1.0.4
  - OpenACC directive
  - modify AVX/AVX512 options


---
- 2019-6-27 Version 1.0.3
  - OpenACC option
  - change SIMD option
  - rm Toolchanin INTEL_SKL

---
- 2019-6-24 Version 1.0.2
  - remove ratain option


---
- 2019-6-23 Version 1.0.1
  - add  trace option
  - IVDEP

---
- 2019-6-23 Version 1.0.0
  - Aurora compiled

---
- 2019-6-23 Version 0.9.9
  - branch aurora

---
- 2019-6-20 Version 0.9.8
  - `lsor_p7` >> `pcr`
  - `lsor_pcr_kij7` >> `pcr`

---
- 2019-6-17 Version 0.9.7
  - confirm OMP behavior
  - change ITO A/B configure script in Readme.md

---
- 2019-6-16 Version 0.9.6
  - fix `lsor_p7_maf`

---
- 2019-6-10 Version 0.9.5
  - scaling for `calc_ax_maf`, `calc_rk_maf`

---
- 2019-6-9 Version 0.9.4
  - pbicgstab_maf


---
- 2019-6-9 Version 0.9.3
- add cz_maf.f90


---
- 2019-6-9 Version 0.9.2
  - add lsor_p7 to preconditioner

---
- 2019-6-9 Version 0.9.1
  - SIMD_width => SIMD_AVX512
  - Halt test6 and 7 due to SIMD intrinsic
  - lsor_pcr_kij7 : remove collapse(2) due to loop  i= , , 2
  - remove extra source
    - jacobi, psor, sor2sma, pbicgstab, lsor_p{1-7}


---
- 2019-4-8 Version 0.9.0
  - lsor_P5, P6, P7

---
- 2019-3-22 Version 0.8.9
  - lsor_pcr_kij4, 5

---
- 2019-3-21 Version 0.8.8
  - lsor_pcr_kij  omp directive modify


---
- 2019-2-18 Version 0.8.7
   - PCR
   - lsor_pcr_kij2

---
- 2018-12-31 Version 0.8.5
   - ベストケースはlsor_k2()
   - スレッドスケジュールはstatic

---
- 2018-12-31 Version 0.8.4
  - 最適化のオプション変更
  - Alignment, SIMD_width

---
- 2018-12-28 Version 0.8.3
  - テストしない関数をobsoleteへ移動

---
- 2018-12-28 Version 0.8.2
  - LSOR_J Best case

---
- 2018-12-25 Version 0.8.1
  - LSOR_J

---
- 2018-12-18 Version 0.8.0
  - tdma6_8() unroll(8)

---
- 2018-12-18 Version 0.7.9
  - tdma7()でIn-palceをやめ、d2ベクトルで出力、スワップ

---
- 2018-12-18 Version 0.7.8
  - tdma7()で受け渡し方法の別実装

---
- 2018-12-18 Version 0.7.7
  - lsor_simd4() 2-color オーダリング

---
- 2018-12-18 Version 0.7.6
  - tdma6()

---
- 2018-12-18 Version 0.7.5
  - relax_256()
  - (i,j)サンプリング > lsor_simd3()

---
- 2018-12-18 Version 0.7.4
  - relax()

---
- 2018-12-17 Version 0.7.3
  - tdma4()

---
- 2018-12-17 Version 0.7.2
  - lsor_simd3() remainder loop の削除

---
- 2018-12-17 Version 0.7.1
  - lsor_simd2() unroll(4)の実装

---
- 2018-12-16 Version 0.7.0
  - lsor_simd() peel/body/remainderのループ長調整

---
- 2018-12-16 Version 0.6.9
  - アラインしたメモリ確保

---
- 2018-12-16 Version 0.6.8
  - LSOR_SIMD()のベース実装

---
- 2018-12-16 Version 0.6.7
  - lsor_simdにk方向内側ループ実装
  - ENABLE_AVX2/AVX512
  - ケースの整理

---
- 2018-12-14 Version 0.6.6
  - MSG テスト

---
- 2018-12-14 Version 0.6.5
  - PBiCGSTAB前処理にLSOR/LJCBを実装

---
- 2018-12-13 Version 0.6.4
  - SIMD test

---
- 2018-12-10 Version 0.6.3
  - Serial 対応
  - Ljcb分割 MSF

---
- 2018-12-09 Version 0.6.2
  - SKL オプション

---
- 2018-11-23 Version 0.6.1
  - LSORのk-loopをTDMA

---
- 2018-11-23 Version 0.6.0
  - PCR

---
- 2018-11-13 Version 0.5.3
  - LSOR
  - test2 マスク付き、D, N境界条件

---
- 2018-11-13 Version 0.5.2
  - PBiCGSTAB bug fix


---
- 2018-11-13 Version 0.5.1
  - example/scripts


---
- 2018-11-13 Version 0.5.0
  - BiCGSTAB
  - 境界条件はディリクレで統一


---
- 2018-11-10 Version 0.1.0
  - Initial setup
