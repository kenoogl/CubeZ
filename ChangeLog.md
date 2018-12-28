# CubeZ

## todo
BiCGSTABの境界条件は係数に織り込むFFVC方式でないと、導入できない
現状はディリクレ条件でのテストのみ



## REVISION HISTORY

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
