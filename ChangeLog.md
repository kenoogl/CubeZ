# CubeZ

## todo
BiCGSTABの境界条件は係数に織り込むFFVC方式でないと、導入できない
現状はディリクレ条件でのテストのみ


## REVISION HISTORY

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
  - examp[le/scripts


---
- 2018-11-13 Version 0.5.0
  - BiCGSTAB
  - 境界条件はディリクレで統一


---
- 2018-11-10 Version 0.1.0
  - Initial setup
