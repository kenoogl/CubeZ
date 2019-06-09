!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################

!> @file   cz_blas.f90
!! @brief  BLAS routine
!! @author RIIT
!! @note based on Onishi version ffv_blas
!<


!> ********************************************************************
!! @brief 要素のマスク
!! @param [in,out] x  スカラー
!! @param [in]     sz 配列長
!! @param [in]     idx         インデクス範囲
!! @param [in]     g  ガイドセル
!<
subroutine imask_k(x, sz, idx, g)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x
!dir$ assume_aligned x:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j=1-g,jx+g
do i=1-g,ix+g
!dir$ vector aligned
!dir$ simd
do k=1-g,kx+g
  x(k, i, j) = 0.0
end do
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  x(k, i, j) = 1.0
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine imask_k


!> ********************************************************************
!! @brief 要素のゼロクリア
!! @param [in,out] x  スカラー
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
subroutine blas_clear(x, sz, g)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x
!dir$ assume_aligned x:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j=1-g,jx+g
do i=1-g,ix+g
  !dir$ vector aligned
  !dir$ simd
do k=1-g,kx+g
  x(k,i,j) = 0.0
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_clear


!> ********************************************************************
!! @brief コピー
!! @param [out]    y  コピー先
!! @param [in]     x  ソース
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
subroutine blas_copy(y, x, sz, g)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  y, x
!dir$ assume_aligned x:64, y:64
ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j=1-g,jx+g
do i=1-g,ix+g
!dir$ vector aligned
!dir$ simd
do k=1-g,kx+g
  y(k,i,j) = x(k,i,j)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_copy


!> ********************************************************************
!! @brief コピー
!! @param [out]    y  コピー先
!! @param [in]     x  ソース
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
subroutine blas_copy_in(y, x, sz, g)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  y, x
!dir$ assume_aligned x:64, y:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
!dir$ vector aligned
!dir$ simd
do k=1,kx
y(k,i,j) = x(k,i,j)
end do
end do
end do
!$OMP END PARALLEL DO


return
end subroutine blas_copy_in

!> ********************************************************************
!! @brief AXPYZ
!! @param [out]    z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine blas_triad(z, x, y, a, sz, idx, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, y, z
double precision                                       ::  flop, a
!dir$ assume_aligned x:64, y:64, z:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 2.0d0    &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  z(k,i,j) = a * x(k,i,j) + y(k,i,j)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_triad


!> ********************************************************************
!! @brief DOT1
!! @param [out] r    内積
!! @param [in]  p    ベクトル
!! @param [in]  sz   配列長
!! @param [in]  idx  インデクス範囲
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
subroutine blas_dot1(r, p, sz, idx, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p
double precision                                       ::  flop, q, r
!dir$ assume_aligned p:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r  = 0.0

flop = flop + 2.0d0    &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:r) &
!$OMP PRIVATE(q)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  q = dble(p(k,i,j))
  r = r + q*q
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_dot1


!> ********************************************************************
!! @brief DOT2
!! @param [out] r    内積
!! @param [in]  p    ベクトル
!! @param [in]  q    ベクトル
!! @param [in]  sz   配列長
!! @param [in]  idx  インデクス範囲
!! @param [in]  g    ガイドセル
!! @param [in,out] flop flop count
!<
subroutine blas_dot2(r, p, q, sz, idx, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, q
double precision                                       ::  flop, r
!dir$ assume_aligned p:64, q:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r  = 0.0

flop = flop + 2.0d0    &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:r)
do j = jst, jed
do i = ist, ied
  !dir$ vector aligned
  !dir$ simd
do k = kst, ked
  r = r + dble(p(k,i,j)) * dble(q(k,i,j))
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_dot2


!> ********************************************************************
!! @brief BiCGstabの部分演算1
!! @param [in,out] p    ベクトル
!! @param [in]     r    ベクトル
!! @param [in]     q    ベクトル
!! @param [in]     beta 係数
!! @param [in]     omg  係数
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine blas_bicg_1(p, r, q, beta, omg, sz, idx, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, r, q
double precision                                       ::  flop, beta, omg
!dir$ assume_aligned p:64, q:64, r:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 4.0d0    &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  p(k,i,j) = r(k,i,j) + beta * ( p(k,i,j) - omg * q(k,i,j) )
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_bicg_1


!> ********************************************************************
!! @brief BiCGstab 2
!! @param [in,out] z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     b    係数
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_bicg_2(z, x, y, a, b, sz, idx, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, y, z
double precision                                       ::  flop, a, b
!dir$ assume_aligned x:64, y:64, z:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 4.0d0    &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  z(k,i,j) = a * x(k,i,j) + b * y(k,i,j) + z(k,i,j)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_bicg_2


!> ********************************************************************
!! @brief AX
!! @param [out] ap   AX
!! @param [in]  p    解ベクトル
!! @param [in]  sz   配列長
!! @param [in]  idx  インデクス範囲
!! @param [in]  g    ガイドセル
!! @param [in]  cf   係数
!! @param [in,out] flop flop count
!<
subroutine blas_calc_ax(ap, p, sz, idx, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  ap, p
double precision                                       ::  flop
real, dimension(7)                                     ::  cf
!dir$ assume_aligned ap:64, p:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + 13.0d0   &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  ss = c1 * p(k  , i+1,j  ) &
     + c2 * p(k  , i-1,j  ) &
     + c3 * p(k  , i  ,j+1) &
     + c4 * p(k  , i  ,j-1) &
     + c5 * p(k+1, i  ,j  ) &
     + c6 * p(k-1, i  ,j  )
  ap(k, i, j) = (ss - dd * p(k, i, j))
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_ax


!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [out]    r    残差ベクトル
!! @param [in]     p    解ベクトル
!! @param [in]     b    定数項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル
!! @param [in]     cf   係数
!! @param [in,out] flop flop count
!<
subroutine blas_calc_rk(r, p, b, sz, idx, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  r, p, b
double precision                                       ::  flop
real, dimension(7)                                     ::  cf
!dir$ assume_aligned r:64, p:64, b:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + 14.0d0   &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  ss = c1 * p(k  , i+1,j  ) &
     + c2 * p(k  , i-1,j  ) &
     + c3 * p(k  , i  ,j+1) &
     + c4 * p(k  , i  ,j-1) &
     + c5 * p(k+1, i  ,j  ) &
     + c6 * p(k-1, i  ,j  )
  r(k, i, j) = (b(k, i, j) - (ss - dd * p(k, i, j)))
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_rk


!> ********************************************************************
!! @brief 残差の自乗和のみ
!! @param [out] res  残差の自乗和
!! @param [in]  p    圧力
!! @param [in]  b    RHS vector
!! @param [in]  sz   配列長
!! @param [in]  idx  インデクス範囲
!! @param [in]  g    ガイドセル長
!! @param [in,out] flop flop count
!<
subroutine blas_calc_r2 (res, p, b, sz, idx, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  dd, ss, dp, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
real, dimension(7)                                     ::  cf
!dir$ assume_aligned p:64, b:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res = 0.0

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + 16.0d0   &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ss, dp)
do j = jst, jed
do i = ist, ied
!dir$ vector aligned
!dir$ simd
do k = kst, ked
  ss = c1 * p(k  , i+1,j  ) &
     + c2 * p(k  , i-1,j  ) &
     + c3 * p(k  , i  ,j+1) &
     + c4 * p(k  , i  ,j-1) &
     + c5 * p(k+1, i  ,j  ) &
     + c6 * p(k-1, i  ,j  )
  dp = ( b(k,i,j) - (ss - dd * p(k,i,j)) )
  res = res + dble(dp*dp)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_r2
