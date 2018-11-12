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
subroutine init_mask(x, sz, idx, g)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x

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
do k=1-g,kx+g
do j=1-g,jx+g
do i=1-g,ix+g
  x(i, j, k) = 0.0
end do
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k = kst, ked
do j = jst, jed
do i = ist, ied
  x(i, j, k) = 1.0
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine init_mask


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
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1-g,kx+g
do j=1-g,jx+g
do i=1-g,ix+g
  x(i, j, k) = 0.0
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
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  y, x

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1-g,kx+g
do j=1-g,jx+g
do i=1-g,ix+g
  y(i, j, k) = x(i, j, k)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_copy


!> ********************************************************************
!! @brief AXPYZ
!! @param [out]    z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine blas_triad(z, x, y, a, sz, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x, y, z
double precision                                       ::  flop, a

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix) * dble(jx) * dble(kx) * 2.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
  z(i, j, k) = a * x(i, j, k) + y(i, j, k)
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
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
subroutine blas_dot1(r, p, sz, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p
double precision                                       ::  flop, q, r

ix = sz(1)
jx = sz(2)
kx = sz(3)

r  = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:r) &
!$OMP PRIVATE(q)
do k=1,kx
do j=1,jx
do i=1,ix
  q = dble(p(i, j, k))
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
!! @param [in]  g    ガイドセル
!! @param [in,out] flop flop count
!<
subroutine blas_dot2(r, p, q, sz, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, q
double precision                                       ::  flop, r

ix = sz(1)
jx = sz(2)
kx = sz(3)

r  = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:r)
do k=1,kx
do j=1,jx
do i=1,ix
  r = r + dble(p(i, j, k)) * dble(q(i, j, k))
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
!! @param [in]     g    ガイドセル
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine blas_bicg_1(p, r, q, beta, omg, sz, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, r, q
double precision                                       ::  flop, beta, omg

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
  p(i,j,k) = r(i,j,k) + beta * ( p(i,j,k) - omg * q(i,j,k) )
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
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_bicg_2(z, x, y, a, b, sz, g, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x, y, z
double precision                                       ::  flop, a, b

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix) * dble(jx) * dble(kx) * 4.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
  z(i, j, k) = a * x(i, j, k) + b * y(i, j, k) + z(i, j, k)
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
!! @param [in]  g    ガイドセル
!! @param [in]  cf   係数
!! @param [in,out] flop flop count
!<
subroutine blas_calc_ax(ap, p, sz, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  ap, p
double precision                                       ::  flop
real, dimension(7)                                     ::  cf

ix = sz(1)
jx = sz(2)
kx = sz(3)

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + dble(ix)*dble(jx)*dble(kx)*13.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
do k=1,kx
do j=1,jx
do i=1,ix
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  ap(i, j, k) = (ss - dd * p(i, j, k))
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_ax


!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [out]    r  残差ベクトル
!! @param [in]     p  解ベクトル
!! @param [in]     b  定数項
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!! @param [in,out] flop flop count
!<
subroutine blas_calc_rk(r, p, b, sz, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  r, p, b
double precision                                       ::  flop
real, dimension(7)                                     ::  cf

ix = sz(1)
jx = sz(2)
kx = sz(3)

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + dble(ix)*dble(jx)*dble(kx)*14.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
do k=1,kx
do j=1,jx
do i=1,ix
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  r(i, j, k) = (b(i, j, k) - (ss - dd * p(i, j, k)))
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
!! @param [in]  g    ガイドセル長
!! @param [in,out] flop flop count
!<
subroutine blas_calc_r2 (res, p, b, sz, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
double precision                                       ::  flop, res
real                                                   ::  dd, ss, dp, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, b
real, dimension(7)                                     ::  cf

ix = sz(1)
jx = sz(2)
kx = sz(3)

res = 0.0
c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + dble(ix)*dble(jx)*dble(kx)*16.0d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ss, dp)
do k=1,kx
do j=1,jx
do i=1,ix
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  dp = ( b(i,j,k) - (ss - dd * p(i,j,k)) )
  res = res + dble(dp*dp)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_r2




!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [in,out]    r  残差ベクトル
!! @param [in]     p  解ベクトル
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!! @param [in,out] flop flop count
!! @note 定数項をrに入れて、上書きする
!<
subroutine blas_calc_rk_wob(r, p, sz, g, cf, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  r, p
double precision                                       ::  flop
real, dimension(7)                                     ::  cf
! !dir$ assume_aligned r:64, p:64

ix = sz(1)
jx = sz(2)
kx = sz(3)

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

flop = flop + dble(ix)*dble(jx)*dble(kx)*14.0d0

! !dir$ vector aligned
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
do k=1,kx
do j=1,jx
do i=1,ix
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  r(i, j, k) = (r(i, j, k) - (ss - dd * p(i, j, k)))
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_rk_wob


!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [in,out]    r  残差ベクトル
!! @param [in]     p  解ベクトル
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!! @param [in,out] flop flop count
!! @note 定数項をrに入れて、上書きする
!<
subroutine blas_calc_rk_1d(r, p, m, sz, g, cf, flop)
implicit none
integer                                                ::  i, ix, jx, kx, g, ist, ied, j1, k1
integer, dimension(3)                                  ::  sz
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension((sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))   ::  r, p, m
double precision                                       ::  flop
real, dimension(7)                                     ::  cf

ix = sz(1)
jx = sz(2)
kx = sz(3)
! ijk = (ix+2*g)*(jx+2*g)*(k+g-1)+(ix+2*g)*(j+g-1)+i+g

ist = (ix+2*g)*(jx+2*g)*g + (ix+2*g)*g + 1+g
ied = (ix+2*g)*(jx+2*g)*(kx+g-1)+(ix+2*g)*(jx+g-1)+ix+g

c1 = cf(1)
c2 = cf(2)
c3 = cf(3)
c4 = cf(4)
c5 = cf(5)
c6 = cf(6)
dd = cf(7)

j1 = (ix+2*g)*g
k1 = (ix+2*g)*(jx+2*g)*g

flop = flop + dble(ied-ist+1)*14.0d0

!$OMP PARALLEL DO SCHEDULE(static) PRIVATE(ss)
do i=ist,ied
  ss = c1 * p(i+1) &
     + c2 * p(i-1) &
     + c3 * p(i+j1) &
     + c4 * p(i-j1) &
     + c5 * p(i+k1) &
     + c6 * p(i-k1)
  r(i) = (r(i) - (ss - dd * p(i))*m(i) )
end do
!$OMP END PARALLEL DO

return
end subroutine blas_calc_rk_1d
