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

! スレッド同期のオーバーヘッド抑制のため，単一のparallel regionとする
#ifndef _OPENACC
!$OMP PARALLEL
#endif


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
do j=1-g,jx+g
do i=1-g,ix+g
do k=1-g,kx+g
  x(k, i, j) = 0.0
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
#endif


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  x(k, i, j) = 1.0
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
#endif

#ifndef _OPENACC
!$OMP END PARALLEL
#endif

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


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
#endif
do j=1-g,jx+g
do i=1-g,ix+g
do k=1-g,kx+g
  x(k,i,j) = 0.0
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


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


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
#endif
do j=1-g,jx+g
do i=1-g,ix+g
do k=1-g,kx+g
  y(k,i,j) = x(k,i,j)
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


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


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
#endif
do j=1,jx
do i=1,ix
do k=1,kx
y(k,i,j) = x(k,i,j)
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


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
double precision                                       ::  flop
real                                                   ::  a
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


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  z(k,i,j) = a * x(k,i,j) + y(k,i,j)
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif

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
integer                                                ::  i, j, k, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p
double precision                                       ::  flop
real                                                   ::  q, r
!dir$ assume_aligned p:64

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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3) private(q) reduction(+:r)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:r) PRIVATE(q)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  q = p(k,i,j)
  r = r + q*q
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif

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
integer                                                ::  i, j, k, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, q
double precision                                       ::  flop
real                                                   ::  r
!dir$ assume_aligned p:64, q:64

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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3) reduction(+:r)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) REDUCTION(+:r)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  r = r + p(k,i,j) * q(k,i,j)
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif

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
integer                                                ::  i, j, k, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, r, q
double precision                                       ::  flop
real                                                   ::  beta, omg
!dir$ assume_aligned p:64, q:64, r:64

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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  p(k,i,j) = r(k,i,j) + beta * ( p(k,i,j) - omg * q(k,i,j) )
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


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
integer                                                ::  i, j, k, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, y, z
double precision                                       ::  flop
real                                                   ::  a, b
!dir$ assume_aligned x:64, y:64, z:64

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

#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  z(k,i,j) = a * x(k,i,j) + b * y(k,i,j) + z(k,i,j)
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


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
integer                                                ::  i, j, k, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  ap, p
double precision                                       ::  flop
real, dimension(7)                                     ::  cf
!dir$ assume_aligned ap:64, p:64

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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3) private(ss)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
#endif
do j = jst, jed
do i = ist, ied
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
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


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
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real                                                   ::  dd, ss, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  r, p, b
double precision                                       ::  flop
real, dimension(7)                                     ::  cf
!dir$ assume_aligned r:64, p:64, b:64

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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3) private(ss)
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(ss)
#endif
do j = jst, jed
do i = ist, ied
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
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


return
end subroutine blas_calc_rk


!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [out]    r    残差ベクトル
!! @param [in]     p    解ベクトル
!! @param [in]     b    定数項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル
!! @param [in]     X,Y,Z  座標
!! @param [in]     pvt  行の最大係数
!! @param [in,out] flop flop count
!<
subroutine calc_rk_maf(r, p, b, sz, idx, g, X, Y, Z, pvt, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real                                                   ::  GX, EY, TZ, YJA, YJAI
real                                                   ::  XG, YE, ZT, XGG, YEE, ZTT
real                                                   ::  C1, C2, C3, C7, C8, C9
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  r, p, b, pvt
double precision                                       ::  flop
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z
!dir$ assume_aligned r:64, p:64, b:64, X:64, Y:64, Z:64, pvt:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 63.0d0   &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3)
#else
!$OMP PARALLEL DO Collapse(2) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked

XG = 0.5 * (X(i+1) - X(i-1))
YE = 0.5 * (Y(j+1) - Y(j-1))
ZT = 0.5 * (Z(k+1) - Z(k-1)) ! 6

XGG= X(i+1) - 2.0*X(i) + X(i-1)
YEE= Y(j+1) - 2.0*Y(j) + Y(j-1)
ZTT= Z(k+1) - 2.0*Z(k) + Z(k-1) ! 9

! Jacobian
YJA  = XG * YE * ZT
YJAI = 1.0 / YJA    ! 3

! 1st order
GX =  YE * ZT * YJAI
EY =  XG * ZT * YJAI
TZ =  XG * YE * YJAI ! 6

! Laplacian
C1 =  GX * GX
C2 =  EY * EY
C3 =  TZ * TZ
C7 = -XGG * C1 * GX
C8 = -YEE * C2 * EY
C9 = -ZTT * C3 * TZ ! 9

r(k,i,j) = ( b(k,i,j)                           &
         + 2.0 * (C1 + C2 + C3) * P(k,i,j)      &
         - (C1 + 0.5 * C7) * P(k  , i+1, j  )   &
         - (C1 - 0.5 * C7) * P(k  , i-1, j  )   &
         - (C2 + 0.5 * C8) * P(k  , i  , j+1)   &
         - (C2 - 0.5 * C8) * P(k  , i  , j-1)   &
         - (C3 + 0.5 * C9) * P(k+1, i  , j  )   &
         - (C3 - 0.5 * C9) * P(k-1, i  , j  ) ) &
         * pvt(k,i,j) ! 30
enddo
enddo
enddo
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif



return
end subroutine calc_rk_maf


!> ********************************************************************
!! @brief AX
!! @param [out] ap   AX
!! @param [in]  p    解ベクトル
!! @param [in]  sz   配列長
!! @param [in]  idx  インデクス範囲
!! @param [in]  g    ガイドセル
!! @param [in]     X,Y,Z  座標
!! @param [in]     pvt  行の最大係数
!! @param [in,out] flop flop count
!<
subroutine calc_ax_maf(ap, p, sz, idx, g, X, Y, Z, pvt, flop)
implicit none
integer                                                ::  i, j, k, g
integer, dimension(3)                                  ::  sz
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(0:5)                                ::  idx
real                                                   ::  GX, EY, TZ, YJA, YJAI
real                                                   ::  XG, YE, ZT, XGG, YEE, ZTT
real                                                   ::  C1, C2, C3, C7, C8, C9
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  ap, p, pvt
double precision                                       ::  flop
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 63.0d0   &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3)
#else
!$OMP PARALLEL DO Collapse(2) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked

XG = 0.5 * (X(i+1) - X(i-1))
YE = 0.5 * (Y(j+1) - Y(j-1))
ZT = 0.5 * (Z(k+1) - Z(k-1)) ! 6

XGG= X(i+1) - 2.0*X(i) + X(i-1)
YEE= Y(j+1) - 2.0*Y(j) + Y(j-1)
ZTT= Z(k+1) - 2.0*Z(k) + Z(k-1) ! 9

! Jacobian
YJA  = XG * YE * ZT
YJAI = 1.0 / YJA    ! 3

! 1st order
GX =  YE * ZT * YJAI
EY =  XG * ZT * YJAI
TZ =  XG * YE * YJAI ! 6

! Laplacian
C1 =  GX * GX
C2 =  EY * EY
C3 =  TZ * TZ
C7 = -XGG * C1 * GX
C8 = -YEE * C2 * EY
C9 = -ZTT * C3 * TZ ! 9

ap(k,i,j) = (                                  &
          + (C1 + 0.5 * C7) * P(k  , i+1, j  ) &
          + (C1 - 0.5 * C7) * P(k  , i-1, j  ) &
          + (C2 + 0.5 * C8) * P(k  , i  , j+1) &
          + (C2 - 0.5 * C8) * P(k  , i  , j-1) &
          + (C3 + 0.5 * C9) * P(k+1, i  , j  ) &
          + (C3 - 0.5 * C9) * P(k-1, i  , j  ) &
          - 2.0 * (C1 + C2 + C3) * P(k,i,j)  ) &
          * pvt(k,i,j) ! 30
enddo
enddo
enddo
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


return
end subroutine calc_ax_maf


!> ********************************************************************
!! @brief 最大要素探索
!! @param [out]    pvt  行の最大要素
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル
!! @param [in]     X,Y,Z  座標
!<
subroutine search_pivot(pvt, sz, idx, g, X, Y, Z)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
real                                                   ::  GX, EY, TZ, YJA, YJAI
real                                                   ::  XG, YE, ZT, XGG, YEE, ZTT
real                                                   ::  C1, C2, C3, C7, C8, C9
real                                                   ::  s1, s2, s3, s4, s5, s6, s7, ss
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  pvt
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z
!dir$ assume_aligned pvt:64, X:64, Y:64, Z:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3)
#else
!$OMP PARALLEL DO Collapse(2) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9) &
!$OMP PRIVATE(s1, s2, s3, s4, s5, s6, s7, ss)
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked

XG = 0.5 * (X(i+1) - X(i-1))
YE = 0.5 * (Y(j+1) - Y(j-1))
ZT = 0.5 * (Z(k+1) - Z(k-1)) ! 6

XGG= X(i+1) - 2.0*X(i) + X(i-1)
YEE= Y(j+1) - 2.0*Y(j) + Y(j-1)
ZTT= Z(k+1) - 2.0*Z(k) + Z(k-1) ! 9

! Jacobian
YJA  = XG * YE * ZT
YJAI = 1.0 / YJA    ! 3

! 1st order
GX =  YE * ZT * YJAI
EY =  XG * ZT * YJAI
TZ =  XG * YE * YJAI ! 6

! Laplacian
C1 =  GX * GX
C2 =  EY * EY
C3 =  TZ * TZ
C7 = -XGG * C1 * GX
C8 = -YEE * C2 * EY
C9 = -ZTT * C3 * TZ ! 9

s1 = abs(C1 + 0.5 * C7)
s2 = abs(C1 - 0.5 * C7)
s3 = abs(C2 + 0.5 * C8)
s4 = abs(C2 - 0.5 * C8)
s5 = abs(C3 + 0.5 * C9)
s6 = abs(C3 - 0.5 * C9)
s7 = abs(2.0 * (C1 + C2 + C3))
ss = max(s1, s2, s3, s4, s5, s6, s7)

pvt(k,i,j) = 1.0/ss 

enddo
enddo
enddo
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END PARALLEL DO
#endif


return
end subroutine search_pivot
