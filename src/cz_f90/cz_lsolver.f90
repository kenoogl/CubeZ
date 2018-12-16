!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################

!> ********************************************************************
!! @brief Boundary condition
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     p     soution vector
!! @param [out]    dh    grid width
!! @param [in]     org   起点
!! @param [in]     nID   隣接ランクテーブル
!! @note resは積算
!!       側面の境界条件はディリクレ。pBiCGSTABの係数がディリクレを想定している実装のため。
!<
subroutine bc (sz, g, p, dh, org, nID)
implicit none
include 'cz_fparam.fi'
integer                                                :: i, j, k, ix, jx, kx, g
integer, dimension(3)                                  :: sz
integer, dimension(0:5)                                :: nID
real, dimension(3)                                     :: org
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: p
real                                                   :: pi, x, y, dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

pi = 2.0*asin(1.0)

! ZMINUS Dirichlet
if( nID(K_MINUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
do j=1,jx
do i=1,ix
  x = org(1) + dh*real(i-1)
  y = org(2) + dh*real(j-1)
  p(i,j,1) = sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END PARALLEL DO
endif


! ZPLUS Dirichlet
if( nID(K_PLUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
do j=1,jx
do i=1,ix
  x = org(1) + dh*real(i-1)
  y = org(2) + dh*real(j-1)
  p(i,j,kx) = sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END PARALLEL DO
endif


! XMINUS
if( nID(I_MINUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
  p(1,j,k) = 0.0 !p(2, j,k)
end do
end do
!$OMP END PARALLEL DO
endif


! XPLUS
if( nID(I_PLUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
  p(ix,j,k) = 0.0 !p(ix-1,j,k)
end do
end do
!$OMP END PARALLEL DO
endif


! YMINUS
if( nID(J_MINUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
  p(i,1, k) = 0.0 !p(i,2, k)
end do
end do
!$OMP END PARALLEL DO
endif


! YPLUS
if( nID(J_PLUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
  p(i,jx,k) = 0.0 !p(i,jx-1,k)
end do
end do
!$OMP END PARALLEL DO
endif

return
end subroutine bc


!> ********************************************************************
!! @brief point SOR法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     idx         インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     cf   係数
!! @param [in]     omg  加速係数
!! @param [in]     b    RHS vector
!! @param [out]    res  residual
!! @param [in,out] flop flop count
!<
subroutine psor (p, sz, idx, g, cf, omg, b, res, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  res
double precision                                       ::  flop, de
real                                                   ::  omg, dd, ss, dp, pp, bb, pn, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, b
real, dimension(7)                                     ::  cf

ix = sz(1)
jx = sz(2)
kx = sz(3)

res = 0.0

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

flop = flop + dble(ix*jx*kx)*28.0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(pp, bb, ss, dp, de, pn)
do k = kst, ked
do j = jst, jed
do i = ist, ied
  pp = p(i,j,k)
  bb = b(i,j,k)
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  p(i,j,k) = pn

  de  = dble( bb - (ss - pn * dd) )
  res = res + de*de
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine psor


!> **********************************************************************
!! @brief 緩和Jacobi法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     idx         インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     omg  加速係数
!! @param [in]     b    RHS vector
!! @param [in,out] res  residual
!! @param [out]    wk2  ワーク用配列
!! @param [in,out] flop flop count
!<
subroutine jacobi (p, sz, idx, g, cf, omg, b, res, wk2, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  res
double precision                                       ::  flop, de
real                                                   ::  omg, dd, ss, dp, pp, bb, pn, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, b, wk2
real, dimension(7)                                     ::  cf

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

flop = flop + dble((kx*jx*ix)*28.0)

!$OMP PARALLEL REDUCTION(+:res) PRIVATE(pp, bb, ss, dp, de, pn)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = kst, ked
do j = jst, jed
do i = ist, ied
  pp = p(i,j,k)
  bb = b(i,j,k)
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  wk2(i,j,k) = pn
  de  = dble( bb - (ss - pn * dd) )
  res = res + de*de
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = kst, ked
do j = jst, jed
do i = ist, ied
  p(i,j,k)=wk2(i,j,k)
end do
end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine jacobi


!> ********************************************************************
!! @brief 2-colored SOR法 stride memory access
!! @param [in,out] p     圧力
!! @param [in]     sz    配列長
!! @param [in]     idx   インデクス範囲
!! @param [in]     g     ガイドセル長
!! @param [in]     ofst  開始点オフセット
!! @param [in]     color グループ番号
!! @param [in]     omg   加速係数
!! @param [in]     b     RHS vector
!! @param [out]    res  residual
!! @param [in,out] flop  浮動小数演算数
!! @note resは積算
!<
subroutine psor2sma_core (p, sz, idx, g, cf, ofst, color, omg, b, res, flop)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, de
double precision                                       ::  res
real                                                   ::  omg, dd, ss, dp, pp, bb, pn, c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, b
integer                                                ::  ip, color, ofst
real, dimension(7)                                     ::  cf

ix = sz(1)
jx = sz(2)
kx = sz(3)
ip = ofst+color

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

flop = flop + dble(ix*jx*kx)*28.0*0.5d0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(pp, bb, ss, dp, de, pn)
do k=kst,ked
do j=jst,jed
do i=ist+mod(k+j+ip,2), ied, 2
  pp = p(i,j,k)
  bb = b(i,j,k)
  ss = c1 * p(i+1,j  ,k  ) &
     + c2 * p(i-1,j  ,k  ) &
     + c3 * p(i  ,j+1,k  ) &
     + c4 * p(i  ,j-1,k  ) &
     + c5 * p(i  ,j  ,k+1) &
     + c6 * p(i  ,j  ,k-1)
  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  p(i,j,k) = pn

  de  = dble( bb - (ss - pn * dd) )
  res = res + de*de
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine psor2sma_core


!> ********************************************************************
!! @brief Dirichlet source
!! @param [in]     b     RHS vector
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [out]    dh    grid width
!! @param [in]     nID   隣接ランクテーブル
!! @note resは積算
!<
subroutine src_dirichlet (b, sz, g, dh, nID)
implicit none
include 'cz_fparam.fi'
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  nID
real                                                   ::  dh, x, y, pi, dh2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  b

ix = sz(1)
jx = sz(2)
kx = sz(3)

pi = 2.0*asin(1.0)
dh2 = dh*dh

if( nID(K_MINUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
do j = 1, jx
do i = 1, ix
  x = dh*real(i)
  y = dh*real(j)
  !b(i,j,1 ) = - sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END PARALLEL DO
endif

if( nID(K_PLUS) < 0 ) then
  !$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
  do j = 1, jx
  do i = 1, ix
    x = dh*real(i)
    y = dh*real(j)
    !b(i,j,kx) = - sin(pi*x)*sin(pi*y)
  end do
  end do
  !$OMP END PARALLEL DO
endif

return
end subroutine src_dirichlet


!> ********************************************************************
!! @brief Line SOR法
!! @param [in,out] x    解ベクトル
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     omg  加速係数
!! @param [in]     q    TDMAの解
!! @param [out]    res  residual
!! @param [in,out] flop flop count
!<
subroutine tdma_sor (x, sz, idx, g, omg, q, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, de, res
real                                                   ::  omg, dd, ss, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x, q

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)
dd = 6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*26.0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(pp, ss, dp, de, pn)
do k = kst, ked
do j = jst, jed
do i = ist, ied
  pp = q(i,j,k) ! Jacobi的なので最新の値がよい
  ss = q(i+1,j  ,k  ) &
     + q(i-1,j  ,k  ) &
     + x(i  ,j+1,k  ) &
     + x(i  ,j-1,k  ) &
     + x(i  ,j  ,k+1) &
     + x(i  ,j  ,k-1)
  dp = ( ss/dd - pp ) * omg
  pn = pp + dp
  x(i,j,k) = pn

  de  = dble( - (ss - pn * dd) )
  res = res + de*de
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine tdma_sor


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     a    L_1 vector
!! @param [in]     b    D vector
!! @param [in]     c    U_1 vector
!! @param [in]     w    work vector (U_1)
!! @param [in,out] flop flop count
!!
!! Algorithm
!!
!! w_1 = c_1 / b_1
!! g_1 = d_1 / b_1
!!
!! ==> Forward elimination
!! for i=2, N do
!!  e = 1 / (b_i - a_i * w_{i-1})
!!  w_i = c_i * e                    // N-1までが意味がある
!!  g_i = (d_i - a_i * g_{i-1}) * e
!!
!!  x_n = g_n
!!
!! ==> Backward substitution
!! for i=N-1, 1 do
!!  x_i = g_i - w_i * x_{i+1}
!!
!! 上記のアルゴリズムで In-place にすると、
!! g と x を d で置き換えることができる
!!
!!  A = (a, b, c) の三重対角
!!  Ax = d を解く
!!  L = (a, e^{-1}) aが下副要素(L_1)、e^{-1}が対角要素
!!  U = (1, w) 対角要素1、wが上副要素(U_1)
!!
!!  LU分解: Lg=d, Ux=g
!<
subroutine tdma_0 (nx, d, a, b, c, w)
implicit none
integer                          ::  i, nx
real                             ::  e
real, dimension(nx)              ::  d, w, a, b, c


d(1) = d(1)/b(1);
w(1) = c(1)/b(1);

do i=2, nx
  e = 1.0 / (b(i) - a(i) * w(i-1))
  w(i) = e * c(i)
  d(i) = (d(i) - a(i) * d(i-1)) * e
end do

do i=nx-1, 1, -1
  d(i) = d(i) - w(i) * d(i+1)
end do

return
end subroutine tdma_0


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine tdma_1 (nx, d, cf, w, flop)
implicit none
integer                   ::  i, nx
double precision          ::  flop
real                      ::  a, b, c, e
real, dimension(nx)       ::  d, w
real, dimension(3)        ::  cf

a = cf(1)
b = cf(2)
c = cf(3)

flop = flop + dble(nx-1)*(14.0 + 2.0) + 8.0*2.0

d(1) = d(1)/b;
w(1) = c/b;

do i=2, nx
  e = 1.0 / (b - a * w(i-1))
  w(i) = e * c
  d(i) = (d(i) - a * d(i-1)) * e
end do

do i=nx-1, 1, -1
  d(i) = d(i) - w(i) * d(i+1)
end do

return
end subroutine tdma_1


!> ********************************************************************
!! @brief TDMAの右辺項
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in,out] flop flop count
!<
subroutine tdma_rhs (d, sz, idx, g, x, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x
real                                                   ::  r

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)
r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*4.0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k = kst, ked
do j = jst, jed
do i = ist, ied
  d(i,j,k) = ( x(i  ,j+1,k  ) &
           +   x(i  ,j-1,k  ) &
           +   x(i  ,j  ,k+1) &
           +   x(i  ,j  ,k-1) ) * r
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine tdma_rhs


!> ********************************************************************
!! @brief Line SOR法
!! @param [in,out] d    RHS / X
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     cf   coef
!! @param [in,out] flop flop count
!<
subroutine tdma_wrap (d, sz, idx, g, w, cf, flop)
implicit none
integer                                                ::  i, j, k, g, nn
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, dummy
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, w
real, dimension(3)                                     ::  cf

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)
nn = ied - ist + 1

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*16.0

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k = kst, ked
do j = jst, jed
  call tdma_1(nn, d(ist,j,k), cf, w(ist,j,k), dummy)
end do
end do
!$OMP END PARALLEL DO

return
end subroutine tdma_wrap

!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     w    work
!! @param [in]     cf   coef
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in]     ofst 開始点補正
!! @param [in]     color グループ
!! @param [in,out] flop flop count
!<
subroutine tdma_lsor2sma (d, sz, idx, g, x, w, cf, omg, res, ofst, color, flop)
implicit none
integer                                                ::  i, j, k, g, nn, ip, color, ofst
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, dummy, res, de
real                                                   ::  omg, dd, ss, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, w
real                                                   ::  r
real, dimension(3)                                     ::  cf

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)
r = 1.0/6.0
nn = ied - ist + 1
dd = 6.0
ip = ofst+color

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(4.0+16.0+13.0)

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(pp, ss, dp, de, pn, dummy)
do k = kst, ked
do j = jst+mod(k+ip,2), jed, 2
  do i = ist, ied
    d(i,j,k) = ( x(i  ,j+1,k  ) &
             +   x(i  ,j-1,k  ) &
             +   x(i  ,j  ,k+1) &
             +   x(i  ,j  ,k-1) ) * r
  end do

  call tdma_1(nn, d(ist,j,k), cf, w(ist,j,k), dummy)

  do i = ist, ied
    pp =  x(i,j,k)  ! ここはxでないと収束しない?
    ss =  d(i+1,j  ,k  ) &
       +  d(i-1,j  ,k  ) &
       +  x(i  ,j+1,k  ) &
       +  x(i  ,j-1,k  ) &
       +  x(i  ,j  ,k+1) &
       +  x(i  ,j  ,k-1)
    dp = ( ss*r - pp ) * omg
    pn = pp + dp
    x(i,j,k) = pn
    de  = dble( - (ss - pn * dd) )
    res = res + de*de
  end do
end do
end do
!$OMP END PARALLEL DO


return
end subroutine tdma_lsor2sma
