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
subroutine bc_k (sz, g, p, dh, org, nID)
implicit none
include 'cz_fparam.fi'
integer                                                :: i, j, k, ix, jx, kx, g
integer, dimension(3)                                  :: sz
integer, dimension(0:5)                                :: nID
real, dimension(3)                                     :: org
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) :: p
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
  p(1,i,j) = sin(pi*x)*sin(pi*y)
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
  p(kx,i,j) = sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END PARALLEL DO
endif


! XMINUS
if( nID(I_MINUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
  p(k,1,j) = 0.0 !p(2, j,k)
end do
end do
!$OMP END PARALLEL DO
endif


! XPLUS
if( nID(I_PLUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
  p(k,ix,j) = 0.0 !p(ix-1,j,k)
end do
end do
!$OMP END PARALLEL DO
endif


! YMINUS
if( nID(J_MINUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
  p(k,i,1) = 0.0 !p(i,2, k)
end do
end do
!$OMP END PARALLEL DO
endif


! YPLUS
if( nID(J_PLUS) < 0 ) then
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
  p(k,i,jx) = 0.0 !p(i,jx-1,k)
end do
end do
!$OMP END PARALLEL DO
endif

return
end subroutine bc_k


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
double precision                                       ::  flop
real                                                   ::  omg, dd, ss, dp, pp, bb, pn
real                                                   ::  c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
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

flop = flop + 25.0     &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(pp, bb, ss, dp, pn)
do j = jst, jed
do i = ist, ied
do k = kst, ked
  pp = p(k,i,j)
  bb = b(k,i,j)
  ss = c1 * p(k  , i+1,j  ) &
     + c2 * p(k  , i-1,j  ) &
     + c3 * p(k  , i  ,j+1) &
     + c4 * p(k  , i  ,j-1) &
     + c5 * p(k+1, i  ,j  ) &
     + c6 * p(k-1, i  ,j  )
  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  p(k,i,j) = pn
  res = res + dp*dp
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
double precision                                       ::  flop
real                                                   ::  omg, dd, ss, dp, pp, bb, pn
real                                                   ::  c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b, wk2
real, dimension(7)                                     ::  cf
!dir$ assume_aligned p:64, b:64, wk2:64

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

flop = flop + 25.0  &
            * dble(ied-ist+1) &
            * dble(jed-jst+1) &
            * dble(ked-kst+1)

!$OMP PARALLEL PRIVATE(pp, bb, ss, dp, pn) &
!$OMP REDUCTION(+:res)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j = jst, jed
do i = ist, ied

!dir$ vector aligned
!dir$ simd
do k = kst, ked
  pp = p(k,i,j)
  bb = b(k,i,j)
  ss = c1 * p(k  , i+1,j  ) &
     + c2 * p(k  , i-1,j  ) &
     + c3 * p(k  , i  ,j+1) &
     + c4 * p(k  , i  ,j-1) &
     + c5 * p(k+1, i  ,j  ) &
     + c6 * p(k-1, i  ,j  )
  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  wk2(k,i,j) = pn
  res = res + dp*dp
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j = jst, jed
do i = ist, ied

!dir$ vector aligned
!dir$ simd
do k = kst, ked
  p(k,i,j)=wk2(k,i,j)
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
double precision                                       ::  flop
double precision                                       ::  res
real                                                   ::  omg, dd, ss, dp, pp, bb, pn
real                                                   ::  c1, c2, c3, c4, c5, c6
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
integer                                                ::  kp, color, ofst
real, dimension(7)                                     ::  cf
!dir$ assume_aligned p:64, b:64

ix = sz(1)
jx = sz(2)
kx = sz(3)
kp = ofst+color

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

flop = flop + 25.0d0*0.5d0  &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(pp, bb, ss, dp, pn)
do j=jst,jed
do i=ist,ied

!dir$ vector aligned
!dir$ simd
do k=kst+mod(i+j+kp,2), ked, 2
  pp = p(k,i,j)
  bb = b(k,i,j)
  ss = c1 * p(k  , i+1,j  ) &
     + c2 * p(k  , i-1,j  ) &
     + c3 * p(k  , i  ,j+1) &
     + c4 * p(k  , i  ,j-1) &
     + c5 * p(k+1, i  ,j  ) &
     + c6 * p(k-1, i  ,j  )
  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  p(k,i,j) = pn
  res = res + dp*dp
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine psor2sma_core
