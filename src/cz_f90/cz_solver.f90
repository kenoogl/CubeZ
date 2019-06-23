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
integer                                                ::  i, j, k, g
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

flop = flop + 18.0     &
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
integer                                                ::  i, j, k, g
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

flop = flop + 18.0  &
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
integer                                                ::  i, j, k, g
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

flop = flop + 18.0d0*0.5d0  &
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


!********************************************************************************
subroutine PCR_RB (sz, idx, g, pn, ofst, color, x, msk, rhs, a, c, d, a1, c1, d1, omg, res, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, kl, kr, s, p, color, ip, ofst
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(-1:sz(3)+2)             ::  a, c, d, a1, c1, d1
real                                     ::  r, ap, cp, e, pp, dp
real                                     ::  jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
s = 2**(pn-1)

flop = flop + dble(          &
  (jed-jst+1)*(ied-ist+1)* ( &
     (ked-kst+1)* 6.0        &  ! Source
   + (ked-kst+1)*(pn-1)*14.0 &  ! PCR
   + 2*s*9.0                 &
   + (ked-kst-2*s+1)*25.0    &
   + (ked-kst+1)*6.0         &  ! Relaxation
     + 6.0 )                 &  ! BC
  ) * 0.5


ip = ofst + color

!$OMP PARALLEL reduction(+:res) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, a1, c1, d1)

!$OMP DO SCHEDULE(static)
do j=jst, jed
do i=ist+mod(j+ip,2), ied, 2

! Reflesh coef. due to override
a(kst) = 0.0
do k=kst+1, ked
a(k) = -r
end do

do k=kst, ked-1
c(k) = -r
end do
c(ked) = 0.0

! Source
!dir$ vector aligned
!dir$ simd
do k = kst, ked
d(k) = (   ( x(k, i  , j-1)        &
       +     x(k, i  , j+1)        &
       +     x(k, i-1, j  )        &
       +     x(k, i+1, j  ) - rhs(k, i, j) ) * r ) &
       *   msk(k, i, j)
end do ! 6 flops

! BC  6 flops
d(kst) = ( d(kst) + x(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked) = ( d(ked) + x(ked+1, i, j) * r ) * msk(ked, i, j)

!d(kst) = ( d(kst) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
!d(ked) = ( d(ked) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(kl) - cp * a(kr) )
a1(k) =  -e * ap * a(kl)
c1(k) =  -e * cp * c(kr)
d1(k) =   e * ( d(k) - ap * d(kl) - cp * d(kr) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = a1(k)
c(k) = c1(k)
d(k) = d1(k)
end do

end do ! p反復


! 最終段の反転
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
do k = kst, kst+s-1
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = c(k)
aa2 = a(kr)
f1  = d(k)
f2  = d(kr)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k ) = dd1
d1(kr) = dd2
end do


!dir$ vector aligned
!dir$ simd
do k = kst+s, ked-s
kl  = max(k-s, kst-1)
kr  = min(k+s, ked+1)
cc1 = c(kr)
aa2 = a(k)
cc2 = c(k)
aa3 = a(kl)
f1  = d(kr)
f2  = d(k)
f3  = d(kl)
jj = 1.0 / (1.0 - cc2 * aa3 - cc1 * aa2)
dd1 = ( f1 * (3.0-cc2*aa3) - cc1*f2 ) * jj
dd2 = (1.0 - f1*aa2 + 2.0*f2 - cc2*f3) * jj
dd3 = (1.0 + 2.0*f3 - aa3*f2 - aa2*cc1) * jj
d1(kr) = dd1
d1(k ) = dd2
d1(kl) = dd3
end do


!dir$ vector aligned
!dir$ simd
do k = ked-s+1, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = c(kl)
aa2 = a(k)
f1  = d(kl)
f2  = d(k)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(kl) = dd1
d1(k ) = dd2
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res = res + real(dp*dp, kind=8)
end do

end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine PCR_RB


!********************************************************************************
subroutine pcr (sz, idx, g, pn, x, msk, rhs, a, c, d, a1, c1, d1, omg, res, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, kl, kr, s, p
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(1-g:sz(3)+g)             ::  a, c, d, a1, c1, d1
real                                     ::  r, ap, cp, e, pp, dp
real                                     ::  jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
s = 2**(pn-1)

flop = flop + dble(          &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ 2*s*9.0                 &
+ (ked-kst-2*s+1)*25.0    &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)


!$OMP PARALLEL reduction(+:res) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, a1, c1, d1)

!$OMP DO SCHEDULE(static) Collapse(2)
do j=jst, jed
do i=ist, ied

! Reflesh coef. due to override
a(kst) = 0.0
do k=kst+1, ked
a(k) = -r
end do

do k=kst, ked-1
c(k) = -r
end do
c(ked) = 0.0

! Source
!dir$ vector aligned
!dir$ simd
do k = kst, ked
d(k) = (   ( x(k, i  , j-1)        &
+     x(k, i  , j+1)        &
+     x(k, i-1, j  )        &
+     x(k, i+1, j  ) - rhs(k, i, j) ) * r ) &
*   msk(k, i, j)
end do ! 6 flops

! BC  6 flops
d(kst) = ( d(kst) + x(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked) = ( d(ked) + x(ked+1, i, j) * r ) * msk(ked, i, j)

!d(kst) = ( d(kst) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
!d(ked) = ( d(ked) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(kl) - cp * a(kr) )
a1(k) =  -e * ap * a(kl)
c1(k) =  -e * cp * c(kr)
d1(k) =   e * ( d(k) - ap * d(kl) - cp * d(kr) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = a1(k)
c(k) = c1(k)
d(k) = d1(k)
end do

end do ! p反復


! 最終段の反転
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
do k = kst, kst+s-1
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = c(k)
aa2 = a(kr)
f1  = d(k)
f2  = d(kr)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k ) = dd1
d1(kr) = dd2
end do


!dir$ vector aligned
!dir$ simd
do k = kst+s, ked-s
kl  = max(k-s, kst-1)
kr  = min(k+s, ked+1)
cc1 = c(kr)
aa2 = a(k)
cc2 = c(k)
aa3 = a(kl)
f1  = d(kr)
f2  = d(k)
f3  = d(kl)
jj = 1.0 / (1.0 - cc2 * aa3 - cc1 * aa2)
dd1 = ( f1 * (3.0-cc2*aa3) - cc1*f2 ) * jj
dd2 = (1.0 - f1*aa2 + 2.0*f2 - cc2*f3) * jj
dd3 = (1.0 + 2.0*f3 - aa3*f2 - aa2*cc1) * jj
d1(kr) = dd1
d1(k ) = dd2
d1(kl) = dd3
end do


!dir$ vector aligned
!dir$ simd
do k = ked-s+1, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = c(kl)
aa2 = a(k)
f1  = d(kl)
f2  = d(k)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(kl) = dd1
d1(k ) = dd2
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res = res + real(dp*dp, kind=8)
end do

end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine pcr
