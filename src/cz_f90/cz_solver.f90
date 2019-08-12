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

! スレッド同期のオーバーヘッド抑制のため，単一のparallel regionとする
!$OMP PARALLEL

! ZMINUS Dirichlet
if( nID(K_MINUS) < 0 ) then
!$OMP DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
#ifdef _OPENACC
!$acc kernels
#endif
do j=1,jx
do i=1,ix
  x = org(1) + dh*real(i-1)
  y = org(2) + dh*real(j-1)
  p(1,i,j) = sin(pi*x)*sin(pi*y)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO NOWAIT
endif


! ZPLUS Dirichlet
if( nID(K_PLUS) < 0 ) then
!$OMP DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
#ifdef _OPENACC
!$acc kernels
#endif
do j=1,jx
do i=1,ix
  x = org(1) + dh*real(i-1)
  y = org(2) + dh*real(j-1)
  p(kx,i,j) = sin(pi*x)*sin(pi*y)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO
endif


! XMINUS
if( nID(I_MINUS) < 0 ) then
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#ifdef _OPENACC
!$acc kernels
#endif
do k=1,kx
do j=1,jx
  p(k,1,j) = 0.0 !p(2, j,k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO NOWAIT
endif


! XPLUS
if( nID(I_PLUS) < 0 ) then
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#ifdef _OPENACC
!$acc kernels
#endif
do k=1,kx
do j=1,jx
  p(k,ix,j) = 0.0 !p(ix-1,j,k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO
endif


! YMINUS
if( nID(J_MINUS) < 0 ) then
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#ifdef _OPENACC
!$acc kernels
#endif
do k=1,kx
do i=1,ix
  p(k,i,1) = 0.0 !p(i,2, k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO NOWAIT
endif


! YPLUS
if( nID(J_PLUS) < 0 ) then
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#ifdef _OPENACC
!$acc kernels
#endif
do k=1,kx
do i=1,ix
  p(k,i,jx) = 0.0 !p(i,jx-1,k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO
endif

!$OMP END PARALLEL

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
real                                                   ::  c1, c2, c3, c4, c5, c6, res1
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
real, dimension(7)                                     ::  cf

res1 = 0.0

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
!$OMP REDUCTION(+:res1) &
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
  res1 = res1 + dp*dp
end do
end do
end do
!$OMP END PARALLEL DO

res = res + real(res1, kind=8)

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
real                                                   ::  c1, c2, c3, c4, c5, c6, res1
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b, wk2
real, dimension(7)                                     ::  cf

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0

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
!$OMP REDUCTION(+:res1)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#ifdef _OPENACC
!$acc kernels
#endif
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
  wk2(k,i,j) = pn
  res1 = res1 + dp*dp
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO

!$OMP DO SCHEDULE(static) COLLAPSE(2)
#ifdef _OPENACC
!$acc kernels
#endif
do j = jst, jed
do i = ist, ied
do k = kst, ked
  p(k,i,j)=wk2(k,i,j)
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO

!$OMP END PARALLEL

res = res + real(res1, kind=8)

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
real                                                   ::  c1, c2, c3, c4, c5, c6, res1
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

res1 = 0.0

flop = flop + 18.0d0*0.5d0  &
     * dble(ied-ist+1) &
     * dble(jed-jst+1) &
     * dble(ked-kst+1)

!$OMP PARALLEL REDUCTION(+:res1) &
!$OMP PRIVATE(pp, bb, ss, dp, pn)
!$OMP DO SCHEDULE(static) COLLAPSE(2) 
#ifdef _OPENACC
!$acc kernels
!$acc loop independent gang reduction(+:res1)
do j=jst,jed
!$acc loop independent gang reduction(+:res1)
do i=ist,ied
!$acc loop independent vector(128) reduction(+:res1)
do k=kst+mod(i+j+kp,2), ked, 2
#else
!pgi$ ivdep
do j=jst,jed
!pgi$ novector
do i=ist,ied
!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ vector
do k=kst+mod(i+j+kp,2), ked, 2
#endif
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
  res1 = res1 + dp*dp
end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO
!$OMP END PARALLEL

res = res + real(res1, kind=8)

return
end subroutine psor2sma_core


!********************************************************************************
subroutine pcr_rb (sz, idx, g, pn, ofst, color, x, msk, rhs, a, c, d, a1, c1, d1, omg, res, flop)
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
real, dimension(-1:sz(3)+2)              ::  a, c, d, a1, c1, d1
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, cc1, cc2, f1, f2

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0
r = 1.0/6.0

flop = flop + dble(          &
  (jed-jst+1)*(ied-ist+1)* ( &
     (ked-kst+1)* 6.0        &  ! Source
   + (ked-kst+1)*(pn-1)*14.0 &  ! PCR
   + 2**(pn-1)*9.0                 &
   + (ked-kst+1)*6.0         &  ! Relaxation
     + 6.0 )                 &  ! BC
  ) * 0.5


ip = ofst + color


! #ifdef _OPENACC
! !$acc kernels
! !$acc loop independent gang reduction(+:res1)
! do j=jst, jed
! !$acc loop independent vector(128) reduction(+:res1)
! do i=ist+mod(j+ip,2), ied, 2
! #else
! !$OMP PARALLEL reduction(+:res1) &
! !$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
! !$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
! !$OMP private(a, c, d, a1, c1, d1)
! !$OMP DO SCHEDULE(static)
! !pgi$ ivdep
! do j=jst, jed
! do i=ist+mod(j+ip,2), ied, 2
! #endif

#ifdef _OPENACC
!$acc kernels
#else
!$OMP PARALLEL &
!$OMP reduction(+:res) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a, c, d, a1, c1, d1)
!$OMP DO SCHEDULE(static)
#endif
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
!NEC$ IVDEP
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
!res1 = res1 + dp*dp
res = res + dp*dp
end do

end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO
!$OMP END PARALLEL

!res = res + real(res1, kind=8)

return
end subroutine pcr_rb


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
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, cc1, cc2, f1, f2

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0
r = 1.0/6.0

flop = flop + dble(          &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ 2**(pn-1)*9.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)


#ifdef _OPENACC
!$acc kernels
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a, c, d, a1, c1, d1)
!$OMP DO SCHEDULE(static) Collapse(2)
#endif
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
!pgi$ ivdep
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
!NEC$ IVDEP
!pgi$ ivdep
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



! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res1 = res1 + dp*dp
end do

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif


res = res + real(res1, kind=8)

return
end subroutine pcr


!********************************************************************************
! pcr for vector (Aurora and GPU)
subroutine pcr_eda (sz, idx, g, pn, x, msk, rhs, a1, c1, d1, omg, res, flop)
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
real, dimension(1-g:sz(3)+g)             ::  a1, c1, d1
real, dimension(:), allocatable          ::  a, c, d
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, cc1, cc2, f1, f2

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0
r = 1.0/6.0
s = 2**(pn-1)

flop = flop + dble(          &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ s*9.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)

s = 2**(pn-2)
kl = kst-s
kr = ked+s
allocate(a(kl:kr))
allocate(c(kl:kr))
allocate(d(kl:kr))
do k=kl, kr
a(k) = 0.0
c(k) = 0.0
d(k) = 0.0
end do


#ifdef _OPENACC
!$acc kernels
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a1, c1, d1) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static) Collapse(2)
#endif
do j=jst, jed
do i=ist, ied

! Reflesh coef. due to override
!a(kst) = 0.0
do k=kst+1, ked
a(k) = -r
end do

do k=kst, ked-1
c(k) = -r
end do
!c(ked) = 0.0

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


! PCR  最終段の一つ手前で停止
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-s) - cp * a(k+s) )
a1(k) =  -e * ap * a(k-s)
c1(k) =  -e * cp * c(k+s)
d1(k) =   e * ( d(k) - ap * d(k-s) - cp * d(k+s) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = a1(k)
c(k) = c1(k)
d(k) = d1(k)
end do

end do ! p反復


! 最終段の反転 512のとき pn=9, s=256
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
do k = kst, kst+s-1 ! 2, 2+256-1=257
cc1 = c(k)
aa2 = a(k+s)
f1  = d(k)
f2  = d(k+s)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k ) = dd1
d1(k+s)= dd2
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res1 = res1 + dp*dp
end do

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

deallocate(a)
deallocate(c)
deallocate(d)

res = res + real(res1, kind=8)

return
end subroutine pcr_eda


!********************************************************************************
! pcr for vector (Aurora and GPU)
subroutine pcr_esa (sz, idx, g, pn, s, x, msk, rhs, a, c, d, a1, c1, d1, omg, res, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, p, ss, s
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(1-g:sz(3)+g)             ::  a1, c1, d1
real, dimension(idx(4)-s:idx(5)+s)       ::  a, c, d
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, cc1, cc2, f1, f2

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0
r = 1.0/6.0

flop = flop + dble(          &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ 2**(pn-1)*9.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)

#ifdef _OPENACC
!$acc kernels
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(ap, cp, e, ss, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a1, c1, d1) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static) Collapse(2)
#endif
do j=jst, jed
do i=ist, ied

! Reflesh coef. due to override
!a(kst) = 0.0
do k=kst+1, ked
a(k) = -r
end do

do k=kst, ked-1
c(k) = -r
end do
!c(ked) = 0.0

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


! PCR  最終段の一つ手前で停止
do p=1, pn-1
ss = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-ss) - cp * a(k+ss) )
a1(k) =   -e * ap * a(k-ss)
c1(k) =   -e * cp * c(k+ss)
d1(k) =   e * ( d(k) - ap * d(k-ss) - cp * d(k+ss) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = a1(k)
c(k) = c1(k)
d(k) = d1(k)
end do

end do ! p反復


! 最終段の反転 512のとき pn=9, ss=256
ss = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
do k = kst, kst+ss-1 ! 2, 2+256-1=257
cc1 = c(k)
aa2 = a(k+ss)
f1  = d(k)
f2  = d(k+ss)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k   ) = dd1
d1(k+ss) = dd2
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res1 = res1 + dp*dp
end do

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

res = res + real(res1, kind=8)

return
end subroutine pcr_esa


!********************************************************************************
subroutine pcr_rb_esa (sz, idx, g, pn, ofst, color, s, x, msk, rhs, a, c, d, a1, c1, d1, omg, res, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, p, color, ip, ofst, ss, s
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(-1:sz(3)+2)              ::  a1, c1, d1
real, dimension(idx(4)-s:idx(5)+s)       ::  a, c, d
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, cc1, cc2, f1, f2

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0
r = 1.0/6.0

flop = flop + dble(          &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ 2**(pn-1)*9.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
) * 0.5

ip = ofst + color


#ifdef _OPENACC
!$acc kernels
!$acc loop independent gang reduction(+:res1)
do j=jst, jed
!$acc loop independent vector(128) reduction(+:res1)
do i=ist+mod(j+ip,2), ied, 2
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(ap, cp, e, ss, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a1, c1, d1) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static)
!pgi$ ivdep
do j=jst, jed
do i=ist+mod(j+ip,2), ied, 2
#endif

! Reflesh coef. due to override
!a(kst) = 0.0
do k=kst+1, ked
a(k) = -r
end do

do k=kst, ked-1
c(k) = -r
end do
!c(ked) = 0.0

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


! PCR  最終段の一つ手前で停止
do p=1, pn-1
ss = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-ss) - cp * a(k+ss) )
a1(k) =   -e * ap * a(k-ss)
c1(k) =   -e * cp * c(k+ss)
d1(k) =   e * ( d(k) - ap * d(k-ss) - cp * d(k+ss) )
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
ss = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
do k = kst, kst+s-1
cc1 = c(k)
aa2 = a(k+ss)
f1  = d(k)
f2  = d(k+ss)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k   ) = dd1
d1(k+ss) = dd2
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
res1 = res1 + dp*dp
end do

end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END DO

!$OMP END PARALLEL

res = res + real(res1, kind=8)

return
end subroutine pcr_rb_esa
