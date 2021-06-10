!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University.
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
#ifdef _OPENACC
!$acc kernels
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static) PRIVATE(x, y)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
#endif
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
#else
!$OMP END DO NOWAIT
#endif
endif


! ZPLUS Dirichlet
if( nID(K_PLUS) < 0 ) then
#ifdef _OPENACC
!$acc kernels
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static) PRIVATE(x, y)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y)
#endif
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
#else
!$OMP END DO
#endif
endif


! XMINUS
if( nID(I_MINUS) < 0 ) then
#ifdef _OPENACC
!$acc kernels
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
do k=1,kx
do j=1,jx
p(k,1,j) = 0.0 !p(2, j,k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO NOWAIT
#endif
endif


! XPLUS
if( nID(I_PLUS) < 0 ) then
#ifdef _OPENACC
!$acc kernels
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
do k=1,kx
do j=1,jx
p(k,ix,j) = 0.0 !p(ix-1,j,k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
#endif
endif


! YMINUS
if( nID(J_MINUS) < 0 ) then
#ifdef _OPENACC
!$acc kernels
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
do k=1,kx
do i=1,ix
p(k,i,1) = 0.0 !p(i,2, k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO NOWAIT
#endif
endif


! YPLUS
if( nID(J_PLUS) < 0 ) then
#ifdef _OPENACC
!$acc kernels
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
do k=1,kx
do i=1,ix
p(k,i,jx) = 0.0 !p(i,jx-1,k)
end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
#endif
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


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3) reduction(+:res1)
#else
!$OMP PARALLEL PRIVATE(pp, bb, ss, dp, pn) &
!$OMP REDUCTION(+:res1)
! auroraはここにcollapseを入れると完全に並列化しない
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
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
#else
!$OMP END DO NOWAIT
#endif


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3)
#else
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
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
#else
!$OMP END DO
!$OMP END PARALLEL
#endif


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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent gang reduction(+:res1)
do j=jst,jed
!$acc loop independent gang reduction(+:res1)
do i=ist,ied
!$acc loop independent vector(128) reduction(+:res1)
do k=kst+mod(i+j+kp,2), ked, 2
#else
!$OMP PARALLEL REDUCTION(+:res1) &
!$OMP PRIVATE(pp, bb, ss, dp, pn)
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
do j=jst,jed
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
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

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


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) gang private(a, c, d, a1, c1, d1) reduction(+:res)
#else
!$OMP PARALLEL &
!$OMP reduction(+:res) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a, c, d, a1, c1, d1)
!$OMP DO SCHEDULE(static) collapse(2)
#endif
do j=jst, jed
do i=ist, ied
if(mod(i+j,2) /= color) cycle

! do i=ist+mod(j+ip,2), ied, 2



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
!$acc loop seq
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
!$acc loop independent
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
!$acc loop reduction(+:res)
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
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

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
integer                                  ::  i, j, k, kl, km, kr, s, p
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(1-g:sz(3)+g)             ::  a, c, d, a1, c1, d1
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, dd3, dd4, aa2, aa3, aa4, cc1, cc2, cc3
real                                     ::  inv_detA, detA1, detA2, detA3, detA4


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
+ (ked-kst+1)*(pn-2)*14.0 & ! PCR4x4
+ 2**(pn-2)*74.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) reduction(+:res1) &
!$acc& private(a, c, d, a1, c1, d1) &
!$acc& private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2)
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(kl, km, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, dd4, aa2, aa3, aa4, cc1, cc2, cc3) &
!$OMP private(a, c, d, a1, c1, d1) &
!$OMP private(inv_detA, detA1, detA2, detA3, detA4)
!$OMP DO SCHEDULE(static) Collapse(2)
#endif
do j=jst, jed
do i=ist, ied

! Reflesh coef. due to override
! 係数行列a[1-g:sz(3)+g]を初期化
! r = 1.0 / 6.0
a(kst-1) = 0.0
a(kst) = 0.0
do k=kst+1, ked
a(k) = -r
end do
a(ked+1) = 0.0

c(kst-1) = 0.0
do k=kst, ked-1
c(k) = -r
end do
c(ked) = 0.0
c(ked+1) = 0.0

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
!$acc loop seq
do p=1, pn-2
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
s = 2**(pn-2)


!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
!$acc loop independent
do k = kst, kst+s-1
kl = min(k+  s, ked+1)
km = min(k+2*s, ked+1)
kr = min(k+3*s, ked+1)

! A = \\
! \begin{pmatrix}
! 1   & cc1 & 0   & 0   \\
! aa2 & 1   & cc2 & 0   \\
! 0   & aa3 & 1   & cc3 \\
! 0   & 0   & aa4 & 1   \\
! \end{pmatrix}
cc1 = c(k)
cc2 = c(kl)
cc3 = c(km)
aa2 = a(kl)
aa3 = a(km)
aa4 = a(kr)

! (dd1, dd2, dd3, dd4 ) = ( d(k) & d(kl) & d(km) & d(kr) )
dd1 = d(k)
dd2 = d(kl)
dd3 = d(km)
dd4 = d(kr)

! 1.0 / |A|
inv_detA= 1.0 / (1.0 - aa4*cc3 - aa3*cc2 - aa2*cc1*(1.0 - cc3*aa4))

! |A_i|
! A_i = A の第 i-列 (i = 1, 2, …, n) を系の右辺である d で置き換えて得られる行列
detA1 = -cc3*(aa4*dd1 + cc1*cc2*dd4 - aa4*cc1*dd2) &
+       dd1 + cc1*cc2*dd3 - aa3*cc2*dd1 - cc1*dd2

detA2 = dd2 + cc2*cc3*dd4 - aa4*cc3*dd2 - cc2*dd3 &
-       aa2*(dd1 - aa4*cc3*dd1)

detA3 = dd3 - cc3*dd4 - aa3*dd2 &
-       aa2*(cc1*dd3 - cc1*cc3*dd4 - aa3*dd1)

detA4 = dd4 + aa3*aa4*dd2 - aa4*dd3 - aa3*cc2*dd4 &
-       aa2*(cc1*dd4 + aa3*aa4*dd1 - aa4*cc1*dd3)


! Cramer's rule
d1(k)  = detA1 * inv_detA
d1(kl) = detA2 * inv_detA
d1(km) = detA3 * inv_detA
d1(kr) = detA4 * inv_detA
end do



! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
!$acc loop reduction(+:res1)
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

flop = flop + dble(          &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ 2**(pn-1)*9.0                 &
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



!$OMP PARALLEL reduction(+:res1) &
!$OMP private(ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a1, c1, d1) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static) Collapse(2)
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
!$OMP END DO
!$OMP END PARALLEL


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
integer                                  ::  i, j, k, kl, km, kr, p, sq, s
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(1-g:sz(3)+g)             ::  a1, c1, d1
real, dimension(idx(4)-s:idx(5)+s)       ::  a, c, d
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, dd3, dd4
real                                     ::  aa2, aa3, aa4, cc1, cc2, cc3
real                                     ::  inv_detA, detA1, detA2, detA3, detA4

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
+ (ked-kst+1)*(pn-2)*14.0 &  ! PCR
+ 2**(pn-2)*78.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)

#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) gang reduction(+:res1) &
!$acc& private(a, c, d, a1, c1, d1) &
!$acc& private(ap, cp, e, sq, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2)
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(ap, cp, e, sq, p, k, km, kl, kr, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, dd4, aa2, aa3, aa4, cc1, cc2, cc3) &
!$OMP private(a1, c1, d1) &
!$OMP private(inv_detA, detA1, detA2, detA3, detA4) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static) Collapse(2)
#endif
do j=jst, jed
do i=ist, ied

! Reflesh coef. due to override

do sq=0, s
a(kst-sq) = 0.0
a(ked+sq) = 0.0
end do
! override a(ked) = -r

do k=kst+1, ked
a(k) = -r
end do
do k=ked+1, ked+s
a(k) = 0.0
end do


do sq=0, s
c(kst-sq) = 0.0
c(ked+sq) = 0.0
end do
do k=kst, ked-1
c(k) = -r
end do
! over write c(kst)


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


! PCR  最終段の2つ手前で停止
!$acc loop seq
do p=1, pn-2
sq = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-sq) - cp * a(k+sq) )
a1(k) =   -e * ap * a(k-sq)
c1(k) =   -e * cp * c(k+sq)
d1(k) =   e * ( d(k) - ap * d(k-sq) - cp * d(k+sq) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = a1(k)
c(k) = c1(k)
d(k) = d1(k)
end do

end do ! p反復


! 最終段の反転 512のとき pn=9, sq=256
sq = 2**(pn-2)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
!$acc loop independent
do k = kst, kst+sq-1 ! 2, 2+256-1=257
kl = k + sq
km = k + 2*sq
kr = k + 3*sq

cc1 = c(k)
cc2 = c(kl)
cc3 = c(km)
aa2 = a(kl)
aa3 = a(km)
aa4 = a(kr)

! (dd1, dd2, dd3, dd4 ) = ( d(k) & d(kl) & d(km) & d(kr) )
dd1 = d(k)
dd2 = d(kl)
dd3 = d(km)
dd4 = d(kr)

! 1.0 / |A|
inv_detA= 1.0 / (1.0 - aa4*cc3 - aa3*cc2 - aa2*cc1*(1.0 - cc3*aa4))

! |A_i|
! A_i = A の第 i-列 (i = 1, 2, …, n) を系の右辺である d で置き換えて得られる行列
detA1 = -cc3*(aa4*dd1 + cc1*cc2*dd4 - aa4*cc1*dd2) &
+       dd1 + cc1*cc2*dd3 - aa3*cc2*dd1 - cc1*dd2

detA2 = dd2 + cc2*cc3*dd4 - aa4*cc3*dd2 - cc2*dd3 &
-       aa2*(dd1 - aa4*cc3*dd1)

detA3 = dd3 - cc3*dd4 - aa3*dd2 &
-       aa2*(cc1*dd3 - cc1*cc3*dd4 - aa3*dd1)

detA4 = dd4 + aa3*aa4*dd2 - aa4*dd3 - aa3*cc2*dd4 &
-       aa2*(cc1*dd4 + aa3*aa4*dd1 - aa4*cc1*dd3)


! Cramer's rule
d1(k)  = detA1 * inv_detA
d1(kl) = detA2 * inv_detA
d1(km) = detA3 * inv_detA
d1(kr) = detA4 * inv_detA
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
!$acc loop reduction(+:res1)
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
integer                                  ::  i, j, k, kl, km, kr, p, color, ip, ofst, sq, s
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(-1:sz(3)+2)              ::  a1, c1, d1
real, dimension(idx(4)-s:idx(5)+s)       ::  a, c, d
real                                     ::  r, ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, dd3, dd4
real                                     ::  aa2, aa3, aa4, cc1, cc2, cc3
real                                     ::  inv_detA, detA1, detA2, detA3, detA4


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
+ (ked-kst+1)*(pn-2)*14.0 &  ! PCR
+ 2**(pn-2)*78.0                 &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
) * 0.5

! unused variable
ip = ofst + color


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) gang reduction(+:res1) &
!$acc& private(a, c, d, a1, c1, d1) &
!$acc& private(ap, cp, e, sq, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2)
#else
!$OMP PARALLEL reduction(+:res1) &
!$OMP private(ap, cp, e, sq, p, k, km, kl, kr, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, dd4, aa2, aa3, aa4, cc1, cc2, cc3) &
!$OMP private(a1, c1, d1) &
!$OMP private(inv_detA, detA1, detA2, detA3, detA4) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static) collapse(2)
#endif
do j=jst, jed
do i=ist, ied
if(mod(i+j,2) /= color) cycle

!do i=ist+mod(j+ip,2), ied, 2


! Reflesh coef. due to override
do k=kst-s, kst
a(k) = 0.0
end do
do k=kst+1, ked
a(k) = -r
end do
do k=ked+1, ked+s
a(k) = 0.0
end do

do k=kst-s, kst-1
c(k) = 0.0
end do
do k=kst, ked-1
c(k) = -r
end do
do k=ked, ked+s
c(k) = 0.0
end do

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


! PCR  最終段の2つ手前で停止
!$acc loop seq
do p=1, pn-2
sq = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-sq) - cp * a(k+sq) )
a1(k) =   -e * ap * a(k-sq)
c1(k) =   -e * cp * c(k+sq)
d1(k) =   e * ( d(k) - ap * d(k-sq) - cp * d(k+sq) )
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
sq = 2**(pn-2)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!$acc loop independent
do k = kst, kst+sq-1
kl = k + sq
km = k + 2*sq
kr = k + 3*sq

cc1 = c(k)
cc2 = c(kl)
cc3 = c(km)
aa2 = a(kl)
aa3 = a(km)
aa4 = a(kr)

! (dd1, dd2, dd3, dd4 ) = ( d(k) & d(kl) & d(km) & d(kr) )
dd1 = d(k)
dd2 = d(kl)
dd3 = d(km)
dd4 = d(kr)

! 1.0 / |A|
inv_detA= 1.0 / (1.0 - aa4*cc3 - aa3*cc2 - aa2*cc1*(1.0 - cc3*aa4))

! |A_i|
! A_i = A の第 i-列 (i = 1, 2, …, n) を系の右辺である d で置き換えて得られる行列
detA1 = -cc3*(aa4*dd1 + cc1*cc2*dd4 - aa4*cc1*dd2) &
+       dd1 + cc1*cc2*dd3 - aa3*cc2*dd1 - cc1*dd2

detA2 = dd2 + cc2*cc3*dd4 - aa4*cc3*dd2 - cc2*dd3 &
-       aa2*(dd1 - aa4*cc3*dd1)

detA3 = dd3 - cc3*dd4 - aa3*dd2 &
-       aa2*(cc1*dd3 - cc1*cc3*dd4 - aa3*dd1)

detA4 = dd4 + aa3*aa4*dd2 - aa4*dd3 - aa3*cc2*dd4 &
-       aa2*(cc1*dd4 + aa3*aa4*dd1 - aa4*cc1*dd3)


! Cramer's rule
d1(k)  = detA1 * inv_detA
d1(kl) = detA2 * inv_detA
d1(km) = detA3 * inv_detA
d1(kr) = detA4 * inv_detA
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!$acc loop reduction(+:res1)
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
end subroutine pcr_rb_esa

!********************************************************************************
! pcr for vector (Aurora and GPU)
subroutine pcr_j_esa (sz, idx, g, pn, s, x, msk, rhs, a, c, d, a1, c1, d1, src, wrk, omg, res, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs, src, wrk
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, p, sq, s
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
!$acc loop independent collapse(3)
#else
!$OMP PARALLEL
! auroraはここにcollapseを入れると完全に並列化しない
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
do j=jst, jed
do i=ist, ied
do k=kst, ked
  src(k,i,j) = ( ( x(k, i  , j-1)        &
             +     x(k, i  , j+1)        &
             +     x(k, i-1, j  )        &
             +     x(k, i+1, j  ) - rhs(k, i, j) ) * r ) &
             *   msk(k, i, j)
end do
end do
end do ! 6 flops
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
#endif



#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) gang reduction(+:res1) &
!$acc& private(a, c, d, a1, c1, d1) &
!$acc& private(ap, cp, e, sq, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2)
#else
!$OMP DO SCHEDULE(static) Collapse(2) reduction(+:res1) &
!$OMP private(ap, cp, e, sq, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, cc1, cc2, f1, f2) &
!$OMP private(a1, c1, d1) &
!$OMP firstprivate(a, c, d)
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
do k = kst, ked
  d(k) = src(k, i, j)
end do ! 6 flops

! BC  6 flops
d(kst) = ( d(kst) + x(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked) = ( d(ked) + x(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
!$acc loop seq
do p=1, pn-1
sq = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
!$acc loop independent
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-sq) - cp * a(k+sq) )
a1(k) =   -e * ap * a(k-sq)
c1(k) =   -e * cp * c(k+sq)
d1(k) =   e * ( d(k) - ap * d(k-sq) - cp * d(k+sq) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = a1(k)
c(k) = c1(k)
d(k) = d1(k)
end do

end do ! p反復


! 最終段の反転 512のとき pn=9, sq=256
sq = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
!$acc loop independent
do k = kst, kst+sq-1 ! 2, 2+256-1=257
cc1 = c(k)
aa2 = a(k+sq)
f1  = d(k)
f2  = d(k+sq)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k   ) = dd1
d1(k+sq) = dd2
end do


! Relaxation
!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
!$acc loop reduction(+:res1)
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k) - pp ) * omg * msk(k, i, j)
wrk(k, i, j) = pp + dp
res1 = res1 + dp*dp
end do

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO NOWAIT
#endif


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(3)
#else
! auroraはここにcollapseを入れると完全に並列化しない
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
do j=jst, jed
do i=ist, ied
do k=kst, ked
x(k,i,j) = wrk(k,i,j)
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
end subroutine pcr_j_esa
