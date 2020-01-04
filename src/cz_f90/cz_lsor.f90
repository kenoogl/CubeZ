!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################

!> ********************************************************************
!! @brief 2x2の行列反転
!! @param [in,out] d    RHS vector(d) -> 解ベクトル(x) (in-place)
!! @param [in]     a    係数
!! @param [in]     c    係数
!! @note Ax = d    9 fp
!!       A=|1  c1|   x={x1, x2}, d={d1, d2}
!!         |a2  1|
!<
subroutine matx2(d, a, c)
implicit none
double precision, dimension(3)       ::  d, a, c
double precision                     ::  a2, c1, j, d1, d2
  a2 = a(2)
  c1 = c(1)
  j = 1.0 / (1.0 - a2 * c1)
  d1 = d(1)
  d2 = d(2)
  d(1) = (d1 - c1 * d2) * j
  d(2) = (d2 - a2 * d1) * j

return
end subroutine matx2


!> ********************************************************************
!! @brief 3x3の行列反転
!! @param [in,out] d    RHS vector(d) -> 解ベクトル(x) (in-place)
!! @param [in]     a    係数
!! @param [in]     c    係数
!! @note Ax = d    25 fp
!!       A=|  1 c1  0 |
!!         | a2  1 c2 |
!!         |  0 a3  1 |
!<
subroutine matx3(d, a, c)
implicit none
double precision, dimension(3)       ::  d, a, c
double precision                     ::  j, d1, d2, d3, a2, a3, c1, c2

  a2 = a(2)
  a3 = a(3)
  c1 = c(1)
  c2 = c(2)
  d1 = d(1)
  d2 = d(2)
  d3 = d(3)
  j = 1.0 / (1.0 - c2 * a3 - c1 * a2)
  d(1) = ( d1 * (3.0-c2*a3) - c1*d2 ) * j
  d(2) = (1.0 - d1*a2 + 2.0*d2 - c2*d3) * j
  d(3) = (1.0 + 2.0*d3 - a3*d2 - a2*c1) * j

return
end subroutine matx3



!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_pcr_kij (sz, idx, g, pn, x, a, c, d, a1, c1, d1, msk, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s, p, pn
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  a, c, d, a1, c1, d1
real                                                   ::  r, ap, cp, e, omg, pp, dp
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble(  &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)*( 6.0       &  ! Source
+ pn * 14.0 &  ! PCR
+ 6.0 )     &  ! Relaxation
+ 6.0                   &  ! BC
) )


!$OMP PARALLEL

! Reflesh coef. due to override

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst+1, ked
a(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst, ked-1
c(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
a(kst,i,j) = 0.0
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
c(ked,i,j) = 0.0
end do
end do
!$OMP END DO


res = 0.0


!$OMP DO SCHEDULE(static) collapse(2) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP reduction(+:res)
do j=jst, jed
do i=ist, ied

! Source
!dir$ vector aligned
!dir$ simd
do k = kst, ked
d(k, i, j) = ( ( x(k, i  , j-1)        &
+     x(k, i  , j+1)        &
+     x(k, i-1, j  )        &
+     x(k, i+1, j  ) ) * r + rhs(k, i, j) ) &
*   msk(k, i, j)
end do

! BC
d(kst, i, j) = ( d(kst, i, j) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked, i, j) = ( d(ked, i, j) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR
do p=1, pn
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
ap = a(k,i,j)
cp = c(k,i,j)
e = 1.0 / ( 1.0 - ap * c(kl,i,j) - cp * a(kr,i,j) )
a1(k,i,j) =  -e * ap * a(kl,i,j)
c1(k,i,j) =  -e * cp * c(kr,i,j)
d1(k,i,j) =   e * ( d(k,i,j) - ap * d(kl,i,j) - cp * d(kr,i,j) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k,i,j) = a1(k,i,j)
c(k,i,j) = c1(k,i,j)
d(k,i,j) = d1(k,i,j)
end do

end do

! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked
pp =   x(k, i, j)
dp = ( d(k, i, j) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res = res + real(dp*dp, kind=8)
end do

end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine lsor_pcr_kij




!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note lsor_pcr_kij()からの変更 最終段を直接反転
!<
subroutine lsor_pcr_kij2 (sz, idx, g, pn, x, a, c, d, a1, c1, d1, msk, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s, p, pn
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  a, c, d, a1, c1, d1
real                                                   ::  r, ap, cp, e, omg, pp, dp
double precision, dimension(3)                         ::  aa, cc, dd
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64
!DIR$ ATTRIBUTES FORCEINLINE::matx2, matx3

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


!$OMP PARALLEL

! Reflesh coef. due to override

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst+1, ked
a(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst, ked-1
c(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
a(kst,i,j) = 0.0
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
c(ked,i,j) = 0.0
end do
end do
!$OMP END DO


res = 0.0



!$OMP DO SCHEDULE(static) collapse(2) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp, aa, cc, dd) &
!$OMP reduction(+:res)
do j=jst, jed
do i=ist, ied

! Source
!dir$ vector aligned
!dir$ simd
do k = kst, ked
d(k, i, j) = ( ( x(k, i  , j-1)        &
+     x(k, i  , j+1)        &
+     x(k, i-1, j  )        &
+     x(k, i+1, j  ) ) * r + rhs(k, i, j) ) &
*   msk(k, i, j)
end do

! BC
d(kst, i, j) = ( d(kst, i, j) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked, i, j) = ( d(ked, i, j) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
ap = a(k,i,j)
cp = c(k,i,j)
e = 1.0 / ( 1.0 - ap * c(kl,i,j) - cp * a(kr,i,j) )
a1(k,i,j) =  -e * ap * a(kl,i,j)
c1(k,i,j) =  -e * cp * c(kr,i,j)
d1(k,i,j) =   e * ( d(k,i,j) - ap * d(kl,i,j) - cp * d(kr,i,j) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k,i,j) = a1(k,i,j)
c(k,i,j) = c1(k,i,j)
d(k,i,j) = d1(k,i,j)
end do

end do ! p反復


! 最終段の反転
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)

if (k<kst+s) then ! 2 eqations
cc(1) = real( c(k ,i,j), kind=8)
aa(2) = real( a(kr,i,j), kind=8)
dd(1) = real( d(k ,i,j), kind=8)
dd(2) = real( d(kr,i,j), kind=8)
call matx2(dd, aa, cc)
d1(k ,i,j) = real( dd(1), kind=4)
d1(kr,i,j) = real( dd(2), kind=4)
else if (k<=ked-s) then ! 3 equations
cc(1) = real( c(kr,i,j), kind=8)
aa(2) = real( a(k ,i,j), kind=8)
cc(2) = real( c(k ,i,j), kind=8)
aa(3) = real( a(kl,i,j), kind=8)
dd(1) = real( d(kr,i,j), kind=8)
dd(2) = real( d(k ,i,j), kind=8)
dd(3) = real( d(kl,i,j), kind=8)
call matx3(dd, aa, cc)
d1(kr,i,j) = real( dd(1), kind=4)
d1(k ,i,j) = real( dd(2), kind=4)
d1(kl,i,j) = real( dd(3), kind=4)
else ! 2 equations
cc(1) = real( c(kl,i,j), kind=8)
aa(2) = real( a(k ,i,j), kind=8)
dd(1) = real( d(kl,i,j), kind=8)
dd(2) = real( d(k ,i,j), kind=8)
call matx2(dd, aa, cc)
d1(kl,i,j) = real( dd(1), kind=4)
d1(k ,i,j) = real( dd(2), kind=4)
endif
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k, i, j) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res = res + real(dp*dp, kind=8)
end do

end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine lsor_pcr_kij2



!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note lsor_pcr_kij2()からの変更 分割
!<
subroutine lsor_pcr_kij3 (sz, idx, g, pn, x, a, c, d, a1, c1, d1, msk, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s, p, pn, t
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  a, c, d, a1, c1, d1
real                                                   ::  r, ap, cp, e, omg, pp, dp
double precision, dimension(3)                         ::  aa, cc, dd
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64
!DIR$ ATTRIBUTES FORCEINLINE::matx2, matx3

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
t = 2**(pn-1)

flop = flop + dble(  &
(jed-jst+1)*(ied-ist+1)* ( &
(ked-kst+1)* 6.0        &  ! Source
+ (ked-kst+1)*(pn-1)*14.0 &  ! PCR
+ 2*t*9.0                 &
+ (ked-kst-2*t+1)*25.0    &
+ (ked-kst+1)*6.0         &  ! Relaxation
+ 6.0 )                 &  ! BC
)


!$OMP PARALLEL

! Refresh coef. due to override

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst+1, ked
a(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst, ked-1
c(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
a(kst,i,j) = 0.0
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
c(ked,i,j) = 0.0
end do
end do
!$OMP END DO


res = 0.0



!$OMP DO SCHEDULE(static) collapse(2) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp, aa, cc, dd) &
!$OMP reduction(+:res)
do j=jst, jed
do i=ist, ied

! Source
!dir$ vector aligned
!dir$ simd
do k = kst, ked ! 6 fp
d(k, i, j) = ( ( x(k, i  , j-1)        &
+     x(k, i  , j+1)        &
+     x(k, i-1, j  )        &
+     x(k, i+1, j  ) ) * r + rhs(k, i, j) ) &
*   msk(k, i, j)
end do

! BC
d(kst, i, j) = ( d(kst, i, j) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked, i, j) = ( d(ked, i, j) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked ! 14 fp
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
ap = a(k,i,j)
cp = c(k,i,j)
e = 1.0 / ( 1.0 - ap * c(kl,i,j) - cp * a(kr,i,j) )
a1(k,i,j) =  -e * ap * a(kl,i,j)
c1(k,i,j) =  -e * cp * c(kr,i,j)
d1(k,i,j) =   e * ( d(k,i,j) - ap * d(kl,i,j) - cp * d(kr,i,j) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k,i,j) = a1(k,i,j)
c(k,i,j) = c1(k,i,j)
d(k,i,j) = d1(k,i,j)
end do

end do ! p反復


! 最終段の反転
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
do k = kst, kst+s-1
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc(1) = real( c(k ,i,j), kind=8)
aa(2) = real( a(kr,i,j), kind=8)
dd(1) = real( d(k ,i,j), kind=8)
dd(2) = real( d(kr,i,j), kind=8)
call matx2(dd, aa, cc) ! 9 fp
d1(k ,i,j) = real( dd(1), kind=4)
d1(kr,i,j) = real( dd(2), kind=4)
end do


!dir$ vector aligned
!dir$ simd
do k = kst+s, ked-s
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc(1) = real( c(kr,i,j), kind=8)
aa(2) = real( a(k ,i,j), kind=8)
cc(2) = real( c(k ,i,j), kind=8)
aa(3) = real( a(kl,i,j), kind=8)
dd(1) = real( d(kr,i,j), kind=8)
dd(2) = real( d(k ,i,j), kind=8)
dd(3) = real( d(kl,i,j), kind=8)
call matx3(dd, aa, cc) ! 25 fp
d1(kr,i,j) = real( dd(1), kind=4)
d1(k ,i,j) = real( dd(2), kind=4)
d1(kl,i,j) = real( dd(3), kind=4)
end do


!dir$ vector aligned
!dir$ simd
do k = ked-s+1, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc(1) = real( c(kl,i,j), kind=8)
aa(2) = real( a(k ,i,j), kind=8)
dd(1) = real( d(kl,i,j), kind=8)
dd(2) = real( d(k ,i,j), kind=8)
call matx2(dd, aa, cc) ! 9 fp
d1(kl,i,j) = real( dd(1), kind=4)
d1(k ,i,j) = real( dd(2), kind=4)
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked ! 6 fp
pp =   x(k, i, j)
dp = ( d1(k, i, j) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res = res + real(dp*dp, kind=8)
end do

end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine lsor_pcr_kij3



!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note lsor_pcr_kij4()からの変更 matx2, matx3を手動展開
!<
subroutine lsor_pcr_kij4 (sz, idx, g, pn, x, a, c, d, a1, c1, d1, msk, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s, p, pn
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  a, c, d, a1, c1, d1
real                                                   ::  r, ap, cp, e, omg, pp, dp
double precision                                       ::  jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64

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


!$OMP PARALLEL

! Reflesh coef. due to override

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst+1, ked
a(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst, ked-1
c(k,i,j) = -r
end do
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
a(kst,i,j) = 0.0
end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
c(ked,i,j) = 0.0
end do
end do
!$OMP END DO


res = 0.0



!$OMP DO SCHEDULE(static) collapse(2) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP reduction(+:res)
do j=jst, jed
do i=ist, ied

! Source
!dir$ vector aligned
!dir$ simd
do k = kst, ked
d(k, i, j) = ( ( x(k, i  , j-1)        &
+     x(k, i  , j+1)        &
+     x(k, i-1, j  )        &
+     x(k, i+1, j  ) ) * r + rhs(k, i, j) ) &
*   msk(k, i, j)
end do

! BC
d(kst, i, j) = ( d(kst, i, j) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked, i, j) = ( d(ked, i, j) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
do k = kst, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
ap = a(k,i,j)
cp = c(k,i,j)
e = 1.0 / ( 1.0 - ap * c(kl,i,j) - cp * a(kr,i,j) )
a1(k,i,j) =  -e * ap * a(kl,i,j)
c1(k,i,j) =  -e * cp * c(kr,i,j)
d1(k,i,j) =   e * ( d(k,i,j) - ap * d(kl,i,j) - cp * d(kr,i,j) )
end do

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k,i,j) = a1(k,i,j)
c(k,i,j) = c1(k,i,j)
d(k,i,j) = d1(k,i,j)
end do

end do ! p反復


! 最終段の反転
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
do k = kst, kst+s-1
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = real( c(k ,i,j), kind=8)
aa2 = real( a(kr,i,j), kind=8)
f1  = real( d(k ,i,j), kind=8)
f2  = real( d(kr,i,j), kind=8)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k ,i,j) = real( dd1, kind=4)
d1(kr,i,j) = real( dd2, kind=4)
end do


!dir$ vector aligned
!dir$ simd
do k = kst+s, ked-s
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = real( c(kr,i,j), kind=8)
aa2 = real( a(k ,i,j), kind=8)
cc2 = real( c(k ,i,j), kind=8)
aa3 = real( a(kl,i,j), kind=8)
f1  = real( d(kr,i,j), kind=8)
f2  = real( d(k ,i,j), kind=8)
f3  = real( d(kl,i,j), kind=8)
jj = 1.0 / (1.0 - cc2 * aa3 - cc1 * aa2)
dd1 = ( f1 * (3.0-cc2*aa3) - cc1*f2 ) * jj
dd2 = (1.0 - f1*aa2 + 2.0*f2 - cc2*f3) * jj
dd3 = (1.0 + 2.0*f3 - aa3*f2 - aa2*cc1) * jj
d1(kr,i,j) = real( dd1, kind=4)
d1(k ,i,j) = real( dd2, kind=4)
d1(kl,i,j) = real( dd3, kind=4)
end do


!dir$ vector aligned
!dir$ simd
do k = ked-s+1, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = real( c(kl,i,j), kind=8)
aa2 = real( a(k ,i,j), kind=8)
f1  = real( d(kl,i,j), kind=8)
f2  = real( d(k ,i,j), kind=8)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(kl,i,j) = real( dd1, kind=4)
d1(k ,i,j) = real( dd2, kind=4)
end do


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked
pp =   x(k, i, j)
dp = ( d1(k, i, j) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp
res = res + real(dp*dp, kind=8)
end do

end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine lsor_pcr_kij4



!********************************************************************************
subroutine lsor_pcr_kij5 (sz, idx, g, pn, x, msk, rhs, omg, res, flop)
implicit none
! arguments
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                      ::  i, j, k, kl, kr, s, p
integer                      ::  ist, ied, jst, jed, kst, ked
real, dimension(1-g:sz(3)+g) ::  a, c, d, a1, c1, d1
real                         ::  r, ap, cp, e, pp, dp
double precision             ::  jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64

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


!$OMP PARALLEL

res = 0.0

!$OMP DO SCHEDULE(static) collapse(2) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, a1, c1, d1) &
!$OMP reduction(+:res)
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
+     x(k, i+1, j  ) ) * r + rhs(k, i, j) ) &
*   msk(k, i, j)
end do

! BC
d(kst) = ( d(kst) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked) = ( d(ked) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


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
cc1 = real( c(k ), kind=8)
aa2 = real( a(kr), kind=8)
f1  = real( d(k ), kind=8)
f2  = real( d(kr), kind=8)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(k ) = real( dd1, kind=4)
d1(kr) = real( dd2, kind=4)
end do


!dir$ vector aligned
!dir$ simd
do k = kst+s, ked-s
kl  = max(k-s, kst-1)
kr  = min(k+s, ked+1)
cc1 = real( c(kr), kind=8)
aa2 = real( a(k ), kind=8)
cc2 = real( c(k ), kind=8)
aa3 = real( a(kl), kind=8)
f1  = real( d(kr), kind=8)
f2  = real( d(k ), kind=8)
f3  = real( d(kl), kind=8)
jj = 1.0 / (1.0 - cc2 * aa3 - cc1 * aa2)
dd1 = ( f1 * (3.0-cc2*aa3) - cc1*f2 ) * jj
dd2 = (1.0 - f1*aa2 + 2.0*f2 - cc2*f3) * jj
dd3 = (1.0 + 2.0*f3 - aa3*f2 - aa2*cc1) * jj
d1(kr) = real( dd1, kind=4)
d1(k ) = real( dd2, kind=4)
d1(kl) = real( dd3, kind=4)
end do


!dir$ vector aligned
!dir$ simd
do k = ked-s+1, ked
kl = max(k-s, kst-1)
kr = min(k+s, ked+1)
cc1 = real( c(kl), kind=8)
aa2 = real( a(k ), kind=8)
f1  = real( d(kl), kind=8)
f2  = real( d(k ), kind=8)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
d1(kl) = real( dd1, kind=4)
d1(k ) = real( dd2, kind=4)
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
end subroutine lsor_pcr_kij5



!********************************************************************************
subroutine lsor_pcr_kij6 (sz, idx, g, pn, x, msk, rhs, omg, res, flop)
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
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64

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


!$OMP PARALLEL

res = 0.0

!$OMP DO SCHEDULE(static) collapse(2) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, a1, c1, d1) &
!$OMP reduction(+:res)
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
+     x(k, i+1, j  ) ) * r + rhs(k, i, j) ) &
*   msk(k, i, j)
end do

! BC
d(kst) = ( d(kst) + rhs(kst-1, i, j) * r ) * msk(kst, i, j)
d(ked) = ( d(ked) + rhs(ked+1, i, j) * r ) * msk(ked, i, j)


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
end subroutine lsor_pcr_kij6

