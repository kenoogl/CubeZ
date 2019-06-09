!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################


!> ********************************************************************
!! @brief point SOR法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     X,Y,Z  座標
!! @param [in]     omg  加速係数
!! @param [in]     b    RHS vector
!! @param [out]    res  residual
!! @param [in,out] flop flop count
!<
subroutine psor_maf (p, sz, idx, g, X, Y, Z, omg, b, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  res
double precision                                       ::  flop
real                                                   ::  omg, dd, dp, pp, bb, pn, rp
real                                                   ::  GX, EY, TZ, YJA, YJAI
real                                                   ::  XG, YE, ZT, XGG, YEE, ZTT
real                                                   ::  C1, C2, C3, C7, C8, C9
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z


res = 0.0

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)


flop = flop + 76.0d0      &
* dble(ied-ist+1) &
* dble(jed-jst+1) &
* dble(ked-kst+1)


!$OMP PARALLEL DO Collapse(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(rp, pp, bb, dd, dp, pn) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
do j = jst, jed
do i = ist, ied
do k = kst, ked
bb = b(k,i,j)
pp = p(k,i,j)

XG = 0.5 * (X(i+1) - X(i-1))
YE = 0.5 * (Y(j+1) - Y(j-1))
ZT = 0.5 * (Z(k+1) - Z(k-1)) ! 6

XGG= X(i+1) - 2.0*X(i) + X(i-1)
YEE= Y(j+1) - 2.0*Y(j) + Y(j-1)
ZTT= Z(k+1) - 2.0*Z(k) + Z(k-1) ! 9

! Jacobian
YJA  = XG * YE * ZT
YJAI = 1.0 / YJA    ! 8 + 2

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


dd = C1 + C2 + C3
rp = 0.5 *((C1 * (P(k  , i+1, j  ) + P(k  , i-1, j  ))        &
   +        C2 * (P(k  , i  , j+1) + P(k  , i  , j-1))        &
   +        C3 * (P(k+1, i  , j  ) + P(k-1, i  , j  ))        &
   +        C7 * (P(k  , i+1, j  ) - P(k  , i-1, j  )) * 0.5  &
   +        C8 * (P(k  , i  , j+1) - P(k  , i  , j-1)) * 0.5  &
   +        C9 * (P(k+1, i  , j  ) - P(k-1, i  , j  )) * 0.5) &
   + bb)
dp = ( rp / dd - pp ) * omg
pn = pp + dp
P(k,i,j) = pn
res = res + real(dp * dp, kind=8) ! 30
enddo
enddo
enddo
!$OMP END PARALLEL DO

return
end subroutine psor_maf


!> **********************************************************************
!! @brief 緩和Jacobi法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     idx         インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     X,Y,Z  座標
!! @param [in]     omg  加速係数
!! @param [in]     b    RHS vector
!! @param [in,out] res  residual
!! @param [out]    wk2  ワーク用配列
!! @param [in,out] flop flop count
!<
subroutine jacobi_maf (p, sz, idx, g, X, Y, Z, omg, b, res, wk2, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  res
double precision                                       ::  flop
real                                                   ::  omg, dd, dp, pp, bb, pn, rp
real                                                   ::  GX, EY, TZ, YJA, YJAI
real                                                   ::  XG, YE, ZT, XGG, YEE, ZTT
real                                                   ::  C1, C2, C3, C7, C8, C9
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b, wk2
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z
!dir$ assume_aligned p:64, b:64, wk2:64, X:64, Y:64, Z:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res = 0.0


flop = flop + 76.0d0  &
* dble(ied-ist+1) &
* dble(jed-jst+1) &
* dble(ked-kst+1)

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(rp, pp, bb, dd, dp, pn) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
!$OMP DO SCHEDULE(static) Collapse(2)
do j = jst, jed
do i = ist, ied
do k = kst, ked
bb = b(k,i,j)
pp = p(k,i,j)

XG = 0.5 * (X(i+1) - X(i-1))
YE = 0.5 * (Y(j+1) - Y(j-1))
ZT = 0.5 * (Z(k+1) - Z(k-1)) ! 6

XGG= X(i+1) - 2.0*X(i) + X(i-1)
YEE= Y(j+1) - 2.0*Y(j) + Y(j-1)
ZTT= Z(k+1) - 2.0*Z(k) + Z(k-1) ! 9

! Jacobian
YJA  = XG * YE * ZT
YJAI = 1.0 / YJA    ! 8 + 2

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


dd = C1 + C2 + C3
rp = 0.5 *((C1 * (P(k  , i+1, j  ) + P(k  , i-1, j  ))        &
+        C2 * (P(k  , i  , j+1) + P(k  , i  , j-1))        &
+        C3 * (P(k+1, i  , j  ) + P(k-1, i  , j  ))        &
+        C7 * (P(k  , i+1, j  ) - P(k  , i-1, j  )) * 0.5  &
+        C8 * (P(k  , i  , j+1) - P(k  , i  , j-1)) * 0.5  &
+        C9 * (P(k+1, i  , j  ) - P(k-1, i  , j  )) * 0.5) &
+ bb)
dp = ( rp / dd - pp ) * omg
pn = pp + dp
wk2(k,i,j) = pn
res = res + real(dp * dp, kind=8) ! 30
enddo
enddo
enddo
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
end subroutine jacobi_maf


!> ********************************************************************
!! @brief 2-colored SOR法 stride memory access
!! @param [in,out] p     圧力
!! @param [in]     sz    配列長
!! @param [in]     idx   インデクス範囲
!! @param [in]     g     ガイドセル長
!! @param [in]     X,Y,Z  座標
!! @param [in]     ofst  開始点オフセット
!! @param [in]     color グループ番号
!! @param [in]     omg   加速係数
!! @param [in]     b     RHS vector
!! @param [out]    res  residual
!! @param [in,out] flop  浮動小数演算数
!! @note resは積算
!<
subroutine psor2sma_core_maf (p, sz, idx, g, X, Y, Z, ofst, color, omg, b, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
double precision                                       ::  res
real                                                   ::  omg, dd, dp, pp, bb, pn, rp
real                                                   ::  GX, EY, TZ, YJA, YJAI
real                                                   ::  XG, YE, ZT, XGG, YEE, ZTT
real                                                   ::  C1, C2, C3, C7, C8, C9
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
integer                                                ::  kp, color, ofst
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z
!dir$ assume_aligned p:64, b:64, X:64, Y:64, Z:64

kp = ofst+color

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 76.0d0*0.5d0  &
* dble(ied-ist+1) &
* dble(jed-jst+1) &
* dble(ked-kst+1)


!$OMP PARALLEL DO Collapse(2) &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(rp, pp, bb, dd, dp, pn) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
do j=jst,jed
do i=ist,ied

!dir$ vector aligned
!dir$ simd
do k=kst+mod(i+j+kp,2), ked, 2
pp = p(k,i,j)
bb = b(k,i,j)

XG = 0.5 * (X(i+1) - X(i-1))
YE = 0.5 * (Y(j+1) - Y(j-1))
ZT = 0.5 * (Z(k+1) - Z(k-1)) ! 6

XGG= X(i+1) - 2.0*X(i) + X(i-1)
YEE= Y(j+1) - 2.0*Y(j) + Y(j-1)
ZTT= Z(k+1) - 2.0*Z(k) + Z(k-1) ! 9

! Jacobian
YJA  = XG * YE * ZT
YJAI = 1.0 / YJA    ! 8 + 2

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


dd = C1 + C2 + C3
rp = 0.5 *((C1 * (P(k  , i+1, j  ) + P(k  , i-1, j  ))        &
+        C2 * (P(k  , i  , j+1) + P(k  , i  , j-1))        &
+        C3 * (P(k+1, i  , j  ) + P(k-1, i  , j  ))        &
+        C7 * (P(k  , i+1, j  ) - P(k  , i-1, j  )) * 0.5  &
+        C8 * (P(k  , i  , j+1) - P(k  , i  , j-1)) * 0.5  &
+        C9 * (P(k+1, i  , j  ) - P(k-1, i  , j  )) * 0.5) &
+ bb)
dp = ( rp / dd - pp ) * omg
pn = pp + dp
P(k,i,j) = pn
res = res + real(dp * dp, kind=8) ! 30
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine psor2sma_core_maf


!********************************************************************************
subroutine lsor_pcr_kij7_maf (sz, idx, g, pn, ofst, color, x, msk, rhs, XX, YY, ZZ, omg, res, flop)
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
real, dimension(1-g:sz(3)+g)             ::  a, c, d, a1, c1, d1
real                                     ::  r, ap, cp, e, pp, dp
real                                     ::  jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3
real, dimension(-1:sz(1)+2)                            ::  XX
real, dimension(-1:sz(2)+2)                            ::  YY
real, dimension(-1:sz(3)+2)                            ::  ZZ
!dir$ assume_aligned x:64, msk:64, rhs:64, a:64, c:64, d:64, a1:64, c1:64, d1:64, XX:64, YY:64, ZZ:64

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

!$OMP PARALLEL

res = 0.0

!$OMP DO SCHEDULE(static) &
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, dd3, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, a1, c1, d1) &
!$OMP reduction(+:res)
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
end subroutine lsor_pcr_kij7_maf

