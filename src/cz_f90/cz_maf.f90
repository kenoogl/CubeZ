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
real                                                   ::  C1, C2, C3, C7, C8, C9, res1
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z

res1 = 0.0

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 66.0d0      &
* dble(ied-ist+1) &
* dble(jed-jst+1) &
* dble(ked-kst+1)


!$OMP PARALLEL DO Collapse(2) &
!$OMP REDUCTION(+:res1) &
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


dd = 2.0 * (C1 + C2 + C3)
rp = (C1 + 0.5 * C7) * P(k  , i+1, j  ) &
   + (C1 - 0.5 * C7) * P(k  , i-1, j  ) &
   + (C2 + 0.5 * C8) * P(k  , i  , j+1) &
   + (C2 - 0.5 * C8) * P(k  , i  , j-1) &
   + (C3 + 0.5 * C9) * P(k+1, i  , j  ) &
   + (C3 - 0.5 * C9) * P(k-1, i  , j  ) &
   + bb
dp = ( rp / dd - pp ) * omg
pn = pp + dp
P(k,i,j) = pn
res1 = res1 + dp * dp ! 30
enddo
enddo
enddo
!$OMP END PARALLEL DO

res = res + real(res1, kind=8)

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
!! @param [in]     tmp  ワーク
!! @param [in,out] flop flop count
!<
subroutine jacobi_maf (p, sz, idx, g, X, Y, Z, omg, b, res, wk2, tmp, flop)
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
real                                                   ::  C1, C2, C3, C7, C8, C9, res1
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b, wk2
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z, tmp

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0

#ifdef _SVR
tmp = 0.0
#endif

flop = flop + 66.0d0  &
* dble(ied-ist+1) &
* dble(jed-jst+1) &
* dble(ked-kst+1)


#ifdef _OPENACC
!$acc kernels
!$acc loop collapse(3) reduction(+:res1)
#else
!$OMP PARALLEL &
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP PRIVATE(rp, pp, bb, dd, dp, pn) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
#ifdef __NEC__
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(static) COLLAPSE(2)
#endif
#endif
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

dd = 2.0 * (C1 + C2 + C3)
rp = (C1 + 0.5 * C7) * P(k  , i+1, j  ) &
   + (C1 - 0.5 * C7) * P(k  , i-1, j  ) &
   + (C2 + 0.5 * C8) * P(k  , i  , j+1) &
   + (C2 - 0.5 * C8) * P(k  , i  , j-1) &
   + (C3 + 0.5 * C9) * P(k+1, i  , j  ) &
   + (C3 - 0.5 * C9) * P(k-1, i  , j  ) &
   + bb
dp = ( rp / dd - pp ) * omg
pn = pp + dp
wk2(k,i,j) = pn

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

enddo
enddo
enddo
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


#ifdef _SVR
do k = kst, ked
  res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

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
!! @param [in]     tmp  ワーク
!! @param [in,out] flop  浮動小数演算数
!! @note resは積算
!<
subroutine psor2sma_core_maf (p, sz, idx, g, X, Y, Z, ofst, color, omg, b, res, tmp, flop)
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
real                                                   ::  C1, C2, C3, C7, C8, C9, res1
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  p, b
integer                                                ::  kp, color, ofst
real, dimension(-1:sz(1)+2)                            ::  X
real, dimension(-1:sz(2)+2)                            ::  Y
real, dimension(-1:sz(3)+2)                            ::  Z, tmp

kp = ofst+color
res1 = 0.0

#ifdef _SVR
tmp = 0.0
#endif

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + 66.0d0*0.5d0  &
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
#ifdef __NEC__
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) & ! ここはcollapseを入れた方がよい
#else
!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
#endif
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP PRIVATE(rp, pp, bb, dd, dp, pn) &
!$OMP PRIVATE(XG, YE, ZT, XGG, YEE, ZTT) &
!$OMP PRIVATE(GX, EY, TZ, YJA, YJAI) &
!$OMP PRIVATE(C1, C2, C3, C7, C8, C9)
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

dd = 2.0 * (C1 + C2 + C3)
rp = (C1 + 0.5 * C7) * P(k  , i+1, j  ) &
   + (C1 - 0.5 * C7) * P(k  , i-1, j  ) &
   + (C2 + 0.5 * C8) * P(k  , i  , j+1) &
   + (C2 - 0.5 * C8) * P(k  , i  , j-1) &
   + (C3 + 0.5 * C9) * P(k+1, i  , j  ) &
   + (C3 - 0.5 * C9) * P(k-1, i  , j  ) &
   + bb

dp = ( rp / dd - pp ) * omg
pn = pp + dp
P(k,i,j) = pn

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

end do
end do
end do
#ifdef _OPENACC
!$acc end kernels
#endif
!$OMP END PARALLEL DO

#ifdef _SVR
do k = kst, ked
res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

return
end subroutine psor2sma_core_maf


!********************************************************************************
subroutine pcr_rb_maf (sz, idx, g, pn, ofst, color, x, msk, rhs, XX, YY, ZZ, &
                       a, c, d, aw, cw, dw, omg, res, tmp, flop)
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
real, dimension(-1:sz(3)+2)              ::  a, c, d, aw, cw, dw
real                                     ::  ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3
real                                                   ::  C1, C2, C7, C8, GX, EY, TZ, ZTT
real, dimension(-1:sz(1)+2)                            ::  XX
real, dimension(-1:sz(2)+2)                            ::  YY
real, dimension(-1:sz(3)+2)                            ::  ZZ, tmp

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)


flop = flop + dble(           &
  (jed-jst+1)*(ied-ist+1)* (  &
    ( 24.0d0 +                & ! metrics
     + 3.0 * 2.0 + 12.0       & ! coef BC
    )                         &
   + (ked-kst+1)* (11.0+10.0) & ! coef + Source
   + (ked-kst-1)* 6.0d0       & ! coef
   + (ked-kst+1)*(pn-1)*16.0  & ! PCR
   + 2**(pn-1) * 11.0           & ! 2x2
   + (ked-kst+1)*6.0          & ! Relaxation
   )                          &
 ) * 0.5


ip = ofst + color
res1 = 0.0

#ifdef _SVR
tmp = 0.0
#endif


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) gang reduction(+:res1) &
!$acc& private(a, c, d, aw, cw, dw) &
!$acc& private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$acc& private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
#else
!$OMP PARALLEL &
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, aw, cw, dw) &
!$OMP private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
!$OMP DO SCHEDULE(static) collapse(2)
#endif
do j=jst, jed
do i=ist, ied
if(mod(i+j,2) /= color) cycle

! do i=ist+mod(j+ip,2), ied, 2

GX =  2.0 / (XX(i+1) - XX(i-1))
EY =  2.0 / (YY(j+1) - YY(j-1))
C1 =  GX * GX
C2 =  EY * EY
C7 = -(XX(i+1) - 2.0*XX(i) + XX(i-1)) * C1 * GX
C8 = -(YY(j+1) - 2.0*YY(j) + YY(j-1)) * C2 * EY
dd1= C1 + 0.5 * C7 ! R1
dd2= C1 - 0.5 * C7 ! R2
cc1= C2 + 0.5 * C8 ! R3
cc2= C2 - 0.5 * C8 ! R4      >>  24 flops

! Reflesh coef. due to override
do k=kst, ked
  f1 = ZZ(k+1)
  f2 = ZZ(k-1)
  TZ = 2.0 / (f1 - f2)
  ZTT= f1 - 2.0*ZZ(k) + f2
  f3 = TZ * TZ
  aw(k) = f3              ! C3
  cw(k) = -ZTT * f3 * TZ  ! C9
  dw(k) = 0.5 / (C1 + C2 + f3) ! 1/R7
end do   ! >>  11 flops


a(kst) = 0.0
c(kst) = -(aw(kst) + 0.5 * cw(kst)) * dw(kst) ! -R5/R7 = -(C3+0.5*C9) / R7      >> 3 flops

do k=kst+1, ked-1
  f1 = aw(k)  ! C3
  f2 = cw(k)  ! C9
  aa3 = dw(k)
  a(k) = -(f1 - 0.5 * f2) * aa3
  c(k) = -(f1 + 0.5 * f2) * aa3
end do    !  >> 6 flops

a(ked) = -(aw(ked) - 0.5 * cw(ked)) * dw(ked) ! -R6/R7 = -(C3-0.5*C9) / R7     >> 3 flops
c(ked) = 0.0

! Source
do k = kst, ked
  d(k) = (                         &
            dd1 * x(k, i+1, j  )   &
          + dd2 * x(k, i-1, j  )   &
          + cc1 * x(k, i  , j+1)   &
          + cc2 * x(k, i  , j-1)   &
          - rhs(k, i, j)           &
        ) * dw(k) * msk(k, i, j)
end do   !  >>  10 flops
! ここまで、dd1, dd2, cc1, cc2は再利用しない
! a, c, dを計算したので、aw, cw, dwは再利用可能

! BC   >>  12 flops
d(kst) = ( d(kst) + (aw(kst) - 0.5 * cw(kst)) * dw(kst) * x(kst-1, i, j) ) * msk(kst, i, j)
d(ked) = ( d(ked) + (aw(ked) + 0.5 * cw(ked)) * dw(ked) * x(ked+1, i, j) ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
!$acc loop seq
do p=1, pn-1
s = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
do k = kst, ked
  kl = max(k-s, kst-1)
  kr = min(k+s, ked+1)
  ap = a(k)
  cp = c(k)
  e = 1.0 / ( 1.0 - ap * c(kl) - cp * a(kr) )
  aw(k) =  -e * ap * a(kl)
  cw(k) =  -e * cp * c(kr)
  dw(k) =   e * ( d(k) - ap * d(kl) - cp * d(kr) )
end do   !  >> 16 flops

!dir$ vector aligned
!dir$ simd
do k = kst, ked
  a(k) = aw(k)
  c(k) = cw(k)
  d(k) = dw(k)
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
  dw(k ) = dd1
  dw(kr) = dd2
end do  ! >>  11 flops



! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!$acc loop reduction(+:res1)
do k = kst, ked
  pp =   x(k, i, j)
  dp = ( dw(k) - pp ) * omg * msk(k, i, j)
  x(k, i, j) = pp + dp

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

end do  !  >> 6 flops

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef _SVR
do k = kst, ked
res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

return
end subroutine pcr_rb_maf


!********************************************************************************
subroutine pcr_maf (sz, idx, g, pn, x, msk, rhs, XX, YY, ZZ, &
                    a, c, d, aw, cw, dw, omg, res, tmp, flop)
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
real, dimension(-1:sz(3)+2)              ::  a, c, d, aw, cw, dw
real                                     ::  ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3
real                                                   ::  C1, C2, C7, C8, GX, EY, TZ, ZTT
real, dimension(-1:sz(1)+2)                            ::  XX
real, dimension(-1:sz(2)+2)                            ::  YY
real, dimension(-1:sz(3)+2)                            ::  ZZ, tmp

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0

s = 2**(pn-1)

flop = flop + dble(           &
(jed-jst+1)*(ied-ist+1)* (  &
( 24.0d0 +                & ! metrics
+ 3.0 * 2.0 + 12.0       & ! coef BC
)                         &
+ (ked-kst+1)* (11.0+10.0) & ! coef + Source
+ (ked-kst-1)* 6.0d0       & ! coef
+ (ked-kst+1)*(pn-1)*16.0  & ! PCR
+ s * 11.0           & ! 2x2
+ (ked-kst+1)*6.0          & ! Relaxation
)                          &
)

#ifdef _SVR
tmp = 0.0
#endif


#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) reduction(+:res1) &
!$acc& private(a, c, d, aw, cw, dw) &
!$acc& private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$acc& private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
#else
!$OMP PARALLEL &
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(a, c, d, aw, cw, dw) &
!$OMP private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
!$OMP DO SCHEDULE(static) collapse(2)
#endif
do j=jst, jed
do i=ist, ied

GX =  2.0 / (XX(i+1) - XX(i-1))
EY =  2.0 / (YY(j+1) - YY(j-1))
C1 =  GX * GX
C2 =  EY * EY
C7 = -(XX(i+1) - 2.0*XX(i) + XX(i-1)) * C1 * GX
C8 = -(YY(j+1) - 2.0*YY(j) + YY(j-1)) * C2 * EY
dd1= C1 + 0.5 * C7 ! R1
dd2= C1 - 0.5 * C7 ! R2
cc1= C2 + 0.5 * C8 ! R3
cc2= C2 - 0.5 * C8 ! R4      >>  24 flops

! Reflesh coef. due to override
do k=kst, ked
f1 = ZZ(k+1)
f2 = ZZ(k-1)
TZ = 2.0 / (f1 - f2)
ZTT= f1 - 2.0*ZZ(k) + f2
f3 = TZ * TZ
aw(k) = f3              ! C3
cw(k) = -ZTT * f3 * TZ  ! C9
dw(k) = 0.5 / (C1 + C2 + f3) ! 1/R7
end do   ! >>  11 flops


a(kst) = 0.0
c(kst) = -(aw(kst) + 0.5 * cw(kst)) * dw(kst) ! -R5/R7 = -(C3+0.5*C9) / R7      >> 3 flops

do k=kst+1, ked-1
f1 = aw(k)  ! C3
f2 = cw(k)  ! C9
aa3 = dw(k)
a(k) = -(f1 - 0.5 * f2) * aa3
c(k) = -(f1 + 0.5 * f2) * aa3
end do    !  >> 6 flops

a(ked) = -(aw(ked) - 0.5 * cw(ked)) * dw(ked) ! -R6/R7 = -(C3-0.5*C9) / R7     >> 3 flops
c(ked) = 0.0

! Source
do k = kst, ked
d(k) = (                         &
dd1 * x(k, i+1, j  )   &
+ dd2 * x(k, i-1, j  )   &
+ cc1 * x(k, i  , j+1)   &
+ cc2 * x(k, i  , j-1)   &
- rhs(k, i, j)           &
) * dw(k) * msk(k, i, j)
end do   !  >>  10 flops
! ここまで、dd1, dd2, cc1, cc2は再利用しない
! a, c, dを計算したので、aw, cw, dwは再利用可能

! BC   >>  12 flops
d(kst) = ( d(kst) + (aw(kst) - 0.5 * cw(kst)) * dw(kst) * x(kst-1, i, j) ) * msk(kst, i, j)
d(ked) = ( d(ked) + (aw(ked) + 0.5 * cw(ked)) * dw(ked) * x(ked+1, i, j) ) * msk(ked, i, j)


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
aw(k) =  -e * ap * a(kl)
cw(k) =  -e * cp * c(kr)
dw(k) =   e * ( d(k) - ap * d(kl) - cp * d(kr) )
end do   !  >> 16 flops

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = aw(k)
c(k) = cw(k)
d(k) = dw(k)
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
dw(k ) = dd1
dw(kr) = dd2
end do  ! >>  11 flops


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!$acc loop reduction(+:res1)
do k = kst, ked
pp =   x(k, i, j)
dp = ( dw(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

end do  !  >> 6 flops

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef _SVR
do k = kst, ked
res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

return
end subroutine pcr_maf


!********************************************************************************
subroutine pcr_eda_maf (sz, idx, g, pn, x, msk, rhs, XX, YY, ZZ, &
aw, cw, dw, omg, res, tmp, flop)
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
real, dimension(-1:sz(3)+2)              ::  aw, cw, dw
real, dimension(:), allocatable          ::  a, c, d
real                                     ::  ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3
real                                                   ::  C1, C2, C7, C8, GX, EY, TZ, ZTT
real, dimension(-1:sz(1)+2)                            ::  XX
real, dimension(-1:sz(2)+2)                            ::  YY
real, dimension(-1:sz(3)+2)                            ::  ZZ, tmp

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0

#ifdef _SVR
tmp = 0.0
#endif

s = 2**(pn-1)

flop = flop + dble(           &
(jed-jst+1)*(ied-ist+1)* (  &
( 24.0d0 +                & ! metrics
+ 3.0 * 2.0 + 12.0       & ! coef BC
)                         &
+ (ked-kst+1)* (11.0+10.0) & ! coef + Source
+ (ked-kst-1)* 6.0d0       & ! coef
+ (ked-kst+1)*(pn-1)*16.0  & ! PCR
+ s * 9.0           & ! 2x2
+ (ked-kst+1)*6.0          & ! Relaxation
)                          &
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


!$OMP PARALLEL &
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP private(ap, cp, e, s, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(aw, cw, dw) &
!$OMP firstprivate(a, c, d) &
!$OMP private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
!$OMP DO SCHEDULE(static) Collapse(2)
do j=jst, jed
do i=ist, ied

GX =  2.0 / (XX(i+1) - XX(i-1))
EY =  2.0 / (YY(j+1) - YY(j-1))
C1 =  GX * GX
C2 =  EY * EY
C7 = -(XX(i+1) - 2.0*XX(i) + XX(i-1)) * C1 * GX
C8 = -(YY(j+1) - 2.0*YY(j) + YY(j-1)) * C2 * EY
dd1= C1 + 0.5 * C7 ! R1
dd2= C1 - 0.5 * C7 ! R2
cc1= C2 + 0.5 * C8 ! R3
cc2= C2 - 0.5 * C8 ! R4      >>  24 flops

! Reflesh coef. due to override
do k=kst, ked
f1 = ZZ(k+1)
f2 = ZZ(k-1)
TZ = 2.0 / (f1 - f2)
ZTT= f1 - 2.0*ZZ(k) + f2
f3 = TZ * TZ
aw(k) = f3              ! C3
cw(k) = -ZTT * f3 * TZ  ! C9
dw(k) = 0.5 / (C1 + C2 + f3) ! 1/R7
end do   ! >>  11 flops


a(kst) = 0.0
c(kst) = -(aw(kst) + 0.5 * cw(kst)) * dw(kst) ! -R5/R7 = -(C3+0.5*C9) / R7      >> 3 flops

do k=kst+1, ked-1
f1 = aw(k)  ! C3
f2 = cw(k)  ! C9
aa3 = dw(k)
a(k) = -(f1 - 0.5 * f2) * aa3
c(k) = -(f1 + 0.5 * f2) * aa3
end do    !  >> 6 flops

a(ked) = -(aw(ked) - 0.5 * cw(ked)) * dw(ked) ! -R6/R7 = -(C3-0.5*C9) / R7     >> 3 flops
c(ked) = 0.0

! Source
do k = kst, ked
d(k) = (                         &
dd1 * x(k, i+1, j  )   &
+ dd2 * x(k, i-1, j  )   &
+ cc1 * x(k, i  , j+1)   &
+ cc2 * x(k, i  , j-1)   &
- rhs(k, i, j)           &
) * dw(k) * msk(k, i, j)
end do   !  >>  10 flops
! ここまで、dd1, dd2, cc1, cc2は再利用しない
! a, c, dを計算したので、aw, cw, dwは再利用可能

! BC   >>  12 flops
d(kst) = ( d(kst) + (aw(kst) - 0.5 * cw(kst)) * dw(kst) * x(kst-1, i, j) ) * msk(kst, i, j)
d(ked) = ( d(ked) + (aw(ked) + 0.5 * cw(ked)) * dw(ked) * x(ked+1, i, j) ) * msk(ked, i, j)


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
aw(k) =  -e * ap * a(k-s)
cw(k) =  -e * cp * c(k+s)
dw(k) =   e * ( d(k) - ap * d(k-s) - cp * d(k+s) )
end do   !  >> 16 flops

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = aw(k)
c(k) = cw(k)
d(k) = dw(k)
end do

end do ! p反復


! 最終段の反転
s = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
do k = kst, kst+s-1
cc1 = c(k)
aa2 = a(k+s)
f1  = d(k)
f2  = d(k+s)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
dw(k ) = dd1
dw(k+s)= dd2
end do  ! >>  9 flops



! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
do k = kst, ked
pp =   x(k, i, j)
dp = ( dw(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

end do  !  >> 6 flops

end do
end do
!$OMP END DO
!$OMP END PARALLEL

#ifdef _SVR
do k = kst, ked
res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

return
end subroutine pcr_eda_maf


!********************************************************************************
subroutine pcr_esa_maf (sz, idx, g, pn, s, x, msk, rhs, XX, YY, ZZ, &
a, c, d, aw, cw, dw, omg, res, tmp, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, p, sq, s
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(-1:sz(3)+2)              ::  aw, cw, dw
real, dimension(idx(4)-s:idx(5)+s)       ::  a, c, d
real                                     ::  ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3
real                                                   ::  C1, C2, C7, C8, GX, EY, TZ, ZTT
real, dimension(-1:sz(1)+2)                            ::  XX
real, dimension(-1:sz(2)+2)                            ::  YY
real, dimension(-1:sz(3)+2)                            ::  ZZ, tmp

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

res1 = 0.0

#ifdef _SVR
tmp = 0.0
#endif

flop = flop + dble(           &
(jed-jst+1)*(ied-ist+1)* (  &
( 24.0d0 +                & ! metrics
+ 3.0 * 2.0 + 12.0       & ! coef BC
)                         &
+ (ked-kst+1)* (11.0+10.0) & ! coef + Source
+ (ked-kst-1)* 6.0d0       & ! coef
+ (ked-kst+1)*(pn-1)*16.0  & ! PCR
+ 2**(pn-1) * 9.0           & ! 2x2
+ (ked-kst+1)*6.0          & ! Relaxation
)                          &
)

#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) reduction(+:res1) &
!$acc& private(a, c, d, aw, cw, dw) &
!$acc& private(ap, cp, e, sq, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$acc& private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
#else
!$OMP PARALLEL &
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP private(ap, cp, e, sq, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(C1, C2, C7, C8, GX, EY, TZ, ZTT) &
!$OMP private(aw, cw, dw) &
!$OMP firstprivate(a, c, d)
!$OMP DO SCHEDULE(static) Collapse(2)
#endif
do j=jst, jed
do i=ist, ied

GX =  2.0 / (XX(i+1) - XX(i-1))
EY =  2.0 / (YY(j+1) - YY(j-1))
C1 =  GX * GX
C2 =  EY * EY
C7 = -(XX(i+1) - 2.0*XX(i) + XX(i-1)) * C1 * GX
C8 = -(YY(j+1) - 2.0*YY(j) + YY(j-1)) * C2 * EY
dd1= C1 + 0.5 * C7 ! R1
dd2= C1 - 0.5 * C7 ! R2
cc1= C2 + 0.5 * C8 ! R3
cc2= C2 - 0.5 * C8 ! R4      >>  24 flops

! Reflesh coef. due to override
do k=kst, ked
f1 = ZZ(k+1)
f2 = ZZ(k-1)
TZ = 2.0 / (f1 - f2)
ZTT= f1 - 2.0*ZZ(k) + f2
f3 = TZ * TZ
aw(k) = f3              ! C3
cw(k) = -ZTT * f3 * TZ  ! C9
dw(k) = 0.5 / (C1 + C2 + f3) ! 1/R7
end do   ! >>  11 flops


a(kst) = 0.0
c(kst) = -(aw(kst) + 0.5 * cw(kst)) * dw(kst) ! -R5/R7 = -(C3+0.5*C9) / R7      >> 3 flops

do k=kst+1, ked-1
f1 = aw(k)  ! C3
f2 = cw(k)  ! C9
aa3 = dw(k)
a(k) = -(f1 - 0.5 * f2) * aa3
c(k) = -(f1 + 0.5 * f2) * aa3
end do    !  >> 6 flops

a(ked) = -(aw(ked) - 0.5 * cw(ked)) * dw(ked) ! -R6/R7 = -(C3-0.5*C9) / R7     >> 3 flops
c(ked) = 0.0

! Source
do k = kst, ked
d(k) = (                         &
  dd1 * x(k, i+1, j  )   &
+ dd2 * x(k, i-1, j  )   &
+ cc1 * x(k, i  , j+1)   &
+ cc2 * x(k, i  , j-1)   &
- rhs(k, i, j)           &
) * dw(k) * msk(k, i, j)
end do   !  >>  10 flops
! ここまで、dd1, dd2, cc1, cc2は再利用しない
! a, c, dを計算したので、aw, cw, dwは再利用可能

! BC   >>  12 flops
d(kst) = ( d(kst) + (aw(kst) - 0.5 * cw(kst)) * dw(kst) * x(kst-1, i, j) ) * msk(kst, i, j)
d(ked) = ( d(ked) + (aw(ked) + 0.5 * cw(ked)) * dw(ked) * x(ked+1, i, j) ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
!$acc loop seq
do p=1, pn-1
sq = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!pgi$ ivdep
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-sq) - cp * a(k+sq) )
aw(k) =   -e * ap * a(k-sq)
cw(k) =   -e * cp * c(k+sq)
dw(k) =   e * ( d(k) - ap * d(k-sq) - cp * d(k+sq) )
end do   !  >> 16 flops

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = aw(k)
c(k) = cw(k)
d(k) = dw(k)
end do

end do ! p反復


! 最終段の反転
sq = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!pgi$ ivdep
!$acc loop independent
do k = kst, kst+sq-1
cc1 = c(k)
aa2 = a(k+sq)
f1  = d(k)
f2  = d(k+sq)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
dw(k   ) = dd1
dw(k+sq) = dd2
end do  ! >>  9 flops


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!$acc loop reduction(+:res1)
do k = kst, ked
pp =   x(k, i, j)
dp = ( dw(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

end do  !  >> 6 flops

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif


#ifdef _SVR
do k = kst, ked
res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

return
end subroutine pcr_esa_maf


!********************************************************************************
subroutine pcr_rb_esa_maf (sz, idx, g, pn, ofst, color, s, x, msk, rhs, XX, YY, ZZ, &
a, c, d, aw, cw, dw, omg, res, tmp, flop)
implicit none
!args
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer                                                ::  g, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs
real                                                   ::  omg
double precision                                       ::  res, flop
! work
integer                                  ::  i, j, k, p, color, ip, ofst, sq, s
integer                                  ::  ist, ied, jst, jed, kst, ked
real, dimension(-1:sz(3)+2)              ::  aw, cw, dw
real, dimension(idx(4)-s:idx(5)+s)       ::  a, c, d
real                                     ::  ap, cp, e, pp, dp, res1
real                                     ::  jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3
real                                                   ::  C1, C2, C7, C8, GX, EY, TZ, ZTT
real, dimension(-1:sz(1)+2)                            ::  XX
real, dimension(-1:sz(2)+2)                            ::  YY
real, dimension(-1:sz(3)+2)                            ::  ZZ, tmp

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble(           &
(jed-jst+1)*(ied-ist+1)* (  &
( 24.0d0 +                & ! metrics
+ 3.0 * 2.0 + 12.0       & ! coef BC
)                         &
+ (ked-kst+1)* (11.0+10.0) & ! coef + Source
+ (ked-kst-1)* 6.0d0       & ! coef
+ (ked-kst+1)*(pn-1)*16.0  & ! PCR
+ 2**(pn-1) * 11.0           & ! 2x2
+ (ked-kst+1)*6.0          & ! Relaxation
)                          &
) * 0.5


ip = ofst + color
res1 = 0.0

#ifdef _SVR
tmp = 0.0
#endif

#ifdef _OPENACC
!$acc kernels
!$acc loop independent collapse(2) reduction(+:res1) &
!$acc& private(a, c, d, aw, cw, dw) &
!$acc& private(ap, cp, e, sq, p, k, pp, dp) &
!$acc& private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$acc& private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
#else
!$OMP PARALLEL &
#ifdef _SVR
!$OMP REDUCTION(+:tmp) &
#else
!$OMP REDUCTION(+:res1) &
#endif
!$OMP private(ap, cp, e, sq, p, k, pp, dp) &
!$OMP private(jj, dd1, dd2, aa2, aa3, cc1, cc2, f1, f2, f3) &
!$OMP private(aw, cw, dw) &
!$OMP firstprivate(a, c, d) &
!$OMP private(C1, C2, C7, C8, GX, EY, TZ, ZTT)
!$OMP DO SCHEDULE(static) collapse(2)
#endif
do j=jst, jed
do i=ist, ied
if(mod(i+j,2) /= color) cycle

!do i=ist+mod(j+ip,2), ied, 2

GX =  2.0 / (XX(i+1) - XX(i-1))
EY =  2.0 / (YY(j+1) - YY(j-1))
C1 =  GX * GX
C2 =  EY * EY
C7 = -(XX(i+1) - 2.0*XX(i) + XX(i-1)) * C1 * GX
C8 = -(YY(j+1) - 2.0*YY(j) + YY(j-1)) * C2 * EY
dd1= C1 + 0.5 * C7 ! R1
dd2= C1 - 0.5 * C7 ! R2
cc1= C2 + 0.5 * C8 ! R3
cc2= C2 - 0.5 * C8 ! R4      >>  24 flops

! Reflesh coef. due to override
do k=kst, ked
f1 = ZZ(k+1)
f2 = ZZ(k-1)
TZ = 2.0 / (f1 - f2)
ZTT= f1 - 2.0*ZZ(k) + f2
f3 = TZ * TZ
aw(k) = f3              ! C3
cw(k) = -ZTT * f3 * TZ  ! C9
dw(k) = 0.5 / (C1 + C2 + f3) ! 1/R7
end do   ! >>  11 flops


a(kst) = 0.0
c(kst) = -(aw(kst) + 0.5 * cw(kst)) * dw(kst) ! -R5/R7 = -(C3+0.5*C9) / R7      >> 3 flops

do k=kst+1, ked-1
f1 = aw(k)  ! C3
f2 = cw(k)  ! C9
aa3 = dw(k)
a(k) = -(f1 - 0.5 * f2) * aa3
c(k) = -(f1 + 0.5 * f2) * aa3
end do    !  >> 6 flops

a(ked) = -(aw(ked) - 0.5 * cw(ked)) * dw(ked) ! -R6/R7 = -(C3-0.5*C9) / R7     >> 3 flops
c(ked) = 0.0

! Source
do k = kst, ked
d(k) = (                         &
dd1 * x(k, i+1, j  )   &
+ dd2 * x(k, i-1, j  )   &
+ cc1 * x(k, i  , j+1)   &
+ cc2 * x(k, i  , j-1)   &
- rhs(k, i, j)           &
) * dw(k) * msk(k, i, j)
end do   !  >>  10 flops
! ここまで、dd1, dd2, cc1, cc2は再利用しない
! a, c, dを計算したので、aw, cw, dwは再利用可能

! BC   >>  12 flops
d(kst) = ( d(kst) + (aw(kst) - 0.5 * cw(kst)) * dw(kst) * x(kst-1, i, j) ) * msk(kst, i, j)
d(ked) = ( d(ked) + (aw(ked) + 0.5 * cw(ked)) * dw(ked) * x(ked+1, i, j) ) * msk(ked, i, j)


! PCR  最終段の一つ手前で停止
!$acc loop seq
do p=1, pn-1
sq = 2**(p-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
do k = kst, ked
ap = a(k)
cp = c(k)
e = 1.0 / ( 1.0 - ap * c(k-sq) - cp * a(k+sq) )
aw(k) =   -e * ap * a(k-sq)
cw(k) =   -e * cp * c(k+sq)
dw(k) =   e * ( d(k) - ap * d(k-sq) - cp * d(k+sq) )
end do   !  >> 16 flops

!dir$ vector aligned
!dir$ simd
do k = kst, ked
a(k) = aw(k)
c(k) = cw(k)
d(k) = dw(k)
end do

end do ! p反復


! 最終段の反転
sq = 2**(pn-1)

!dir$ vector aligned
!dir$ simd
!NEC$ IVDEP
!$acc loop independent
do k = kst, kst+sq-1
cc1 = c(k)
aa2 = a(k+sq)
f1  = d(k)
f2  = d(k+sq)
jj  = 1.0 / (1.0 - aa2 * cc1)
dd1 = (f1 - cc1 * f2) * jj
dd2 = (f2 - aa2 * f1) * jj
dw(k   ) = dd1
dw(k+sq) = dd2
end do  ! >>  11 flops


! a_{i-1} x_{i-2} + x_{i-1} + c_{i-1} x_i     = d_{i-1}
! a_{i}   x_{i-1} + x_{i}   + c_{i}   x_{i+1} = d_{i}
! a_{i+1} x_{i}   + x_{i+1} + c_{i+1} x_{i+2} = d_{i+1}


! Relaxation
!dir$ vector aligned
!dir$ simd
!$acc loop reduction(+:res1)
do k = kst, ked
pp =   x(k, i, j)
dp = ( dw(k) - pp ) * omg * msk(k, i, j)
x(k, i, j) = pp + dp

#ifdef _SVR
tmp(k) = tmp(k) + dp * dp ! 30
#else
res1 = res1 + dp * dp ! 30
#endif

end do  !  >> 6 flops

end do
end do
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef _SVR
do k = kst, ked
res1 = res1 + tmp(k)
end do
#endif

res = res + real(res1, kind=8)

return
end subroutine pcr_rb_esa_maf
