

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
subroutine tdma_p (nx, d, a, c, w)
implicit none
integer                          ::  i, nx
real                             ::  e
real, dimension(nx)              ::  d, w, a, c

w(1) = c(1);

do i=2, nx
  e = 1.0 / (1.0 - a(i) * w(i-1))
  w(i) = e * c(i)
  d(i) = (d(i) - a(i) * d(i-1)) * e
end do

do i=nx-1, 1, -1
  d(i) = d(i) - w(i) * d(i+1)
end do

return
end subroutine tdma_p


!> ********************************************************************
subroutine tdma_mp (nx, ms, d, a, c, w)
implicit none
integer                          ::  i, m, nx, ms
real                             ::  e
real, dimension(ms, 0:nx+1)      ::  d, w, a, c
!dir$ assume_aligned d:64, w:64, a:64, c:64

!$OMP PARALLEL DO SCHEDULE(static)
!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do m=1, ms
  w(m,1) = c(m,1)
end do

!$OMP PARALLEL DO SCHEDULE(static) private(e)
do i=2, nx
!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do m=1, ms
  e = 1.0 / (1.0 - a(m,i) * w(m,i-1))
  w(m,i) = e * c(m,i)
  d(m,i) = (d(m,i) - a(m,i) * d(m,i-1)) * e
end do
end do

!$OMP PARALLEL DO SCHEDULE(static)
do i=nx-1, 1, -1
!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do m=1, ms
  d(m,i) = d(m,i) - w(m,i) * d(m,i+1)
end do
end do

return
end subroutine tdma_mp


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
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in,out] flop flop count
!! @note あとでlsor_lu_rhs1でマスクをかける
!<
subroutine lsor_lu_rhs_jd (d, sz, idx, g, j, x, rhs, msk, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, kst
integer                                                ::  ied, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, rhs, msk
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,rhs:64,msk:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(ked-kst+1))*5.0

do k = kst, ked

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i = ist, ied
  d(i  ,k,j) = ( x(i-1, k  , j  ) &
             +   x(i+1, k  , j  ) &
             +   x(i  , k-1, j  ) &
             +   x(i  , k+1, j  ) &
             +   x(i  , k  , j-1) &
             +   x(i  , k  , j+1) ) * r &
             + rhs(i  ,k,j)
end do
end do

return
end subroutine lsor_lu_rhs_jd


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine lsor_relax_d (d, sz, idx, g, j, x, msk, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg
real                                                   ::  dp1, pp1, pn1
real                                                   ::  dp2, pp2, pn2
real                                                   ::  dp3, pp3, pn3
real                                                   ::  dp4, pp4, pn4
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk
!dir$ assume_aligned d:64,x:64,msk:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ked-kst+1)*(ied-ist+1))*5.0

do k = kst, ked

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i = ist, ied
  dp1 = d(i  , k  , j  ) * omg * msk(i  , k  , j  )
  x(i  , k  , j  ) = x(i  , k  , j  ) + dp1

  res = res + dp1 * dp1
end do
end do

return
end subroutine lsor_relax_d




!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine lsor_lu_rhs_k4 (d, sz, idx, g, j, k, x, d2, msk, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, ff
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk, d2, rhs
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,msk:64,d2:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
ff = dble(ied-ist+1)

flop = flop + ff * 4.0

!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i  ,k, j) = (d2(i  , k, j  ) &
              + ( x(i-1, k, j  ) &
              +   x(i+1, k, j  ) ) * r ) &
              * msk(i  , k, j  )
end do

if (k == kst) then
  !dir$ vector aligned
  !dir$ simd
  do i=ist, ied
    d(i  , kst, j  ) = ( d(i  , kst  , j  ) &
                     + rhs(i  , kst-1, j  ) * r ) &
                     * msk(i  , kst  , j  )
  end do
  flop = flop + ff*3.0
endif

if (k == ked) then
  !dir$ vector aligned
  !dir$ simd
  do i=ist, ied
    d(i  , ked, j  ) = ( d(i  , ked  , j  ) &
                     + rhs(i  , ked+1, j  ) * r ) &
                     * msk(i  , ked  , j  )
  end do
  flop = flop + ff*3.0
endif

return
end subroutine lsor_lu_rhs_k4

!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_tdma_f4 (d, sz, idx, g, j, d2, x, a, e, msk, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, ff
real                                                   ::  aa, ee, r
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e
!dir$ assume_aligned d:64, a:64, e:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
ff = dble( (ied-ist+1)*(ked-kst+1) )

flop = flop + ff * 4.0

do k = kst, ked, 4
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i, k  , j) = (d2(i  , k  , j) &
               + ( x(i-1, k  , j) &
               +   x(i+1, k  , j) ) * r ) &
               * msk(i  , k  , j)

  d(i, k+1, j) = (d2(i  , k+1, j) &
               + ( x(i-1, k+1, j) &
               +   x(i+1, k+1, j) ) * r ) &
               * msk(i  , k+1, j)

  d(i, k+2, j) = (d2(i  , k+2, j) &
               + ( x(i-1, k+2, j) &
               +   x(i+1, k+2, j) ) * r ) &
               * msk(i  , k+2, j)

  d(i, k+3, j) = (d2(i  , k+3, j) &
               + ( x(i-1, k+3, j) &
               +   x(i+1, k+3, j) ) * r ) &
               * msk(i  , k+3, j)
end do
end do


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, kst, j) = ( d(i, kst  , j) &
               + rhs(i, kst-1, j) * r ) &
               * msk(i, kst  , j)
end do

flop = flop + ff*3.0


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, ked, j) = ( d(i, ked  , j) &
               + rhs(i, ked+1, j) * r ) &
               * msk(i, ked  , j)
end do

flop = flop + ff*3.0



do k = 3, ked, 4

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
  d(i, k+1, j) = (d(i, k+1, j) - a(k+1) * d(i, k  , j)) * e(k+1)
  d(i, k+2, j) = (d(i, k+2, j) - a(k+2) * d(i, k+1, j)) * e(k+2)
  d(i, k+3, j) = (d(i, k+3, j) - a(k+3) * d(i, k+2, j)) * e(k+3)
end do
end do

flop = flop + ff*3.0

return
end subroutine lsor_tdma_f4


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ms_rhs8 (d, sz, idx, g, ia, ja, x, rhs, msk, flop)
implicit none
integer                                                ::  k, g
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer, dimension(8)                                  ::  ia, ja
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x, rhs, msk
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,rhs:64,msk:64

kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble(ked-kst+1)*48.0

!dir$ vector aligned
!dir$ ivdep
do k = kst, ked
  d(k,ia(1),ja(1)) = (( x(k, ia(1)-1, ja(1)  ) &
                   +    x(k, ia(1)+1, ja(1)  ) &
                   +    x(k, ia(1)  , ja(1)-1) &
                   +    x(k, ia(1)  , ja(1)+1) ) * r &
                   +  rhs(k, ia(1)  , ja(1)) ) &
                   *  msk(k, ia(1)  , ja(1))

  d(k,ia(2),ja(2)) = (( x(k, ia(2)-1, ja(2)  ) &
                   +    x(k, ia(2)+1, ja(2)  ) &
                   +    x(k, ia(2)  , ja(2)-1) &
                   +    x(k, ia(2)  , ja(2)+1) ) * r &
                   +  rhs(k, ia(2)  , ja(2)) ) &
                   *  msk(k, ia(2)  , ja(2))

  d(k,ia(3),ja(3)) = (( x(k, ia(3)-1, ja(3)  ) &
                   +    x(k, ia(3)+1, ja(3)  ) &
                   +    x(k, ia(3)  , ja(3)-1) &
                   +    x(k, ia(3)  , ja(3)+1) ) * r &
                   +  rhs(k, ia(3)  , ja(3)) ) &
                   *  msk(k, ia(3)  , ja(3))

  d(k,ia(4),ja(4)) = (( x(k, ia(4)-1, ja(4)  ) &
                   +    x(k, ia(4)+1, ja(4)  ) &
                   +    x(k, ia(4)  , ja(4)-1) &
                   +    x(k, ia(4)  , ja(4)+1) ) * r &
                   +  rhs(k, ia(4)  , ja(4)) ) &
                   *  msk(k, ia(4)  , ja(4))

  d(k,ia(5),ja(5)) = (( x(k, ia(5)-1, ja(5)  ) &
                   +    x(k, ia(5)+1, ja(5)  ) &
                   +    x(k, ia(5)  , ja(5)-1) &
                   +    x(k, ia(5)  , ja(5)+1) ) * r &
                   +  rhs(k, ia(5)  , ja(5)) ) &
                   *  msk(k, ia(5)  , ja(5))

  d(k,ia(6),ja(6)) = (( x(k, ia(6)-1, ja(6)  ) &
                   +    x(k, ia(6)+1, ja(6)  ) &
                   +    x(k, ia(6)  , ja(6)-1) &
                   +    x(k, ia(6)  , ja(6)+1) ) * r &
                   +  rhs(k, ia(6)  , ja(6)) ) &
                   *  msk(k, ia(6)  , ja(6))

  d(k,ia(7),ja(7)) = (( x(k, ia(7)-1, ja(7)  ) &
                   +    x(k, ia(7)+1, ja(7)  ) &
                   +    x(k, ia(7)  , ja(7)-1) &
                   +    x(k, ia(7)  , ja(7)+1) ) * r &
                   +  rhs(k, ia(7)  , ja(7)) ) &
                   *  msk(k, ia(7)  , ja(7))

  d(k,ia(8),ja(8)) = (( x(k, ia(8)-1, ja(8)  ) &
                   +    x(k, ia(8)+1, ja(8)  ) &
                   +    x(k, ia(8)  , ja(8)-1) &
                   +    x(k, ia(8)  , ja(8)+1) ) * r &
                   +  rhs(k, ia(8)  , ja(8)) ) &
                   *  msk(k, ia(8)  , ja(8))
end do

return
end subroutine ms_rhs8


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine ms_relax8 (d, sz, idx, g, ia, ja, x, msk, omg, res, flop)
implicit none
integer                                                ::  k, g
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer, dimension(8)                                  ::  ia, ja
double precision                                       ::  flop, res
real                                                   ::  omg
real                                                   ::  dp1, pp1, pn1
real                                                   ::  dp2, pp2, pn2
real                                                   ::  dp3, pp3, pn3
real                                                   ::  dp4, pp4, pn4
real                                                   ::  dp5, pp5, pn5
real                                                   ::  dp6, pp6, pn6
real                                                   ::  dp7, pp7, pn7
real                                                   ::  dp8, pp8, pn8
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x, msk
!dir$ assume_aligned d:64,x:64,msk:64

kst = idx(4)
ked = idx(5)

flop = flop + dble(ked-kst+1)*48.0

!dir$ vector aligned
!dir$ ivdep
do k = kst, ked
  pp1 =   x(k, ia(1), ja(1))
  dp1 = ( d(k, ia(1), ja(1)) - pp1 ) * omg * msk(k, ia(1), ja(1))
  pn1 = pp1 + dp1
  x(k, ia(1), ja(1)) = pn1

  pp2 =   x(k, ia(2), ja(2))
  dp2 = ( d(k, ia(2), ja(2)) - pp2 ) * omg * msk(k, ia(2), ja(2))
  pn2 = pp2 + dp2
  x(k, ia(2), ja(2)) = pn2

  pp3 =   x(k, ia(3), ja(3))
  dp3 = ( d(k, ia(3), ja(3)) - pp3 ) * omg * msk(k, ia(3), ja(3))
  pn3 = pp3 + dp3
  x(k, ia(3), ja(3)) = pn3

  pp4 =   x(k, ia(4), ja(4))
  dp4 = ( d(k, ia(4), ja(4)) - pp4 ) * omg * msk(k, ia(4), ja(4))
  pn4 = pp4 + dp4
  x(k, ia(4), ja(4)) = pn4

  pp5 =   x(k, ia(5), ja(5))
  dp5 = ( d(k, ia(5), ja(5)) - pp5 ) * omg * msk(k, ia(5), ja(5))
  pn5 = pp5 + dp5
  x(k, ia(5), ja(5)) = pn5

  pp6 =   x(k, ia(6), ja(6))
  dp6 = ( d(k, ia(6), ja(6)) - pp6 ) * omg * msk(k, ia(6), ja(6))
  pn6 = pp6 + dp6
  x(k, ia(6), ja(6)) = pn6

  pp7 =   x(k, ia(7), ja(7))
  dp7 = ( d(k, ia(7), ja(7)) - pp7 ) * omg * msk(k, ia(7), ja(7))
  pn7 = pp7 + dp7
  x(k, ia(7), ja(7)) = pn7

  pp8 =   x(k, ia(8), ja(8))
  dp8 = ( d(k, ia(8), ja(8)) - pp8 ) * omg * msk(k, ia(8), ja(8))
  pn8 = pp8 + dp8
  x(k, ia(8), ja(8)) = pn8

  res = res + dp1 * dp1 &
            + dp2 * dp2 &
            + dp3 * dp3 &
            + dp4 * dp4 &
            + dp5 * dp5 &
            + dp6 * dp6 &
            + dp7 * dp7 &
            + dp8 * dp8
end do

return
end subroutine ms_relax8



!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ms_bc (d, sz, idx, g, ia, ja, rhs, msk, flop)
implicit none
integer                                                ::  k, g, l
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer, dimension(8)                                  ::  ia, ja
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, rhs, msk
real                                                   ::  r
!dir$ assume_aligned d:64,rhs:64,msk:64

kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + 48.0;

!dir$ vector aligned
!dir$ ivdep
!dir$ unroll(8)
do l = 1, 8
  d(kst, ia(l), ja(l)) = ( d(kst  , ia(l), ja(l)) &
                       + rhs(kst-1, ia(l), ja(l)) * r ) &
                       * msk(kst  , ia(l), ja(l))
  d(ked, ia(l), ja(l)) = ( d(ked  , ia(l), ja(l)) &
                       + rhs(ked+1, ia(l), ja(l)) * r ) &
                       * msk(ked  , ia(l), ja(l))
end do

return
end subroutine ms_bc


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine ms_tdma (d, sz, idx, g, ia, ja, a, e, w, flop)
implicit none
integer                                                ::  k, g
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
integer, dimension(8)                                  ::  ia, ja
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
!dir$ assume_aligned d:64,a:64,e:64,w:64

ked = idx(5)

flop = flop + dble(ked-2)*40.0

!dir$ vector aligned
!dir$ ivdep
!dir$ unroll(8)
do k=3, ked
  d(k,ia(1),ja(1)) = (d(k,ia(1),ja(1)) - a(k) * d(k-1,ia(1),ja(1))) * e(k)
  d(k,ia(2),ja(2)) = (d(k,ia(2),ja(2)) - a(k) * d(k-1,ia(2),ja(2))) * e(k)
  d(k,ia(3),ja(3)) = (d(k,ia(3),ja(3)) - a(k) * d(k-1,ia(3),ja(3))) * e(k)
  d(k,ia(4),ja(4)) = (d(k,ia(4),ja(4)) - a(k) * d(k-1,ia(4),ja(4))) * e(k)
  d(k,ia(5),ja(5)) = (d(k,ia(5),ja(5)) - a(k) * d(k-1,ia(5),ja(5))) * e(k)
  d(k,ia(6),ja(6)) = (d(k,ia(6),ja(6)) - a(k) * d(k-1,ia(6),ja(6))) * e(k)
  d(k,ia(7),ja(7)) = (d(k,ia(7),ja(7)) - a(k) * d(k-1,ia(7),ja(7))) * e(k)
  d(k,ia(8),ja(8)) = (d(k,ia(8),ja(8)) - a(k) * d(k-1,ia(8),ja(8))) * e(k)
end do


!dir$ vector aligned
!dir$ ivdep
!dir$ unroll(8)
do k=ked-1, 2, -1
  d(k,ia(1),ja(1)) = d(k,ia(1),ja(1)) - w(k) * d(k+1,ia(1),ja(1))
  d(k,ia(2),ja(2)) = d(k,ia(2),ja(2)) - w(k) * d(k+1,ia(2),ja(2))
  d(k,ia(3),ja(3)) = d(k,ia(3),ja(3)) - w(k) * d(k+1,ia(3),ja(3))
  d(k,ia(4),ja(4)) = d(k,ia(4),ja(4)) - w(k) * d(k+1,ia(4),ja(4))
  d(k,ia(5),ja(5)) = d(k,ia(5),ja(5)) - w(k) * d(k+1,ia(5),ja(5))
  d(k,ia(6),ja(6)) = d(k,ia(6),ja(6)) - w(k) * d(k+1,ia(6),ja(6))
  d(k,ia(7),ja(7)) = d(k,ia(7),ja(7)) - w(k) * d(k+1,ia(7),ja(7))
  d(k,ia(8),ja(8)) = d(k,ia(8),ja(8)) - w(k) * d(k+1,ia(8),ja(8))
end do

return
end subroutine ms_tdma



!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine tdma_lsor_a (d, sz, idx, g, x, w, a, b, c, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g, nn
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, dummy, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, w, rhs
real                                                   ::  r
real, dimension(1-g:sz(1)+g)                           ::  a, b, c

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)
r = 1.0/6.0
nn = ied - ist + 1


flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(4.0+16.0+5.0)

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn)
do k = kst, ked
do j = jst, jed

  !$DIR IVDEP
  do i = ist, ied
    d(i,j,k) = ( x(i  ,j+1,k  ) &
             +   x(i  ,j-1,k  ) &
             +   x(i  ,j  ,k+1) &
             +   x(i  ,j  ,k-1) ) * r + rhs(i,j,k)
  end do

  d(ist,j,k) = d(ist,j,k) + rhs(ist-1,j,k)*r
  d(ied,j,k) = d(ied,j,k) + rhs(ied+1,j,k)*r

  call tdma_0(nn, d(ist,j,k), a(ist), b(ist), c(ist), w(ist,j,k))

  !$DIR IVDEP
  do i = ist, ied
    pp = x(i,j,k)
    dp = ( d(i,j,k) - pp ) * omg
    pn = pp + dp
    x(i,j,k) = pn
    res = res + dp*dp
  end do

end do
end do
!$OMP END PARALLEL DO


return
end subroutine tdma_lsor_a


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine tdma_lsor_d (d, sz, idx, g, x, w, a, b, c, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, w, rhs
real                                                   ::  r, e, y
real, dimension(1-g:sz(3)+g)                           ::  a, b, c

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(4.0+16.0+5.0) &
            + dble((ied-ist+1)*(jed-jst+1)) * 2.0;

y = 1.0 / b(kst)

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn, e)
do j = jst, jed
do i = ist, ied

  !$DIR IVDEP
  do k = kst, ked
     d(i,j,k) = ( x(i-1,j  ,k  ) &
              +   x(i+1,j  ,k  ) &
              +   x(i  ,j-1,k  ) &
              +   x(i  ,j+1,k  ) ) * r + rhs(i,j,k)
  end do
  d(i,j,kst) = d(i,j,kst) + rhs(i,j,kst-1)*r
  d(i,j,ked) = d(i,j,ked) + rhs(i,j,ked+1)*r

  ! TDMA
  d(i,j,kst) = d(i,j,kst) * y
  w(i,j,kst) = c(kst) * y

  !$DIR VECTOR
  do k=kst+1, ked
    e = 1.0 / (b(k) - a(k) * w(i,j,k-1))
    w(i,j,k) = e * c(k)
    d(i,j,k) = (d(i,j,k) - a(k) * d(i,j,k-1)) * e
  end do

  do k=ked-1, kst, -1
    d(i,j,k) = d(i,j,k) - w(i,j,k) * d(i,j,k+1)
  end do


  !$DIR IVDEP
  do k = kst, ked
    pp = x(i,j,k)
    dp = ( d(i,j,k) - pp ) * omg
    pn = pp + dp
    x(i,j,k) = pn
    res = res + dp*dp
  end do

end do
end do
!$OMP END PARALLEL DO


return
end subroutine tdma_lsor_d


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ljcb_f04 (d, sz, idx, g, x, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, rhs
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*5.0

!$OMP PARALLEL DO SCHEDULE(static)
do k = kst, ked
do j = jst, jed, 4
do i = ist, ied
  d(i,j  ,k) = ( x(i-1,j  ,k) &
             +   x(i+1,j  ,k) &
             +   x(i  ,j-1,k) &
             +   x(i  ,j+1,k) ) * r + rhs(i,j,k)

  d(i,j+1,k) = ( x(i-1,j+1,k) &
             +   x(i+1,j+1,k) &
             +   x(i  ,j  ,k) &
             +   x(i  ,j+2,k) ) * r + rhs(i,j+1,k)

  d(i,j+2,k) = ( x(i-1,j+2,k) &
             +   x(i+1,j+2,k) &
             +   x(i  ,j+1,k) &
             +   x(i  ,j+3,k) ) * r + rhs(i,j+2,k)

  d(i,j+3,k) = ( x(i-1,j+3,k) &
             +   x(i+1,j+3,k) &
             +   x(i  ,j+2,k) &
             +   x(i  ,j+4,k) ) * r + rhs(i,j+3,k)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f04

!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ljcb_f0t (d, sz, idx, g, x, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x, rhs
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*5.0

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
do k = kst, ked
  d(k,i,j) = ( x(k,i-1,j  ) &
           +   x(k,i+1,j  ) &
           +   x(k,i  ,j-1) &
           +   x(k,i  ,j+1) ) * r + rhs(k,i,j)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f0t


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ljcb_f1 (d, sz, idx, g, w, b, c, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, w, rhs
real                                                   ::  r, z, y
real, dimension(1-g:sz(3)+g)                           ::  b, c
!dir$ assume_aligned d:64,w:64,b:64,c:64,rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)) * 5.0;

y = 1.0 / b(kst)
z = c(kst) * y

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  d(i,j,kst) = ( d(i,j,kst) + rhs(i,j,kst-1)*r ) * y
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  d(i,j,ked) = d(i,j,ked) + rhs(i,j,ked+1)*r
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  w(i,j,kst) = z
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f1


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ljcb_f1t (d, sz, idx, g, w, b, c, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, w, rhs
real                                                   ::  r, z, y
real, dimension(1-g:sz(3)+g)                           ::  b, c
!dir$ assume_aligned d:64,w:64,b:64,c:64,rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)) * 5.0;

y = 1.0 / b(kst)
z = c(kst) * y

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  d(kst,i,j) = ( d(kst,i,j) + rhs(kst-1,i,j)*r ) * y
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  d(ked,i,j) = d(ked,i,j) + rhs(ked+1,i,j)*r
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  w(kst,i,j) = z
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f1t


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in,out] flop flop count
!<
subroutine ljcb_f2 (d, sz, idx, g, w, a, b, c, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, w
real                                                   ::  e, aa, bb, cc
real, dimension(1-g:sz(3)+g)                           ::  a, b, c
!dir$ assume_aligned d:64,w:64,a:64,b:64,c:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*14.0


do k=kst+1, ked

aa = a(k)
bb = b(k)
cc = c(k)

!$OMP PARALLEL DO SCHEDULE(static) PRIVATE(e)
do j = jst, jed
do i = ist, ied
  e = 1.0 / (bb - aa * w(i,j,k-1))
  w(i,j,k) = e * cc
  d(i,j,k) = (d(i,j,k) - aa * d(i,j,k-1)) * e
end do
end do
!$OMP END PARALLEL DO

end do

return
end subroutine ljcb_f2

!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in,out] flop flop count
!<
subroutine ljcb_f2t (d, sz, idx, g, w, a, b, c, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, w
real                                                   ::  e, aa, bb, cc
real, dimension(1-g:sz(3)+g)                           ::  a, b, c
!dir$ assume_aligned d:64,w:64,a:64,b:64,c:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*14.0


!$OMP PARALLEL DO SCHEDULE(static) PRIVATE(e)
do j = jst, jed
do i = ist, ied
do k = kst+1, ked
  aa = a(k)
  bb = b(k)
  cc = c(k)
  e = 1.0 / (bb - aa * w(k-1,i,j))
  w(k,i,j) = e * cc
  d(k,i,j) = (d(k,i,j) - aa * d(k-1,i,j)) * e
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f2t


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in,out] flop flop count
!<
subroutine ljcb_f3 (d, sz, idx, g, w, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, w
!dir$ assume_aligned d:64,w:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*2.0


do k=ked-1, kst, -1

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
  d(i,j,k) = d(i,j,k) - w(i,j,k) * d(i,j,k+1)
end do
end do
!$OMP END PARALLEL DO

end do

return
end subroutine ljcb_f3


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in,out] flop flop count
!<
subroutine ljcb_f3t (d, sz, idx, g, w, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, w
!dir$ assume_aligned d:64,w:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*2.0


!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
do i = ist, ied
do k = ked-1, kst, -1
  d(k,i,j) = d(k,i,j) - w(k,i,j) * d(k+1,i,j)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f3t


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine ljcb_f4t (d, sz, idx, g, x, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x
!dir$ assume_aligned d:64,x:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*5.0

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn)
do j = jst, jed
do i = ist, ied
do k = kst, ked
  pp = x(k,i,j)
  dp = ( d(k,i,j) - pp ) * omg
  pn = pp + dp
  x(k,i,j) = pn
  res = res + dp*dp
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f4t


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine tdma_ljcb_a (d, sz, idx, g, x, w, a, b, c, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, w, rhs
real                                                   ::  r, e, y
real, dimension(1-g:sz(3)+g)                           ::  a, b, c
!dir$ assume_aligned d:64,x:64,w:64,a:64,b:64,c:64,rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(4.0+16.0+5.0) &
            + dble((ied-ist+1)*(jed-jst+1)) * 6.0;

y = 1.0 / b(kst)

!$OMP PARALLEL DO SCHEDULE(static)
do k = kst, ked
do j = jst, jed
!dir$ simd vectorlength(16)
do i = ist, ied
  d(i,j,k) = ( x(i-1,j  ,k  ) &
           +   x(i+1,j  ,k  ) &
           +   x(i  ,j-1,k  ) &
           +   x(i  ,j+1,k  ) ) * r + rhs(i,j,k)
end do
end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
!dir$ simd vectorlength(16)
do i = ist, ied
  d(i,j,kst) = d(i,j,kst) + rhs(i,j,kst-1)*r
  d(i,j,ked) = d(i,j,ked) + rhs(i,j,ked+1)*r

  ! TDMA
  d(i,j,kst) = d(i,j,kst) * y
  w(i,j,kst) = c(kst) * y
end do
end do
!$OMP END PARALLEL DO


do k=kst+1, ked

!$OMP PARALLEL DO SCHEDULE(static) PRIVATE(e)
do j = jst, jed
!dir$ simd vectorlength(16)
do i = ist, ied
  e = 1.0 / (b(k) - a(k) * w(i,j,k-1))
  w(i,j,k) = e * c(k)
  d(i,j,k) = (d(i,j,k) - a(k) * d(i,j,k-1)) * e
end do
end do
!$OMP END PARALLEL DO

end do



do k=ked-1, kst, -1

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
!dir$ simd vectorlength(16)
do i = ist, ied
  d(i,j,k) = d(i,j,k) - w(i,j,k) * d(i,j,k+1)
end do
end do
!$OMP END PARALLEL DO

end do


!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn)
do k = kst, ked
do j = jst, jed
!dir$ simd vectorlength(16)
do i = ist, ied
  pp = x(i,j,k)
  dp = ( d(i,j,k) - pp ) * omg
  pn = pp + dp
  x(i,j,k) = pn
  res = res + dp*dp
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine tdma_ljcb_a


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in]     m    mask
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine tdma_ljcb_b(d, sz, idx, g, x, w, a, b, c, rhs, m, omg, res, flop)
implicit none
integer                                                ::  i, k, g, js
integer                                                ::  ist, jst, kst, st
integer                                                ::  ied, jed, ked, ed
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension((sz(1)+2*g)*(sz(2)+2*g), 1-g:sz(3)+g)  ::  d, x, w, m, rhs
real                                                   ::  r, e, y ,z
real, dimension(1-g:sz(3)+g)                           ::  a, b, c
!dir$ assume_aligned d:64,x:64,w:64,m:64,a:64,b:64,c:64,rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

js  = sz(1) + 2 * g
st  = ist+g + js * (jst+g-1)
ed  = ied+g + js * (jed+g-1)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(5.0+16.0+5.0) &
            + dble((ied-ist+1)*(jed-jst+1)) * 6.0;

y = 1.0 / b(kst)
z = c(kst) * y

!$OMP PARALLEL DO
do k = kst, ked
! !dir$ simd vectorlength(16)
do i = st, ed
  d(i,k) =(( x(i-1 ,k  ) &
         +   x(i+1 ,k  ) &
         +   x(i-js,k  ) &
         +   x(i+js,k  ) ) * r + rhs(i,k) )* m(i,k)
end do
end do
!$OMP END PARALLEL DO



!$OMP PARALLEL DO
do i = st, ed
  d(i,kst) = ( (d(i,kst) + rhs(i,kst-1)*r) * m(i,kst) ) * y
  d(i,ked) = (  d(i,ked) + rhs(i,ked+1)*r) * m(i,ked)

  ! TDMA
  !d(i,kst) = d(i,kst) * y
  w(i,kst) = z
end do
!$OMP END PARALLEL DO



do k = kst+1, ked
!$OMP PARALLEL DO PRIVATE(e)
!dir$ simd vectorlength(16)
do i = st, ed
  e = 1.0 / (b(k) - a(k) * w(i,k-1))
  w(i,k) = e * c(k)
  d(i,k) = (d(i,k) - a(k) * d(i,k-1)) * e
end do
!$OMP END PARALLEL DO
end do



do k=ked-1, kst, -1
!$OMP PARALLEL DO
!dir$ simd vectorlength(16)
do i = st, ed
  d(i,k) = d(i,k) - w(i,k) * d(i,k+1)
end do
!$OMP END PARALLEL DO
end do



do k = kst, ked
!$OMP PARALLEL DO REDUCTION(+:res) PRIVATE(pp, dp, pn)
!dir$ simd vectorlength(16)
do i = st, ed
  pp = x(i,k)
  dp = ( d(i,k) - pp ) * omg * m(i,k)
  pn = pp + dp
  x(i,k) = pn
  res = res + dp*dp
end do
!$OMP END PARALLEL DO
end do


return
end subroutine tdma_ljcb_b


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine lsor_lu_rhs_k (d, sz, idx, g, j, k, x, d2, msk, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk, d2
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,msk:64,d2:64

ist = idx(0)
ied = idx(1)

r = 1.0/6.0

flop = flop + dble(ied-ist+1)*4.0

!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i  ,k, j) = (d2(i  , k, j  ) &
              + ( x(i-1, k, j  ) &
              +   x(i+1, k, j  ) ) * r ) &
              * msk(i  , k, j  )
end do

return
end subroutine lsor_lu_rhs_k




!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine lsor_lu_bc_kst (d, sz, idx, g, j, rhs, msk, flop)
implicit none
integer                                                ::  i, j, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, rhs, msk
real                                                   ::  r
!dir$ assume_aligned d:64,rhs:64,msk:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)
r = 1.0/6.0
flop = flop + dble(ied-ist+1)*3.0;

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i  , kst, j  ) = ( d(i  , kst  , j  ) &
                   + rhs(i  , kst-1, j  ) * r ) &
                   * msk(i  , kst  , j  )
end do

return
end subroutine lsor_lu_bc_kst


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine lsor_lu_bc_ked (d, sz, idx, g, j, rhs, msk, flop)
implicit none
integer                                                ::  i, j, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, rhs, msk
real                                                   ::  r
!dir$ assume_aligned d:64,rhs:64,msk:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)
r = 1.0/6.0
flop = flop + dble(ied-ist+1)*3.0;

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i  , ked, j  ) = ( d(i  , ked  , j  ) &
                   + rhs(i  , ked+1, j  ) * r ) &
                   * msk(i  , ked  , j  )
end do

return
end subroutine lsor_lu_bc_ked


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_tdma_f (d, sz, idx, g, j, k, a, e, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real                                                   ::  aa, ee
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d
real, dimension(1-g:sz(3)+g)                           ::  a, e
!dir$ assume_aligned d:64, a:64, e:64

ist = idx(0)
ied = idx(1)

flop = flop + dble(ied-ist+1)*3.0

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k, j  ) = (d(i, k, j  ) - a(k) * d(i, k-1, j  )) * e(k)
end do

return
end subroutine lsor_tdma_f


!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     x    解ベクトル
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine tdma_lsor_b (d, sz, idx, g, x, w, a, c, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, w, rhs
real                                                   ::  r, e
real, dimension(1-g:sz(1)+g)                           ::  a, c
!dir$ assume_aligned d:64, x:64, w:64, rhs:64, a:64, c:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(5.0+16.0+5.0) &
            + dble((ked-kst+1)*(jed-jst+1)) * 10.0;

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn, e)
do k = kst, ked
do j = jst, jed

  !dir$ vector aligned
  !dir$ simd
  do i = ist, ied
    d(i,j,k) = ( x(i  ,j-1,k  ) &
             +   x(i  ,j+1,k  ) &
             +   x(i  ,j  ,k-1) &
             +   x(i  ,j  ,k+1) ) * r + rhs(i,j,k)
  end do

  d(ist,j,k) = d(ist,j,k) + rhs(ist-1,j,k)*r
  d(ied,j,k) = d(ied,j,k) + rhs(ied+1,j,k)*r

  ! TDMA
  w(ist,j,k) = c(ist)

  !dir$ vector aligned
  do i = ist, ied
    e = 1.0 / (1.0 - a(i) * w(i-1,j,k))
    w(i,j,k) = e * c(i)
    d(i,j,k) = (d(i,j,k) - a(i) * d(i-1,j,k)) * e
  end do

  !dir$ vector aligned
  do i=ied-1, ist, -1
    d(i,j,k) = d(i,j,k) - w(i,j,k) * d(i+1,j,k)
  end do


  !dir$ vector aligned
  !dir$ simd
  do i = ist, ied
    pp = x(i,j,k)
    dp = ( d(i,j,k) - pp ) * omg
    pn = pp + dp
    x(i,j,k) = pn
    res = res + dp*dp
  end do

end do
end do
!$OMP END PARALLEL DO


return
end subroutine tdma_lsor_b


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_tdma_fwd (d, sz, idx, g, j, d2, x, a, e, msk, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, ff
real                                                   ::  aa, ee, r
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e
!dir$ assume_aligned d:64, a:64, e:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
ff = dble( (ied-ist+1)*(ked-kst+1) )

flop = flop + ff * 4.0

do k = kst, ked
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i, k  , j) = (d2(i  , k  , j) &
               + ( x(i-1, k  , j) &
               +   x(i+1, k  , j) ) * r ) &
               * msk(i  , k  , j)
end do
end do


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, kst, j) = ( d(i, kst  , j) &
               + rhs(i, kst-1, j) * r ) &
               * msk(i, kst  , j)
end do

flop = flop + ff*3.0


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, ked, j) = ( d(i, ked  , j) &
               + rhs(i, ked+1, j) * r ) &
               * msk(i, ked  , j)
end do

flop = flop + ff*3.0



do k = 3, ked

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
end do
end do

flop = flop + ff*3.0

return
end subroutine lsor_tdma_fwd
