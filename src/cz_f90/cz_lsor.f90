!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################

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
subroutine lsor_lu_rhs_j (d, sz, idx, g, j, x, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, kst
integer                                                ::  ied, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, rhs
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,rhs:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(ked-kst+1))*12.0

do k = kst, ked

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i = ist, ied, 4
  d(i  ,k,j) = ( x(i  , k  , j-1) &
             +   x(i  , k  , j+1) ) * r + rhs(i  ,k,j)

  d(i+1,k,j) = ( x(i+1, k  , j-1) &
             +   x(i+1, k  , j+1) ) * r + rhs(i+1,k,j)

  d(i+2,k,j) = ( x(i+2, k  , j-1) &
             +   x(i+2, k  , j+1) ) * r + rhs(i+2,k,j)

  d(i+3,k,j) = ( x(i+3, k  , j-1) &
             +   x(i+3, k  , j+1) ) * r + rhs(i+3,k,j)
end do
end do

return
end subroutine lsor_lu_rhs_j


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
subroutine lsor_lu_rhs_k (d, sz, idx, g, j, k, x, msk, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk
real                                                   ::  r
!dir$ assume_aligned d:64,x:64,msk:64

ist = idx(0)
ied = idx(1)

r = 1.0/6.0

flop = flop + dble(ied-ist+1)*16.0

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i = ist, ied, 4
  d(i  ,k,j) = ( d(i  , k  , j  ) &
             + ( x(i-1, k  , j  ) &
             +   x(i+1, k  , j  ) ) * r ) &
             * msk(i  , k  , j  )

  d(i+1,k,j) = ( d(i+1, k  , j  ) &
             + ( x(i  , k  , j  ) &
             +   x(i+2, k  , j  ) ) * r ) &
             * msk(i+1, k  , j  )

  d(i+2,k,j) = ( d(i+2, k  , j  ) &
             + ( x(i+1, k  , j  ) &
             +   x(i+3, k  , j  ) ) * r ) &
             * msk(i+2, k  , j  )

  d(i+3,k,j) = ( d(i+3, k  , j  ) &
             + ( x(i+2, k  , j  ) &
             +   x(i+4, k  , j  ) ) * r ) &
             * msk(i+3, k  , j  )
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
flop = flop + dble(ied-ist+1)*12.0;

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i=ist, ied, 4
  d(i  , kst, j  ) = ( d(i  , kst  , j  ) &
                   + rhs(i  , kst-1, j  ) * r ) &
                   * msk(i  , kst  , j  )

  d(i+1, kst, j  ) = ( d(i+1, kst  , j  ) &
                   + rhs(i+1, kst-1, j  ) * r ) &
                   * msk(i+1, kst  , j  )

  d(i+2, kst, j  ) = ( d(i+2, kst  , j  ) &
                   + rhs(i+2, kst-1, j  ) * r ) &
                   * msk(i+2, kst  , j  )

  d(i+3, kst, j  ) = ( d(i+3, kst  , j  ) &
                   + rhs(i+3, kst-1, j  ) * r ) &
                   * msk(i+3, kst  , j  )
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
flop = flop + dble(ied-ist+1)*12.0;

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i=ist, ied, 4
  d(i  , ked, j  ) = ( d(i  , ked  , j  ) &
                   + rhs(i  , ked+1, j  ) * r ) &
                   * msk(i  , ked  , j  )

  d(i+1, ked, j  ) = ( d(i+1, ked  , j  ) &
                   + rhs(i+1, ked+1, j  ) * r ) &
                   * msk(i+1, ked  , j  )

  d(i+2, ked, j  ) = ( d(i+2, ked  , j  ) &
                   + rhs(i+2, ked+1, j  ) * r ) &
                   * msk(i+2, ked  , j  )

  d(i+3, ked, j  ) = ( d(i+3, ked  , j  ) &
                   + rhs(i+3, ked+1, j  ) * r ) &
                   * msk(i+3, ked  , j  )
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
!dir$ assume_aligned d:64,a:64,e:64

ist = idx(0)
ied = idx(1)

flop = flop + dble(ied-ist+1)*12.0

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i=ist, ied, 4
  aa = a(k)
  ee = e(k)
  d(i  ,k  ,j  ) = (d(i  ,k  ,j  ) - aa * d(i  , k-1,j  )) * ee
  d(i+1,k  ,j  ) = (d(i+1,k  ,j  ) - aa * d(i+1, k-1,j  )) * ee
  d(i+2,k  ,j  ) = (d(i+2,k  ,j  ) - aa * d(i+2, k-1,j  )) * ee
  d(i+3,k  ,j  ) = (d(i+3,k  ,j  ) - aa * d(i+3, k-1,j  )) * ee
end do

return
end subroutine lsor_tdma_f


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_tdma_b (d, sz, idx, g, j, w, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real                                                   ::  ww
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d
real, dimension(1-g:sz(3)+g)                           ::  w
!dir$ assume_aligned d:64,w:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)
flop = flop + dble((ied-ist+1)*(ked-2))*8.0

do k=ked-1, 2, -1

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i=ist, ied, 4
  ww = w(k)
  d(i  , k  , j  ) = d(i  , k  , j  ) - ww * d(i  , k+1, j  )
  d(i+1, k  , j  ) = d(i+1, k  , j  ) - ww * d(i+1, k+1, j  )
  d(i+2, k  , j  ) = d(i+2, k  , j  ) - ww * d(i+2, k+1, j  )
  d(i+3, k  , j  ) = d(i+3, k  , j  ) - ww * d(i+3, k+1, j  )
end do
end do

return
end subroutine lsor_tdma_b


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
subroutine lsor_relax (d, sz, idx, g, j, x, msk, omg, res, flop)
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

flop = flop + dble((ked-kst+1)*(ied-ist+1))*24.0

do k = kst, ked

!dir$ vector always
!dir$ vector aligned
!dir$ ivdep
do i = ist, ied, 4
  pp1 =   x(i  , k  , j  )
  dp1 = ( d(i  , k  , j  ) - pp1 ) * omg * msk(i  , k  , j  )
  pn1 = pp1 + dp1
  x(i  , k  , j  ) = pn1

  pp2 =   x(i+1, k  , j  )
  dp2 = ( d(i+1, k  , j  ) - pp2 ) * omg * msk(i+1, k  , j  )
  pn2 = pp2 + dp2
  x(i+1, k  , j  ) = pn2

  pp3 =   x(i+2, k  , j  )
  dp3 = ( d(i+2, k  , j  ) - pp3 ) * omg * msk(i+2, k  , j  )
  pn3 = pp3 + dp3
  x(i+2, k  , j  ) = pn3

  pp4 =   x(i+3, k  , j  )
  dp4 = ( d(i+3, k  , j  ) - pp4 ) * omg * msk(i+3, k  , j  )
  pn4 = pp4 + dp4
  x(i+3, k  , j  ) = pn4

  res = res + dp1 * dp1 &
            + dp2 * dp2 &
            + dp3 * dp3 &
            + dp4 * dp4
end do
end do

return
end subroutine lsor_relax


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
subroutine tdma_lsor_b (d, sz, idx, g, x, w, a, b, c, rhs, omg, res, flop)
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
real, dimension(1-g:sz(1)+g)                           ::  a, b, c

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
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn, e, y)
do k = kst, ked
do j = jst, jed

  !$DIR IVDEP
  do i = ist, ied
    d(i,j,k) = ( x(i  ,j-1,k  ) &
             +   x(i  ,j+1,k  ) &
             +   x(i  ,j  ,k-1) &
             +   x(i  ,j  ,k+1) ) * r + rhs(i,j,k)
  end do

  d(ist,j,k) = d(ist,j,k) + rhs(ist-1,j,k)*r
  d(ied,j,k) = d(ied,j,k) + rhs(ied+1,j,k)*r

  ! TDMA
  y = 1.0 / b(ist)
  d(ist,j,k) = d(ist,j,k) * y
  w(ist,j,k) = c(ist) * y

  !$DIR VECTOR
  do i = ist, ied
    e = 1.0 / (b(i) - a(i) * w(i-1,j,k))
    w(i,j,k) = e * c(i)
    d(i,j,k) = (d(i,j,k) - a(i) * d(i-1,j,k)) * e
  end do

  do i=ied-1, ist, -1
    d(i,j,k) = d(i,j,k) - w(i,j,k) * d(i+1,j,k)
  end do


  !$DIR IVDEP
  !do k = kst, ked
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
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine tdma_lsor_e (d, sz, idx, g, x, w, a, b, c, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x, w, rhs
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
     d(k,i,j) = ( x(k,i-1,j  ) &
              +   x(k,i+1,j  ) &
              +   x(k,i  ,j-1) &
              +   x(k,i  ,j+1) ) * r + rhs(k,i,j)
  end do
  d(kst,i,j) = d(kst,i,j) + rhs(kst-1,i,j)*r
  d(ked,i,j) = d(ked,i,j) + rhs(ked+1,i,j)*r

  ! TDMA
  d(kst,i,j) = d(kst,i,j) * y
  w(kst,i,j) = c(kst) * y

  do k=kst+1, ked
    e = 1.0 / (b(k) - a(k) * w(k-1,i,j))
    w(k,i,j) = e * c(k)
    d(k,i,j) = (d(k,i,j) - a(k) * d(k-1,i,j)) * e
  end do

  do k=ked-1, kst, -1
    d(k,i,j) = d(k,i,j) - w(k,i,j) * d(k+1,i,j)
  end do


  !$DIR IVDEP
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
end subroutine tdma_lsor_e

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
subroutine tdma_lsor_f (d, sz, idx, g, x, w, e, a, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x, rhs
real                                                   ::  r, y
real, dimension(1-g:sz(3)+g)                           ::  a, w, e

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(5.0+5.0+5.0) &
            + dble((ied-ist+1)*(jed-jst+1)) * 5.0;

y = 1.0 / 1.0 ! 1.0/b(kst)

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn)
do j = jst, jed
do i = ist, ied

  !$DIR IVDEP
  do k = kst, ked
     d(k,i,j) = ( x(k,i-1,j  ) &
              +   x(k,i+1,j  ) &
              +   x(k,i  ,j-1) &
              +   x(k,i  ,j+1) ) * r + rhs(k,i,j)
  end do
  d(kst,i,j) = d(kst,i,j) + rhs(kst-1,i,j)*r
  d(ked,i,j) = d(ked,i,j) + rhs(ked+1,i,j)*r

  ! TDMA
  d(kst,i,j) = d(kst,i,j) * y

  !$DIR VECTOR
  do k=kst+1, ked
    d(k,i,j) = (d(k,i,j) - a(k) * d(k-1,i,j)) * e(k)
  end do

  do k=ked-1, kst, -1
    d(k,i,j) = d(k,i,j) - w(k) * d(k+1,i,j)
  end do


  !$DIR IVDEP
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
end subroutine tdma_lsor_f
