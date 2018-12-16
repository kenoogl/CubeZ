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
