
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
subroutine ljcb_f0 (d, sz, idx, g, x, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x, rhs
real                                                   ::  r
!dir$ assume_aligned d:64, x:64, rhs:64

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
do j = jst, jed
!dir$ vector aligned
!$OMP simd
do i = ist, ied
  d(i,j,k) = ( x(i-1,j  ,k  ) &
           +   x(i+1,j  ,k  ) &
           +   x(i  ,j-1,k  ) &
           +   x(i  ,j+1,k  ) ) * r + rhs(i,j,k)
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f0



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
subroutine ljcb_f4 (d, sz, idx, g, x, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, x
!dir$ assume_aligned d:64, x:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*5.0

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp)
do k = kst, ked
do j = jst, jed
!dir$ vector aligned
!$OMP simd
do i = ist, ied
  pp = x(i,j,k)
  dp = ( d(i,j,k) - pp ) * omg
  x(i,j,k) = pp + dp
  res = res + dp*dp
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine ljcb_f4



!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     b    coef
!! @param [in]     rhs  オリジナルの線形方程式の右辺項
!! @param [in,out] flop flop count
!<
subroutine ljcb_g1 (d, sz, idx, g, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, rhs
real                                                   ::  r, y
!dir$ assume_aligned d:64, rhs:64


ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)) * 5.0;

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
!dir$ vector aligned
!$OMP simd
do i = ist, ied
  d(i,j,kst) = d(i,j,kst) + rhs(i,j,kst-1)*r
  d(i,j,ked) = d(i,j,ked) + rhs(i,j,ked+1)*r
end do
end do
!$OMP END PARALLEL DO


return
end subroutine ljcb_g1

!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     e    work
!! @param [in,out] flop flop count
!<
subroutine ljcb_g2 (d, sz, idx, g, e, a, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d
real, dimension(1-g:sz(3)+g)                           ::  e, a
!dir$ assume_aligned d:64, a:64, e:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*3.0


do k=kst+1, ked

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
!dir$ vector aligned
!$OMP simd
do i = ist, ied
  d(i,j,k) = (d(i,j,k) - a(k) * d(i,j,k-1)) * e(k)
end do
end do
!$OMP END PARALLEL DO

end do

return
end subroutine ljcb_g2

!> ********************************************************************
!! @brief lsor
!! @param [in,out] d    ソース項
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in,out] flop flop count
!<
subroutine ljcb_g3 (d, sz, idx, g, w, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d
real, dimension(1-g:sz(3)+g)                           ::  w
!dir$ assume_aligned d:64, w:64

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
!dir$ vector aligned
do i = ist, ied
  d(i,j,k) = d(i,j,k) - w(k) * d(i,j,k+1)
end do
end do
!$OMP END PARALLEL DO

end do

return
end subroutine ljcb_g3
