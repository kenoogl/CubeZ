


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
do j = jst, jed
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
!dir$ assume_aligned d:64,x:64

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
!! @param [in,out] e    work
!! @param [in]     sz   配列長
!! @param [in]     idx  インデクス範囲
!! @param [in]     g    ガイドセル長
!! @param [in]     w    work
!! @param [in]     a    coef
!! @param [in]     b    coef
!! @param [in]     c    coef
!<
subroutine ljcb_g0 (e, sz, idx, g, w, a, b, c, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g)                           ::  a, b, c, w, e
real                                                   ::  f


ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble(ked-kst+1)*11.0

w(kst) = c(kst)/b(kst)

do k=kst+1, ked
  f = 1.0 / (b(k) - a(k) * w(k-1))
  w(k) = f * c(k)
  e(k) = f
end do

return
end subroutine ljcb_g0

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
subroutine ljcb_g1 (d, sz, idx, g, b, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  d, rhs
real                                                   ::  r, y
real, dimension(1-g:sz(3)+g)                           ::  b

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)) * 5.0;

y = 1.0 / b(kst)

!$OMP PARALLEL DO SCHEDULE(static)
do j = jst, jed
!$OMP simd
do i = ist, ied
  d(i,j,kst) = ( d(i,j,kst) + rhs(i,j,kst-1)*r ) * y
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
!$OMP simd
do i = ist, ied
  d(i,j,k) = d(i,j,k) - w(k) * d(i,j,k+1)
end do
end do
!$OMP END PARALLEL DO

end do

return
end subroutine ljcb_g3


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
