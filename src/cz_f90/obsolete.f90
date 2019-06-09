

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
