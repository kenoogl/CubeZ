!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################

!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_pcr_kij2 (sz, idx, g, pn, x, a, c, d, a1, c1, d1, msk, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s, p, pn
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs, d, d1
real, dimension(1-g:sz(3)+g)                           ::  a, c, a1, c1
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
+ pn * 21.0 &  ! PCR
+ 6.0 )     &  ! Relaxation
+ 6.0                   &   ! BC
) )


res = 0.0

!$OMP PARALLEL DO SCHEDULE(static) collapse(2) &
!$OMP          private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP          reduction(+:res)
do j=jst, jed
do i=ist, ied

  ! Reflesh coef. due to override

  !dir$ vector aligned
  !dir$ simd
  do k=kst+1, ked
    a(k) = -r
  end do

  !dir$ vector aligned
  !dir$ simd
  do k=kst, ked-1
    c(k) = -r
  end do

  a(kst) = 0.0
  c(ked) = 0.0


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
      ap = a(k)
      cp = c(k)
      e = 1.0 / ( 1.0 - ap * c(kl) - cp * a(kr) )
      a1(k) =  -e * ap * a(kl)
      c1(k) =  -e * cp * c(kr)
      d1(k,i,j) =   e * ( d(k,i,j) - ap * d(kl,i,j) - cp * d(kr,i,j) )
    end do

    !dir$ vector aligned
    !dir$ simd
    do k = kst, ked
      a(k) = a1(k)
      c(k) = c1(k)
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
    res = res + dp * dp
  end do

end do
end do

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
                             + pn * 21.0 &  ! PCR
                             + 6.0 )     &  ! Relaxation
                 + 6.0                   &  ! BC
              ) )

! Reflesh coef. due to override
!$OMP PARALLEL DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst+1, ked
  a(k,i,j) = -r
end do
end do
end do

!$OMP PARALLEL DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
do k=kst, ked-1
  c(k,i,j) = -r
end do
end do
end do

!$OMP PARALLEL DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
  a(kst,i,j) = 0.0
end do
end do

!$OMP PARALLEL DO SCHEDULE(static) collapse(2)
do j=jst, jed
do i=ist, ied
  c(ked,i,j) = 0.0
end do
end do


res = 0.0


!$OMP PARALLEL DO SCHEDULE(static) collapse(2) &
!$OMP          private(kl, kr, ap, cp, e, s, p, k, pp, dp) &
!$OMP          reduction(+:res)
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
    res = res + dp * dp
  end do

end do
end do

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
!! @note resDのFalse sharering回避
!<
subroutine lsor_pcr_src_q (sz, idx, g, i, j, x, d, msk, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, rhs, d
real                                                   ::  r
!dir$ assume_aligned x:64, msk:64, rhs:64, d:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble( (ked-kst+1)*6.0 + 6.0 )


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

return
end subroutine lsor_pcr_src_q


!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note resDのFalse sharering回避
!<
subroutine lsor_pcr_q (sz, idx, g, s, i, j, a, c, d, a1, c1, d1, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  a, c, d, a1, c1, d1
real                                                   ::  r, ap, cp, e
!dir$ assume_aligned a:64, c:64, d:64, a1:64, c1:64, d1:64

!ist = idx(0)
!ied = idx(1)
!jst = idx(2)
!jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble( (ked-kst+1)*21.0 )


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

return
end subroutine lsor_pcr_q


!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_pcr_relax_q (sz, idx, g, i, j, x, d, msk, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  x, msk, d
real                                                   ::  omg, dp, pp
!dir$ assume_aligned x:64, msk:64, d:64

!ist = idx(0)
!ied = idx(1)
!jst = idx(2)
!jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble( (ked-kst+1)*6.0 )


! !dir$ vector aligned
! !dir$ simd

!$OMP SIMD reduction(+:res)
do k = kst, ked
  pp =   x(k, i, j)
  dp = ( d(k, i, j) - pp ) * omg * msk(k, i, j)
  x(k, i, j) = pp + dp
  res = res + dp * dp
end do

return
end subroutine lsor_pcr_relax_q


!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note resDのFalse sharering回避
!<
subroutine lsor_pcr_src (sz, idx, g, x, d, msk, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x, msk, rhs, d
real                                                   ::  r
!dir$ assume_aligned x:64, msk:64, rhs:64, d:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble( (ked-kst+1)*(jed-jst+1)*(ied-ist+1)*6.0 &
                  + (jed-jst+1)*(ied-ist+1)*6.0 )


!$OMP PARALLEL DO SCHEDULE(static)
do k = kst, ked
do j = jst, jed
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i, j, k) = ( ( x(i  , j-1, k)        &
             +     x(i  , j+1, k)        &
             +     x(i-1, j  , k)        &
             +     x(i+1, j  , k) ) * r + rhs(i, j, k) ) &
             *   msk(i, j, k)
end do
end do
end do


! BC

!$OMP PARALLEL DO SCHEDULE(static)
do j=jst, jed
!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, j, kst) = ( d(i, j, kst) + rhs(i, j, kst-1) * r ) * msk(i, j, kst)
  d(i, j, ked) = ( d(i, j, ked) + rhs(i, j, ked+1) * r ) * msk(i, j, ked)
end do
end do

return
end subroutine lsor_pcr_src


!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note resDのFalse sharering回避
!<
subroutine lsor_pcr (sz, idx, g, s, a, c, d, a1, c1, d1, flop)
implicit none
integer                                                ::  i, j, k, g, kl, kr, s
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  a, c, d, a1, c1, d1
real                                                   ::  r, ap, cp, e
!dir$ assume_aligned a:64, c:64, d:64, a1:64, c1:64, d1:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble( (ked-kst+1)*(jed-jst+1)*(ied-ist+1)*21.0 )


!$OMP PARALLEL DO SCHEDULE(static) private(kl, kr, ap, cp, e)
do k = kst, ked
do j = jst, jed
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  kl = max(k-s, kst-1)
  kr = min(k+s, ked+1)
  ap = a(i,j,k)
  cp = c(i,j,k)
  e = 1.0 / ( 1.0 - ap * c(i,j,kl) - cp * a(i,j,kr) )
  a1(i,j,k) =  -e * ap * a(i,j,kl)
  c1(i,j,k) =  -e * cp * c(i,j,kr)
  d1(i,j,k) =   e * ( d(i,j,k) - ap * d(i,j,kl) - cp * d(i,j,kr) )
  !write (*,*) i,j,k, s, kl,kr, d(i,j,kl),d(i,j,kr)
end do
end do
end do

return
end subroutine lsor_pcr


!> ********************************************************************
!! @brief PCR
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_pcr_relax (sz, idx, g, x, d, msk, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, jst, jed, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  x, msk, d
real                                                   ::  omg, dp, pp
!dir$ assume_aligned x:64, msk:64, d:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

flop = flop + dble( (ked-kst+1)*(jed-jst+1)*(ied-ist+1)*6.0 )


!$OMP PARALLEL DO SCHEDULE(static) private(dp,pp) reduction(+:res)
do k = kst, ked
do j = jst, jed
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  pp =   x(i, j, k)
  dp = ( d(i, j, k) - pp ) * omg * msk(i, j, k)
  x(i, j, k) = pp + dp
  res = res + dp * dp
  !write (*,*) i,j,k, x(i,j,k), msk(i,j,k)
end do
end do
end do

return
end subroutine lsor_pcr_relax


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note resDのFalse sharering回避
!<
subroutine lsor_inner_c (d, sz, idx, g, d2, x, a, e, w, msk, rhs, omg, &
                         clsz, nt, itri, rd, flop)
!$ use omp_lib
implicit none
integer                                                ::  i, j, k, g, clsz, l, nt, itri
integer                                                ::  ist, ied, kst, ked, jst, jed, id
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
double precision, dimension(clsz, nt)                  ::  rd
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
real                                                   ::  omg, r
real                                                   ::  dp, pp
!dir$ assume_aligned d:64, a:64, e:64, w:64, d2:64, x:64, msk:64, rhs:64, rd:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble(jed-jst+1) * ( &
       dble((ied-ist+1)*(ked-kst+1))*3.0 &
     + dble(itri)* ( &
       dble((ied-ist+1)*(ked-kst+1))*(4.0+6.0) &
      +dble(ied-ist+1)*(3.0+3.0) &
      +dble((ied-ist+1)*(ked-2))*(3.0+2.0) ) )

!$OMP PARALLEL DO SCHEDULE(dynamic,1) private(dp,pp,id) reduction(+:rd)
do j=jst, jed
  id = omp_get_thread_num() + 1

  ! Source term I for AX=b

  do k = kst, ked
    !dir$ vector aligned
    !dir$ simd
    do i = ist, ied
      d2(i, k, j) = ( x(i, k, j-1) &
                  +   x(i, k, j+1) ) * r + rhs(i, k, j)
    end do
  end do


  ! Inner Iteration
  do l=1, itri

    ! Source term II for AX=b

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


    ! BC

    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, kst, j) = ( d(i, kst  , j) &
                   + rhs(i, kst-1, j) * r ) &
                   * msk(i, kst  , j)
    end do


    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, ked, j) = ( d(i, ked  , j) &
                   + rhs(i, ked+1, j) * r ) &
                   * msk(i, ked  , j)
    end do

    ! Forward

    do k = 3, ked
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
      end do
    end do


    ! Backwad

    do k=ked-1, 2 , -1
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
      end do
    end do


    ! Relax

    do k = kst, ked
      !dir$ vector aligned
      !dir$ simd
      do i = ist, ied
        pp =   x(i, k  , j)
        dp = ( d(i, k  , j) - pp ) * omg * msk(i, k  , j)
        x(i, k, j) = pp + dp
        rd(l,id) = rd(l,id) + dp * dp
      end do
    end do

  end do ! ItrInner
end do ! j-loop

!do l=1,itri
!  do i=1,nt
!    write(*,*) l, i, rd(l,i)
!  end do
!end do

return
end subroutine lsor_inner_c

!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_inner_rb (d, sz, idx, g, d2, x, a, e, w, &
                         msk, rhs, omg, itri, rd, flop)
implicit none
integer                                                ::  i, j, k, g, itri, l, col
integer                                                ::  ist, ied, kst, ked, jst, jed
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
double precision, dimension(itri)                      ::  rd
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
real                                                   ::  omg, r
real                                                   ::  dp1, pp1
!dir$ assume_aligned d:64, a:64, e:64, w:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble(jed-jst+1) * ( &
       dble((ied-ist+1)*(ked-kst+1))*3.0 &
     + dble(itri)* ( &
       dble((ied-ist+1)*(ked-kst+1))*(4.0+6.0) &
      +dble(ied-ist+1)*(3.0+3.0) &
      +dble((ied-ist+1)*(ked-2))*(3.0+2.0) ) )

do col=0, 1

!$OMP PARALLEL DO SCHEDULE(static) reduction(+:rd) &
!$OMP private(dp1,pp1)
do j=jst+col, jed, 2

  ! Source term I for AX=b

  ! Intel 19.0 compiler can not unroll
  do k = kst, ked
    !dir$ vector aligned
    !dir$ simd
    do i = ist, ied
      d2(i, k  , j) = ( x(i, k  , j-1) &
                    +   x(i, k  , j+1) ) * r + rhs(i, k  , j)
    end do
  end do


  ! Inner Iteration
  do l=1, itri

    ! Source term II for AX=b

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


    ! BC

    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, kst, j) = ( d(i, kst  , j) &
                   + rhs(i, kst-1, j) * r ) &
                   * msk(i, kst  , j)
    end do


    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, ked, j) = ( d(i, ked  , j) &
                   + rhs(i, ked+1, j) * r ) &
                   * msk(i, ked  , j)
    end do

    ! Forward

    do k = 3, ked
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
      end do
    end do


    ! Backwad

    do k=ked-1, 2 , -1
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
      end do
    end do


    ! Relax

    do k = kst, ked
      !dir$ vector aligned
      !dir$ simd
      do i = ist, ied
        pp1 =   x(i, k  , j)
        dp1 = ( d(i, k  , j) - pp1 ) * omg * msk(i, k  , j)
        x(i, k  , j) = pp1 + dp1
        rd(l) = rd(l) + dp1 * dp1
      end do
    end do

  end do ! ItrInner
end do ! j-loop

end do ! col

return
end subroutine lsor_inner_rb

!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_inner_cb (d, sz, idx, g, d2, x, a, e, w, &
                         msk, rhs, omg, itri, rd, flop)
implicit none
integer                                                ::  i, j, k, g, itri, l
integer                                                ::  ti, tk, ldvi, ldvk, ib, kb
integer                                                ::  ist, ied, kst, ked, jst, jed
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
double precision, dimension(itri)                      ::  rd
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
real                                                   ::  omg, r
real                                                   ::  dp, pp
!dir$ assume_aligned d:64, a:64, e:64, w:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

ldvi = 2
ldvk = 80
ti = (ied-ist+1)/ldvi
if ( ied-ist+1 /= ti*ldvi ) ti = ti + 1
tk = (ked-kst+1)/ldvk
if ( ked-kst+1 /= tk*ldvk ) tk = tk + 1
!write(*,*) ti, tk, 5*ti*tk*3*4

r = 1.0/6.0

flop = flop + dble(jed-jst+1) * ( &
       dble((ied-ist+1)*(ked-kst+1))*3.0 &
     + dble(itri)* ( &
       dble((ied-ist+1)*(ked-kst+1))*(4.0+6.0) &
      +dble(ied-ist+1)*(3.0+3.0) &
      +dble((ied-ist+1)*(ked-2))*(3.0+2.0) ) )

! Cache blocking, I:128, J:3, K:4 x 5 vars
!       128 * 3 * 4 * 5 * 4B = 30 kB

!$OMP PARALLEL DO SCHEDULE(static,3) private(dp, pp) reduction(+:rd)
do j=jst, jed

  ! Source term I for Ax=b
  do kb=kst, ked, tk
  do ib=ist, ied, ti

    do k = kb, min(kb+tk-1, ked)
      !dir$ vector aligned
      !dir$ simd
      do i = ib, min(ib+ti-1, ied)
        d2(i, k, j) = ( x(i, k, j-1) &
                    +   x(i, k, j+1) ) * r + rhs(i, k, j)
      end do
    end do
  end do
  end do


  ! Inner Iteration
  do l=1, itri

    ! Source term II for AX=b
    do kb=kst, ked, tk
    do ib=ist, ied, ti

      do k = kb, min(kb+tk-1, ked)
        !dir$ vector aligned
        !dir$ simd
        do i = ib, min(ib+ti-1, ied)
          d(i, k  , j) = (d2(i  , k  , j) &
                       + ( x(i-1, k  , j) &
                       +   x(i+1, k  , j) ) * r ) &
                       * msk(i  , k  , j)
        end do
      end do
    end do
    end do



    ! BC

    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, kst, j) = ( d(i, kst  , j) &
                   + rhs(i, kst-1, j) * r ) &
                   * msk(i, kst  , j)
    end do


    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, ked, j) = ( d(i, ked  , j) &
                   + rhs(i, ked+1, j) * r ) &
                   * msk(i, ked  , j)
    end do

    ! Forward

    do kb=3, ked, tk
    do ib=ist, ied, ti

    do k = kb, min(kb+tk-1, ked)
      !dir$ vector aligned
      !dir$ simd
      do i = ib, min(ib+ti-1, ied)
        d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
      end do
    end do

    end do
    end do


    ! Backwad
    do kb=ked-1, 2, -tk
    do ib=ist, ied, ti

    do k=kb, max(kb-tk+1, 2), -1
      !dir$ vector aligned
      !dir$ simd
      do i = ib, min(ib+ti-1, ied)
        d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
      end do
    end do

    end do
    end do

    ! Relax
    do kb=kst, ked, tk
    do ib=ist, ied, ti

    do k = kb, min(kb+tk-1, ked)
      !dir$ vector aligned
      !dir$ simd
      do i = ib, min(ib+ti-1, ied)
        pp =   x(i, k  , j)
        dp = ( d(i, k  , j) - pp ) * omg * msk(i, k  , j)
        x(i, k, j) = pp + dp
        rd(l) = rd(l) + dp * dp
      end do
    end do

    end do
    end do

  end do ! ItrInner
end do ! j-loop

return
end subroutine lsor_inner_cb


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_inner_b (d, sz, idx, g, d2, x, a, e, w, &
                         msk, rhs, omg, itri, rd, flop)
implicit none
integer                                                ::  i, j, k, g, itri, l
integer                                                ::  ist, ied, kst, ked, jst, jed
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
double precision, dimension(itri)                      ::  rd
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
real                                                   ::  omg, r
real                                                   ::  dp, pp
!dir$ assume_aligned d:64, a:64, e:64, w:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble(jed-jst+1) * ( &
       dble((ied-ist+1)*(ked-kst+1))*3.0 &
     + dble(itri)* ( &
       dble((ied-ist+1)*(ked-kst+1))*(4.0+6.0) &
      +dble(ied-ist+1)*(3.0+3.0) &
      +dble((ied-ist+1)*(ked-2))*(3.0+2.0) ) )

!$OMP PARALLEL DO SCHEDULE(static,1) private(dp, pp) reduction(+:rd)
do j=jst, jed

  ! Source term I for AX=b

  do k = kst, ked
    !dir$ vector aligned
    !dir$ simd
    do i = ist, ied
      d2(i, k, j) = ( x(i, k, j-1) &
                  +   x(i, k, j+1) ) * r + rhs(i, k, j)
    end do
  end do


  ! Inner Iteration
  do l=1, itri

    ! Source term II for AX=b

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


    ! BC

    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, kst, j) = ( d(i, kst  , j) &
                   + rhs(i, kst-1, j) * r ) &
                   * msk(i, kst  , j)
    end do


    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, ked, j) = ( d(i, ked  , j) &
                   + rhs(i, ked+1, j) * r ) &
                   * msk(i, ked  , j)
    end do

    ! Forward

    do k = 3, ked
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
      end do
    end do


    ! Backwad

    do k=ked-1, 2 , -1
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
      end do
    end do


    ! Relax

    do k = kst, ked
      !dir$ vector aligned
      !dir$ simd
      do i = ist, ied
        pp =   x(i, k  , j)
        dp = ( d(i, k  , j) - pp ) * omg * msk(i, k  , j)
        x(i, k, j) = pp + dp
        rd(l) = rd(l) + dp * dp
      end do
    end do

  end do ! ItrInner
end do ! j-loop

return
end subroutine lsor_inner_b


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!! @note 予測したdのループを繰り返す
!<
subroutine lsor_inner_d (d, sz, idx, g, d2, x, a, e, w, &
                         msk, rhs, omg, itri, rd, flop)
implicit none
integer                                                ::  i, j, k, g, itri, l
integer                                                ::  ist, ied, kst, ked, jst, jed
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
double precision, dimension(itri)                      ::  rd
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
real                                                   ::  omg, r
real                                                   ::  dp, pp
!dir$ assume_aligned d:64, a:64, e:64, w:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble(jed-jst+1) * ( &
       dble((ied-ist+1)*(ked-kst+1))*3.0 &
     + dble(itri)* ( &
       dble((ied-ist+1)*(ked-kst+1))*(4.0+6.0) &
      +dble(ied-ist+1)*(3.0+3.0) &
      +dble((ied-ist+1)*(ked-2))*(3.0+2.0) ) )

!$OMP PARALLEL DO SCHEDULE(static,1) private(dp, pp) reduction(+:rd)
do j=jst, jed

  ! Source term I for AX=b

  do k = kst, ked
    !dir$ vector aligned
    !dir$ simd
    do i = ist, ied
      d2(i, k, j) = ( x(i, k, j-1) &
                  +   x(i, k, j+1) ) * r + rhs(i, k, j)
    end do
  end do


  ! Inner Iteration
  do l=1, itri

    ! Source term II for AX=b

    if (l == 1) then
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
    else
      do k = kst, ked
        !dir$ vector aligned
        !dir$ simd
        do i = ist, ied
          d(i, k  , j) = (d2(i  , k  , j) &
                       + ( d(i-1, k  , j) &
                       +   d(i+1, k  , j) ) * r ) &
                       * msk(i  , k  , j)
        end do
      end do
    endif


    ! BC

    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, kst, j) = ( d(i, kst  , j) &
                   + rhs(i, kst-1, j) * r ) &
                   * msk(i, kst  , j)
    end do


    !dir$ vector aligned
    !dir$ simd
    do i=ist, ied
      d(i, ked, j) = ( d(i, ked  , j) &
                   + rhs(i, ked+1, j) * r ) &
                   * msk(i, ked  , j)
    end do

    ! Forward

    do k = 3, ked
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = (d(i, k  , j) - a(k) * d(i, k-1, j)) * e(k)
      end do
    end do


    ! Backwad

    do k=ked-1, 2 , -1
      !dir$ vector aligned
      !dir$ simd
      do i=ist, ied
        d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
      end do
    end do

  end do ! ItrInner

    ! Relax

  do k = kst, ked
    !dir$ vector aligned
    !dir$ simd
    do i = ist, ied
      pp =   x(i, k  , j)
      dp = ( d(i, k  , j) - pp ) * omg * msk(i, k  , j)
      x(i, k, j) = pp + dp
      rd(1) = rd(1) + dp * dp
    end do
  end do

end do ! j-loop

return
end subroutine lsor_inner_d


!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_inner (d, sz, idx, g, j, d2, x, a, e, w, msk, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, f1, f2, res
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, d2, x, msk, rhs
real, dimension(1-g:sz(3)+g)                           ::  a, e, w
real                                                   ::  omg, r
real                                                   ::  dp1, pp1
real                                                   ::  dp2, pp2
real                                                   ::  dp3, pp3
real                                                   ::  dp4, pp4
!dir$ assume_aligned d:64, a:64, e:64, w:64, d2:64, x:64, msk:64, rhs:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0
f1 = dble( (ied-ist+1)*(ked-kst+1) )
f2 = dble( (ied-ist+1)*(ked-2) )

flop = flop + f1 * 4.0

! Source term for AX=b

!dir$ unroll(4)
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


! BC

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, kst, j) = ( d(i, kst  , j) &
               + rhs(i, kst-1, j) * r ) &
               * msk(i, kst  , j)
end do

flop = flop + f1 * 3.0


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, ked, j) = ( d(i, ked  , j) &
               + rhs(i, ked+1, j) * r ) &
               * msk(i, ked  , j)
end do

flop = flop + f1 * 3.0



! Forward

do k = 3, ked
!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k  , j) = (d(i, k  , j) - a(k  ) * d(i, k-1, j)) * e(k  )
end do
end do

flop = flop + f2 * 3.0


! Backwad

do k=ked-1, 2 , -1
!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
end do
end do

flop = flop + f2 * 2.0


! Relax

do k = kst, ked
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  pp1 =   x(i, k  , j)
  dp1 = ( d(i, k  , j) - pp1 ) * omg * msk(i, k  , j)
  x(i, k, j) = pp1 + dp1
  res = res + dp1 * dp1
end do
end do

flop = flop + f1 * 6.0

return
end subroutine lsor_inner




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
!dir$ assume_aligned d:64, x:64, rhs:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(ked-kst+1))*3.0

!dir$ unroll(4)
do k = kst, ked
!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i, k, j) = ( x(i, k, j-1) &
             +   x(i, k, j+1) ) * r + rhs(i, k, j)
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
subroutine lsor_lu_fwd (d, sz, idx, g, j, x, a, e, d2, msk, rhs, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied, kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk, d2, rhs
real                                                   ::  r
real, dimension(1-g:sz(3)+g)                           ::  a, e
!dir$ assume_aligned d:64, x:64, msk:64, d2:64, rhs:64, a:64, e:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble( (ied-ist+1)*(ked-kst+1) )*4.0

do k = kst, ked

!dir$ vector aligned
!dir$ simd
do i = ist, ied
  d(i  ,k, j) = (d2(i  , k, j  ) &
              + ( x(i-1, k, j  ) &
              +   x(i+1, k, j  ) ) * r ) &
              * msk(i  , k, j  )
end do
end do


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, kst, j) = ( d(i, kst  , j) &
               + rhs(i, kst-1, j) * r ) &
               * msk(i, kst  , j)
end do

flop = flop + dble( (ied-ist+1) )*3.0


!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, ked, j) = ( d(i, ked  , j) &
               + rhs(i, ked+1, j) * r ) &
               * msk(i, ked  , j)
end do

flop = flop + dble( (ied-ist+1) )*3.0


do k = 3, ked

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k, j  ) = (d(i, k, j  ) - a(k) * d(i, k-1, j  )) * e(k)
end do
end do

flop = flop + dble( (ied-ist+1)*(ked-2) )*3.0


return
end subroutine lsor_lu_fwd



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
!dir$ assume_aligned d:64, w:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)
flop = flop + dble((ied-ist+1)*(ked-2))*2.0

do k=ked-1, 2, -1

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i  , k  , j  ) = d(i  , k  , j  ) - w(k) * d(i  , k+1, j  )
end do
end do

return
end subroutine lsor_tdma_b

!> ********************************************************************
!! @brief TDMA
!! @param [in]     nx   配列長
!! @param [in]     g    ガイドセル長
!! @param [in,out] d    RHS vector -> 解ベクトル (in-place)
!! @param [in]     cf   係数
!! @param [in]     w    U_1 vector
!! @param [in,out] flop flop count
!<
subroutine lsor_tdma_b4 (d, sz, idx, g, j, w, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d
real, dimension(1-g:sz(3)+g)                           ::  w
!dir$ assume_aligned d:64, w:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)
flop = flop + dble((ied-ist+1)*(ked-2))*2.0

do k=ked-1, 2, -4

!dir$ vector aligned
!dir$ simd
do i=ist, ied
  d(i, k  , j) = d(i, k  , j) - w(k  ) * d(i, k+1, j)
  d(i, k-1, j) = d(i, k-1, j) - w(k-1) * d(i, k  , j)
  d(i, k-2, j) = d(i, k-2, j) - w(k-2) * d(i, k-1, j)
  d(i, k-3, j) = d(i, k-3, j) - w(k-3) * d(i, k-2, j)
end do
end do

return
end subroutine lsor_tdma_b4


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
real                                                   ::  dp1, pp1
real                                                   ::  dp2, pp2
real                                                   ::  dp3, pp3
real                                                   ::  dp4, pp4
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk
!dir$ assume_aligned d:64, x:64, msk:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ked-kst+1)*(ied-ist+1))*6.0

do k = kst, ked

!dir$ vector aligned
!dir$ simd
do i = ist, ied
  pp1 =   x(i  , k  , j  )
  dp1 = ( d(i  , k  , j  ) - pp1 ) * omg * msk(i  , k  , j  )
  x(i  , k  , j  ) = pp1 + dp1

  res = res + dp1 * dp1
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
!! @param [in]     omg  加速係数
!! @param [out]    res  残差
!! @param [in,out] flop flop count
!<
subroutine lsor_relax4 (d, sz, idx, g, j, x, msk, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, ied
integer                                                ::  kst, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg
real                                                   ::  dp1, pp1
real                                                   ::  dp2, pp2
real                                                   ::  dp3, pp3
real                                                   ::  dp4, pp4
real, dimension(1-g:sz(1)+g, 1-g:sz(3)+g, 1-g:sz(2)+g) ::  d, x, msk
!dir$ assume_aligned d:64, x:64, msk:64

ist = idx(0)
ied = idx(1)
kst = idx(4)
ked = idx(5)

flop = flop + dble((ked-kst+1)*(ied-ist+1))*6.0

do k = kst, ked, 4

!dir$ vector aligned
!dir$ simd
do i = ist, ied
  pp1 =   x(i, k  , j)
  dp1 = ( d(i, k  , j) - pp1 ) * omg * msk(i, k  , j)
  x(i, k, j) = pp1 + dp1

  pp2 =   x(i, k+1, j)
  dp2 = ( d(i, k+1, j) - pp2 ) * omg * msk(i, k+1, j)
  x(i, k+1, j) = pp2 + dp2

  pp3 =   x(i, k+2, j)
  dp3 = ( d(i, k+2, j) - pp3 ) * omg * msk(i, k+2, j)
  x(i, k+2, j) = pp3 + dp3

  pp4 =   x(i, k+3, j)
  dp4 = ( d(i, k+3, j) - pp4 ) * omg * msk(i, k+3, j)
  x(i, k+3, j) = pp4 + dp4

  res = res + dp1 * dp1 &
            + dp2 * dp2 &
            + dp3 * dp3 &
            + dp4 * dp4
end do
end do

return
end subroutine lsor_relax4




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
subroutine tdma_lsor_e (d, sz, idx, g, x, w, a, c, rhs, omg, res, flop)
implicit none
integer                                                ::  i, j, k, g
integer                                                ::  ist, jst, kst
integer                                                ::  ied, jed, ked
integer, dimension(3)                                  ::  sz
integer, dimension(0:5)                                ::  idx
double precision                                       ::  flop, res
real                                                   ::  omg, dp, pp, pn
real, dimension(1-g:sz(3)+g, 1-g:sz(1)+g, 1-g:sz(2)+g) ::  d, x, w, rhs
real                                                   ::  r, e
real, dimension(1-g:sz(3)+g)                           ::  a, c
!dir$ assume_aligned d:64, x:64, w:64, rhs:64, a:64, c:64


ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(4.0+16.0+5.0) &
            + dble((ied-ist+1)*(jed-jst+1)) * 2.0;

!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn, e)
do j = jst, jed
do i = ist, ied

  !dir$ vector aligned
  !dir$ simd
  do k = kst, ked
     d(k,i,j) = ( x(k,i-1,j  ) &
              +   x(k,i+1,j  ) &
              +   x(k,i  ,j-1) &
              +   x(k,i  ,j+1) ) * r + rhs(k,i,j)
  end do
  d(kst,i,j) = d(kst,i,j) + rhs(kst-1,i,j)*r
  d(ked,i,j) = d(ked,i,j) + rhs(ked+1,i,j)*r

  ! TDMA
  w(kst,i,j) = c(kst)

  !dir$ vector aligned
  do k=kst+1, ked
    e = 1.0 / (1.0 - a(k) * w(k-1,i,j))
    w(k,i,j) = e * c(k)
    d(k,i,j) = (d(k,i,j) - a(k) * d(k-1,i,j)) * e
  end do

  !dir$ vector aligned
  do k=ked-1, kst, -1
    d(k,i,j) = d(k,i,j) - w(k,i,j) * d(k+1,i,j)
  end do


  !dir$ vector aligned
  !dir$ simd
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
real                                                   ::  r
real, dimension(1-g:sz(3)+g)                           ::  a, w, e
!dir$ assume_aligned d:64, x:64, w:64, rhs:64, a:64, e:64

ist = idx(0)
ied = idx(1)
jst = idx(2)
jed = idx(3)
kst = idx(4)
ked = idx(5)

r = 1.0/6.0

flop = flop + dble((ied-ist+1)*(jed-jst+1)*(ked-kst+1))*(5.0+5.0+5.0) &
            + dble((ied-ist+1)*(jed-jst+1)) * 5.0;


!$OMP PARALLEL DO SCHEDULE(static) &
!$OMP REDUCTION(+:res) PRIVATE(pp, dp, pn)
do j = jst, jed
do i = ist, ied

  !dir$ vector aligned
  !dir$ simd
  do k = kst, ked
     d(k,i,j) = ( x(k,i-1,j  ) &
              +   x(k,i+1,j  ) &
              +   x(k,i  ,j-1) &
              +   x(k,i  ,j+1) ) * r + rhs(k,i,j)
  end do
  d(kst,i,j) = d(kst,i,j) + rhs(kst-1,i,j)*r
  d(ked,i,j) = d(ked,i,j) + rhs(ked+1,i,j)*r

  ! TDMA
  d(kst,i,j) = d(kst,i,j)

  !dir$ vector aligned
  do k=kst+1, ked
    d(k,i,j) = (d(k,i,j) - a(k) * d(k-1,i,j)) * e(k)
  end do

  !dir$ vector aligned
  do k=ked-1, kst, -1
    d(k,i,j) = d(k,i,j) - w(k) * d(k+1,i,j)
  end do


  !dir$ vector aligned
  !dir$ simd
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
