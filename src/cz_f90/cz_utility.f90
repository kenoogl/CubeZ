!###################################################################################
!#
!# CubeZ
!#
!# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
!# All rights reserved.
!#
!###################################################################################

!> @file   cz_utility.f90
!! @brief  Utilitiy functions
!! @author aics
!<

subroutine fileout (sz, g, s, dh, org, fname)
implicit none
integer                                                :: nn, ix, jx, kx, i, j, k, g
integer, dimension(3)                                  :: sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: s
real                                                   :: dh, rtime
real, dimension(3)                                     :: org
character*20                                           :: fname

ix = sz(1)
jx = sz(2)
kx = sz(3)

nn = 0
rtime = 0.0

open (unit=22, file=fname, form='unformatted')
write (22) 1, 1
write (22) ix, jx, kx
write (22) org(1), org(2), org(3)
write (22) dh, dh, dh
write (22) nn, rtime
write (22) (((s(i,j,k),i=1,ix),j=1,jx),k=1,kx)
close (unit=22)

return
end subroutine fileout


!> *******************************
subroutine exact (sz, g, e, dh, org)
implicit none
include 'cz_fparam.fi'
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz
real, dimension(3)                                     ::  org
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  e
real                                                   ::  dh, pi, r2, x, y, z

ix = sz(1)
jx = sz(2)
kx = sz(3)

r2 = sqrt(2.0)
pi = 2.0*asin(1.0)

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) PRIVATE(x, y, z)
do k=1, kx
do j=1, jx
do i=1, ix
  x = org(1) + dh*real(i-1)
  y = org(2) + dh*real(j-1)
  z = org(3) + dh*real(k-1)
  e(i,j,k) = sin(pi*x)*sin(pi*y) / sinh(r2*pi) * ( sinh(r2*pi*z)-sinh(r2*pi*(z-1.0)) )
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine exact


!> *******************************
subroutine err (sz, idx, g, d, p, e, loc)
implicit none
integer                                                ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                  ::  sz, loc
integer, dimension(0:5)                                ::  idx
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, e
real                                                   ::  r
double precision                                       ::  d, q

ix = sz(1)
jx = sz(2)
kx = sz(3)

d=0.0
loc(1)=-1
loc(2)=-1
loc(3)=-1

!$OMP PARALLEL DO SCHEDULE(static) COLLAPSE(2) &
!$OMP PRIVATE(r, q)
do k = 1,kx
do j = 1,jx
do i = 1,ix
  r = p(i,j,k) - e(i,j,k)
  e(i,j,k) = r
  q = abs(dble(r))
  if (d < q) then
    d = q
    loc(1) = i
    loc(2) = j
    loc(3) = k
  endif
end do
end do
end do
!$OMP END PARALLEL DO

return
end subroutine err
