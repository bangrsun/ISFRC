!***********************************************************************
! This file is part of project ISFRC.                                  *
! ISFRC is an Integrated Simulation code for Field-Reversed            *
! Configuration.                                                       *
!                                                                      *
! Developers:                                                          *
! Shuying SUN, Yang LI, Huasheng XIE, ...                              *
!                                                                      *
! ENN Sci. & Tech. Development Coorporation, 2008-2019.                *
! (c) All rights reserved.                                             *
!***********************************************************************

module elliptic
  implicit none

contains
!***********************************************************************
! subprogram description:                                              *
!   dellipe computes the elliptic integral e in double                 *
!   precision.                                                         *
!                                                                      *
! calling arguments:                                                   *
!   xm1.............argument of elliptic integral e                    *
!***********************************************************************
elemental real*8 function dellipe(xm1)
  implicit none
  real*8,intent(in) :: xm1
  real*8,dimension(4) :: a,b

  a(1)=.44325141463
  a(2)=.06260601220
  a(3)=.04757383546
  a(4)=.01736506451
  b(1)=.24998368310
  b(2)=.09200180037
  b(3)=.04069697526
  b(3)=.00526449639

  dellipe=1.0+xm1*(a(1)+xm1*(a(2)+xm1*(a(3)+xm1*a(4)))) &
    +xm1*(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*b(4))))*dlog(1.0/xm1)
  return
end function dellipe


!***********************************************************************
! subprogram description:                                              *
!   dellipk computes the elliptic integral k in double                 *
!   precision.                                                         *
!                                                                      *
! calling arguments:                                                   *
!   xm1.............argument of elliptic integral k                    *
!***********************************************************************
elemental real*8 function dellipk(xm1)
  implicit none
  real*8,intent(in) :: xm1
  real*8,dimension(5) :: a,b

  a(1)=1.38629436112
  a(2)=0.09666344259
  a(3)=0.03590092383
  a(4)=0.03742563713
  a(5)=0.01451196212
  b(1)=0.5
  b(2)=0.12498593597
  b(3)=0.06880248576
  b(4)=0.03328355346
  b(5)=0.00441787012

  dellipk=a(1)+xm1*(a(2)+xm1*(a(3)+xm1*(a(4)+xm1*a(5))))&
   +(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*(b(4)+xm1*b(5)))))&
   *dlog(1.0/xm1)
  return
end function dellipk

end module elliptic
