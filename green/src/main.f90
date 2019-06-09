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

!***********************************************************************
! subprogram description:                                              *
!      gfun generates the necessary response                           *
!      functions used by ISFRC for reconstruction of the               *
!      magnetic surfaces and current density profile.                  *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
program gfun
  implicit none
  real*8 :: start, finish

  call cpu_time(start)

  call getset
  call calcgfunc
  call writegfunc
  !call calcrfunc
  !call writerfunc

  call cpu_time(finish)
  write(*,*) "CPU time used = ",finish - start,'s'

  stop 'GREEN TABLE GENERATED!'
end program gfun

