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
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
program main
  use readin_params,only:Izeta,nt
  use global_params
  use nio
  implicit none
  !real*8 :: start, finish
  integer*4 :: it

  open(unit=fu_output,status='unknown',file='output.dat')
!----------------------------------------------------------------------
!-- read settings and allocate memories                              --
!----------------------------------------------------------------------
  write(fu_output,*) "=================================================="
  write(fu_output,*) "SETUP STARTED!"
  call setup
!----------------------------------------------------------------------
!-- load initial equilibrium                                         --
!----------------------------------------------------------------------
  write(fu_output,*) "=================================================="
  write(fu_output,*) 'LOAD STARTED!'
  call load
!----------------------------------------------------------------------
!-- main iteration loops                                             --
!----------------------------------------------------------------------
  write(fu_output,*) "=================================================="
  write(fu_output,*) 'MAIN ITERATION LOOP STARTED!'
  do it=1,nt
  enddo

  close(fu_output)
  stop
end program main


