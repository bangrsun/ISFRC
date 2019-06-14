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
program main
  use consta,only:mu0,pi
  use readin_params,only:nitermax,tol
  use global_params
  use nio
  implicit none
  !real*8 :: start, finish
  integer*4 :: it
  real*8 :: dpsi

  !call cpu_time(start)
!----------------------------------------------------------------------
!-- read settings and allocate memories                              --
!----------------------------------------------------------------------
  write(*,*) 'SETUP STARTED!'
  call setup
!----------------------------------------------------------------------
!-- calculate the Green's function for plasma and fcoil              --
!----------------------------------------------------------------------
  write(*,*) 'CALCULATE GREEN TABLE STARTED!'
  call calcgfunc
  call writegfunc
!----------------------------------------------------------------------
!-- calculate psi by fcoils                                          --
!----------------------------------------------------------------------
  call calcpsif
!----------------------------------------------------------------------
!-- main iteration loops                                             --
!----------------------------------------------------------------------
  write(*,*) 'MAIN ITERATION LOOP STARTED!'
  it=0
  dpsi=10*tol
  psip_rz=0.0d0
  psi_rz(:,:)=psip_rz(:,:)+psif_rz(:,:)
  do while(it<nitermax .and. dpsi>tol)
    it=it+1
    write(*,*) "it = ",it
    !-- calculate Jzeta from psi---------------------------------------
    call calcJzeta
    !-- update psip from Jzeta ----------------------------------------
    call calcpsip
    psinew_rz(:,:)=psip_rz(:,:)+psif_rz(:,:)
    !-- calculate delta psi after update ------------------------------
    dpsi_rz(:,:)=psinew_rz(:,:)-psi_rz(:,:)
    dpsi=maxval(abs(dpsi_rz))

    psi_rz(:,:)=psinew_rz(:,:)
  enddo
!----------------------------------------------------------------------
!-- output psi and Jzeta                                             --
!----------------------------------------------------------------------
  open(unit=fu_snap,status='unknown',file='snap.dat')
  write(fu_snap,*) rgrid_rz
  write(fu_snap,*) zgrid_rz
  write(fu_snap,*) psi_rz
  write(fu_snap,*) Jzeta_rz
  flush(fu_snap)
  close(fu_snap)

  stop
end program main


