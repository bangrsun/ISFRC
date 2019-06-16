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
  use readin_params,only:Izeta,nitermax,tol
  use global_params
  use nio
  implicit none
  !real*8 :: start, finish
  integer*4 :: it
  real*8 :: dpsip,dpsi

  open(unit=fu_output,status='unknown',file='output.dat')
!----------------------------------------------------------------------
!-- read settings and allocate memories                              --
!----------------------------------------------------------------------
  write(fu_output,*) "=================================================="
  write(fu_output,*) "SETUP STARTED!"
  call setup
!----------------------------------------------------------------------
!-- load initial psi by calculate psif from fcoil                    --
!----------------------------------------------------------------------
  write(fu_output,*) "=================================================="
  write(fu_output,*) 'LOAD STARTED!'
  call load
!----------------------------------------------------------------------
!-- main iteration loops                                             --
!----------------------------------------------------------------------
  write(fu_output,*) "=================================================="
  write(fu_output,*) 'MAIN ITERATION LOOP STARTED!'
  it=0
  dpsi=10.0d0*tol
  !write(*,*) "Izeta = ",Izeta
  do while(it<nitermax .and. dpsi>tol)
    it=it+1
    !-- psi_rz -> pprim_rz & Jzeta_rz through fix ---------------------
    call calcJzeta
    !-- pprim_rz -> psip_rz through SOR iteration ---------------------
    dpsip=10*tol
    psip_rz(:,:)=psi_rz(:,:)
    do while(dpsip>tol)
      call calcpsipnew
      dpsip_rz(:,:)=psipnew_rz(:,:)-psip_rz(:,:)
      dpsip=maxval(abs(dpsip_rz))
      psip_rz(:,:)=psipnew_rz(:,:)
    enddo
    !-- psip_rz -> psinew_rz by adding psif_rz ------------------------
    psinew_rz(:,:)=psip_rz(:,:)+psif_rz(:,:)
    dpsi_rz(:,:)=psinew_rz(:,:)-psi_rz(:,:)
    dpsi=maxval(abs(dpsi_rz))
    psi_rz(:,:)=psinew_rz(:,:)
    write(*,*) "it = ",it,", dpsi = ",dpsi,", Jztot = ",Jztot
    !write(*,*) "Izeta-Jztot = ",Izeta-Jztot
  enddo
!----------------------------------------------------------------------
!-- output psi and Jzeta                                             --
!----------------------------------------------------------------------
  open(unit=fu_snap,status='unknown',file='snap.dat')
  write(fu_snap,*) rgrid_rz
  write(fu_snap,*) zgrid_rz
  write(fu_snap,*) psi_rz
  write(fu_snap,*) psif_rz
  write(fu_snap,*) psip_rz
  write(fu_snap,*) pprim_rz
  flush(fu_snap)
  close(fu_snap)

  close(fu_output)
  stop
end program main


