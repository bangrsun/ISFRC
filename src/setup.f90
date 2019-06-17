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
!      setup performs inputing and initialization.                     *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine setup
  use nio
  use readin_params
  use global_params
  implicit none
  logical file_exist

!----------------------------------------------------------------------
!-- read experimental settings                                       --
!----------------------------------------------------------------------
  inquire(file='input.dat',exist=file_exist)
  if(.not. file_exist) then
    write(*,*) 'Cannot find file input.dat !!!'
    stop
  endif
  open(unit=fu_input,status='old',file='input.dat')
  read(fu_input,nml=simuset)
  if(ngr <= 0 .or. ngz <= 0) then
    write(*,*) "Wrong grid settings for ngr and ngz !!!"
    stop
  endif
!----------------------------------------------------------------------
!-- read fcoil settings                                              --
!----------------------------------------------------------------------
  if(nfcoil <= 0) then
    write(*,*) "Wrong fcoil settings : nfcoil = ",nfcoil
    stop
  endif
  allocate(r_f(nfcoil),z_f(nfcoil),w_f(nfcoil),h_f(nfcoil), &
    ar_f(nfcoil),az_f(nfcoil),nsr_f(nfcoil),nsz_f(nfcoil),J_f(nfcoil))
  r_f=0.0d0
  z_f=0.0d0
  w_f=0.0d0
  h_f=0.0d0
  ar_f=0.0d0
  az_f=0.0d0
  nsr_f=1
  nsz_f=1
  J_f=0.0d0
  read(fu_input,nml=fcoilset)
  if(w_f(nfcoil) == 0.0d0 .or. h_f(nfcoil) == 0.0d0) then
    write(*,*) "Warning: length of fcoil coordinate < nfcoil !"
  endif
  close(fu_input)
!----------------------------------------------------------------------
!-- allocate memory for global parameters                            --
!----------------------------------------------------------------------
  ngr1=ngr+1
  ngz1=ngz+1
  allocate(rgrid_rz(ngr1,ngz1),zgrid_rz(ngr1,ngz1))
  allocate(psif_rz(ngr1,ngz1))
  allocate(psi_rz(ngr1,ngz1),psinew_rz(ngr1,ngz1),dpsi_rz(ngr1,ngz1))
  allocate(psip_rz(ngr1,ngz1),psipnew_rz(ngr1,ngz1),dpsip_rz(ngr1,ngz1))
  allocate(pprim_rz(ngr1,ngz1),Jzeta_rz(ngr1,ngz1))
  rgrid_rz=0.0d0
  zgrid_rz=0.0d0
  psif_rz=0.0d0
  psi_rz=0.0d0
  psinew_rz=0.0d0
  dpsi_rz=0.0d0
  psip_rz=0.0d0
  psipnew_rz=0.0d0
  dpsip_rz=0.0d0
  pprim_rz=0.0d0
  Jzeta_rz=0.0d0

  return
end subroutine setup

