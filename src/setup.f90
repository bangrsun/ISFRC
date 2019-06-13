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
  integer i,j
  real*8 dr,dz

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
  allocate(rgrid_rz(ngr,ngz),zgrid_rz(ngr,ngz))
  allocate(dr_rz(ngr-1,ngz-1),dz_rz(ngr-1,ngz-1))
  allocate(gffcoil_rzf(ngr,ngz,nfcoil),gfplas_rzrz(ngr,ngz,ngr,ngz))
  allocate(psi_rz(ngr,ngz),psip_rz(ngr,ngz),psif_rz(ngr,ngz))
  allocate(psinew_rz(ngr,ngz),dpsi_rz(ngr,ngz))
  allocate(pres_rz(ngr,ngz),dp_rz(ngr-1,ngz-1),Jzeta_rz(ngr,ngz))
  rgrid_rz=0.0d0
  zgrid_rz=0.0d0
  dr_rz=0.0d0
  dz_rz=0.0d0
  gfplas_rzrz=0.0d0
  gffcoil_rzf=0.0d0
  psi_rz=0.0d0
  psip_rz=0.0d0
  psif_rz=0.0d0
  psinew_rz=0.0d0
  dpsi_rz=0.0d0
  pres_rz=0.0d0
  dp_rz=0.0d0
  Jzeta_rz=0.0d0
!----------------------------------------------------------------------
!-- read grid parameters or make grid meshes                         --
!----------------------------------------------------------------------
  if(igrid == 0) then
    dr=(rmax-rmin)/float(ngr-1)
    dz=(zmax-zmin)/float(ngz-1)
    dr_rz=dr
    dz_rz=dz
    do j=1,ngz
      do i=1,ngr
        rgrid_rz(i,j)=rmin+(i-1)*dr
        zgrid_rz(i,j)=zmin+(j-1)*dz
      enddo
    enddo
    open(unit=fu_grid,status='unknown',file='grid.dat')
    write(fu_grid,*) ((rgrid_rz(i,j), i=1,ngr), j=1,ngz)
    write(fu_grid,*) ((zgrid_rz(i,j), i=1,ngr), j=1,ngz)
    close(fu_grid)
  elseif(igrid == 1) then
    inquire(file='grid.dat',exist=file_exist)
    if(.not. file_exist) then
      write(*,*) 'Cannot find file grid.dat !!!'
      stop
    endif
    open(unit=fu_grid,status='old',file='grid.dat')
    read(fu_grid,*) ((rgrid_rz(i,j), i=1,ngr), j=1,ngz)
    read(fu_grid,*) ((zgrid_rz(i,j), i=1,ngr), j=1,ngz)
    do j=1,ngz-1
      do i=1,ngr-1
        dr_rz(i,j)=rgrid_rz(i,j+1)-rgrid_rz(i,j)
        dz_rz(i,j)=zgrid_rz(i,j+1)-zgrid_rz(i,j)
      enddo
    enddo
    close(fu_grid)
  endif

  return
end subroutine setup

