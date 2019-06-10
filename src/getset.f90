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
!      getset performs inputing and initialization.                    *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine getset
  use exparam
  use nio
  use domain
  use vesel
  use fcoil
  use acoil
  use ecoil
  !use fldiag
  !use mpdiag
  !use rogdiag
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
  read(fu_input,nml=expset)
  if(nrgrid <= 0 .or. nzgrid <= 0) then
    write(*,*) "Wrong grid settings for nrgrid and nzgrid !!!"
    stop
  endif
  !-- if vesel enabled ------------------------------------------------
  if(nvesel > 0) then
    allocate(r_v(nvesel),z_v(nvesel),w_v(nvesel),h_v(nvesel), &
      ar_v(nvesel),az_v(nvesel))
    allocate(gfvesel_rzv(nrgrid,nzgrid,nvesel))
    r_v=0.0d0
    z_v=0.0d0
    w_v=0.0d0
    h_v=0.0d0
    ar_v=0.0d0
    az_v=0.0d0
    read(fu_input,nml=veselset)
    if(w_v(nvesel) == 0.0d0 .or. h_v(nvesel) == 0.0d0) then
      write(*,*) "Waring: length of vesel segments coordinate < nfcoil !"
    endif
  endif
  !-- if fcoil enabled ------------------------------------------------
  if(nfcoil > 0) then
    allocate(r_f(nfcoil),z_f(nfcoil),w_f(nfcoil),h_f(nfcoil), &
      ar_f(nfcoil),az_f(nfcoil),nsr_f(nfcoil),nsz_f(nfcoil))
    allocate(gffcoil_rzf(nrgrid,nzgrid,nfcoil))
    r_f=0.0d0
    z_f=0.0d0
    w_f=0.0d0
    h_f=0.0d0
    ar_f=0.0d0
    az_f=0.0d0
    nsr_f=1
    nsz_f=1
    read(fu_input,nml=fcoilset)
    if(w_f(nfcoil) == 0.0d0 .or. h_f(nfcoil) == 0.0d0) then
      write(*,*) "Waring: length of fcoil coordinate < nfcoil !"
    endif
  endif
  !-- if acoil enabled ------------------------------------------------
  if(nacoil > 0) then

  endif
  !-- if ecoil enabled ------------------------------------------------
  if(necoil > 0) then

  endif
  !-- if fldiag enabled -----------------------------------------------
  if(nfldiag > 0) then

  endif
  !-- if mpdiag enabled -----------------------------------------------
  if(nmpdiag > 0) then

  endif
  !-- if rogdiag enabled ----------------------------------------------
  if(nrogdiag > 0) then

  endif

  close(fu_input)
!----------------------------------------------------------------------
!-- read grid parameters or make grid meshes                         --
!----------------------------------------------------------------------
  allocate(rgrid_rz(nrgrid,nzgrid),zgrid_rz(nrgrid,nzgrid))
  allocate(gfplas_rzrz(nrgrid,nzgrid,nrgrid,nzgrid))
  rgrid_rz=0.0d0
  zgrid_rz=0.0d0
  if(igrid == 0) then
    dr=(rmax-rmin)/float(nrgrid-1)
    dz=(zmax-zmin)/float(nzgrid-1)
    do j=1,nzgrid
      do i=1,nrgrid
        rgrid_rz(i,j)=rmin+(i-1)*dr
        zgrid_rz(i,j)=zmin+(j-1)*dz
      enddo
    enddo
    open(unit=fu_grid,status='unknown',file='grid.dat')
    write(fu_grid,*) ((rgrid_rz(i,j), i=1,nrgrid), j=1,nzgrid)
    write(fu_grid,*) ((zgrid_rz(i,j), i=1,nrgrid), j=1,nzgrid)
    close(fu_grid)
  elseif(igrid == 1) then
    inquire(file='grid.dat',exist=file_exist)
    if(.not. file_exist) then
      write(*,*) 'Cannot find file grid.dat !!!'
      stop
    endif
    open(unit=fu_grid,status='old',file='grid.dat')
    read(fu_grid,*) ((rgrid_rz(i,j), i=1,nrgrid), j=1,nzgrid)
    read(fu_grid,*) ((zgrid_rz(i,j), i=1,nrgrid), j=1,nzgrid)
    close(fu_grid)
  endif

  return
end subroutine getset

