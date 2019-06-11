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
! module for constants                                                 *
!***********************************************************************
module consta
  implicit none
  public

  real*8,parameter :: pi = 3.1415926535897932
  real*8,parameter :: tmu = 2.0e-07
  real*8,parameter :: tole = 1.0d-10
end module consta


!***********************************************************************
! module for experimental parameters                                   *
!***********************************************************************
module exparam
  implicit none
  public

!> grids and active coil settings --------------------------------------
  !> flag for define grids, =0 uniform, =1 customed
  integer*4 :: igrid
  !> number of grids in R direction
  integer*4 :: nrgrid
  !> number of grids in Z direction
  integer*4 :: nzgrid
  !> border of calculation box, 
  !  must be defined if igrid=0, no meaning if igrid=1
  real*8 :: rmin,rmax,zmin,zmax
  !> number of vessel segments, <=0 disabled, >0 enabled
  integer*4 :: nvesel
  !> number of p.f. coils, <=0 disabled, >0 enabled
  integer*4 :: nfcoil

  namelist/expset/ igrid,nrgrid,nzgrid,rmin,rmax,zmin,zmax &
    ,nvesel,nfcoil

end module exparam


!***********************************************************************
! module for input and output file unit                                *
!***********************************************************************
module nio
  implicit none
  public

  !> expset namelist are defined in this file
  integer*4,parameter :: fu_input = 100
  !> rgrid & zgrid are defined in this file
  integer*4,parameter :: fu_grid = 101
  !> expset namelist are outputed in this file
  integer*4,parameter :: fu_output = 200
  !> green functions are outputed in this file
  integer*4,parameter :: fu_gfunc = 201
  !> response functions are outputed in this file
  integer*4,parameter :: fu_rfunc = 202

end module nio


!***********************************************************************
! module for domain                                                    *
!***********************************************************************
module domain
  implicit none
  public

  !> R coordinate of grid points, must be defined if igrid=1
  real*8,dimension(:,:),allocatable :: rgrid_rz
  !> Z coordinate of grid points, must be defined if igrid=1
  real*8,dimension(:,:),allocatable :: zgrid_rz
  !> Green's function for plasma itself
  real*8,dimension(:,:,:,:),allocatable :: gfplas_rzrz

end module domain


!***********************************************************************
! module for vessel                                                    *
!***********************************************************************
module vesel
  implicit none
  public

  !> R coordinate of vessel segments
  real*8,dimension(:),allocatable :: r_v
  !> Z coordinate of vessel segments
  real*8,dimension(:),allocatable :: z_v
  !> width of vessel segments
  real*8,dimension(:),allocatable :: w_v
  !> height of vessel segments
  real*8,dimension(:),allocatable :: h_v
  !> angle to R direction
  real*8,dimension(:),allocatable :: ar_v
  !> angle to Z direction
  real*8,dimension(:),allocatable :: az_v
  !> Green's function for p.f. coils
  real*8,dimension(:,:,:),allocatable :: gfvesel_rzv

  namelist/veselset/ r_v,z_v,w_v,h_v,ar_v,az_v

end module vesel


!***********************************************************************
! module for poloidal field coil                                       *
!***********************************************************************
module fcoil
  implicit none
  !> R coordinate of poloidal field coil
  real*8,dimension(:),allocatable :: r_f
  !> Z coordinate of poloidal field coil
  real*8,dimension(:),allocatable :: z_f
  !> width of poloidal field coil
  real*8,dimension(:),allocatable :: w_f
  !> height of poloidal field coil
  real*8,dimension(:),allocatable :: h_f
  !> angle to R direction
  real*8,dimension(:),allocatable :: ar_f
  !> angle to Z direction
  real*8,dimension(:),allocatable :: az_f
  !> splits in R direction
  integer*4,dimension(:),allocatable :: nsr_f
  !> splits in Z direction
  integer*4,dimension(:),allocatable :: nsz_f
  !> Green's function for p.f. coils
  real*8,dimension(:,:,:),allocatable :: gffcoil_rzf

  namelist/fcoilset/ r_f,z_f,w_f,h_f,ar_f,az_f,nsr_f,nsz_f

end module fcoil


