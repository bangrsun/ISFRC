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
  real*8,parameter :: mu0 = pi*4.0d-07
  !> tmu=mu0/2.0d0/pi
  real*8,parameter :: tmu = 2.0d-07
  real*8,parameter :: small = 1.0d-10
end module consta


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
  !> snapshot are outputed in this file
  integer*4,parameter :: fu_snap = 202

end module nio


!***********************************************************************
! module for read in parameters                                        *
!***********************************************************************
module readin_params
  implicit none
  public

  !> number of grids in R & Z direction, without outer boundary
  integer*4 :: ngr,ngz
  !> border of calculation box, 
  real*8 :: rmin,rmax,zmin,zmax
  !> number of p.f. coils, <=0 disabled, >0 enabled
  integer*4 :: nfcoil
!----------------------------------------------------------------------
!-- fcoil parameters                                                 --
!----------------------------------------------------------------------
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
  !> current at each ploidal field coil
  real*8,dimension(:),allocatable :: J_f


  integer*4 :: ipres
  real*8 :: Izeta
  real*8,dimension(3) :: cn
  !> maximum iteration steps
  integer*4 :: nitermax
  !> output tolerance for delta psi
  real*8 :: tol

  namelist/simuset/ ngr,ngz,rmin,rmax,zmin,zmax &
    ,nfcoil &
    ,ipres,Izeta,cn,nitermax,tol

  namelist/fcoilset/ r_f,z_f,w_f,h_f,ar_f,az_f,nsr_f,nsz_f,J_f

end module readin_params


!***********************************************************************
! module for global_params                                             *
!***********************************************************************
module global_params
  implicit none
  public

  !> number of grids in R & Z direction, with outer boundary,
  ! ngr1=ngr+1, ngz1=ngz+1
  integer*4 :: ngr1,ngz1
  !> delta r & delta z
  real*8 :: dr,dz
  !> R coordinate of grid points
  real*8,dimension(:,:),allocatable :: rgrid_rz
  !> Z coordinate of grid points
  real*8,dimension(:,:),allocatable :: zgrid_rz
  !> psi at grid points
  real*8,dimension(:,:),allocatable :: psi_rz
  !> psi induced by fcoil at grid points
  real*8,dimension(:,:),allocatable :: psif_rz
  !> new psi after every iteration
  real*8,dimension(:,:),allocatable :: psinew_rz
  !> delta psi between two successive iteration
  real*8,dimension(:,:),allocatable :: dpsi_rz
  !> pressure at grid points
  real*8,dimension(:,:),allocatable :: pres_rz
  !> dp/dpsi at grid points
  real*8,dimension(:,:),allocatable :: pprim_rz
  !> current density at grid points
  real*8,dimension(:,:),allocatable :: Jzeta_rz
  real*8 :: Jztot

end module global_params


