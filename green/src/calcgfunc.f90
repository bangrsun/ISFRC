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
!      calc_gf calculate the Green's function of all active coils,     *
!      including the fcoil, plasma, vessel, acoil, ecoil.              *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcgfunc
  use exparam,only:nrgrid,nzgrid,nfcoil,nvesel
  use domain
  use fcoil
  use vesel
  implicit none
  integer i,j
  real*8,dimension(nrgrid,nzgrid) :: gftmp

  if(nfcoil > 0) then
    do i=1,nfcoil
      call greenfunc(r_f(i),z_f(i),w_f(i),h_f(i),ar_f(i),az_f(i), &
        nsr_f(i),nsz_f(i),nrgrid,nzgrid,rgrid_rz,zgrid_rz,gftmp)
      gffcoil_rzf(:,:,i)=gftmp
    enddo
  endif
  if(nvesel > 0) then
    do i=1,nvesel
      call greenfunc(r_v(i),z_v(i),w_v(i),h_v(i),ar_v(i),az_v(i), &
        1,1,nrgrid,nzgrid,rgrid_rz,zgrid_rz,gftmp)
      gfvesel_rzv(:,:,i)=gftmp
    enddo
  endif
 
  return
end subroutine calcgfunc


!***********************************************************************
! subprogram description:                                              *
!      greenfunc computes the Green's functions at (r,z)               *
!      due to current loops.                                           *
!***********************************************************************
subroutine greenfunc(rc,zc,wc,hc,acr,acz,nsr,nsz, &
    ngr,ngz,rgrid_rz,zgrid_rz,gfgrid_rz)
  implicit none
  !> R & Z coordinate of current loop
  real*8,intent(in) :: rc,zc
  !> width & height of current loop
  real*8,intent(in) :: wc,hc
  !> angle with R & Z direction of current loop
  real*8,intent(in) :: acr,acz
  !> number of splits in R & Z direction of current loop
  integer*4,intent(in) :: nsr,nsz
  !> number of grid points in R & Z direction
  integer*4,intent(in) :: ngr,ngz
  !> coordinate of grid points in R & Z direction
  real*8,dimension(ngr,ngz),intent(in ) :: rgrid_rz,zgrid_rz
  !> output Green's function at grid points
  real*8,dimension(ngr,ngz),intent(out) :: gfgrid_rz
  !> coordinate of split coil in R & Z direction
  real*8,dimension(:,:),allocatable :: rsplit,zsplit
  !> Green's function of each split coil at grid points
  real*8,dimension(ngr,ngz) :: gsplit_rz
  real*8,external :: mutpsi
  integer*4 :: ic,isz,isr,igz,igr

  gfgrid_rz=0.0d0
  gsplit_rz=0.0d0
  allocate(rsplit(nsr,nsz),zsplit(nsr,nsz))
  call splitcoil(rc,zc,wc,hc,acr,acz,nsr,nsz,rsplit,zsplit)
  do isz=1,nsz
    do isr=1,nsr
!----------------------------------------------------------------------
!--  compute the Green's functions at (r,z) due to split coil        --
!----------------------------------------------------------------------
      do igz=1,ngz
        do igr=1,ngr
          gsplit_rz(igr,igz)=mutpsi(rsplit(isr,isz), &
            rgrid_rz(igr,igz), zgrid_rz(igr,igz)-zsplit(isr,isz))
        enddo
      enddo
      gfgrid_rz=gfgrid_rz+gsplit_rz
    enddo
  enddo
  deallocate(rsplit,zsplit)

  return
end subroutine greenfunc

