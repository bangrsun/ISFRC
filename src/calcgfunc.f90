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
!      including the fcoil, vessel, acoil, ecoil. This subroutine      *
!      does not calculate the Green's function from plasma itself.     *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcgfunc
  use readin_params
  use global_params
  implicit none
  real*8,dimension(ngr,ngz) :: gftmp
  real*8,external :: mutpsi
  integer i,j,jg,ig

!----------------------------------------------------------------------
!-- calculate Green's function of plasma itself                      --
!----------------------------------------------------------------------
  do jg=1,ngz
    do ig=1,ngr
      if(rgrid_rz(ig,jg)<=0.0d0) then
        gfplas_rzrz(:,:,ig,jg)=0.0d0
      else
        do j=1,ngz
          do i=1,ngr
            if(i==ig .and. j==jg) then
              gfplas_rzrz(i,j,ig,jg)=0.0d0
            else
              gfplas_rzrz(i,j,ig,jg)=mutpsi(rgrid_rz(ig,jg), &
                rgrid_rz(i,j),zgrid_rz(i,j)-zgrid_rz(ig,jg))
            endif
          enddo
        enddo
      endif
    enddo
  enddo
!----------------------------------------------------------------------
!-- calculate Green's function of fcoil                              --
!----------------------------------------------------------------------
  if(nfcoil > 0) then
    do i=1,nfcoil
      call greenfunc(r_f(i),z_f(i),w_f(i),h_f(i),ar_f(i),az_f(i), &
        nsr_f(i),nsz_f(i),ngr,ngz,rgrid_rz,zgrid_rz,gftmp)
      gffcoil_rzf(:,:,i)=gftmp
    enddo
  endif
 
  return
end subroutine calcgfunc


!***********************************************************************
! subprogram description:                                              *
!      greenfunc computes the Green's functions at (r,z)               *
!      due to coil with width and height.                              *
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
  real*8,dimension(nsr,nsz) :: rsplit,zsplit
  !> Green's function of each split coil at grid points
  real*8,dimension(ngr,ngz) :: gsplit_rz
  real*8,external :: mutpsi
  integer*4 :: isz,isr,igz,igr

  gfgrid_rz=0.0d0
  gsplit_rz=0.0d0
  call splitcoil(rc,zc,wc,hc,acr,acz,nsr,nsz,rsplit,zsplit)
  do isz=1,nsz
    do isr=1,nsr
!----------------------------------------------------------------------
!--  compute the Green's functions at (r,z) due to split coil        --
!----------------------------------------------------------------------
      do igz=1,ngz
        do igr=1,ngr
          gsplit_rz(igr,igz)=mutpsi(rsplit(isr,isz), &
            rgrid_rz(igr,igz),zgrid_rz(igr,igz)-zsplit(isr,isz))
        enddo
      enddo
      gfgrid_rz(:,:)=gfgrid_rz(:,:)+gsplit_rz(:,:)
    enddo
  enddo
  gfgrid_rz=gfgrid_rz/nsr/nsz

  return
end subroutine greenfunc

