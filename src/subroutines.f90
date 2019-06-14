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
!      calcJzeta calculate Jzeta at grid points from psi according to  *
!      Grad-Shafranov equation.                                        *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcJzeta
  use consta,only:mu0
  use readin_params,only:ngr,ngz
  use global_params,only:rgrid_rz,dr_rz,dz_rz,Jzeta_rz,psi_rz
  implicit none
  integer*4 :: i,j

  do j=2,ngz-1
    do i=2,ngr-1
      Jzeta_rz(i,j)=-1.0d0/(mu0*rgrid_rz(i,j))*( &
        (psi_rz(i+1,j)-2.0*psi_rz(i,j)+psi_rz(i-1,j))/(dr_rz(i,j)*dr_rz(i,j)) &
        -(1.0d0/rgrid_rz(i,j))*(psi_rz(i+1,j)-psi_rz(i-1,j))/(2*dr_rz(i,j)) &
        +(psi_rz(i,j+1)-2.0*psi_rz(i,j)+psi_rz(i,j-1))/(dz_rz(i,j)*dz_rz(i,j)) &
        )
    enddo
  enddo
  ! solid boundary condition in r direction
  Jzeta_rz(1,:)=0.0d0
  Jzeta_rz(ngr,:)=0.0d0
  ! solid boundary condition in z direction
  Jzeta_rz(:,1)=0.0d0
  Jzeta_rz(:,ngz)=0.0d0
  ! periodic boundary condition in z direction
  !do i=2,ngr-1
  !  Jzeta_rz(i,1)=-1.0d0/(mu0*rgrid_rz(i,1))*( &
  !    (psi_rz(i+1,1)-2.0*psi_rz(i,1)+psi_rz(i-1,1))/(dr_rz(i,1)*dr_rz(i,1)) &
  !    -(1.0d0/rgrid_rz(i,1))*(psi_rz(i+1,1)-psi_rz(i-1,1))/(2*dr_rz(i,1)) &
  !    +(psi_rz(i,2)-2.0*psi_rz(i,1)+psi_rz(i,ngz))/(dz_rz(i,1)*dz_rz(i,1)) &
  !    )
  !enddo
  !Jzeta_rz(:,ngz)=Jzeta_rz(:,1)
  
  return
end subroutine calcJzeta


!***********************************************************************
! subprogram description:                                              *
!      calcpsip calculate the psi according to the plasma current      *
!      Jzeta and the Green's function by plasma itself.                *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcpsip
  use consta,only:mu0,pi
  use readin_params,only:ngr,ngz
  use global_params,only:dr_rz,dz_rz,Jzeta_rz,psip_rz,gfplas_rzrz
  implicit none
  integer*4 :: ig,jg

  psip_rz=0.0d0
  do jg=1,ngz-1
    do ig=1,ngr-1
      psip_rz(:,:)=psip_rz(:,:) &
        +mu0/2.0d0/pi*Jzeta_rz(ig,jg)*dr_rz(ig,jg)*dz_rz(ig,jg) &
        *gfplas_rzrz(:,:,ig,jg)
    enddo
  enddo

end subroutine calcpsip

