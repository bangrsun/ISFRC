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
  use readin_params,only:ngr,ngz,ipres,cn
  use global_params,only:ngr1,ngz1,rgrid_rz,dr_rz,dz_rz, &
    Jzeta_rz,pprim_rz,psi_rz,Jztot
  implicit none
  integer*4 :: i,j
  real*8 :: xtmp,sechx

  Jzeta_rz=0.0d0
  do j=1,ngz
    do i=1,ngr
      if(ipres==1) then
        pprim_rz(i,j)=0.0
        if(psi_rz(i,j)>0.0) then
          pprim_rz(i,j)=cn(1)
        endif
      elseif(ipres==2) then ! assume P=c1*sech(c2*(psi-c3))
        xtmp=cn(2)*(psi_rz(i,j)-cn(3))
        sechx=2.0/(exp(xtmp)+exp(-xtmp))
        pprim_rz(i,j)=-cn(2)*cn(1)*tanh(xtmp)*sechx/mu0
      else
        xtmp=cn(2)*(psi_rz(i,j)-cn(3))
        pprim_rz(i,j)=cn(1)*(1.0-tanh(xtmp))/mu0
      endif
      Jzeta_rz(i,j)=rgrid_rz(i,j)*pprim_rz(i,j)
    enddo
  enddo
  ! solid boundary condition in r direction
  Jzeta_rz(1,:)=0.0d0
  Jzeta_rz(ngr1,:)=0.0d0
  ! solid boundary condition in z direction
  Jzeta_rz(:,1)=0.0d0
  Jzeta_rz(:,ngz1)=0.0d0
  Jztot=0.0d0
  do j=2,ngz
    do i=2,ngr
      Jztot=Jztot+Jzeta_rz(i,j)*dr_rz(i,j)*dz_rz(i,j)
    enddo
  enddo
  
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
  do jg=1,ngz
    do ig=1,ngr
      psip_rz(:,:)=psip_rz(:,:) &
        +mu0/2.0d0/pi*Jzeta_rz(ig,jg)*dr_rz(ig,jg)*dz_rz(ig,jg) &
        *gfplas_rzrz(:,:,ig,jg)
    enddo
  enddo

end subroutine calcpsip

