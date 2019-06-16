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
  use global_params,only:ngr1,ngz1,dr,dz,rgrid_rz &
    ,Jzeta_rz,pprim_rz,psi_rz,Jztot
  implicit none
  integer*4 :: i,j
  real*8 :: xtmp,sechx

  Jzeta_rz=0.0d0
  do j=1,ngz1
    do i=2,ngr1
      if(ipres==1) then
        pprim_rz(i,j)=0.0
        if(psi_rz(i,j)>0.0) then
          pprim_rz(i,j)=cn(1)
        endif
      elseif(ipres==2) then ! assume P=c1*sech(c2*(psi-c3))
        xtmp=cn(2)*(psi_rz(i,j)-cn(3))
        sechx=2.0/(dexp(xtmp)+dexp(-xtmp))
        pprim_rz(i,j)=-cn(2)*cn(1)*dtanh(xtmp)*sechx/mu0
      else
        xtmp=cn(2)*(psi_rz(i,j)-cn(3))
        pprim_rz(i,j)=cn(1)*(1.0-dtanh(xtmp))/mu0
      endif
    enddo
  enddo
  Jzeta_rz(:,:)=rgrid_rz(:,:)*pprim_rz(:,:)
  ! solid boundary condition in r direction
  Jzeta_rz(1,:)=0.0d0
  Jzeta_rz(ngr1,:)=0.0d0
  ! solid boundary condition in z direction
  !Jzeta_rz(:,1)=0.0d0
  !Jzeta_rz(:,ngz1)=0.0d0
  Jzeta_rz(:,ngz1)=Jzeta_rz(:,1)
  Jztot=0.0d0
  do j=2,ngz
    do i=2,ngr
      Jztot=Jztot+Jzeta_rz(i,j)*dr*dz
    enddo
  enddo
  
  return
end subroutine calcJzeta


!***********************************************************************
! subprogram description:                                              *
!      calcpsi calculate the psi according to the plasma current       *
!      Jzeta and the Green's function by plasma itself.                *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcpsipnew
  use consta,only:mu0
  use readin_params,only:ngr,ngz
  use global_params,only:ngr1,ngz1,dr,dz,rgrid_rz &
    ,pprim_rz,psip_rz,psipnew_rz
  implicit none
  integer*4 :: i,j
  real*8 :: w=0.8
  real*8 :: r,cdiag,cr1,cr2,cz

  do j=2,ngz
    do i=2,ngr
      r=rgrid_rz(i,j)
      cr1=r/(dr*dr*(r+0.5d0*dr))
      cr2=r/(dr*dr*(r-0.5d0*dr))
      cz=1.0d0/(dz*dz)
      cdiag=1.0d0/(cr1+cr2+2.0d0*cz)
      psipnew_rz(i,j)=(1-w)*psip_rz(i,j) &
        +w*cdiag*( &
                   cr1*psip_rz(i+1,j)+cr2*psip_rz(i-1,j) &
                  +cz*(psip_rz(i,j+1)+psip_rz(i,j-1)) &
                  +mu0*r*r*pprim_rz(i,j) &
                 )
    enddo
  enddo
  psipnew_rz(1,:)=0.0d0
  psipnew_rz(ngr1,:)=0.0d0
  !psipnew_rz(:,1)=0.0d0
  !psipnew_rz(:,ngz1)=0.0d0

end subroutine calcpsipnew

