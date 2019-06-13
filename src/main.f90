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
!      gfun generates the necessary response                           *
!      functions used by ISFRC for reconstruction of the               *
!      magnetic surfaces and current density profile.                  *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
program main
  use consta,only:mu0,pi
  use readin_params,only:nitermax,tol
  use global_params
  use nio
  implicit none
  !real*8 :: start, finish
  integer*4 :: it
  real*8 :: dpsi

  !call cpu_time(start)
!----------------------------------------------------------------------
!-- read settings and allocate memories                              --
!----------------------------------------------------------------------
  write(*,*) 'SETUP STARTED!'
  call setup
!----------------------------------------------------------------------
!-- calculate the Green's function for plasma and fcoil              --
!----------------------------------------------------------------------
  write(*,*) 'CALCULATE GREEN TABLE STARTED!'
  call calcgfunc
  call writegfunc
!----------------------------------------------------------------------
!-- calculate psi by fcoils                                          --
!----------------------------------------------------------------------
  call calcpsif
!----------------------------------------------------------------------
!-- main iteration loops                                             --
!----------------------------------------------------------------------
  write(*,*) 'MAIN ITERATION LOOP STARTED!'
  it=0
  dpsi=10*tol
  psi_rz(:,:)=psip_rz(:,:)+psif_rz(:,:)
  do while(it<nitermax .and. dpsi>tol)
    it=it+1
    write(*,*) "it = ",it
    !-- calculate Jzeta from psi---------------------------------------
    call calcJzeta
    !-- update psip from Jzeta ----------------------------------------
    call calcpsip
    psinew_rz(:,:)=psip_rz(:,:)+psif_rz(:,:)
    !-- calculate delta psi after update ------------------------------
    dpsi_rz(:,:)=psinew_rz(:,:)-psi_rz(:,:)
    dpsi=maxval(abs(dpsi_rz))

    psi_rz(:,:)=psinew_rz(:,:)
  enddo
!----------------------------------------------------------------------
!-- output psi and Jzeta                                             --
!----------------------------------------------------------------------
  open(unit=fu_snap,status='unknown',file='snap.dat')
  write(fu_snap,*) rgrid_rz
  write(fu_snap,*) zgrid_rz
  write(fu_snap,*) psi_rz
  write(fu_snap,*) Jzeta_rz
  flush(fu_snap)
  close(fu_snap)

  stop
end program main


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
!      calcpsif calculate the psi according to the ploidal field       *
!      coil current J_f and it's Green's function.                     *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcpsif
  use consta, only:mu0,pi
  use readin_params,only:nfcoil,J_f
  use global_params,only:psif_rz,gffcoil_rzf
  implicit none
  integer*4 :: i

  psif_rz=0.0d0
  do i=1,nfcoil
    psif_rz(:,:)=psif_rz(:,:)+mu0/2.0d0/pi*J_f(i)*gffcoil_rzf(:,:,i)
  enddo
end subroutine calcpsif

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

