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
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine load
  use nio,only:fu_equ,fu_snap
  use readin_params,only:iequ,ngr,ngz,rmax,rmin,zmax,zmin,tol
  use global_params,only:dr,dz,ngr1,ngz1,rgrid_rz,zgrid_rz,psif_rz &
    ,psi_rz,psinew_rz,dpsi_rz,psip_rz,psipnew_rz,dpsip_rz,Jzeta_rz,Jztot
  implicit none
  real*8 :: dpsi,dpsip
  integer*4 :: i,j
  logical file_exist

  if(iequ==0) then
!----------------------------------------------------------------------
!-- for calculated equilibirum                                       --
!----------------------------------------------------------------------
    !-- make grid meshes ----------------------------------------------
    dr=(rmax-rmin)/float(ngr)
    dz=(zmax-zmin)/float(ngz)
    do j=1,ngz1
      do i=1,ngr1
        rgrid_rz(i,j)=rmin+(i-1)*dr
        zgrid_rz(i,j)=zmin+(j-1)*dz
      enddo
    enddo
    !-- calculate psi by fcoils ---------------------------------------
    call calcpsif
    psi_rz(:,:)=psif_rz(:,:)

    dpsi=10.0d0*tol
    !write(*,*) "Izeta = ",Izeta
    do while(dpsi>tol)
      !-- psi_rz -> pprim_rz & Jzeta_rz through fix -------------------
      call calcJzeta
      psip_rz(:,:)=psi_rz(:,:)
      !-- pprim_rz -> psip_rz through SOR iteration -------------------
      dpsip=10*tol
      psipnew_rz=0.0d0
      do while(dpsip>tol)
        call calcpsipnew
        dpsip_rz(:,:)=psipnew_rz(:,:)-psip_rz(:,:)
        dpsip=maxval(abs(dpsip_rz))
        psip_rz(:,:)=psipnew_rz(:,:)
      enddo
      !-- psip_rz -> psinew_rz by adding psif_rz ----------------------
      psinew_rz(:,:)=psip_rz(:,:)+psif_rz(:,:)
      dpsi_rz(:,:)=psinew_rz(:,:)-psi_rz(:,:)
      dpsi=maxval(abs(dpsi_rz))
      psi_rz(:,:)=psinew_rz(:,:)
      write(*,*) "dpsi = ",dpsi,", Jztot = ",Jztot
      !write(*,*) "Izeta-Jztot = ",Izeta-Jztot
    enddo
  else
!----------------------------------------------------------------------
!-- for input equilibirum                                            --
!----------------------------------------------------------------------
    inquire(file='equilibrium.dat',exist=file_exist)
    if(.not. file_exist) then
      write(*,*) 'Cannot file file equilibirum.dat !!!'
      stop
    endif
    open(unit=fu_equ,status='old',file='equilibiurm.dat')
    read(fu_equ,*) ((rgrid_rz(i,j),i=1,ngr1),j=1,ngz1)
    read(fu_equ,*) ((zgrid_rz(i,j),i=1,ngr1),j=1,ngz1)
    read(fu_equ,*) ((psi_rz(i,j),i=1,ngr1),j=1,ngz1)
    dr=rgrid_rz(2,1)-rgrid_rz(1,1)
    dz=zgrid_rz(1,2)-zgrid_rz(1,1)
    !-- calculate psi by fcoils ---------------------------------------
    call calcpsif
  endif ! iequ==0
!----------------------------------------------------------------------
!-- output psi and Jzeta                                             --
!----------------------------------------------------------------------
  open(unit=fu_snap,status='unknown',file='snap.dat')
  write(fu_snap,*) rgrid_rz
  write(fu_snap,*) zgrid_rz
  write(fu_snap,*) psi_rz
  write(fu_snap,*) psif_rz
  write(fu_snap,*) psip_rz
  write(fu_snap,*) Jzeta_rz
  flush(fu_snap)
  close(fu_snap)

  return
end subroutine load


!***********************************************************************
! subprogram description:                                              *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine calcpsif
  use consta,only:tmu
  use readin_params,only:nfcoil,r_f,z_f,w_f,h_f,ar_f,az_f,nsr_f,nsz_f,J_f
  use global_params,only:ngr1,ngz1,rgrid_rz,zgrid_rz,psif_rz
  implicit none
  real*8,external :: mutpsi
  real*8 :: r1,z1,jsplit,gsplit
  real*8,dimension(:,:),allocatable :: rsplit,zsplit
  integer i,j,jf,isr,isz
  psif_rz=0.0d0
  if(nfcoil > 0) then
    do j=1,ngz1
      do i=1,ngr1
        r1=rgrid_rz(i,j)
        z1=zgrid_rz(i,j)
        do jf=1,nfcoil
          allocate(rsplit(nsr_f(jf),nsz_f(jf)),zsplit(nsr_f(jf),nsz_f(jf)))
          call splitcoil(r_f(jf),z_f(jf),w_f(jf),h_f(jf),ar_f(jf),az_f(jf), &
            nsr_f(jf),nsz_f(jf),rsplit,zsplit)
          Jsplit=J_f(jf)/nsr_f(jf)/nsz_f(jf)
          do isz=1,nsz_f(jf)
            do isr=1,nsr_f(jf)
              gsplit=mutpsi(rsplit(isr,isz),r1,z1-zsplit(isr,isz))
              psif_rz(i,j)=psif_rz(i,j)+tmu*Jsplit*gsplit
            enddo
          enddo
          deallocate(rsplit,zsplit)
        enddo
      enddo
    enddo
  endif

  return
end subroutine calcpsif


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

