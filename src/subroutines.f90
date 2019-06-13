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
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine splitcoil(rcoil,zcoil,wcoil,hcoil,a1coil,a2coil, &
    nsr,nsz,rsplit_rz,zsplit_rz)
  use consta,only:pi
  implicit none
  real*8,   intent(in) :: rcoil,zcoil,wcoil,hcoil,a1coil,a2coil
  integer*4,intent(in) :: nsr,nsz
  real*8,dimension(nsr,nsz),intent(out) :: rsplit_rz,zsplit_rz
  real*8 :: frd,dw,dh,dr1,dz1,dr2,dz2,rstrt,zstrt
  integer :: ii,jj

  frd=pi/180.0d0
  dw=wcoil/nsr
  dh=hcoil/nsz
  dr1=dw*dcos(frd*a1coil)
  dz1=dh*dcos(frd*a2coil)
  dr2=dh*dsin(frd*a2coil)
  dz2=dw*dsin(frd*a1coil)
  rstrt=rcoil-0.5d0*(wcoil*dcos(frd*a1coil)+hcoil*dsin(frd*a2coil)) &
    +0.5d0*(dr1+dr2)
  zstrt=zcoil-0.5d0*(hcoil*dcos(frd*a2coil)+wcoil*dsin(frd*a1coil)) &
    +0.5d0*(dz1+dz2)
  do jj=1,nsz
    do ii=1,nsr
      rsplit_rz(ii,jj)=rstrt+(ii-1)*dr1+(jj-1)*dr2
      zsplit_rz(ii,jj)=zstrt+(jj-1)*dz1+(ii-1)*dz2
    enddo 
  enddo

  return
end subroutine splitcoil
!      ^ R coordinate                                                  
!      |              /<                                               
!      |      |a1coil/`---_                                            
!      |      |     /      ``---_hcoil                                 
!      |      |    /             ``---_                                
!      |      |   /                    ``---_ >/                       
!      |      |  /                           `/`--_                    
! rcoil| _ _ _|_/ _ _ _ _ _ _ _              /                         
!      |      |/______________|____         /                          
!      |        `---_         |            /                           
!      |             ``---_   | a2coil    / wcoil                      
!      |                   ``---_        /                             
!      |                      |  ``---_ /                              
!      |                      |        ``--_                           
!      |                      |                                        
!------|-------------------------------------------------------------->
!     O|                     zcoil                         Z coordinate

!***********************************************************************
! subprogram description:                                              *
!      mutpsi computes mutual inductance/2/pi between two              *
!      circular filaments of radii a1 and r1 and                       *
!      separation of z1, for mks units multiply returned               *
!      value by 2.0e-07.                                               *
!                                                                      *
! calling arguments:                                                   *
!   a1..............filament radius                                    *
!   r1..............radius of calcuated points                         *
!   z1..............vertical separation                                *
!***********************************************************************
real*8 function mutpsi(a1,r1,z1)
  use consta,only:small
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,tmp1,ksq,kk,ee

  a=a1
  r=r1
  z=z1
  tmp1=(a+r)*(a+r)+z*z
  !tmp2=(a-r)*(a-r)+z*z
  ksq=4.0d0*a*r/tmp1
  if(ksq<small) ksq=small
  call comelp(ksq,kk,ee)

  mutpsi=dsqrt(a*r)*((2.0d0-ksq)*kk-2.0d0*ee)/dsqrt(ksq)
  return
end function mutpsi


!***********************************************************************
! subprogram description:                                              *
!      mutbr computes mutual inductance/2/pi between two               *
!      circular filaments of radii a1 and r1 and                       *
!      separation of z1, for mks units multiply returned               *
!      value by 2.0e-07.                                               *
!                                                                      *
! calling arguments:                                                   *
!   a1..............filament radius                                    *
!   r1..............radius of calcuated points                         *
!   z1..............vertical separation                                *
!***********************************************************************
real*8 function mutbr(a1,r1,z1)
  use consta,only:small
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,tmp1,tmp2,ksq,kk,ee

  a=a1
  r=r1
  z=z1
  tmp1=(a+r)*(a+r)+z*z
  tmp2=(a-r)*(a-r)+z*z
  if(tmp2<small) tmp2=small
  ksq=4.0d0*a*r/tmp1
  call comelp(ksq,kk,ee)

  mutbr=z/(r*dsqrt(tmp1))*((a*a+r*r+z*z)/tmp2*ee-kk)
  return
end function mutbr


!***********************************************************************
! subprogram description:                                              *
!      mutbz computes mutual inductance/2/pi between two               *
!      circular filaments of radii a1 and r1 and                       *
!      separation of z1, for mks units multiply returned               *
!      value by 2.0e-07.                                               *
!                                                                      *
! calling arguments:                                                   *
!   a1..............filament radius                                    *
!   r1..............radius of calcuated points                         *
!   z1..............vertical separation                                *
!***********************************************************************
real*8 function mutbz(a1,r1,z1)
  use consta,only:small
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,tmp1,tmp2,ksq,kk,ee

  a=a1
  r=r1
  z=z1
  tmp1=(a+r)*(a+r)+z*z
  tmp2=(a-r)*(a-r)+z*z
  if(tmp2<small) tmp2=small
  ksq=4.0d0*a*r/tmp1
  call comelp(ksq,kk,ee)

  mutbz=((a*a-r*r-z*z)/tmp2*ee+kk)/dsqrt(tmp1)
  return
end function mutbz

