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
!      |                                     __\_                      
!      |                               __---`  / ``-                   
!      |                   hcoil __---`       /                        
!      |                   __---`            /                         
!      |             \_---`                 /wcoil                     
!      |      |a1coil/                     /                           
! rcoil| - - -|- - -/- - - - - -          /                            
!      |      |    /           |       __/                             
!      |      |   /            | __---`   ``-                          
!      |      |  /         __---`                                      
!      |      | /    __---`    | a2coil                                
!      |      |/_---`__________|___                                    
!      |                       |                                       
!------|-------------------------------------------------------------->
!     O|                       zcoil                       Z coordinate


!***********************************************************************
! subprogram description:                                              *
!      mutpsi computes mutual inductance/2/pi between two              *
!      circular filaments of radii a1 and r1 and                       *
!      separation of z1, for mks units multiply returned               *
!      value by 2.0e-07.                                               *
!                                                                      *
! calling arguments:                                                   *
!   a1..............first filament radius                              *
!   r1..............second filament radius                             *
!   z1..............vertical separation                                *
!***********************************************************************
real*8 function mutpsi(a1,r1,z1)
  use consta,only:tole
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,tmp1,tmp2,xk,x1,cay,ee
  real*8,external :: dellipk,dellipe

  a=a1
  r=r1
  z=z1
  tmp1=(a+r)*(a+r)+z*z
  tmp2=(a-r)*(a-r)+z*z
  xk=4.0d0*a*r/tmp1
  x1=tmp2/tmp1
  if(x1 < 1.0d-10) x1=1.0d-10
  cay=dellipk(x1)
  ee=dellipe(x1)

  mutpsi=dsqrt(tmp1)*((1.0d0-0.5d0*xk)*cay-ee)
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
!   a1..............first filament radius                              *
!   r1..............second filament radius                             *
!   z1..............vertical separation                                *
!***********************************************************************
real*8 function mutbr(a1,r1,z1)
  use consta,only:tole
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,tmp1,tmp2,xk,cay,ee
  real*8,external :: dellipk,dellipe

  a=a1
  r=r1
  z=z1
  tmp1=(a+r)*(a+r)+z*z
  tmp2=(a-r)*(a-r)+z*z
  if(tmp2 < tole) tmp2=tole
  xk=4.0d0*a*r/tmp1
  cay=dellipk(xk)
  ee=dellipe(xk)

  mutbr=z/(r*dsqrt(tmp1))*((a*a+r*r+z*z)/tmp2*ee-cay)
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
!   a1..............first filament radius                              *
!   r1..............second filament radius                             *
!   z1..............vertical separation                                *
!***********************************************************************
real*8 function mutbz(a1,r1,z1)
  use consta,only:tole
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,tmp1,tmp2,xk,cay,ee
  real*8,external :: dellipk,dellipe

  a=a1
  r=r1
  z=z1
  tmp1=(a+r)*(a+r)+z*z
  tmp2=(a-r)*(a-r)+z*z
  if(tmp2 < tole) tmp2=tole
  xk=4.0d0*a*r/tmp1
  cay=dellipk(xk)
  ee=dellipe(xk)

  mutbz=((a*a-r*r-z*z)/tmp2*ee+cay)/dsqrt(tmp1)
  return
end function mutbz


!***********************************************************************
! subprogram description:                                              *
!   dellipe computes the elliptic integral e in double                 *
!   precision.                                                         *
!                                                                      *
! calling arguments:                                                   *
!   xm1.............argument of elliptic integral e                    *
!***********************************************************************
real*8 function dellipe(xm1)
  implicit none
  real*8,intent(in) :: xm1
  real*8,dimension(4) :: a,b

  data a(1),a(2),a(3),a(4)/.44325141463,.06260601220,&
    .04757383546,.01736506451/
  data b(1),b(2),b(3),b(4)/.24998368310,.09200180037,&
    .04069697526,.00526449639/

  dellipe=1.0+xm1*(a(1)+xm1*(a(2)+xm1*(a(3)+xm1*a(4))))&
   +xm1*(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*b(4))))*dlog(1.0/xm1)
  return
end function dellipe


!***********************************************************************
! subprogram description:                                              *
!   dellipk computes the elliptic integral k in double                 *
!   precision.                                                         *
!                                                                      *
! calling arguments:                                                   *
!   xm1.............argument of elliptic integral k                    *
!***********************************************************************
real*8 function dellipk(xm1)
  implicit none
  real*8,intent(in) :: xm1
  real*8,dimension(5) :: a,b

  data a(1),a(2),a(3),a(4),a(5)/1.38629436112,.09666344259,&
    .03590092383,.03742563713,.01451196212/
  data b(1),b(2),b(3),b(4),b(5)/.5,.12498593597,.06880248576,&
    .03328355346,.00441787012/

  dellipk=a(1)+xm1*(a(2)+xm1*(a(3)+xm1*(a(4)+xm1*a(5))))&
   +(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*(b(4)+xm1*b(5)))))&
   *dlog(1.0/xm1)
  return
end function dellipk

