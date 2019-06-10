!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          psical computes mutual inductance/2/pi between two      **
!**          circular filaments of radii a1 and r1 and               **
!**          separation of z1, for mks units multiply returned       **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     calling arguments:                                           **
!**       a1..............first filament radius                      **
!**       r1..............second filament radius                     **
!**       z1..............vertical separation                        **
!**                                                                  **
!**     references:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
real*8 function psical(a1,r1,z1)
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,den,xk,x1,cay,ee
  real*8,external :: xmdelk,xmdele

  a=a1
  r=r1
  z=z1
  den=a*a+r*r+z*z+2.0d0*a*r
  xk=4.*a*r/den
  x1=(a*a+r*r+z*z-2.*a*r)/den
  if(x1 < 1.0d-10) x1=1.0d-10
  cay=xmdelk(x1)
  ee=xmdele(x1)
!----------------------------------------------------------------------
!--   psi computation                                                --
!----------------------------------------------------------------------
  psical= sqrt(den)*((1.0d+00-0.5d+00*xk)*cay-ee)
  return
end function psical


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          br     computes mutual inductance/2/pi between two      **
!**          circular filaments of radii a1 and r1 and               **
!**          separation of z1, for mks units multiply returned       **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     calling arguments:                                           **
!**       a1..............first filament radius                      **
!**       r1..............second filament radius                     **
!**       z1..............vertical separation                        **
!**                                                                  **
!**     references:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
real*8 function br(a1,r1,z1)
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,den,xk,x1,cay,ee
  real*8,external :: xmdelk,xmdele

  a=a1
  r=r1
  z=z1
  den=a*a+r*r+z*z+2.0d0*a*r
  xk=4.0d0*a*r/den
  x1=(a*a+r*r+z*z-2.0d0*a*r)/den
  if(x1 < 1.0d-10) x1=1.0d-10
  cay=xmdelk(x1)
  ee=xmdele(x1)
!----------------------------------------------------------------------
!--   br  computation                                                --
!----------------------------------------------------------------------
  br=z/(r* sqrt(den))*(-cay+(a*a+r*r+z*z)/((a-r)*(a-r)+z*z)*ee)
  return
end function br


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          bz     computes mutual inductance/2/pi between two      **
!**          circular filaments of radii a1 and r1 and               **
!**          separation of z1, for mks units multiply returned       **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     calling arguments:                                           **
!**       a1..............first filament radius                      **
!**       r1..............second filament radius                     **
!**       z1..............vertical separation                        **
!**                                                                  **
!**     references:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
real*8 function bz(a1,r1,z1)
  implicit none
  real*8,intent(in) :: a1,r1,z1
  real*8 :: a,r,z,den,xk,x1,cay,ee
  real*8,external :: xmdelk,xmdele

  a=a1
  r=r1
  z=z1
  den=a*a+r*r+z*z+2.0d0*a*r
  xk=4.0d0*a*r/den
  x1=(a*a+r*r+z*z-2.0d0*a*r)/den
  if(x1 < 1.0d-10) x1=1.0d-10
  cay=xmdelk(x1)
  ee=xmdele(x1)
!----------------------------------------------------------------------
!--   bz  computation                                                --
!----------------------------------------------------------------------
  bz=(cay+(a*a-r*r-z*z)/((a-r)*(a-r)+z*z)*ee)/ sqrt(den)
  return
end function bz


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine splitc(isplit,rs,zs,cs,rcoil,zcoil,wcoil,hcoil,a1coil,a2coil,cc)
  use consta,only:pi
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  integer*4,intent(in) :: isplit
  real*8,dimension(isplit),intent(out) :: rs,zs,cs
  real*8,intent(in) :: rcoil,zcoil,wcoil,hcoil,a1coil,a2coil,cc

  frd=pi/180.0d0
  if(a1coil+a2coil==0.0d0) then
!----------------------------------------------------------------------
!-- rectangle                                                        --
!----------------------------------------------------------------------
    wdelt=wcoil/isplit
    hdelt=hcoil/isplit
    rstrt=rcoil-wcoil/2.0d0+wdelt/2.0d0
    zstrt=zcoil-hcoil/2.0d0+hdelt/2.0d0
    zz=zstrt
    ic=0
    c=cc/(isplit*isplit)
    do ii=1,isplit
      rr=rstrt
      do jj=1,isplit
        ic=ic+1
        zs(ic)=zz
        rs(ic)=rr
        cs(ic)=c
        rr=rr+wdelt
      enddo 
      zz=zz+hdelt
    enddo
  elseif(a1coil ~= 0.0d0) then
!----------------------------------------------------------------------
!-- a1coil ~= 0                                                      --
!----------------------------------------------------------------------
    side=tan(frd*a1coil)*wcoil
    hdelt=hcoil/isplit
    wdelt=wcoil/isplit
    zdelt=tan(frd*a1coil)*wdelt
    rstrt=rcoil-wcoil/2.+wdelt/2.
    tsid=hcoil+side
    zstrt =zcoil-tsid/2.+tsid/2.*1./isplit
    rr=rstrt
    ic=0
    c=cc/(isplit*isplit)
    do ii=1,isplit
      zz=zstrt+(ii-1)*zdelt
      do jj=1,isplit
        ic=ic+1
        zs(ic)=zz
        rs(ic)=rr
        cs(ic)=c
        zz=zz+hdelt
      enddo 
      rr=rr+wdelt
    enddo
  elseif(a2coil ~= 0.0d0) then
!----------------------------------------------------------------------
!-- a2coil ~= 0                                                      --
!----------------------------------------------------------------------
    side=hcoil/tan(frd*a2coil)
    hdelt=hcoil/isplit
    wdelt=wcoil/isplit
    zstrt=zcoil-hcoil/2.+hdelt/2.
    rdelt=hdelt/tan(frd*a2coil)
    rstrt=rcoil-side/2.-wcoil/2.+rdelt/2.+wdelt/2.
    side=hcoil/tan(frd*a2coil)
    wtot=side+wcoil
    whaf=(side+wcoil)/2.
    rcorn=rcoil-whaf
    rcornr=rcoil+whaf
    rcorn2=rcorn+wtot/isplit
    rstrt=(rcorn+rcorn2)/2.
    zz=zstrt
    ic=0
    c=cc/(isplit*isplit)
    do ii=1,isplit
      rr=rstrt+(ii-1)*rdelt
      do jj=1,isplit
        ic=ic+1
        zs(ic)=zz
        rs(ic)=rr
        cs(ic)=c
        rr=rr+wdelt
      enddo 
      zz=zz+hdelt
    enddo
  endif

  return
end subroutine splitc

