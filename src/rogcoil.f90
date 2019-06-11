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
      subroutine rogowc(rr, nrr, zz, nz, coef, nr, nc)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      use rogowl
      use coilsp
      use consta
      use fcoil
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rogpth(nrogow)
			real*8,dimension(nr) :: rr
			real*8,dimension(nz) :: zz
      dimension coef(nr,nc)
!
!
      ngrid=25
      isplit=17
      itot=isplit*isplit
      fitot=itot
      dels = 0.
      mm=1
!
      do m=1,nrogow
         if (nz > 0) then
            do inn=1,nrr
               do imm=1,nz
                  ikk=(inn-1)*nz+imm
                  coef(m,ikk)=0.0
               enddo 
            enddo 
         else
            do k=1,nfcoil
               coef(m,k)=0.0
            enddo 
         endif
         k=m
         call rogrid(ngrid,mm,k,dels)
         mm=mm+narc(m)+1
         do i=1,ngrid
            iii=i
            if(i==ngrid) then
               zl=zpg(i)-zpg(i-1)
               rl=rpg(i)-rpg(i-1)
            else
               zl=zpg(i+1)-zpg(i)
               rl=rpg(i+1)-rpg(i)
            endif
            hl=sqrt(zl*zl+rl*rl)
            sint=zl/hl
            cost=rl/hl
!
            if (nz <= 0) then
               do k=1,nfcoil
                  call splitc(isplit,rsplt,zsplt,csplt, &
                           rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
                  do l=1,itot
                     a=rsplt(l)
                     r1=rpg(i)
                     z1=zpg(i)-zsplt(l)
                     brc=br(a,r1,z1)*tmu/fitot
                     bzc=bz(a,r1,z1)*tmu/fitot
                     part=brc*cost+bzc*sint
                     call simpf(iii,fact)
                     coef(m,k)=coef(m,k)+fact*part*dels/rogpth(m)
                  enddo 
               enddo
            else
               do inn=1,nrr
                  do imm=1,nz
                     ikk=(inn-1)*nz+imm
                     a=rr(inn)
                     r1=rpg(i)
                     z1=zpg(i)-zz(imm)
                     brg=br(a,r1,z1)*tmu
                     bzg=bz(a,r1,z1)*tmu
                     part=brg*cost+bzg*sint
                     call simpf(iii,fact)
                     coef(m,ikk)=coef(m,ikk)+fact*part*dels/rogpth(m)
                  enddo 
               enddo 
            endif
         enddo
      enddo
!
      return
      end subroutine rogowc


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          rogrid calculates grid points along the given arc       **
!**          made up of at most six straight line segments.          **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**       ngrid...........                                           **
!**       mm..............                                           **
!**       m...............                                           **
!**       dels............                                           **
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
subroutine rogrid(ngrid,mm,m,dels)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
  use rogowl
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension sl(6)

  s = 0.
  mm1 = mm+narc(m)-1
  j = 1
  do i = mm,mm1
    sl(j) = sqrt((rp(i+1)-rp(i))**2 + (zp(i+1)-zp(i))**2)
    s = s+sl(j)
    j = j+1
  enddo 
  dels = s/float(ngrid-1)
  ds = 0.
  i1 = 1
  do j = 1,narc(m)
    dr = rp(mm+j)-rp(mm+j-1)
    dz = zp(mm+j)-zp(mm+j-1)
    rpg(i1) = rp(mm+j-1)+dr*ds/sl(j)
    zpg(i1) = zp(mm+j-1)+dz*ds/sl(j)
    dd = sl(j)-ds
    n1 = int(dd/dels)+i1
    i2 = i1+1
    dr = dr*dels/sl(j)
    dz = dz*dels/sl(j)
    do i = i2,n1
      rpg(i) = rpg(i-1)+dr
      zpg(i) = zpg(i-1)+dz
    enddo 
    ds = dels-(dd-float(n1-i1)*dels)
    i1 = n1+1
  enddo
  rpg(ngrid) = rp(mm+narc(m))
  zpg(ngrid) = zp(mm+narc(m))
  return
end subroutine rogrid


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**                                                                  **
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
subroutine simpf(i,f)
  !implicit integer*4 (i-n), real*8 (a-h, o-z)
  implicit none
  integer*4,intent(in) :: i
  real*8,intent(out) :: f

  if (i==1 .or. i==25) then
     f = 1./3.
  elseif ((i/2)*2==i) then
     f = 4./3.
  else
     f = 2./3.
  endif
  return
end subroutine simpf

