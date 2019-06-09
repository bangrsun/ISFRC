!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          gvesel gets the Green's functions for the vessel        **
!**          segments.                                               **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine gvesel(rsilvs,rmp2vs,gridvs,rgrid,mw, &
                        zgrid,mh,rfcvs,rvsfc,rvsec)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      use consta
      use fcoil
      use cecoil
      use cvesel
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rsilvs(nsilop,nvesel),rmp2vs(magpr2,nvesel), &
                rgrid(1),zgrid(1),rvsec(nvesel,nesum), &
                rfcvs(nfcoil,nvesel),rvsfc(nvesel,nfcoil)
      dimension gridvs(1,nvesel)
      dimension taf(nfcoil),taf2(nfcoil)
      dimension tas(nvesel),tas2(nvesel)
      do j=1,nsilop
         jj=j
         do i=1,nvesel
            ii=i
            call v1coef(work,jj,ii)
            rsilvs(j,i)=work
         enddo 
      enddo 
!
      do j=1,magpr2
         jj=j
         do i=1,nvesel
            ii=i
            call v2coef(work,jj,ii)
            rmp2vs(j,i)=work
         enddo 
      enddo 
!
      do i=1,nw
         nr=i
         do j=1,nh
            nz=j
            kk=(i-1)*nh+j
            do n=1,nvesel
               nn=n
               call vgrid(work,rgrid,nr,zgrid,nz,nn)
               gridvs(kk,n)=work
            enddo 
         enddo 
		  enddo
!
      do j=1,nfcoil
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         do i=1,nvesel
            tas(i)=tan(avs(i)*pi/180.)
            tas2(i)=tan(avs2(i)*pi/180.)
            call flux(rvs(i),zvs(i),wvs(i),hvs(i),tas(i),tas2(i), &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            rfcvs(j,i)=work
         enddo 
      enddo 
!
      aaa=0.0
      do j=1,nvesel
         tas(j)=tan(avs(j)*pi/180.)
         tas2(j)=tan(avs2(j)*pi/180.)
         do i=1,nesum
            rvsec(j,i)=0.0
         enddo 
         do i=1,necoil
            call flux(re(i),ze(i),we(i),he(i),aaa,aaa, &
                      rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j),work)
            work=work*0.5/pi
            kkm=ecid(i)
            kk=kkm+1
            rvsec(j,kkm)=rvsec(j,kkm)+work*(kk-ecid(i))
         enddo
      enddo 
!
      do j=1,nvesel
         tas(j)=tan(avs(j)*pi/180.)
         tas2(j)=tan(avs2(j)*pi/180.)
         do i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
            call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                      rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j),work)
            work=work*0.5/pi
            rvsfc(j,i)=work
         enddo 
      enddo 
      return
      end subroutine gvesel


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          vgrid computes the response functions due to            **
!**          the grid points and the vessel segments.                **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine vgrid(coef, rgrid, nr, zgrid, nz, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                  nesum,nfsum,nvsum,nvesel,nacoil
  use coilsp
  use siloop
  use consta
  use nio
  use cvesel
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  real*8,dimension(nr) :: rgrid
  real*8,dimension(nz) :: zgrid
  data init/0/

  radeg=pi/180.
  isplit=17
  itot=isplit*isplit
  fitot=itot

  k=ne
  psict=0
  call splitc(isplit,rsplt,zsplt,csplt, &
              rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
  do l=1,itot
    a=rsplt(l)
    r1=rgrid(nr)
    z1=zgrid(nz)-zsplt(l)
    psic=psical(a,r1,z1)*tmu
    psict=psict+psic/fitot
  enddo 
  coef=psict

  return
end subroutine vgrid


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          v1coef computes the response functions due to           **
!**          the thin flux loops and the vessel segments.            **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          15/07/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine v1coef(coef,  nl, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                  nesum,nfsum,nvsum,nvesel,nacoil
  use coilsp
  use consta
  use nio
  use siloop
  use cvesel
  implicit integer*4 (i-n), real*8 (a-h, o-z)

  radeg=pi/180.
  isplit=17
  itot=isplit*isplit
  fitot=itot

  m=nl
  k=ne
  psict=0
  call splitc(isplit,rsplt,zsplt,csplt, &
              rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
  do l=1,itot
     a=rsplt(l)
     r1=rsi(m)
     z1=zsi(m)-zsplt(l)
     psic=psical(a,r1,z1)*tmu
     psict=psict+psic/fitot
  enddo 
  coef=psict

  return
end subroutine v1coef


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          v2coef computes the response functions due to           **
!**          the magnetic probes and the vessel segments.            **
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
subroutine v2coef(coef, mp, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                  nesum,nfsum,nvsum,nvesel,nacoil
  use coilsp
  use mprobe
  use consta
  use cvesel
  implicit integer*4 (i-n), real*8 (a-h, o-z)

  radeg=pi/180.
  isplit=17
  itot=isplit*isplit
  fitot=itot

  m=mp
  if (smp2(m) > 0.0) then
    sinm=sin(radeg*amp2(m))
    cosm=cos(radeg*amp2(m))
    delsx=smp2(m)/nsmp2*cosm
    delsy=smp2(m)/nsmp2*sinm
  else
!------------------------------------------------------------------------------
!--  perpendicular probes    96/02/04                                        --
!------------------------------------------------------------------------------
    sinm=sin(radeg*amp2(m))
    cosm=cos(radeg*amp2(m))
    sinms=sin(radeg*(amp2(m)+90.))
    cosms=cos(radeg*(amp2(m)+90.))
    delsx=abs(smp2(m))/nsmp2*cosms
    delsy=abs(smp2(m))/nsmp2*sinms
  endif
  xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
  ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
  k=ne
  brct=0
  bzct=0
  call splitc(isplit,rsplt,zsplt,csplt, &
              rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
  do l=1,itot
     a=rsplt(l)
     do mmm=1,nsmp2
        r1=xmp20+(mmm-1)*delsx
        z1=ymp20+(mmm-1)*delsy-zsplt(l)
        brc=br(a,r1,z1)*tmu
        bzc=bz(a,r1,z1)*tmu
        brct=brct+brc/fitot
        bzct=bzct+bzc/fitot
     enddo 
  enddo 
  coef=(brct*cosm+bzct*sinm)/nsmp2

  return
end subroutine v2coef

