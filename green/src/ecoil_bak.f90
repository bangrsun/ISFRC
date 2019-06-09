!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          gecoil gets the Green's functions for Ohmic heating     **
!**          coils.                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          30/01/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine gecoil(rsilec,rmp2ec,gridec,rgrid,mw,zgrid,mh, &
                  rfcec,recec,rsisec)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
  use consta
  use fcoil
  use cecoil
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension rsilec(nsilop,nesum),rmp2ec(magpr2,nesum), &
            rfcec(nfcoil,nesum),recec(nesum,nesum), &
            rsisec(nesum)
	real*8,dimension(mw) :: rgrid
	real*8,dimension(mh) :: zgrid
	real*8,dimension(nwnh,nesum) :: gridec
  dimension taf(nfcoil),taf2(nfcoil)
	real*8 :: zetaec = 3.5e-08
  do j=1,nsilop
    do i=1,nesum
      rsilec(j,i)=0.0
    enddo 
    jj=j
    do i=1,necoil
      ii=i
      call e1coef(work,jj,ii)
      kkm=ecid(i)
      rsilec(j,kkm)=rsilec(j,kkm)+work*ecturn(i)
    enddo 
  enddo

  do j=1,magpr2
    do i=1,nesum
      rmp2ec(j,i)=0.0
    enddo 
    jj=j
    do i=1,necoil
      ii=i
      call e2coef(work,jj,ii)
      kkm=ecid(i)
      rmp2ec(j,kkm)=rmp2ec(j,kkm)+work*ecturn(i)
    enddo 
  enddo

  do i=1,nw
    nr=i
    do j=1,nh
      nz=j
      kk=(i-1)*nh+j
      do m=1,nesum
        gridec(kk,m)=0.0
      enddo 
      do n=1,necoil
        nn=n
        call egrid(work,rgrid,nr,zgrid,nz,nn)
        kkkm=ecid(n)
        gridec(kk,kkkm)=gridec(kk,kkkm)+work*ecturn(n)
      enddo 
    enddo
  enddo

  aaa=0.0
  bbb=0.0
  do j=1,nfcoil
    do i=1,nesum
      rfcec(j,i)=0.0
    enddo 
    taf(j)=tan(af(j)*pi/180.)
    taf2(j)=tan(af2(j)*pi/180.)
    do i=1,necoil
      call flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
        rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
      work=work*0.5/pi
      kkm=ecid(i)
      kk=kkm+1
      rfcec(j,kkm)=rfcec(j,kkm)+work*(kk-ecid(i))
    enddo 
  enddo

  do j=1,nesum
    rsisec(j)=0.0
    do i=1,nesum
      recec(j,i)=0.0
    enddo 
  enddo 
  do j=1,necoil
    jjm=ecid(j)
    do i=1,necoil
      call flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
        re(j),ze(j),we(j),he(j),aaa,bbb,work)
      work=work*0.5/pi
      kkm=ecid(i)
      recec(jjm,kkm)=recec(jjm,kkm)+work
    enddo 
    rsisec(jjm)=rsisec(jjm)+2.*pi*re(j)/we(j)/he(j)*zetaec
  enddo
  return
end subroutine gecoil


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          egrid computes the response functions due to            **
!**          the grid points and Ohmic heating coils.                **
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
subroutine egrid(coef, rgrid, nr, zgrid, nz, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil
  use cecoil
  use coilsp
  use consta
  use nio
  use mprobe
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  real*8,dimension(nr) :: rgrid
  real*8,dimension(nz) :: zgrid


  radeg=pi/180.
  isplit=17
  itot=isplit*isplit
  fitot=itot

  k=ne
  psict=0
  aaa=0.0
  bbb=0.0
  call splitc(isplit,rsplt,zsplt,csplt, &
              re(k),ze(k),we(k),he(k),aaa,bbb,cdum)
  do l=1,itot
    a=rsplt(l)
    r1=rgrid(nr)
    z1=zgrid(nz)-zsplt(l)
    psic=psical(a,r1,z1)*tmu
    psict=psict+psic/fitot
  enddo 
  coef=psict


  return
end subroutine egrid


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          e1coef computes the response functions due to           **
!**          the thin flux loops and Ohmic heating coils.            **
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
subroutine e1coef(coef,  nl, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil
  use cecoil
  use coilsp
  use consta
  use nio
  use siloop
  implicit integer*4 (i-n), real*8 (a-h, o-z)

  radeg=pi/180.
  isplit=17
  itot=isplit*isplit
  fitot=itot

  m=nl
  k=ne
  psict=0
  aaa=0.0
  bbb=0.0
  call splitc(isplit,rsplt,zsplt,csplt, &
              re(k),ze(k),we(k),he(k),aaa,bbb,cdum)
  do l=1,itot
     a=rsplt(l)
     r1=rsi(m)
     z1=zsi(m)-zsplt(l)
     psic=psical(a,r1,z1)*tmu
     psict=psict+psic/fitot
  enddo 
  coef=psict

  return
end subroutine e1coef


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          e2coef computes the response functions due to           **
!**          the magnetic probes and Ohmic heating coils.            **
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
subroutine e2coef(coef, mp, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil
  use cecoil
  use coilsp
  use consta
  use mprobe
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
  aaa=0.0
  bbb=0.0
  call splitc(isplit,rsplt,zsplt,csplt, &
              re(k),ze(k),we(k),he(k),aaa,bbb,cdum)
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
end subroutine e2coef

