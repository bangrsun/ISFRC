!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          m2coef computes the response functions due to           **
!**          the magnetic probes.                                    **
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
      subroutine m2coef(rr, nr, zz, nz, coef,  mp, nc)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nw,nh,nwnh
      use fcoil
      use coilsp
      use mprobe
      use pmodel
      use consta
      use fshift
      use bfgrid
      implicit integer*4 (i-n), real*8 (a-h, o-z)
			real*8,dimension(nr) :: rr
			real*8,dimension(nz) :: zz
      dimension coef(mp,nc)
!
      if (.not. allocated(brgridfc)) then
        allocate(brgridfc(nwnh,nfcoil))
        brgridfc(:,:) = 0.0
      endif
      if (.not. allocated(bzgridfc)) then
        allocate(bzgridfc(nwnh,nfcoil))
        bzgridfc(:,:) = 0.0
      endif
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      do m=1,magpr2
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
         if (nz > 0) then
            do ii=1,nr
               do jj=1,nz
                  kk=(ii-1)*nz+jj
                  a=rr(ii)
                  brct=0.0
                  bzct=0.0
                  do mmm=1,nsmp2
                     r=xmp20+(mmm-1)*delsx
                     z=ymp20+(mmm-1)*delsy-zz(jj)
                     brtmp=br(a,r,z)*tmu
                     bztmp=bz(a,r,z)*tmu
                     brct=brtmp+brct
                     bzct=bztmp+bzct
                  enddo 
                  cmp2=(brct*cosm+bzct*sinm)/nsmp2
                  coef(m,kk)=cmp2
               enddo
            enddo
         else
            do k=1,nfcoil
!---------------------------------------------------------------
!-- Shifted F-coil                                            --
!---------------------------------------------------------------
               if (k==nshiftrz(k)) then
                  pmnow=radeg*pmprobe(m)
                  pfnow=radeg*pshift(k)
               endif
!
               brct=0
               bzct=0
               call splitc(isplit,rsplt,zsplt,csplt, &
                           rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
               do l=1,itot
                  a=rsplt(l)
                  do mmm=1,nsmp2
                     r1=xmp20+(mmm-1)*delsx
                     z1=ymp20+(mmm-1)*delsy-zsplt(l)
!---------------------------------------------------------------
!-- Shifted F-coil                                            --
!---------------------------------------------------------------
                     if (k==nshiftrz(k)) then
                        rcos=r1*cos(pmnow)-rshift(k)*cos(pfnow)
                        rsin=r1*sin(pmnow)-rshift(k)*sin(pfnow)
                        r1=sqrt(rcos**2+rsin**2)
                        z1=z1-zshift(k)
                     endif
!
                     brc=br(a,r1,z1)*tmu
                     bzc=bz(a,r1,z1)*tmu
                     bzct=bzct+bzc/fitot
!---------------------------------------------------------------
!-- Shifted F-coil ?                                          --
!---------------------------------------------------------------
                     if (k==nshiftrz(k)) then
                        cospp=rcos/r1
                        sinpp=rsin/r1
                        cfactor=cos(pmnow)*cospp+sin(pmnow)*sinpp
                        brct=brct+brc/fitot*cfactor
                     else
                        brct=brct+brc/fitot
                     endif
                  enddo
               enddo
!
               coef(m,k)=(brct*cosm+bzct*sinm)/nsmp2
            enddo
         endif
      enddo
!----------------------------------------------------------------
!-- BR, BZ at grid due to F-coils                              --
!----------------------------------------------------------------
      if (nz > 0) then
         return
      else
         do ii=1,nw
            do jj=1,nh
               kk=(ii-1)*nh+jj
               do k=1,nfcoil
                  brct=0
                  bzct=0
                  call splitc(isplit,rsplt,zsplt,csplt, &
                              rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
                  do l=1,itot
                     a=rsplt(l)
                     r1=rgrid(ii)
                     z1=zgrid(jj)-zsplt(l)
                     brc=br(a,r1,z1)*tmu
                     bzc=bz(a,r1,z1)*tmu
                     bzct=bzct+bzc/fitot
                     brct=brct+brc/fitot
                  enddo 
                  brgridfc(kk,k)=brct
                  bzgridfc(kk,k)=bzct
               enddo
            enddo
         enddo
      endif
!
      return
      end subroutine m2coef

