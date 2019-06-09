!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          gsilop computes the Green's function at si loops due    **
!**          to filament currents flowing in r(n) and z(n).          **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**       rr..............r coordinates                              **
!**       zz..............z coordinates                              **
!**       nr..............dimension of r                             **
!**       rspfun..........computed green's functions values          **
!**       nz..............dimension of z                             **
!**                                                                  **
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
      subroutine gsilop(rr, nr, zz, nz, rspfun, ns, rsi, zsi, wsi &
           , hsi, as, as2, ndim)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nw,nh,nwnh
      use consta
      implicit integer*4 (i-n), real*8 (a-h, o-z)
			real*8,dimension(nr) :: rr
			real*8,dimension(nz) :: zz
			real*8,dimension(ns) :: rsi,zsi,wsi,hsi,as,as2
			real*8,dimension(ndim,nwnh) :: rspfun
!     dimension rsi(1),zsi(1),wsi(1),hsi(1),as(1),as2(1)
!     dimension rr(1),zz(1),rspfun(ndim,1)
      dimension taf(nfcoil),taf2(nfcoil)
      data isplit/17/,tole/1.0e-10/

      ndr = isplit
      ndz = isplit
      in=ns
        taf(in) = tan(as(in)*pi/180.0)
        taf2(in) = tan(as2(in)*pi/180.0)
        do 1500 nc=1,nr
        do 1500 nd=1,nz
        nk=(nc-1)*nz+nd
        rsum = 0.
        dwc = wsi(in)/ndr
        dhc = hsi(in)/ndz
        if (as2(in) ~= 0.) go to 640
        z1 = zsi(in) - taf(in)*(wsi(in)-dwc)/2. - .5*hsi(in) + .5*dhc
        ab = rsi(in) - .5*wsi(in) + .5*dwc
        do 620 iw = 1, isplit
          drc = ab + (iw-1)*dwc + iw*tole
          z2 = z1 + (iw-1)*taf(in)*dwc
          do 610 ih = 1, isplit
            dzc = z2 + (ih-1)*dhc
            rtmp = psical(drc,rr(nc),zz(nd)-dzc)
            rsum = rsum + rtmp
  610     continue
  620   continue
        go to 660
  640   continue
        do 655 ih = 1, ndz
          dzc = zsi(in) - .5*hsi(in) + .5*dhc + dhc*(ih-1)
          do 645 iw = 1, ndr
            drc = rsi(in) - .5*wsi(in) - .5*hsi(in)/taf2(in)&
                   + .5*dwc + .5*dhc/taf2(in)&
                   + dhc/taf2(in)*(ih-1) + dwc*(iw-1)
            rtmp = psical(drc,rr(nc),zz(nd)-dzc)
            rsum = rsum + rtmp
  645     continue
  655   continue
  660   continue
        cmut = rsum*2.e-07/(isplit*isplit)
        rspfun(in,nk)=cmut
 1500   continue

      return
      end subroutine gsilop


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          m1coef computes the response functions due to           **
!**          the thin flux loops.                                    **
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
      subroutine m1coef(rr, zz, nr, nz, coef,  nl, nf)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      use fcoil
      use coilsp
      use consta
      use nio
      use siloop
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(nr) :: rr
      real*8,dimension(nz) :: zz
      real*8,dimension(nsilop,nr*nz) :: coef
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      if (nf <= 0) then
         do ii=1,nr
            do jj=1,nz
               kk=(ii-1)*nz+jj
               a=rr(ii)
               r=rsi(m)
               z=zsi(m)-zz(jj)
               cmp2=psical(a,r,z)*tmu
               coef(m,kk)=cmp2
            enddo 
         enddo 
      else
         k=nf
         psict=0
         call splitc(isplit,rsplt,zsplt,csplt, &
                     rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
         do l=1,itot
            a=rsplt(l)
            r1=rsi(m)
            z1=zsi(m)-zsplt(l)
            psic=psical(a,r1,z1)*tmu
            psict=psict+psic/fitot
         enddo 
         coef(m,k)=psict
      endif
!
      return
      end subroutine m1coef

