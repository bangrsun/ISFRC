!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          efun generates the necessary response                   **
!**          functions used by efit for reconstruction of the        **
!**          magnetic surfaces and current density profile.          **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1) d.w. swain and g.h. neilson, nucl. fusion           **
!**              22 (1982) 1015.                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**          28/01/85..........modified for D3D                      **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!efund module
!Revisions:
!$Log: efundud.f90,v $
!Revision 1.3.2.1  2008/11/18 22:24:13  radhakri
!*** empty log message ***
!
!Revision 1.3  2007/06/12 04:04:04  renq
!*** empty log message ***
!
!Revision 1.2  2007/06/01 05:29:02  renq
!
!efundud.f90: added by Qilong Ren
!  main EFUND program. Allocatable arrays are used to accomadate
!  different grid sizes such as 65x65,129x129,257x257 and 513x513.
!  Some memory problem with grid size 1025x1025 is to be fixed.
!
!Revision 1.3  2007/04/24 18:26:34  renq
!dummy arguments can not be allocatable arrays.
!
!Revision 1.2  2007/04/24 05:43:46  renq
!unnecessary write statements removed
!
!Revision 1.1  2007/04/24 05:16:17  renq
!Initial revision
!
program efund
  integer :: time_cc, time_cr, time_cm
  real :: elapsed_time
  ! integer :: istart, ifinish
  real :: start, finish

  ! call system_clock(istart)
  call cpu_time(start)

  call efund_getset
  call efund_matrix
  call efund_grid

  call cpu_time(finish)
  write(*,*) "CPU time used = ",finish - start,'s'
  ! call system_clock(ifinish,time_cr,time_cm)
  ! write(*,*) "CPU time used = ",(ifinish - istart)/time_cr*1.0, " s"

  stop 'GREEN TABLE GENERATED!'
end program efund


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          getset performs inputing and initialization.            **
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
subroutine efund_getset
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
  use siloop
  use cecoil
  use fcoil
  use pmodel
  use consta
  use input
  use cacoil
  use nio
  use mprobe
  use cvesel
  use fshift
!vas
  use var_filech
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension patmp2(magpr2)
  character*10 mpnam2(magpr2),lpname(nsilop),vsname(nvesel)
  NAMELIST/in3/igrid,rleft,rright,zbotto,ztop,ifcoil &
      ,islpfc,iecoil,mpnam2,xmp2,ymp2,amp2,smp2,isize,rsi,zsi,wsi &
      ,hsi,as,as2,lpname,nsmp2,ivesel,rsisvs,vsname,turnfc,patmp2 &
      ,iacoil,racoil,zacoil,wacoil,hacoil &
      ,rf,zf,fcid,wf,hf,wvs,hvs,avs,avs2,af,af2,fcturn &
      ,re,ze,ecid,ecturn,vsid,rvs,zvs,we,he &
      ,nshiftrz,rshift,zshift,pshift,pmprobe &
      ,nw,nh
!
  open(unit=nin,status='old',file='mhdin.dat')
  open(unit=nout,status='unknown',file='mhdout.dat')
!
  nw = 65
  nh = 65
!
  ifcoil=0
  igrid=0
  isize=0
  islpfc=0
  iecoil=0
  ivesel=0
  nsmp2=1
  rsi(1)=-1
  rf(1)=-1.
  re(1)=-1.
  rvs(1)=-1.
  wvs(1)=-1.
  do i=1,nfcoil
    nshiftrz(i)=0
  enddo
!---------------------------------------------------------------------
!--  isize=0      no finite size correction for flux loops          --
!--        1         finite size correction for flux loops          --
!--  islpfc=1     flux loops at F coils                             --
!---------------------------------------------------------------------
  read (nin,in3)
  write (nout,in3)
!
  if (.not. allocated(rgrid)) then
    allocate(rgrid(nw))
    rgrid(:) = 0.0
  endif
  if (.not. allocated(zgrid)) then
    allocate(zgrid(nh))
    zgrid(:) = 0.0
  endif
!
  nwnh = nw * nh
!make the file names for green-table
!vas
  call inp_file_ch(nw,nh,ch1,ch2)
!      print*,'file name : ', 'ep'//trim(ch1)// &
!                         trim(ch2)//'.ddd'
!
!----------------------------------------------------------------------
!-- read p.f. coil and flux loop dimensions                          --
!----------------------------------------------------------------------
  if (rf(1) < 0.0) then
    read (nin,*) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i),i=1,nfcoil)
  endif
  if (rsi(1) < 0.0) then
    read (nin,*) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i),i=1,nsilop)
  endif
!----------------------------------------------------------------------
!-- read ohmic heating coil and vessel segments                      --
!----------------------------------------------------------------------
  if ((iecoil > 0).or.(ivesel > 0)) then
    if (re(1) < 0.0) then
      read (nin,*) (re(i),ze(i),we(i),he(i),ecid(i),i=1,necoil)
    endif
    if ((ivesel > 0).and.(rvs(1) < 0.0)) then
      if (wvs(1) < 0.0) then
        read (nin,*) (rvs(i),zvs(i),wvs(i),hvs(i),avs(i),avs2(i), &
              i=1,nvesel)
      else
        do i=1,nvesel
          read (nin,*) rvs(i),zvs(i)
        enddo
      endif
    endif
  endif
! 200 continue
!----------------------------------------------------------------------
!--  compute r and z arrays                                          --
!----------------------------------------------------------------------
  dr=(rright-rleft)/float(nw-1)
  dz=(ztop-zbotto)/float(nh-1)
  do i=1,nw
    rgrid(i)=rleft+dr*(i-1)
  enddo 
  do i=1,nh
    zgrid(i)=zbotto+dz*(i-1)
  enddo 
! 300 continue
!
!vas
  close(nin)
  return
!10000 FORMAT (6e12.6)
!10010 FORMAT (4e12.6)
!10020 FORMAT (5e10.4)
end subroutine efund_getset


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          matrix calculates the appropriate response functions.   **
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
subroutine efund_matrix
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
  use consta
  use nio
  use cvesel
  use input
  use cacoil
  use pmodel
  use siloop
  use fcoil
  use fshift
  use bfgrid
!vas
  use var_filech
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension rfcfc(nfcoil,nfcoil)
  dimension rsilfc(nsilop,nfcoil),rmp2fc(magpr2,nfcoil), &
            rgowfc(nrogow,nfcoil)
  dimension gsilfc(nsilop,nfsum),gmp2fc(magpr2,nfsum)
  dimension rsilec(nsilop,nesum),rmp2ec(magpr2,nesum), &
            rfcec(nfcoil,nesum), &
            recec(nesum,nesum),rsisec(nesum)
  dimension rsilvs(nsilop,nvesel),rmp2vs(magpr2,nvesel), &
            rfcvs(nfcoil,nvesel), &
            rvsec(nvesel,nesum),rvsfc(nvesel,nfcoil), &
            rvsvs(nvesel,nvesel),tas(nvesel),tas2(nvesel)
  dimension gsilvs(nsilop,nvsum),gmp2vs(magpr2,nvsum)
  dimension taf(nfcoil),taf2(nfcoil)
  dimension rsilac(nsilop,nacoil),rmp2ac(magpr2,nacoil)
  real*8,dimension(:,:),allocatable :: rfcpc,brgrfc,bzgrfc, &
                                     rsilpc,rmp2pc,rgowpc, &
                                     gridec,gridvs,ggridvs, &
                                     gridac

  if (.not. allocated(rfcpc)) then
    allocate(rfcpc(nfcoil,nwnh))
    rfcpc(:,:) = 0.0
  endif
  if (.not. allocated(brgrfc)) then
    allocate(brgrfc(nwnh,nfsum))
    brgrfc(:,:) = 0.0
  endif
  if (.not. allocated(bzgrfc)) then
    allocate(bzgrfc(nwnh,nfsum))
    bzgrfc(:,:) = 0.0
  endif
  if (.not. allocated(rsilpc)) then
    allocate(rsilpc(nsilop,nwnh))
    rsilpc(:,:) = 0.0
  endif
  if (.not. allocated(rmp2pc)) then
    allocate(rmp2pc(magpr2,nwnh))
    rmp2pc(:,:) = 0.0
  endif
  if (.not. allocated(rgowpc)) then
    allocate(rgowpc(nrogow,nwnh))
    rgowpc(:,:) = 0.0
  endif
  if (.not. allocated(gridec)) then
    allocate(gridec(nwnh,nesum))
    gridec(:,:) = 0.0
  endif
  if (.not. allocated(gridvs)) then
    allocate(gridvs(nwnh,nvesel))
    gridvs(:,:) = 0.0
  endif
  if (.not. allocated(ggridvs)) then
    allocate(ggridvs(nwnh,nvsum))
    ggridvs(:,:) = 0.0
  endif
  if (.not. allocated(gridac)) then
    allocate(gridac(nwnh,nacoil))
    gridac(:,:) = 0.0
  endif
!
  n000=0
  if (ifcoil <= 0) go to 1100
!----------------------------------------------------------------------
!-- calculate the response function of psi loops due to f coils      --
!----------------------------------------------------------------------
  do i=1,nfcoil
    taf(i)=tan(af(i)*pi/180.)
    taf2(i)=tan(af2(i)*pi/180.)
  enddo 
  if (islpfc <= 0) go to 300
  do i=1,nfcoil
    do j=1,nfcoil
      call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j), &
                rfcfc(j,i))
      rfcfc(j,i)=rfcfc(j,i)*0.5/pi
    enddo 
    ii=i
    call gsilop(rgrid,nw,zgrid,nh,rfcpc,ii,rf,zf,wf,hf,af,af2,nfcoil)
  enddo
!vas
  print*,'file name : ','fc'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrspfc,status='unknown',file='fcfcpc.dat', &
  open(unit=nrspfc,status='unknown',file='fc'//trim(ch1)// & 
                         trim(ch2)//'.ddd',form='unformatted')
  write (nrspfc) rfcfc
  write (nrspfc) rfcpc
  close(unit=nrspfc)
 300  continue  
  msilop=nsilop
  if (msilop <= 1) go to 520
  if (isize <= 0) go to 420
!---------------------------------------------------------------------
!-- finite size flux loops                                          --
!---------------------------------------------------------------------
  do i=1,nfcoil
    do j=1,isize
      taz=tan(as(j)*pi/180.)
      taz2=tan(as2(j)*pi/180.)
      call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                rsi(j),zsi(j),wsi(j),hsi(j),taz,taz2, &
                rsilfc(j,i))
      rsilfc(j,i)=rsilfc(j,i)*0.5/pi
    enddo 
  enddo
  420 continue
  if (isize >= nsilop) go to 520
!---------------------------------------------------------------------
!-- thin flux loops                                                 --
!---------------------------------------------------------------------
  do i=1,nfcoil
    ii=i
    do j=isize+1,nsilop
      jj=j
      call m1coef(xdum,xdum,n000,n000,rsilfc,jj,ii)
    enddo 
  enddo 
  520 continue
!----------------------------------------------------------------------
!-- compute the response function of magnetic probes due to f coils  --
!----------------------------------------------------------------------
  magprr=magpr2
  if (magprr > 1) then
    call m2coef(xdum,n000,ydum,n000,rmp2fc,magpr2,nfcoil)
  endif
!----------------------------------------------------------------------
!-- compute the response function of partial rogowski loops due to   --
!-- f coils                                                          --
!----------------------------------------------------------------------
  mrogow=nrogow
  if (mrogow > 1) then
     call rogowc(xdum,n000,ydum,n000,rgowfc,nrogow,nfcoil)
  endif
!----------------------------------------------------------------------
!-- write f coil response functions                                  --
!----------------------------------------------------------------------
  do i=1,nfsum
    do j=1,nsilop
      gsilfc(j,i)=0.0
    enddo
    do j=1,magpr2
      gmp2fc(j,i)=0.0
    enddo
  enddo
!
  do i=1,nfcoil
    k=abs(fcid(i))
    do j=1,nsilop
      gsilfc(j,k)=gsilfc(j,k)+fcturn(i)*rsilfc(j,i)
    enddo
    do j=1,magpr2
      gmp2fc(j,k)=gmp2fc(j,k)+fcturn(i)*rmp2fc(j,i)
    enddo
  enddo
!
!vas
  print*,'file name : ','rfcoil.ddd' 
!vasorg      open(unit=nrspfc,status='unknown',file='rfcoil.dat', &
  open(unit=nrspfc,status='unknown',file='rfcoil.ddd', &
       form='unformatted')
  write (nrspfc) gsilfc
  write (nrspfc) gmp2fc
  close(unit=nrspfc)
!
  do i=1,nfsum
    do j=1,nwnh
      brgrfc(j,i)=0.0
      bzgrfc(j,i)=0.0
    enddo
  enddo
  do i=1,nfcoil
    k=abs(fcid(i))
    do j=1,nwnh
      brgrfc(j,k)=brgrfc(j,k)+fcturn(i)*brgridfc(j,i)
      bzgrfc(j,k)=bzgrfc(j,k)+fcturn(i)*bzgridfc(j,i)
    enddo
  enddo
!
! -- Disable brzgfc.dat output by lzp 2015/11/12
! open(unit=nrspfc,status='unknown',file='brzgfc.dat', &
!      form='unformatted')
! write (nrspfc) brgrfc
! write (nrspfc) bzgrfc
! close(unit=nrspfc)
 1100 continue
!----------------------------------------------------------------------
!-- plasma response functions                                        --
!----------------------------------------------------------------------
  if (igrid <= 0) go to 3200
  msilop=nsilop
  if (msilop <= 1) go to 1220
!----------------------------------------------------------------------
!-- filament plasma current model                                    --
!----------------------------------------------------------------------
  if (isize <= 0) go to 1160
  do j=1,isize
    jj=j
    call gsilop(rgrid,nw,zgrid,nh,rsilpc,jj,rsi,zsi,wsi,hsi,as,as2 &
      ,nsilop)
  enddo 
 1160 continue
  f (isize >= nsilop) go to 1220
  do j=isize+1,nsilop
    jj=j
    call m1coef(rgrid,zgrid,nw,nh,rsilpc,jj,n000)
  enddo 
 1220 continue
  magprr=magpr2
  if (magprr > 1) then
    call m2coef(rgrid,nw,zgrid,nh,rmp2pc,magpr2,nwnh)
  endif
  mrogow=nrogow
  if (mrogow > 1) then
    call rogowc(rgrid,nw,zgrid,nh,rgowpc,nrogow,nwnh)
  endif
!----------------------------------------------------------------------
!-- write the plasma response function                               --
!----------------------------------------------------------------------
!vas
      print*,'file name : ','ep'//trim(ch1)// & 
                        trim(ch2)//'.ddd' 
!vasorg      open(unit=nrsppc,status='unknown',file='eplasm.dat', &
      open(unit=nrsppc,status='unknown',file='ep'//trim(ch1)// & 
                         trim(ch2)//'.ddd' , &
           form='unformatted')
      write (nrsppc) rsilpc
      write (nrsppc) rmp2pc
      close(unit=nrsppc)
!
 3200 continue
      if (iecoil > 0) then
         call gecoil(rsilec,rmp2ec,gridec,rgrid,nw, &
                     zgrid,nh,rfcec,recec,rsisec)
      endif
!vas
      print*,'file name : ','re'//trim(ch1)// & 
                        trim(ch2)//'.ddd' 
!vasorg      open(unit=nrsppc,status='unknown',file='recoil.dat', &
      open(unit=nrsppc,status='unknown',file='re'//trim(ch1)// & 
                     trim(ch2)//'.ddd', &
           form='unformatted')
      write (nrsppc) rsilec
      write (nrsppc) rmp2ec
      write (nrsppc) gridec
      close(unit=nrsppc)
!
      if (ivesel > 0) then
         call gvesel(rsilvs,rmp2vs,gridvs,rgrid,nw, &
                     zgrid,nh,rfcvs,rvsfc,rvsec)
         do i=1,nvesel
           tas(i)=tan(avs(i)*pi/180.)
           tas2(i)=tan(avs2(i)*pi/180.)
         enddo 
         do i=1,nvesel
            do j=1,nvesel
               call flux(rvs(i),zvs(i),wvs(i),hvs(i),tas(i),tas2(i), &
                         rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j), &
                         rvsvs(j,i))
               rvsvs(j,i)=rvsvs(j,i)*0.5/pi
            enddo 
         enddo 
!
         do i=1,nvsum
            do j=1,nsilop
               gsilvs(j,i)=0.0
            enddo
            do j=1,magpr2
               gmp2vs(j,i)=0.0
            enddo
            do j=1,nwnh
               ggridvs(j,i)=0.0
            enddo
         enddo
!
         do i=1,nvesel
            k=abs(vsid(i))
            do j=1,nsilop
               gsilvs(j,k)=gsilvs(j,k)+rsilvs(j,i)
            enddo
            do j=1,magpr2
               gmp2vs(j,k)=gmp2vs(j,k)+rmp2vs(j,i)
            enddo
            do j=1,nwnh
               ggridvs(j,k)=ggridvs(j,k)+gridvs(j,i)
            enddo
         enddo
!
!vas
      print*,'file name : ','rv'//trim(ch1)// & 
                        trim(ch2)//'.ddd' 
!vasorg         open(unit=nrsppc,status='unknown',file='rvesel.dat', &
         open(unit=nrsppc,status='unknown',file='rv'//trim(ch1)// & 
                          trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) gsilvs
         write (nrsppc) gmp2vs
         write (nrsppc) ggridvs
         close(unit=nrsppc)
      endif
!---------------------------------------------------------------------
!-- advance divertor coil                                           --
!---------------------------------------------------------------------
      if (iacoil > 0) then
         call gacoil(rsilac,rmp2ac,gridac,rgrid,nw, &
                     zgrid,nh)
!vas
      print*,'file name : ','ra'//trim(ch1)// & 
                        trim(ch2)//'.ddd' 
!vasorg         open(unit=nrsppc,status='unknown',file='racoil.dat', &
         open(unit=nrsppc,status='unknown',file='ra'//trim(ch1)// & 
                         trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) gridac
         write (nrsppc) rsilac
         write (nrsppc) rmp2ac
         close(unit=nrsppc)
      endif
!
      if ( allocated(rfcpc)) deallocate(rfcpc)
      if ( allocated(brgrfc)) deallocate(brgrfc)
      if ( allocated(bzgrfc)) deallocate(bzgrfc)
      if ( allocated(rsilpc)) deallocate(rsilpc)
      if ( allocated(rmp2pc)) deallocate(rmp2pc)
      if ( allocated(rgowpc)) deallocate(rgowpc)
      if ( allocated(gridec)) deallocate(gridec)
      if ( allocated(gridvs)) deallocate(gridvs)
      if ( allocated(ggridvs)) deallocate(ggridvs)
      if ( allocated(gridac)) deallocate(gridac)
!
      return
      end subroutine efund_matrix


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          grid computes the green's functions at (r,z)            **
!**          due to plasma currents and f coils.                     **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          09/06/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine efund_grid
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
  use consta
  use pmodel
  use input
  use fcoil
  use nio

  use var_filech
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  real*8,dimension(:,:),allocatable :: gridfc,gridpc,ggridfc
  dimension taf(nfcoil),taf2(nfcoil)
  data isplit/17/,tole/1.0e-10/

  if (.not. allocated(gridfc)) then
    allocate(gridfc(nwnh,nfcoil))
    gridfc(:,:) = 0.0
  endif
  if (.not. allocated(gridpc)) then
    allocate(gridpc(nwnh,nw))
    gridpc(:,:) = 0.0
  endif
  if (.not. allocated(ggridfc)) then
    allocate(ggridfc(nwnh,nfsum))
    ggridfc(:,:) = 0.0
  endif
!
  if (igrid <= 0) return
!----------------------------------------------------------------------
!--  compute the green's functions at (r,z) due to f coils           --
!----------------------------------------------------------------------
  ndr = isplit
  ndz = isplit
  do in=1,nfcoil
    taf(in) = tan(af(in)*pi/180.0)
    taf2(in) = tan(af2(in)*pi/180.0)
    do ni=1,nw
      do nj=1,nh
        kk = (ni-1)*nh + nj
        rsum = 0.
        dwc = wf(in)/ndr
        dhc = hf(in)/ndz
        if (af2(in) ~= 0.) go to 640
        z1 = zf(in) - taf(in)*(wf(in)-dwc)/2. - .5*hf(in) + .5*dhc
        ab = rf(in) - .5*wf(in) + .5*dwc
        do iw = 1, isplit
          drc = ab + (iw-1)*dwc + iw*tole
          z2 = z1 + (iw-1)*taf(in)*dwc
          do ih = 1, isplit
            dzc = z2 + (ih-1)*dhc
            rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
            rsum = rsum + rtmp
          enddo 
        enddo
        go to 660
  640 continue
        do ih = 1, ndz
          dzc = zf(in) - .5*hf(in) + .5*dhc + dhc*(ih-1)
          do iw = 1, ndr
            drc = rf(in) - .5*wf(in) - .5*hf(in)/taf2(in) &
                  + .5*dwc + .5*dhc/taf2(in) &
                  + dhc/taf2(in)*(ih-1) + dwc*(iw-1)
            rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
            rsum = rsum + rtmp
          enddo 
        enddo
  660 continue
        cmut = rsum*2.e-07/(isplit*isplit)
        gridfc(kk,in) = cmut
      enddo
    enddo
  enddo
!----------------------------------------------------------------------
!--  compute the green's functions at (r,z) due to itself            --
!----------------------------------------------------------------------
  aaa=0.0
  do i=1,nw
    do j=1,nh
      kk=(i-1)*nh+j
      do ni=1,nw
        if ((j > 1).or.(i ~= ni)) then
          zdif=(j-1)*dz
          gridpc(kk,ni)=psical(rgrid(i),rgrid(ni),zdif)*tmu
        else
          call flux(rgrid(ni),aaa,dr,dz,aaa,aaa,rgrid(ni),aaa,&
                    dr,dz,aaa,aaa,fridpc)
          gridpc(kk,ni)=fridpc*0.5/pi
        endif
      enddo 
    enddo
  enddo
!----------------------------------------------------------------------
!--  store green's function table                                    --
!----------------------------------------------------------------------
!
  do i=1,nfsum
    do j=1,nwnh
      ggridfc(j,i)=0.0
    enddo
  enddo
!
  do i=1,nfcoil
    k=abs(fcid(i))
    do j=1,nwnh
      ggridfc(j,k)=ggridfc(j,k)+fcturn(i)*gridfc(j,i)
    enddo
  enddo
!vas
  print*,'file name : ','ec'//trim(ch1)//trim(ch2)//'.ddd' 
!
!vasorg      open(unit=ncontr,status='unknown',file='econto.dat', &
  open(unit=ncontr,status='unknown',file='ec'//trim(ch1)// &
       trim(ch2)//'.ddd',form='unformatted')
  mw=nw
  mh=nh
  write (ncontr) mw,mh
  write (ncontr) rgrid,zgrid
  write (ncontr) ggridfc
  write (ncontr) gridpc
!vas ... just for testing
!      open(35,file='test-ec1.dat',status='new')
!      write (35,*) mw,mh
!      write (35,1009) rgrid,zgrid
!      close(35)
!      open(35,file='test-ec2.dat',status='new')
!      write (35,1009) ggridfc
!      close(35)
!      open(35,file='test-ec3.dat',status='new')
!      write (35,1009) gridpc
!      close(35)
!1009  format(3(1x,e14.8))
!vas ...
  close(unit=ncontr)
!
  if (allocated(gridfc)) deallocate(gridfc)
  if (allocated(ggridfc)) deallocate(ggridfc)
  if (allocated(gridpc)) deallocate(gridpc)
!
  return
end subroutine efund_grid


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          for SCCS control revision information.                  **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          11/07/95..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine efundu_rev(i)
  character*100 opt
  character*10 s
  if(i==0)  &
  s='@(#)efund.for,v 2.3 1996/10/17 15:53:28 lao Exp\000'
  return
end subroutine efundu_rev
