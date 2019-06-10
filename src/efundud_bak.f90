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
!**          gacoil gets the Green's functions for the advance       **
!**          divertor coil.                                          **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine gacoil(rsilac,rmp2ac,gridac,rgrid,mw,zgrid,mh)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                  nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
  use consta
  use cacoil
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension rsilac(nsilop,nacoil),rmp2ac(magpr2,nacoil)
  real*8,dimension(mw) ::  rgrid
  real*8,dimension(mh) ::  zgrid
  real*8,dimension(nwnh,nacoil) :: gridac
  do j=1,nsilop
    jj=j
    do i=1,nacoil
      ii=i
      call a1coef(work,jj,ii)
      rsilac(j,i)=work
    enddo 
  enddo 

  do j=1,magpr2
    jj=j
    do i=1,nacoil
      ii=i
      call a2coef(work,jj,ii)
      rmp2ac(j,i)=work
    enddo 
  enddo 

  do i=1,nw
    nr=i
    do j=1,nh
      nz=j
      kk=(i-1)*nh+j
      do n=1,nacoil
        nn=n
        call agrid(work,rgrid,nr,zgrid,nz,nn)
        gridac(kk,n)=work
      enddo 
    enddo
  enddo 
  return
end subroutine gacoil


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
!**          flux computes the mutual inductance/2/pi between        **
!**          two circulars of rectangular cross section.             **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**       r1,r2...........radius of first and second coil            **
!**       z1,z2...........elevation                                  **
!**       w1,w2...........width                                      **
!**       h1,h2...........height                                     **
!**                                                                  **
!**                                                                  **
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
subroutine flux(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,fuxx)
  use exparm,only:mgaus1,mgaus2
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension post1(mgaus1),wght1(mgaus1),post2(mgaus2) &
           ,wght2(mgaus2)
  data init/0/

  if (init > 0) go to 100
!vas introduced to make it ok in hydra
!org      ngaus1=mgaus1
!org      ngaus2=mgaus2
!org      call lgauss(post1,wght1,ngaus1,ner)
!org      call lgauss(post2,wght2,ngaus2,ner)
!      print*,'initializing',init,mgaus1,mgaus2
  call lgauss(post1,wght1,mgaus1,ner)
  call lgauss(post2,wght2,mgaus2,ner)
  init=1
  100 continue
!
  x = 0.
  fuxx = 0.
  zf = z1
  zs = z2
  hf2 = h1*.5
  hs2 = h2*.5

  if (t12==0) go to 200
  xl1 = r1-0.5*w1-0.5*h1/abs(t12)
  xr1 = r1+0.5*w1+0.5*h1/abs(t12)
  xbm1 = xl1+w1
  xtm1 = xr1-w1
  if (t12 < 0.) xbm1 = xr1 - w1
  if (t12 < 0.) xtm1 = xl1 + w1
!
  200 if (t22==0) go to 300
  xl2 = r2-0.5*w2-0.5*h2/abs(t22)
  xr2 = r2+0.5*w2+0.5*h2/abs(t22)
  xbm2 = xl2+w2
  xtm2 = xr2-w2
  if (t22 < 0.) xbm2 = xr2 - w2
  if (t22 < 0.) xtm2 = xl2 + w2
!
  300 continue
!      print*,'vasan : ',init,mgaus1,mgaus2
  do i = 1,mgaus1
!org      do i = 1,ngaus1
    rf = r1+.5*w1*post1(i)
    if (t12 ~= 0) rf = r1+(0.5*w1+0.5*h1/abs(t12))*post1(i)
    do j = 1,mgaus2
!org         do j = 1,ngaus2
      rs = r2+0.5*w2*post2(j)
      if (t22 ~= 0) rs = r2+(0.5*w2+0.5*h2/abs(t22))*post2(j)
      call soleno(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2,&
                  xtm1,xtm2,hfa,hsa,rf,rs,solx)
      fuxx = fuxx+wght1(i)*wght2(j)/hfa/hsa*solx
    enddo 
  enddo 
!
return
end subroutine flux


      subroutine gsilop(rr, nr, zz, nz, rspfun, ns, rsi, zsi, wsi &
           , hsi, as, as2, ndim)
!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          gsilop computes the green's function at si loops due    **
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
      subroutine gvesel(rsilvs,rmp2vs,gridvs,rgrid,mw, &
                        zgrid,mh,rfcvs,rvsfc,rvsec)
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
      subroutine lgauss(x,w,n,nn)
!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          lgauss computes the zeroes of the legendre polynomial   **
!**          and their associated weights for a gaussian quadrature. **
!**                                                                  **
!**     calling arguments:                                           **
!**       x...............zeroes of legendre polynomial              **
!**       w...............weights                                    **
!**       n...............order of legendre polynomial               **
!**       nn..............error flag                                 **
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
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(n) ::  x,w
!
      nn = 0
      if (n-1) 10,20,30
   10 continue
      nn = 1
      return
!----------------------------------------------------------------------
!-- request for a zero point formula is meaningless                  --
!----------------------------------------------------------------------
   20 x(1) = 0.
      w(1) = 2.
      return
!----------------------------------------------------------------------
!-- for a one point formula, send back                               --
!-- results without computing.                                       --
!----------------------------------------------------------------------
   30 r = n
      g = -1.
!----------------------------------------------------------------------
!-- the initial guess for the smallest root                          --
!-- of p(n) is taken as -1.                                          --
!----------------------------------------------------------------------
      do i = 1,n
      test = -2.
      ic = n+1-i
!----------------------------------------------------------------------
!-- whenever we find a root of the                                   --
!-- polynomial, its negative is also a root.                         --
!-- the index ic tells WHERE to store the other root                 --
!----------------------------------------------------------------------
      if (ic < i) go to 950
   40 s = g
      t = 1.
      u = 1.
      v = 0.
!----------------------------------------------------------------------
!-- evaluation of the n-th legendre polynomial                       --
!-- and its first derivative.                                        --
!-- WHERE   u = ds/dx                                                --
!--         v = dt/dx                                                --
!--         dp=dp/dx                                                 --
!----------------------------------------------------------------------
      do k = 2,n
      a = k
      p = ((2.0*a-1.0)*s*g-(a-1.0)*t)/a
      dp = ((2.0*a-1.0)*(s+g*u)-(a-1.0)*v)/a
      v = u
      u = dp
      t = s
      s = p
      enddo 
      if (abs((test-g)/(test+g)) < 0.0000005) go to 100
      sum = 0.
      if (i==1) go to 70
!----------------------------------------------------------------------
!-- the following computes the reduced                               --
!-- legendre polynomial and its derivative.                          --
!----------------------------------------------------------------------
      do k = 2,i
         sum = sum+1./(g-x(k-1))
      enddo 
   70 test = g
      g = g-p/(dp-p*sum)
      go to 40
  100 x(ic) = -g
      x(i) = g
      w(i) = 2./(r*t*dp)
      w(ic) = w(i)
  900 g = g-r*t/((r+2.)*g*dp+r*v-2.*r*t*sum)
	    enddo
  950 return
      end subroutine lgauss
      subroutine m1coef(rr, zz, nr, nz, coef,  nl, nf)
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
      subroutine m2coef(rr, nr, zz, nz, coef,  mp, nc)
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


      function psical(a1,r1,z1)
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
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8 x1,cay,ee,xmdelk,xmdele
!
      isw=1
      go to 10
      ENTRY br(a1,r1,z1)
      isw=2
      go to 10
      ENTRY bz(a1,r1,z1)
      isw=3
!
   10 continue
      a=a1
      r=r1
      z=z1
      den=a*a+r*r+z*z+2.*a*r
      xk=4.*a*r/den
      x1=(a*a+r*r+z*z-2.*a*r)/den
      if (x1 < 1.0e-10) x1=1.0e-10
      cay=xmdelk(x1)
      ee=xmdele(x1)
      go to (20,30,40),isw
!----------------------------------------------------------------------
!--   psi computation                                                --
!----------------------------------------------------------------------
   20 psical= sqrt(den)*((1.e+00-0.5e+00*xk)*cay-ee)
      return
!----------------------------------------------------------------------
!--   br  computation                                                --
!----------------------------------------------------------------------
   30 psical=z/(r* sqrt(den))*(-cay+(a*a+r*r+z*z)/((a-r)*(a-r)+z*z)*ee)
      return
!----------------------------------------------------------------------
!--   bz  computation                                                --
!----------------------------------------------------------------------
   40 psical=(cay+(a*a-r*r-z*z)/((a-r)*(a-r)+z*z)*ee)/ sqrt(den)
      return
      end function psical
      subroutine rogowc(rr, nrr, zz, nz, coef, nr, nc)
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
      subroutine rogrid(ngrid,mm,m,dels)
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
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      use rogowl
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension sl(6)
!
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
      subroutine simpf(i,f)
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
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      if (i==1 .or. i==25) then
         f = 1./3.
      elseif ((i/2)*2==i) then
         f = 4./3.
      else
         f = 2./3.
      endif
      return
      end subroutine simpf
      subroutine soleno(ra,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,sol)
!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          soleno computes the inductance for a solenoid           **
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
      use consta,only:pi
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension z(2,2)
      real*8 ksq,kpsq,msl,mut
      data init/0/
!
      if (init > 0) go to 5
      rpi=pi
      rpi2=rpi*0.5
      rh=rpi*2.0e-07
      ut=rh/3.0
      err=1.0e-05
      init=1
    5 continue
!
      dr = rf-rs
      drsq = dr*dr
      tr = rf+rs
      trsq = tr*tr
      fr = 4.*rf*rs
      sol = 0.
      csq = fr/trsq
      cpsq = 1.-csq
!
   10 if (t12 ~= 0.) go to 20
      r = rf-ra+0.5*w1
      zc1 = z1-0.5*h1-0.5*w1*t1
      zb1 = zc1+t1*r
      zt1 = zb1+h1
      go to 30
!
   20 zb1 = z1-h1/2.
      if (t12 < 0.) go to 25
      if (rf > xbm1) zb1 = zb1+t12*(rf-xbm1)
      zt1 = z1+h1/2.
      if (rf < xtm1) zt1 = zt1-t12*(xtm1-rf)
      go to 30
!
   25 if (rf < xbm1) zb1 = zb1 + t12*(rf-xbm1)
      zt1 = z1 + h1/2.
      if (rf > xtm1) zt1 = zt1 - t12*(xtm1-rf)
!
   30 if (t22 ~= 0.) go to 40
      r = rs-r2+0.5*w2
      zc2 = z2-0.5*h2-0.5*w2*t2
      zb2 = zc2+t2*r
      zt2 = zb2+h2
      go to 50
!
   40 zb2 = z2-h2/2.
      if (t22 < 0.) go to 45
      if (rs > xbm2) zb2 = zb2+t22*(rs-xbm2)
      zt2 = z2+h2/2.
      if (rs < xtm2) zt2 = zt2-t22*(xtm2-rs)
      go to 50
!
   45 if (rs < xbm2) zb2 = zb2 + t22*(rs-xbm2)
      zt2 = z2 + h2/2.
      if (rs > xtm2) zt2 = zt2 - t22*(xtm2-rs)
!
   50 z(1,1) = zb1
      z(2,1) = zb2
      z(1,2) = zt1
      z(2,2) = zt2
      hfa = zt1-zb1
      hsa = zt2-zb2
!
      do i = 1,2
      do j = 1,2
      sign = -.25
      if (i ~= j) sign = .25
      dz = z(1,i)-z(2,j)
      dzsq = dz*dz
      r2sq = dzsq+drsq
      r1sq = dzsq+trsq
      r1 = sqrt(r1sq)
      ksq = fr/r1sq
      t = 2./ksq-1.
      kpsq = 1.-ksq
      alpha = 1.
!--------------------------------------------------------------------------
!--  To avoid numerical truncation                                       --
!--------------------------------------------------------------------------
      if (kpsq < 1.0e-30) kpsq = 1.0e-30
      beta = sqrt(kpsq)
      if (beta < 1.0e-30) beta = 1.0e-10
      if (cpsq < 1.0e-30) cpsq = 1.0e-10
      delta = cpsq/beta
      epsi = csq/cpsq
      zeta = 0.
      sinf = 0.
      sa = .25
!
  100 continue
      sa = 2.*sa
      ambsq = (alpha-beta)*(alpha-beta)
      sinf = sinf+sa*ambsq
      alphat = alpha
      epsit = epsi
      alpha = .5*(alpha+beta)
      beta = sqrt(alphat*beta)
      epsi = (delta*epsi+zeta)/(1.+delta)
      delta = beta/4./alpha*(2.+delta+1./delta)
      zeta = .5*(epsit+zeta)
      if (abs(delta-1.) > err) go to 100
      if (ambsq > 1.e-14) go to 100
      cay = rpi2/alpha
      pik = cay*zeta
      ek = .5*cay*(ksq+sinf)
      msl = rh*dzsq*(r1*ek-drsq*pik/r1)
      if (csq-1.) 290,190,290
  190 msl = rh*dzsq*(r1*ek-cay*fr/r1*.5)
  290 continue
!
      mut = msl+ut*fr*r1*(cay-t*ek)
      sol = sol+sign*mut
      enddo 
      enddo 
!
      return
      end subroutine soleno
      subroutine splitc(is,rs,zs,cs,rc,zc,wc,hc,ac,ac2,cc)
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
      use consta
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(is) :: rs,zs,cs
!
      frd=pi/180.
      if(ac+ac2==0.) go to 100
      if(ac ~= 0.) go to 200
      if(ac2 ~= 0.) go to 300
!----------------------------------------------------------------------
!-- rectangle                                                        --
!----------------------------------------------------------------------
  100 continue
      wdelt=wc/is
      hdelt=hc/is
      rstrt=rc-wc/2.+wdelt/2.
      zstrt=zc-hc/2.+hdelt/2.
      zz=zstrt
      ic=0
      c=cc/(is*is)
      do ii=1,is
         rr=rstrt
         do jj=1,is
            ic=ic+1
            zs(ic)=zz
            rs(ic)=rr
            cs(ic)=c
            rr=rr+wdelt
         enddo 
         zz=zz+hdelt
      enddo
      go to 900
!----------------------------------------------------------------------
!-- ac ~= 0                                                        --
!----------------------------------------------------------------------
  200 continue
      side=tan(frd*ac)*wc
      hdelt=hc/is
      wdelt=wc/is
      zdelt=tan(frd*ac)*wdelt
      rstrt=rc-wc/2.+wdelt/2.
      tsid=hc+side
      zstrt =zc-tsid/2.+tsid/2.*1./is
      rr=rstrt
      ic=0
      c=cc/(is*is)
      do ii=1,is
         zz=zstrt+(ii-1)*zdelt
         do jj=1,is
            ic=ic+1
            zs(ic)=zz
            rs(ic)=rr
            cs(ic)=c
            zz=zz+hdelt
         enddo 
         rr=rr+wdelt
      enddo
      go to 900
!----------------------------------------------------------------------
!-- ac2 ~= 0                                                         --
!----------------------------------------------------------------------
  300 continue
!
  340 continue
      side=hc/tan(frd*ac2)
      hdelt=hc/is
      wdelt=wc/is
      zstrt=zc-hc/2.+hdelt/2.
      rdelt=hdelt/tan(frd*ac2)
      rstrt=rc-side/2.-wc/2.+rdelt/2.+wdelt/2.
      side=hc/tan(frd*ac2)
      wtot=side+wc
      whaf=(side+wc)/2.
      rcorn=rc-whaf
      rcornr=rc+whaf
      rcorn2=rcorn+wtot/is
      rstrt=(rcorn+rcorn2)/2.
      zz=zstrt
      ic=0
      c=cc/(is*is)
      do ii=1,is
         rr=rstrt+(ii-1)*rdelt
         do jj=1,is
            ic=ic+1
            zs(ic)=zz
            rs(ic)=rr
            cs(ic)=c
            rr=rr+wdelt
         enddo 
         zz=zz+hdelt
      enddo
      go to 900
!
  900 continue
      return
      end subroutine splitc
      subroutine v1coef(coef,  nl, ne)
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
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use consta
      use nio
      use siloop
      use cvesel
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
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
!
      return
      end subroutine v1coef


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          e2coef computes the response functions due to           **
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


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          egrid computes the response functions due to            **
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
!**          a1coef computes the response functions due to           **
!**          the thin flux loops and the advance divertor coil.      **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine a1coef(coef, nl, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                  nesum,nfsum,nvsum,nvesel,nacoil
  use coilsp
  use consta
  use cacoil
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
              racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb,cdum)
  do l=1,itot
     a=rsplt(l)
     r1=rsi(m)
     z1=zsi(m)-zsplt(l)
     psic=psical(a,r1,z1)*tmu
     psict=psict+psic/fitot
  enddo 
  coef=psict

  return
end subroutine a1coef


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          a2coef computes the response functions due to           **
!**          the magnetic probes and the advance divertor coil.      **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine a2coef(coef, mp, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                  nesum,nfsum,nvsum,nvesel,nacoil
  use coilsp
  use mprobe
  use consta
  use cacoil
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
              racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb,cdum)
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
end subroutine a2coef


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          agrid computes the response functions due to            **
!**          the grid points and the advance divertor coil.          **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine agrid(coef, rgrid, nr, zgrid, nz, ne)
  use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                  nesum,nfsum,nvsum,nvesel,nacoil
  use coilsp
  use consta
  use cacoil
  use nio
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  real*8,dimension(nr) :: rgrid
  real*8,dimension(nr) :: zgrid
  data init/0/

  radeg=pi/180.
  isplit=17
  itot=isplit*isplit
  fitot=itot

  k=ne
  psict=0
  aaa=0.0
  bbb=0.0
  call splitc(isplit,rsplt,zsplt,csplt, &
              racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb,cdum)
  do l=1,itot
    a=rsplt(l)
    r1=rgrid(nr)
    z1=zgrid(nz)-zsplt(l)
    psic=psical(a,r1,z1)*tmu
    psict=psict+psic/fitot
  enddo 
  coef=psict

  return
end subroutine agrid


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**       xmdele computes the elliptic integral e.                   **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**       xm1.............argument of elliptic integral e            **
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
function xmdele(xm1)
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension a(4),b(4)
  real*8 a,b,xm1,xmdele
  data a(1),a(2),a(3),a(4)/.44325141463,.06260601220,&
    .04757383546,.01736506451/
  data b(1),b(2),b(3),b(4)/.24998368310,.09200180037,&
    .04069697526,.00526449639/

  xmdele=1.0+xm1*(a(1)+xm1*(a(2)+xm1*(a(3)+xm1*a(4))))&
   +xm1*(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*b(4))))*log(1.0/xm1)
  return
end function xmdele


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**       xmdelk computes the elliptic integral k.                   **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**       xm1.............argument of elliptic integral k            **
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
function xmdelk(xm1)
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension a(5),b(5)
  real*8  a,b,xm1,xmdelk
  data a(1),a(2),a(3),a(4),a(5)/1.38629436112,.09666344259,&
    .03590092383,.03742563713,.01451196212/
  data b(1),b(2),b(3),b(4),b(5)/.5,.12498593597,.06880248576,&
    .03328355346,.00441787012/

  xmdelk=a(1)+xm1*(a(2)+xm1*(a(3)+xm1*(a(4)+xm1*a(5))))&
   +(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*(b(4)+xm1*b(5)))))&
   *log(1.0/xm1)
  return
end function xmdelk


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
