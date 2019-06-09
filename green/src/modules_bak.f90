!***********************************************************************
! module for constants                                                 *
!***********************************************************************
module consta
  implicit none
  public

  real*8,parameter :: pi = 3.1415926535897932
  real*8,parameter :: tmu = 2.0e-07

end module consta


!***********************************************************************
! module for experimental parameters                                   *
!***********************************************************************
module exparm
  implicit none
  public
  
  !> number of p.f. coils                                      
  integer*4,parameter :: nfcoil = 6
  !> number of flux loops                                      
  integer*4,parameter :: nsilop = 41
  !> total number of magnetic probes (magpri in EFIT)           
  integer*4,parameter :: magpr2 = 63
  !> number of partial rogowski loops                            
  integer*4,parameter :: nrogow = 1
  !> number of ohmic heating coils                             
  integer*4,parameter :: necoil = 1
  !> number of o.h. coil groups                                
  integer*4,parameter :: nesum = 1
  !> number of p.f. coil groups                                
  integer*4,parameter :: nfsum = 6
  !> number of vessel segement groups                          
  integer*4,parameter :: nvsum = 24
  !> number of vessel segements     
  integer*4,parameter :: nvesel = 24
  !> number of advance divertor coils                          
  integer*4,parameter :: nacoil = 1
  integer*4,parameter :: mgaus1 = 8
  integer*4,parameter :: mgaus2 = 10
  integer*4 :: nw,nh,nwnh
end module exparm


!***********************************************************************
! module for input                                                     *
!***********************************************************************
module input
  implicit none
  public

  !> flag for shapping coil calculation, =1 enabled, =0 disabled
  integer*4 :: ifcoil
  integer*4 :: islpfc
  !> flag for Ohmic coil calculation, =1 enabled, =0 disabled
  integer*4 :: iecoil
  integer*4 :: isize
  !> flag for vessel calculation, =1 enabled, =0 disabled
  integer*4 :: ivesel
  integer*4 :: igrid
  real*8    :: rleft,rright,zbotto,ztop

end module input


!***********************************************************************
! module for bfgrid                                                    *
!***********************************************************************
module bfgrid
  implicit none
  public

  real*8,dimension(:,:),allocatable,save :: brgridfc,bzgridfc

end module bfgrid


!***********************************************************************
! module for poloidal field coil                                       *
!***********************************************************************
module fcoil
  use exparm,only:nfcoil
  implicit none
  public

  !> R coordinate of poloidal field coil
  real*8,dimension(nfcoil) :: rf
  !> Z coordinate of poloidal field coil
  real*8,dimension(nfcoil) :: zf
  !> width of poloidal field coil
  real*8,dimension(nfcoil) :: wf
  !> height of poloidal field coil
  real*8,dimension(nfcoil) :: hf
  !> angle 1 of poloidal field coil
  real*8,dimension(nfcoil) :: af
  !> angle 2 of poloidal field coil
  real*8,dimension(nfcoil) :: af2
  real*8,dimension(nfcoil) :: turnfc
  real*8,dimension(nfcoil) :: fcid
  real*8,dimension(nfcoil) :: fcturn

end module fcoil


!***********************************************************************
! module for advanced divertor coil                                    *
!***********************************************************************
module cacoil
  use exparm,only: nacoil
  implicit none
  public

  !> flag for advance divertor coil calculation, =1 enabled, =0 disabled
  integer*4 :: iacoil
  !> R coordinate of advanced divertor coils
  real*8,dimension(nacoil) :: racoil
  !> Z coordinate of advanced divertor coils
  real*8,dimension(nacoil) :: zacoil
  !> width of advanced divertor coils
  real*8,dimension(nacoil) :: wacoil
  !> height of advanced divertor coils
  real*8,dimension(nacoil) :: hacoil

end module cacoil


!***********************************************************************
! module for Ohmic heating coil                                        *
!***********************************************************************
module cecoil
  use exparm,only:necoil
  implicit none
  public

  !> R coordinate of Ohmic heating coils
  real*8,dimension(necoil):: re
  !> Z coordinate of Ohmic heating coils
  real*8,dimension(necoil):: ze
  !> width of Ohmic heating coils
  real*8,dimension(necoil):: we
  !> height of Ohmic heating coils
  real*8,dimension(necoil):: he
  real*8,dimension(necoil):: ecid
  real*8,dimension(necoil):: ecturn

end module cecoil


!***********************************************************************
! module for coilsp                                                    *
!***********************************************************************
module coilsp
  implicit none
  public
  real*8,dimension(300) :: rsplt,zsplt,csplt

end module coilsp


!***********************************************************************
! module for cvesel                                                    *
!***********************************************************************
module cvesel
  use exparm,only:nvesel
  implicit none
  public

  real*8,dimension(nvesel) :: rvs,zvs,wvs,hvs,avs,avs2,&
                            rsisvs,vsid

end module cvesel


!***********************************************************************
! module for var_filech                                                *
!***********************************************************************
module var_filech
        
  character*4 :: ch1, ch2
        
end module var_filech


!***********************************************************************
! module for fshift                                                    *
!***********************************************************************
module fshift
  use exparm,only:nfcoil,magpr2
  implicit none
  public

  integer*4,dimension(nfcoil) :: nshiftrz
  real*8,dimension(nfcoil) :: rshift,zshift,pshift
  real*8,dimension(magpr2) :: pmprobe

end module fshift


!***********************************************************************
! module for nio                                                       *
!***********************************************************************
module nio
  implicit none
  public

  !> file unit of input file
  integer*4,parameter :: nin = 11
  !> file unit of output file
  integer*4,parameter :: nout = 10
  !> file unit for command window
  integer*4,parameter :: ntty = 5
  integer*4,parameter :: nrsppc = 25
  integer*4,parameter :: nrspfc = 26
  integer*4,parameter :: ncontr = 35

end module nio


!***********************************************************************
! module for pmodel                                                    *
!***********************************************************************
module pmodel
  implicit none
  public

  real*8,dimension(:),allocatable,save :: rgrid
  real*8,dimension(:),allocatable,save :: zgrid
  real*8 :: dr,dz

end module pmodel


!***********************************************************************
! module for magnetic probe diagnostics                                *
!***********************************************************************
module mprobe
  use exparm,only:magpr2
  implicit none
  public

  !> R coordinate of magnetic probe in meter
  real*8,dimension(magpr2) :: xmp2
  !> Z coordinate of magnetic probe in meter
  real*8,dimension(magpr2) :: ymp2
  !> angle of magnetic probe in degree
  real*8,dimension(magpr2) :: amp2
  !> length of magnetic probe in meter
  real*8,dimension(magpr2) :: smp2
  integer*4 :: nsmp2

end module mprobe


!***********************************************************************
! module for flux loop diagnostics                                     *
!***********************************************************************
module siloop
  use exparm,only:nsilop
  implicit none
  public

  !> R coordinate of flux loop diagnostics
  real*8,dimension(nsilop) :: rsi
  !> Z coordinate of flux loop diagnostics
  real*8,dimension(nsilop) :: zsi
  !> width of flux loop diagnostics
  real*8,dimension(nsilop) :: wsi
  !> height of flux loop diagnostics
  real*8,dimension(nsilop) :: hsi
  real*8,dimension(nsilop) :: as
  real*8,dimension(nsilop) :: as2

end module siloop


!***********************************************************************
! module for Rogowski coil diagnostics                                 *
!***********************************************************************
module rogowl
  use exparm,only:nrogow
  implicit none
  public

  integer*4,dimension(nrogow) :: narc
  real*8,dimension(nrogow) :: prname
  real*8,dimension(36) :: rp,zp
  real*8,dimension(101) :: rpg,zpg

end module rogowl

