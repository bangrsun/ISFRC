      module exparm
! module for experimental parameters
! Revisions:
! $Log: exparm.f90,v $
! Revision 1.1.2.1  2008/11/18 22:24:35  radhakri
! *** empty log message ***
!
! Revision 1.1  2007/06/01 05:16:13  renq
! *** empty log message ***
!nesum changed from 18 to 6 by srini to match with public domain ver 09/25/2008
!
!     magpr2    total number of magnetic probes (magpri in EFIT)           
!     nacoil    number of advance divertor coils                          
!     necoil    number of ohmic heating coils                             
!     nesum     number of o.h. coil groups                                
!     nfcoil    number of p.f. coils                                      
!     nfsum     number of p.f. coil groups                                
!     nvsum     number of vessel segement groups                          
!     nrogow	  number of partial rogowski loops                            
!     nsilop    number of flux loops                                      
!     nvesel    number of vessel segements     

      implicit none
      public

      integer*4,parameter :: nfcoil = 6,nsilop = 41,magpr2 = 63,&
                           nrogow = 1,necoil = 1,nesum = 1,&
                           nfsum = 6,nvsum = 24,nvesel = 24,&
                           nacoil = 1,mgaus1 = 8,mgaus2 = 10
      integer*4 :: nw,nh,nwnh
      end module exparm
