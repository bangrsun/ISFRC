      module rogowl
!module for rogowl
!Revisions:
!$Log: rogowl.f90,v $
!Revision 1.1  2007/06/01 05:16:35  renq
!*** empty log message ***
!
!
      use exparm,only:nrogow
      implicit none
      public

      integer*4,dimension(nrogow) :: narc
      real*8,dimension(nrogow) :: prname
      real*8,dimension(36) :: rp,zp
      real*8,dimension(101) :: rpg,zpg

      end module rogowl
