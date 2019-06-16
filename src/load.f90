!***********************************************************************
! This file is part of project ISFRC.                                  *
! ISFRC is an Integrated Simulation code for Field-Reversed            *
! Configuration.                                                       *
!                                                                      *
! Developers:                                                          *
! Shuying SUN, Yang LI, Huasheng XIE, ...                              *
!                                                                      *
! ENN Sci. & Tech. Development Coorporation, 2008-2019.                *
! (c) All rights reserved.                                             *
!***********************************************************************

!***********************************************************************
! subprogram description:                                              *
!      calc_gf calculate the Green's function of all active coils,     *
!      including the fcoil, vessel, acoil, ecoil. This subroutine      *
!      does not calculate the Green's function from plasma itself.     *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine load
  use consta,only:tmu
  use readin_params
  use global_params
  implicit none
  real*8,external :: mutpsi
  real*8 :: r1,z1,jsplit,gsplit
  real*8,dimension(:,:),allocatable :: rsplit,zsplit
  integer i,j,jf,isr,isz

!----------------------------------------------------------------------
!-- calculate psi by fcoils                                          --
!----------------------------------------------------------------------
  psif_rz=0.0d0
  if(nfcoil > 0) then
    do j=1,ngz1
      do i=1,ngr1
        r1=rgrid_rz(i,j)
        z1=zgrid_rz(i,j)
        do jf=1,nfcoil
          allocate(rsplit(nsr_f(jf),nsz_f(jf)),zsplit(nsr_f(jf),nsz_f(jf)))
          call splitcoil(r_f(jf),z_f(jf),w_f(jf),h_f(jf),ar_f(jf),az_f(jf), &
            nsr_f(jf),nsz_f(jf),rsplit,zsplit)
          Jsplit=J_f(jf)/nsr_f(jf)/nsz_f(jf)
          do isz=1,nsz_f(jf)
            do isr=1,nsr_f(jf)
              gsplit=mutpsi(rsplit(isr,isz),r1,z1-zsplit(isr,isz))
              psif_rz(i,j)=psif_rz(i,j)+tmu*Jsplit*gsplit
            enddo
          enddo
          deallocate(rsplit,zsplit)
        enddo
      enddo
    enddo
  endif
!----------------------------------------------------------------------
!-- load initial psi at grid points                                  --
!----------------------------------------------------------------------
  psi_rz(:,:)=psif_rz(:,:)
 
  return
end subroutine load

