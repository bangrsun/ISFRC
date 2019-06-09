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
!      writegfunc writes the Green's function into output files.       *
!                                                                      *
! calling arguments:                                                   *
!                                                                      *
!***********************************************************************
subroutine writegfunc
  use exparam
  use nio
  use fcoil
  use vesel
  !use ecoil
  !use acoil
  !use fldiag
  !use mpdiag
  !use rogdiag
  implicit none
  logical file_exist
  integer i,j
  real*8 dr,dz

  open(unit=fu_gfunc,status='unknown',file='gfunc.dat')
  if(nfcoil > 0) then
    write(fu_gfunc,*) gffcoil_rzf
  endif
  if(nvesel > 0) then
    write(fu_gfunc,*) gfvesel_rzv
  endif

  flush(fu_gfunc)
  close(fu_gfunc)
  return
end subroutine writegfunc

