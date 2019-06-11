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
  use domain
  use vesel
  use fcoil
  implicit none
  logical file_exist
  integer i,j
  real*8 dr,dz

  open(unit=fu_gfunc,status='unknown',file='gfunc.dat')
  if(nvesel > 0) then
    write(fu_gfunc,*) gfvesel_rzv
  endif
  if(nfcoil > 0) then
    write(fu_gfunc,*) gffcoil_rzf
  endif

  flush(fu_gfunc)
  close(fu_gfunc)
  return
end subroutine writegfunc

