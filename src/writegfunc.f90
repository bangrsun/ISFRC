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
  use nio,only:fu_gfunc
  use global_params,only:gffcoil_rzf
  implicit none

  open(unit=fu_gfunc,status='unknown',file='gfunc.dat')
  write(fu_gfunc,*) gffcoil_rzf

  flush(fu_gfunc)
  close(fu_gfunc)
  return
end subroutine writegfunc

