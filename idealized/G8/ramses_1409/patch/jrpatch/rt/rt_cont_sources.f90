!************************************************************************
SUBROUTINE add_cont_sources(ilevel,dt)

! 'Interface' for injecting any continuous sources (RT and non-RT) into 
! the grid.
! This routine can be overridden by any patch to include continuous 
! sources.
!------------------------------------------------------------------------
  use amr_parameters,only:dp
  implicit none
  integer::ilevel
  real(dp)::dt
!------------------------------------------------------------------------

END SUBROUTINE add_cont_sources




