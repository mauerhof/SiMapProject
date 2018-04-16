!************************************************************************
SUBROUTINE rt_read_params(nml_ok)

! 'Interface' for reading any additional parameters from .nml file.
! This routine can be overridden by any patch to include more parameters.
!------------------------------------------------------------------------
  use ideal_SN_module !-----------------------------------------!ideal_SN
  implicit none
  logical::nml_ok
!------------------------------------------------------------------------
  call read_rt_params(nml_ok)
  call read_ideal_SN_params(nml_ok) !---------------------------!ideal_SN

END SUBROUTINE rt_read_params




