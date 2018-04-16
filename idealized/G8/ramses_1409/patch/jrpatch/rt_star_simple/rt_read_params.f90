!RT_STAR_SIMPLE patch (SS):
! Added parameters that specify how stellar intensity is divided between
! photon packages.
!************************************************************************

!************************************************************************
SUBROUTINE read_simple_star_params(nml_ok)

! Read simple_star_params namelist
!------------------------------------------------------------------------
  !use amr_commons
  !use rt_parameters
  use SED_module
  implicit none
  logical::nml_ok
!------------------------------------------------------------------------
  namelist/SIMPLE_STAR_PARAMS/ Ipac_fractions
  ! Read namelist file
  rewind(1)
  read(1,NML=SIMPLE_STAR_PARAMS,END=101)
101 continue                                  ! No harm if no rt namelist

  if(sum(Ipac_fractions) .ne. 1.) then
     write(*,*) 'WARNING: The photon Ipac_fractions do not sum up to 1'
  endif
END SUBROUTINE read_simple_star_params

!************************************************************************
SUBROUTINE rt_read_params(nml_ok)

! 'Interface' for reading any additional parameters from .nml file.
! This routine can be overridden by any patch to include more parameters.
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  logical::nml_ok
!------------------------------------------------------------------------
  call read_simple_star_params(nml_ok)
  call read_rt_params(nml_ok)
END SUBROUTINE rt_read_params

