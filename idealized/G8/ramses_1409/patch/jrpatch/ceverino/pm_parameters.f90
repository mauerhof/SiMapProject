module pm_parameters
  use amr_parameters, ONLY: dp
  integer::npartmax=0               ! Maximum number of particles
  integer::nsinkmax=10000           ! Maximum number of sinks
  integer::npart=0                  ! Actual number of particles
  integer::nsink=0                  ! Actual number of sinks
  integer::iseed=0                  ! Seed for stochastic star formation
  integer::nstar_tot=0              ! Total number of star particle
  integer::ir_cloud=4               ! Radius of cloud region in unit of grid spacing
  real(dp)::mstar_tot=0             ! Total star mass
  real(dp)::mstar_lost=0            ! Missing star mass
  real(dp)::T2thres_SF=1d4          ! Temperature threshold for FS            !cev patch
  real(dp)::SN_delay_Myr=5d0        ! Time delay for SN feedback              !cev patch
end module pm_parameters
