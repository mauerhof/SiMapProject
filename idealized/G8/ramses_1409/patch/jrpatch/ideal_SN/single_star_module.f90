MODULE single_star_module

! Obtains values for stellar winds based on:
! Stellar mass (Msolar)
! Stellar age (Myears)
! Stellar metallicty (log([Fe/H])
! Outputs:
! Wind speed (cm/s)
! Mass lost in the given timestep (Msolar)
! Sam Geen, June 2010

! Use lookup table module for speed, mass-loss arrays
use amr_parameters
use lookup_table_module

implicit none

public

type(lookup_table) :: speed_table
type(lookup_table) :: massloss_table

CONTAINS

!*************************************************************************
! Sets up the tables, and then clears them (e.g. on program exit)
SUBROUTINE setup_single_star(speedfile,massfile)
  use amr_commons,only:myid
  character(len=128)::speedfile,massfile
  logical::ok
!-------------------------------------------------------------------------
  inquire(file=TRIM(speedfile), exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then 
        write(*,*)'Cannot fine starwind_speed file...'
        write(*,*)'File '//TRIM(speedfile)//' not found'
        write(*,*)'You need to set the SN_speedfile ' //   &
                  ' parameter to the correct path'
     endif
     call clean_stop
  end if

  inquire(file=TRIM(massfile), exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then 
        write(*,*)'Cannot fine starwind_massloss file...'
        write(*,*)'File '//TRIM(massfile)//' not found'
        write(*,*)'You need to set the SN_massfile to the correct path'
     endif
     call clean_stop
  end if

  call setup_table(speed_table, speedfile)
  call setup_table(massloss_table, massfile)
  !call setup_table(speed_table, 'starwind_speedtable.dat')
  !call setup_table(massloss_table, 'starwind_masslosstable.dat')
END SUBROUTINE setup_single_star

!*************************************************************************
SUBROUTINE cleanup_single_star()
!-------------------------------------------------------------------------
  call clear_table(speed_table)
  call clear_table(massloss_table)
END SUBROUTINE cleanup_single_star

!************************************************************************
! Finds the wind speed and massloss for a given star (mass, age, metallicity)
! Note: mass is *initial* mass, *not* current mass
! Values from Padova model (e.g. Girardi et al 2008)
SUBROUTINE star_wind_values(mass, age, metal, speedout, masslossout)
  real(dp)::mass, age, metal,logage
  real(dp),intent(out)::speedout, masslossout
!-------------------------------------------------------------------------
  ! Lookup table contains time data based on log(years)
  ! So we need to convert that shit right up
  ! Convert Myr --> log(yr)
  logage = DLOG10(age)+6d0
  ! Find values in the tables
  call find_value3(speed_table, mass, logage, metal, speedout)
  call find_value3(massloss_table, mass, logage, metal, masslossout)
  ! Uh, yeah. That's basically it.

END SUBROUTINE star_wind_values

!************************************************************************
! How long before a star croaks it, and what is its final mass and energy?
! Note: "mass" is *initial* mass, *not* current mass
! Mass in Msolar, time in Myr, energy in 10^51 ergs
! Values from Smartt et al 2009, Kovetz et al 2009
SUBROUTINE star_snparams(mass, metal, fmassout, lifetimeout, snenergyout)
!-------------------------------------------------------------------------
  ! HACK! HACK! HORRIBLE HACK! ONLY VALUES FOR 15Msolar
  real(dp)::mass, metal
  real(dp),intent(out)::fmassout, lifetimeout, snenergyout
  ! HACK! HARD-CODED VALUES FOR STAR OF 15 Msolar
  fmassout = 1.5d0
  lifetimeout = 14.125d0
  snenergyout = 1.2d0 ! This is a total guess
END SUBROUTINE star_snparams

END MODULE
