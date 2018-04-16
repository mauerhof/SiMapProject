!************************************************************************
subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

! Conversion factors from user units into cgs units
! For gravity runs, make sure that G=1 in user units.
!
! scale_l   <=  length scale [cm]
! scale_t   <=  time scale [t]
! scale_d   <=  mass density scale [g cm-3]
! scale_v   <=  speed scale [cm s-1]
! scale_nH  <=  number density scale [cm-3]
! scale_T2  <=  Temperature scale [Kelvin]
!------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
!------------------------------------------------------------------------
  ! scale_d converts mass density from user units into g/cc
  scale_d = units_density
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  scale_t = units_time
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_l converts distance from user units into cm
  ! since the box always has length 1 in UU, scale_l is the 
  ! length of the box in cgs
  scale_l = units_length
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin,
  ! using ideal gas equation
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units

!BEGIN !RT---------------------------------------------------------------
!************************************************************************
subroutine rt_units(scale_Np, scale_Fp)

! Conversion factors from user units into cgs units for RT variables.
!
! scale_np   <=  photon number density scale [# cm-3]
! scale_pf   <=  photon flux density scale [# cm-2 s-1]
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  implicit none
  real(dp)::scale_Np, scale_Fp
  real(dp)::scale_nH, scale_T2, scale_t, scale_v, scale_d, scale_l
!------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  scale_Np = 1.0
  if(cosmo) scale_Np = 1./aexp**3         ! or maybe it should be aexp**4

  scale_Fp = scale_Np  * scale_v

end subroutine rt_units
!END-- !RT----------------------------------------------------------------
