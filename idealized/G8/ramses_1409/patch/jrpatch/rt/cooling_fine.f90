!RT patch: A lot of change here, since we're not only evolving the 
!temperature but also ionization states and photon densities and fluxes.
!To see the RT changes, do a diff.
!*************************************************************************
SUBROUTINE cooling_fine(ilevel, dt)

! Compute ionization, recombination and cooling for leaf cells
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use cooling_module
  use rt_parameters, ONLY: rt_freeflow
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  real(dp)::dt  !---------------------------------------------------------
  integer:: ncache,i,igrid,ngrid,info
  integer,dimension(1:nvector),save:: ind_grid
!-------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Compute sink accretion rates
  if(sink)call compute_accretion_rate(0)

  if(rt_freeflow) return            ! ATT: use rt_smooth=.false. with this

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid, ngrid, ilevel, dt)
  end do

  if(ilevel==levelmin.and.cosmo)then
     !if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if
111 format('   Entering cooling_fine for level',i2)
END SUBROUTINE cooling_fine

!*************************************************************************
SUBROUTINE coolfine1(ind_grid, ngrid, ilevel, dt)

! Vector sweep for ionization and cooling of grids
! ind_grid => Indexes of grids/octs to cool in
! ngrid    => Number of valid indexes in ind_grid (i.e. number of grids)
! ilevel   => Level at which we're cooling
! dt       => timestep size in code units
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
  use rt_parameters
  use cooling_module
  implicit none
  integer::ngrid,ilevel
  integer,dimension(1:nvector)::ind_grid
  real(dp)::dt !----------------------------------------------------------
  integer::i, ind, iskip, idim, ncell, ii, ip
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(dp)::dtcool,nISM,nCOM,damp_factor,t_blast
  real(dp)::polytropic_constant
  integer, dimension(1:nvector),save::ind_cell
  real(dp),dimension(1:nvector),save::nH,T2,ekk
  real(dp),dimension(1:nvector),save::T2min,Zsolar
  logical,dimension(1:nvector),save::cooling_switch=.true. !-----------!RT
  real(dp)::scale_Np, scale_Fp !---------------------------------------!RT
  real(dp)::Fpnew !----------------------------------------------------!RT
  real(dp),dimension(1:nvector, n_U),save::U !-------------------------!RT
  real(dp),dimension(1:nvector,nPacs),save::Fp, Fp_precool !-----------!RT
  real(dp),dimension(1:nvector,nPacs),save::dNpdt=0., dFpdt=0. !-------!RT
!-------------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo) nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)  &
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax) & 
            *scale_l/aexp)**2/(twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2

  ! Loop over cells in the octs
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     ind_cell(1:ngrid)=iskip+ind_grid(1:ngrid)

     ! Gather leaf cells
     ncell=0
     do i=1,ngrid
        if(son(ind_cell(i))==0) then
           ncell = ncell + 1
           ind_cell(ncell) = ind_cell(i)
        end if
     end do
     if(ncell .eq. 0) cycle                        !    No leaf cells here

     ! Get gas density (mass density in code units)-----------------------
     do i=1,ncell
        nH(i) = MAX(uold(ind_cell(i),1),smallr)
     enddo
     ! Get the ionization fractions---------------------------------------
     do i=0,nIons-1
        U(1:ncell,2+i) = uold(ind_cell(1:ncell),iIons+i)/nH(1:ncell)
     end do
     ! Get metallicity----------------------------------------------------
     if(metal) then
        Zsolar(1:ncell)= uold(ind_cell(1:ncell),imetal) /nH(1:ncell) /0.02
     else                                          !    Zsolar is fraction 
           Zsolar(1:ncell) = z_ave                 !  of solar metallicity
     endif
     ! Get T/mu-----------------------------------------------------------
     T2(1:ncell) = uold(ind_cell(1:ncell),ndim+2)  !        e_kin + e_heat
     ekk(1:ncell) = 0.0d0                          !
     do idim=1,ndim                                !     Calculate kinetic
        ekk(1:ncell) = ekk(1:ncell) &              !       energy density.
                      +0.5*uold(ind_cell(1:ncell),idim+1)**2/nH(1:ncell)
     end do
     T2(1:ncell) = (gamma-1.0)*(T2(1:ncell)-ekk(1:ncell))
     ! now T2 is pressure (code units)...convert that to T/mu in Kelvin
     T2(1:ncell) = T2(1:ncell)/nH(1:ncell)*scale_T2 !       ideal gas eq.

     ! Compute nH in H/cc (number of H nuclei per cubic centimeter)
     nH(1:ncell) = nH(1:ncell)*scale_nH

     !====================================================================
     ! Compute temperature from polytrope EOS
     !====================================================================
     if(jeans_ncells>0)then
        T2min(1:ncell) = nH(1:ncell)*polytropic_constant*scale_T2
     else
        T2min(1:ncell) = T2_star*(nH(1:ncell)/nISM)**(g_star-1.0)
     endif
     !====================================================================
     ! You can put your own polytrope EOS here
     !====================================================================
     ! Compute "thermal" temperature by subtracting polytrope
     do i=1,ncell
        U(i,1) = MAX(T2(i)-T2min(i),T2_min_fix) 
     enddo

     ! Delayed cooling: Turn off cooling in blast wave regions
     if(cooling .and. delayed_cooling) then
        cooling_switch(1:ncell)=.true.
        do i=1,ncell
           if(uold(ind_cell(i),idelay)/uold(ind_cell(i),1) .gt. 1d-3)    &
                cooling_switch(i)=.false.
        end do
     endif
     if(isothermal) cooling_switch(1:ncell)=.false.

     ! Compute cooling time step in seconds
     dtcool = dt*scale_t

     ! Get photon densities and flux magnitudes---------------------------
     do ip=1,nPacs
        do i=1,ncell
           U(i,iNpU(ip)) = scale_Np * uold(ind_cell(i),iPac(ip))
           U(i,iFpU(ip)) = scale_Fp *                                    &
                sqrt(sum((uold(ind_cell(i),iPac(ip)+1:iPac(ip)+ndim))**2))
        enddo
        if(rt_smooth) then                              ! Smooth RT update 
           do i=1,ncell   ! Calc addition per sec to Np, Fp for current dt
              dNpdt(i,ip) = scale_Np / dtcool *                          &
                   (unew(ind_cell(i),iPac(ip))-uold(ind_cell(i),iPac(ip)))
              Fpnew = scale_Fp * sqrt(                                   &
                   sum((unew(ind_cell(i),iPac(ip)+1:iPac(ip)+ndim))**2))
              dFpdt(i,ip) = (Fpnew-U(i,iFpU(ip)))/dtcool!Change in magnit.
              ! Update flux vector to get the right direction
              uold(ind_cell(i),iPac(ip)+1:iPac(ip)+ndim) =               &
                                unew(ind_cell(i),iPac(ip)+1:iPac(ip)+ndim)
              Fp_precool(i,ip)=Fpnew      ! For update after solve_cooling
           enddo
        else
           Fp_precool(:,ip)=U(:,iFpU(ip)) ! For update after solve_cooling
        endif
     end do
     ! Now uold contains the post-hydro Fp vectors and Fp_precool contains
     ! their magnitudes, for update after solve_cooling.------------------

     !====================================================================
     call solve_cooling(U, dNpdt,  dFpdt,  nH,   cooling_switch, Zsolar  &
                         , dtcool, aexp,   ncell)
     !====================================================================

     nH(1:ncell) = nH(1:ncell)/scale_nH   ! Convert rho back to code units

     ! a) Update energy density-------------------------------------------
     uold(ind_cell(1:ncell),ndim+2) = (U(1:ncell,1)+T2min(1:ncell))      &
                     * nH(1:ncell) / scale_T2 / (gamma-1.0) + ekk(1:ncell)

     do ii=0,nIons-1 ! b) Update ionization fractions---------------------
        uold(ind_cell(1:ncell),iIons+ii) = U(1:ncell,2+ii)*nH(1:ncell)
     end do

     ! c) ...and update photon densities and flux magnitudes--------------
     do ip=1,nPacs
        do i=1,ncell
           uold(ind_cell(i),iPac(ip)) = U(i,iNpU(ip)) / scale_Np        !&
                !max(U(i,iNpU(ip)) / scale_Np,  smallNp)
           if(Fp_precool(i,ip) .gt. 0.d0)                                &
                uold(ind_cell(i),iPac(ip)+1:iPac(ip)+ndim) =             &
                   U(i,iFpU(ip))/Fp_precool(i,ip)                        &
                   *uold(ind_cell(i),iPac(ip)+1:iPac(ip)+ndim)
        enddo
     end do

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=20d0*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,ncell
           uold(ind_cell(i),idelay)=uold(ind_cell(i),idelay)*damp_factor
        end do
     endif

  end do
  ! End loop over cells
END SUBROUTINE coolfine1



