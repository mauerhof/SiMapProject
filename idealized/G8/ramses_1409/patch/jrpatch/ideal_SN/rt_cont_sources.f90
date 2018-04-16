MODULE ideal_SN_module
  use amr_parameters,only:dp,ndim,MAXREGION
  use rt_parameters,only:rt
  use single_star_module
  implicit none

  public

  ! SN namelist parameters------------------------------------------------
  integer                         ::n_SN_source=0
  real(dp),dimension(1:MAXREGION) ::SN_posx=0.
  real(dp),dimension(1:MAXREGION) ::SN_posy=0.
  real(dp),dimension(1:MAXREGION) ::SN_posz=0.
  real(dp),dimension(1:MAXREGION) ::SN_stellar_mass_solar=0.
  logical,dimension(1:MAXREGION)  ::SN_wind=.false.
  real(dp),dimension(1:MAXREGION) ::SN_metallicity=-2.  ! Log(metallicity)
  character(len=128)              ::SN_tabledir=''
  integer                         ::SN_r_wind      ! Wind rad [fine cells]
  logical::source_odd=.true.   ! Odd means the source is at a cell center,
                               ! even=intersection of cells
  real(dp)::w_speed_cgs=0., w_mdot_cgs=0., w_Tk=0.
  ! ----------------------------------------------------------------------
  ! Other variables
  real(dp)::w_dm_cgs=0., w_dp_cgs=0., w_dEk_cgs=0., w_dET_cgs=0.
  logical,dimension(1:MAXREGION)  ::SN_done=.false.
  real(dp),dimension(1:MAXREGION,1:ndim)::SN_pos=0.
  real(dp),dimension(1:MAXREGION) ::SN_E_cgs=0.
  real(dp),dimension(1:MAXREGION) ::SN_time_myr=-1.
  real(dp),dimension(1:MAXREGION) ::SN_time_uu=-1.
  real(dp),dimension(1:MAXREGION) ::SN_mremnant_solar=0.
  real(dp)::dx_levelmax, r_star
  logical::first_add_source=.true.
  logical::is_feedbk_file=.false.
  real(dp),parameter::courant_wind = 0.5     ! Max % of mass to add per dt
  real(dp),parameter::m_sun=1.9891d33        ! Solar mass in grams
  real(dp),parameter:: Myr = 3.15569d+13
  real(dp),parameter ::mH    = 1.6600000d-24    !          H atom mass [g]
  real(dp),parameter ::kB    = 1.38062d-16      !  Boltzm.const. [erg K-1]
  real(dp)::mlost_sol_tbl=0.,mlost_sol_rl=0.    !          For bookkeeping
  logical::SN_GO=.false.     ! If true then inject SN at first opportunity
  logical::SN_wait=.false.   !       Keep zero timestep while this is true

  real(dp),allocatable,dimension(:,:,:)::Wgrid,Wgrid_all ! SN wind weights
  ! Wgrid represents a sphere around the source of radius SN_r_wind
  ! ----------------------------------------------------------------------

CONTAINS

!*************************************************************************
SUBROUTINE read_ideal_SN_params(nml_ok)

! Read ideal_SN_params namelist
!-------------------------------------------------------------------------
  use amr_commons
  use sed_module
  implicit none
  logical::nml_ok
  integer::i
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::fmassout
!-------------------------------------------------------------------------
  namelist/IDEAL_SN_PARAMS/ n_SN_source, SN_posx, SN_posy, SN_posz       &
       & ,SN_stellar_mass_solar, SN_wind, SN_metallicity, SN_tabledir    &
       & ,SN_r_wind, source_odd, w_speed_cgs, w_mdot_cgs, w_TK
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rewind(1)
  read(1,NML=IDEAL_SN_PARAMS,END=101)
  do i=1,n_SN_source
     SN_pos(i,1)=SN_posx(i)
     if(ndim .gt. 1) SN_pos(i,2)=SN_posy(i)
     if(ndim .gt. 2) SN_pos(i,3)=SN_posz(i)
  end do
  dx_levelmax=boxlen*0.5D0**nlevelmax
  if(SN_wind(1)) call prepare_Wind_Weights(1000000)
101 continue
  if(myid==1) write(*,*) 'Initializing SED table for stellar source'
  if(rt) call init_SED_table()
  call update_SN_PacProps()
  if(SN_tabledir .ne. '') is_feedbk_file=.true. 
  call setup_single_star(SN_tabledir)
  do i=1,n_SN_source
     call star_snparams(SN_stellar_mass_solar(i), SN_metallicity(i)      &
          ,SN_mremnant_solar(i), SN_time_myr(i),  SN_E_cgs(i))
  enddo
  SN_E_cgs(1:n_SN_source)=SN_E_cgs(1:n_SN_source)*1.d51
  SN_time_uu(1:n_SN_source) = &
       SN_time_myr(1:n_SN_source)*3.1556926d13/scale_t

  if(myid==1) then
     write(*,*),'Ideal supernova properties:-----------------------------'
     write(*,'(a, 20(1pe16.6))')                                         &
          '  Stellar mass [Solar] = ',SN_stellar_mass_solar(1:n_SN_source) 
     write(*,'(a, 20(1pe16.6))')                                         &
          '  SN time [Myr]        = ',SN_time_myr(1:n_SN_source) 
     write(*,'(a, 20(1pe16.6))')                                         &
          '  SN E [erg]           = ',SN_e_cgs(1:n_SN_source) 
     write(*,'(a, 20(1pe16.6))')                                         &
          '  Remnant mass [Solar] = ',SN_mremnant_solar(1:n_SN_source) 
     write(*,'(a, 20(1pe16.6))')                                         &
          '  SN Z [log(Fe/H)]     = ',SN_metallicity(1:n_SN_source)
     if(.not. is_feedbk_file) then
        write(*,'(a, 1pe16.6)') '  Wind speed [cm/s] = ',w_speed_cgs
        write(*,'(a, 1pe16.6)') '  Mass loss [g s]   = ',w_mdot_cgs
        write(*,'(a, 1pe16.6)') '  Wind T [K]        = ',w_TK
     endif
     write(*,*),'--------------------------------------------------------'
  endif

END SUBROUTINE read_ideal_SN_params

!*************************************************************************
SUBROUTINE SN_courant_fine(ilevel,dt)

! Set a minimum timestep length, as recommended by this module.
! We don't want to overstep a SN event, so the timestep recommended is
! up until the next supernova event.
! Also we constrain the timestep so that the wind doesn't traverse more
! than one cell (this is to get a small timestep at t=0 as the initial
! conditions in an ideal_SN experiment are static conditions leading to 
! a very large default timestep).
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_parameters
  use rt_parameters
  implicit none
  real(dp)::dt,SN_dtmin,dx,dx_loc,scale
  integer::ilevel,i,nx_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::dt_ok=.false.
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  if(ilevel.lt.nlevelmax) return       ! Only constrain dt if finest level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  SN_dtmin=dt
 
  if(any(SN_wind(1:n_SN_source))) then     ! Wind mass injection condition
     ! The mass condition must be checked before the velocity condition  !
     ! because w_speed_cgs and w_TK are set in stellar_feedback, which   !
     ! is called in the mass condition check.                            !
     i=0
     do                   ! Perform fake feedback to see if timestep is ok
        call stellar_feedback(ilevel, t+SN_dtmin, SN_dtmin, dt_ok, .true.)
        i=i+1
        if(i.gt.1) print*,'dt loop ',myid,i,SN_dtmin
        if(dt_ok) exit                                 !    Timestep is ok
        SN_dtmin = SN_dtmin/2.                         !  Timestep too big
     enddo
  endif
  do i=1,n_SN_source
     if(.not. SN_done(i) .and. SN_time_uu(i)-t.ge.0.) then
        SN_dtmin=min(SN_dtmin, SN_time_uu(i)-t)        ! Don't overstep SN
     endif
     if(SN_wind(i)) then
        if(nstep_coarse .lt. 10) then
           ! Velocity Courant:
           SN_dtmin=min(SN_dtmin,dx_loc/(w_speed_cgs/scale_v))
           ! Sound speed Courant:
           SN_dtmin=min(SN_dtmin                                        &
                          ,dx_loc / (sqrt(gamma*kb*w_TK/mH)/scale_v*100.))
        endif
     endif
  end do
  if(doDump .or. SN_GO .or. SN_wait) SN_dtmin=0.d0  ! Superova event!
  SN_dtmin=max(0., SN_dtmin)
  dt=min(dt, SN_dtmin)

END SUBROUTINE SN_courant_fine

!*************************************************************************
SUBROUTINE ideal_SN_refine(xx,ok,ncell,ilevel)

! This routine flags cells immediately around SN sources to the finest
! level of refinement. The criteria for refinement at a point are:
! a) The point is less than one ilevel cell width from an SN source.
! b) The point is within SN_r_wind finest level cell widths from
!    the SN source.
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell,ilevel,i,k,nx_loc
  real(dp),dimension(1:nvector,1:ndim)::xx
  logical ,dimension(1:nvector)::ok
  real(dp)::dx_loc,rvec(ndim),w,rmag,rSN
!-------------------------------------------------------------------------
  nx_loc=(icoarse_max-icoarse_min+1)
  dx_loc = boxlen*0.5D0**ilevel/dble(nx_loc)
  rSN = SN_r_wind*dx_levelmax
  ! Loop over regions
  do k=1,n_SN_source
     do i=1,ncell
        rvec(:)=xx(i,:)-SN_pos(k,:)
        rmag=sqrt(sum(rvec**2))
        if(rmag .le. 2*rSN+dx_loc) then
           ok(i)=.true.
        endif
     end do
  end do
END SUBROUTINE ideal_SN_refine

!*************************************************************************
SUBROUTINE update_SN_PacProps()

! Compute and assign to all SED photon packages the quantities pac_csn, 
! pac_cse and pac_egy corresponding to the age and metallicity of the 
! first SN source.
!-------------------------------------------------------------------------
  use amr_commons,only:myid
  use rt_parameters
  use sed_module
  integer :: ip, ii
  real(dp):: egy(nSEDPacs), age, csn(nSEDPacs,nIons), cse(nSEDPacs,nIons)
!--------------------------------------------------------------------------
  if(.not. rt) return
  call getAgeGyr(0.d0, age)                                !  Age = [Gyrs]
  age = log10(max(age, SED_minAge))                        ! In case age=0
  call inp_SED_table(age, SN_metallicity(1), 2, egy)       !          [eV]
  do ii=1,nIons
     call inp_SED_table(age,SN_metallicity(1),1+2*ii,csn(:,ii))   ! [cm^2]
     call inp_SED_table(age,SN_metallicity(1),2+2*ii,cse(:,ii))   ! [cm^2]
  end do
  do ip=1,nSEDpacs
     pac_egy(ip)   = egy(ip)
     pac_csn(ip,:) = csn(ip,:)
     pac_cse(ip,:) = cse(ip,:)
  end do
  call updateRTPac_CoolConstants
  if(myid==1) write(*,*) 'SED photon packages updated for SN source'
  call write_PacProps(.true.,6)
END SUBROUTINE update_SN_PacProps


!*************************************************************************
SUBROUTINE update_wind_props(tf, dt)

! Read the current stellar wind properties into w_s_cgs, w_mdot_cgs, w_TK
! and w_Sr, for a given timestep.
! tf => End time of timestep to be used [uu]
! dt => Timestep length [uu]
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_parameters
  real(dp)::tf, dt, ti_sec, tf_sec, wm_solar_i, wm_solar_f
  real(dp)::wp_cgs_i, wp_cgs_f, wEk_cgs_i, wEk_cgs_f, wET_cgs_i, wET_cgs_f
  real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp),parameter:: Myr = 3.15569d+13
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call getProperTime(tf-dt, ti_sec)  ! [uu]
  ti_sec = ti_sec*scale_t
  call getProperTime(tf,    tf_sec)  ! [uu]
  tf_sec = tf_sec*scale_t

  if(.not. is_feedbk_file) then      ! Just using constant wind properties
     w_dm_cgs  = w_mdot_cgs * (tf_sec-ti_sec)
     w_dp_cgs  = w_speed_cgs*w_dm_cgs
     w_dEk_cgs = 0.5d0*w_dp_cgs**2/w_dm_cgs
     w_dET_cgs = w_TK * kB/mH * 1.7/(gamma-1.0) * w_dm_cgs    ! 1/mu = 1.7
     return
  endif

  call star_wind_values(SN_stellar_mass_solar(1), ti_sec/Myr             &
                       ,10.d0**SN_metallicity(1)                         &
                       ,wm_solar_i, wp_cgs_i, wEk_cgs_i, wET_cgs_i)
  call star_wind_values(SN_stellar_mass_solar(1), tf_sec/Myr             &
                       ,10.d0**SN_metallicity(1)                         &
                       ,wm_solar_f, wp_cgs_f, wEk_cgs_f, wET_cgs_f)
  w_dm_cgs    = ( wm_solar_f-wm_solar_i ) * m_sun
  w_dp_cgs    = wp_cgs_f  - wp_cgs_i
  w_dEk_cgs   = wEk_cgs_f - wEk_cgs_i
  w_dET_cgs   = wET_cgs_f - wET_cgs_i

  w_mdot_cgs=0. ; w_speed_cgs=0. ; w_TK=0.
  if((tf_sec-ti_sec) .gt. 0.) w_mdot_cgs  = w_dm_cgs / (tf_sec-ti_sec)
  if(w_dm_cgs .gt. 0.) then
     w_speed_cgs = w_dp_cgs / w_dm_cgs
     w_TK = w_dET_cgs * mH/kB * (gamma-1.0)/1.7 / w_dm_cgs
    endif

  mlost_sol_tbl = wm_solar_f          !                    For bookkeeping

END SUBROUTINE update_wind_props

!*************************************************************************
SUBROUTINE stellar_feedback(ilevel, tf, dt, dt_ok, fakeIt)

! Inject SN regions (from the IDEAL_SN_PARAMS namelist). 
! AN SN region inputs ionizing radiation into a point and mass+velocity 
! around that point until the SN_time, at which a SN explosion is 
! triggered (at which a certain amount of energy is put into that cell).
!
! ilevel  => amr level at which to inject the radiation
! tf      => End time for the timestep to be used for feedback [uu]
! dt      => Timestep length [uu]. (The feedback then starts at t-dt)
! dt_ok  <=  Returns true if the used timestep doesn't break the 
!            recommended timestep constraint of not adding more than a
!            certain percent of mass to cells per timestep.
! fakeIt  => If true, don't actually do any feedback. This is useful
!            for determining the appropriate timestep for the feedback.
!------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
  integer::ilevel
  real(dp)::tf, dt
  logical::dt_ok, fakeIt, anySources
  integer::i,igrid,ncache,iskip,ngrid,ind,idim,ivar,ix,iy,iz,nx_loc
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::dx, dx_loc, scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
!------------------------------------------------------------------------
  dt_ok=.true.
  if(first_add_source) then         !              First time we're here:
     SN_done(1:n_SN_source)=.true.  !           Then we need to go though 
     do i=1,n_SN_source             !              all the SN sources and
        if(SN_time_uu(i).ge.0 .and. SN_time_uu(i).ge.t) & !   check which
             SN_done(i)=.false.     !              are already triggered,
     end do                         !                 in case of restart.
     first_add_source=.false.      
     if(any(SN_done(1:n_SN_source))) then   !      Turn off RT if SN done
        if(rt) rt_cool_only=.true.  
        rt=.false. 
     endif
  endif                             
  if(ilevel.eq.levelmin) call update_SN_PacProps()

  ! Renew the wind weights grid every coarse timestep
  if(ilevel.eq.levelmin .and. SN_wind(1))call prepare_Wind_Weights(100000)
  if(ilevel.lt.nlevelmax) return    !  Source should be refined to the max
  if(numbtot(1,ilevel)==0)return    !               No grids at this level

  if(SN_wind(1)) call update_wind_props(tf,dt)  !   Update wind properties
  if(.not. fakeIt) mlost_sol_rl=mlost_sol_rl+w_dm_cgs/m_sun !  Bookkeeping

  if(is_feedbk_file .and. .not.fakeIt .and. myid==1 .and. SN_wind(1)) then
     write(*,'(" Wind: t=",1pe9.3,  ", dm=",1pe9.3, " msun"              &
               ", mdot=",1pe9.3,    " g/s, s=",1pe9.3, " cm/s"           &
               ", T=",1pe9.3,       " K, Ekin/Ep=",1pe9.3                &
               ", dm_accu=",1pe9.3, " msun, dm_tbl/dm_real=",1pe9.3      &
               )') &
               tf,           w_dm_cgs/m_sun,                             &
               w_mdot_cgs,   w_dp_cgs/w_dm_cgs,                          &
               w_TK,         w_dEk_cgs/ (0.5*w_dp_cgs**2/w_dm_cgs),      &
               mlost_sol_rl, mlost_sol_tbl/mlost_sol_rl
     
  endif


  anysources=.false.
  if(SN_GO)                             anySources=.true.
  if(any(SN_wind(1:n_SN_source)))       anySources=.true.
  if(rt)                                anySources=.true.
  if(.not. anySources) return
  dx=0.5D0**ilevel       !  Mesh size at level ilevel in coarse cell units
  do ind=1,twotondim     ! Set pos of cell centers relative to grid center
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid
  ! dx (and dx_loc=dx) are just equal to 1/nx (where 1 is the boxlength)
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)             ! Gather cell indices
        end do
        ! Gather cell center positions and rescale to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do
        do ivar=1,nvar
           do i=1,ngrid
              uu(i,ivar)=uold(ind_cell(i),ivar)       ! Gather cell values
           end do
        end do
        ! find new values in injection cells
        call stellar_feedback_vsweep(xx,uu,dx_loc,tf,dt,dt_ok,ngrid      &
             ,fakeIt)
        ! ...and update those cells
        if( .not. fakeIt) then
           do ivar=1,nvar
              do i=1,ngrid
                 uold(ind_cell(i),ivar)=uu(i,ivar)
              end do
           end do
        endif
     end do ! End loop over cells
  end do ! End loop over grids
  if(.not. fakeIt) SN_GO=.false.

END SUBROUTINE stellar_feedback

!*************************************************************************
SUBROUTINE stellar_feedback_vsweep(x, uu, dx, tf, dt, dt_ok, nn, fakeIt)

! Do a vector sweep, injecting SN regions into cells, that is if
! the cells are in any of these regions.
!
! x      =>  ncells*ndim: positions of grid cells
! uu    <=   ncells*nvars: injected variables in each cell
! dx     =>  Cell width in code units
! tf     =>  End time of timestep [uu]
! dt     =>  Timestep length [uu]
! dt_ok <=   True if timestep ok, otherwise false
! nn     =>  int number of cells
!------------------------------------------------------------------------
  use amr_commons
  use hydro_parameters
  use rt_parameters
  use sed_module
  implicit none
  integer ::nn,i,j,k,ip
  real(dp)::dx,tf,dt,dt_rec,dx_cgs,dt_cgs
  logical::dt_ok,fakeIt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:ndim)::rCICvec,rvec
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(dp)::scale_np, scale_fp, vol_cgs, rCIC
  real(dp)::meject_solar, meject_cgs, Z
  real(dp)::sed_I(nSEDpacs), r, w, pmag, age, drho
#ifdef SOLVERhydro
  integer ::imetal=6
#endif
#ifdef SOLVERmhd
  integer ::imetal=9
#endif
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np,scale_fp)
  dx_cgs=dx*scale_l  ;  dt_cgs=dt*scale_t
  vol_cgs=dx_cgs**ndim

  do k=1,n_SN_source                                   ! Loop over regions
     do i=1,nn                                         !   Loop over cells
        ! Compute CIC weights relative to SN source
        if(source_odd) then 
           rCICvec(:)=max(1.0-abs(x(i,:)-SN_pos(k,:)+dx/2.)/dx, 0.0_dp)
        else
           rCICvec(:)=max(1.0-abs(x(i,:)-SN_pos(k,:))/dx, 0.0_dp)
        endif
        ! Compute distance, in user units, from SN source
        rvec(:)=x(i,:)-SN_pos(k,:)            !  Distance [uu] from SN src
        if(source_odd) rvec(:)=rvec(:)+dx/2.
        rCIC=product(rCICvec)                 !  Zero if outside CIC cloud

        ! A: Light injection:---------------------------------------------
        if(rCIC.gt.0. .and. rt) then            ! Add if cell in CIC cloud
           call getAgeGyr(0.d0, age)               !                [Gyrs]
           age=log10(max(age, sed_minAge))         !         In case age=0
           call inp_SED_table(age, SN_metallicity(k), 1, sed_I)
           sed_I = 10.d0**sed_I * SN_stellar_mass_solar(k) ! [# ph s-1]
           if(.not. fakeIt) write(*,*) &
                'Photons/sec: ',sed_I(1)*rCIC,dt_cgs/31536000.d0
           do ip=1,nSEDPacs
              uu(i,iPac(ip)) = uu(i,iPac(ip))  &
                   + sed_I(ip)*rCIC/vol_cgs*dt_cgs/scale_np
           end do
        endif

        ! B: SN injection:------------------------------------------------
        !if((rCIC.gt.0.) .and. (.not. SN_done(k))                        &
        !     .and. (SN_time_uu(k).ge.0.) .and. (SN_time_uu(k).le.t)) then
        if(rCIC.gt.0. .and. SN_GO) then
           ! Inject remaining stellar mass
           meject_solar = &
                SN_stellar_mass_solar(k)-mlost_sol_rl-SN_mremnant_solar(k)
           meject_cgs = meject_solar * m_sun
           uu(i,1) = uu(i,1) + meject_cgs * rCIC / vol_cgs / scale_d
           uu(i,ndim+2) = uu(i,ndim+2)  &
                        + SN_E_cgs(k)*rCIC/vol_cgs /(scale_d*scale_v**2)
           if(.not. fakeIt) then
              print *,'***************************************'
              print *,'SN BLAST at t = ',t,' !!!!! '
              print *,'SN ENERGY     = ',SN_E_cgs(k), rCIC, ' erg'
              print *,'SN MASS EJECT = ',meject_solar, ' msun'
              print *,'***************************************'
           endif
           ! Probably need something better here..put SN as kinetic E
        endif

        ! C: Wind injection:----------------------------------------------
        if(SN_wind(k) .and. is_in_wind_sphere(rvec,w))  then
           drho = w_dm_cgs**(ndim/3.) * w/vol_cgs/scale_d ! Injection [uu]
           if(drho/uu(i,1) .gt. courant_wind ) dt_ok=.false.
           ! Inject wind into wind domain:
           uu(i,1)=uu(i,1)+drho                            !   Inject mass
           if(metal) &                                     ! Inject metals
                uu(i,imetal)=uu(i,imetal)+drho*10.d0**SN_metallicity(k)
           uu(i,ndim+2) = uu(i,ndim+2) & ! Inject thermal + kinetic energy
                + w * (w_dET_cgs+w_dEk_cgs) /vol_cgs /(scale_d*scale_v**2) 
           ! Inject momentum:             
           r=sqrt(sum(rvec**2))              !  uu distance from SN source
           if(r .gt. dx_levelmax/100.) then  ! Dont inject in central cell
              pmag=w*w_dp_cgs/vol_cgs/scale_v!      Momentum (p) magnitude
              rvec=rvec/r                    !        Momentum unit vector
              uu(i,2:ndim+1) = uu(i,2:ndim+1) + rvec*pmag ! Mom. injection
           endif
        endif
     end do
     
  end do

END SUBROUTINE stellar_feedback_vsweep

!*************************************************************************
SUBROUTINE prepare_Wind_Weights(nrays)

! Prepare wind weights for an odd number of cells, i.e. where the source
! is at the center of a cell 
! nrays => Number of rays shot from center into sphere per cpu to generate
!          sphere of wind weights. One million rays per cpu is acceptable
!          but a bit time consuming.
!-------------------------------------------------------------------------
  use amr_commons
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::nRays
  real(dp)::twoPi    
  integer,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer,save,dimension(IRandNumSize)::localseed=-1
  integer::iseed=1,i,nix,niy,niz,ix,iy,iz,ierr
  real(dp)::RandNum,theta,phi,r,rx,ry,rz
!-------------------------------------------------------------------------
  twoPi=4.*asin(1.d0)                                  !              2*pi
  nix=2*SN_r_wind ; niy=2*SN_r_wind ; niz=2*SN_r_wind
  if(source_odd) then
     nix=nix+1 ; niy=niy+1 ; niz=niz+1
  endif
  if(ndim .lt. 3) niz=1
  if(ndim .lt. 2) niy=1
  if(.not. allocated(Wgrid)) then
     allocate(Wgrid(nix, niy, niz))
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  endif
  allocate(Wgrid_all(nix, niy, niz))
  Wgrid=0.
  phi=twoPi/4.   ;    theta=0
  do i=1,nRays
     call ranf(localseed,RandNum)                 
     r=RandNum                                         !             [0,1]
     if(ndim .gt. 1) then
        call ranf(localseed,RandNum)
        theta=RandNum*twoPi                            !           [0,2pi]
     endif
     if(ndim .gt. 2) then
        call ranf(localseed,RandNum)
        phi=RandNum*twoPi                              !           [0,2pi]
     endif
     rx = r * sin(phi) * cos(theta)
     ry=-1 ; rz=-1
     if(ndim .gt. 1) ry = r * sin(phi) * sin(theta)
     if(ndim .gt. 2) rz = r * cos(phi)
     if(source_odd) then
        ix = floor((SN_r_wind+0.5)*(rx + 1))+1
        iy = floor((SN_r_wind+0.5)*(ry + 1))+1
        iz = floor((SN_r_wind+0.5)*(rz + 1))+1
     else
        ix = floor(SN_r_wind*(rx + 1))+1
        iy = floor(SN_r_wind*(ry + 1))+1
        iz = floor(SN_r_wind*(rz + 1))+1
     endif
     Wgrid(ix,iy,iz) = Wgrid(ix,iy,iz)+1.
  end do
#ifndef WITHOUTMPI
  Wgrid_all = 0.
  call MPI_ALLREDUCE(Wgrid, Wgrid_all, nix**ndim,  &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  Wgrid = Wgrid_all
  deallocate(Wgrid_all)
#endif     
  Wgrid=Wgrid/(nRays*ncpu)

  !if(myid==1)then                            ! Just for debugging
  !   open(10, file='wind_weights.dat', form='unformatted')
  !   
  !   write(10) Wgrid
  !   print*,'Weight grid:===================='
  !   write(*,111)Wgrid(:,13,1)
  !   write(*,111)Wgrid(:,12,1)
  !   write(*,111)Wgrid(:,11,1)
  !   write(*,111)Wgrid(:,10,1)
  !   write(*,111)Wgrid(:,9,1)
  !   write(*,111)Wgrid(:,8,1)
  !   write(*,111)Wgrid(:,7,1)
  !   write(*,111)Wgrid(:,6,1)
  !   write(*,111)Wgrid(:,5,1)
  !   write(*,111)Wgrid(:,4,1)
  !   write(*,111)Wgrid(:,3,1)
  !   write(*,111)Wgrid(:,2,1)
  !   write(*,111)Wgrid(:,1,1)
  !   print*,'================================'
  !   print*,'Total weight of grid is ',sum(Wgrid)
  !endif

111 format(40f10.6)
END SUBROUTINE prepare_Wind_Weights

!************************************************************************
FUNCTION is_in_wind_sphere(location, weight)

! Returns true if cell is within the 'input region' of a SN source.
! location => Location of cell center relative to SN source (uu)
! weight  <=  Weight of cell within region
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::location(ndim)
  real(dp),intent(out)::weight
  logical::is_in_wind_sphere
  integer::i,ix(3)=1
!------------------------------------------------------------------------
  is_in_wind_sphere=.false.
  weight=0.d0
  if(.not. SN_wind(1)) return
  if(source_odd) then
     ix(1:ndim)=nint(location(1:ndim)/dx_levelmax)
     if(maxval(abs(ix)) .gt. SN_r_wind) RETURN   ! Outside wind domain
     ix(1:ndim)=ix(1:ndim)+SN_r_wind+1
  else
     ix(1:ndim)=location(1:ndim)/dx_levelmax + SN_r_wind + 1
     if(maxval(ix) .gt. SN_r_wind * 2) RETURN   ! Outside wind domain
     if(minval(ix) .lt. 1)             RETURN   ! Outside wind domain
  endif
  is_in_wind_sphere=.true.
  weight=Wgrid(ix(1),ix(2),ix(3))
END FUNCTION is_in_wind_sphere

END MODULE ideal_SN_module

!*************************************************************************
SUBROUTINE add_cont_sources(ilevel,dt)

! Inject any continuous sources into the grid.
!-------------------------------------------------------------------------
  use amr_parameters,only:dp
  use amr_commons,only:t,myid
  use rt_parameters
  use ideal_SN_module
  implicit none
  integer::ilevel,i
  real(dp)::dt
  logical::dt_ok
  logical,save::first_step_in_main=.true.
!-------------------------------------------------------------------------
  if( .not. SN_done(1) .and. SN_time_uu(1).ge.0.                         &
       .and. SN_time_uu(1).le.t                                          &
       .and. .not. SN_wait) then
     ! Start triggering a SN event. Before injecting the SN energy, we
     ! freeze the simulation until the next main timestep to allow for 
     ! dumping of output. Then we inject the SN energy but still keep
     ! zero timestep to again allow for output. After that things go back
     ! to normal.
     doDump=.true.  ! Tell amr-step to to output when next at coarse step
     SN_wait=.true. ! Keep zero timestep
  else if(SN_wait .and. first_step_in_main) then
     !We're at the beginning of a main step, so ok to inject the SN
     SN_wait=.false. ! Stop keeping zero timestep
     SN_GO=.true.    ! Inject SN energy ASAP
     doDump=.true.   ! Tell amr_step to do output (and keep dt=0)
     SN_done(1)=.true.                  
     if(rt) rt_cool_only=.true.
     rt=.false.      ! Turn off RT after SN event
  endif

  call stellar_feedback(ilevel, t, dt, dt_ok, .false.)

  first_step_in_main=.false.
  if(ilevel==levelmin) first_step_in_main=.true. ! Next time we're here  !
                                                 ! will be first step in !
                                                 ! the main step         !
END SUBROUTINE add_cont_sources



