!RT_STAR_SIMPLE patch (SS):
! Simple model for stellar emission. 
! Each star emits I=9.1691d13 photons per sec per gram, times an escape
! fraction, until it is 10Myrs old at which time it completely stops
! shining.
!************************************************************************
! These are modules for reading, integrating and interpolating RT-relevant
! values from spectral tables, specifically SED (spectral energy 
! distrbution) tables for stellar particle sources, as functions of age 
! and metallicity, and UV-spectrum tables, as functions of redshift.
! The modules are:
! spectrum_integrator_module
!    For dealing with the integration of spectra in general
! SED_module
!    For dealing with the SED data
! UV_module
!    For dealing with the UV data
!_________________________________________________________________________

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for integrating a wavelength-dependent spectrum
!_________________________________________________________________________
!
MODULE spectrum_integrator_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  implicit none

  PUBLIC integrateSpectrum, f1, fLambda, fdivLambda, fSig, fSigLambda,   &
         fSigdivLambda

  PRIVATE   ! default

CONTAINS

!************************************************************************
FUNCTION integrateSpectrum(X, Y, N, e0, e1, species, func, doPrint)

! Integrate spectral weighted function in energy interval [e0,e1]
! X      => Wavelengths [angstrom]
! Y      => Spectral intensity per angstrom at wavelenghts [XX A-1] 
! N      => Length of X and Y
! e0,e1  => Integrated interval [ev]
! species=> ion species, used as an argument in fx
! func   => Function which is integrated (of X, Y, species)
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs,eV_to_erg, hp
  use amr_commons,only:myid
  real(dp):: integrateSpectrum, X(N), Y(N), e0, e1
  integer :: N, species
  interface
     real(dp) function func(wavelength,intensity,species)
       use amr_parameters,only:dp
       real(dp)::wavelength,intensity
       integer::species
     end function func
  end interface!----------------------------------------------------------
  real(dp),dimension(:),allocatable:: xx, yy, f
  real(dp):: la0, la1
  integer :: i
  logical,optional::doPrint
!-------------------------------------------------------------------------
  integrateSpectrum=0.
  if(N .le. 2) RETURN
  ! Convert energy interval to wavelength interval
  la0 = X(1) ; la1 = X(N)
  if(e1.gt.0) la0 = max(la0, 1.d8 * hp * c_cgs / e1 / eV_to_erg)                
  if(e0.gt.0) la1 = min(la1, 1.d8 * hp * c_cgs / e0 / eV_to_erg)
  if(la0 .ge. la1) RETURN         
  ! If we get here, the [la0, la1] inverval is completely within X
  allocate(xx(N)) ; allocate(yy(N)) ; allocate(f(N))
  xx =  la0   ;   yy =  0.   ;   f = 0.
  i=2
  do while ( i.lt.N .and. X(i).le.la0 )
     i = i+1                      !              Below wavelength interval
  enddo                           !   X(i) is now the first entry .gt. la0 
  ! Interpolate to value at la0
  yy(i-1) = Y(i-1) + (xx(i-1)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i-1)  = func(xx(i-1), yy(i-1), species)
  do while ( i.lt.N .and. X(i).le.la1 )              ! Now within interval
     xx(i) = X(i) ; yy(i) = Y(i) ; f(i) = func(xx(i),yy(i),species)
     i = i+1
  enddo                          ! i=N or X(i) is the first entry .gt. la1
  xx(i:) = la1                   !             Interpolate to value at la1
  yy(i) = Y(i-1) + (xx(i)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i)  = func(xx(i),yy(i),species)

  !if(present(doPrint)) then
  !   if(doprint) then
  !      write(*,*) e0,e1,la0,la1,N
  !      write(*,*) '***************'
  !      do i=1,N
  !         write(*,*) xx(i),f(i),yy(i)
  !      end do
  !      stop
  !   endif
  !endif

  integrateSpectrum = trapz1(xx,f,i)
  deallocate(xx) ; deallocate(yy) ; deallocate(f)

END FUNCTION integrateSpectrum

!************************************************************************
! FUNCTIONS FOR USE WITH integrateSpectrum:
! lambda  => wavelengths in Angstrom
! f       => function of wavelength (a spectrum in some units)
! species => 1=HI, 2=HeI or 3=HeII
!_________________________________________________________________________
FUNCTION f1(lambda, f, species)
  real(dp):: f1, lambda, f
  integer :: species
  f1 = f
END FUNCTION f1

FUNCTION fLambda(lambda, f, species)
  real(dp):: fLambda, lambda, f
  integer :: species
  fLambda = f * lambda
END FUNCTION fLambda

FUNCTION fdivLambda(lambda, f, species)
  real(dp):: fdivlambda, lambda, f
  integer :: species
  fdivLambda = f / lambda
END FUNCTION fdivLambda

FUNCTION fSig(lambda, f, species)
  real(dp):: fSig, lambda, f
  integer :: species
  fSig = f * getCrosssection_Hui(lambda,species)
END FUNCTION fSig

FUNCTION fSigLambda(lambda, f, species)
  real(dp):: fSigLambda, lambda, f
  integer :: species
  fSigLambda = f * lambda * getCrosssection_Hui(lambda,species)
END FUNCTION fSigLambda

FUNCTION fSigdivLambda(lambda, f, species)
  real(dp):: fSigdivLambda, lambda, f
  integer :: species
  fSigdivLambda = f / lambda * getCrosssection_Hui(lambda,species)
END FUNCTION fSigdivLambda
!_________________________________________________________________________

!*************************************************************************
FUNCTION trapz1(X,Y,N)

! Integrates function Y(X) along the whole interval 1..N, using a very 
! simple staircase method and returns the result.
!-------------------------------------------------------------------------
  integer :: N,j
  real(dp):: trapz1
  real(dp):: X(N),Y(N)
!-------------------------------------------------------------------------
  trapz1=0.
  if (N.le.1) RETURN
  do j=2,N
     trapz1= trapz1 + abs(X(j)-X(j-1)) * (Y(j)+Y(j-1)) / 2.
  end do
END FUNCTION trapz1

!************************************************************************
FUNCTION getCrosssection_Hui(lambda, species)

! Gives an atom-photon cross-section of given species at given wavelength, 
! as given by Hui and Gnedin (1997).
! lambda  => Wavelength in angstrom
! species => 1=HI, 2=HeI or 3=HeII
! returns :  photoionization cross-section in cm^2
!------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
  real(dp)          :: lambda, getCrosssection_Hui
  integer           :: species
  real(dp)          :: E, E0, cs0, P, ya, yw, y0, y1, x, y
!------------------------------------------------------------------------
  E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev
  if ( E .lt. ionEvs(species) ) then            ! below ionization energy
     getCrosssection_Hui=0.
     RETURN
  endif
  select case (species)
  case (1) ! HI
     E0 = 4.298d-1 ; cs0 = 5.475d-14    ; P  = 2.963 
     ya = 32.88    ; yw  = 0            ; y0 = 0         ; y1 = 0
  case (2) ! HeI
     E0 = 1.361d1    ; cs0 = 9.492d-16  ; P  = 3.188 
     ya = 1.469      ; yw  = 2.039      ; y0 = 0.4434    ; y1 = 2.136
  case (3) ! HeII
     E0 = 1.720      ; cs0 = 1.369d-14  ; P  = 2.963 
     ya = 32.88      ; yw  = 0          ; y0 = 0         ; y1 = 0
  end select

  x = E/E0 - y0
  y = sqrt(x**2+y1**2)

  getCrosssection_Hui = &
       cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P

END FUNCTION getCrosssection_Hui


END MODULE spectrum_integrator_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for Stellar Energy Distribution table. 
!_________________________________________________________________________
!
MODULE SED_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  use rt_parameters,only:nPacs
  implicit none

  PUBLIC nSEDpacs, SED_minage, Ipac_fractions                            &
       ,init_SED_table,update_SED_PacProps,star_RT_feedback

  PRIVATE   ! default

  ! Light properties for different spectral energy distributions----------
  integer::nSEDpacs=nPacs               ! All pacs are SED pacs by default 
  real(dp)::SED_minAge=0.  
  real(dp)::I_star=9.1691d13            ! Photons per sec per gram
  real(dp),dimension(nPacs)::Ipac_fractions=1./nPacs! I per package
  real(dp)::deathAgeSec=3.15569d+14              ! Turn off at 10 Myr
  ! ---------------------------------------------------------------------

CONTAINS

!************************************************************************
SUBROUTINE init_SED_table()

! Do nothing
!------------------------------------------------------------------------
END SUBROUTINE init_SED_table

!************************************************************************
SUBROUTINE update_SED_PacProps()

! Do nothing.
!------------------------------------------------------------------------
END SUBROUTINE update_SED_PacProps

!************************************************************************
SUBROUTINE star_RT_feedback(ilevel, dt)

! This routine adds photons from radiating stars to appropriate cells in 
! the  hydro grid.
! ilegel =>  grid level in which to perform the feedback
! ti     =>  initial time for the timestep (code units)
! dt     =>  real timestep length in code units
!------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use rt_parameters
  integer:: ilevel
  real(dp):: dt
  integer:: igrid, jgrid, ipart, jpart, next_part
  integer:: i, ig, ip, npart1, npart2, icpu
  integer,dimension(1:nvector),save:: ind_grid, ind_part, ind_grid_part
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return ! number of grids in the level
  if(nstar_tot .le. 0 ) return
  if(verbose)write(*,111)ilevel
  ! Gather star particles only.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel) ! grid index
     ig=0                     ! index of grid with stars (in ind_grid)
     ip=0                     ! index of star particle   (in ind_part)
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)   ! Number of particles in the grid
        npart2=0              ! number of selected (i.e. star) particles
        
        ! Count star particles in the grid
        if(npart1 > 0)then
          ipart = headp(igrid)        ! particle index
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle       <--- Very important !!!
              next_part = nextp(ipart)
              if(idp(ipart) .gt. 0 .and. tp(ipart) .ne. 0.d0 &
                   .and. tpp(ipart) .le. t_proper) then 
                 npart2 = npart2+1     ! only stars
              endif
              ipart = next_part        ! Go to next particle
           end do
        endif

        ! Gather star particles within the grid
        if(npart2 > 0)then        
           ig = ig+1
           ind_grid(ig) = igrid
           ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle      <--- Very important !!!
              next_part = nextp(ipart)
              if(ig == 0)then
                 ig = 1
                 ind_grid(ig) = igrid
              end if
              ! Select only star particles
              if(idp(ipart) .gt. 0 .and. tp(ipart) .ne. 0.d0 &
                   .and. tpp(ipart) .le. t_proper) then
                 ip = ip+1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig ! points to grid a star is in  
              endif
              if(ip == nvector)then
                 if(rt_pp .or. rt_coupling_pp) then
                    call star_RT_vsweep_pp( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 else
                    call star_RT_vsweep( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 endif
                 ip = 0
                 ig = 0
              end if
              ipart = next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid = next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip > 0) then
        if(rt_pp .or. rt_coupling_pp) then
           call star_RT_vsweep_pp( &
                ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
        else
           call star_RT_vsweep( &
                ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
        endif
     endif
  end do 
  ! End loop over cpus

111 format('   Entering star_rt_feedback for level ',I2)

END SUBROUTINE star_RT_feedback

!*************************************************************************
!*************************************************************************
! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

!************************************************************************
SUBROUTINE star_RT_vsweep(ind_grid, ind_part, ind_grid_part, ng, np,    &
                          dt, ilevel)
! This routine is called by subroutine star_rt_feedback.
! Each star particle dumps a number of photons into the nearest grid cell
! using array unew.
! Radiation is injected into cells at level ilevel, but it is important
! to know that ilevel-1 cells may also get some radiation. This is due
! to star particles that have just crossed to a coarser level. 
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  star indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of stars
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use hydro_commons
  use rt_parameters
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,nPacs),save::irad
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0,igd=0,icd=0
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  & 
       ,scale_fp, scale_inp, scale_nphot, ageUU, deathAgeUU, overageUU   &
       ,dt_loc
  real(dp),parameter::vol_factor=2**ndim  ! Vol factor for ilevel-1 cells
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)
  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = I_star*rt_esc_frac * scale_d*scale_t / scale_np / vol_loc
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50

  ! Lower left corners of 3x3x3 grid-cubes (with given grid in center)
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - 3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(&
          ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)
  ! now nbors_father cells are a cube of 27 cells with ind_cell in the 
  ! middle and nbors_father_grids are the 8 grids that contain these 27
  ! cells (though only one of those is fully included in the cube)

  ! Rescale position of stars to positions within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx 
        ! so 0<x<2 is bottom cell, ...., 4<x<6 is top cell 
     end do
  end do
  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! So id=0-5 is the cell (in the 
     end do                    ! 3x3x3 supercube) containing the star
  end do
  ! Compute parent grids
  do idim = 1, ndim
     do j = 1, np
        igd(j,idim) = id(j,idim)/2 ! must be 0, 1 or 2
     end do
  end do
  do j = 1, np
     kg(j) = 1 + igd(j,1) + 3*igd(j,2) + 9*igd(j,3) ! 1 to 27
  end do
  do j = 1, np
     igrid(j) = son(nbors_father_cells(ind_grid_part(j),kg(j))) 
     ! grid (not cell) containing the star 
  end do

  ! Check if particles are entirely in level ilevel.
  ! This should always be the case in post-processing
  ok(1:np) = .true.
  do j = 1, np
     ok(j) = ok(j) .and. igrid(j) > 0
  end do ! if ok(j) is true then particle j's cell contains a grid. 
  ! Otherwise it is a leaf, and we need to fill it with radiation.

  ! Compute parent cell position within it's grid
  do idim = 1, ndim
     do j = 1, np
        if( ok(j) ) then
           icd(j,idim) = id(j,idim) - 2*igd(j,idim) ! 0 or 1
        end if
     end do
  end do
  do j = 1, np
     if( ok(j) ) then
        icell(j) = 1 + icd(j,1) + 2*icd(j,2) + 4*icd(j,3) ! 1 to 8
     end if
  end do

  ! Compute parent cell adress and particle radiation contribution
  do j = 1, np
     ageUU = t - tp(ind_part(j))              !     Age at end of timestep
     overAgeUU = ageUU-deathAgeUU             ! 'Surplus' age at end of dt
     !Possibilities: born before and dead after dt, born after dt, 
     !               dead within dt, born and dead within dt
     dt_loc = max(0.d0,min(dt,ageUU,dt-overAgeUU,deathAgeUU))
     irad(j,1:nSEDpacs)= &             !             [#/cc, in Code Units]
          Ipac_fractions * mp(ind_part(j)) * dt_loc * scale_inp

     if(showSEDstats) then
        step_nPhot = step_nPhot+irad(j,1)*scale_nPhot
        step_nStar = step_nStar+dt_loc
     endif

     if( ok(j) )then
        indp(j) = ncoarse + (icell(j)-1)*ngridmax + igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
     end if
  end do
  ! Increase photon density in cell due to stellar radiation
  do j=1,np
     if( ok(j) ) then                                      !   ilevel cell
        do ip=1,nSEDPacs
           unew(indp(j),iPac(ip)) = unew(indp(j),iPac(ip)) + irad(j,ip)
        end do
     else                                                  ! ilevel-1 cell
        do ip=1,nSEDPacs
           unew(indp(j),iPac(ip)) = unew(indp(j),iPac(ip)) + irad(j,ip)  &
                / vol_factor
        end do
     endif
  end do

END SUBROUTINE star_RT_vsweep

!*************************************************************************
SUBROUTINE star_RT_vsweep_pp(ind_grid, ind_part, ind_grid_part, ng, np,  &
                             dt, ilevel)
! This routine is called by subroutine star_rt_feedback,
! in the case of rt-post-processing. The difference from doing coupled 
! RT is that in the post-processing case, the pointer from a grid to a 
! particle is always true, so we can trust completely that 
! ind_grid(ind_grid_part(j)) refers to the grid that particle ind_part(j)
! resides in, and we don't need to bother searching for the particle's 
! grid.  Each star particle dumps a number of photons into the nearest 
! grid cell using array unew.
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  star indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of stars
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use hydro_commons
  use rt_parameters
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !---------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  real(dp),dimension(1:nvector,nPacs),save::irad
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0
  integer ,dimension(1:nvector),save::icell,indp
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_np   & 
       ,scale_fp, scale_inp, scale_nphot, ageUU, deathAgeUU, overageUU   &
       ,dt_loc
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)

  deathAgeUU = deathAgeSec/scale_t             ! Maximum age in code units

  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = I_star*rt_esc_frac * scale_d*scale_t / scale_np / vol_loc
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50

  ! Lower left corners of parent grid
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - dx
     end do
  end do

  ! Rescale position of star to position within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx 
        ! so 0<x<1 is bottom cell, ...., 1<x<2 is top cell
     end do
  end do

  ! Compute position of parent cell within parent grid
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! id=0-1 is cell in parent containing star
     end do
  end do
  do j = 1, np
     icell(j) = 1 + id(j,1) + 2*id(j,2) + 4*id(j,3) ! 1 to 8
  end do

  do j = 1, np ! Calculate emissivity of particle into cell
     ageUU = t - tp(ind_part(j))              !     Age at end of timestep
     overAgeUU = ageUU-deathAgeUU             ! 'Surplus' age at end of dt
     !Possibilities: born before and dead after dt, born after dt, 
     !               dead within dt, born and dead within dt
     dt_loc = max(0.d0,min(dt,ageUU,dt-overAgeUU,deathAgeUU))
     irad(j,1:nSEDpacs)= &             !             [#/cc, in Code Units]
          Ipac_fractions * mp(ind_part(j)) * dt_loc * scale_inp

     if(showSEDstats) then
        step_nPhot = step_nPhot+irad(j,1)*scale_nPhot
        step_nStar = step_nStar+dt_loc
     endif

     indp(j)= ncoarse + (icell(j)-1)*ngridmax + ind_grid(ind_grid_part(j))

  end do

  ! Update RT variables due to feedback
  do j=1,np
     do ip=1,nSEDPacs
        unew(indp(j),iPac(ip)) = unew(indp(j),iPac(ip)) + irad(j,ip)
     end do
  end do

END SUBROUTINE star_RT_vsweep_pp

END MODULE SED_module


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for UV table of redshift dependent photon fluxes, cross sections
! and photon energies, per photon package. This is mainly useful for 
! calculating heating rates due to the homogeneous UV background in the 
! cooling-module, but can also be used in an implementation of radiative
! transfer of the UV background, where we then assign intensities, 
! average cross sections and energies to UV photon packages.
! 
! There are two tables of UV properties:
! UV_rates_table: Integrated ionization and heating rates for a 
!                 homogeneous UV background.
! UV_pacs_table:  UV properties integrated over each photon package 
!                 wavelength interval: photon flux, average energy, 
!                 average cross section, energy weighted cross section.
!_________________________________________________________________________
!
MODULE UV_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  use hydro_commons,only:nvar
  use spectrum_integrator_module
  use rt_parameters,only:c_cgs, eV_to_erg, hp, nIons, ionEVs, nPacs

  implicit none

  PUBLIC nUVpacs, iUVpacs, iUVvars                                       &
       , init_UV_tables, inp_UV_rates_table, inp_UV_pacs_table           &
       , UV_minz, UV_maxz

  PRIVATE   ! default

  ! UV properties for different redshifts---------------------------------
  integer::UV_nz
  real(dp),allocatable,dimension(:)::UV_zeds
  ! Redshift interval
  real(dp)::UV_minz, UV_maxz  
  ! Table of integrated ionization and heating rates per ion species
  ! (UV_nz, nIons, 2):(z, nIon, ionization rate [# s-1]: Hrate [erg s-1]) 
  real(dp),allocatable,dimension(:,:,:)::UV_rates_table
  integer::nUVpacs=0             !  # of UV photon packages
  ! UV pac indexes among nPac pacs, UV Np indexes among nVar vars.
  integer,allocatable::iUVpacs(:),iUVvars(:) 
  ! Table of photon package props (UV_nz, nUVpacs, 2+2*nIons): 
  !(z, pac, flux [#/cm2/s]: egy[eV]: nIons*csn [cm-2]: nIons*cse [cm-2])
  real(dp),allocatable,dimension(:,:,:)::UV_pacs_table
  ! The tables do not have constant redshift intervals, so we need to
  ! first locate the correct interval when doing interpolation.
  real(dp),parameter::fourPi=12.566371
  logical::is_UV_initiated=.false.
  ! ----------------------------------------------------------------------

CONTAINS

!*************************************************************************
SUBROUTINE init_UV_tables()
! Initiate UV props table, which gives photon fluxes and average photon 
! cross sections and energies of each photon package as a function of z.
! The data is read from a file specified by uv_file containing: 
! a) number of redshifts
! b) redshifts (increasing order)
! c) number of wavelengths
! d) wavelengths (increasing order) [Angstrom]
! e) fluxes per (redshift,wavelength) [photons cm-2 s-1 A-1 sr-1]
!-------------------------------------------------------------------------
  use amr_commons,only:myid
  use rt_parameters
  use SED_module
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer:: nLs                                 ! # of bins of wavelength
  real(dp),allocatable  :: Ls(:)                ! Wavelengths
  real(dp),allocatable  :: UV(:,:)              ! UV f(lambda,z)
  real(dp),allocatable  :: tbl(:,:), tbl2(:,:)
  integer::i,ia,iz,ip,ii,dum,locid,ncpu,ierr
  logical::ok
  real(dp)::da, dz, pL0,pL1
!------------------------------------------------------------------------
  if(is_UV_initiated) RETURN                 ! Table is already initiated 
  is_UV_initiated=.true.
  if(UV_FILE .eq. '') &
       call get_environment_variable('RAMSES_UV_FILE', UV_FILE)
  inquire(file=TRIM(uv_file), exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then 
        write(*,*)'Cannot access UV file...',TRIM(uv_file)
        write(*,*)'File '//TRIM(uv_file)//' not found'
        write(*,*)'You need to set the RAMSES_UV_FILE envvar' // &
                  ' to the correct path, or use the namelist'
     endif
     call clean_stop
  end if
  ! Read redshifts, wavelengths and spectra
  open(unit=10,file=TRIM(uv_file),status='old',form='unformatted')
  read(10) UV_nz ; allocate(UV_zeds(UV_nz)) ; read(10) UV_zeds(:)
  read(10) nLs   ; allocate(Ls(nLs))        ; read(10) Ls(:)
  allocate(UV(nLs,UV_nz))                   ; read(10) UV(:,:)
  close(10)
  ! If mpi then share the UV integration between the cpus:
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,locid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  locid=0 ; ncpu=1
#endif

  ! Shift the highest z in the table (10) to reionization epoch,
  ! so that we start injecting at z_reion
  if(z_reion .gt. UV_zeds(UV_nz-1))  UV_zeds(UV_nz) = z_reion
 
  ! Initialize rates table------------------------------------------------
  if(rt_UV_hom) then
     allocate(UV_rates_table(UV_nz, nIons, 2))
     allocate(tbl(UV_nz, 2))
     do ii = 1, nIons
        tbl=0.
        do iz = locid+1,UV_nz,ncpu
           tbl(iz,1)= getUV_Irate(Ls,UV(:,iz),nLs,ii)
           tbl(iz,2)= getUV_Hrate(Ls,UV(:,iz),nLs,ii)
        end do
#ifndef WITHOUTMPI
        allocate(tbl2(UV_nz,2))
        call MPI_ALLREDUCE(tbl, tbl2, UV_nz*2,  &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        tbl = tbl2
        deallocate(tbl2)
#endif     
        UV_rates_table(:,ii,:)=tbl 
     end do
     deallocate(tbl)
     if (myid==1) call write_UVrates_table
  endif

  ! Note: Need to take special care of the highest redshift. Some UV     !
  !       spectra (read: Faucher-Giguere) have zero fluxes at the highest!
  !       redshift, and integrations of energies and cross sections of   !
  !       these fluxes results in a division by zero and NaNs.           !
  !       The approach taken is then to keep zero flux at the highest    !
  !       redshift but otherwise use the same values of energies and     !
  !       cross sections as the next-highest redshift.                   !

  ! Initialize photon packages table--------------------------------------
  allocate(UV_pacs_table(UV_nz, nUVpacs, 2+2*nIons))
  allocate(tbl(UV_nz, 2+2*nIons))
  do ip = 1,nUVpacs             !                     Loop photon packages
     tbl=0.
     pL0 = pacL0(nSEDpacs+ip)   !     energy interval of photon package ip
     pL1 = pacL1(nSEDpacs+ip)   !
     do iz = locid+1,UV_nz,ncpu
        tbl(iz,1) =         getUVFlux(Ls,UV(:,iz),nLs,pL0,pL1)
        if(tbl(iz,1) .eq. 0.d0) cycle        ! Can't integrate zero fluxes
        tbl(iz,2) = getUVphotonEnergy(Ls,UV(:,iz),nLs,pL0,pL1)
        do ii = 1,nIons
           tbl(iz,1+ii*2)= getUVcsn(Ls,UV(:,iz),nLs,pL0,pL1,ii)
           tbl(iz,2+ii*2)= getUVcse(Ls,UV(:,iz),nLs,pL0,pL1,ii)
        end do
     end do
#ifndef WITHOUTMPI
     allocate(tbl2(UV_nz,2+2*nIons))
     call MPI_ALLREDUCE(tbl,tbl2,UV_nz*(2+2*nIons),&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     tbl = tbl2
     deallocate(tbl2)
#endif     
     if(tbl(UV_nz,1) .eq. 0.d0) &            !                   Zero flux
        tbl(UV_nz,2:)=tbl(UV_nz-1,2:)
     UV_pacs_table(:,ip,:)=tbl 
  end do
  UV_minz = UV_zeds(1) ; UV_maxz=UV_zeds(UV_nz)
  deallocate(tbl) ; deallocate(Ls) ; deallocate(UV)

  if (myid==1) call write_UVpacs_tables

END SUBROUTINE init_UV_tables

!************************************************************************
SUBROUTINE inp_UV_rates_table(z, ret)
! Compute UV heating and ionization rates by interpolation from table.
! z     => Redshift
! ret  <=  (nIons,2) interpolated values of all UV rates.
!          (1=ionization rate [# s-1], 2=heating rate [erg s-1]) 
!------------------------------------------------------------------------
  real(dp), intent(in):: z
  real(dp):: ret(nIons,2), dz0, dz1
  integer:: iz0, iz1
!------------------------------------------------------------------------
  ret=0. ; if(z .gt. UV_maxz) RETURN
  call inp_1d(UV_zeds, UV_nz, z, iz0, iz1, dz0, dz1)
  ret = dz0*UV_rates_table(iz1, :, :) + dz1*UV_rates_table(iz0, :, :)
END SUBROUTINE inp_UV_rates_table

!************************************************************************
SUBROUTINE inp_UV_pacs_table(z, ret)
! Compute UV properties by interpolation from table.
! z     => Redshift
! ret  <=  (nPac,2+2*nIons) interpolated values of all UV properties.
!          (1=photon flux [#/cm2/s]), 2=pac_egy [ev],3+2*i=pac_csn[cm-2], 
!          4+2*i=pac_cse [cm-2])
!------------------------------------------------------------------------
  real(dp), intent(in):: z
  real(dp):: ret(nUVpacs,2+2*nIons), dz0, dz1
  integer:: iz0, iz1
!------------------------------------------------------------------------
  ret=0. ; if(z .gt. UV_maxz) RETURN
  call inp_1d(UV_zeds, UV_nz, z, iz0, iz1, dz0, dz1)
  ret = dz0 * UV_pacs_table(iz1, :, :) + dz1 * UV_pacs_table(iz0, :, :)
END SUBROUTINE inp_UV_pacs_table

!*************************************************************************
! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

!*************************************************************************
FUNCTION getUV_Irate(X, Y, N, species)
! Compute and return photoionization rate per ion, given the UV background
! Y(X). Assumes X is in Angstroms and Y in [# cm-2 s-1 sr-1 A-1].
! returns: Photoionization rate in # s-1
!-------------------------------------------------------------------------
  real(dp):: getUV_Irate, X(N), Y(N)
  integer :: N, species
!-------------------------------------------------------------------------
  getUV_Irate = & 
       fourPi*integrateSpectrum(X,Y,N, ionEvs(species), 0.d0, species, fsig)
END FUNCTION getUV_Irate

!*************************************************************************
FUNCTION getUV_Hrate(X, Y, N, species)
! Compute and return heating rate per ion, given the UV background   
! Y(X). Assumes X is in Angstroms and Y in [# cm-2 s-1 sr-1 A-1].
! returns: Heating rate in erg s-1
!-------------------------------------------------------------------------
  real(dp):: getUV_Hrate, X(N), Y(N), e0
  integer :: N, species
  real(dp),parameter :: const1=fourPi*1.d8*hp*c_cgs
  real(dp),parameter :: const2=fourPi*eV_to_erg
!-------------------------------------------------------------------------
  e0=ionEvs(species)
  getUV_Hrate = &
        const1*integrateSpectrum(X,Y,N, e0, 0.d0, species, fsigDivLambda)&
       -const2*ionEvs(species) *                                         &
               integrateSpectrum(X,Y,N, e0, 0.d0, species, fsig) 
END FUNCTION getUV_Hrate

!*************************************************************************
FUNCTION getUVFlux(X, Y, N, e0, e1)
! Compute and return UV photon flux in energy interval (e0,e1) [eV]   
! in UV spectrum Y(X). Assumes X is in [A] and Y in # cm-2 s-1 sr-1 A-1.
! returns: Photon flux in # cm-2 s-1
!-------------------------------------------------------------------------
  real(dp):: getUVflux, X(N), Y(N), e0, e1
  integer :: N, species
!-------------------------------------------------------------------------
  species          = 1                   ! irrelevant but must be included
  getUVflux = fourPi*integrateSpectrum(X, Y, N, e0, e1, species, f1)
END FUNCTION getUVflux

!*************************************************************************
FUNCTION getUVphotonEnergy(X, Y, N, e0, e1)
! Compute average photon energy, in eV, in energy interval (e0,e1) [eV] in
! UV spectrum Y(X). Assumes X is in [A] and Y is # cm-2 s-1 sr-1 A-1.
!-------------------------------------------------------------------------
  real(dp):: getUVphotonEnergy, X(N), Y(N), e0, e1, norm
  integer :: N,species
  real(dp),parameter :: const=1.d8*hp*c_cgs/eV_to_erg    ! unit conversion
!-------------------------------------------------------------------------
  species      = 1                       ! irrelevant but must be included
  norm         = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getUVphotonEnergy = const * &
            integrateSpectrum(X, Y, N, e0, e1, species, fdivLambda) / norm
END FUNCTION getUVphotonEnergy

!*************************************************************************
FUNCTION getUVcsn(X, Y, N, e0, e1, species)
! Compute and return average photoionization cross-section [cm2] for given 
! energy interval (e0,e1) [eV] in UV spectrum Y. Assumes X is in Angstroms
! and that Y is photon flux per angstrom.
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  real(dp):: getUVcsn, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getUVcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getUVcsn = integrateSpectrum(X, Y, N, e0, e1, species, fSig)/norm
END FUNCTION getUVcsn

!************************************************************************
FUNCTION getUVcse(X, Y, N, e0, e1, species)
! Compute average energy weighted photoionization cross-section [cm2] for
! given energy interval (e0,e1) [eV] in UV spectrum Y. Assumes X is in 
! Angstroms and that Y is energy intensity per angstrom.
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  real(dp):: getUVcse, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getUVcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, fdivLambda)
  getUVcse = integrateSpectrum(X, Y, N, e0, e1, species, fSigdivLambda)  &
           / norm
END FUNCTION getUVcse

!************************************************************************
SUBROUTINE write_UVrates_table()
! Write the UV rates to a file (this is just in 
! debugging, to check if the UV spectra are being read correctly).
!------------------------------------------------------------------------
  character(len=128)::filename
  integer::i, j
!------------------------------------------------------------------------
  write(filename,'(A, I1, A)') 'UVrates.list'
  open(10, file=filename, status='unknown')
  write(10,*) UV_nz
  
  do i = 1,UV_nz
     write(10,900) UV_zeds(i), UV_rates_table(i,:,:)
  end do
  close(10)
900 format (f21.6, 20(1pe21.6))
END SUBROUTINE write_UVrates_table

!************************************************************************
SUBROUTINE write_UVpacs_tables()

! Write the UV photon package properties to files (this is just in 
! debugging, to check if the UV spectra are being read correctly).
!------------------------------------------------------------------------
  use rt_parameters,only:nPacs
  character(len=128)::filename
  integer::ip, i, j
!------------------------------------------------------------------------
  do ip=1,nUVPacs
     write(filename,'(A, I1, A)') 'UVtable', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) UV_nz
     
     do i = 1,UV_nz
        write(10,901)                                                   &
                UV_zeds(i)           ,                                  &
                UV_pacs_table(i,ip,1),    UV_pacs_table(i,ip,2),        &
                UV_pacs_table(i,ip,3),    UV_pacs_table(i,ip,4),        &
                UV_pacs_table(i,ip,5),    UV_pacs_table(i,ip,6),        &
                UV_pacs_table(i,ip,7),    UV_pacs_table(i,ip,8)
     end do
     close(10)
  end do
901 format (f21.6,   f21.6,   f21.6,   1pe21.6, &
          & 1pe21.6, 1pe21.6, 1pe21.6, 1pe21.6, 1pe21.6  )
END SUBROUTINE write_UVpacs_tables

END MODULE UV_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!************************************************************************
SUBROUTINE locate(xx,n,x,j)
! Locates position j of a value x in an ordered array xx of n elements
! After: xx(j) <= x <= xx(j+1) (assuming increasing order)
! j is lower bound, so it can be zero but not larger than n
!------------------------------------------------------------------------
  use amr_commons,only:dp
  integer ::  n,j,jl,ju,jm
  real(dp)::  xx(n),x
!------------------------------------------------------------------------
  jl = 0
  ju = n+1
  do while (ju-jl > 1) 
     jm = (ju+jl)/2
     if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
        jl = jm
     else
        ju = jm
     endif
  enddo
  j = jl
END SUBROUTINE locate

!************************************************************************
SUBROUTINE inp_1d(xax,nx,x,ix0,ix1,dx0,dx1)
! Compute variables by interpolation from table.
! xax      => Axis of x-values in table
! nx       => Length of x-axis
! x        => x-value to interpolate to
! ix0,ix1 <=  Lower and upper boundaries of x in xax
! dx0,dx1 <=  Weights of ix0 and ix1 indexes
!------------------------------------------------------------------------
  use amr_commons,only:dp
  real(dp), intent(in)::xax(nx), x
  real(dp):: x_step, dx0, dx1
  integer:: nx, ix0, ix1
!------------------------------------------------------------------------
  call locate(xax, nx, x, ix0)
  if(ix0 < 1) ix0=1
  if (ix0 < nx) then
     ix1  = ix0+1
     x_step = xax(ix1) - xax(ix0)
     dx0  = max(           x - xax(ix0), 0.0d0 ) / x_step
     dx1  = min(xax(ix1) - x           , x_step) / x_step
  else
     ix1  = ix0
     dx0  = 0.0d0 ;  dx1  = 1.0d0
  end if

  if (abs(dx0+dx1-1.0d0) .gt. 1.0d-5) then
     write(*,*) 'Screwed up the 1d interpolation ... '
     write(*,*) dx0+dx1
     call clean_stop
  end if
  !ret = dx0 * table(ix1, :) + dx1 * table(ix0, :)

END SUBROUTINE inp_1d

