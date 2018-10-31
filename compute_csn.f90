program star_luminosities

  ! generate photons emitted by star particles within a given domain

  ! use module_photon
  !use module_utils
  use module_domain
  !use module_random
  use module_constants
  use module_ramses

  implicit none

  type(domain)    :: emission_domain
  character(2000) :: parameter_file
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), star_L(:,:)
  integer(kind=4) :: i,nstars,narg

  real(kind=8),allocatable,dimension(:,:,:,:)::SED_table
  real(kind=8),allocatable,dimension(:)::SED_ages, SED_zeds
  real(kind=8),parameter::SED_dlgA=0.02d0
  real(kind=8)::SED_dlgZ
  real(kind=8)::SED_lgA0, SED_lgZ0
  integer::SED_nA, SED_nZ=8           ! Number of age bins and Z bins

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [PhotonsFromStars] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: sed_dir = './seds/'        ! File to the SED library, in the same format as for Ramses
  character(2000)           :: outputfile = 'PhotICs.dat' ! file to which outputs will be written
  character(2000)           :: repository = './'          ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1                ! ramses output number to use
  integer(kind=4)           :: nSEDgroups = 1             ! number of photon groups
  real(kind=8),allocatable  :: group_decomp(:)            ! energy intervals of the photon groups, in eV

  ! --- domain whithin which star particles will be selected (should be within computational domain used for RT). 
  character(10)             :: star_dom_type      = 'sphere'         ! shape type of domain  // default is sphere.
  real(kind=8),dimension(3) :: star_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: star_dom_rsp       = 0.3              ! radius of spher [code units]
  real(kind=8)              :: star_dom_size      = 0.3              ! size of cube [code units]
  real(kind=8)              :: star_dom_rin       = 0.0              ! inner radius of shell [code units]
  real(kind=8)              :: star_dom_rout      = 0.3              ! outer radius of shell [code units]
  real(kind=8)              :: star_dom_thickness = 0.1              ! thickness of slab [code units]

  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------



  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: PhotonsFromStars path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_PhotonsFromStars_params(parameter_file)
  if (verbose) call print_PhotonsFromStars_params
  ! --------------------------------------------------------------------------------------
  ! define domain within which stars may shine
  ! --------------------------------------------------------------------------------------
  select case(star_dom_type)
  case('sphere')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, &
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),r=star_dom_rsp)
  case('shell')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, &
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),r_inbound=star_dom_rin,r_outbound=star_dom_rout)
  case('cube')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, & 
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),size=star_dom_size)
  case('slab')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, &
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),thickness=star_dom_thickness)
  end select
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! read star particles within domain
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> reading star particles'
  call ramses_read_stars_in_domain(repository,snapnum,emission_domain,star_pos,star_age,star_mass,star_vel,star_met)
  ! --------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------
  ! Initialize SED properties and compute luminosities of stellar particles
  ! --------------------------------------------------------------------------------------
  call init_SED_table(SED_table, SED_ages, SED_zeds)
  ! --------------------------------------------------------------------------------------

  nstars = size(star_age)
  allocate(star_L(nSEDgroups,nstars))
  do i=1,nstars
     call inp_sed_table(star_age(i), star_met(i), 1, .false., star_L(:,i))   !1 for luminosity, 2 for cumulative luminosty, 3 for Egy,  2+2*i for csn,  3+2*i for cse,  i=1,7 for HI,HeI,HeII,SiI,SiII,SiIII,SiIV
     star_L(:,i) = star_L(:,i)*star_mass(i)/msun
  end do

  print*, nSEDgroups


  ! --------------------------------------------------------------------------------------
  ! Write on file
  ! --------------------------------------------------------------------------------------
  open(unit=10, file=trim(outputfile), status='replace', form='unformatted', action='write')
  write(10) nSEDgroups
  write(10) nstars
  write(10) star_pos
  write(10) star_L
  close(10)



contains


  SUBROUTINE init_SED_table(SED_table, SED_ages, SED_zeds)

    ! Initiate SED properties table, which gives photon luminosities, 
    ! integrated luminosities, average photon cross sections and energies of 
    ! each photon group as a function of stellar population age and
    ! metallicity.  The SED is read from a directory specified by sed_dir.
    !-------------------------------------------------------------------------

    integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
    real(kind=8),allocatable::ages(:), Zs(:), Ls(:), rebAges(:)
    real(kind=8),allocatable::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
    real(kind=8),allocatable::tbl(:,:,:), tbl2(:,:,:), reb_tbl(:,:,:)
    integer::i,ia,iz,ip,ii,dum
    character(len=128)::fZs, fAges, fSEDs                        ! Filenames
    logical::ok,okAge,okZ
    real(kind=8)::dlgA, pL0, pL1, tmp
    integer::locid,ncpu2,ierr
    integer::nv=3+2*7  ! # vars in SED table: L,Lacc,egy,nions*(csn,cse)
    integer,parameter::tag=1132
    integer::dummy_io,info2
    real(kind=8),allocatable,dimension(:,:,:,:)::SED_table
    real(kind=8),allocatable,dimension(:)::SED_ages, SED_zeds
    real(kind=8), parameter:: Gyr2sec = 3.15569d+16       !       Gyr to sec conversion constant

    !-------------------------------------------------------------------------
    inquire(FILE=TRIM(sed_dir)//'/all_seds.dat', exist=ok)
    if(.not. ok)then 
       write(*,*)'Cannot access SED directory ',TRIM(sed_dir)
       write(*,*)'Directory '//TRIM(sed_dir)//' not found'
       write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
            ' to the correct path, or use the namelist.'
       stop
    end if
    write(fZs,'(a,a)')   trim(sed_dir),"/metallicity_bins.dat"
    write(fAges,'(a,a)') trim(sed_dir),"/age_bins.dat"
    write(fSEDs,'(a,a)') trim(sed_dir),"/all_seds.dat"
    inquire(file=fZs, exist=okZ)
    inquire(file=fAges, exist=okAge)
    inquire(file=fSEDs, exist=ok)
    if(.not. ok .or. .not. okAge .or. .not. okZ) then
       write(*,*) 'Cannot read SED files...'
       write(*,*) 'Check if SED-directory contains the files ',  &
            'metallicity_bins.dat, age_bins.dat, all_seds.dat'
       stop
    end if


    ! READ METALLICITY BINS-------------------------------------------------
    open(unit=10,file=fZs,status='old',form='formatted')
    read(10,'(i8)') nzs
    allocate(zs(nzs))
    do i = 1, nzs
       read(10,'(e14.6)') zs(i)
    end do
    close(10)
    ! READ AGE BINS---------------------------------------------------------
    open(unit=10,file=fAges,status='old',form='formatted')
    read(10,'(i8)') nAges
    allocate(ages(nAges))
    do i = 1, nAges
       read(10,'(e14.6)') ages(i)
    end do
    close(10)
    ages = ages*1.e-9                       !         Convert from yr to Gyr
    if(ages(1) .ne. 0.) ages(1) = 0.
    ! READ SEDS-------------------------------------------------------------
    open(unit=10,file=fSEDs,status='old',form='unformatted')
    read(10) nLs, dum
    allocate(Ls(nLs))
    read(10) Ls(:)
    allocate(SEDs(nLs,nAges,nzs))
    do iz = 1, nzs
       do ia = 1, nAges
          read(10) SEDs(:,ia,iz)
       end do
    end do
    close(10)

    ! Perform SED integration of luminosity, csn and egy per (age,Z) bin----
    allocate(tbl(nAges,nZs,nv))
    do ip = 1,nSEDgroups                                ! Loop photon groups
       tbl=0.
       pL0 = group_decomp(ip) ; pL1 = group_decomp(ip+1)! eV interval of photon group ip
       do iz = 1, nzs                                     ! Loop metallicity
          do ia = 1,nAges                                ! Loop age
             tbl(ia,iz,1) = getSEDLuminosity(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
             tbl(ia,iz,3) = getSEDEgy(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
             do ii=1,7
                tbl(ia,iz,2+ii*2) = getSEDcsn(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
                tbl(ia,iz,3+ii*2) = getSEDcse(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
             end do !End species loop
          end do ! End age loop
       end do ! End Z loop

       dlgA = SED_dlgA ; SED_dlgZ = -SED_nz
       call rebin_log(dlgA, SED_dlgZ                                       &
            , tbl(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, nv        &
            , reb_tbl, SED_nA, SED_nZ, rebAges, SED_Zeds)
       SED_nA=SED_nA+1                              ! Make room for zero age
       if(ip .eq. 1 ) allocate(SED_table(SED_nA, SED_nZ, nSEDgroups, nv))
       SED_table(1, :,ip,:) = reb_tbl(1,:,:)            ! Zero age properties
       SED_table(1, :,ip,2) = 0.                        !  Lacc=0 at zero age
       SED_table(2:,:,ip,:) = reb_tbl
       deallocate(reb_tbl)

       if(ip==1)then
          SED_lgZ0 = log10(SED_Zeds(1))                  ! Interpolation intervals
          SED_lgA0 = log10(rebAges(1))
          allocate(SED_ages(SED_nA))
          SED_ages(1)=0.d0 ; SED_ages(2:)=rebAges ;    ! Must have zero initial age
       endif

       ! Integrate the cumulative luminosities:
       SED_table(:,:,ip,2)=0.d0
       do iz = 1, SED_nZ ! Loop metallicity
          tmp = trapz1( SED_ages, SED_table(:,iz,ip,1), SED_nA, SED_table(:,iz,ip,2) )
          SED_table(:,iz,ip,2) = SED_table(:,iz,ip,2) * Gyr2sec
       end do

    end do ! End photon group loop

    if(verbose) print*, 'Initialization of SED_table : done !'

    deallocate(SEDs) ; deallocate(tbl)
    deallocate(ages) ; deallocate(rebAges) 
    deallocate(zs)
    deallocate(Ls)


  END SUBROUTINE init_SED_table



  !*************************************************************************
  SUBROUTINE inp_SED_table(age, Z, nProp, same, ret)

    ! Compute SED property by interpolation from table.
    ! input/output:
    ! age   => Star population age [Gyrs]
    ! Z     => Star population metallicity [m_metals/m_tot] 
    ! nprop => Number of property to fetch 
    !          1=log(photon # intensity [# Msun-1 s-1]),
    !          2=log(cumulative photon # intensity [# Msun-1]),
    !          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse 
    ! same  => If true then assume same age and Z as used in last call. 
    !          In this case the interpolation indexes can be recycled.
    ! ret   => The interpolated values of the sed property for every photon 
    !          group
    !-------------------------------------------------------------------------
    real(kind=8), intent(in):: age, Z
    real(kind=8):: lgAge, lgZ
    integer:: nProp
    logical:: same
    real(kind=8) :: ret(:)
    integer,save:: ia, iz
    real(kind=8),save:: da, da0, da1, dz, dz0, dz1
    !-------------------------------------------------------------------------
    if(.not. same) then
       if(age.le.0d0) then
          lgAge=-4d0 
       else
          lgAge = log10(age)
       endif
       lgZ=log10(Z)
       ia = min(max(floor((lgAge-SED_lgA0)/SED_dlgA ) + 2, 1  ),  SED_nA-1 )
       da = SED_ages(ia+1)-SED_ages(ia)
       da0= min( max(   (age-SED_ages(ia)) /da,       0. ), 1.          )
       da1= min( max(  (SED_ages(ia+1)-age)/da,       0. ), 1.          )

       iz = min(max(floor((lgZ-SED_lgZ0)/SED_dlgZ ) + 1,   1  ),  SED_nZ-1 )
       dz = sed_Zeds(iz+1)-SED_Zeds(iz)
       dz0= min( max(   (Z-SED_zeds(iz)) /dz,         0. ),  1.         )
       dz1= min( max(  (SED_Zeds(iz+1)-Z)/dz,         0. ),  1.         )

       if (abs(da0+da1-1.0d0) > 1.0d-5 .or. abs(dz0+dz1-1.0d0) > 1.0d-5) then
          write(*,*) 'Screwed up the sed interpolation ... '
          write(*,*) da0+da1,dz0+dz1
          stop
       end if
    endif

    ret = da0 * dz0 * SED_table(ia+1, iz+1, :, nProp) + &
         da1 * dz0 * SED_table(ia,   iz+1, :, nProp) + &
         da0 * dz1 * SED_table(ia+1, iz,   :, nProp) + &
         da1 * dz1 * SED_table(ia,   iz,   :, nProp)

  END SUBROUTINE inp_SED_table


  !*************************************************************************
  FUNCTION getCrosssection(lambda, species)

    ! Gives an atom-photon cross-section of given species at given wavelength,
    ! as given by Hui and Gnedin (1997) for HI and He, and Abel97 for H2.
    ! lambda  => Wavelength in angstrom
    ! species => ixHI (H2), ixHII (HI), ixHeII (HeI) or ixHeIII (HeII)
    ! returns :  photoionization or photodissociation cross-section in cm^2
    !------------------------------------------------------------------------
    use module_constants,only:clight, hp, evtoerg,ionEvs
    real(kind=8)      :: lambda, getCrosssection
    integer           :: species
    real(kind=8)      :: E0=1., cs0=0., P=1., ya=1., yw=0., y0=0., y1=1.
    real(kind=8)      :: E, x, y
    real(kind=8)      :: ionEvs(7)
    !------------------------------------------------------------------------

    !ionEvs = (/ 13.598, 24.587, 54.418, 8.152, 16.346, 33.493, 45.142 /)

    E = hp * clight/(lambda*1.d-8) / evtoerg      ! photon energy in ev
    if ( E .lt. ionEvs(species) ) then            ! below ionization energy
       getCrosssection=0.
       RETURN
    endif

    if(species .eq. 1) then ! HI
       E0 = 4.298d-1 ; cs0 = 5.475d-14  ; P  = 2.963
       ya = 32.88    ; yw  = 0          ; y0 = 0         ; y1 = 0
    endif
    if(species .eq. 2) then ! HeI
       E0 = 1.361d1  ; cs0 = 9.492d-16  ; P  = 3.188
       ya = 1.469    ; yw  = 2.039      ; y0 = 0.4434    ; y1 = 2.136
    endif
    if(species .eq. 3) then ! HeII
       E0 = 1.720    ; cs0 = 1.369d-14  ; P  = 2.963
       ya = 32.88    ; yw  = 0          ; y0 = 0         ; y1 = 0
    endif
    if(species .eq. 4) then ! SiI			
       E0 = 2.317d1    ; cs0 = 2.506d-17  ; P  = 3.546 
       ya = 20.57      ; yw  = 2.837d-1   ; y0 = 1.627d-5  ; y1 = 4.207d-1
    endif
    if(species .eq. 5) then ! SiII
       E0 = 25.56      ; cs0 = 4.14d-18   ; P  = 11.91 
       ya = 13.37      ; yw  = 1.57       ; y0 = 6.634     ; y1 = 1.272d-1
    endif
    if(species .eq. 6) then ! SiIII
       E0 = 1.659d-1   ; cs0 = 5.79d-22   ; P  = 13.36 
       ya = 1.474d2    ; yw  = 8.626d-1   ; y0 = 9.613d1   ; y1 = 6.442d-1
    endif
    if(species .eq. 7) then ! SiIV
       E0 = 12.88      ; cs0 = 6.083d-18  ; P  = 3.353 
       ya = 1.356d6    ; yw  = 0          ; y0 = 0         ; y1 = 0
    endif
    x = E/E0 - y0
    y = sqrt(x**2+y1**2)

    getCrosssection = &
         cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P
  END FUNCTION getCrosssection



  !*************************************************************************
  SUBROUTINE rebin_log(xint_log, yint_log,                                 &
       data,       nx,       ny,     x,     y,     nz,           & 
       new_data, new_nx, new_ny, new_x, new_y      )

    ! Rebin the given 2d data into constant logarithmic intervals, using 
    ! linear interpolation.
    ! xint_log,  => x and y intervals in the rebinned data. If negative,
    ! yint_log      these values represent the new number of bins.
    ! data       => The 2d data to be rebinned
    ! nx,ny      => Number of points in x and y in the original data
    ! nz         => Number of values in the data
    ! x,y        => x and y values for the data
    ! new_data  <=  The rebinned 2d data
    ! new_nx,ny <=  Number of points in x and y in the rebinned data
    ! new_x,y   <=  x and y point values for the rebinned data
    !-------------------------------------------------------------------------
    real(kind=8):: xint_log, yint_log
    integer,intent(in):: nx, ny, nz
    integer::new_nx, new_ny
    real(kind=8),intent(in):: x(nx),y(ny)
    real(kind=8),intent(in):: data(nx,ny,nz)
    real(kind=8),dimension(:,:,:),allocatable:: new_data
    real(kind=8),dimension(:),allocatable:: new_x, new_lgx, new_y, new_lgy
    real(kind=8):: dx0, dx1, dy0, dy1, x_step, y_step
    real(kind=8):: x0lg, x1lg, y0lg, y1lg
    integer :: i, j, ix, iy, ix1, iy1
    !-------------------------------------------------------------------------
    if(allocated(new_x)) deallocate(new_x)
    if(allocated(new_y)) deallocate(new_y)
    if(allocated(new_data)) deallocate(new_data)

    ! Find dimensions of the new_data and initialize it
    x0lg = log10(x(1));   x1lg = log10(x(nx)) 
    y0lg = log10(y(1));   y1lg = log10(y(ny))

    if(xint_log .lt. 0 .and. nx .gt. 1) then
       new_nx=-xint_log                             ! xint represents wanted
       xint_log = (x1lg-x0lg)/(new_nx-1)            !     number of new bins
    else
       new_nx = (x1lg-x0lg)/xint_log + 1
    endif
    allocate(new_x(new_nx)) ; allocate(new_lgx(new_nx))
    do i = 0, new_nx-1                              !  initialize the x-axis
       new_lgx(i+1) = x0lg + i*xint_log 
    end do
    new_x=10.d0**new_lgx

    if(yint_log .lt. 0 .and. ny .gt. 1) then        ! yint represents wanted
       new_ny=-yint_log                             !     number of new bins
       yint_log = (y1lg-y0lg)/(new_ny-1)
    else
       new_ny = (y1lg-y0lg)/yint_log + 1
    endif
    allocate(new_y(new_ny)) ; allocate(new_lgy(new_ny))
    do j = 0, new_ny-1                              !      ...and the y-axis
       new_lgy(j+1) = y0lg + j*yint_log
    end do
    new_y=10.d0**new_lgy

    ! Initialize new_data and find values for each point in it
    allocate(new_data(new_nx, new_ny, nz))
    do j = 1, new_ny  
       call locate(y, ny, new_y(j), iy) 
       ! y(iy) <= new_y(j) <= y(iy+1)
       ! iy is lower bound, so it can be zero but not larger than ny
       if(iy < 1) iy=1
       if (iy < ny) then
          iy1  = iy + 1
          y_step = y(iy1) - y(iy)
          dy0  = max(new_y(j) - y(iy),    0.0d0)  / y_step
          dy1  = min(y(iy1)   - new_y(j), y_step) / y_step
       else
          iy1  = iy
          dy0  = 0.0d0 ;  dy1  = 1.0d0
       end if

       do i = 1, new_nx
          call locate(x, nx, new_x(i), ix)
          if(ix < 1) ix=1
          if (ix < nx) then
             ix1  = ix+1
             x_step = x(ix1)-x(ix)
             dx0  = max(new_x(i) - x(ix),    0.0d0)  / x_step
             dx1  = min(x(ix1)   - new_x(i), x_step) / x_step
          else
             ix1  = ix
             dx0  = 0.0d0 ;  dx1  = 1.0d0
          end if

          if (abs(dx0+dx1-1.0d0) .gt. 1.0d-6 .or.                          &
               abs(dy0+dy1-1.0d0) > 1.0d-6) then
             write(*,*) 'Screwed up the rebin interpolation ... '
             write(*,*) dx0+dx1,dy0+dy1
             stop
          end if

          new_data(i,j,:) =                                                &
               dx0 * dy0 * data(ix1,iy1,:) + dx1 * dy0 * data(ix, iy1,:) + &
               dx0 * dy1 * data(ix1,iy, :) + dx1 * dy1 * data(ix, iy, :)
       end do
    end do

    deallocate(new_lgx)
    deallocate(new_lgy)

  END SUBROUTINE rebin_log


  SUBROUTINE locate(xx,n,x,j)
    ! Locates position j of a value x in an ordered array xx of n elements
    ! After: xx(j) <= x <= xx(j+1) (assuming increasing order)
    ! j is lower bound, so it can be zero but not larger than n
    !-------------------------------------------------------------------------
    !use amr_commons,only:kind=8
    integer ::  n,j,jl,ju,jm
    real(kind=8)::  xx(n),x
    !-------------------------------------------------------------------------
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




  FUNCTION getSEDLuminosity(X, Y, N, e0, e1)

    ! Compute and return luminosity in energy interval (e0,e1) [eV]   
    ! in SED Y(X). Assumes X is in Angstroms and Y in Lo/Angstroms/Msun. 
    ! (Lo=[Lo_sun], Lo_sun=[erg s-1]. total solar luminosity is 
    ! Lo_sun=10^33.58 erg/s)
    ! returns: Photon luminosity in, [# s-1 Msun-1], 
    !                             or [eV s-1 Msun-1] if SED_isEgy=true
    !-------------------------------------------------------------------------
    use module_constants,only:clight, hp, evtoerg
    !use spectrum_integrator_module
    real(kind=8):: getSEDLuminosity, X(n), Y(n), e0, e1
    integer :: N, species
    real(kind=8),parameter :: const=1.0e-8/hp/clight
    real(kind=8),parameter :: Lsun=3.8256d33 ! Solar luminosity [erg/s]
    ! const is a div by ph energy => ph count.  1e-8 is a conversion into 
    ! cgs, since wly=[angstrom] h=[erg s-1], c=[cm s-1]
    !-------------------------------------------------------------------------
    species          = 1                   ! irrelevant but must be included
    !if(.not. SED_isEgy) then               !  Photon number per sec per Msun
    getSEDLuminosity = const &
         * integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
    getSEDLuminosity = getSEDLuminosity*Lsun  ! Scale by solar luminosity
    !else                             ! SED_isEgy=true -> eV per sec per Msun
    !  getSEDLuminosity = integrateSpectrum(X, Y, N, e0, e1, species, f1)
    ! Scale by solar lum and convert to eV (bc group energies are in eV)
    ! getSEDLuminosity = getSEDLuminosity/evtoerg*Lsun 
    ! endif
  END FUNCTION getSEDLuminosity


  !*************************************************************************
  FUNCTION getSEDEgy(X, Y, N, e0, e1)

    ! Compute average energy, in eV, in energy interval (e0,e1) [eV] in SED
    ! Y(X). Assumes X is in Angstroms and Y is energy weight per angstrom
    ! (not photon count).
    !-------------------------------------------------------------------------
    use module_constants,only:clight, hp, evtoerg
    use spectrum_integrator_module
    real(dp):: getSEDEgy, X(N), Y(N), e0, e1, norm
    integer :: N,species
    real(dp),parameter :: const=1.d8*hp*clight/evtoerg! energy conversion
    !-------------------------------------------------------------------------
    species      = 1                       ! irrelevant but must be included
    norm         = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
    getSEDEgy    = const * &
         integrateSpectrum(X, Y, N, e0, e1, species, f1) / norm
  END FUNCTION getSEDEgy


  !*************************************************************************
  FUNCTION getSEDcsn(X, Y, N, e0, e1, species)

    ! Compute and return average photoionization
    ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
    ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
    ! angstrom (not photon #).
    ! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
    !-------------------------------------------------------------------------
    use spectrum_integrator_module
    use module_constants,only:ionEvs
    real(kind=8):: getSEDcsn, X(N), Y(N), e0, e1, norm
    integer :: N, species
    !-------------------------------------------------------------------------
    if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
       getSEDcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
    endif
    norm     = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
    getSEDcsn= integrateSpectrum(X, Y, N, e0, e1, species, fSigLambda)/norm
  END FUNCTION getSEDcsn


  !************************************************************************
  FUNCTION getSEDcse(X, Y, N, e0, e1, species)

    ! Compute and return average energy weighted photoionization
    ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
    ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
    ! angstrom (not photon #).
    ! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
    !-------------------------------------------------------------------------
    use spectrum_integrator_module
    use module_constants,only:ionEvs
    real(dp):: getSEDcse, X(N), Y(N), e0, e1, norm
    integer :: N, species
    !-------------------------------------------------------------------------
    if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
       getSEDcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
    endif
    norm      = integrateSpectrum(X, Y, N, e0, e1, species, f1)
    getSEDcse = integrateSpectrum(X, Y, N, e0, e1, species, fSig) / norm
  END FUNCTION getSEDcse



  FUNCTION integrateSpectrum(X, Y, N, e0, e1, species, func, doPrint)

    ! Integrate spectral weighted function in energy interval [e0,e1]
    ! X      => Wavelengths [angstrom]
    ! Y      => Spectral luminosity per angstrom at wavelenghts [XX A-1] 
    ! N      => Length of X and Y
    ! e0,e1  => Integrated interval [ev]
    ! species=> ion species, used as an argument in fx
    ! func   => Function which is integrated (of X, Y, species)
    !-------------------------------------------------------------------------
    use module_constants,only:clight,evtoerg, hp
    real(kind=8):: integrateSpectrum, X(N), Y(N), e0, e1
    integer :: N, species
    interface
       real(kind=8) function func(wavelength,intensity,species)
         ! use amr_parameters,only:kind=8
         real(kind=8)::wavelength,intensity
         integer::species
       end function func
    end interface!----------------------------------------------------------
    real(kind=8),dimension(:),allocatable:: xx, yy, f
    real(kind=8):: la0, la1
    integer :: i
    logical,optional::doPrint
    !-------------------------------------------------------------------------
    integrateSpectrum=0.
    if(N .le. 2) RETURN
    ! Convert energy interval to wavelength interval
    la0 = X(1) ; la1 = X(N)
    if(e1.gt.0) la0 = max(la0, 1.d8 * hp * clight / e1 / evtoerg)                
    if(e0.gt.0) la1 = min(la1, 1.d8 * hp * clight / e0 / evtoerg)
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

    integrateSpectrum = trapz1(xx,f,i)
    deallocate(xx) ; deallocate(yy) ; deallocate(f)

  END FUNCTION integrateSpectrum



  FUNCTION fLambda(lambda, f, species)
    real(kind=8):: fLambda, lambda, f
    integer :: species
    fLambda = f * lambda
  END FUNCTION fLambda

  FUNCTION f1(lambda, f, species)
    real(kind=8):: f1, lambda, f
    integer :: species
    f1 = f
  END FUNCTION f1

  FUNCTION fdivLambda(lambda, f, species)
    real(kind=8):: fdivlambda, lambda, f
    integer :: species
    fdivLambda = f / lambda
  END FUNCTION fdivLambda

  FUNCTION fSig(lambda, f, species)
    real(kind=8):: fSig, lambda, f
    integer :: species
    fSig = f * getCrosssection(lambda,species)
  END FUNCTION fSig

  FUNCTION fSigLambda(lambda, f, species)
    real(kind=8):: fSigLambda, lambda, f
    integer :: species
    fSigLambda = f * lambda * getCrosssection(lambda,species)
  END FUNCTION fSigLambda

  FUNCTION fSigdivLambda(lambda, f, species)
    real(kind=8):: fSigdivLambda, lambda, f
    integer :: species
    fSigdivLambda = f / lambda * getCrosssection(lambda,species)
  END FUNCTION fSigdivLambda


  !*************************************************************************
  FUNCTION trapz1(X,Y,N,cum)

    ! Integrates function Y(X) along the whole interval 1..N, using a very 
    ! simple staircase method and returns the result.
    ! Optionally, the culumative integral is returned in the cum argument.
    !-------------------------------------------------------------------------
    integer :: N,i
    real(kind=8):: trapz1
    real(kind=8):: X(N),Y(N)
    real(kind=8),optional::cum(N)
    real(kind=8),allocatable::cumInt(:)
    !-------------------------------------------------------------------------
    allocate(cumInt(N))
    cumInt(:)=0.d0
    if (N.le.1) RETURN
    do i=2,N
       cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2.d0
    end do
    trapz1 = cumInt(N)
    if(present(cum)) cum=cumInt
    deallocate(cumInt)
  END FUNCTION trapz1

  !*************************************************************************




  subroutine read_PhotonsFromStars_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present,ok

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:18) == '[PhotonsFromStars]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case('sed_dir')
             write(sed_dir,'(a)') trim(value)
          case ('outputfile')
             write(outputfile,'(a)') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case('nSEDgroups')
             read(value,*) nSEDgroups
          allocate(group_decomp(nSEDgroups+1))
          case('group_decomp')
             read(value,*) group_decomp
          case ('star_dom_type')
             write(star_dom_type,'(a)') trim(value)
          case ('star_dom_pos')
             read(value,*) star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
          case ('star_dom_rsp')
             read(value,*) star_dom_rsp
          case ('star_dom_size')
             read(value,*) star_dom_size
          case ('star_dom_rin')
             read(value,*) star_dom_rin
          case ('star_dom_rout')
             read(value,*) star_dom_rout
          case ('star_dom_thickness')
             read(value,*) star_dom_thickness
          case ('verbose')
             read(value,*) verbose
          case default
             write(*,'(a,a,a)') '> WARNING: parameter ',trim(name),' unknown '
          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)


    return

  end subroutine read_PhotonsFromStars_params


  subroutine print_PhotonsFromStars_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[PhotonsFromStars]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  sed_dir         = ',trim(sed_dir)
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a,a)')           '  repository      = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum         = ',snapnum
       write(unit,'(a,i5)')          '  nSEDgroups      = ',nSEDgroups
       write(unit,2000)              '  group_decomp    = ',group_decomp
       write(unit,'(a)')             '# computational domain parameters'
       write(unit,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(unit,'(a,3(ES10.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       select case (trim(star_dom_type))
       case ('sphere')
          write(unit,'(a,ES10.3)')       '  star_dom_rsp       = ',star_dom_rsp
       case ('shell')
          write(unit,'(a,ES10.3)')       '  star_dom_rin       = ',star_dom_rin
          write(unit,'(a,ES10.3)')       '  star_dom_rout      = ',star_dom_rout
       case('cube')
          write(unit,'(a,ES10.3)')       '  star_dom_size      = ',star_dom_size
       case('slab')
          write(unit,'(a,ES10.3)')       '  star_dom_thickness = ',star_dom_thickness
       end select
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
    else
       write(*,'(a,a,a)')         '[PhotonsFromStars]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  sed_dir         = ',trim(sed_dir)
       write(*,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(*,'(a,a)')           '  repository      = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum         = ',snapnum
       write(*,'(a,i5)')          '  nSEDgroups      = ',nSEDgroups
       write(*,2000)              '  group_decomp    = ',group_decomp
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(*,'(a,3(ES10.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       select case (trim(star_dom_type))
       case ('sphere')
          write(*,'(a,ES10.3)')       '  star_dom_rsp       = ',star_dom_rsp
       case ('shell')
          write(*,'(a,ES10.3)')       '  star_dom_rin       = ',star_dom_rin
          write(*,'(a,ES10.3)')       '  star_dom_rout      = ',star_dom_rout
       case('cube')
          write(*,'(a,ES10.3)')       '  star_dom_size      = ',star_dom_size
       case('slab')
          write(*,'(a,ES10.3)')       '  star_dom_thickness = ',star_dom_thickness
       end select
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
    end if
    2000 format (a,1000(ES10.3,1x))

    return

  end subroutine print_PhotonsFromStars_params


end program star_luminosities
