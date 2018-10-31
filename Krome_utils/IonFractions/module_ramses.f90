module module_ramses

  use module_constants, only : kB, mp, XH, mSi, mMg
  use module_domain

  implicit none

  private 

  ! stuff read from AMR files
  integer(kind=4)                  :: ncell,ncoarse,ngridmax
  real(kind=8),allocatable         :: xg(:,:)      ! grids position
  integer,allocatable              :: nbor(:,:)    ! neighboring father cells
  integer,allocatable              :: next(:)      ! next grid in list
  integer,allocatable              :: son(:)       ! sons grids
  integer,allocatable              :: cpu_map(:)  ! domain decomposition
  integer,allocatable              :: headl(:,:),taill(:,:),numbl(:,:),numbtot(:,:)
  integer,allocatable              :: headb(:,:),tailb(:,:),numbb(:,:)
  real(KIND=8),dimension(1:3)      :: xbound=(/0d0,0d0,0d0/)  

  ! Stop pretending this would work in 2D 
  integer(kind=4),parameter :: ndim = 3
  integer(kind=4),parameter :: twondim = 6
  integer(kind=4),parameter :: twotondim= 8 

  ! stuff read from the HYDRO files
  real(kind=8),allocatable         :: var(:,:)
  real(kind=8),allocatable         :: cell_x(:),cell_y(:),cell_z(:)
  integer(kind=4),allocatable      :: cell_level(:)

  integer(kind=4)                  :: ncpu

  ! conversion factors (units)
  logical                        :: conversion_scales_are_known = .False. 
  real(kind=8)                   :: dp_scale_l,dp_scale_d,dp_scale_t,dp_scale_T2,dp_scale_zsun,dp_scale_nh,dp_scale_v,dp_scale_m

  ! cooling-related stuff -------------------------------------------------------------
  type cooling_table
     integer(kind=4)          :: n11
     integer(kind=4)          :: n22
     real(kind=8),allocatable :: nH(:)    
     real(kind=8),allocatable :: T2(:)    
     real(kind=8),allocatable :: metal(:,:)  
     real(kind=8),allocatable :: cool(:,:)  
     real(kind=8),allocatable :: heat(:,:)  
     real(kind=8),allocatable :: cool_com(:,:)  
     real(kind=8),allocatable :: heat_com(:,:)  
     real(kind=8),allocatable :: cool_com_prime(:,:)  
     real(kind=8),allocatable :: heat_com_prime(:,:)  
     real(kind=8),allocatable :: metal_prime(:,:)  
     real(kind=8),allocatable :: cool_prime(:,:)  
     real(kind=8),allocatable :: heat_prime(:,:)  
     real(kind=8),allocatable :: mu(:,:)  
     real(kind=8),allocatable :: spec(:,:,:)  ! last dimension (6) is n_e, n_HI, n_HII, n_HeI, n_HeII, n_HeIII
  end type cooling_table
  type(cooling_table) :: cooling
  type cool_interp
     integer(kind=4)  :: n_nH
     real(kind=8)     :: nH_start,nH_step
     integer(kind=4)  :: n_T2
     real(kind=8)     :: T2_start,T2_step
  end type cool_interp
  type(cool_interp)  :: cool_int
  logical            :: cooling_is_read = .False. 
  ! ----------------------------------------------------------------------------------

  ! particle-related stuff -----------------------------------------------------------
  character(30) :: ParticleFields(20)  ! array of particle fields (e.g. (/'pos','vel','mass','iord','level'/) for a DM-only run)
  ! conformal time things
  integer(kind=4),parameter             :: n_frw = 1000
  real(KIND=8),dimension(:),allocatable :: aexp_frw,hexp_frw,tau_frw,t_frw
  ! ----------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ramses] of the parameter file
  ! --------------------------------------------------------------------------
  ! ramses options (not guessable from outputs)
  !Val--
  character(20)            :: star_reading_method  = 'default'  ! "default" : like originally in rascas,  "teyssier" : like in Romain's github Ramses version, and changes minit,  "fromlist", from an external list
  character(200)           :: star_info_file       = '.'        ! Where the file with star lists is, in case star_reading_method = 'fromlist'
  real(kind=8)             :: t_sne_Myr            = 1d1        ! Time before supernovae in Romain's version of Ramses
  real(kind=8)             :: eta_sn               = 2d-1       ! In Teyssier version :  if(t > t_sne_myr), minit = m/(1 - eta_sn)
  !--laV
  logical                  :: self_shielding       = .true.     ! if true, reproduce self-shielding approx made in ramses to compute nHI. 
  logical                  :: ramses_rt            = .false.    ! if true, read ramses-RT output and compute nHI and T accordingly.
  logical                  :: read_rt_variables    = .false.    ! if true, read RT variables (e.g. to compute heating terms)
  logical                  :: use_initial_mass     = .false.    ! if true, use initial masses of star particles instead of mass at output time
  logical                  :: cosmo                = .true.     ! if false, assume idealised simulation
  logical                  :: use_proper_time      = .false.    ! if true, use proper time instead of conformal time for cosmo runs. 
  ! miscelaneous
  logical                  :: verbose        = .false. ! display some run-time info on this module
  ! RT variable indices
  integer(kind=4) :: itemp  = 5 ! index of thermal pressure
  integer(kind=4) :: imetal = 6 ! index of metallicity 
  integer(kind=4) :: ihii   = 7 ! index of HII fraction 
  integer(kind=4) :: iheii  = 8 ! index of HeII fraction 
  integer(kind=4) :: iheiii = 9 ! index of HeIII fraction 
  ! Solar abundance ratio
  ! Si
  ! abundance_Si_mass == abundance_Si_number * 28.085
  real(kind=8),parameter    :: abundance_Si_mass = 9.1d-4 ! From Scarlata (private comm.): abundance_Si_number = 3.24d-5
  ! Mg
  ! abundance_Mg_mass == abundance_Mg_number * 24.305
  real(kind=8),parameter    :: abundance_Mg_mass = 8.24d-4 ! From Scarlata (private comm.): abundance_Mg_number = 3.39d-5
  ! --------------------------------------------------------------------------
  
  
  public :: read_leaf_cells, read_leaf_cells_omp, read_leaf_cells_in_domain
  public :: get_ngridtot, ramses_get_box_size_cm, get_cpu_list, get_cpu_list_periodic, get_ncpu
  public :: ramses_get_velocity_cgs, ramses_get_T_nhi_cgs, ramses_get_metallicity,  ramses_get_nh_cgs, ramses_get_T_nSiII_cgs, ramses_get_T_nMgII_cgs
  public :: ramses_read_stars_in_domain
  public :: read_ramses_params, print_ramses_params, dump_ramses_info
  
  !==================================================================================
contains

  ! ----------------
  ! public functions 
  ! ----------------

  subroutine read_leaf_cells(repository, snapnum, nleaftot, nvar, &
       & xleaf, ramses_var, leaf_level)

    ! read all leaf cell from a simulation snapshot. Return standard 
    ! ramses variables through ramses_var(nvar,nleaftot) and
    ! positions (xleaf(3,nleaftot)) and levels (leaf_level).

    implicit none 

    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum
    integer(kind=4),intent(inout)             :: nleaftot, nvar
    real(kind=8),allocatable, intent(inout)   :: ramses_var(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level(:)

    logical                                   :: do_allocs
    integer(kind=4)                           :: icpu, ileaf, icell, ivar

    if(verbose)then
       print *,' '
       print *,'...reading RAMSES cells...'
       print *,' '
    endif

    nleaftot = get_nleaf(repository,snapnum)  ! sets ncpu too 
    nvar     = get_nvar(repository,snapnum)
    allocate(ramses_var(nvar,nleaftot), xleaf(nleaftot,3), leaf_level(nleaftot))
    ncpu = get_ncpu(repository,snapnum)

    if(verbose)print *,'-- read_leaf_cells --> nleaftot, nvar, ncpu =',nleaftot,nvar,ncpu

    do_allocs = .true.
    ileaf = 0
    do icpu = 1,ncpu
       call read_amr(repository,snapnum,icpu,do_allocs)
       call read_hydro(repository,snapnum,icpu,do_allocs)
       do_allocs = .false.
       ! collect leaf cells
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             ileaf = ileaf + 1
             do ivar = 1,nvar
                ramses_var(ivar,ileaf) = var(icell,ivar)
             end do
             xleaf(ileaf,1)    = cell_x(icell)
             xleaf(ileaf,2)    = cell_y(icell)
             xleaf(ileaf,3)    = cell_z(icell)
             leaf_level(ileaf) = cell_level(icell)
          end if
       end do
    end do
    call clear_amr

    return

  end subroutine read_leaf_cells


  subroutine read_leaf_cells_omp(repository, snapnum, ncpu_read, cpu_list, &
       & nleaftot, nvar, xleaf_all, ramses_var_all, leaf_level_all)
    ! read all leaf cell from a selection of cpu files in a simulation snapshot.
    ! Return standard ramses variables through ramses_var(nvar,nleaftot) and
    ! positions (xleaf(3,nleaftot)) and levels (leaf_level).
    implicit none 
    
    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum, ncpu_read
    integer(kind=4),dimension(:),allocatable,intent(in) :: cpu_list
    integer(kind=4),intent(inout)             :: nleaftot, nvar
    real(kind=8),allocatable, intent(inout)   :: ramses_var_all(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf_all(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level_all(:)
    real(kind=8),allocatable                  :: ramses_var(:,:)
    real(kind=8),allocatable                  :: xleaf(:,:)
    integer(kind=4),allocatable               :: leaf_level(:)
    integer(kind=4)                           :: k, icpu, ileaf, icell, ivar, ilast, iloop
    logical                                   :: do_allocs
    
    if(verbose)then
       print *,' '
       print *,'...reading RAMSES cells...'
       print *,' '
    endif

    nleaftot = get_nleaf(repository,snapnum)  ! sets ncpu too 
    nvar     = get_nvar(repository,snapnum)
    allocate(ramses_var_all(nvar,nleaftot), xleaf_all(nleaftot,3), leaf_level_all(nleaftot))
    !ncpu = get_ncpu(repository,snapnum)

    if(verbose) print *,'-- read_leaf_cells_omp --> nleaftot, nvar, ncpu =',nleaftot,nvar,ncpu_read

    ileaf = 0
    iloop = 0
    ilast = 1
!$OMP PARALLEL &
!$OMP DEFAULT(private) &
!$OMP SHARED(iloop, ilast, xleaf_all, leaf_level_all, ramses_var_all, repository, snapnum, nvar, nleaftot, ncpu_read, cpu_list)
    do_allocs = .true.
!$OMP DO
    do k=1,ncpu_read
    !do icpu = 1, ncpu
       icpu=cpu_list(k)
       call read_amr_hydro(repository,snapnum,icpu,&
            & son,cpu_map,var,cell_x,cell_y,cell_z,cell_level,ncell)
       
       if (do_allocs) allocate(ramses_var(nvar,ncell), xleaf(ncell,3), leaf_level(ncell))
       do_allocs = .false.
       ! collect leaf cells
       ileaf = 0
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             ileaf = ileaf + 1
             do ivar = 1,nvar
                ramses_var(ivar,ileaf) = var(icell,ivar)
             end do
             xleaf(ileaf,1)    = cell_x(icell)
             xleaf(ileaf,2)    = cell_y(icell)
             xleaf(ileaf,3)    = cell_z(icell)
             leaf_level(ileaf) = cell_level(icell)
          end if
       end do

!$OMP CRITICAL
       ! only one CRITICAL zone
       !write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
       !     ' Reading leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
       !iloop=iloop+1

       ! ileaf is now the number of leaves on local cpu
       if(ileaf .gt. 0) then
          ! save leaf cells to return arrays
          xleaf_all(ilast:ilast-1+ileaf,1:3)  = xleaf(1:ileaf,1:3)
          leaf_level_all(ilast:ilast-1+ileaf) = leaf_level(1:ileaf)
          ramses_var_all(1:nvar,ilast:ilast-1+ileaf) = ramses_var(1:nvar,1:ileaf)
       endif
       ilast=ilast+ileaf
!$OMP END CRITICAL
    end do
!$OMP END DO
    if(.not. do_allocs) deallocate(ramses_var,xleaf,leaf_level)
!$OMP END PARALLEL
    if(verbose)then
       print*,' '
       print*,'--> Nleaves read =',ilast-1
       print*,' '
    end if

    return

  end subroutine read_leaf_cells_omp


  subroutine read_leaf_cells_in_domain(repository, snapnum, selection_domain, &
       & ncpu_read, cpu_list, &
       & nleaftot_all, nvar, xleaf_all, ramses_var_all, leaf_level_all)

    ! Read leaf cells from a selection of cpu files in a simulation snapshot
    ! and select cells inside selection_domain.
    !
    ! Return standard ramses variables through ramses_var(nvar,nleaftot) and
    ! positions (xleaf(3,nleaftot)) and levels (leaf_level).
    !
    ! Slower method since it needs to read 2 times each file (one to count, one to extract)
    
    implicit none 

    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum, ncpu_read
    type(domain),intent(in)                   :: selection_domain
    integer(kind=4),intent(inout)             :: nleaftot_all, nvar
    integer(kind=4)                           :: nleaftot
    real(kind=8),allocatable                  :: ramses_var(:,:)
    real(kind=8),allocatable                  :: xleaf(:,:)
    integer(kind=4),allocatable               :: leaf_level(:)
    real(kind=8),allocatable, intent(inout)   :: ramses_var_all(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf_all(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level_all(:)

    integer(kind=4),dimension(:),allocatable,intent(in) :: cpu_list
    
    logical                                   :: do_allocs
    integer(kind=4)                           :: icpu, ileaf, icell, ivar, nleaf_in_domain, k
    integer(kind=4)                           :: ilast
    real(kind=8),dimension(3)                 :: temp
    
    if(verbose)then
       print *,' '
       print *,'...reading RAMSES cells...'
       print *,' '
    endif

    nvar = get_nvar(repository,snapnum)
    ncpu = get_ncpu(repository,snapnum)

    ! first count leaf cells in domain...
    nleaftot = 0 ; nleaf_in_domain = 0
!$OMP PARALLEL &
!$OMP REDUCTION(+:nleaftot,nleaf_in_domain) &
!$OMP DEFAULT(private) &
!$OMP SHARED(repository, snapnum, ncpu_read, cpu_list, selection_domain)
!$OMP DO
    do k=1,ncpu_read
       icpu=cpu_list(k)
       call read_amr_hydro(repository,snapnum,icpu,&
            & son,cpu_map,var,cell_x,cell_y,cell_z,cell_level,ncell)
       !call read_amr(repository,snapnum,icpu,do_allocs)
       !call read_hydro(repository,snapnum,icpu,do_allocs)
       ! collect leaf cells
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             temp(:) = (/cell_x(icell), cell_y(icell), cell_z(icell)/)
             nleaftot = nleaftot+1
             if (domain_contains_point(temp,selection_domain)) then
                nleaf_in_domain = nleaf_in_domain + 1
             end if
          end if
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    if(verbose)print *,'-- read_leaf_cells_in_domain --> nleaftot, nleaf_in_domain, nvar, ncpu =',nleaftot, nleaf_in_domain, nvar, ncpu_read
    
    allocate(ramses_var_all(nvar,nleaf_in_domain), xleaf_all(nleaf_in_domain,3), leaf_level_all(nleaf_in_domain))
    ilast = 1
!$OMP PARALLEL &
!$OMP DEFAULT(private) &
!$OMP SHARED(ilast, xleaf_all, leaf_level_all, ramses_var_all, repository, snapnum, nvar, nleaftot, ncpu_read, cpu_list, selection_domain)
    do_allocs=.true.
!$OMP DO
    do k=1,ncpu_read
       icpu=cpu_list(k)
       call read_amr_hydro(repository,snapnum,icpu,&
            & son,cpu_map,var,cell_x,cell_y,cell_z,cell_level,ncell)
       if (do_allocs) allocate(ramses_var(nvar,ncell), xleaf(ncell,3), leaf_level(ncell))
       do_allocs=.false.
       ! collect leaf cells
       ileaf = 0
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             temp(:) = (/cell_x(icell), cell_y(icell), cell_z(icell)/)
             if (domain_contains_point(temp,selection_domain)) then
                ileaf = ileaf + 1
                do ivar = 1,nvar
                   ramses_var(ivar,ileaf) = var(icell,ivar)
                end do
                xleaf(ileaf,1)    = cell_x(icell)
                xleaf(ileaf,2)    = cell_y(icell)
                xleaf(ileaf,3)    = cell_z(icell)
                leaf_level(ileaf) = cell_level(icell)
             end if
          end if
       end do
!$OMP CRITICAL
       ! only one CRITICAL zone
       !write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
       !     ' Reading leaves ',dble(iloop) / ncpu_read * 100,' % ',char(13)
       !iloop=iloop+1

       ! ileaf is now the number of leaves on local cpu
       if(ileaf .gt. 0) then
          ! save leaf cells to return arrays
          xleaf_all(ilast:ilast-1+ileaf,1:3)  = xleaf(1:ileaf,1:3)
          leaf_level_all(ilast:ilast-1+ileaf) = leaf_level(1:ileaf)
          ramses_var_all(1:nvar,ilast:ilast-1+ileaf) = ramses_var(1:nvar,1:ileaf)
       endif
       ilast=ilast+ileaf
!$OMP END CRITICAL
    end do
!$OMP END DO
    if(allocated(ramses_var)) deallocate(ramses_var,xleaf,leaf_level)
!$OMP END PARALLEL
    
    !call clear_amr
    nleaftot_all = nleaf_in_domain
    
    return

  end subroutine read_leaf_cells_in_domain


  subroutine get_cpu_list(repository, snapnum, xmin,xmax,ymin,ymax,zmin,zmax, ncpu_read, cpu_list)
    
    implicit none
    
    character(2000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    real(KIND=8),intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    
    integer(kind=4) :: i,j,ndom
    integer(kind=4) :: ilevel, lmax, lmin
    integer(kind=4) :: ncpu_read, ncpu
    integer(kind=4) :: imin,imax,jmin,jmax,kmin,kmax
    integer(kind=4) :: impi,bit_length,maxdom
    integer(kind=4),dimension(1:8):: idom,jdom,kdom,cpu_min,cpu_max
    real(KIND=8),dimension(1:8):: bounding_min,bounding_max, order_min
    real(KIND=8)::dkey,dmax
    real(KIND=8)::dx
    
    real(kind=8),dimension(:),allocatable :: bound_key
    logical,dimension(:),allocatable      :: cpu_read
    integer(kind=4),dimension(:),allocatable,intent(out)      :: cpu_list
    
    lmax = nint(get_param_real(repository,snapnum,'levelmax'))
    ncpu = get_ncpu(repository,snapnum)
    
    dmax = max(xmax-xmin,ymax-ymin,zmax-zmin)
    
    allocate(cpu_list(1:ncpu))
    allocate(bound_key(0:ncpu))
    allocate(cpu_read(1:ncpu))
    cpu_read=.false.
    cpu_list=0
    
    call read_hilbert_keys(repository,snapnum,ncpu,bound_key)
    
    do ilevel=1,lmax
       dx=0.5d0**ilevel
       if(dx.lt.dmax)exit
    end do
    lmin=ilevel
    bit_length=lmin-1
    maxdom=2**bit_length
    imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
    if(bit_length>0)then
       imin=int(xmin*dble(maxdom))
       imax=imin+1
       jmin=int(ymin*dble(maxdom))
       jmax=jmin+1
       kmin=int(zmin*dble(maxdom))
       kmax=kmin+1
    endif
    
    dkey=(dble(2**(lmax+1)/dble(maxdom)))**ndim
    ndom=1
    if(bit_length>0)ndom=8
    idom(1)=imin; idom(2)=imax
    idom(3)=imin; idom(4)=imax
    idom(5)=imin; idom(6)=imax
    idom(7)=imin; idom(8)=imax
    jdom(1)=jmin; jdom(2)=jmin
    jdom(3)=jmax; jdom(4)=jmax
    jdom(5)=jmin; jdom(6)=jmin
    jdom(7)=jmax; jdom(8)=jmax
    kdom(1)=kmin; kdom(2)=kmin
    kdom(3)=kmin; kdom(4)=kmin
    kdom(5)=kmax; kdom(6)=kmax
    kdom(7)=kmax; kdom(8)=kmax
    
    do i=1,ndom
       if(bit_length>0)then
          call hilbert3d(idom(i),jdom(i),kdom(i),order_min(i),bit_length,1)
       else
          order_min(i)=0.0d0
       endif
       bounding_min(i)=(order_min(i))*dkey
       bounding_max(i)=(order_min(i)+1.0D0)*dkey
    end do
    
    cpu_min=0; cpu_max=0
    do impi=1,ncpu
       do i=1,ndom
          if (   bound_key(impi-1).le.bounding_min(i).and.&
               & bound_key(impi  ).gt.bounding_min(i))then
             cpu_min(i)=impi
          endif
          if (   bound_key(impi-1).lt.bounding_max(i).and.&
               & bound_key(impi  ).ge.bounding_max(i))then
             cpu_max(i)=impi
          endif
       end do
    end do
    
    ncpu_read=0
    do i=1,ndom
       do j=cpu_min(i),cpu_max(i)
          if(.not. cpu_read(j))then
             ncpu_read=ncpu_read+1
             cpu_list(ncpu_read)=j
             cpu_read(j)=.true.
          endif
       enddo
    enddo
    
    deallocate(bound_key,cpu_read)
    
    print*,'--> Ncpu to read = ',ncpu_read
    
    return
    
  end subroutine get_cpu_list
  
  
  subroutine get_cpu_list_periodic(repository, snapnum, xmin, xmax, ymin, &
       & ymax, zmin, zmax, ncpu_read, cpu_list)
    
    implicit none
    
    character(2000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    real(KIND=8),intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    
    ! Can have dom_min<0 and dom_max>1
    real(KIND=8),dimension(3)::     dom_min,dom_max
    ! Can NOT have dom_pmin<0 and dom_pmax>1
    ! (therefore extra domains across boundaires)
    real(KIND=8),dimension(3,1:2):: dom_pmin,dom_pmax
    integer(kind=4) :: i,j,ix,iy,iz,ndom
    integer(kind=4) :: ilevel, lmax, lmin
    integer(kind=4) :: ncpu_read, ncpu
    integer(kind=4) :: imin,imax,jmin,jmax,kmin,kmax
    integer(kind=4) :: impi,bit_length,maxdom
    integer(kind=4),dimension(1:8):: idom,jdom,kdom,cpu_min,cpu_max
    real(KIND=8),dimension(1:8):: bounding_min,bounding_max, order_min
    real(KIND=8)::dkey,dmin,dmax
    real(KIND=8)::dx
    
    real(kind=8),dimension(:),allocatable :: bound_key
    logical,dimension(:),allocatable      :: cpu_read
    integer(kind=4),dimension(:),allocatable,intent(out)      :: cpu_list
    
    lmax = nint(get_param_real(repository,snapnum,'levelmax'))
    ncpu = get_ncpu(repository,snapnum)
    
    allocate(cpu_list(1:ncpu))
    allocate(bound_key(0:ncpu))
    allocate(cpu_read(1:ncpu))
    cpu_read=.false.
    cpu_list=0
    
    print*,'...getting CPU list...'
    
    call read_hilbert_keys(repository,snapnum,ncpu,bound_key)
    
    ! Set up the periodic domains
    dom_min(1) = xmin;   dom_max(1) = xmax
    dom_min(2) = ymin;   dom_max(2) = ymax
    dom_min(3) = zmin;   dom_max(3) = zmax
    dom_pmin(:,1) = dom_min ! Main domain ...
    dom_pmax(:,1) = dom_max ! ...always used
    dom_pmin(:,2) = 1.
    dom_pmax(:,2) = 0.
    do i=1,3 ! loop x,y,z
       if (dom_min(i) < 0.) then ! Lower limit left of boundary
          dom_pmin(i,1) = 0.     ! Adjust main domain
          dom_pmin(i,2) = 1. + dom_min(i) ! Extra domain
          dom_pmax(i,2) = 1.
       endif
       if (dom_max(i) > 1.) then ! Upper limit right of boundary
          dom_pmax(i,1) = 1.     ! Adjust main domain
          dom_pmax(i,2) = dom_min(i) - 1. ! Extra domain
          dom_pmin(i,2) = 0.
       endif
    end do
    
    ncpu_read=0
    do ix=1,2
       do iy=1,2
          do iz=1,2
             
             dmin = min(dom_pmax(1,ix)-dom_pmin(1,ix)  &
                  ,dom_pmax(2,iy)-dom_pmin(2,iy)  &
                  ,dom_pmax(3,iz)-dom_pmin(3,iz) )
             if(dmin<=0.) cycle
             !print*,'Extracting from hilbert domain '
             !print*,xmin,xmax
             !print*,ymin,ymax
             !print*,zmin,zmax
             !print*,dom_pmin(1,ix),dom_pmax(1,ix)
             !print*,dom_pmin(2,iy),dom_pmax(2,iy)
             !print*,dom_pmin(3,iz),dom_pmax(3,iz)
             !dmax = max(xmax-xmin,ymax-ymin,zmax-zmin)
             dmax = max(dom_pmax(1,ix)-dom_pmin(1,ix)  &
                  ,dom_pmax(2,iy)-dom_pmin(2,iy)  &
                  ,dom_pmax(3,iz)-dom_pmin(3,iz) )
             
             ! Set up at most 8 sub-domains in periodic box
             do ilevel=1,lmax
                dx=0.5d0**ilevel
                if(dx.lt.dmax)exit
             end do
             lmin=ilevel
             bit_length=lmin-1
             maxdom=2**bit_length
             imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
             if(bit_length>0)then
                imin=int(dom_pmin(1,ix)*dble(maxdom))
                imax=imin+1
                jmin=int(dom_pmin(2,iy)*dble(maxdom))
                jmax=jmin+1
                kmin=int(dom_pmin(3,iz)*dble(maxdom))
                kmax=kmin+1
             endif
             
             dkey=(dble(2**(lmax+1)/dble(maxdom)))**ndim
             ndom=1
             if(bit_length>0)ndom=8
             idom(1)=imin; idom(2)=imax
             idom(3)=imin; idom(4)=imax
             idom(5)=imin; idom(6)=imax
             idom(7)=imin; idom(8)=imax
             jdom(1)=jmin; jdom(2)=jmin
             jdom(3)=jmax; jdom(4)=jmax
             jdom(5)=jmin; jdom(6)=jmin
             jdom(7)=jmax; jdom(8)=jmax
             kdom(1)=kmin; kdom(2)=kmin
             kdom(3)=kmin; kdom(4)=kmin
             kdom(5)=kmax; kdom(6)=kmax
             kdom(7)=kmax; kdom(8)=kmax
             
             do i=1,ndom
                if(bit_length>0)then
                   call hilbert3d(idom(i),jdom(i),kdom(i),order_min(i),bit_length,1)
                else
                   order_min(i)=0.0d0
                endif
                bounding_min(i)=(order_min(i))*dkey
                bounding_max(i)=(order_min(i)+1.0D0)*dkey
             end do
             
             cpu_min=0; cpu_max=0
             do impi=1,ncpu
                do i=1,ndom
                   if (   bound_key(impi-1).le.bounding_min(i).and.&
                        & bound_key(impi  ).gt.bounding_min(i))then
                      cpu_min(i)=impi
                   endif
                   if (   bound_key(impi-1).lt.bounding_max(i).and.&
                        & bound_key(impi  ).ge.bounding_max(i))then
                      cpu_max(i)=impi
                   endif
                end do
             end do
             
             do i=1,ndom
                do j=cpu_min(i),cpu_max(i)
                   if(.not. cpu_read(j))then
                      ncpu_read=ncpu_read+1
                      cpu_list(ncpu_read)=j
                      cpu_read(j)=.true.
                   endif
                enddo
             enddo
             
          end do !ix=1,2
       end do !iy=1,2
    end do !iz=1,2
    
    deallocate(bound_key,cpu_read)
    
    print*,'--> Ncpu to read = ',ncpu_read
    
    return
    
  end subroutine get_cpu_list_periodic
  
  
  subroutine read_hilbert_keys(repository,snapnum,ncpu,bound_key)
    
    implicit none
    
    character(2000),intent(in)                   :: repository
    integer(kind=4),intent(in)                   :: snapnum, ncpu
    real(kind=8),dimension(0:ncpu),intent(inout) :: bound_key
    
    logical(kind=4)            :: not_ok
    character(512)             :: nomfich
    character(512)             :: line,name,value,orderingtype
    integer(kind=4)            :: i, impi
    integer(kind=4),parameter  :: param_unit = 13
    
    not_ok = .true.
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,'/info_',snapnum,'.txt'
    open(unit=param_unit,file=nomfich,status='old',form='formatted')
    do
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       if (trim(name) .eq. 'ordering type') then
          read(value,*) orderingtype
          if (orderingtype=='hilbert') then
             ! start reading keys
             ! skip one line
             read(param_unit,*)
             do impi=1,ncpu
                read(param_unit,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
             end do
          end if
       end if
    end do
2   close (param_unit)
    
    return
    
  end subroutine read_hilbert_keys
  
  
  subroutine hilbert3d(x,y,z,order,bit_length,npoint)
    implicit none
    
    integer     ,INTENT(IN)                     ::bit_length,npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
    real(kind=8),INTENT(OUT),dimension(1:npoint)::order
    
    logical,dimension(0:3*bit_length-1)::i_bit_mask
    logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
    integer,dimension(0:7,0:1,0:11)::state_diagram
    integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit
    
    if(bit_length>bit_size(bit_length))then
       write(*,*)'Maximum bit length=',bit_size(bit_length)
       write(*,*)'stop in hilbert3d'
       stop
    endif
    
    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
         &   0, 1, 3, 2, 7, 6, 4, 5,&
         &   2, 6, 0, 7, 8, 8, 0, 7,&
         &   0, 7, 1, 6, 3, 4, 2, 5,&
         &   0, 9,10, 9, 1, 1,11,11,&
         &   0, 3, 7, 4, 1, 2, 6, 5,&
         &   6, 0, 6,11, 9, 0, 9, 8,&
         &   2, 3, 1, 0, 5, 4, 6, 7,&
         &  11,11, 0, 7, 5, 9, 0, 7,&
         &   4, 3, 5, 2, 7, 0, 6, 1,&
         &   4, 4, 8, 8, 0, 6,10, 6,&
         &   6, 5, 1, 2, 7, 4, 0, 3,&
         &   5, 7, 5, 3, 1, 1,11,11,&
         &   4, 7, 3, 0, 5, 6, 2, 1,&
         &   6, 1, 6,10, 9, 4, 9,10,&
         &   6, 7, 5, 4, 1, 0, 2, 3,&
         &  10, 3, 1, 1,10, 3, 5, 9,&
         &   2, 5, 3, 4, 1, 6, 0, 7,&
         &   4, 4, 8, 8, 2, 7, 2, 3,&
         &   2, 1, 5, 6, 3, 0, 4, 7,&
         &   7, 2,11, 2, 7, 5, 8, 5,&
         &   4, 5, 7, 6, 3, 2, 0, 1,&
         &  10, 3, 2, 6,10, 3, 4, 4,&
         &   6, 1, 7, 0, 5, 2, 4, 3 /), &
         & (/8 ,2, 12 /) )
    
    do ip=1,npoint
       
       ! convert to binary
       do i=0,bit_length-1
          x_bit_mask(i)=btest(x(ip),i)
          y_bit_mask(i)=btest(y(ip),i)
          z_bit_mask(i)=btest(z(ip),i)
       enddo
       
       ! interleave bits
       do i=0,bit_length-1
          i_bit_mask(3*i+2)=x_bit_mask(i)
          i_bit_mask(3*i+1)=y_bit_mask(i)
          i_bit_mask(3*i  )=z_bit_mask(i)
       end do
       
       ! build Hilbert ordering using state diagram
       cstate=0
       do i=bit_length-1,0,-1
          b2=0 ; if(i_bit_mask(3*i+2))b2=1
          b1=0 ; if(i_bit_mask(3*i+1))b1=1
          b0=0 ; if(i_bit_mask(3*i  ))b0=1
          sdigit=b2*4+b1*2+b0
          nstate=state_diagram(sdigit,0,cstate)
          hdigit=state_diagram(sdigit,1,cstate)
          i_bit_mask(3*i+2)=btest(hdigit,2)
          i_bit_mask(3*i+1)=btest(hdigit,1)
          i_bit_mask(3*i  )=btest(hdigit,0)
          cstate=nstate
       enddo
       
       ! save Hilbert key as double precision real
       order(ip)=0.
       do i=0,3*bit_length-1
          b0=0 ; if(i_bit_mask(i))b0=1
          order(ip)=order(ip)+dble(b0)*dble(2)**i
       end do
       
    end do
    
  end subroutine hilbert3d
  
  
  function get_nGridTot(repository,snapnum)

    ! get total number of grids in the simulation 

    implicit none 

    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    integer(kind=4)             :: get_nGridTot
    character(1000)             :: nomfich
    logical                     :: ok
    integer(kind=4)             :: icpu,ngrid_current

    ncpu = get_ncpu(repository,snapnum)
    get_nGridTot = 0
    do icpu = 1,ncpu
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       inquire(file=nomfich, exist=ok)
       if(.not. ok)then
          write(*,*)'File '//TRIM(nomfich)//' not found'    
          stop
       end if
       open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)ngrid_current
       close(10)
       get_nGridTot = get_nGridTot + ngrid_current
    end do

    return

  end function get_nGridTot



  subroutine ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nhi)

    implicit none 

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nhi(nleaf), temp(nleaf)
    real(kind=8),allocatable    :: boost(:)
    real(kind=8)                :: xhi
    integer(kind=4) :: ihx,ihy,i
    real(kind=8)    :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4) :: if1,if2,jf1,jf2

    real(kind=8),allocatable,dimension(:)    :: mu

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       nhi  = ramses_var(1,:) * dp_scale_nh  * (1.d0 - ramses_var(ihii,:))   ! nb of H atoms per cm^3
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu   
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2            ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard

       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if

       nhi  = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]

       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nhi(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if

    
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nhi(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! neutral H density 
          f = dxx1 * dyy1 * cooling%spec(if2,jf2,2) + dxx2 * dyy1 * cooling%spec(if1,jf2,2) &
               & + dxx1 * dyy2 * cooling%spec(if2,jf1,2) + dxx2 * dyy2 * cooling%spec(if1,jf1,2)
          xhi = 10.0d0**(f-xx)  ! this is xHI = nHI/nH where the n's include the boost. 
          nhi(i) = xhi * nhi(i) ! nHI (cm^-3) (boost-free)
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do

       deallocate(boost)

    endif
    
    return

  end subroutine ramses_get_T_nhi_cgs

  subroutine ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,velocity_cgs)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: velocity_cgs(3,nleaf)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    velocity_cgs = ramses_var(2:4,:) * dp_scale_v ! [ cm / s ]

    return

  end subroutine ramses_get_velocity_cgs

  subroutine ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nh)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nh(nleaf)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    nh = ramses_var(1,:) * dp_scale_nh ! [ H / cm^3 ]

    return

  end subroutine ramses_get_nh_cgs


  subroutine ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)

    implicit none

    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: metallicity(nleaf)

    if (nvar < 6) then
       print*,'No metals !!! '
       stop
    end if
    metallicity = ramses_var(imetal,:) 

    return

  end subroutine ramses_get_metallicity

  !Val--
  subroutine ramses_get_he_fractions(nleaf,nvar,ramses_var,xheii,xheiii)

    implicit none

    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: xheii(nleaf), xheiii(nleaf)

    xheii = ramses_var(iheii,:)
    xheiii = ramses_var(iheiii,:) 

    return

  end subroutine ramses_get_he_fractions
  !--laV

  
  subroutine ramses_get_T_nSiII_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nSiII)

    implicit none 

    character(1000),intent(in)            :: repository
    integer(kind=4),intent(in)            :: snapnum
    integer(kind=4),intent(in)            :: nleaf, nvar
    real(kind=8),intent(in)               :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)            :: nSiII(nleaf), temp(nleaf)
    real(kind=8),allocatable              :: boost(:)
    integer(kind=4)                       :: ihx,ihy,i
    real(kind=8)                          :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4)                       :: if1,if2,jf1,jf2
    real(kind=8),allocatable,dimension(:) :: mu, nh

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if


    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2            ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard
       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if
       allocate(nh(nleaf))
       nh   = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]
       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nh(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nh(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do
       deallocate(boost,nh)
    endif       
    
    ! from T, we can compute the SiII fraction as 100% between 6.306571e+04 K and 1.264600e+05 K, and 0% elsewhere.
    ! (These limits correspond to the following ionisation energies:
    ! - Si-Si+   : 8.15169 ev
    ! - Si+-Si++ : 16.34585 eV
    do i = 1,nleaf
       if (temp(i) >= 6.306571d4 .and. temp(i) <= 1.2646d5) then
          nSiII(i) = ramses_var(1,i) * dp_scale_d * ramses_var(imetal,i) * dp_scale_zsun * abundance_Si_mass / mSi   ! [#/cm3]
       else
          nSiII(i) = 0.0d0
       end if
    end do
    !print*, 'test ramses_T : ', minval(temp(:)), maxval(temp(:))
    
    return

  end subroutine ramses_get_T_nSiII_cgs



  subroutine ramses_get_T_nMgII_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nMgII)

    implicit none 

    character(1000),intent(in)            :: repository
    integer(kind=4),intent(in)            :: snapnum
    integer(kind=4),intent(in)            :: nleaf, nvar
    real(kind=8),intent(in)               :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)            :: nMgII(nleaf), temp(nleaf)
    real(kind=8),allocatable              :: boost(:)
    integer(kind=4)                       :: ihx,ihy,i
    real(kind=8)                          :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4)                       :: if1,if2,jf1,jf2
    real(kind=8),allocatable,dimension(:) :: mu, nh

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2                ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard
       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if
       allocate(nh(nleaf))
       nh   = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]
       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nh(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nh(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do
       deallocate(boost,nh)
    endif       
    
    ! from T, we can compute the MgII fraction as 100% between 5.915393e+04 K and 1.163181e+05 K, and 0% elsewhere.
    ! (These limits correspond to the following ionisation energies:
    ! - Mg  - Mg+  :  7.646235 ev
    ! - Mg+ - Mg++ : 15.03527 eV
    do i = 1,nleaf
       if (temp(i) >= 5.915393d4 .and. temp(i) <= 1.163181d5) then
          nMgII(i) = ramses_var(1,i) * dp_scale_d * ramses_var(imetal,i) * dp_scale_zsun * abundance_Mg_mass / mMg   ! [#/cm3]
       else
          nMgII(i) = 0.0d0
       end if
    end do
    
    return

  end subroutine ramses_get_T_nMgII_cgs



  function ramses_get_box_size_cm(repository,snapnum)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    real(kind=8)                :: ramses_get_box_size_cm 

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    ramses_get_box_size_cm = get_param_real(repository,snapnum,'boxlen') * dp_scale_l  ! [ cm ] 

    return

  end function ramses_get_box_size_cm


  subroutine dump_ramses_info(repository,snapnum,ilun)

    integer(kind=4),intent(in) :: ilun
    integer(kind=4),intent(in) :: snapnum
    character(1000),intent(in) :: repository
    real(kind=8) :: unit_t, unit_l, unit_d, boxlen, time, aexp, H0, omega_m, omega_l, omega_k, omega_b
    
    ! get data
    unit_l  = get_param_real(repository,snapnum,'unit_l')
    unit_d  = get_param_real(repository,snapnum,'unit_d')
    unit_t  = get_param_real(repository,snapnum,'unit_t')
    boxlen  = get_param_real(repository,snapnum,'boxlen')
    time    = get_param_real(repository,snapnum,'time')
    aexp    = get_param_real(repository,snapnum,'aexp')
    H0      = get_param_real(repository,snapnum,'H0')
    omega_m = get_param_real(repository,snapnum,'omega_m')
    omega_l = get_param_real(repository,snapnum,'omega_l')
    omega_k = get_param_real(repository,snapnum,'omega_k')
    omega_b = get_param_real(repository,snapnum,'omega_b')

    ! dump data in open file using unit
    ! Write physical parameters as in RAMSES
    write(ilun,'(" boxlen      =",E23.15)')boxlen
    write(ilun,'(" time        =",E23.15)')time
    write(ilun,'(" aexp        =",E23.15)')aexp
    write(ilun,'(" H0          =",E23.15)')H0
    write(ilun,'(" omega_m     =",E23.15)')omega_m
    write(ilun,'(" omega_l     =",E23.15)')omega_l
    write(ilun,'(" omega_k     =",E23.15)')omega_k
    write(ilun,'(" omega_b     =",E23.15)')omega_b
    write(ilun,'(" unit_l      =",E23.15)')unit_l
    write(ilun,'(" unit_d      =",E23.15)')unit_d
    write(ilun,'(" unit_t      =",E23.15)')unit_t
    write(ilun,*)
  
  end subroutine dump_ramses_info


  !==================================================================================
  ! ----------------
  ! private functions 
  ! ----------------

  subroutine read_hydro(repository,snapnum,icpu,do_allocs)

    implicit none

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository
    logical,intent(in)          :: do_allocs
    real(kind=8)                :: dx
    character(1000)             :: nomfich
    integer(kind=4)             :: i,nlevelmax,nboundary,ix,iy,iz,ind,ilevel,ibound,ncache,istart,ivar,iskip,igrid,nvarH,nvarRT
    real(kind=8),allocatable    :: xc(:,:),xx(:)
    integer(kind=4),allocatable :: ind_grid(:)

    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    read(10)
    read(10)nvarH
    read(10)
    read(10)nlevelmax
    read(10)nboundary
    read(10)
    
    if (read_rt_variables) then
       ! Open RT file and get nvarRT
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
       open(unit=12,file=nomfich,status='old',form='unformatted')
       read(12)
       read(12)nvarRT
       read(12)
       read(12)
       read(12)
       read(12)
    else
       nvarRT = 0
    end if
    
    if (do_allocs) then
       allocate(var(1:ncell,1:nvarH+nvarRT))
       allocate(cell_x(1:ncell),cell_y(1:ncell),cell_z(1:ncell))
       allocate(cell_level(1:ncell))
    end if
    allocate(xc(1:twotondim,1:ndim))

    cell_level = -1

    do ilevel=1,nlevelmax

       dx=0.5d0**ilevel
       do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc(ind,1)=(dble(ix)-0.5D0)*dx
          xc(ind,2)=(dble(iy)-0.5D0)*dx
          xc(ind,3)=(dble(iz)-0.5D0)*dx
       end do

       do ibound=1,nboundary+ncpu
          if(ibound<=ncpu)then  ! in the box 
             ncache=numbl(ibound,ilevel)   ! nb of grids in the simulated box. 
             istart=headl(ibound,ilevel)   ! head of grid list of simulated box
          else                  ! boundaries of simulated volume (aka useless)
             ncache=numbb(ibound-ncpu,ilevel)
             istart=headb(ibound-ncpu,ilevel)
          end if
          read(10)!ilevel2
          read(10)!numbl2
          if (read_rt_variables) then
             read(12)
             read(12)
          end if
          if(ncache>0)then
             allocate(ind_grid(1:ncache))
             allocate(xx(1:ncache))
             ! Loop over level grids
             igrid=istart
             do i=1,ncache
                ind_grid(i)=igrid
                igrid=next(igrid)
             end do
             ! Loop over cells
             do ind=1,twotondim
                iskip=ncoarse+(ind-1)*ngridmax
                ! Loop over conservative variables
                do ivar=1,nvarH
                   read(10) xx
                   if (ibound > ncpu) cycle  ! dont bother with boundaries
                   do i = 1, ncache
                      var(ind_grid(i)+iskip,ivar) = xx(i)
                   end do
                end do
                if (read_rt_variables) then 
                   do ivar=1,nvarRT
                      read(12) xx
                      if (ibound > ncpu) cycle  ! dont bother with boundaries
                      do i = 1, ncache
                         var(ind_grid(i)+iskip,ivar+nvarH) = xx(i)
                      end do
                   end do
                end if
                do i = 1,ncache
                   cell_x(ind_grid(i)+iskip) = xc(ind,1) + xg(ind_grid(i),1) -xbound(1)
                   cell_y(ind_grid(i)+iskip) = xc(ind,2) + xg(ind_grid(i),2) -xbound(2)
                   cell_z(ind_grid(i)+iskip) = xc(ind,3) + xg(ind_grid(i),3) -xbound(3)
                   cell_level(ind_grid(i)+iskip)      = ilevel
                end do
             end do
             deallocate(ind_grid,xx)
          end if
       end do
    end do
    deallocate(xc)
    close(10)
    if (read_rt_variables) close(12)
    
    return

  end subroutine read_hydro


  subroutine read_amr(repository,snapnum,icpu,do_allocs)

    implicit none 

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository
    logical,intent(in)          :: do_allocs
    character(1000)             :: nomfich 
    integer,allocatable         :: ind_grid(:),iig(:),grid(:)
    real(kind=8),allocatable    :: xxg(:)
    logical                     :: ok
    integer(kind=4)             :: i,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)             :: ilevel,ncache,ibound,idim,ind,iskip

    ! Vérification de l'existence des fichiers AMR
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    inquire(file=nomfich, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(nomfich)//' not found'    
       stop
    end if
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    ! Read grid variables
    read(10)
    read(10)
    read(10)nx,ny,nz
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

    ! Critical parameter: define the root level of the tree
    ncoarse=nx*ny*nz
    read(10)nlevelmax
    read(10)ngridmax
    read(10)nboundary
    read(10)!ngrid_current
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    if (do_allocs) allocate( headl(1:ncpu,1:nlevelmax),taill(1:ncpu,1:nlevelmax), &
         & numbl(1:ncpu,1:nlevelmax),numbtot(1:10,1:nlevelmax), &
         & headb(1:nboundary,1:nlevelmax),tailb(1:nboundary,1:nlevelmax), &
         & numbb(1:nboundary,1:nlevelmax) )
    headl=0;taill=0;numbl=0;numbtot=0;headb=0;tailb=0;numbb=0
    ! Allocate tree arrays
    if (do_allocs) then
       allocate(next(1:ngridmax))
       allocate(nbor(1:ngridmax,1:twondim))
    end if
    nbor=0; next=0
    ! Allocate grid center coordinates
    if (do_allocs) allocate(xg(1:ngridmax,1:ndim))
    xg=0.0D0
    ! Read levels variables
    read(10)headl(1:ncpu,1:nlevelmax)
    read(10)taill(1:ncpu,1:nlevelmax)
    read(10)numbl(1:ncpu,1:nlevelmax)
    read(10)numbtot(1:10,1:nlevelmax)
    ! Read boundary linked list
    if(nboundary>0)then
       read(10)headb(1:nboundary,1:nlevelmax)
       read(10)tailb(1:nboundary,1:nlevelmax)
       read(10)numbb(1:nboundary,1:nlevelmax)
    end if
    !  Read free memory
    read(10)
    next(ngridmax) = 0
    ! Read cpu boundaries
    read(10)
    read(10)
    ncell=ncoarse+twotondim*ngridmax
    if (do_allocs) allocate(son(1:ncell),cpu_map(1:ncell))
    son=0; cpu_map=0
    ! Read coarse level
    read(10)son(1:ncoarse)       
    read(10)
    read(10)cpu_map(1:ncoarse)
    do ilevel=1,nlevelmax
       do ibound=1,nboundary+ncpu
          if(ibound<=ncpu)then
             ncache=numbl(ibound,ilevel)
          else
             ncache=numbb(ibound-ncpu,ilevel)
          end if
          if(ncache>0)then
             allocate(ind_grid(1:ncache))
             allocate(xxg(1:ncache))
             allocate(iig(1:ncache))
             allocate(grid(1:ncache))
             ! Read grid index
             read(10)ind_grid
             ! Read next index
             read(10)iig
             do i=1,ncache
                next(ind_grid(i))=iig(i)
             end do
             ! Read prev index (skip)
             read(10)iig
             ! Read grid center
             do idim=1,ndim
                read(10)xxg
                do i=1,ncache
                   xg(ind_grid(i),idim)=xxg(i)
                end do
             end do
             ! Read father index (skip)
             read(10)iig
             ! Read nbor index
             do ind=1,twondim
                read(10)iig
                do i=1,ncache
                   nbor(ind_grid(i),ind)=iig(i)
                end do
             end do
             ! Read son index
             do ind=1,twotondim
                iskip=ncoarse+(ind-1)*ngridmax
                read(10)iig
                do i=1,ncache
                   son(ind_grid(i)+iskip)=iig(i)
                end do
             end do
             ! Read cpu map
             do ind=1,twotondim
                iskip=ncoarse+(ind-1)*ngridmax
                read(10)iig
                do i=1,ncache
                   cpu_map(ind_grid(i)+iskip)=iig(i)
                end do
             end do
             ! Read refinement map (skip)
             do ind=1,twotondim
                read(10)!iig 
             end do
             deallocate(xxg,iig,grid,ind_grid)
          end if
       end do
    end do
    close(10)
    return

  end subroutine read_amr


  subroutine clear_amr

    implicit none

    deallocate(son,cpu_map,xg,nbor,next)
    deallocate(headl,taill,numbl,numbtot,headb,tailb,numbb)
    deallocate(var,cell_x,cell_y,cell_z,cell_level)

    return

  end subroutine clear_amr


  subroutine read_amr_hydro(repository,snapnum,icpu,&
       & son_l,cpu_map_l,var_l,cell_x_l,cell_y_l,cell_z_l,cell_level_l,ncell_l)
    ! purpose: use only local variables for OMP

    !$ use OMP_LIB
    implicit none 

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository

    character(1000)             :: nomfich 
    integer(kind=4),allocatable :: ind_grid(:),iig(:),grid(:)
    real(kind=8),allocatable    :: xxg(:)
    logical                     :: ok
    integer(kind=4)             :: i,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)             :: ilevel,ncache,ibound,idim,ind,iskip,iunit,iu2,rank

    ! stuff read from AMR files
    integer(kind=4),intent(out)      :: ncell_l
    integer(kind=4)                  :: ncoarse_l,ngridmax_l
    real(kind=8),allocatable         :: xg_l(:,:)      ! grids position
    integer,allocatable              :: nbor_l(:,:)    ! neighboring father cells
    integer,allocatable              :: next_l(:)      ! next grid in list
    integer,allocatable,intent(out)  :: son_l(:)       ! sons grids
    integer,allocatable,intent(out)  :: cpu_map_l(:)   ! domain decomposition
    integer,allocatable              :: headl_l(:,:),taill_l(:,:),numbl_l(:,:),numbtot_l(:,:)
    integer,allocatable              :: headb_l(:,:),tailb_l(:,:),numbb_l(:,:)
    real(KIND=8),dimension(1:3)      :: xbound_l=(/0d0,0d0,0d0/)  

    real(kind=8)                :: dx
    integer(kind=4)             :: ix,iy,iz,istart,ivar,igrid,nvarH,nvarRT
    real(kind=8),allocatable    :: xc(:,:),xx(:)

    ! stuff read from the HYDRO files
    real(kind=8),allocatable,intent(out)     :: var_l(:,:)
    real(kind=8),allocatable,intent(out)     :: cell_x_l(:),cell_y_l(:),cell_z_l(:)
    integer(kind=4),allocatable,intent(out)  :: cell_level_l(:)
    
    rank = 1
    !$ rank = OMP_GET_THREAD_NUM()
    iunit=10+rank*2
    iu2 = 10+rank*2+1

    ! Vérification de l'existence des fichiers AMR
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    inquire(file=nomfich, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(nomfich)//' not found'    
       stop
    end if
    open(unit=iunit,file=nomfich,form='unformatted',status='old',action='read')
    ! Read grid variables
    read(iunit)
    read(iunit)
    read(iunit)nx,ny,nz
    xbound_l=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

    ! Critical parameter: define the root level of the tree
    ncoarse_l=nx*ny*nz
    read(iunit)nlevelmax
    read(iunit)ngridmax_l
    read(iunit)nboundary
    read(iunit)!ngrid_current
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    allocate( headl_l(1:ncpu,1:nlevelmax),taill_l(1:ncpu,1:nlevelmax), &
         & numbl_l(1:ncpu,1:nlevelmax),numbtot_l(1:10,1:nlevelmax), &
         & headb_l(1:nboundary,1:nlevelmax),tailb_l(1:nboundary,1:nlevelmax), &
         & numbb_l(1:nboundary,1:nlevelmax) )
    headl_l=0;taill_l=0;numbl_l=0;numbtot_l=0;headb_l=0;tailb_l=0;numbb_l=0
    ! Allocate tree arrays
    allocate(next_l(1:ngridmax_l))
    allocate(nbor_l(1:ngridmax_l,1:twondim))
    nbor_l=0; next_l=0
    ! Allocate grid center coordinates
    allocate(xg_l(1:ngridmax_l,1:ndim))
    xg_l=0.0D0
    ! Read levels variables
    read(iunit)headl_l(1:ncpu,1:nlevelmax)
    read(iunit)taill_l(1:ncpu,1:nlevelmax)
    read(iunit)numbl_l(1:ncpu,1:nlevelmax)
    read(iunit)numbtot_l(1:10,1:nlevelmax)
    ! Read boundary linked list
    if(nboundary>0)then
       read(iunit)headb_l(1:nboundary,1:nlevelmax)
       read(iunit)tailb_l(1:nboundary,1:nlevelmax)
       read(iunit)numbb_l(1:nboundary,1:nlevelmax)
    end if
    !  Read free memory
    read(iunit)
    next_l(ngridmax_l) = 0
    ! Read cpu boundaries
    read(iunit)
    read(iunit)
    ncell_l=ncoarse_l+twotondim*ngridmax_l
    allocate(son_l(1:ncell_l),cpu_map_l(1:ncell_l))
    son_l=0; cpu_map_l=0
    ! Read coarse level
    read(iunit)son_l(1:ncoarse_l)       
    read(iunit)
    read(iunit)cpu_map_l(1:ncoarse_l)
    do ilevel=1,nlevelmax
       do ibound=1,nboundary+ncpu
          if(ibound<=ncpu)then
             ncache=numbl_l(ibound,ilevel)
          else
             ncache=numbb_l(ibound-ncpu,ilevel)
          end if
          if(ncache>0)then
             allocate(ind_grid(1:ncache))
             allocate(xxg(1:ncache))
             allocate(iig(1:ncache))
             allocate(grid(1:ncache))
             ! Read grid index
             read(iunit)ind_grid
             ! Read next index
             read(iunit)iig
             do i=1,ncache
                next_l(ind_grid(i))=iig(i)
             end do
             ! Read prev index (skip)
             read(iunit)iig
             ! Read grid center
             do idim=1,ndim
                read(iunit)xxg
                do i=1,ncache
                   xg_l(ind_grid(i),idim)=xxg(i)
                end do
             end do
             ! Read father index (skip)
             read(iunit)iig
             ! Read nbor index
             do ind=1,twondim
                read(iunit)iig
                do i=1,ncache
                   nbor_l(ind_grid(i),ind)=iig(i)
                end do
             end do
             ! Read son index
             do ind=1,twotondim
                iskip=ncoarse_l+(ind-1)*ngridmax_l
                read(iunit)iig
                do i=1,ncache
                   son_l(ind_grid(i)+iskip)=iig(i)
                end do
             end do
             ! Read cpu map
             do ind=1,twotondim
                iskip=ncoarse_l+(ind-1)*ngridmax_l
                read(iunit)iig
                do i=1,ncache
                   cpu_map_l(ind_grid(i)+iskip)=iig(i)
                end do
             end do
             ! Read refinement map (skip)
             do ind=1,twotondim
                read(iunit)!iig 
             end do
             deallocate(xxg,iig,grid,ind_grid)
          end if
       end do
    end do
    close(iunit)
    ! => can return cpu_map_l & son_l

    !print*,'in module_ramses.read_amr_hydro -> ',ncell_l,nvar,icpu,snapnum
    ! and then the hydro file
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=iunit,file=nomfich,form='unformatted',status='old',action='read')
    read(iunit)
    read(iunit)nvarH
    read(iunit)
    read(iunit)nlevelmax
    read(iunit)nboundary
    read(iunit)

    if (read_rt_variables) then
       ! Open RT file and get nvarRT
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
       open(unit=iu2,file=nomfich,status='old',form='unformatted')
       read(iu2)
       read(iu2)nvarRT
       read(iu2)
       read(iu2)
       read(iu2)
       read(iu2)
    else
       nvarRT = 0
    end if
    
    allocate(var_l(1:ncell_l,1:nvarH+nvarRT))
    allocate(cell_x_l(1:ncell_l),cell_y_l(1:ncell_l),cell_z_l(1:ncell_l))
    allocate(cell_level_l(1:ncell_l))
    allocate(xc(1:twotondim,1:ndim))

    cell_level_l = -1

    do ilevel=1,nlevelmax

       dx=0.5d0**ilevel
       do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc(ind,1)=(dble(ix)-0.5D0)*dx
          xc(ind,2)=(dble(iy)-0.5D0)*dx
          xc(ind,3)=(dble(iz)-0.5D0)*dx
       end do

       do ibound=1,nboundary+ncpu
          if(ibound<=ncpu)then  ! in the box 
             ncache=numbl_l(ibound,ilevel)   ! nb of grids in the simulated box. 
             istart=headl_l(ibound,ilevel)   ! head of grid list of simulated box
          else                  ! boundaries of simulated volume (aka useless)
             ncache=numbb_l(ibound-ncpu,ilevel)
             istart=headb_l(ibound-ncpu,ilevel)
          end if
          read(iunit)!ilevel2
          read(iunit)!numbl2
          if (read_rt_variables) then
             read(iu2)
             read(iu2)
          end if
          if(ncache>0)then
             allocate(ind_grid(1:ncache))
             allocate(xx(1:ncache))
             ! Loop over level grids
             igrid=istart
             do i=1,ncache
                ind_grid(i)=igrid
                igrid=next_l(igrid)
             end do
             ! Loop over cells
             do ind=1,twotondim
                iskip=ncoarse_l+(ind-1)*ngridmax_l
                ! Loop over conservative variables
                do ivar=1,nvarH
                   read(iunit) xx
                   if (ibound > ncpu) cycle  ! dont bother with boundaries
                   do i = 1, ncache
                      var_l(ind_grid(i)+iskip,ivar) = xx(i)
                   end do
                end do
                if (read_rt_variables) then 
                   do ivar=1,nvarRT
                      read(iu2) xx
                      if (ibound > ncpu) cycle  ! dont bother with boundaries
                      do i = 1, ncache
                         var_l(ind_grid(i)+iskip,ivar+nvarH) = xx(i)
                      end do
                   end do
                end if
                do i = 1,ncache
                   cell_x_l(ind_grid(i)+iskip) = xc(ind,1) + xg_l(ind_grid(i),1) -xbound_l(1)
                   cell_y_l(ind_grid(i)+iskip) = xc(ind,2) + xg_l(ind_grid(i),2) -xbound_l(2)
                   cell_z_l(ind_grid(i)+iskip) = xc(ind,3) + xg_l(ind_grid(i),3) -xbound_l(3)
                   cell_level_l(ind_grid(i)+iskip)      = ilevel
                end do
             end do
             deallocate(ind_grid,xx)
          end if
       end do
    end do
    deallocate(xc)
    close(iunit)
    if (read_rt_variables) close(iu2)
    ! => can return var_l, cell_x_l, cell_y_l, cell_z_l, cell_level_l

    deallocate(headl_l, taill_l, numbl_l, numbtot_l, headb_l, tailb_l, numbb_l)
    deallocate(next_l, nbor_l, xg_l) 
    
    return
    
  end subroutine read_amr_hydro


  function get_nleaf(repository,snapnum)

    implicit none

    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    integer(kind=4)             :: get_nleaf
    integer(kind=4)             :: icpu,icell
    logical                     :: do_allocs
    character(1000)             :: nomfich 
    integer,allocatable         :: ind_grid(:),iig(:)
    integer,allocatable         :: son(:)
    integer,allocatable         :: cpu_map(:)
    integer,allocatable         :: numbl(:,:)
    integer,allocatable         :: numbb(:,:)
    logical                     :: ok
    integer(kind=4)             :: i,nx,ny,nz,nlevelmax,nboundary,ncell,ncoarse,ngridmax
    integer(kind=4)             :: ilevel,ncache,ibound,idim,ind,iskip
    integer(kind=4)             :: ndim,twondim,twotondim
    integer(kind=4)             :: nleaf, nleaftot, iloop, iunit

    ncpu = get_ncpu(repository,snapnum)
    get_nleaf = 0
    ndim = 3
    twondim = 2*ndim
    twotondim = 2**ndim

    nleaftot = 0
    iloop = 0
    
!$OMP PARALLEL &
!$OMP DEFAULT(private) &
!$OMP SHARED(iloop, nleaftot, ncpu, repository, snapnum, ndim, twondim, twotondim)
    do_allocs = .true. 
!$OMP DO SCHEDULE(DYNAMIC, 10) 
    do icpu = 1, ncpu
!!!!$OMP CRITICAL
!!!       write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
!!!            ' Reading nleaftot ',dble(iloop) / ncpu * 100,' % ',char(13)
!!!       iloop=iloop+1
!!!!$OMP END CRITICAL

       iunit = icpu+10
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       inquire(file=nomfich, exist=ok)
       if(.not. ok)then
          write(*,*)'File '//TRIM(nomfich)//' not found'    
          stop
       end if
       open(unit=iunit,file=nomfich,form='unformatted',status='old',action='read')
       read(iunit)
       read(iunit)
       read(iunit)nx,ny,nz
       ncoarse=nx*ny*nz ! Critical parameter: define the root level of the tree
       read(iunit)nlevelmax
       read(iunit)ngridmax
       read(iunit)nboundary
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       read(iunit)
       if (do_allocs) allocate(numbl(1:ncpu,1:nlevelmax),numbb(1:nboundary,1:nlevelmax))
       numbl=0;numbb=0
       ! Read levels variables
       read(iunit)!headl(1:ncpu,1:nlevelmax)
       read(iunit)!taill(1:ncpu,1:nlevelmax)
       read(iunit)numbl(1:ncpu,1:nlevelmax)
       read(iunit)!numbtot(1:10,1:nlevelmax)
       ! Read boundary linked list
       if(nboundary>0)then
          read(iunit)!headb(1:nboundary,1:nlevelmax)
          read(iunit)!tailb(1:nboundary,1:nlevelmax)
          read(iunit)numbb(1:nboundary,1:nlevelmax)
       end if
       read(iunit) ! Read free memory
       read(iunit) ! Read cpu boundaries
       read(iunit)
       ncell=ncoarse+twotondim*ngridmax
       if (do_allocs) then
          allocate(son(1:ncell),cpu_map(1:ncell))
          do_allocs = .false. 
       end if
       son=0; cpu_map=0
       ! Read coarse level
       read(iunit)son(1:ncoarse)       
       read(iunit)
       read(iunit)cpu_map(1:ncoarse)
       do ilevel=1,nlevelmax
          do ibound=1,nboundary+ncpu
             if(ibound<=ncpu)then
                ncache=numbl(ibound,ilevel)
             else
                ncache=numbb(ibound-ncpu,ilevel)
             end if
             if(ncache>0)then
!!!if (ilevel < fg_levelmin) print*,'lev < levmin',ilevel
                allocate(ind_grid(1:ncache))
                allocate(iig(1:ncache))
                read(iunit)ind_grid ! Read grid index
                read(iunit) ! Read next index
                read(iunit) ! Read prev index
                do idim=1,ndim
                   read(iunit) ! Read grid center
                end do
                read(iunit) ! Read father index
                do ind=1,twondim
                   read(iunit) ! Read nbor index
                end do
                do ind=1,twotondim
                   iskip=ncoarse+(ind-1)*ngridmax
                   read(iunit)iig  ! Read son index 
                   do i=1,ncache
                      son(ind_grid(i)+iskip)=iig(i) 
                   end do
                end do
                do ind=1,twotondim
                   iskip=ncoarse+(ind-1)*ngridmax
                   read(iunit)iig ! Read cpu map
                   do i=1,ncache
                      cpu_map(ind_grid(i)+iskip)=iig(i)
                   end do
                end do

                do ind=1,twotondim
                   read(iunit) ! Read refinement map (skip)
                end do
                deallocate(iig,ind_grid)
             end if
          end do
       end do
       close(iunit)
       ! count leaf cells
       nleaf = 0
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             nleaf = nleaf + 1
          end if
       end do

       ! nleaf is now the number of leaves on local cpu
!$OMP ATOMIC
       nleaftot = nleaftot + nleaf
    end do
!$OMP END DO
    if(.not. do_allocs) deallocate(son,cpu_map,numbl,numbb)
!$OMP END PARALLEL

    get_nleaf = nleaftot
    
    return
  end function get_nleaf


  function get_nvar(repository,snapnum)
    implicit none 
    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    character(1000)             :: nomfich
    integer(kind=4)             :: get_nvar,icpu,nvarRT

    icpu = 1
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    read(10)
    read(10)get_nvar
    close(10)
    if (read_rt_variables) then
       ! Open RT file and get nvarRT
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
       open(unit=12,file=nomfich,status='old',form='unformatted')
       read(12)
       read(12)nvarRT
       get_nvar = get_nvar + nvarRT
       close(12)
    end if

    return
  end function get_nvar


  function get_ncpu(repository,snapnum)

    implicit none

    integer(kind=4)            :: get_ncpu
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    get_ncpu = nint(get_param_real(repository,snapnum,'ncpu'))

    return
  end function get_ncpu



  function get_param_real(repository,snapnum,param)

    implicit none 

    real(kind=8)               :: get_param_real
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    character(*),intent(in)    :: param
    logical(kind=4)            :: not_ok
    character(1000)             :: nomfich
    character(512)             :: line,name,value
    integer(kind=4)            :: i
    integer(kind=4),parameter  :: param_unit = 13

    not_ok = .true.
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,'/info_',snapnum,'.txt'
    open(unit=param_unit,file=nomfich,status='old',form='formatted')
    do 
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))

       if (trim(name) .eq. trim(param)) then 
          read(value,*) get_param_real
          not_ok = .false.
       end if

    end do
2   close (param_unit)  

    if (not_ok) then 
       write(6,*) '> parameter not found in infoxxx.txt :',trim(param)
       stop
    end if

    return

  end function get_param_real

  subroutine read_conversion_scales(repository,snapnum)

    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    ! set global variables 
    dp_scale_l    = get_param_real(repository,snapnum,'unit_l')
    dp_scale_d    = get_param_real(repository,snapnum,'unit_d')
    dp_scale_t    = get_param_real(repository,snapnum,'unit_t')
    dp_scale_nH   = XH/mp * dp_scale_d      ! convert mass density (code units) to numerical density of H atoms [/cm3]
    dp_scale_v    = dp_scale_l/dp_scale_t   ! -> converts velocities into cm/s
    dp_scale_T2   = mp/kB * dp_scale_v**2   ! -> converts P/rho to T/mu, in K
    dp_scale_zsun = 1.d0/0.0127     
    dp_scale_m    = dp_scale_d * dp_scale_l**3 ! convert mass in code units to cgs. 

    return

  end subroutine read_conversion_scales

  subroutine read_cosmo_params(repository,snapnum,omega_0,lambda_0,little_h)

    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    real(kind=8),intent(out)   :: omega_0,lambda_0,little_h

    omega_0  = get_param_real(repository,snapnum,'omega_m')
    lambda_0 = get_param_real(repository,snapnum,'omega_l')
    little_h = get_param_real(repository,snapnum,'H0') / 100.0d0

    return

  end subroutine read_cosmo_params

  subroutine get_fields_from_header(dir,ts,nfields)

    implicit none

    character(1000),intent(in)  :: dir
    integer(kind=4),intent(in)  :: ts
    integer(kind=4),intent(out) :: nfields
    character(2000)             :: nomfich,line
    integer(kind=4) :: i

    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(dir),'/output_',ts,'/header_',ts,'.txt'
    open(unit=50,file=nomfich,status='old',action='read',form='formatted')
    read(50,*) ! total nb of particles
    read(50,*)
    read(50,*) ! nb of DM particles
    read(50,*)
    read(50,*) ! nb of star particles
    read(50,*)
    read(50,*) ! nb of sinks
    read(50,*)
    read(50,*) ! Field list
    read(50,'(a)') line
    close(50)

    ! parse the Field list ...
    nfields = 0
    do
       i    = scan(line,' ') ! find a blank
       nfields = nfields + 1
       ParticleFields(nfields) = trim(adjustl(line(:i)))
       line = trim(adjustl(line(i:)))
       if (len_trim(line) == 0) exit
    end do

    return

  end subroutine get_fields_from_header

  !*****************************************************************************************************************

  subroutine read_cooling(repository,snapnum)

    implicit none

    character(1000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    character(1024)            :: filename
    integer(kind=4)            :: n1,n2

    ! initialize cooling variables
    call clear_cooling

    ! read cooling variables from cooling.out file
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,"/cooling_",snapnum,".out"
    open(unit=44,file=filename,form='unformatted')
    read(44) n1,n2
    cooling%n11 = n1
    cooling%n22 = n2
    allocate(cooling%nH(n1),cooling%T2(n2))
    allocate(cooling%cool_com(n1,n2),cooling%heat_com(n1,n2),cooling%cool_com_prime(n1,n2),cooling%heat_com_prime(n1,n2))
    allocate(cooling%cool(n1,n2),cooling%heat(n1,n2),cooling%mu(n1,n2))
    allocate(cooling%cool_prime(n1,n2),cooling%heat_prime(n1,n2),cooling%metal_prime(n1,n2))
    allocate(cooling%metal(n1,n2),cooling%spec(n1,n2,6))
    read(44)cooling%nH
    read(44)cooling%T2
    read(44)cooling%cool
    read(44)cooling%heat
    read(44)cooling%cool_com
    read(44)cooling%heat_com
    read(44)cooling%metal
    read(44)cooling%cool_prime
    read(44)cooling%heat_prime
    read(44)cooling%cool_com_prime
    read(44)cooling%heat_com_prime
    read(44)cooling%metal_prime
    read(44)cooling%mu
    read(44)cooling%spec
    close(44)

    ! define useful quantities for interpolation 
    cool_int%n_nh     = n1
    cool_int%nh_start = minval(cooling%nh)
    cool_int%nh_step  = cooling%nh(2) - cooling%nh(1)
    cool_int%n_t2     = n2
    cool_int%t2_start = minval(cooling%t2)
    cool_int%t2_step  = cooling%t2(2) - cooling%t2(1)

    return

  end subroutine read_cooling

  !*****************************************************************************************************************

  subroutine clear_cooling

    implicit none

    if (cooling%n11 > 0 .or. cooling%n22 > 0) then 
       cooling%n11 = 0
       cooling%n22 = 0
       deallocate(cooling%nH,cooling%T2,cooling%metal)
       deallocate(cooling%heat_com,cooling%cool_com,cooling%heat_com_prime,cooling%cool_com_prime)
       deallocate(cooling%cool,cooling%heat,cooling%metal_prime,cooling%cool_prime)
       deallocate(cooling%heat_prime,cooling%mu,cooling%spec)
    end if

    cool_int%n_nh = 0
    cool_int%nh_start = 0.0d0
    cool_int%nh_step  = 0.0d0
    cool_int%n_t2 = 0
    cool_int%t2_start = 0.0d0
    cool_int%t2_step  = 0.0d0

    return

  end subroutine clear_cooling



  !==================================================================================
  ! STARS utilities 

  subroutine ramses_read_stars_in_domain(repository,snapnum,selection_domain,star_pos,star_age,star_mass,star_vel,star_met)

    implicit none

    character(1000),intent(in)             :: repository
    integer(kind=4),intent(in)             :: snapnum
    type(domain),intent(in)                :: selection_domain
    real(kind=8),allocatable,intent(inout) :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:)
    integer(kind=4)                        :: nstars
    real(kind=8)                           :: omega_0,lambda_0,little_h,omega_k,H0
    real(kind=8)                           :: aexp,stime,time_cu,boxsize
    integer(kind=4)                        :: ncpu,ilast,icpu,npart,i,ifield,nfields
    character(1000)                        :: filename
    integer(kind=4),allocatable            :: id(:)
    real(kind=8),allocatable               :: age(:),m(:),x(:,:),v(:,:),mets(:),skipy(:),imass(:)
    real(kind=8)                           :: temp(3)


    !Val--
    if(star_reading_method == 'fromlist') then
       
       open(unit=20, file=star_info_file, status='old', form='unformatted')

       read(20) nstars

       allocate(star_pos(3,nstars), star_age(nstars), star_mass(nstars), star_vel(3,nstars), star_met(nstars))

       read(20) star_pos(1,:) !code units
       read(20) star_pos(2,:)
       read(20) star_pos(3,:)
       read(20) star_age(:)   !Myr
       
       if(use_initial_mass) then
          read(20) star_mass(:)  !g
          read(20)
       else
          read(20)
          read(20) star_mass(:)  !g
       end if

       read(20) star_vel(1,:) ! cm/s
       read(20) star_vel(2,:)
       read(20) star_vel(3,:)
       read(20) star_met(:)   !Mass fraction of metals

       close(20)

    !--laV

    else

       ! get cosmological parameters to convert conformal time into ages
       call read_cosmo_params(repository,snapnum,omega_0,lambda_0,little_h)
       omega_k = 0.0d0
       h0      = little_h * 100.0d0
       call ct_init_cosmo(omega_0,lambda_0,omega_k,h0)
       ! compute cosmic time of simulation output (Myr)
       aexp  = get_param_real(repository,snapnum,'aexp') ! exp. factor of output
       stime = ct_aexp2time(aexp) ! cosmic time
       ! read units
       if (.not. conversion_scales_are_known) then 
          call read_conversion_scales(repository,snapnum)
          conversion_scales_are_known = .True.
       end if

       if(.not.cosmo)then
          ! read time
          time_cu = get_param_real(repository,snapnum,'time') ! code unit
          write(*,*)'Time simu [Myr] =',time_cu, time_cu*dp_scale_t/(365.*24.*3600.*1d6)
          boxsize = get_param_real(repository,snapnum,'boxlen') !!!* dp_scale_l  ! [ cm ]
          write(*,*)'boxlen =',boxsize
       endif

       ! read stars 
       nstars = get_tot_nstars(repository,snapnum)
       if (nstars == 0) then
          write(*,*) 'ERROR : no star particles in output '
          stop
       end if
       allocate(star_pos(3,nstars),star_age(nstars),star_mass(nstars),star_vel(3,nstars),star_met(nstars))
       ! get list of particle fields in outputs 
       call get_fields_from_header(repository,snapnum,nfields)
       ncpu  = get_ncpu(repository,snapnum)
       ilast = 1
       do icpu = 1, ncpu
          write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository), '/output_', snapnum, '/part_', snapnum, '.out', icpu
          open(unit=11,file=filename,status='old',form='unformatted')


          !Val--
          if(star_reading_method == 'Teyssier') then

             read(11) !ncpu2
             read(11) !ndim2
             read(11) npart
             read(11) !localseed
             read(11) !nstar_tot
             read(11) !mstar_tot
             read(11) !mstar_lost
             read(11) !nsink

             allocate(age(1:npart))
             allocate(x(1:npart,1:ndim),m(npart),imass(npart))
             allocate(id(1:npart))
             allocate(mets(1:npart))
             allocate(v(1:npart,1:ndim))
             allocate(skipy(1:npart))


             ! Read position
             do i=1,ndim
                read(11) x(1:npart,i)
             end do
             ! Read velocity
             do i=1,ndim
                read(11) v(1:ndim,i)
             end do
             ! Read mass
             read(11) m(1:npart)
             ! Read identity
             read(11) id(1:npart)
             ! Read level
             read(11)
             ! Read family
             read(11)
             ! Read tag
             read(11)
             ! Read birth epoch
             read(11) age(1:npart)
             ! Read metallicity
             read(11) mets(1:npart)

          !--laV

          !By default,  like originally in Rascas
          else
             read(11)
             read(11)
             read(11)npart
             read(11)
             read(11)
             read(11)
             read(11)
             read(11)
             allocate(age(1:npart))
             allocate(x(1:npart,1:ndim),m(npart),imass(npart))
             allocate(id(1:npart))
             allocate(mets(1:npart))
             allocate(v(1:npart,1:ndim))
             allocate(skipy(1:npart))
             do ifield = 1,nfields
                select case(trim(ParticleFields(ifield)))
                case('pos')
                   do i = 1,ndim
                      read(11) x(1:npart,i)
                   end do
                case('vel')
                   do i = 1,ndim 
                      read(11) v(1:npart,i)
                   end do
                case('mass')
                   read(11) m(1:npart)
                case('iord') 
                   read(11) id(1:npart)
                case('level')
                   read(11)
                case('tform')
                   read(11) age(1:npart)
                case('metal')
                   read(11) mets(1:npart)
                case('imass')
                   read(11) imass(1:npart)
                case default
                   ! Note: we presume here that the unknown field is an 1d array of size 1:npart
                   read(11) skipy(1:npart)
                   print*,'Error, Field unknown: ',trim(ParticleFields(ifield))
                end select
             end do

          end if
          close(11)


          if(.not.cosmo)then
             x=x/boxsize
          endif

          ! save star particles within selection region
          do i = 1,npart
             if (age(i).ne.0.0d0) then ! This is a star
                temp(:) = x(i,:)
                if (domain_contains_point(temp,selection_domain)) then ! it is inside the domain
                   if(cosmo)then
                      if (use_proper_time) then
                         star_age(ilast) = (stime - ct_proptime2time(age(i),h0))*1.d-6 ! Myr
                      else
                         ! Convert from conformal time to age in Myr
                         star_age(ilast) = (stime - ct_conftime2time(age(i)))*1.d-6 ! Myr
                      end if
                   else
                      ! convert from tborn to age in Myr
                      star_age(ilast)   = max(0.d0, (time_cu - age(i)) * dp_scale_t / (365.d0*24.d0*3600.d0*1.d6))
                   endif
                   if (use_initial_mass) then
                      !Val--  From Romain Teyssier version : have to divide mass by 1-eta_sn
                      if(star_reading_method == 'Teyssier') then
                         star_mass(ilast) = m(i) * dp_scale_m ! [g]
                         if(star_age(ilast) > t_sne_Myr) star_mass(ilast) = m(i) * dp_scale_m / (1d0 - eta_sn)
                      !laV--
                      !Normal way, like initially in Rascas
                      else
                         star_mass(ilast) = imass(i) * dp_scale_m ! [g]
                      end if
                   else
                      star_mass(ilast) = m(i)     * dp_scale_m ! [g]
                   end if
                   star_pos(:,ilast) = x(i,:)              ! [code units]
                   star_vel(:,ilast) = v(i,:) * dp_scale_v ! [cm/s]
                   star_met(ilast) = mets(i) 
                   ilast = ilast + 1
                end if
             end if
          end do

          deallocate(age,m,x,id,mets,v,skipy,imass)

       end do

       ! resize star arrays
       nstars = ilast-1
       ! ages
       allocate(age(nstars))
       age = star_age(1:nstars)
       deallocate(star_age)
       allocate(star_age(nstars))
       star_age = age
       deallocate(age)
       ! masses
       allocate(m(nstars))
       m = star_mass(1:nstars)
       deallocate(star_mass)
       allocate(star_mass(nstars))
       star_mass = m
       deallocate(m)
       ! positions
       allocate(x(3,nstars))
       do i = 1,nstars 
          x(:,i) = star_pos(:,i)
       end do
       deallocate(star_pos)
       allocate(star_pos(3,nstars))
       star_pos = x
       deallocate(x)
       ! velocities
       allocate(v(3,nstars))
       do i = 1,nstars 
          v(:,i) = star_vel(:,i)
       end do
       deallocate(star_vel)
       allocate(star_vel(3,nstars))
       star_vel = v
       deallocate(v)
       ! metals
       allocate(mets(nstars))
       mets = star_met(1:nstars)
       deallocate(star_met)
       allocate(star_met(nstars))
       star_met = mets
       deallocate(mets)

    end if
    
    return
  end subroutine ramses_read_stars_in_domain


  function get_tot_nstars(dir,ts)

    implicit none 

    integer(kind=4),intent(in) :: ts
    character(1000),intent(in) :: dir
    !Val--
    character(2000)            :: nomfich, useless
    !--Val
    integer(kind=4)            :: get_tot_nstars

    get_tot_nstars = 0
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(dir),'/output_',ts,'/header_',ts,'.txt'
    open(unit=50,file=nomfich,status='old',action='read',form='formatted')
    !Val--
    if(star_reading_method == 'Teyssier') then
       read(50,*) ! total nb of particles
       read(50,*)
       read(50,*) ! nb of DM particles
       read(50,*)
       read(50,*) ! nb of star particles 
       read(50,*) ! nb of DM particles
       read(50,*)
       read(50,*) ! nb of star particles 
       read(50,*) useless, get_tot_nstars
    !--laV
    else
       read(50,*) ! total nb of particles
       read(50,*)
       read(50,*) ! nb of DM particles
       read(50,*)
       read(50,*) ! nb of star particles 
       read(50,*) get_tot_nstars
    end if

    close(50)

    return

  end function get_tot_nstars

  ! conformal time utils :
    function ct_conftime2time(tau)

    ! return look-back time in yr
    
    implicit none 
    
    real(kind=8),intent(in) :: tau
    real(kind=8)            :: ct_conftime2time
    integer(kind=4)         :: i
    
    
    ! locate bracketing conf. times
    i = 1
    do while(tau_frw(i) > tau .and. i < n_frw)
       i = i + 1
    end do
    ! Interploate time
    ct_conftime2time = t_frw(i) * (tau-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
         & t_frw(i-1)        * (tau-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
    
    return

  end function ct_conftime2time
  
  
  function ct_proptime2time(tau,h0)
    ! return look-back time in yr
    
    implicit none 
    real(kind=8),intent(in) :: tau,h0
    real(kind=8)            :: ct_proptime2time
    
    ct_proptime2time = tau / (h0 / 3.08d19) / (365.25*24.*3600.)
    return
  end function ct_proptime2time
  
  
  function ct_aexp2time(aexp)

    ! return look-back time in yr

    implicit none

    real(kind=8),intent(in) :: aexp
    real(kind=8)            :: ct_aexp2time
    integer(kind=4)         :: i

    ! find bracketting aexp's 
     i = 1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i = i + 1
     end do
     ! Interploate time
     ct_aexp2time = t_frw(i) * (aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)    * (aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))

    return
    
  end function ct_aexp2time

  
  subroutine ct_init_cosmo(omega_m,omega_l,omega_k,h0)
    
    ! h0 is in km/s/Mpc

    implicit none 
    real(kind=8),intent(in) :: omega_m,omega_l,omega_k,h0
    real(kind=8)            :: time_tot

    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
    call ct_friedman(omega_m,omega_l,omega_k,1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
    ! convert time to yr
    t_frw = t_frw / (h0 / 3.08d19) / (365.25*24.*3600.)

    return
    
  end subroutine ct_init_cosmo


  subroutine ct_clear_cosmo
    
    implicit none
    
    deallocate(aexp_frw,hexp_frw,tau_frw,t_frw)

    return

  end subroutine ct_clear_cosmo


  subroutine ct_friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
       & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

    implicit none
    integer::ntable
    real(kind=8)::O_mat_0, O_vac_0, O_k_0
    real(kind=8)::alpha,axp_min,age_tot
    real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
    ! ######################################################!
    ! This subroutine assumes that axp = 1 at z = 0 (today) !
    ! and that t and tau = 0 at z = 0 (today).              !
    ! axp is the expansion factor, hexp the Hubble constant !
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
    ! time, and t the look-back time, both in unit of 1/H0. !
    ! alpha is the required accuracy and axp_min is the     !
    ! starting expansion factor of the look-up table.       !
    ! ntable is the required size of the look-up table.     !
    ! ######################################################!
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dtau,dt
    real(kind=8)::tau,t
    integer::nstep,nout,nskip

    !  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
    !     write(*,*)'Error: non-physical cosmological constants'
    !     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
    !     write(*,*)'The sum must be equal to 1.0, but '
    !     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
    !     stop
    !  end if

    axp_tau = 1.0D0
    axp_t = 1.0D0
    tau = 0.0D0
    t = 0.0D0
    nstep = 0

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

    end do

    age_tot=-t
!!$    write(*,666)-t
!!$666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

    nskip=nstep/ntable

    axp_t = 1.d0
    t = 0.d0
    axp_tau = 1.d0
    tau = 0.d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

       if(mod(nstep,nskip)==0)then
          nout=nout+1
          t_out(nout)=t
          tau_out(nout)=tau
          axp_out(nout)=axp_tau
          hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
       end if

    end do
    t_out(ntable)=t
    tau_out(ntable)=tau
    axp_out(ntable)=axp_tau
    hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  contains
    function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
      implicit none 
      real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
      dadtau = axp_tau*axp_tau*axp_tau *  &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
           &     O_k_0   * axp_tau )
      dadtau = sqrt(dadtau)
      return
    end function dadtau
    
    function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
      implicit none
      real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
      dadt   = (1.0D0/axp_t)* &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_t*axp_t*axp_t + &
           &     O_k_0   * axp_t )
      dadt = sqrt(dadt)
      return
    end function dadt
    
  end subroutine ct_friedman


  subroutine read_ramses_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:8) == '[ramses]') then
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
          !Val--
          case('star_reading_method')
             read(value,*) star_reading_method
          case('star_info_file')
            write(star_info_file,'(a)') trim(value)
          case('t_sne_Myr')
             read(value,*) t_sne_Myr
          case('eta_sn')
             read(value,*) eta_sn
          !--laV
          case ('self_shielding')
             read(value,*) self_shielding
          case ('ramses_rt')
             read(value,*) ramses_rt
          case ('read_rt_variables')
             read(value,*) read_rt_variables
          case ('verbose')
             read(value,*) verbose
          case ('use_initial_mass')
             read(value,*) use_initial_mass
          case ('cosmo')
             read(value,*) cosmo
          case ('use_proper_time')
             read(value,*) use_proper_time
          case('itemp') ! index of thermal pressure
             read(value,*) itemp
          case('imetal')! index of metallicity  
             read(value,*) imetal
          case('ihii') ! index of HII fraction 
             read(value,*) ihii
          case ('iheii') ! index of HeII fraction 
             read(value,*) iheii
          case('iheiii') ! index of HeIII fraction 
             read(value,*) iheiii
          end select
       end do
    end if
    close(10)
    return
  end subroutine read_ramses_params


  
  subroutine print_ramses_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[ramses]'
       !Val--
       write(unit,'(a,a)')  '  star_reading_method  = ', star_reading_method
       write(unit,'(a,a)')  '  star_info_file       = ', trim(star_info_file)
       write(unit,'(a,ES10.3)')  '  t_sne_Myr       = ', t_sne_Myr
       write(unit,'(a,ES10.3)')  '  eta_sn          = ', eta_sn
       !--laV
       write(unit,'(a,L1)') '  self_shielding       = ', self_shielding
       write(unit,'(a,L1)') '  ramses_rt            = ', ramses_rt
       write(unit,'(a,L1)') '  read_rt_variables    = ', read_rt_variables
       write(unit,'(a,L1)') '  use_initial_mass     = ', use_initial_mass
       write(unit,'(a,L1)') '  cosmo                = ', cosmo
       write(unit,'(a,L1)') '  use_proper_time      = ', use_proper_time
       write(unit,'(a,L1)') '  verbose              = ', verbose
       write(unit,'(a,i2)') '  itemp                = ', itemp
       write(unit,'(a,i2)') '  imetal               = ', imetal
       write(unit,'(a,i2)') '  ihii                 = ', ihii
       write(unit,'(a,i2)') '  iheii                = ', iheii
       write(unit,'(a,i2)') '  iheiii               = ', iheiii
    else
       write(*,'(a,a,a)') '[ramses]'
       !Val--
       write(*,'(a,a)')  '  star_reading_method  = ', star_reading_method
       write(*,'(a,a)')  '  star_info_file       = ', trim(star_info_file)
       write(*,'(a,ES10.3)')    't_sne_Myr       = ', t_sne_Myr
       write(*,'(a,ES10.3)')    'eta_sn          = ', eta_sn
       !--laV
       write(*,'(a,L1)') '  self_shielding       = ', self_shielding
       write(*,'(a,L1)') '  ramses_rt            = ', ramses_rt
       write(*,'(a,L1)') '  read_rt_variables    = ', read_rt_variables
       write(*,'(a,L1)') '  use_initial_mass     = ', use_initial_mass
       write(*,'(a,L1)') '  cosmo                = ', cosmo
       write(*,'(a,L1)') '  use_proper_time      = ', use_proper_time
       write(*,'(a,L1)') '  verbose              = ', verbose
       write(*,'(a,i2)') '  itemp                = ', itemp
       write(*,'(a,i2)') '  imetal               = ', imetal
       write(*,'(a,i2)') '  ihii                 = ', ihii
       write(*,'(a,i2)') '  iheii                = ', iheii
       write(*,'(a,i2)') '  iheiii               = ', iheiii
    end if
    
    return
  end subroutine print_ramses_params




end module module_ramses
!==================================================================================
!==================================================================================
