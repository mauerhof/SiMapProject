module star_formation_module

  ! ISSUES/QUESTIONS 
  ! 
  ! - Should we force use_proper_time=.true.  (i assume yes)
  ! 
  ! - (for feedback) how do we get the cells around the current cell
  ! (voir pm/move_fine.f90). Attention: les cellules voisines peuvent
  ! etre dans d'autres cpus ... on peut les lire dans le buffer mais
  ! pas les modifier... (cf kinetic feedback de Yohan).
  ! 
  ! - We use a hard-coded omega_b here ... 
  ! 
  ! - initialization of random number generator seeds -> do elsewhere ? 
  ! 
  ! - solar masses used here and there with different hardcoded values (grep 33 in *f90).
  ! 

  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  use random
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  public 

  ! global parameters that define the star formation model (TAKEN FROM AMR_PARAMETERS)
!!  real(dp) ::m_star =-1.0      ! Star particle mass in units of mass_sph
!!  real(dp) ::n_star =0.1D0     ! Star formation density threshold in H/cc
!!  real(dp) ::t_star =0.0D0     ! Star formation time scale in Gyr
!!  real(dp) ::eps_star=0.0D0    ! Star formation efficiency (0.02 at n_star=0.1 gives t_star=8 Gyr)
!!  real(dp) ::del_star=2.D2     ! Minimum overdensity to define ISM
!!  integer  ::iseed=0           ! Seed for stochastic star formation (TAKEN FROM PM_PARAMETERS) 
                               ! WARNING: THIS SHOULD BE SHARED SOMEHOW ACCROSS CPUs... 
  ! new parameters ... 
  real(dp),parameter :: star_particle_mass = 1d3 ! mass of a star particle in Msun
  real(dp),parameter :: Msun2g = 1.98892d33
  real(dp),parameter :: Gyr2s = 1d9 * 365. * 24. * 3600d0 

  logical,parameter  :: jejose = .true. ! local verbose ... 
  
contains 

  subroutine star_formation(ilevel)

    implicit none 

    integer(kind=4),intent(in) :: ilevel
    real(dp)                   :: scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
    integer                    :: nx_loc
    real(dp)                   :: dx, dx_loc, vol_loc, scale 
    real(dp)                   :: t0,nISM,nCOM,d0,d00,mstar


    integer ,dimension(1:nvector),save::ind_grid,ind_cell,nstar  !JB: why the save attribute ? 
    logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.  !JB: why the save attribute ? 
    real(dp),dimension(1:3)::skip_loc ! JB: this could be a global defined where icoarse_min is ... 


    integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  

    ! checks for compatibility ...  (to clean up later)
    if (.not. use_proper_time) then 
       write(*,*) 'ABORTING : Star formation requires use_proper_time==.true. '
       call clean_stop
    end if
    if(numbtot(1,ilevel)==0) return ! JB: this should go in amr-step before call to star_formation ? 
    if(.not. hydro)return ! JB: this should go in amr-step before call to star_formation
    if(ndim.ne.3)return   ! JB: this should go in amr-step before call to star_formation
    if(static)return      ! JB: this should go in amr-step before call to star_formation

    ! If necessary, initialize random number generator  
    if(localseed(1)==-1)then
       call rans(ncpu,iseed,allseed)
       localseed=allseed(myid,1:IRandNumSize)
    end if
   
    ! ----------------------------------------------------------
    ! compute values of model parameters with current user units
    ! ----------------------------------------------------------
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2) ! conversion factors to cgs units at current expansion factor
    t0    = t_star*Gyr2s/scale_t   ! timescale for star formation in code units
    nISM  = n_star                                ! ISM density threshold (H/cc)
    nCOM  = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
    nISM  = MAX(nCOM,nISM) 
    d0    = nISM/scale_nH                         ! SF threshhold in code units
    d00   = n_star/scale_nH                       ! rho_star (efficiency parameter) in code units
    mstar = star_particle_mass * Msun2g / (scale_d * scale_l**3)  ! particle mass in code units. 
    ! compute volume of of a cell at level ilevel
    nx_loc  = (icoarse_max-icoarse_min+1)
    scale   = boxlen/dble(nx_loc)
    dx      = 0.5D0**ilevel          ! cell size in units of boxlen=1
    dx_loc  = dx*scale               ! cell size in user units
    vol_loc = dx_loc**ndim           ! volume of an ilevel cell (user units)
    
    if (jejose) then 
       write(*,*) 'Parameters for star formation at level ', ilevel
       write(*,*) '...'
       write(*,*) 'nx_loc = ',nx_loc
    end if

    ! ----------------------------------------------
    ! Convert hydro variables to primitive variables
    ! ----------------------------------------------
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector ! loop on ncache grids, by sweeps of nvector elements
       ngrid=MIN(nvector,ncache-igrid+1)
       do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
       end do
       do ind=1,twotondim  ! loop on cells 
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
             ind_cell(i)=iskip+ind_grid(i)
          end do
          do i=1,ngrid
             d=uold(ind_cell(i),1)     ! density 
             u=uold(ind_cell(i),2)/d   ! x-velocity
             v=uold(ind_cell(i),3)/d   ! 
             w=uold(ind_cell(i),4)/d   !
             e=uold(ind_cell(i),5)/d   ! total energy
#ifdef SOLVERmhd
             bx1=uold(ind_cell(i),6)
             by1=uold(ind_cell(i),7)
             bz1=uold(ind_cell(i),8)
             bx2=uold(ind_cell(i),nvar+1)
             by2=uold(ind_cell(i),nvar+2)
             bz2=uold(ind_cell(i),nvar+3)
             e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d  
#endif
             e=e-0.5d0*(u**2+v**2+w**2)
             uold(ind_cell(i),1)=d
             uold(ind_cell(i),2)=u
             uold(ind_cell(i),3)=v
             uold(ind_cell(i),4)=w
             uold(ind_cell(i),5)=e     ! specific thermal energy
          end do
          do ivar=imetal,nvar
             do i=1,ngrid
                d=uold(ind_cell(i),1)
                w=uold(ind_cell(i),ivar)/d
                uold(ind_cell(i),ivar)=w
             end do
          end do
       end do
    end do

    
    
    !------------------------------------------------
    ! Compute number of new stars in each cell
    !------------------------------------------------
    ntot=0
    ! Loop over grids
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
       ngrid=MIN(nvector,ncache-igrid+1)
       do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
       end do
       ! Star formation criterion ---> logical array ok(i)
       do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
             ind_cell(i)=iskip+ind_grid(i)
          end do
          ! Flag leaf cells
          do i=1,ngrid
             ok(i)=son(ind_cell(i))==0
          end do
          ! Density criterion
          do i=1,ngrid
             d=uold(ind_cell(i),1)
             if(d<=d0)ok(i)=.false. 
          end do
          ! Geometrical criterion
          if(ivar_refine>0)then
             do i=1,ngrid
                d=uold(ind_cell(i),ivar_refine)
                if(d<=var_cut_refine)ok(i)=.false.
             end do
          endif
          ! Calculate number of new stars in each cell using Poisson statistics
          do i=1,ngrid
             nstar(i)=0
             if(ok(i))then
                ! Compute mean number of events
                d=uold(ind_cell(i),1)
                mcell=d*vol_loc
                tstar=t0*sqrt(d00/d)
                PoissMean=dtnew(ilevel)/tstar*mcell/mstar
                ! Compute Poisson realisation
                call poissdev(localseed,PoissMean,nstar(i))
                ! Compute depleted gas mass
                mgas=nstar(i)*mstar
                ! Security to prevent more than 90% of gas depletion
                if (mgas > 0.9*mcell) then
                   nstar_corrected=int(0.9*mcell/mstar)
                   mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
                   nstar(i)=nstar_corrected
                endif
                ! Compute new stars local statistics
                mstar_tot=mstar_tot+nstar(i)*mstar
                if(nstar(i)>0)then
                   ntot=ntot+1   
                endif
             endif
          enddo
          ! Store nstar in array flag2
          do i=1,ngrid
             flag2(ind_cell(i))=nstar(i)
          end do
       end do
    end do


    !---------------------------------
    ! Check for free particle memory
    !---------------------------------
    ok_free=(numbp_free-ntot-ndebris_tot)>=0
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#else
    numbp_free_tot=numbp_free
#endif
    if(.not. ok_free)then
       write(*,*)'No more free memory for particles'
       write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
       call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
       stop
#endif
    end if


    ! compute position of each cell and skip_loc thing 
    skip_loc=(/0.0d0,0.0d0,0.0d0/)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    do ind=1,twotondim  
       iz=(ind-1)/4
       iy=(ind-1-4*iz)/2
       ix=(ind-1-2*iy-4*iz)
       xc(ind,1)=(dble(ix)-0.5D0)*dx   ! Cells center position relative to grid center position
       xc(ind,2)=(dble(iy)-0.5D0)*dx
       xc(ind,3)=(dble(iz)-0.5D0)*dx
    end do


    ! --------------------------------------------------------------
    ! Create new star particles (spawn'em and update gas properties) 
    ! --------------------------------------------------------------
    ! Starting identity number
    if(myid==1)then
       index_star=nstar_tot-ntot_all
    else
       index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
    end if
    ! Loop over grids
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
       ngrid=MIN(nvector,ncache-igrid+1)
       do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
       end do
       ! Loop over cells
       do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
             ind_cell(i)=iskip+ind_grid(i)
          end do
          ! Flag cells with at least one new star
          do i=1,ngrid
             ok(i)=flag2(ind_cell(i))>0
          end do
          ! Gather new star arrays
          nnew=0
          do i=1,ngrid
             if (ok(i))then
                nnew=nnew+1
                ind_grid_new(nnew)=ind_grid(i)
                ind_cell_new(nnew)=ind_cell(i)
             end if
          end do
          ! Update linked list for stars
          call remove_free(ind_part,nnew)
          call add_list(ind_part,ind_grid_new,ok_new,nnew)
          ! Update linked list for debris
          if(f_w>0)then
             call remove_free(ind_debris,nnew)
             call add_list(ind_debris,ind_grid_new,ok_new,nnew)
          endif

          ! Calculate new star particle and modify gas density
          do i=1,nnew
             index_star=index_star+1
             ! Get gas variables
             n=flag2(ind_cell_new(i))
             d=uold(ind_cell_new(i),1)
             u=uold(ind_cell_new(i),2)
             v=uold(ind_cell_new(i),3)
             w=uold(ind_cell_new(i),4)
             x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
             y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
             z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
             if(metal)zg=uold(ind_cell_new(i),imetal)
             ! Set star particle variables
             tp(ind_part(i))=birth_epoch  ! Birth epoch
             mp(ind_part(i))=n*mstar      ! Mass
             levelp(ind_part(i))=ilevel   ! Level
             idp(ind_part(i))=index_star  ! Star identity
             xp(ind_part(i),1)=x
             xp(ind_part(i),2)=y
             xp(ind_part(i),3)=z
             vp(ind_part(i),1)=u
             vp(ind_part(i),2)=v
             vp(ind_part(i),3)=w
             if(metal)zp(ind_part(i))=zg  ! Initial star metallicity
             ! Set GMC particle variables
             if(f_w>0)then
                ! Compute GMC mass without more than 50% of gas depletion
                mdebris=min(f_w*n*mstar,0.5*d*vol_loc-n*mstar)
                ! Add supernova ejecta
                mdebris=mdebris+eta_sn*n*mstar
                ! Remove ejecta from the long lived star mass
                mp(ind_part(i))=n*mstar-eta_sn*n*mstar
                ! Set GMC particle variables
                tp(ind_debris(i))=birth_epoch  ! Birth epoch
                mp(ind_debris(i))=mdebris      ! Mass
                levelp(ind_debris(i))=ilevel   ! Level
                idp(ind_debris(i))=-n          ! Number of individual stars
                xp(ind_debris(i),1)=x
                xp(ind_debris(i),2)=y
                xp(ind_debris(i),3)=z
                vp(ind_debris(i),1)=u
                vp(ind_debris(i),2)=v
                vp(ind_debris(i),3)=w
                ! GMC metallicity + yield from ejecta 
                if(metal)zp(ind_debris(i))=zg+eta_sn*yield*(1-zg)*n*mstar/mdebris
             endif
          end do  ! End loop over new star particles
          ! Modify gas density according to mass depletion
          do i=1,ngrid
             if(flag2(ind_cell(i))>0)then
                n=flag2(ind_cell(i))
                d=uold(ind_cell(i),1)
                uold(ind_cell(i),1)=max(d-n*dstar*(1.0+f_w),0.5*d)
             endif
          end do
       end do ! End loop over cells
    end do    ! End loop over grids
    
    ! Convert hydro variables back to conservative variables


    return
    
  end subroutine star_formation




end module star_formation_module


