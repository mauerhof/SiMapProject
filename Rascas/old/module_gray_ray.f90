module module_gray_ray

  ! This is a trimmed and slightly modified version of module_photon.
  !
  ! THIS IMPLENMENTATION HAS NO WAVELENGTH / VELOCITY INFO (used for escape fraction measurements). 
  
  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain
  use module_utils, only: path, isotropic_direction

  implicit none

  ! todonext define accuracy
  real(kind=8),parameter :: accuracy=1.d-15
  integer(kind=4) :: iran = -10
  logical::identical_ray_distribution=.false. ! Same ray directions for all sources?
  logical:: use_halos = .false.


  type ray_type
     integer(kind=4)           :: ID       ! a positive unique ID 
     real(kind=8)              :: dist     ! distance traveled along ray (box units)
     real(kind=8)              :: tau      ! integrated opacity along ray
     real(kind=8),dimension(3) :: x_em     ! emission location (box units)
     integer(kind=4)           :: halo_ID  ! a positive unique ID 
     real(kind=8),dimension(3) :: k_em     ! emission direction == propagation direction (normalised vector)
  end type ray_type

  type halo_type
     integer(kind=4)           :: ID       ! a positive unique ID
     type(domain)              :: domain   ! rvir and center (box units)
  end type halo_type

  type(halo_type),dimension(:),allocatable :: halos

  public  :: init_rays_from_file, init_halos_from_file, dump_rays, ray_advance !column_density
  private :: path

contains

  
  subroutine ray_advance(ray,domesh,domaine_calcul,maxdist,n,flux)

    type(ray_type),intent(inout)   :: ray            ! a ray 
    type(mesh),intent(in)          :: domesh         ! mesh
    type(domain),intent(in)        :: domaine_calcul ! domaine dans lequel on propage les photons...
    integer(kind=4),intent(in)     :: n
    real(kind=8),intent(in)        :: maxdist      ! stop propagation after a distance of maxdist
    real(kind=8),intent(inout)     :: flux(1,n)        ! stop propagation after a distance of maxdist
    type(gas)                      :: cell_gas       ! gas in the current cell 
    integer(kind=4)                :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                   :: cell_size, cell_size_cm, scalar, nu_cell, maxdist_cm
    real(kind=8),dimension(3)      :: ppos,ppos_cell ! working coordinates of photon (in box and in cell units)
    real(kind=8)                   :: distance_to_border,distance_to_border_cm
    real(kind=8)                   :: dist, tau_cell, tau, dborder, d
    integer(kind=4)                :: i, icellnew, npush
    real(kind=8),dimension(3)      :: vgas, kray, cell_corner, posoct
    logical                        :: flagoutvol, in_domain
    real(kind=8)                   :: excess, tmp
    
    ! initialise ray tracing 
    ppos  = ray%x_em   ! emission position 
    kray  = ray%k_em   ! propagation direction
    dist  = 0.0d0      ! distance covered
    tau   = 0.0d0      ! corresponding optical depth
    maxdist_cm = maxdist ! dist_tot is now provided in cm ... 


    ! check that the ray starts in the domain
    if (.not. domain_contains_point(ppos,domaine_calcul)) then
       print * ,'Ray outside domain at start '
       stop
    end if
    
    ! find the (leaf) cell in which the photon is, and define all its indices
    icell = in_cell_finder(domesh,ppos)
    if(domesh%son(icell)>=0)then
       print*,'ERROR: not a leaf cell ',ppos
       stop
    endif
    ileaf = - domesh%son(icell)
    ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct

    ! advance ray until it escapes the computational domain ... 
    ray_propagation : do
       
       ! gather properties properties of current cell
       cell_level   = domesh%octlevel(ioct)      ! level of current cell
       cell_size    = 0.5d0**cell_level          ! size of current cell in box units
       cell_size_cm = cell_size * box_size_cm    ! size of the current cell in cm
       cell_gas     = domesh%gas(ileaf)
       ! compute position of photon in current-cell units
       posoct(:)    = domesh%xoct(ioct,:)
       cell_corner  = get_cell_corner(posoct,ind,cell_level)   ! position of cell corner, in box units.
       ppos_cell    = (ppos - cell_corner) / cell_size         ! position of photon in cell units (x,y,z in [0,1] within cell)
       if((ppos_cell(1)>1.0d0).or.(ppos_cell(2)>1.0d0).or.(ppos_cell(3)>1.0d0).or. &
            (ppos_cell(1)<0.0d0).or.(ppos_cell(2)<0.0d0).or.(ppos_cell(3)<0.0d0))then
          print*,"ERROR: problem in computing ppos_cell ",ppos_cell(1),ppos_cell(2),ppos_cell(3)
          stop
       endif
       
       ! compute distance of photon to border of cell or domain along propagation direction
       distance_to_border    = path(ppos_cell,kray) * cell_size            ! in box units
       distance_to_border    = min(distance_to_border, &
            & domain_distance_to_border_along_k(ppos,kray,domaine_calcul)) ! in box units
       distance_to_border_cm = distance_to_border * box_size_cm ! cm
       
       
       ! compute (total) optical depth along ray in cell 
       tau_cell = gas_get_tau_gray(cell_gas, distance_to_border_cm)

       ! update traveled distance and optical depth
       !if(ray%ID==30327)print*,distance_to_border_cm
       dist = dist + distance_to_border_cm/2
       tau  = tau + tau_cell/2

       flux(1,ileaf) = exp(-tau)/dist**2

       dist = dist + distance_to_border_cm/2
       tau  = tau + tau_cell/2

       if (dist > maxdist_cm) then
          excess   = maxdist_cm - dist
          ray%dist = dist
          ray%tau  = tau - (excess / distance_to_border_cm)*tau_cell
          exit ray_propagation
       endif


       ! update head of ray position
       ppos = ppos + kray * distance_to_border *(1.0d0 + epsilon(1.0d0))
       ! correct for periodicity
       do i=1,3
          if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
          if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
       enddo
       
       ! check if photon still in computational domain after position update
       in_domain = domain_contains_point(ppos,domaine_calcul)
       if (.not.(in_domain)) then
          ray%dist = dist
          ray%tau  = tau
          ray%halo_ID = 0
          !print*,'WARNING: escaping domain before dist_tot or maxtau is reached... 1'
          !print*,'initial distance to border of domain [cm] : ',domain_distance_to_border(ray%x_em,domaine_calcul)*box_size_cm
          !print*,'dist_tot [cm]                              : ',dist_tot_cm
          exit ray_propagation  ! ray is done 
       end if
       ! Ray moves to next cell : find it
       call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
       ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
       ! -> give it an extra push untill it is out. 
       npush = 0
       do while (icell==icellnew)
          npush = npush + 1
          if (npush>100) then
             print*,'Too many pushes, npush>100 '
             print*,'ray id = ',ray%ID
             print*,'ray initial position = ',ray%x_em
             print*,'ray direction = ',ray%k_em
             print*,'halo domain type = ',domaine_calcul%type
             print*,'halo center = ',domaine_calcul%sp%center
             print*,'halo radius = ',domaine_calcul%sp%radius
             print*,'ray position = ',ppos
             print*,'distance to border = ',distance_to_border,distance_to_border_cm
             tmp = sqrt(sum((ppos-domaine_calcul%sp%center)**2))
             print*,'dist of ray from center = ',tmp
             print*,'dist - rvir = ',tmp - domaine_calcul%sp%radius
             stop
          endif

          ! hack
          do i=1,3
             ppos(i) = ppos(i) + kray(i) * 1d7 * epsilon(ppos(i))
          end do
          !!$ ppos(1) = ppos(1) + merge(-1.0d0,1.0d0,kray(1)<0.0d0) * epsilon(ppos(1))
          !!$ ppos(2) = ppos(2) + merge(-1.0d0,1.0d0,kray(2)<0.0d0) * epsilon(ppos(2))
          !!$ ppos(3) = ppos(3) + merge(-1.0d0,1.0d0,kray(3)<0.0d0) * epsilon(ppos(3))
          ! kcah
          
          ! correct for periodicity
          do i=1,3
             if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
             if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
          enddo
          ! test that we are still in domain before calling WhereIsPhotonGoing... 
          in_domain = domain_contains_point(ppos,domaine_calcul)
          if (.not.(in_domain)) then
             ray%dist = dist
             ray%tau  = tau
             ray%halo_ID = 0
             !print*,'WARNING: escaping domain before dist_tot or maxtau is reached... 2'
             !print*,'initial distance to border of domain [cm] : ',domain_distance_to_border(ray%x_em,domaine_calcul)*box_size_cm
             !print*,'dist_tot [cm]                              : ',dist_tot_cm
             exit ray_propagation  ! ray is done 
          end if
          call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
       end do
       if (npush > 1) print*,'WARNING : npush > 1 needed in module_gray_ray:propagate.'
       
       ! check if photon outside of cpu domain (flagoutvol)
       if(flagoutvol)then
          ! photon out of cpu domain -> we have a problem... 
          print*,'ERROR: photon out of CPU domain when it should not ... '
          stop
       endif
       
       ! else, new cell is in the cpu domain and photon goes to this new cell
       icell = icellnew
       if(domesh%son(icell)>=0)then
          print*,'ERROR: not a leaf cell',icell,flagoutvol
          print*,'This should not happen in module_ray (on single domains). '
          stop
       endif
       ileaf = - domesh%son(icell)
       ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
       ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
       
    end do ray_propagation
    

  end subroutine ray_advance

  

  ! function column_density(domesh, domaine_calcul, pos_i, k, maxdist_cm)
 
  !   type(mesh),intent(in)          :: domesh         ! mesh
  !   type(domain),intent(in)        :: domaine_calcul ! domaine dans lequel on propage les photons...
  !   real(kind=8),intent(in)        :: pos_i(3), k(3), maxdist_cm
  !   real(kind=8)                   :: column_density
  !   type(gas)                      :: cell_gas       ! gas in the current cell 
  !   integer(kind=4)                :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
  !   real(kind=8)                   :: cell_size, cell_size_cm, scalar, nu_cell, dist_tot_cm
  !   real(kind=8),dimension(3)      :: pos, pos_cell ! working coordinates of photon (in box and in cell units)
  !   real(kind=8)                   :: distance_to_border,distance_to_border_cm
  !   real(kind=8)                   :: dist, tau_cell, tau, dborder, d
  !   integer(kind=4)                :: i, icellnew, npush
  !   real(kind=8),dimension(3)      :: vgas, kray, cell_corner, posoct
  !   logical                        :: flagoutvol, in_domain
  !   real(kind=8)                   :: excess, tmp

    
  !   pos = pos_i             ! emission position
  !   dist  = 0.0d0           ! distance covered
  !   column_density = 0d0


  !   ! check that the ray starts in the domain
  !   if (.not. domain_contains_point(pos,domaine_calcul)) then
  !      print * ,'Ray outside domain at start '
  !      stop
  !   end if
    
  !   ! find the (leaf) cell in which the photon is, and define all its indices
  !   icell = in_cell_finder(domesh,pos)
  !   if(domesh%son(icell)>=0)then
  !      print*,'ERROR: not a leaf cell ',pos
  !      stop
  !   endif
  !   ileaf = - domesh%son(icell)
  !   ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
  !   ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct

  !   ! advance ray until it escapes the computational domain ... 
  !   ray_propagation : do
       
  !      ! gather properties properties of current cell
  !      cell_level   = domesh%octlevel(ioct)      ! level of current cell
  !      cell_size    = 0.5d0**cell_level          ! size of current cell in box units
  !      cell_size_cm = cell_size * box_size_cm    ! size of the current cell in cm
  !      cell_gas     = domesh%gas(ileaf)
  !      ! compute position of photon in current-cell units
  !      posoct(:)    = domesh%xoct(ioct,:)
  !      cell_corner  = get_cell_corner(posoct,ind,cell_level)   ! position of cell corner, in box units.
  !      pos_cell    = (pos - cell_corner) / cell_size         ! position of photon in cell units (x,y,z in [0,1] within cell)
  !      if((pos_cell(1)>1.0d0).or.(pos_cell(2)>1.0d0).or.(pos_cell(3)>1.0d0).or. &
  !           (pos_cell(1)<0.0d0).or.(pos_cell(2)<0.0d0).or.(pos_cell(3)<0.0d0)) then
  !         print*,"ERROR: problem in computing pos_cell ",pos_cell(1),pos_cell(2),pos_cell(3)
  !         stop
  !      endif
       
  !      ! compute distance of photon to border of cell or domain along propagation direction
  !      distance_to_border    = path(pos_cell,kray) * cell_size            ! in box units
  !      distance_to_border    = min(distance_to_border, &
  !           & domain_distance_to_border_along_k(pos,kray,domaine_calcul)) ! in box units
  !      distance_to_border_cm = distance_to_border * box_size_cm ! cm

  !      ! update traveled distance and column density
  !      dist = dist + distance_to_border_cm
  !      column_density = column_density + cell_gas%nHI*distance_to_border_cm

  !      ! check if we reached tau or distance limits
  !      if (dist > maxdist_cm) then
  !         ! dist exceeding boundary -> exit. 
  !         print*, 'Stopping the computation of column density because maxdist is reached'
  !         exit ray_propagation  ! no need to update other unused properties of ray. 
  !      end if

  !      ! update head of ray position
  !      pos = pos + k * distance_to_border *(1.0d0 + epsilon(1.0d0))
  !      ! correct for periodicity
  !      do i=1,3
  !         if (pos(i) < 0.0d0) pos(i)=pos(i)+1.0d0
  !         if (pos(i) > 1.0d0) pos(i)=pos(i)-1.0d0
  !      enddo
       
  !      ! check if photon still in computational domain after position update
  !      in_domain = domain_contains_point(pos,domaine_calcul)
  !      if (.not.(in_domain)) then

  !         print*, 'Stopping the computation of column density because the boundary of the domain is reached'

  !         exit ray_propagation  ! ray is done 
  !      end if
  !      ! Ray moves to next cell : find it
  !      call whereIsPhotonGoing(domesh,icell,pos,icellnew,flagoutvol)
  !      ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
  !      ! -> give it an extra push untill it is out. 
  !      npush = 0
  !      do while (icell==icellnew)
  !         npush = npush + 1
  !         if (npush>100) then
  !            print*,'Too many pushes, npush>100 '
  !            stop
  !         endif

  !         ! hack
  !         do i=1,3
  !            pos(i) = pos(i) + k(i) * 1d7 * epsilon(pos(i))
  !         end do
  !         !!$ pos(1) = pos(1) + merge(-1.0d0,1.0d0,kray(1)<0.0d0) * epsilon(pos(1))
  !         !!$ pos(2) = pos(2) + merge(-1.0d0,1.0d0,kray(2)<0.0d0) * epsilon(pos(2))
  !         !!$ pos(3) = pos(3) + merge(-1.0d0,1.0d0,kray(3)<0.0d0) * epsilon(pos(3))
  !         ! kcah
          
  !         ! correct for periodicity
  !         do i=1,3
  !            if (pos(i) < 0.0d0) pos(i)=pos(i)+1.0d0
  !            if (pos(i) > 1.0d0) pos(i)=pos(i)-1.0d0
  !         enddo
  !         ! test that we are still in domain before calling WhereIsPhotonGoing... 
  !         in_domain = domain_contains_point(pos,domaine_calcul)
  !         if (.not.(in_domain)) then
  !            print*, 'Stopping the computation of column density because the boundary of the domain is reached'
  !            exit ray_propagation  ! ray is done 
  !         end if
  !         call whereIsPhotonGoing(domesh,icell,pos,icellnew,flagoutvol)
  !      end do
  !      if (npush > 1) print*,'WARNING : npush > 1 needed in module_gray_ray:propagate.'
       
  !      ! check if photon outside of cpu domain (flagoutvol)
  !      if(flagoutvol)then
  !         ! photon out of cpu domain -> we have a problem... 
  !         print*,'ERROR: photon out of CPU domain when it should not ... '
  !         stop
  !      endif
       
  !      ! else, new cell is in the cpu domain and photon goes to this new cell
  !      icell = icellnew
  !      if(domesh%son(icell)>=0)then
  !         print*,'ERROR: not a leaf cell',icell,flagoutvol
  !         print*,'This should not happen in module_ray (on single domains). '
  !         stop
  !      endif
  !      ileaf = - domesh%son(icell)
  !      ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
  !      ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
       
  !   end do ray_propagation
    

  ! end function column_density
  

  subroutine init_rays_from_file(file,rays)

    character(2000),intent(in)                           :: file
    type(ray_type),dimension(:),allocatable, intent(out) :: rays
    integer(kind=4)                                      :: i, n_rays

    ! read ICs
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='read')
    read(14) n_rays
    allocate(rays(n_rays))
    if (n_rays==0) return
    read(14) (rays(i)%ID,         i=1,n_rays)
    read(14) (rays(i)%x_em(1),    i=1,n_rays)
    read(14) (rays(i)%x_em(2),    i=1,n_rays)
    read(14) (rays(i)%x_em(3),    i=1,n_rays)
    read(14) (rays(i)%halo_ID,     i=1,n_rays)
    close(14)
    ! initialise other properties. 
    do i=1,n_rays
       rays(i)%dist  = 0.d0
       rays(i)%tau  = 0.d0
    enddo

    print*,minval(rays(:)%x_em(1)),maxval(rays(:)%x_em(1))
    print*,minval(rays(:)%x_em(2)),maxval(rays(:)%x_em(2))
    print*,minval(rays(:)%x_em(3)),maxval(rays(:)%x_em(3))
    
  end subroutine init_rays_from_file

  subroutine init_halos_from_file(file,halos)

    character(2000),intent(in)                            :: file
    type(halo_type),dimension(:),allocatable, intent(out) :: halos
    integer(kind=4)                                       :: i, n_halos

    print*,'Initialising halos from file (in module_gray_ray.f90)'
    ! read Halos (for escape fractions out of virial radii)
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='read')
    read(14) n_halos
    allocate(halos(n_halos))
    halos(1:n_halos)%domain%type='sphere'
    read(14) (halos(i)%ID,               i=1,n_halos)
    read(14) (halos(i)%domain%sp%radius,    i=1,n_halos)
    read(14) (halos(i)%domain%sp%center(1), i=1,n_halos)
    read(14) (halos(i)%domain%sp%center(2), i=1,n_halos)
    read(14) (halos(i)%domain%sp%center(3), i=1,n_halos)
    close(14)
    
  end subroutine init_halos_from_file


  subroutine dump_rays(file,rays)

    character(2000),intent(in)             :: file
    type(ray_type),dimension(:),intent(in) :: rays
    integer(kind=4)                        :: i,np

    np = size(rays)
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='write')
    write(14) np
    if(np.gt.0) then
       !write(14) (rays(i)%fesc,i=1,np)
       write(14) (rays(i)%ID,i=1,np)
    endif
    close(14)

  end subroutine dump_rays

end module module_gray_ray
