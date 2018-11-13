module module_domain

  use module_utils, only: path
  
  implicit none

  public
  
  type shell
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: r_inbound,r_outbound
  end type shell

  type cube
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: size            ! convention: size is the full size of the cube, corners are x+-size/2
  end type cube

  type sphere
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: radius
  end type sphere

  type slab                                   ! infinite slab in xy direction (make use of periodic boundaries)
     real(kind=8)              :: zc          ! position of the slab in the z direction
     real(kind=8)              :: thickness   ! thickness of the slab in the z direction 
  end type slab

  type domain
     character(10) :: type  ! one of the possible shapes ('shell', 'cube', 'sphere', 'slab')
     type(shell)   :: sh
     type(cube)    :: cu
     type(sphere)  :: sp
     type(slab)    :: sl
  end type domain
  

contains
  
  !--------------------------------------------------------------------------------------------------
  ! domain constructors
  !--------------------------------------------------------------------------------------------------

  subroutine domain_constructor_from_scratch(dom,type,xc,yc,zc,r,r_inbound,r_outbound,size,thickness)

    implicit none
    character(10),intent(in)         :: type
    real(kind=8),intent(in),optional :: xc,yc,zc
    real(kind=8),intent(in),optional :: r                     ! parameters for sphere
    real(kind=8),intent(in),optional :: r_inbound,r_outbound  ! parameters for shell
    real(kind=8),intent(in),optional :: size                  ! parameters for cube
    real(kind=8),intent(in),optional :: thickness             ! parameters for slab
    type(domain),intent(out)         :: dom
    logical                          :: ok

    select case(trim(type))

    case('sphere')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a sphere domain are missing...'
          stop
       endif
       dom%type=type
       dom%sp%center(1)=xc
       dom%sp%center(2)=yc
       dom%sp%center(3)=zc
       dom%sp%radius=r

    case('shell')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r_inbound).and.present(r_outbound)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a shell domain are missing...'
          stop
       endif
       dom%type=type
       dom%sh%center(1)=xc
       dom%sh%center(2)=yc
       dom%sh%center(3)=zc
       dom%sh%r_inbound=r_inbound
       dom%sh%r_outbound=r_outbound

    case('cube')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(size)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a cube domain are missing...'
          stop
       endif
       dom%type=type
       dom%cu%center(1)=xc
       dom%cu%center(2)=yc
       dom%cu%center(3)=zc
       dom%cu%size=size


    case('slab')
       ! check if optional argument required for sphere are present
       ok = present(zc).and.present(thickness)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a slab domain are missing...'
          stop
       endif
       dom%type=type
       dom%sl%zc=zc
       dom%sl%thickness=thickness

    case default
       print *,'ERROR: type not defined ',trim(type)
       stop
    end select

    return

  end subroutine domain_constructor_from_scratch


  
  subroutine domain_constructor_from_file(filename,dom)
    ! read a domain file (filename) and initialise domain dom from it

    implicit none
    character(2000),intent(in) :: filename
    type(domain),intent(inout) :: dom

    open(unit=13, file=trim(filename), status='old', form='formatted', action='read')

    call read_domain(unit=13,dom=dom)

    close(13)

  end subroutine domain_constructor_from_file


  !--------------------------------------------------------------------------------------------------
  ! public use
  !--------------------------------------------------------------------------------------------------

  subroutine select_in_domain(dom,n,xp,indsel)

    implicit none
    integer(kind=4),intent(in)                           :: n
    type(domain),intent(in)                              :: dom 
    real(kind=8),dimension(1:n,1:3),intent(in)           :: xp
    integer(kind=4),dimension(:),allocatable,intent(out) :: indsel
    integer(kind=4)                                      :: i,ii,nsel
    integer(kind=4),dimension(:),allocatable             :: tmpi
    real(kind=8)                                         :: dd
    
    select case(trim(dom%type))
    
    case('sphere')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          dd = (xp(i,1)-dom%sp%center(1))**2 + (xp(i,2)-dom%sp%center(2))**2 + (xp(i,3)-dom%sp%center(3))**2
          if(sqrt(dd)<=dom%sp%radius)then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)

    case('shell')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          dd = (xp(i,1)-dom%sh%center(1))**2 + (xp(i,2)-dom%sh%center(2))**2 + (xp(i,3)-dom%sh%center(3))**2
          if((sqrt(dd)>=dom%sh%r_inbound).and.(sqrt(dd)<dom%sh%r_outbound))then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)


    case('cube')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          if((abs(xp(i,1)-dom%cu%center(1)) <= dom%cu%size*0.5d0).and. & 
             (abs(xp(i,2)-dom%cu%center(2)) <= dom%cu%size*0.5d0).and. &
             (abs(xp(i,3)-dom%cu%center(3)) <= dom%cu%size*0.5d0))then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)


    case('slab')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          if(abs(xp(i,3)-dom%sl%zc) <= dom%sl%thickness*0.5d0)then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)


    case default
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine select_in_domain



  subroutine read_domain(unit,dom)

    integer(kind=4),intent(in)  :: unit
    type(domain),intent(out)    :: dom 

    read(unit,*) dom%type
    select case(trim(dom%type))
    case('sphere')
       read(unit,*) dom%sp%center(:)
       read(unit,*) dom%sp%radius
    case('shell')
       read(unit,*) dom%sh%center(:)
       read(unit,*) dom%sh%r_inbound
       read(unit,*) dom%sh%r_outbound
    case('cube')
       read(unit,*) dom%cu%center(:)
       read(unit,*) dom%cu%size
    case('slab')
       read(unit,*) dom%sl%zc
       read(unit,*) dom%sl%thickness
    case default
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine read_domain



  subroutine read_domain_bin(unit,dom)

    integer(kind=4),intent(in)  :: unit
    type(domain),intent(out)    :: dom 

    read(unit) dom%type
    select case(trim(dom%type))
    case('sphere')
       read(unit) dom%sp%center(:)
       read(unit) dom%sp%radius
    case('shell')
       read(unit) dom%sh%center(:)
       read(unit) dom%sh%r_inbound
       read(unit) dom%sh%r_outbound
    case('cube')
       read(unit) dom%cu%center(:)
       read(unit) dom%cu%size
    case('slab')
       read(unit) dom%sl%zc
       read(unit) dom%sl%thickness
    case default
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine read_domain_bin



  subroutine dump_domain_bin(unit,dom)

    integer(kind=4),intent(in) :: unit
    type(domain),intent(in)    :: dom 

    write(unit) dom%type
    select case(trim(dom%type))
    case('sphere')
       write(unit) dom%sp%center(:)
       write(unit) dom%sp%radius
    case('shell')
       write(unit) dom%sh%center(:)
       write(unit) dom%sh%r_inbound
       write(unit) dom%sh%r_outbound
    case('cube')
       write(unit) dom%cu%center(:)
       write(unit) dom%cu%size
    case('slab')
       write(unit) dom%sl%zc
       write(unit) dom%sl%thickness
    case default 
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine dump_domain_bin

  
  
  subroutine dump_domain_form(unit,dom)

    integer(kind=4),intent(in) :: unit
    type(domain),intent(in)    :: dom 

    write(unit,'(a)') trim(dom%type)
    select case(trim(dom%type))
    case('sphere')
       write(unit,'(3(f14.10))') dom%sp%center(:)
       write(unit,'(f14.10)') dom%sp%radius
    case('shell')
       write(unit,'(3(f14.10))') dom%sh%center(:)
       write(unit,'(f14.10)') dom%sh%r_inbound
       write(unit,'(f14.10)') dom%sh%r_outbound
    case('cube')
       write(unit,'(3(f14.10))') dom%cu%center(:)
       write(unit,'(f14.10)') dom%cu%size
    case('slab')
       write(unit,'(f14.10)') dom%sl%zc
       write(unit,'(f14.10)') dom%sl%thickness
    case default 
       print *,'ERROR: type not defined ',trim(dom%type)
       stop
    end select

  end subroutine dump_domain_form



  subroutine domain_write_file(file,dom)

    character(*),intent(in) :: file
    type(domain),intent(in) :: dom 

    !!!write(file ,'(a,a)') trim(fichier),'.dom'
    open(unit=14, file=trim(file), status='unknown', form='formatted', action='write')
    call dump_domain_form(unit=14,dom=dom)
    close(14)

  end subroutine domain_write_file



  function domain_contains_point(x,dom)
    ! -> returns T/F if point xyz is in domain dom.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    logical                              :: domain_contains_point
    real(kind=8)                         :: rr
    domain_contains_point=.false.
    select case(trim(dom%type))
    case('sphere')
       rr = (x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2
       if(rr<dom%sp%radius*dom%sp%radius)domain_contains_point=.true.
    case('shell')
       rr = (x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2
       if((rr>dom%sh%r_inbound*dom%sh%r_inbound).and.(rr<dom%sh%r_outbound*dom%sh%r_outbound))domain_contains_point=.true.
    case('cube')
       if((abs(x(1)-dom%cu%center(1)) < dom%cu%size*0.5d0).and. & 
            (abs(x(2)-dom%cu%center(2)) < dom%cu%size*0.5d0).and. &
            (abs(x(3)-dom%cu%center(3)) < dom%cu%size*0.5d0))domain_contains_point=.true.
    case('slab')
       if(abs(x(3)-dom%sl%zc) < dom%sl%thickness*0.5d0)domain_contains_point=.true.
    end select
    return
  end function domain_contains_point

  

  function domain_contains_cell(x,dx,dom)
    ! -> returns T/F if the full cell at x, of size dx, is in domain dom.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8),intent(in)              :: dx
    logical                              :: domain_contains_cell
    real(kind=8)                         :: rr
    real(kind=8),parameter :: sqrt3over2 = sqrt(3.0d0)*0.5d0
    
    domain_contains_cell=.false.

    select case(trim(dom%type))

    case('sphere')
       
       rr = sqrt((x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2)
       rr = rr + dx*sqrt3over2
       if (rr < dom%sp%radius) domain_contains_cell=.true.
      
    case('shell')
       
       rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       if(( (rr-dx*sqrt3over2)>dom%sh%r_inbound) .and. ((rr+dx*sqrt3over2)<dom%sh%r_outbound) ) then
          domain_contains_cell=.true.
       end if

    case('cube')
       
       if((x(1)+dx*0.5d0 < dom%cu%center(1)+dom%cu%size*0.5d0).and. &
          (x(1)-dx*0.5d0 > dom%cu%center(1)-dom%cu%size*0.5d0).and. &
          (x(2)+dx*0.5d0 < dom%cu%center(2)+dom%cu%size*0.5d0).and. &
          (x(2)-dx*0.5d0 > dom%cu%center(2)-dom%cu%size*0.5d0).and. &
          (x(3)+dx*0.5d0 < dom%cu%center(3)+dom%cu%size*0.5d0).and. &
          (x(3)-dx*0.5d0 > dom%cu%center(3)-dom%cu%size*0.5d0)) domain_contains_cell=.true.

    case('slab')
       
       if((x(3)+dx*0.5d0 < dom%sl%zc+dom%sl%thickness*0.5d0).and. &
          (x(3)-dx*0.5d0 > dom%sl%zc-dom%sl%thickness*0.5d0)) domain_contains_cell=.true.

    end select

    return
  end function domain_contains_cell



  function get_my_new_domain(x,liste_domaines)
    ! Return the domain from list of domains “liste_domaines” in which point x is most deeply embedded. 
    ! (i.e. such that the distance of x to domain border is the largest). 
    type(domain),intent(in),dimension(:) :: liste_domaines
    real(kind=8),dimension(3),intent(in) :: x
    integer(kind=4)                      :: get_my_new_domain, ndom, i, imax
    real(kind=8)                         :: dmax,db

    ndom = size(liste_domaines)
    dmax = 0.0d0
    imax = 0
    do i=1,ndom
       db = domain_distance_to_border(x,liste_domaines(i))
       ! if db < 0 point outside domain
       if (db>=dmax)then
          imax = i
          dmax = db
       endif
    enddo
    if(dmax<0.0d0 .or. imax==0)then
       print *,'ERROR: problem with get_my_new_domain'
       stop
    endif
    get_my_new_domain = imax

    return
  end function get_my_new_domain


  
  function domain_distance_to_border(x,dom)
    ! return distance of point x to the closest border of domain dom
    ! convention: negative distance means outside domain
    real(kind=8),dimension(3),intent(in) :: x
    type(domain),intent(in)              :: dom
    real(kind=8)                         :: domain_distance_to_border, rr, ddx, ddy, ddz

    select case(trim(dom%type))

    case('sphere')
       rr = sqrt((x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2)
       domain_distance_to_border = dom%sp%radius - rr

    case('shell')
       rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       domain_distance_to_border = min((rr-dom%sh%r_inbound),(dom%sh%r_outbound-rr))

    case('cube')
       ddx = min((dom%cu%center(1)+dom%cu%size*0.5d0-x(1)), (x(1)-dom%cu%center(1)+dom%cu%size*0.5d0))
       ddy = min((dom%cu%center(2)+dom%cu%size*0.5d0-x(2)), (x(2)-dom%cu%center(2)+dom%cu%size*0.5d0))
       ddz = min((dom%cu%center(3)+dom%cu%size*0.5d0-x(3)), (x(3)-dom%cu%center(3)+dom%cu%size*0.5d0))
       if((ddx>=0.0d0).and.(ddy>=0.0d0).and.(ddz>=0.0d0))then
          ! inside domain
          domain_distance_to_border = min(ddx,ddy,ddz)
       else
          ! outside domain
          domain_distance_to_border = sqrt((min(0.0d0,ddx))**2 + (min(0.0d0,ddy))**2 + (min(0.0d0,ddz))**2)  
       endif

    case('slab')
       domain_distance_to_border = min((dom%sl%zc+dom%sl%thickness*0.5d0-x(3)), (x(3)-dom%sl%zc+dom%sl%thickness*0.5d0)) 

    end select

    return
  end function domain_distance_to_border



  function domain_distance_to_border_along_k(x,k,dom)
    ! return the distance of point x to the closest border of domain dom along propagation vector k.
    ! inputs:
    ! - x : position vector, in simulation-box units
    ! - k : direction vector, normalized
    ! - dom : a domain
    
    real(kind=8),dimension(3),intent(in) :: x, k
    type(domain),intent(in)              :: dom
    real(kind=8)                         :: domain_distance_to_border_along_k
    ! variables for the spherical case
    real(kind=8) :: b, c, delta, dx, dy, dz 
    ! variables for the shell case
    real(kind=8) :: t1,t2,tin,tout
    ! variables for the cube case
    real(kind=8),dimension(3) :: x_dom

    ! point x should be inside domain dom
    if(.not.(domain_contains_point(x,dom)))then
       print*,'ERROR: point x outside domain dom in domain_distance_to_border_along_k...'
       stop
    endif
    
    select case(trim(dom%type))
       
    case('sphere')
       
       dx = x(1) - dom%sp%center(1)
       dy = x(2) - dom%sp%center(2)
       dz = x(3) - dom%sp%center(3)
       b = 2.d0 * ( k(1)*dx + k(2)*dy +  k(3)*dz )
       c = dx*dx + dy*dy + dz*dz - dom%sp%radius*dom%sp%radius
       delta = b*b - 4.0d0*c
       if (delta <= 0.0d0) then
          print*,'ERROR: pb with domain_distance_to_border_along_k' ! This can only mean that x is out of the domain ... 
          stop
       end if
       domain_distance_to_border_along_k = (-b + sqrt(delta))*0.5d0 ! select the only positive solution

    case('shell')

       dx = x(1) - dom%sh%center(1)
       dy = x(2) - dom%sh%center(2)
       dz = x(3) - dom%sh%center(3)
       ! inner shell intersection
       b = 2.0d0 * ( k(1)*dx + k(2)*dy +  k(3)*dz )
       c = dx*dx + dy*dy + dz*dz - dom%sh%r_inbound*dom%sh%r_inbound
       delta = b*b - 4.0d0*c
       tin   = 2.0d0 ! larger than box size 
       if (delta >= 0.0d0) then
          t1 =  (-b + sqrt(delta))*0.5d0
          t2 = (-b - sqrt(delta))*0.5d0
          ! select smallest positive solution
          if (t1 > 0 .and. t1 < tin) tin = t1
          if (t2 > 0 .and. t2 < tin) tin = t2
       end if
       ! outer shell intersection
       b = 2.0d0 * ( k(1)*dx + k(2)*dy +  k(3)*dz )
       c = dx*dx + dy*dy + dz*dz - dom%sh%r_outbound*dom%sh%r_outbound
       delta = b*b - 4.0d0*c
       tout  = 2.0d0 ! larger than box size 
       if (delta >= 0.0d0) then
          t1 =  (-b + sqrt(delta))*0.5d0
          t2 = (-b - sqrt(delta))*0.5d0
          ! select smallest positive solution
          if (t1 > 0 .and. t1 < tout) tout = t1
          if (t2 > 0 .and. t2 < tout) tout = t2
       end if
       ! Select smallest distance
       domain_distance_to_border_along_k = min(tin,tout)
       
    case('slab')

       if(k(3)<0.0d0)then
          domain_distance_to_border_along_k = (dom%sl%zc-dom%sl%thickness*0.5d0-x(3))/k(3)
       else
          domain_distance_to_border_along_k = (dom%sl%zc+dom%sl%thickness*0.5d0-x(3))/k(3)
       endif
       if (domain_distance_to_border_along_k < 0.0d0)then
          print *,'ERROR: pb with distance to border along k, slab case, point outside domain...'
          stop
       endif

    case('cube')  

       ! get position relative to the domain
       x_dom = (x - dom%cu%center)/dom%cu%size + 0.5d0
       if((x_dom(1) < 0.0d0).or.(x_dom(1)>1.0d0).or.&
            (x_dom(2) < 0.0d0).or.(x_dom(2)>1.0d0).or.&
            (x_dom(3) < 0.0d0).or.(x_dom(3)>1.0d0))then
          ! x is outside the domain
          print *,'ERROR: pb with distance to border along k, cube case, point outside domain...'
          stop
       endif
       domain_distance_to_border_along_k = path(x_dom,k) * dom%cu%size
       
    case default
       print *,'ERROR: domain type not defined ',trim(dom%type)
       stop

    end select
       
    return
    
  end function domain_distance_to_border_along_k

  
end module module_domain

