module py_cell_utils

  use hydro_reader

  public
  
contains

  subroutine count_cells(ramDir,ts,lmax,ncells,center,radius)

    implicit none

    character(1000),intent(in)  :: ramDir
    integer(kind=4),intent(in)  :: ts,lmax
    real(kind=8),intent(in)     :: center(3),radius
    integer(kind=4),intent(out) :: ncells
    
    character(1000) :: repo
    real(kind=8)    :: d2,radius2
    integer(kind=4) :: k,ilevel,l,ind,i,icell,lmax_local
    
    
    write(repo,'(a,a,i5.5)') trim(ramDir),'/output_',ts
    call read_info(repo)
    ! convention : if lmax <= 0, use nlevelmax
    if (lmax <= 0) then
       lmax_local = nlevelmax
    else
       lmax_local = lmax
    end if
    ! Loop over processor files and collect gas cells ----------------------
    icell = 0
    radius2 = radius*radius
    xmin = center(1) - radius
    xmax = center(1) + radius
    ymin = center(2) - radius
    ymax = center(2) + radius
    zmin = center(3) - radius
    zmax = center(3) + radius
    call init_reader(repo) ! Hydro reader
    do k=1,ncpu_read
       call read_cpu(k,repo,ncpu_read)
       write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
            ' Progress ',dble(k) / ncpu_read * 100,' % ',char(13)
       ! Loop over levels
       do ilevel=1,lmax_local
          call read_level(ilevel)
          if(ngrida .eq. 0) cycle
          do ind=1,twotondim ! Loop over cells
             do i=1,ngrida ! Loop over grids
                if(ref(i,ind) .and. ilevel < lmax_local) cycle  ! Cell is not a leaf and l < lmax
                d2 = 0
                do l=1,3
                   d2 = d2 + (xp(i,ind,l)-center(l))**2
                end do
                if (d2 < radius2) icell = icell + 1
             end do ! End loop over grids
          end do ! End loop over cells
       end do ! End loop over levels     
       call close_cpu()
    end do ! End loop over cpus

    ncells = icell 
    
    return
     
  end subroutine count_cells

  
  subroutine read_cells_hydro(ramDir,ts,lmax,cells,cell_pos,cell_level,ncells,fields,nfields,center,radius,readRT)

    implicit none

    character(1000),intent(in)  :: ramDir
    integer(kind=4),intent(in)  :: ts,lmax,nfields
    integer(kind=4),intent(in)  :: fields(nfields)
    integer(kind=4),intent(in)  :: ncells
    real(kind=8),intent(in)     :: center(3),radius
    real(kind=8),intent(out)    :: cells(ncells,nfields)
    real(kind=8),intent(out)    :: cell_pos(ncells,3)
    integer(kind=4),intent(out) :: cell_level(ncells)
    logical,intent(in)          :: readRT

    character(1000) :: repo
    real(kind=8)    :: d2,radius2
    integer(kind=4) :: k,ilevel,l,ind,i,icell,lmax_local
    
    write(repo,'(a,a,i5.5)') trim(ramDir),'/output_',ts
    call read_info(repo)
    ! convention : if lmax <= 0, use nlevelmax
    if (lmax <= 0) then
       lmax_local = nlevelmax
    else
       lmax_local = lmax
    end if

    ! Loop over processor files and collect gas cells ----------------------
    icell = 1
    radius2 = radius*radius
    xmin = center(1) - radius
    xmax = center(1) + radius
    ymin = center(2) - radius
    ymax = center(2) + radius
    zmin = center(3) - radius
    zmax = center(3) + radius
    call init_reader(repo,readRT) ! Hydro reader
    print*,'Number of read cpus: ',ncpu_read
    do k=1,ncpu_read
       call read_cpu(k,repo,ncpu_read)
       write (*, "(A, f5.2, A, A)", advance='no') &           ! Progress bar
            ' Progress ',dble(k) / ncpu_read * 100,' % ',char(13)
       ! Loop over levels
       do ilevel=1,lmax_local
          call read_level(ilevel)
          if(ngrida .eq. 0) cycle
          do ind=1,twotondim ! Loop over cells
             do i=1,ngrida ! Loop over grids
                if(ref(i,ind) .and. ilevel < lmax_local) cycle  ! Cell is not a leaf and l < lmax
                ! collect cells if in sphere
                d2 = 0
                do l=1,3
                   d2 = d2 + (xp(i,ind,l)-center(l))**2
                end do
                if (d2 < radius2) then
                   do l = 1,nfields
                      cells(icell,l) = var(i,ind,fields(l))
                   end do
                   do l = 1,3
                      cell_pos(icell,l) = xp(i,ind,l)
                   end do
                   cell_level(icell) = ilevel
                   icell = icell + 1
                end if
             end do ! End loop over grids
          end do ! End loop over cells
       end do ! End loop over levels     
       call close_cpu()
    end do ! End loop over cpus
    
    return
     
  end subroutine read_cells_hydro
  

  subroutine get_map_nxny(lmax,xmin,xmax,ymin,ymax,nx,ny)

    implicit none 
    integer(kind=4),intent(in)  :: lmax
    real(kind=8),intent(in)     :: xmin,xmax,ymin,ymax
    integer(kind=4),intent(out) :: nx,ny
    integer(kind=4)             :: n,imin,imax,jmin,jmax
    n=2**lmax
    imin=int(xmin*dble(n))+1
    imax=int(xmax*dble(n))
    jmin=int(ymin*dble(n))+1
    jmax=int(ymax*dble(n))
    nx = imax - imin + 1
    ny = jmax - jmin + 1

    return

  end subroutine get_map_nxny


    subroutine make_map(lmax,do_max,xmin,xmax,ymin,ymax,zmin,zmax,nleaf,nume,denu,cellx,celly,cellz,celllev,nx,ny,map,weight_map)

    ! project along z axis (change arguments to obtain the other projections). 

    implicit none

    integer(kind=4),intent(in) :: lmax
    logical,intent(in)         :: do_max
    real(kind=8),intent(in)    :: xmin,xmax,ymin,ymax,zmin,zmax
    integer(kind=4),intent(in) :: nleaf
    real(kind=8),intent(in)    :: nume(nleaf),denu(nleaf),cellx(nleaf),celly(nleaf),cellz(nleaf)
    integer(kind=4),intent(in) :: celllev(nleaf)
    integer(kind=4),intent(in)  :: nx,ny
    real(kind=8),intent(out)   :: map(nx,ny),weight_map(nx,ny)
    integer(kind=4)            :: ilevel,n_full,icell,ix,iy,iz,imin,imax,jmin,jmax,ndom,i,j
    real(kind=8)               :: dx, weight, oneovertwotolmax,xxmin,yymin
    type level
       integer::ilevel
       integer::ngrid
       real(KIND=8),dimension(:,:),pointer::map
       real(KIND=8),dimension(:,:),pointer::rho
       integer::imin
       integer::imax
       integer::jmin
       integer::jmax
       integer::kmin
       integer::kmax
    end type level
    type(level),dimension(1:100)::grid

    map = 0.0d0

    print*,minval(nume),maxval(nume)
    
    !-----------------------
    ! Map parameters
    !-----------------------
    do ilevel=1,lmax
       n_full=2**ilevel
       grid(ilevel)%imin=int(xmin*dble(n_full))+1
       grid(ilevel)%imax=int(xmax*dble(n_full))+1
       grid(ilevel)%jmin=int(ymin*dble(n_full))+1
       grid(ilevel)%jmax=int(ymax*dble(n_full))+1
       grid(ilevel)%kmin=int(zmin*dble(n_full))+1
       grid(ilevel)%kmax=int(zmax*dble(n_full))+1
       allocate(grid(ilevel)%map(grid(ilevel)%imin:grid(ilevel)%imax,grid(ilevel)%jmin:grid(ilevel)%jmax))
       if (.not. do_max) allocate(grid(ilevel)%rho(grid(ilevel)%imin:grid(ilevel)%imax,grid(ilevel)%jmin:grid(ilevel)%jmax))
       grid(ilevel)%map(:,:)=0.0
       if (.not. do_max) grid(ilevel)%rho(:,:)=0.0
    end do

    ! loop over cells to construct ilevel maps 
    do icell = 1,nleaf
       ilevel = celllev(icell)
       n_full = 2**ilevel
       dx     = 1./n_full
       ix=int(cellx(icell)*dble(n_full))+1
       iy=int(celly(icell)*dble(n_full))+1
       iz=int(cellz(icell)*dble(n_full))+1 
       if(    ix>=grid(ilevel)%imin.and.iy>=grid(ilevel)%jmin.and.iz>=grid(ilevel)%kmin.and.&            
            & ix<=grid(ilevel)%imax.and.iy<=grid(ilevel)%jmax.and.iz<=grid(ilevel)%kmax) then              
          if(do_max)then
             grid(ilevel)%map(ix,iy)=max(grid(ilevel)%map(ix,iy),nume(icell))
          else
             weight=(min(cellz(icell)+dx/2.,zmax)-max(cellz(icell)-dx/2.,zmin))
             weight=min(1.0d0,max(weight,0.0d0)) / (zmax - zmin)
             grid(ilevel)%map(ix,iy)=grid(ilevel)%map(ix,iy)+nume(icell)*denu(icell)*weight
             grid(ilevel)%rho(ix,iy)=grid(ilevel)%rho(ix,iy)+denu(icell)*weight
          endif
       endif
    end do

    ! flatten maps
    n_full=2**lmax
    n_full=2**lmax
    imin=int(xmin*dble(n_full))+1
    imax=int(xmax*dble(n_full))
    jmin=int(ymin*dble(n_full))+1
    jmax=int(ymax*dble(n_full))
    oneovertwotolmax = 1.0d0/2**lmax
    do ix=imin,imax
       xxmin=(ix-0.5)*oneovertwotolmax
       do iy=jmin,jmax
          yymin=(iy-0.5)*oneovertwotolmax
          do ilevel=1,lmax-1
             ndom=2**ilevel
             i=int(xxmin*ndom)+1
             j=int(yymin*ndom)+1
             if(do_max) then
                grid(lmax)%map(ix,iy)=max(grid(lmax)%map(ix,iy), &
                     & grid(ilevel)%map(i,j))
             else
                grid(lmax)%map(ix,iy)=grid(lmax)%map(ix,iy) + &
                     & grid(ilevel)%map(i,j)
                grid(lmax)%rho(ix,iy)=grid(lmax)%rho(ix,iy) + &
                     & grid(ilevel)%rho(i,j)
             endif
          end do
       end do
    end do
    if(do_max)then
       map=grid(lmax)%map(imin:imax,jmin:jmax)
       weight_map = 0.
    else
       weight_map = grid(lmax)%rho(imin:imax,jmin:jmax)
       map=grid(lmax)%map(imin:imax,jmin:jmax)/weight_map
!!$       do ix=imin,imax
!!$          do iy=jmin,jmax
!!$             if (weight_map(ix,iy) /= 0.) map(ix,iy) = map(ix,iy) / weight_map(ix,iy)
!!$          end do
!!$       end do
    endif


    return

  end subroutine make_map

  !***********************************************************************
  subroutine get_outflow_inflow_rates(ramDir,ts,hCenter,hVel,hRvir &
       ,inRate,outRate)

  ! Calculate rates of inflow and outflow (Msun/yr) for (one) given halo
  ! position and radius (both in box units)
  ! ramDir  => ramses run directory  
  ! ts      => snapshot number
  ! hCenter => Halo center (code units, i.e. 0 to 1)
  ! hRvir   => Distance at which to calculate rates (code units)
  ! hVel    => Halo velocity (code units)
  ! inRate  <= Inflow rate (Msun/yr)
  ! outRate <= Outflow rate (Msun/yr)
  !-----------------------------------------------------------------------
    implicit none
    character(1000),intent(in):: ramDir
    integer(kind=4),intent(in)::ts
    real(kind=8),intent(in):: hCenter(3),hRvir,hVel(3)
    real(kind=8),intent(out):: inRate,outRate
    character(1000) :: repo
    integer(kind=4)::ncells,nFields,ic
    real(kind=8),allocatable:: cells(:,:),cell_pos(:,:)
    integer(kind=4),allocatable:: cell_level(:),fields(:)
    real(kind=8)::rad,rMin,rMax,mass,massOut,massIn,mom,momOut,momIn
    real(kind=8)::vol,volTot,scale,unit_Mrate,area,vsh(3)
    logical::readRT=.false.
  !-----------------------------------------------------------------------
    write(repo,'(a,a,i5.5)') trim(ramDir),'/output_',ts
    call read_info(repo)    
    nFields = 4 ! Density, vx, vy, vz
    rMin=hRvir*0.9
    rMax=hRvir*1.1
    scale=1d0
    call count_cells(ramDir,ts,nlevelmax,ncells,hCenter,rMax)
    allocate(cells(ncells,nfields),cell_pos(ncells,3),cell_level(ncells))
    allocate(fields(nFields))
    fields=(/ 1,2,3,4 /)
    print*,'Fields = ',fields
    print*,'Levelmax = ',nlevelmax
    call read_cells_hydro(ramDir,ts,nlevelmax,cells,cell_pos,cell_level &
         ,ncells,fields,nFields,hCenter,rMax,readRT)
    massOut=0d0 ; momOut=0d0 ; massIn=0d0 ; momIn=0d0 ; volTot=0d0
    print*,hCenter
    do ic=1,ncells
       ! Add momenta from all cells at 0.9 rVir to 1.1 rVir
       rad = sqrt((cell_pos(ic,1)-hCenter(1))**2 &
                 +(cell_pos(ic,2)-hCenter(2))**2 &
                 +(cell_pos(ic,3)-hCenter(3))**2)
       if( rad > rMin .and. rad < rMax ) then ! Cell is in rVir shell
          ! Shift the cell gas velocity into the rest-frame of the halo
          vsh =cells(ic,2:4)-hVel
          vol = (2d0**(cell_level(ic))*scale)**3 ! Cell volume (code units)
          mass = vol * cells(ic,1) ! Gas mass in cell (code units)
          mom  = mass * SUM(vsh*cell_pos(ic,:)) / rad ! Radial momentum (code u)
          volTot  = volTot + vol ! Total volume added to shell
          if(mom .gt.0d0) then
             massOut = massOut + mass
             momOut  = momOut  + mom
          else
             massIn = massIn + mass
             momIn  = momIn  + mom
          endif
       endif
    end do

    unit_Mrate = unit_m/2d33/unit_t*yr2s  ! Msun/yr
    !unit_speed = unit_l/unit_t/1d5 ! km/s
    !unit_rad = unit_l*scale/3.08d21/r1_kpc ! rvir
    area = 4d0*3.14*hRvir**2 ! code units
    inRate   = momIn /volTot*area * scale**2 * unit_Mrate
    outRate  = momOut/volTot*area * scale**2 * unit_Mrate
    !inSpeed  = momIn/massIn * unit_speed 
    !outSpeed = momOut/massOut * unit_speed 

  end subroutine get_outflow_inflow_rates



  
end module py_cell_utils
