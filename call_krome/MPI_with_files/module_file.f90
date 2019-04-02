module module_file

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi_mine
  use module_krome
  use module_ramses
  use module_cell

  implicit none

  type file
     integer(kind=4)               :: ID
     integer(kind=4)               :: status  ! 0 if not done, 1 if done
  end type file

  integer(kind=4)                  :: ncpu, ncell, ncoarse, ngridmax!, nleaf_tot
  real(kind=8),allocatable         :: xg(:,:)      ! grids position
  real(kind=8),allocatable         :: var(:,:)
  integer,allocatable              :: headl(:,:),taill(:,:),numbl(:,:),numbtot(:,:)
  integer,allocatable              :: headb(:,:),tailb(:,:),numbb(:,:)
  integer,allocatable              :: nbor(:,:)    ! neighboring father cells
  integer,allocatable              :: son(:)       ! sons grids
  integer,allocatable              :: cpu_map(:)  ! domain decomposition
  integer,allocatable              :: next(:)      ! next grid in list
  real(KIND=8),dimension(1:3)      :: xbound=(/0d0,0d0,0d0/) 
  integer,parameter                :: twotondim = 8, ndim = 3, twondim = 6
  logical                          :: first=.true., first2=.true., first3=.true.

  public :: init_files, compute_file

contains  


  subroutine init_files(repository, snapnum, filegrid)


    implicit none

    character(2000),intent(in)              :: repository
    integer(kind=4),intent(in)              :: snapnum
    type(file),allocatable,intent(inout)    :: filegrid(:)
    integer(kind=4)                         :: i

    ncpu = get_ncpu(repository, snapnum)

    allocate(filegrid(ncpu))
    do i=1,ncpu
       filegrid(i)%ID = i
       filegrid(i)%status = 0
    end do
    !nleaf_tot = 0

  end subroutine init_files



  subroutine compute_file(repository, snapnum, repo_restart, snap_restart, icpu, output_path)

    implicit none

    character(2000),intent(in)              :: repository, repo_restart, output_path
    integer(kind=4),intent(in)              :: snapnum, snap_restart, icpu
    real(kind=8),allocatable                :: ramses_var(:,:)
    real(kind=8),allocatable,dimension(:,:) :: cells_rt, fractions
    real(kind=8),allocatable,dimension(:)   :: nH, Tgas, mets, nHI
    integer(kind=4)                         :: nSEDgroups, nvar, nleaf
    character(1000)                         :: nomfich, nomfich2
    logical                                 :: exist, exist2


    write(nomfich,'(a,a,a,a,a,i5.5,a,i5.5)') trim(output_path),'/',trim(element_names(elements(n_elements))),trim(roman_num(n_ions(n_elements))),'_',snapnum,'.out',icpu
    inquire(file=nomfich, exist=exist)

    if(exist) then
       print*, 'file ', icpu, ' already computed, passing to the next'
    else

       ncpu = get_ncpu(repository,snapnum)

       write(nomfich2,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repo_restart),'/output_',snap_restart,'/rt_',snap_restart,'.out',icpu
       inquire(file=nomfich2, exist=exist2)
       
       if(exist2) then
          if(first2 .and. rank==1) print*, 'Using restart output'
          first2=.false.
          call read_hydro_mine_restart(repository, snapnum, repo_restart, snap_restart, icpu, nvar, nleaf, ramses_var)
       else
          if(first) print*, 'No restart file found, so continuing without restart file'
          first = .false.
          call read_hydro_mine(repository, snapnum, icpu, nvar, nleaf, ramses_var)
       end if
       
       call init_cells(repository, snapnum, nvar, nleaf, ramses_var)

       call compute_cells()

       call write_ion_files(snapnum, icpu, output_path)
       
    end if

  end subroutine compute_file


  subroutine read_hydro_mine(repository, snapnum, icpu, nvar, nleaf, ramses_var)

    implicit none

    character(2000),intent(in)             :: repository
    integer(kind=4),intent(in)             :: snapnum, icpu
    real(kind=8),intent(out),allocatable   :: ramses_var(:,:)
    integer(kind=4),intent(out)            :: nvar,nleaf
    character(1000)                        :: nomfich
    integer(kind=4)                        :: i,nlevelmax,nboundary,ix,iy,iz,ind,ilevel,ibound,ncache,istart,ivar,iskip,igrid,nvarH,nvarRT,ileaf,icell
    real(kind=8),allocatable               :: xc(:,:),xx(:)
    integer(kind=4),allocatable            :: ind_grid(:)
    real(kind=8)                           :: dx

    call init_amr(repository, snapnum, icpu)

    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    read(10)
    read(10)nvarH
    read(10)
    read(10)nlevelmax
    read(10)nboundary
    read(10)

    ! Open RT file and get nvarRT
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
    open(unit=12,file=nomfich,status='old',form='unformatted')
    read(12)
    read(12)nvarRT
    read(12)
    read(12)
    read(12)
    read(12)

    allocate(ramses_var(1:ncell,1:nvarH+nvarRT))
    nvar = nvarH + nvarRT
    if(rank==1 .and. first3) print*, 'nvar : ', nvar
    first3=.false.

    allocate(xc(1:twotondim,1:ndim))

    do ilevel=1,nlevelmax

       dx=0.5d0**ilevel
       do ind=1,8
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

          read(12)
          read(12)

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
                      ramses_var(ind_grid(i)+iskip,ivar) = xx(i)
                   end do
                end do

                do ivar=1,nvarRT
                   read(12) xx
                   if (ibound > ncpu) cycle  ! dont bother with boundaries
                   do i = 1, ncache
                      ramses_var(ind_grid(i)+iskip,ivar+nvarH) = xx(i)
                   end do
                end do
               
             end do
             deallocate(ind_grid,xx)
          end if
       end do
    end do
    deallocate(xc)
    close(10)
    close(12)

    allocate(var(nvar,ncell))
    var = 0d0
    ileaf=0
    do icell=1,ncell
       if(son(icell)==0 .and. cpu_map(icell)==icpu) then
          ileaf = ileaf+1
          do ivar=1,nvar
             var(ivar,ileaf) = ramses_var(icell,ivar)
          end do
       end if
    end do

    nleaf = ileaf
    deallocate(ramses_var)
    allocate(ramses_var(nvar,nleaf))
    do i=1,nleaf
       ramses_var(:,i) = var(:,i)
    end do
    deallocate(var)

    call clear_amr

    return

  end subroutine read_hydro_mine



  subroutine read_hydro_mine_restart(repository, snapnum, repo_restart, snap_restart, icpu, nvar, nleaf, ramses_var)

    implicit none

    character(2000),intent(in)             :: repository, repo_restart
    integer(kind=4),intent(in)             :: snapnum, icpu, snap_restart
    real(kind=8),intent(out),allocatable   :: ramses_var(:,:)
    integer(kind=4),intent(out)            :: nvar,nleaf
    character(1000)                        :: nomfich
    integer(kind=4)                        :: i,nlevelmax,nboundary,ix,iy,iz,ind,ilevel,ibound,ncache,istart,ivar,iskip,igrid,nvarH,nvarRT,nvarRT_restart,ileaf,icell
    real(kind=8),allocatable               :: xc(:,:),xx(:)
    integer(kind=4),allocatable            :: ind_grid(:)
    real(kind=8)                           :: dx

    call init_amr(repository, snapnum, icpu)

    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    read(10)
    read(10)nvarH
    read(10)
    read(10)nlevelmax
    read(10)nboundary
    read(10)

    ! Open RT file and get nvarRT
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
    open(unit=12,file=nomfich,status='old',form='unformatted')
    read(12)
    read(12)nvarRT
    read(12)
    read(12)
    read(12)
    read(12)

    !Same for restart
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repo_restart),'/output_',snap_restart,'/rt_',snap_restart,'.out',icpu
    open(unit=13,file=nomfich,status='old',form='unformatted')
    read(13)
    read(13)nvarRT_restart
    read(13)
    read(13)
    read(13)
    read(13)

    allocate(ramses_var(1:ncell,1:nvarH+nvarRT+nvarRT_restart))
    nvar = nvarH + nvarRT + nvarRT_restart
    if(rank==1 .and. first3) print*, 'nvar : ', nvar
    first3=.false.

    allocate(xc(1:twotondim,1:ndim))

    do ilevel=1,nlevelmax

       dx=0.5d0**ilevel
       do ind=1,8
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

          read(12)
          read(12)
          read(13)
          read(13)

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
                      ramses_var(ind_grid(i)+iskip,ivar) = xx(i)
                   end do
                end do

                !restart
                do ivar=1,nvarRT_restart
                   read(13) xx
                   if (ibound > ncpu) cycle  ! dont bother with boundaries
                   do i = 1, ncache
                      ramses_var(ind_grid(i)+iskip,ivar+nvarH) = xx(i)
                   end do
                end do

                do ivar=1,nvarRT
                   read(12) xx
                   if (ibound > ncpu) cycle  ! dont bother with boundaries
                   do i = 1, ncache
                      ramses_var(ind_grid(i)+iskip,ivar+nvarH+nvarRT_restart) = xx(i)
                   end do
                end do
                
             end do
             deallocate(ind_grid,xx)
          end if
       end do
    end do
    deallocate(xc)
    close(10)
    close(12)
    close(13)

    allocate(var(nvar,ncell))
    var = 0d0
    ileaf=0
    do icell=1,ncell
       if(son(icell)==0 .and. cpu_map(icell)==icpu) then
          ileaf = ileaf+1
          do ivar=1,nvar
             var(ivar,ileaf) = ramses_var(icell,ivar)
          end do
       end if
    end do

    nleaf = ileaf
    deallocate(ramses_var)
    allocate(ramses_var(nvar,nleaf))
    do i=1,nleaf
       ramses_var(:,i) = var(:,i)
    end do
    deallocate(var)

    call clear_amr

    return

  end subroutine read_hydro_mine_restart



  subroutine init_amr(repository, snapnum, icpu)

    implicit none

    character(2000),intent(in)     :: repository
    integer(kind=4),intent(in)     :: snapnum, icpu
    integer(kind=4)                :: get_ncell
    character(1000)                :: nomfich 
    integer,allocatable            :: ind_grid(:),iig(:),grid(:)
    real(kind=8),allocatable       :: xxg(:)
    logical                        :: ok
    integer(kind=4)                :: i,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)                :: ilevel,ncache,ibound,idim,ind,iskip


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
    allocate( headl(1:ncpu,1:nlevelmax),taill(1:ncpu,1:nlevelmax), &
         & numbl(1:ncpu,1:nlevelmax),numbtot(1:10,1:nlevelmax), &
         & headb(1:nboundary,1:nlevelmax),tailb(1:nboundary,1:nlevelmax), &
         & numbb(1:nboundary,1:nlevelmax) )
    headl=0;taill=0;numbl=0;numbtot=0;headb=0;tailb=0;numbb=0
    ! Allocate tree arrays
    allocate(next(1:ngridmax))
    allocate(nbor(1:ngridmax,1:twondim))
    nbor=0; next=0
    ! Allocate grid center coordinates
    allocate(xg(1:ngridmax,1:ndim))
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
    allocate(son(1:ncell),cpu_map(1:ncell))
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

  end subroutine init_amr


  subroutine clear_amr

    implicit none

    if(allocated(son)) deallocate(son,cpu_map)
    if(allocated(xg))  deallocate(xg,nbor,next)
    if(allocated(headl)) deallocate(headl,taill,numbl,numbtot,headb,tailb,numbb)
    if(allocated(var)) deallocate(var)

    return

  end subroutine clear_amr



end module module_file



 ! subroutine read_amr_hydro_mine(repository,snapnum,icpu,ramses_var,nleaf,nvar)

 !    ! purpose: use only local variables for OMP
 !    !$ use OMP_LIB
 !    implicit none 

 !    integer(kind=4),intent(in)  :: snapnum,icpu
 !    character(1000),intent(in)  :: repository

 !    character(1000)             :: nomfich 
 !    integer(kind=4),allocatable :: ind_grid(:),iig(:),grid(:)
 !    real(kind=8),allocatable    :: xxg(:)
 !    logical                     :: ok
 !    integer(kind=4)             :: i,nx,ny,nz,nlevelmax,nboundary
 !    integer(kind=4)             :: ilevel,ncache,ibound,idim,ind,iskip,iunit,iu2,rank

 !    ! stuff read from AMR files
 !    integer(kind=4),intent(out)      :: nleaf, nvar
 !    integer(kind=4)                  :: ncoarse_l,ngridmax_l
 !    real(kind=8),allocatable         :: xg_l(:,:)      ! grids position
 !    integer,allocatable              :: nbor_l(:,:)    ! neighboring father cells
 !    integer,allocatable              :: next_l(:)      ! next grid in list
 !    integer,allocatable              :: son_l(:)       ! sons grids
 !    integer,allocatable              :: cpu_map_l(:)   ! domain decomposition
 !    integer,allocatable              :: headl_l(:,:),taill_l(:,:),numbl_l(:,:),numbtot_l(:,:)
 !    integer,allocatable              :: headb_l(:,:),tailb_l(:,:),numbb_l(:,:)
 !    real(KIND=8),dimension(1:3)      :: xbound_l=(/0d0,0d0,0d0/)  

 !    real(kind=8)                :: dx
 !    integer(kind=4)             :: ix,iy,iz,istart,ivar,igrid,nvarH,nvarRT,ncell_l,ileaf,icell
 !    real(kind=8),allocatable    :: xc(:,:),xx(:)

 !    ! stuff read from the HYDRO files
 !    real(kind=8),allocatable                 :: var_l(:,:)
 !    real(kind=8),allocatable,intent(out)     :: ramses_var(:,:)

 !    rank = 1
 !    !$ rank = OMP_GET_THREAD_NUM()
 !    iunit=10+rank*2
 !    iu2 = 10+rank*2+1

 !    ! Vérification de l'existence des fichiers AMR
 !    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
 !    inquire(file=nomfich, exist=ok)
 !    if(.not. ok)then
 !       write(*,*)'File '//TRIM(nomfich)//' not found'    
 !       stop
 !    end if

 !    open(unit=iunit,file=nomfich,form='unformatted',status='old',action='read')
 !    ! Read grid variables
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)nx,ny,nz
 !    xbound_l=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

 !    ! Critical parameter: define the root level of the tree
 !    ncoarse_l=nx*ny*nz
 !    read(iunit)nlevelmax
 !    read(iunit)ngridmax_l
 !    read(iunit)nboundary
 !    read(iunit)!ngrid_current
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    read(iunit)
 !    allocate( headl_l(1:ncpu,1:nlevelmax),taill_l(1:ncpu,1:nlevelmax), &
 !         & numbl_l(1:ncpu,1:nlevelmax),numbtot_l(1:10,1:nlevelmax), &
 !         & headb_l(1:nboundary,1:nlevelmax),tailb_l(1:nboundary,1:nlevelmax), &
 !         & numbb_l(1:nboundary,1:nlevelmax) )
 !    headl_l=0;taill_l=0;numbl_l=0;numbtot_l=0;headb_l=0;tailb_l=0;numbb_l=0
 !    ! Allocate tree arrays
 !    allocate(next_l(1:ngridmax_l))
 !    allocate(nbor_l(1:ngridmax_l,1:twondim))
 !    nbor_l=0; next_l=0
 !    ! Allocate grid center coordinates
 !    allocate(xg_l(1:ngridmax_l,1:ndim))
 !    xg_l=0.0D0
 !    ! Read levels variables
 !    read(iunit)headl_l(1:ncpu,1:nlevelmax)
 !    read(iunit)taill_l(1:ncpu,1:nlevelmax)
 !    read(iunit)numbl_l(1:ncpu,1:nlevelmax)
 !    read(iunit)numbtot_l(1:10,1:nlevelmax)
 !    ! Read boundary linked list
 !    if(nboundary>0)then
 !       read(iunit)headb_l(1:nboundary,1:nlevelmax)
 !       read(iunit)tailb_l(1:nboundary,1:nlevelmax)
 !       read(iunit)numbb_l(1:nboundary,1:nlevelmax)
 !    end if
 !    !  Read free memory
 !    read(iunit)
 !    next_l(ngridmax_l) = 0
 !    ! Read cpu boundaries
 !    read(iunit)
 !    read(iunit)
 !    ncell_l=ncoarse_l+twotondim*ngridmax_l
 !    allocate(son_l(1:ncell_l),cpu_map_l(1:ncell_l))
 !    son_l=0; cpu_map_l=0
 !    ! Read coarse level
 !    read(iunit)son_l(1:ncoarse_l)       
 !    read(iunit)
 !    read(iunit)cpu_map_l(1:ncoarse_l)
 !    do ilevel=1,nlevelmax
 !       do ibound=1,nboundary+ncpu
 !          if(ibound<=ncpu)then
 !             ncache=numbl_l(ibound,ilevel)
 !          else
 !             ncache=numbb_l(ibound-ncpu,ilevel)
 !          end if
 !          if(ncache>0)then
 !             allocate(ind_grid(1:ncache))
 !             allocate(xxg(1:ncache))
 !             allocate(iig(1:ncache))
 !             allocate(grid(1:ncache))
 !             ! Read grid index
 !             read(iunit)ind_grid
 !             ! Read next index
 !             read(iunit)iig
 !             do i=1,ncache
 !                next_l(ind_grid(i))=iig(i)
 !             end do
 !             ! Read prev index (skip)
 !             read(iunit)iig
 !             ! Read grid center
 !             do idim=1,ndim
 !                read(iunit)xxg
 !                do i=1,ncache
 !                   xg_l(ind_grid(i),idim)=xxg(i)
 !                end do
 !             end do
 !             ! Read father index (skip)
 !             read(iunit)iig
 !             ! Read nbor index
 !             do ind=1,twondim
 !                read(iunit)iig
 !                do i=1,ncache
 !                   nbor_l(ind_grid(i),ind)=iig(i)
 !                end do
 !             end do
 !             ! Read son index
 !             do ind=1,twotondim
 !                iskip=ncoarse_l+(ind-1)*ngridmax_l
 !                read(iunit)iig
 !                do i=1,ncache
 !                   son_l(ind_grid(i)+iskip)=iig(i)
 !                end do
 !             end do
 !             ! Read cpu map
 !             do ind=1,twotondim
 !                iskip=ncoarse_l+(ind-1)*ngridmax_l
 !                read(iunit)iig
 !                do i=1,ncache
 !                   cpu_map_l(ind_grid(i)+iskip)=iig(i)
 !                end do
 !             end do
 !             ! Read refinement map (skip)
 !             do ind=1,twotondim
 !                read(iunit)!iig 
 !             end do
 !             deallocate(xxg,iig,grid,ind_grid)
 !          end if
 !       end do
 !    end do
 !    close(iunit)
 !    ! => can return cpu_map_l & son_l

 !    !print*,'in module_ramses.read_amr_hydro -> ',ncell_l,nvar,icpu,snapnum
 !    ! and then the hydro file
 !    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
 !    open(unit=iunit,file=nomfich,form='unformatted',status='old',action='read')
 !    read(iunit)
 !    read(iunit)nvarH
 !    read(iunit)
 !    read(iunit)nlevelmax
 !    read(iunit)nboundary
 !    read(iunit)

 !    ! Open RT file and get nvarRT
 !    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
 !    open(unit=iu2,file=nomfich,status='old',form='unformatted')
 !    read(iu2)
 !    read(iu2)nvarRT
 !    read(iu2)
 !    read(iu2)
 !    read(iu2)
 !    read(iu2)
 !    nvar  = nvarH + nvarRT

 !    allocate(ramses_var(1:ncell_l,1:nvarH+nvarRT))
 !    allocate(xc(1:twotondim,1:ndim))

 !    do ilevel=1,nlevelmax

 !       dx=0.5d0**ilevel
 !       do ind=1,twotondim
 !          iz=(ind-1)/4
 !          iy=(ind-1-4*iz)/2
 !          ix=(ind-1-2*iy-4*iz)
 !          xc(ind,1)=(dble(ix)-0.5D0)*dx
 !          xc(ind,2)=(dble(iy)-0.5D0)*dx
 !          xc(ind,3)=(dble(iz)-0.5D0)*dx
 !       end do

 !       do ibound=1,nboundary+ncpu
 !          if(ibound<=ncpu)then  ! in the box 
 !             ncache=numbl_l(ibound,ilevel)   ! nb of grids in the simulated box. 
 !             istart=headl_l(ibound,ilevel)   ! head of grid list of simulated box
 !          else                  ! boundaries of simulated volume (aka useless)
 !             ncache=numbb_l(ibound-ncpu,ilevel)
 !             istart=headb_l(ibound-ncpu,ilevel)
 !          end if
 !          read(iunit)!ilevel2
 !          read(iunit)!numbl2

 !          read(iu2)
 !          read(iu2)

 !          if(ncache>0)then
 !             allocate(ind_grid(1:ncache))
 !             allocate(xx(1:ncache))
 !             ! Loop over level grids
 !             igrid=istart
 !             do i=1,ncache
 !                ind_grid(i)=igrid
 !                igrid=next_l(igrid)
 !             end do
 !             ! Loop over cells
 !             do ind=1,twotondim
 !                iskip=ncoarse_l+(ind-1)*ngridmax_l
 !                ! Loop over conservative variables
 !                do ivar=1,nvarH
 !                   read(iunit) xx
 !                   if (ibound > ncpu) cycle  ! dont bother with boundaries
 !                   do i = 1, ncache
 !                      ramses_var(ind_grid(i)+iskip,ivar) = xx(i)
 !                   end do
 !                end do

 !                do ivar=1,nvarRT
 !                   read(iu2) xx
 !                   if (ibound > ncpu) cycle  ! dont bother with boundaries
 !                   do i = 1, ncache
 !                      ramses_var(ind_grid(i)+iskip,ivar+nvarH) = xx(i)
 !                   end do
 !                end do
 !             end do
 !             deallocate(ind_grid,xx)
 !          end if
 !       end do
 !    end do
 !    deallocate(xc)
 !    close(iunit)
 !    close(iu2)
 !    ! => can return var_l, cell_x_l, cell_y_l, cell_z_l, cell_level_l

 !    deallocate(headl_l, taill_l, numbl_l, numbtot_l, headb_l, tailb_l, numbb_l)
 !    deallocate(next_l, nbor_l, xg_l)

 !    allocate(var_l(nvar,ncell_l))

 !    ileaf = 0
 !    !print*, 'ncell : ',ncell_l
 !       do icell = 1,ncell_l
 !          if (son_l(icell)==0 .and. cpu_map_l(icell) == icpu) then
 !             ileaf = ileaf + 1
 !             do ivar = 1,nvar
 !                var_l(ivar,ileaf) = ramses_var(icell,ivar)
 !             end do
 !          end if
 !       end do
 !       deallocate(ramses_var)
 !       nleaf = ileaf
 !       !nleaf_tot = nleaf_tot + nleaf
 !       !print*, 'nleaf, nleaf_tot : ',nleaf, nleaf_tot
 !       allocate(ramses_var(nvar,nleaf))
 !       do ivar=1,nvar
 !          ramses_var(ivar,:) = var_l(ivar,1:nleaf)
 !       end do
 !       deallocate(var_l)

 !    return

 !  end subroutine read_amr_hydro_mine

