module module_master_mine

  use module_parallel_mpi_mine
  use module_spectra
  use module_ramses
  use module_cell
  use module_krome

  implicit none

  private
  
  !integer(kind=4),dimension(:),allocatable         :: first,last,nqueue
  integer(kind=4),dimension(:,:),allocatable       :: next
  integer(kind=4),dimension(:), allocatable        :: cpu
  real(kind=8), dimension(:), allocatable          :: delta
  type(cell),dimension(:),allocatable              :: cellgrid, cellpacket
  integer(kind=4)                                  :: last

  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [master] of the parameter file
  ! --------------------------------------------------------------------------
  logical                   :: verbose = .false.
  ! checkpoint/restart
  logical                   :: restart = .false.    ! if true, start the run from backup file IonBakFile
  !character(2000)           :: IonBakFile = 'backup.dat'
  real(kind=8)              :: dt_backup = 2.    ! time in seconds between 2 backups, default is 7200
  ! --------------------------------------------------------------------------

  public :: master, read_master_params, print_master_params

contains

  subroutine master(repository, snapnum, reading_method, cell_data_file, nbundle, fileout)

    implicit none

    real(kind=8)                            :: start_init, end_init, time_now, dt_since_last_backup, time_last_backup
    character(2000),intent(in)              :: repository, reading_method, cell_data_file, fileout
    integer(kind=4),intent(in)              :: snapnum, nbundle
    integer(kind=4)                         :: ncells, ncellstodo, ncellsdone, i, j, k, icpu, idcpu, ncpuended, ntest
    real(kind=8)                            :: percentDone, percentBefore
    logical                                 :: everything_not_done

    
    call cpu_time(start_init)

    ! read cells (or restore from backup)
    if(restart) then
       if (verbose) print *,'[master] restoring cells from file: ',trim(fileout)
       call restore_cells(fileout,cellgrid)
       ncells = size(cellgrid)
       do i=1,ncells
          if(cellgrid(i)%den_ions(1) == 0d0) then
             last = i-1
             exit
          end if
       end do
       ncellstodo = ncells - last
       if (verbose)then
          print *,'[master] Ncells =',ncells
          print *,'[master] Ncells to do =',ncellstodo
       endif
    else
       if (verbose) print *,'[master] reading cells'
       call init_cells(repository, snapnum, reading_method, cell_data_file, cellgrid, fileout)
       ncells = size(cellgrid)
       ncellstodo = ncells
       last=0
       if (verbose) print *,'[master] N cells =',ncells
    end if
    percentBefore = real(ncells-ncellstodo)/ncells*100.


    allocate(cellpacket(nbundle))

    call cpu_time(end_init)

    if (verbose) print '(" [master] time to initialize cells in master = ",f12.3," seconds. (omp)")',end_init-start_init
    if (verbose) print*,'[master] send a first bundle of cells to each worker'
    time_last_backup = end_init

    ! send a first bundle of cells to each worker
    do icpu=1,nworker

       ! construct a list of photon packets to send
       call fill_bundle(cellpacket,nbundle)

       j=1
       ! Send something to the worker (for the exit tag)
       call MPI_SEND(j, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)

       do i = 1,nbundle 
          call MPI_SEND(cellpacket(i)%id, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(cellpacket(i)%nHI, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(cellpacket(i)%nHII, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(cellpacket(i)%nHeII, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code) 
          call MPI_SEND(cellpacket(i)%nHeIII, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code) 
          call MPI_SEND(cellpacket(i)%T, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(cellpacket(i)%Z, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(cellpacket(i)%rates, nPhotoRea, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
       end do
    end do

    everything_not_done=.true.

    ncpuended=0

    ! Receive and append to pertinent domain list
    do while(everything_not_done)

       ! First, receive information from a given CPU and identify the CPU
       call MPI_RECV(ntest, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       idcpu = status(MPI_SOURCE)

       do i = 1,nbundle
          call MPI_RECV(cellpacket(i)%id, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status, IERROR)
          call MPI_RECV(cellpacket(i)%den_ions, nmols-natoms-3, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG, MPI_COMM_WORLD, status, IERROR)
       end do

       do i=1,nbundle
          if(cellpacket(i)%id<1) exit  
          call update_grid(i)
       end do

       ncellstodo = ncells - last

       ! check if it is time to back up
       call cpu_time(time_now)
       dt_since_last_backup = time_now - time_last_backup
       if(dt_since_last_backup > dt_backup)then
          if(verbose) print*, '[master] beginning backup'
          call save_cells(fileout,cellgrid,.true.,last)
          time_last_backup = time_now
          call cpu_time(dt_since_last_backup)
          if (verbose) print *,'[master] backup done'
          if (verbose) print '(" [master] time to do the backup = ",f12.3," seconds.")',dt_since_last_backup-time_now
       endif

       if(ncellstodo <= 0)then
          ! no more cells to send

          ! first count ended cpu, to not skip last working cpu...
          ncpuended=ncpuended+1
          if(ncpuended==nworker)then
             everything_not_done=.false.
          endif
          if(verbose) print '(" [master] no more cells to send to worker ",i5," then send exit code.")',idcpu
          call MPI_SEND(idcpu, 1, MPI_INTEGER, idcpu, exi_tag , MPI_COMM_WORLD, code)

       else
          ! keep sending photons

          j=1
          call MPI_SEND(j, 1, MPI_INTEGER, idcpu, tag, MPI_COMM_WORLD, code)
          
          ! Construct a new bundle of photon packets
          call fill_bundle(cellpacket,nbundle)

          ! send it

          do i = 1,nbundle 
             call MPI_SEND(cellpacket(i)%id, 1, MPI_INTEGER, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(cellpacket(i)%nHI, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(cellpacket(i)%nHII, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(cellpacket(i)%nHeII, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code) 
             call MPI_SEND(cellpacket(i)%nHeIII, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code) 
             call MPI_SEND(cellpacket(i)%T, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(cellpacket(i)%Z, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(cellpacket(i)%rates, nPhotoRea, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
          end do

       end if

       ! print progress
       ncellsdone = count(mask=(cellgrid(:)%den_ions(1)>0))
       percentdone = real(ncellsdone)/ncells*100.

       if(percentdone>=percentBefore+1.)then
          print '(" [master] number of cell packets done = ",i8," (",f4.1," %)")',ncellsdone,percentDone
          percentBefore = percentDone + 1.
       endif

    end do

    ! final synchronization, for profiling purposes
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    call save_cells(fileout,cellgrid,.true.,ncells)
    call delete_bak(fileout)

    deallocate(cellgrid,cellpacket)

    if(verbose) print *,'[master] end'
  end subroutine master
  !===================================================================================================
  
  
  subroutine delete_bak(fileout)

    implicit none

    character(2000),intent(in)     :: fileout
    character(2000)                :: filebak, command
    integer(kind=4)                :: i,j

    do i=1,n_elements
       do j=1,n_ions(i)
          filebak = trim(fileout)//'.'//trim(element_names(elements(i)))//trim(roman_num(j))//'.bak'
          command = 'rm '//trim(filebak)
          call system(command)
       end do
    end do

  end subroutine delete_bak

  ! subroutine backup_run(fileout)

  !   character(2000),intent(in) :: fileout
  !   character(1000)            :: filebak, command
  !   logical                    :: file_exists

    

  !   ! first, copy last backup file into filebak
  !   filebak = trim(IonBakFile)//'.bak'
  !   ! check if file exists

  !   command = 'mv '//trim(IonBakFile)//' '//trim(filebak) 
  !   call system(command)

  !   ! then, write a new backup file
  !   if(verbose) print*, '[master] beginning backup'
  !   call save_cells(IonBakFile,cellgrid,.true.)

  !   command = 'rm '//trim(filebak)
  !   call system(command)

  !   if (verbose) print *,'[master] backup done'

  ! end subroutine backup_run



  subroutine update_grid(i)
    ! update cell grid with cellpacket
    ! NB: it would be clearer to pass cellpacket in argument

    implicit none
    integer(kind=4)             :: myID
    integer(kind=4), intent(in) :: i

    myID = cellpacket(i)%ID
    if (cellgrid(myID)%ID /= myID) then
       print*,'ERROR: id mismatch updating the grid'
       call stop_mpi
       return
    endif
    cellgrid(myID)%den_ions = cellpacket(i)%den_ions

  end subroutine update_grid

  
  subroutine fill_bundle(cellpacket,nbundle)
    ! fill bundle of cells packets cellpacket(nbundle)

    implicit none
    integer(kind=4), intent(in)                           :: nbundle
    type(cell), dimension(nbundle), intent(out)           :: cellpacket
    integer(kind=4)                                       :: i,fsave

    i=1
    do
       if(last+i > size(cellgrid)) then
          cellpacket(i)%ID = 0
       else
          cellpacket(i) = cellgrid(last+i)
       end if
       i = i+1
       if(i>nbundle) exit
    end do
    last = last + nbundle

  end subroutine fill_bundle


  subroutine read_master_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present
    
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:8) == '[master]') then
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
          case ('verbose')
             read(value,*) verbose
          case ('restart')
             read(value,*) restart
          !case ('IonBakFile')
             !write(IonBakFile,'(a)') trim(value)
          case ('dt_backup')
             read(value,*) dt_backup
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_master_params


  
  subroutine print_master_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[master]'
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a,L1)')          '  restart        = ',restart
       !write(unit,'(a,a)')           '  IonBakFile  = ',trim(IonBakFile)
       write(unit,'(a,f12.3)')       '  dt_backup      = ',dt_backup
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')             '[master]'
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a,L1)')          '  restart        = ',restart
      ! write(*,'(a,a)')           '  IonBakFile  = ',trim(IonBakFile)
       write(*,'(a,f12.3)')       '  dt_backup      = ',dt_backup
       write(*,'(a)')             ' '       
    end if

    return

  end subroutine print_master_params
  

end module module_master_mine

