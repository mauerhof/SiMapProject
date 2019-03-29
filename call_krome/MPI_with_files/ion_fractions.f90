program main

  use module_parallel_mpi_mine
  use module_master_mine
  use module_worker_mine
  use module_ramses
  use module_spectra
  use module_krome

  implicit none

  real(kind=8)                             :: start,finish
  character(2000)                          :: parameter_file, line, file_compute_dom
  integer(kind=4)                          :: narg, i, j, ndomain

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ion_fractions] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: repository = '.'                 ! Path to the Ramses simulation
  character(2000)           :: repo_restart = '.'               ! Path to the restart
  integer(kind=4)           :: snapnum = 12                     ! Timestep of the simulation
  integer(kind=4)           :: snap_restart = 12                ! Timestep of the restart

  character(2000)           :: output_path = '.'   ! Path to the output files. The files will be written in /output_path/output_snapnum/,  so make sure that this directory exists
  !Miscellaneous
  logical                   :: verbose = .true.
  character(2000)           :: csn_file = 'csn.dat'
  ! --------------------------------------------------------------------------


  call cpu_time(start)

  call start_mpi

  nworker=nb_cpus-1
  if(nworker==0)then
     print*,'You have to run this code with MPI'
     stop
  end if

  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: ion_fractions params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  if(rank==0) print*, ' '
  call read_ion_fractions_params(parameter_file)
  if (verbose .and. rank==0) call print_ion_fractions_params
  ! ------------------------------------------------------------

  if (rank == 0 .and. verbose) then
     print*,'--> Working with Nworker =',nworker
     print*,'--> Starting master/workers pattern'
     print*,' '
  end if


  call MPI_BARRIER(MPI_COMM_WORLD,code)

  ! Master - Worker separation
  if (rank == 0) then
     ! Master section, will dispatch the jobs.
     call master(repository, snapnum, csn_file)
  else
     ! Worker section, will mostly do the Krome calls
     call worker(repository, snapnum, repo_restart, snap_restart, output_path)
  end if


  call finish_mpi
  call cpu_time(finish)
  if(verbose .and. rank==0)then
     print*,' '
     print*,'--> work done, MPI finalized'
     print '(" --> Time = ",f12.3," seconds.")',finish-start
     print*,' '
  endif



contains

  subroutine read_ion_fractions_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    logical         :: ndomain_present 

    section_present = .false.
    ndomain_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:15) == '[ion_fractions]') then
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
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('repo_restart')
             write(repo_restart,'(a)') trim(value)
          case('snapnum')
             read(value,*) snapnum
          case('snap_restart')
             read(value,*) snap_restart
          case ('output_path')
             write(output_path,'(a)') trim(value)
          case ('csn_file')
             write(csn_file,'(a)') trim(value)
          case('verbose')
             read(value,*) verbose

          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    call read_spectra_params(pfile)
    call read_master_params(pfile)
    call read_worker_params(pfile)
    call read_krome_params(pfile)

    return

  end subroutine read_ion_fractions_params


  subroutine print_ion_fractions_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt

    if (present(unit)) then 
       write(unit,'(a,a,a)')     '[ion_fractions]'
       write(unit,'(a,a)')       '  repository           = ',trim(repository)
       write(unit,'(a,i5)')      '  snapnum              = ',snapnum
       write(unit,'(a,a)')       '  repo_restart         = ',trim(repo_restart)
       write(unit,'(a,i5)')      '  snap_restart         = ',snap_restart
       write(unit,'(a,a)')       '  output_path          = ',trim(output_path)
       write(unit,'(a,L1)')      '  verbose              = ',verbose
       call print_ramses_params(unit)
       call print_spectra_params(unit)
       call print_master_params(unit)
       call print_worker_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a,a,a)')     '[ion_fractions]'
       write(*,'(a,a)')       '  repository           = ',trim(repository)
       write(*,'(a,i5)')      '  snapnum              = ',snapnum
       write(*,'(a,a)')       '  repo_restart         = ',trim(repo_restart)
       write(*,'(a,i5)')      '  snap_restart         = ',snap_restart
       write(*,'(a,a)')       '  output_path          = ',trim(output_path)
       write(*,'(a,L1)')      '  verbose              = ',verbose
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       call print_ramses_params
       call print_spectra_params
       call print_master_params
       call print_worker_params
    end if

    return

  end subroutine print_ion_fractions_params

end program main

