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
  character(2000)           :: reading_method = 'fromlist'        ! If 'fromlist', takes all the data from a data_krome file (written with python, FortranFile). Otherwise read directly from simu
  character(2000)           :: repository = '.'                       ! Path to the Ramses simulation
  integer(kind=4)           :: snapnum = 12                      ! Timestep of the simulation
  character(2000)           :: cell_data_file = '../cells_infos'  ! In case reading_method=fromlist, path to the file where the list of cells is

  character(2000)           :: output_path = '../outputs/out'     ! Path and name of file where the outputs will be written,  one per core

  !Miscellaneous
  logical                   :: verbose = .true.
  integer(kind=4)           :: nbundle = 10
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
  !if (verbose .and. rank==0) call print_ion_fractions_params
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
     call master(repository, snapnum, reading_method, cell_data_file, nbundle, output_path)
  else
     ! Worker section, will mostly do the Krome calls
     call worker(nbundle)
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
          case('reading_method')
             write(reading_method,'(a)') trim(value)
          case('cell_data_file')
             write(cell_data_file,'(a)') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case('snapnum')
             read(value,*) snapnum
          case ('output_path')
             write(output_path,'(a)') trim(value)
          case('verbose')
             read(value,*) verbose
          case('nbundle')
             read(value,*) nbundle

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
       write(unit,'(a,a)')       '  reading_method   = ',trim(reading_method)
       if(reading_method=='fromlist') then
          write(unit,'(a,a)')       '  cell_data_file   = ',trim(cell_data_file)
       else
          write(unit,'(a,a)')       '  repository           = ',trim(repository)
          write(unit,'(a,i5)')      '  snapnum         = ',snapnum
       endif
       write(unit,'(a,a)')       '  output_path      = ',trim(output_path)
       write(unit,'(a,L1)')      '  verbose          = ',verbose
       write(unit,'(a,i5)')      '  nbundle          = ',nbundle
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')     '[ion_fractions]'
       write(*,'(a,a)')       '  reading_method   = ',trim(reading_method)
       if(reading_method=='fromlist') then
          write(*,'(a,a)')       '  cell_data_file   = ',trim(cell_data_file)
       else
          write(*,'(a,a)')       '  repository           = ',trim(repository)
          write(*,'(a,i5)')      '  snapnum         = ',snapnum
       endif
       write(*,'(a,a)')       '  output_path      = ',trim(output_path)
       write(*,'(a,L1)')      '  verbose          = ',verbose
       write(*,'(a,i5)')      '  nbundle          = ',nbundle
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_ion_fractions_params

end program main

