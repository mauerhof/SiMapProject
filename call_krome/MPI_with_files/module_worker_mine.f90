module module_worker_mine

  use module_parallel_mpi_mine
  use module_cell
  use module_file

  private 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [worker] of the parameter file
  ! --------------------------------------------------------------------------
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  public :: worker, read_worker_params, print_worker_params
  
contains

  subroutine worker(repository, snapnum, output_path)
    
    implicit none
    
    character(2000),intent(in)     :: repository, output_path
    integer(kind=4),intent(in)     :: snapnum
    integer(kind=4)                :: juseless, i, icpu
    real(kind=8)                   :: start_file, end_file

    allocate(csn(3,2))
    call MPI_RECV(csn, 3*2, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

    do while (status(MPI_TAG) .ne. EXI_TAG)

       ! receive my domain number
       call MPI_RECV(juseless, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       if(status(MPI_TAG) == EXI_TAG) then
          write(*,'(a,i5.5,a)') ' [w',rank,'] exit tagged received'
          exit
       end if
       
       ! receive my list of cells to compute
       call MPI_RECV(icpu, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       if(verbose) write(*,'(a,i5.5,a,i5.5)') ' [w',rank,'] begins computation of file ', icpu
       call cpu_time(start_file)
       call compute_file(repository, snapnum, icpu, output_path)
       call cpu_time(end_file)

       if(verbose)then
          write(*,'(a,i5.5,a,i5.5,a,f12.6,a)') ' [w',rank,'] time to compute file ', icpu, ' = ',end_file-start_file,' seconds.'
       endif

       ! send my results
       call MPI_SEND(juseless, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, code)
       call MPI_SEND(icpu, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, code)

    end do

    if(verbose) write(*,'(a,i4.4,a)') ' [w',rank,'] : exit of loop...'

    ! final synchronization, for profiling purposes
    call MPI_BARRIER(MPI_COMM_WORLD,code)

  end subroutine worker

  

  subroutine read_worker_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! NB: does not call read_params of depdencies (master module does that). 
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
       if (line(1:8) == '[worker]') then
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
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_worker_params


  
  subroutine print_worker_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[worker]'
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')             '[worker]'
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
    end if

    return

  end subroutine print_worker_params

end module module_worker_mine
