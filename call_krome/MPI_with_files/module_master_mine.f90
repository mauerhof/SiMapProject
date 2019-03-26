module module_master_mine

  use module_parallel_mpi_mine
  use module_spectra
  use module_ramses
  use module_cell
  use module_krome
  use module_file

  implicit none

  private

  type(file),dimension(:),allocatable              :: filegrid
  integer(kind=4)                                  :: last, ilast, cpu_computed


  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [master] of the parameter file
  ! --------------------------------------------------------------------------
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  public :: master, read_master_params, print_master_params

contains

  subroutine master(repository, snapnum)

    implicit none

    real(kind=8)                            :: start_init, end_init
    character(2000),intent(in)              :: repository
    integer(kind=4),intent(in)              :: snapnum
    integer(kind=4)                         :: nfiles, nfilestodo, nfilesdone, i, j, k, icpu, idcpu, ncpuended, ntest, nSEDgroups
    real(kind=8)                            :: percentDone, percentBefore
    logical                                 :: everything_not_done


    call cpu_time(start_init)

    !Read the number of files and initialize filegrid
    call init_files(repository, snapnum, filegrid)
    nfiles = size(filegrid)
    nfilestodo = nfiles
    if (verbose) print*,'[master] Number of cpu used in the simulation : ', nfiles

    !Compute the csn in the box
    nSEDgroups = get_nOptBins()
    call init_csn(repository, snapnum)
    do icpu=1,nworker
       call MPI_SEND(csn, nSEDgroups*sum(n_ions), MPI_DOUBLE_PRECISION, icpu, tag, MPI_COMM_WORLD, code)
    end do
    call cpu_time(end_init)

    last = 1

    if (verbose) print '(" [master] time to initialize = ",f12.3," seconds.")',end_init-start_init

    ! send a first bundle of cells to each worker
    do icpu=1,nworker
       j=1
       ! Send something to the worker (for the exit tag)
       call MPI_SEND(j, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)

       !Now send the number of the file to compute
       call MPI_SEND(last, 1, MPI_INTEGER, icpu, tag, MPI_COMM_WORLD, code)
       last = last + 1
    end do

    everything_not_done=.true.

    ncpuended=0

    ! Receive and append to pertinent domain list
    do while(everything_not_done)

       ! First, receive information from a given CPU and identify the CPU
       call MPI_RECV(ntest, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)
       
       idcpu = status(MPI_SOURCE)
       
       call MPI_RECV(cpu_computed, 1, MPI_INTEGER, idcpu, DONE_TAG, MPI_COMM_WORLD, status, IERROR)
       
       filegrid(cpu_computed)%status = 1

       nfilestodo = nfiles - last + 1


       if(nfilestodo <= 0)then
          ! no more cells to send

          ! first count ended cpu, to not skip last working cpu...
          ncpuended=ncpuended+1
          if(ncpuended==nworker)then
             everything_not_done=.false.
          endif
          if(verbose) print '(" [master] no more files to send to worker ",i5," then send exit code.")',idcpu
          call MPI_SEND(idcpu, 1, MPI_INTEGER, idcpu, exi_tag , MPI_COMM_WORLD, code)

       else
          ! keep sending files
          j=1
          call MPI_SEND(j, 1, MPI_INTEGER, idcpu, tag, MPI_COMM_WORLD, code)

          !Now send the number of the file to compute
          call MPI_SEND(last, 1, MPI_INTEGER, idcpu, tag, MPI_COMM_WORLD, code)
          last = last + 1

       end if

       ! print progress
       nfilesdone = count(mask=(filegrid(:)%status)>0)
       percentdone = real(nfilesdone)/nfiles*100.

       if(percentdone>=percentBefore+1.)then
          print '(" [master] number of files done = ",i8," (",f4.1," %)")',nfilesdone,percentDone
          percentBefore = percentDone + 1.
       endif

    end do

    ! final synchronization, for profiling purposes
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    deallocate(filegrid)

    if(verbose) print *,'[master] end'
  end subroutine master
  !===================================================================================================


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
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')             '[master]'
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
    end if

    return

  end subroutine print_master_params


end module module_master_mine

