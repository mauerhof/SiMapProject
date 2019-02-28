module module_parallel_mpi_mine

  use mpi

  implicit none

  ! MPI constants
  integer(kind=4)                                 :: rank          ! cpu number
  integer(kind=4)                                 :: nb_cpus       ! total number of cpus
  integer(kind=4)                                 :: code, error
  integer(kind=4)                                 :: tag, done_tag, ierror, exi_tag=2
  integer(kind=4)                                 :: nworker
  integer(kind=4),dimension(MPI_STATUS_SIZE)      :: status

  ! define a MPI derived type for photons
  integer(kind=4),parameter                       :: nbloc=10
  integer(kind=4),dimension(nbloc)                :: types, longueurs_blocs
  integer(kind=MPI_ADDRESS_KIND),dimension(nbloc) :: deplacements, adresses


contains


  subroutine start_mpi

    ! Initialize the MPI execution environment
    call mpi_init(code)

    ! Determines the rank of the cpu
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

    ! Determines the number of cpus
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_cpus, code)
    
  end subroutine start_mpi



  subroutine finish_mpi

    ! Terminates MPI execution environment
    call mpi_finalize(code)
    
  end subroutine finish_mpi



  subroutine stop_mpi

    call MPI_ABORT(MPI_COMM_WORLD,error,code)

  end subroutine stop_mpi


end module module_parallel_mpi_mine

