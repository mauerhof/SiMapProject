PROGRAM Si_fractions

  use krome_main
  use krome_user
  use mpi

  implicit none

  !Variables for mpi
  integer(kind=4)	 	:: ierr, rank, npsize, status0(MPI_STATUS_SIZE)
  integer(kind=4),allocatable   :: chunksize(:), disp(:)
  integer,parameter		:: master=0
  real(kind=8)	 		:: t1=0d0, t2=0d0

  !Variables for Krome
  integer,parameter		:: nBinK=krome_nPhotoBins, nsp=krome_nmols, nIon=4, nPhoto=krome_nPhotoRates
  real(kind=8)                  :: densities(nsp), binsK(nBinK+1)


  integer(kind=4)               :: ncells, nSEDgroups
  real(kind=8)                  :: rates_Si(nPhoto), nHi, nHe, nSi, Tgas_i, xHII, xHeII, xHeIII
  real(kind=8),allocatable      :: nH(:), Tgas(:), mets(:), cells_rt(:,:), fractions(:,:), SiII(:), Si_csn(:,:)
  real(kind=8),allocatable      :: nH_cpu(:), Tgas_cpu(:), mets_cpu(:), cells_rt_cpu(:,:), fractions_cpu(:,:), SiII_cpu(:)
  integer(kind=4)		:: i, j, k, icells, count, narg
  character(2000)               :: parameter_file

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [Si_fractions] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: cell_data_file = '../cells_infos'       ! Path to the file where the list of cells is
  character(2000)           :: Si_csn_file = '../Si_csn'          ! Path to the file with the mean silicone cross-sections
  character(2000)           :: output_path = '../outputs/out'     !Path and name of file where the outputs will be written,  one per core

  !Miscellaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------


  ! -------------------- Commands related to mpi ---------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  allocate(chunksize(npsize), disp(npsize))
  !Initialize time to display the total time of execution
  t1=MPI_WTIME()				
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: Si_fractions params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_Si_fractions_params(parameter_file)
  if (verbose .and. rank==0) call print_Si_fractions_params
  ! ------------------------------------------------------------


  ! -------------------- Open the data files and distribute among the cores ------------------------
  open(unit=10, file=trim(cell_data_file), form='unformatted', action='read')
  read(10) ncells
  read(10) nSEDgroups

  !The master reads all the data
  if(rank==master) then
     allocate(nH(ncells), Tgas(ncells), mets(ncells), cells_rt(nSEDgroups,ncells), fractions(3,ncells))
     read(10) nH
     read(10) Tgas
     read(10) mets
     read(10) cells_rt
     read(10) fractions
  end if
  close(10)

  !Size of the array for each core.  A little bit smaller for the last core.
  chunksize = (ncells+npsize-1)/npsize
  chunksize(npsize) = ncells - (npsize-1)*chunksize(1)
  disp(:) = (/ ((i-1)*chunksize(1), i=1,npsize) /) !Necessary for scatterv
  
  allocate(nH_cpu(chunksize(rank+1)), Tgas_cpu(chunksize(rank+1)), mets_cpu(chunksize(rank+1)), cells_rt_cpu(nSEDgroups, chunksize(rank+1)), fractions_cpu(3, chunksize(rank+1)))


  ! Disrtibute the data among the processes
  call MPI_SCATTERV(nH, chunksize, disp, MPI_DOUBLE_PRECISION, nH_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(Tgas, chunksize, disp, MPI_DOUBLE_PRECISION, Tgas_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(mets, chunksize, disp, MPI_DOUBLE_PRECISION, mets_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(cells_rt, chunksize*nSEDgroups, disp*nSEDgroups, MPI_DOUBLE_PRECISION, cells_rt_cpu, chunksize(rank+1)*nSEDgroups, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(fractions, chunksize*3, disp*3, MPI_DOUBLE_PRECISION, fractions_cpu, chunksize(rank+1)*3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

  if(rank==master) deallocate(nH, Tgas, mets, cells_rt, fractions)
  ! ------------------------------------------------------------------------------------------------

  

  ! -------------------- Reads the silicone mean cross-sections ------------------------------------
  allocate(Si_csn(nSEDgroups,nIon))
  open(unit=10, file=trim(Si_csn_file), form='unformatted', action='read')
  read(10) Si_csn
  close(10)
  ! ------------------------------------------------------------------------------------------------


  !init krome (mandatory)
  call krome_init()

  !Where the results are saved
  allocate(SiII_cpu(chunksize(rank+1)))
  SiII_cpu = 0d0
  

  ! -------------------- Loop over all the cells ---------------------------------------------------
  do icells=1,chunksize(rank+1)
     
     !Rates for Si are the sum over each radiation group i of the product (flux_i in the cell) *  (mean cross_section_i)
     rates_Si(4:7) = (/ (sum(cells_rt_cpu(:,icells)*Si_csn(:,i)), i=1,nIon) /)
     !Zero rate for hydrogen and helium, we don't want to change their density
     rates_Si(1:3) = 0d0
     !Hack so that all the SiI is photoionized in SiII, have to check the accuracy of this approximation
     rates_Si(4) = 1d-5


     !Sets the Krome inut to have the correct photoionization rates (subroutine added by hand in krome_user.f90)
     call krome_set_photoBin_rates(rates_Si)

     !Initialization of values in the cell
     nHi = nH_cpu(icells)
     nHe = 7.895d-2*nHi !Assuming mass fraction of helium of 0.24 
     nSi = mets_cpu(icells)/0.0134*3.24d-5*nHi !Assuming solar abundances
     xHII = fractions_cpu(1,icells) ; xHeII = fractions_cpu(2,icells) ; xHeIII = fractions_cpu(3,icells)
     Tgas_i = Tgas_cpu(icells)

     
     !Computing and writing the outputs
     
     !If the gas is very hot, no need to compute,  all the silicone is in SiV or more ionized
     if(Tgas_i > 1d6) then
        SiII_cpu(icells) = 1d-18
     else

        !Here, the gas is ~dominated by SiV   
        if(Tgas_i > 8d4) then
           densities(:) = (/ max(nHi*xHII + nHe*(xHeII + 2*xHeIII) + 4*nSi,1d-18), max(nHi*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nHi*xHII,1d-18), max(nHe*xHeII,1d-18), 1d-18, max(nHe*xHeIII,1d-18), 1d-18, 1d-18, max(nSi,1d-18) /)

        !Here, the gas is ~dominated by SiIII   
        elseif(Tgas_i > 2.04d4) then
           densities(:) = (/ max(nHi*xHII + nHe*(xHeII + 2*xHeIII) + 2*nSi,1d-18), max(nHi*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nHi*xHII,1d-18), max(nHe*xHeII,1d-18), 1d-18, max(nHe*xHeIII,1d-18), max(nSi,1d-18), 1d-18, 1d-18 /)

        !Here, the gas is ~dominated by SiII
        else
           densities(:) = (/ max(nHi*xHII + nHe*(xHeII + 2*xHeIII) + nSi,1d-18), max(nHi*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nHi*xHII,1d-18), max(nHe*xHeII,1d-18), max(nSi,1d-18), max(nHe*xHeIII,1d-18), 1d-18, 1d-18, 1d-18 /)
        end if
        
        call krome_equilibrium(densities(:), Tgas_i)
        SiII_cpu(icells) = densities(7)
     end if

     if(modulo(icells,50000)==0 .and. verbose) print*,icells, 'cells computed by rank', rank
  end do
  ! ------------------------------------------------------------------------------------------------
  

  ! -------------------- Gather the outputs and write ----------------------------------------------
  if(rank==master) allocate(SiII(ncells))
  call MPI_GATHERV(SiII_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, SiII, chunksize, disp, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  if(rank==master) then
     open(unit=11,file=trim(output_path),status='replace',form='unformatted',action='write')
     write(11) SiII
     close(11)
  end if
  ! ------------------------------------------------------------------------------------------------
  
  
  ! -------------------- Ending of program ---------------------------------------------------------
  t2=MPI_WTIME()
  write(*,*) rank,t2-t1

  call MPI_FINALIZE(ierr)
  ! ------------------------------------------------------------------------------------------------



contains

  subroutine read_Si_fractions_params(pfile)

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
       if (line(1:14) == '[Si_fractions]') then
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
          case('cell_data_file')
             write(cell_data_file,'(a)') trim(value)
          case ('Si_csn_file')
             write(Si_csn_file,'(a)') trim(value)
          case ('output_path')
             write(output_path,'(a)') trim(value)
          case('verbose')
             read(value,*) verbose

          end select
       end do
    end if
    close(10)

    return

  end subroutine read_Si_fractions_params


  subroutine print_Si_fractions_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt

    if (present(unit)) then 
       write(unit,'(a,a,a)')     '[Si_fractions]'
       write(unit,'(a,a)')       '  cell_data_file   = ',trim(cell_data_file)
       write(unit,'(a,a)')       '  Si_csn_file      = ',trim(Si_csn_file)
       write(unit,'(a,a)')       '  output_path      = ',trim(output_path)
       write(unit,'(a,L1)')      '  verbose          = ',verbose
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')     '[Si_fractions]'
       write(*,'(a,a)')       '  cell_data_file   = ',trim(cell_data_file)
       write(*,'(a,a)')       '  Si_csn_file      = ',trim(Si_csn_file)
       write(*,'(a,a)')       '  output_path      = ',trim(output_path)
       write(*,'(a,L1)')      '  verbose          = ',verbose
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_Si_fractions_params


END PROGRAM Si_fractions
