PROGRAM Si_fractions

  use krome_main
  use krome_user
  use mpi

  implicit none

  !Variables for mpi
  integer(kind=4)	 	:: ierr, rank, npsize, status0(MPI_STATUS_SIZE), chunksize
  integer,parameter		:: master=0
  real(kind=8)	 		:: t1=0d0, t2=0d0
  character*10     		:: ranktxt

  !Variables for Krome
  integer,parameter		:: nBinK=krome_nPhotoBins, nsp=krome_nmols, nIon=4, nPhoto=krome_nPhotoRates
  real(kind=8)                  :: densities(nsp), binsK(nBinK+1)


  integer(kind=4)               :: ncells, nSEDgroups
  real(kind=8)                  :: rates_Si(nPhoto), nHi, nHe, nSi, Tgas_i, xHII, xHeII, xHeIII
  real(kind=8),allocatable      :: nH(:), Tgas(:), mets(:), cells_rt(:,:), fractions(:,:), Si_csn(:,:)
  integer(kind=4)		:: i, j, k, icells, count, narg
  character(2000)               :: parameter_file

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [CreateDomDump] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: cell_data_file = '../cells_infos'       ! Path to the file where the list of cells is
  character(2000)           :: Si_csn_file = '../Si_csn'          ! Path to the file with the mean silicone cross-sections
  character(2000)           :: output_path = '../outputs/out'     !Path and name of file where the outputs will be written,  one per core
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------


  ! -------------------- Commands related to mpi ---------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  !Initialize time to display the total time of execution
  t1=MPI_WTIME()				
  open(unit=14,status='scratch')
  rewind(14)
  write (14,*) rank
  rewind(14)		!5 lines to create a character of the rank of the core,  to put it in names of output files   (more straightforward way ?)
  read  (14,*) ranktxt
  close (14)
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

  !Size of the array for each core.  A little bit smaller for the last core if ncells is not divisible by npsize
  chunksize = (ncells+npsize-1)/npsize			
  if (rank==npsize-1) chunksize = ncells - (npsize-1)*chunksize	

  !The master reads all the data
  if(rank==master) then	
     allocate(nH(ncells), Tgas(ncells), mets(ncells), cells_rt(ncells,nSEDgroups), fractions(ncells,3))
     read(10) nH
     read(10) Tgas
     read(10) mets
     read(10) cells_rt
     read(10) fractions

     !Then the master sends each piece of data to the other threads
     do i=1,npsize-2	
        call MPI_SEND(nH(i*chunksize+1:(i+1)*chunksize),chunksize,MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,ierr)
        call MPI_SEND(Tgas(i*chunksize+1:(i+1)*chunksize),chunksize,MPI_DOUBLE_PRECISION,i,12,MPI_COMM_WORLD,ierr)
        call MPI_SEND(mets(i*chunksize+1:(i+1)*chunksize),chunksize,MPI_DOUBLE_PRECISION,i,13,MPI_COMM_WORLD,ierr)
        call MPI_SEND(cells_rt(i*chunksize+1:(i+1)*chunksize,:),chunksize*nSEDgroups,MPI_DOUBLE_PRECISION,i,14,MPI_COMM_WORLD,ierr)
        call MPI_SEND(fractions(i*chunksize+1:(i+1)*chunksize,:),chunksize*3,MPI_DOUBLE_PRECISION,i,15,MPI_COMM_WORLD,ierr)
     end do
     !The last core's vector have another size
     call MPI_SEND(nH((npsize-1)*chunksize+1:ncells),ncells-(npsize-1)*chunksize,MPI_DOUBLE_PRECISION,npsize-1,11,MPI_COMM_WORLD,ierr) 
     call MPI_SEND(Tgas((npsize-1)*chunksize+1:ncells),ncells-(npsize-1)*chunksize,MPI_DOUBLE_PRECISION,npsize-1,12,MPI_COMM_WORLD,ierr)
     call MPI_SEND(mets((npsize-1)*chunksize+1:ncells),ncells-(npsize-1)*chunksize,MPI_DOUBLE_PRECISION,npsize-1,13,MPI_COMM_WORLD,ierr)
     call MPI_SEND(cells_rt((npsize-1)*chunksize+1:ncells,:),(ncells-(npsize-1)*chunksize)*nSEDgroups,MPI_DOUBLE_PRECISION,npsize-1,14,MPI_COMM_WORLD,ierr)
     call MPI_SEND(fractions((npsize-1)*chunksize+1:ncells,:),(ncells-(npsize-1)*chunksize)*3,MPI_DOUBLE_PRECISION,npsize-1,15,MPI_COMM_WORLD,ierr)
     !The other cores receive their piece of data
  else
     allocate(nH(chunksize), Tgas(chunksize), mets(chunksize), cells_rt(chunksize,nSEDgroups), fractions(chunksize,3))
     call MPI_RECV(nH,chunksize,MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(Tgas,chunksize,MPI_DOUBLE_PRECISION,master,12,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(mets,chunksize,MPI_DOUBLE_PRECISION,master,13,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(cells_rt,chunksize*nSEDgroups,MPI_DOUBLE_PRECISION,master,14,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(fractions,chunksize*3,MPI_DOUBLE_PRECISION,master,15,MPI_COMM_WORLD,status0,ierr)
  end if
  close(10)
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- Reads the silicone mean cross-sections ------------------------------------
  allocate(Si_csn(nSEDgroups,nIon))
  open(unit=10, file=trim(Si_csn_file), form='unformatted', action='read')
  read(10) Si_csn
  close(10)
  ! ------------------------------------------------------------------------------------------------


  !init krome (mandatory)
  call krome_init()

  !To write the outputs
  open(unit=rank,file=trim(output_path)//trim(ranktxt),status='replace',form='unformatted',action='write') 


  ! -------------------- Loop over all the cells ---------------------------------------------------
  do icells=1,chunksize
     
     !Rates for Si are the sum over each radiation group i of the product (flux_i in the cell) *  (mean cross_section_i)
     rates_Si(4:7) = (/ (sum(cells_rt(icells,:)*Si_csn(:,i)), i=1,nIon) /)
     !Zero rate for hydrogen and helium, we don't want to change their density
     rates_Si(1:3) = 0d0
     !Hack so that all the SiI is photoionized in SiII, have to check the accuracy of this approximation
     rates_Si(4) = 1d-5


     !Sets the Krome inut to have the correct photoionization rates (subroutine added by hand in krome_user.f90)
     call krome_set_photoBin_rates(rates_Si)

     !Initialization of values in the cell
     nHi = nH(icells)
     nHe = 7.895d-2*nHi !Assuming mass fraction of helium of 0.24 
     nSi = mets(icells)/0.0134*3.24d-5*nHi !Assuming solar abundances
     xHII = fractions(icells,1) ; xHeII = fractions(icells,2) ; xHeIII = fractions(icells,3)
     Tgas_i = Tgas(icells)


     !Computing and writing the outputs

     !If the gas is very hot, no need to compute,  all the silicone is in SiV or more ionized
     if(Tgas_i > 1d6) then
        write(rank) 1d-18, 1d-18, 1d-18, 1d0  !!SiII, SiIII, SiIV, SiV

     !Here, the gas is ~dominated by SiV   
     elseif(Tgas_i > 8d4) then
        densities(:) = (/ max(nHi*xHII + nHe*(xHeII + 2*xHeIII) + 4*nSi,1d-18), max(nHi*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nHi*xHII,1d-18), max(nHe*xHeII,1d-18), 1d-18, max(nHe*xHeIII,1d-18), 1d-18, 1d-18, max(nSi,1d-18) /)
        call krome_equilibrium(densities(:), Tgas_i)
        write(rank) max(densities(7), 1d-18), max(densities(9), 1d-18), max(densities(10), 1d-18), max(densities(11), 1d-18)   !SiII, SiIII, SiIV, SiV

     !Here, the gas is ~dominated by SiIII   
     elseif(Tgas_i > 2.04d4) then
        densities(:) = (/ max(nHi*xHII + nHe*(xHeII + 2*xHeIII) + 2*nSi,1d-18), max(nHi*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nHi*xHII,1d-18), max(nHe*xHeII,1d-18), 1d-18, max(nHe*xHeIII,1d-18), max(nSi,1d-18), 1d-18, 1d-18 /)
        call krome_equilibrium(densities(:), Tgas_i)
        write(rank) max(densities(7), 1d-18), max(densities(9), 1d-18), max(densities(10), 1d-18), max(densities(11), 1d-18)   !SiII, SiIII, SiIV, SiV

     !Here, the gas is ~dominated by SiII
     else
        densities(:) = (/ max(nHi*xHII + nHe*(xHeII + 2*xHeIII) + nSi,1d-18), max(nHi*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nHi*xHII,1d-18), max(nHe*xHeII,1d-18), max(nSi,1d-18), max(nHe*xHeIII,1d-18), 1d-18, 1d-18, 1d-18 /)
        call krome_equilibrium(densities(:), Tgas_i)
        write(rank) max(densities(7), 1d-18), max(densities(9), 1d-18), max(densities(10), 1d-18), max(densities(11), 1d-18)   !SiII, SiIII, SiIV, SiV
     end if
     
     
     if(modulo(icells,10000)==0) print*,icells, 'cells computed by rank', rank
  end do
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- Ending of program ---------------------------------------------------------
  close(rank)
  deallocate(nH, Tgas, mets, cells_rt, fractions, Si_csn)

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
