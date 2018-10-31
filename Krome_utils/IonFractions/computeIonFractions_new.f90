PROGRAM computeSiFractions

  use module_ramses
  use krome_main
  use krome_user
  use mpi


  implicit none

  !Variables for mpi
  integer(kind=4)	 		:: ierr, rank, npsize, status0(MPI_STATUS_SIZE)
  integer(kind=4),allocatable           :: chunksize(:)
  integer,parameter			:: master=0
  real(kind=8)	 			:: t1=0d0, t2=0d0
  character*10     			:: ranktxt


  !Variables for Krome
  integer,parameter			:: nBinK=krome_nPhotoBins, nsp=krome_nmols, nIon=4, nPhoto=krome_nPhotoRates
  real(kind=8) 				:: binsK(nBinK+1), nH, nHe, Tgas, densities(nsp), kromeCoeffs(nIon,nBinK)


  !Variable for Ramses
  integer(kind=4)			:: nSEDgroups, nfields, ncells
  real(kind=8),allocatable       	:: cells(:,:), cell_pos(:,:), nH(:), nHI(:), Tgas(:), mets(:), xHeII(:), xHeIII(:)       							
  integer(kind=4),allocatable    	:: cell_level(:)
  real(kind=8)				:: center(3), radius, rates_Si(nPhoto), sigmaN(nIon,nBinR)  

  !General variables
  integer(kind=4)			:: i, j, k, icells, l, orderIonE(nIon)!, orderBins(nBins-3)
  real(kind=8)                          :: nHI_i, nHII, Tgas_i, nHeI, nHeII, nHeIII, nSi
  character(2000)                       :: parameter_file

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [Si_fractions] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: ramDir = '.'
  integer(kind=4)           :: timestep = 12
  character(2000)           :: Si_csn_file = '../Si_csn'          ! Path to the file with the mean silicone cross-sections
  character(2000)           :: output_path = '../outputs/out'     !Path and name of file where the outputs will be written,  one per core

  !Miscellaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------


  ! -------------------- Commands related to mpi ---------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  allocate(chunksize(npsize))
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



  ! ! -------------------- Reads the silicone mean cross-sections ------------------------------------
  ! allocate(Si_csn(nSEDgroups,nIon))
  ! open(unit=10, file=trim(Si_csn_file), form='unformatted', action='read')
  ! read(10) Si_csn
  ! close(10)
  ! ! ------------------------------------------------------------------------------------------------


  ! -------------------- Read cells with module_ramses ---------------------------------------------
  !The master reads all the data
  if(rank==master) then
     call read_leaf_cells(ramDir, timestep, ncells, nfields, cell_pos, cells, cell_level)

     call ramses_get_nh_cgs(ramDir,timestep,ncells,nfields,cells,nH)
     call ramses_get_T_nhi_cgs(ramDir,timestep,ncells,nfields,cells,Tgas,nHI)
     call ramses_get_metallicity(ncells,nfields,cells,mets)
     call ramses_get_he_fractions(ncells,nfields,cells,xHeII,xHeIII)


     !The master sends ncells to the others
     do i=1,npsize-1		
        call MPI_SEND(ncells,1,MPI_INTEGER,i,11,MPI_COMM_WORLD,ierr)
     end do
  else
     call MPI_RECV(ncells,1,MPI_INTEGER,master,11,MPI_COMM_WORLD,status0,ierr)
  end if

  !Size of the array for each core.  A little bit smaller for the last core.
  chunksize = (ncells+npsize-1)/npsize							
  chunksize(npsize) = ncells - (npsize-1)*chunksize(1)	

  !Then the master sends each piece of data to the other threads
  if(rank==master) then		
     do i=1,npsize-1		
        call MPI_SEND(nH(i*chunksize(i+1)+1:(i+1)*chunksize(i+1)),chunksize(i+1),MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,ierr)
        call MPI_SEND(nHI(i*chunksize(i+1)+1:(i+1)*chunksize(i+1)),chunksize(i+1),MPI_DOUBLE_PRECISION,i,12,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xHeII(i*chunksize(i+1)+1:(i+1)*chunksize(i+1)),chunksize(i+1),MPI_DOUBLE_PRECISION,i,13,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xHeIII(i*chunksize(i+1)+1:(i+1)*chunksize(i+1)),chunksize(i+1),MPI_DOUBLE_PRECISION,i,14,MPI_COMM_WORLD,ierr)
        call MPI_SEND(Tgas(i*chunksize(i+1)+1:(i+1)*chunksize(i+1)),chunksize(i+1),MPI_DOUBLE_PRECISION,i,15,MPI_COMM_WORLD,ierr)
        call MPI_SEND(mets(i*chunksize(i+1)+1:(i+1)*chunksize(i+1)),chunksize(i+1),MPI_DOUBLE_PRECISION,i,16,MPI_COMM_WORLD,ierr)
     end do
  else
     !The other cores receive their piece of data
     allocate(nH(chunksize(rank+1)), nHI(chunksize(rank+1)), xHeII(chunksize(rank+1)), xHeIII(chunksize(rank+1)), Tgas(chunksize(rank+1)), mets(chunksize(rank+1)))
     call MPI_RECV(nH,chunksize(rank+1),MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(nHI,chunksize(rank+1),MPI_DOUBLE_PRECISION,master,12,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(xHeII,chunksize(rank+1),MPI_DOUBLE_PRECISION,master,13,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(xHeIII,chunksize(rank+1),MPI_DOUBLE_PRECISION,master,14,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(Tgas,chunksize(rank+1),MPI_DOUBLE_PRECISION,master,15,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(mets,chunksize(rank+1),MPI_DOUBLE_PRECISION,master,16,MPI_COMM_WORLD,status0,ierr)
  end if
  ! ------------------------------------------------------------------------------------------------


  !init krome (mandatory)
  call krome_init() 


  !To write the outputs
  open(unit=rank,file=trim(output_path)//trim(ranktxt),status='replace',form='unformatted',action='write') 



  ! -------------------- Loop over all the cells ---------------------------------------------------
  do icells=1,chunksize(rank+1)

     !Photoionization rate of each Si ion, computed from ramses output
     !rates_Si(4:7) = (/ (sum(cells(icells,size(hydroIndices)+1:size(hydroIndices)+nBinR)*Si_csn(:,i))*unit_length/unit_time, i=1,nIon) /)
     !Zero rate for hydrogen and helium, we don't want to change their density
     !rates_Si(1:3) = 0d0
     !Hack so that all the SiI is photoionized in SiII, have to check the accuracy of this approximation
     !rates_Si(4) = 1d-5
     rates_Si = 0d0

     !Sets the Krome inut to have the correct photoionization rates (subroutine added by hand in krome_user.f90)
     call krome_set_photoBin_rates(rates_Si)


     nHI_i = nHI(icells)
     nHII = nH(icells) - nHI_i
     !Assuming mass fraction of helium of 0.24
     nHeI = (1-xHeII(icells)-xHeIII(icells))*7.895d-2*nHi
     nHeII = xHeII(icells)*7.895d-2*nHi
     nHeIII = xHeIII(icells)*7.895d-2*nHi
     nSi = mets(icells)/0.0134*3.24d-5*nHi !Assuming solar abundances
     Tgas_i = Tgas(icells)



     !Computing and writing the outputs

     !If the gas is very hot, no need to compute,  all the silicone is in SiV or more ionized
     if(Tgas_i > 1d6)  then
        write(rank) 1d-18, 1d-18, 1d-18, 1d0  !!SiII, SiIII, SiIV, SiV

     else
        !Here, the gas is ~dominated by SiV   
        if(Tgas_i > 8d4) then
           densities(:) = (/ max(nHII + nHeII+2*nHeIII + 4*nSi,1d-18), max(nHI_,1d-18), max(nHeI, 1d-18), max(nHII,1d-18), max(nHeII,1d-18), 1d-18, max(nHeIII,1d-18), 1d-18, 1d-18, max(nSi,1d-18) /)

           !Here, the gas is ~dominated by SiIII   
        elseif(Tgas_i > 2.04d4) then
           densities(:) = (/ max(nHII + nHeII+2*nHeIII + 2*nSi,1d-18), max(nHI_,1d-18), max(nHeI,1d-18), 1d-18, max(nHII,1d-18), max(nHeII,1d-18), 1d-18, max(nHeIII,1d-18), max(nSi,1d-18), 1d-18, 1d-18 /)

           !Here, the gas is ~dominated by SiII
        else
           densities(:) = (/ max(nHII + nHeII+2*nHeIII + nSi,1d-18), max(nHI_,1d-18), max(nHeI,1d-18), 1d-18, max(nHII,1d-18), max(nHeII,1d-18), max(nSi,1d-18), max(nHeIII,1d-18), 1d-18, 1d-18, 1d-18 /)
        end if

        call krome_equilibrium(densities(:), Tgas_i)
        write(rank) max(densities(7), 1d-18), max(densities(9), 1d-18), max(densities(10), 1d-18), max(densities(11), 1d-18)   !SiII, SiIII, SiIV, SiV
     end if

     if(modulo(icells,10000)==0) print*,icells, 'cells computed by rank', rank
  end do
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- Ending of program ---------------------------------------------------------
  close(rank)

  if (rank==master) then
     deallocate(cells,cell_pos,cell_level, nH, nHI, Tgas, xHeII, xHeIII, mets)
  else
     deallocate(nH, nHI, Tgas, xHeII, xHeIII, mets)
  end if

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
          case('ramDir')
             write(ramDir,'(a)') trim(value)
          case ('Si_csn_file')
             write(Si_csn_file,'(a)') trim(value)
          case ('output_path')
             write(output_path,'(a)') trim(value)
          case('verbose')
             read(value,*) verbose
          case('timestep')
             read(value,*) timestep

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
       write(unit,'(a,a)')       '  ramDir           = ',trim(ramDir)
       write(unit,'(a,i2)')      '  timestep         = ',timestep
       write(unit,'(a,a)')       '  Si_csn_file      = ',trim(Si_csn_file)
       write(unit,'(a,a)')       '  output_path      = ',trim(output_path)
       write(unit,'(a,L1)')      '  verbose          = ',verbose
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')     '[Si_fractions]'
       write(*,'(a,a)')       '  ramDir           = ',trim(ramDir)
       write(*,'(a,i2)')      '  timestep         = ',timestep
       write(*,'(a,a)')       '  Si_csn_file      = ',trim(Si_csn_file)
       write(*,'(a,a)')       '  output_path      = ',trim(output_path)
       write(*,'(a,L1)')      '  verbose          = ',verbose
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_Si_fractions_params


END PROGRAM computeSiFractions
