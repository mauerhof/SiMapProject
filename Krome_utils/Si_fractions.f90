PROGRAM Si_fractions

  use krome_main
  use krome_user
  use mpi
  use module_spectra
  use module_ramses
  use module_domain

  implicit none

  !Variables for mpi
  integer(kind=4)	 	:: ierr, rank, npsize, status0(MPI_STATUS_SIZE)
  integer(kind=4),allocatable   :: chunksize(:), disp(:)
  integer,parameter		:: master=0
  real(kind=8)	 		:: t1=0d0, t2=0d0

  !Variables for Krome
  integer,parameter		:: nBinK=krome_nPhotoBins, nsp=krome_nmols, nIon=4, nPhoto=krome_nPhotoRates
  real(kind=8)                  :: densities(nsp), binsK(nBinK+1)


  integer(kind=4)               :: ncells, nSEDgroups, nfields, nstars
  real(kind=8)                  :: rates_Si(nPhoto), nH_i, nHe, nSi, Tgas_i, xHII, xHeII, xHeIII
  integer(kind=4),allocatable   :: cell_level(:)
  real(kind=8),allocatable      :: nH(:), Tgas(:), mets(:), cells_rt(:,:), fractions(:,:), Si(:,:), cell_pos(:,:), cells(:,:), nHI(:)
  real(kind=8),allocatable      :: nH_cpu(:), Tgas_cpu(:), mets_cpu(:), cells_rt_cpu(:,:), fractions_cpu(:,:), Si_cpu(:,:)
  real(kind=8),allocatable      :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), Si_csn(:,:), Si_csn_help(:,:), star_L(:), star_L_tot(:)
  integer(kind=4)		:: i, j, k, icells, count, narg
  character(2000)               :: parameter_file
  character(10)                 :: type
  type(domain)                  :: compute_dom

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [Si_fractions] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: reading_method = 'fromlist'        ! If 'fromlist', takes all the data from a data_krome file (written with python, FortranFile). Otherwise read directly from simu
  character(2000)           :: ramDir = '.'                       ! Path to the Ramses simulation
  integer(kind=4)           :: timestep = 12                      ! Timestep of the simulation
  character(2000)           :: cell_data_file = '../cells_infos'  ! In case reading_method=fromlist, path to the file where the list of cells is
  
  character(2000)           :: output_path = '../outputs/out'     ! Path and name of file where the outputs will be written,  one per core

  !Miscellaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------


  ! -------------------- Command related to mpi ---------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  allocate(chunksize(npsize), disp(npsize))
  count = 0
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

  
  ! -------------------- master computes the mean-cross-sections -----------------------------------
  nSEDgroups = get_nSEDgroups(ramDir,timestep)
  allocate(Si_csn(nSEDgroups,nIon)) ; Si_csn = 0d0
  if(rank==master) then

     !get star properties in domain
     type = 'sphere'
     call domain_constructor_from_scratch(compute_dom, type, 5d-1, 5d-1, 5d-1, 1d3) 
     call ramses_read_stars_in_domain(ramDir,timestep,compute_dom,star_pos,star_age,star_mass,star_vel,star_met) ; deallocate(star_vel)

     nstars = size(star_age)

     !initialize SED properties
     call init_SED_table()

     !compute SED luminosities,  in #photons/s
     allocate(star_L(nSEDgroups), star_L_tot(nSEDgroups), Si_csn_help(nSEDgroups,nIon)) ; star_L_tot = 0d0 
     do i=1,nstars
        call inp_sed_table(star_age(i), star_met(i), 1, .false., star_L(:))    !Third variable :  1 for L[#photons/s],  3 for mean energy in bin,  2+2*Iion for mean cross-section,  Iion: 1:HI,2:HeI,3:HeII,4:SiI, etc
        do j=1,nIon
           call inp_sed_table(star_age(i), star_met(i), 8+2*j, .false., Si_csn_help(:,j))
        end do

        do j=1,nSEDgroups
           Si_csn(j,:) = Si_csn(j,:) + star_L(j)*Si_csn_help(j,:)
           star_L_tot(j) = star_L_tot(j) + star_L(j)
        end do
     end do
     do j=1,nSEDgroups
        Si_csn(j,:) = Si_csn(j,:)/star_L_tot(j)
     end do
     deallocate(star_age, star_met, star_mass, Si_csn_help)
  end if
  call MPI_BCAST(Si_csn, nSEDgroups*nIon, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)          !Every process needs Si_csn
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- read/compute the data and distribute among the cores ------------------------
  if(reading_method=='fromlist') then
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
  else
     if(rank==master) then
        call read_leaf_cells(ramDir, timestep, ncells, nfields, cell_pos, cells, cell_level)   ; deallocate(cell_pos,cell_level)
        allocate(nH(ncells), nHI(ncells), Tgas(ncells), mets(ncells), cells_rt(nSEDgroups,ncells), fractions(3,ncells))
        call ramses_get_nh_cgs(ramDir,timestep,ncells,nfields,cells,nH)
        call ramses_get_T_nhi_cgs(ramDir,timestep,ncells,nfields,cells,Tgas,nHI) ; deallocate(nHI)
        call ramses_get_metallicity(ncells,nfields,cells,mets)
        call ramses_get_fractions(ncells,nfields,cells,fractions)
        call ramses_get_flux(ramDir,timestep,ncells,nfields,nSEDgroups,cells,cells_rt)
     endif
     call MPI_BCAST(ncells, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)          !Every process needs ncells
  endif

  !Size of the array for each core.  A little bit smaller for the last core.
  chunksize = (ncells+npsize-1)/npsize
  chunksize(npsize) = ncells - (npsize-1)*chunksize(1)
  disp(:) = (/ ((i-1)*chunksize(1), i=1,npsize) /) !Necessary for scatterv

  allocate(nH_cpu(chunksize(rank+1)), Tgas_cpu(chunksize(rank+1)), mets_cpu(chunksize(rank+1)), cells_rt_cpu(nSEDgroups, chunksize(rank+1)), fractions_cpu(3, chunksize(rank+1)))


  ! Distribute the data among the processes
  call MPI_SCATTERV(nH, chunksize, disp, MPI_DOUBLE_PRECISION, nH_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(Tgas, chunksize, disp, MPI_DOUBLE_PRECISION, Tgas_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(mets, chunksize, disp, MPI_DOUBLE_PRECISION, mets_cpu, chunksize(rank+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(cells_rt, chunksize*nSEDgroups, disp*nSEDgroups, MPI_DOUBLE_PRECISION, cells_rt_cpu, chunksize(rank+1)*nSEDgroups, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(fractions, chunksize*3, disp*3, MPI_DOUBLE_PRECISION, fractions_cpu, chunksize(rank+1)*3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

  if(rank==master) deallocate(nH, Tgas, mets, cells_rt, fractions)
  ! ------------------------------------------------------------------------------------------------

  

  !init krome (mandatory)
  call krome_init()

  !Where the results are saved
  allocate(Si_cpu(nIon, chunksize(rank+1)))
  Si_cpu = 0d0

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  
  t2=MPI_WTIME()
  if(rank==master) write(*,*) 'Initilization finished,  beginning computations after time ', t2-t1
  
  ! -------------------- Loop over all the cells ---------------------------------------------------
  do icells=1,chunksize(rank+1)

     !Initialization of values in the cell
     nH_i = nH_cpu(icells)
     nHe = 7.895d-2*nH_i !Assuming mass fraction of helium of 0.24 
     nSi = mets_cpu(icells)/0.0134*3.24d-5*nH_i !Assuming solar abundances
     xHII = fractions_cpu(1,icells) ; xHeII = fractions_cpu(2,icells) ; xHeIII = fractions_cpu(3,icells)
     Tgas_i = Tgas_cpu(icells)

     !Rates for Si are the sum over each radiation group i of the product (flux_i in the cell) *  (mean cross_section_i)
     rates_Si(:) = (/ (sum(cells_rt_cpu(:,icells)*Si_csn(:,i)), i=1,nIon) /)
     !Hack so that all the SiI is photoionized in SiII, have to check the accuracy of this approximation
     rates_Si(1) = min(2*nSi, 1d-9)


     !Sets the Krome inut to have the correct photoionization rates (subroutine added by hand in krome_user.f90)
     call krome_set_photoBin_rates(rates_Si)

     
     !Computing and writing the outputs
     
     !If the gas is very hot, no need to compute,  all the silicone is in SiV or more ionized
     if(Tgas_i > 1d6) then
        Si_cpu(:,icells) = (/ 1d-18, 1d-18, 1d-18, 1d-18 /)
     else
        !Here, the gas is ~dominated by SiV   
        if(Tgas_i > 8d4) then
           densities(:) = (/ max(nH_i*xHII + nHe*(xHeII + 2*xHeIII) + 4*nSi,1d-18), max(nH_i*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nH_i*xHII,1d-18), max(nHe*xHeII,1d-18), 1d-18, max(nHe*xHeIII,1d-18), 1d-18, 1d-18, max(nSi,1d-18) /)
        !Here, the gas is ~dominated by SiIII   
        elseif(Tgas_i > 2.04d4) then
           densities(:) = (/ max(nH_i*xHII + nHe*(xHeII + 2*xHeIII) + 2*nSi,1d-18), max(nH_i*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nH_i*xHII,1d-18), max(nHe*xHeII,1d-18), 1d-18, max(nHe*xHeIII,1d-18), max(nSi,1d-18), 1d-18, 1d-18 /)
        !Here, the gas is ~dominated by SiII
        else
           densities(:) = (/ max(nH_i*xHII + nHe*(xHeII + 2*xHeIII) + nSi,1d-18), max(nH_i*(1-xHII),1d-18), max(nHe*(1-xHeII-xHeIII),1d-18), 1d-18, max(nH_i*xHII,1d-18), max(nHe*xHeII,1d-18), max(nSi,1d-18), max(nHe*xHeIII,1d-18), 1d-18, 1d-18, 1d-18 /)
        end if
        
        call krome_equilibrium(densities(:), Tgas_i)

        if(max(densities(4), densities(7), densities(9), densities(10), densities(11))/nSi > 1.0001) then  !cell bugs,  < 1/100000,   no real solution
           count = count + 1
           if(Tgas_i < 2d4) then
              Si_cpu(:,icells) = (/ 1d-18, nSi, 1d-18, 1d-18 /)
           elseif(Tgas_i < 1d5) then
              Si_cpu(:,icells) = (/ 1d-18, 1d-18, nSi, 1d-18 /)
           else
              Si_cpu(:,icells) = (/ 1d-18, 1d-18, 1d-18, 1d-18 /)
           endif
        else
           Si_cpu(:,icells) = (/ max(densities(4),1d-18), max(densities(7),1d-18), max(densities(9),1d-18), max(densities(10),1d-18) /)
        endif
       ! if(nH_i > 100) write(*,*) rank, icells, Tgas_i, nH_i, densities(4)/nSi

     end if

     if(modulo(icells,100000)==0 .and. verbose) print*,icells, 'cells computed by rank', rank
  end do
  ! ------------------------------------------------------------------------------------------------
  

  ! -------------------- Gather the outputs and write ----------------------------------------------
  if(rank==master .and. verbose) write(*,*) count, 'cells had a bug'
  if(rank==master) allocate(Si(nIon,ncells))
  call MPI_GATHERV(Si_cpu, chunksize(rank+1)*nIon, MPI_DOUBLE_PRECISION, Si*nIon, chunksize*nIon, disp*nIon, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  if(rank==master) then
     open(unit=11,file=trim(output_path),status='replace',form='unformatted',action='write')
     write(11) Si
     close(11)
  end if
  ! ------------------------------------------------------------------------------------------------
  
  
  ! -------------------- Ending of program ---------------------------------------------------------
  t1=MPI_WTIME()
  write(*,*) rank, 'computation time :', t1-t2

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
          case('reading_method')
             write(reading_method,'(a)') trim(value)
          case('cell_data_file')
             write(cell_data_file,'(a)') trim(value)
          case ('ramDir')
             write(ramDir,'(a)') trim(value)
          case('timestep')
             read(value,*) timestep
          case ('output_path')
             write(output_path,'(a)') trim(value)
          case('verbose')
             read(value,*) verbose

          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    call read_spectra_params(pfile)

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
       write(unit,'(a,a)')       '  reading_method   = ',trim(reading_method)
       if(reading_method=='fromlist') then
          write(unit,'(a,a)')       '  cell_data_file   = ',trim(cell_data_file)
       else
          write(unit,'(a,a)')       '  ramDir           = ',trim(ramDir)
          write(unit,'(a,i5)')      '  timestep         = ',timestep
       endif
       write(unit,'(a,a)')       '  output_path      = ',trim(output_path)
       write(unit,'(a,L1)')      '  verbose          = ',verbose
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')     '[Si_fractions]'
       write(*,'(a,a,a)')     '[Si_fractions]'
       write(*,'(a,a)')       '  reading_method   = ',trim(reading_method)
       if(reading_method=='fromlist') then
          write(*,'(a,a)')       '  cell_data_file   = ',trim(cell_data_file)
       else
          write(*,'(a,a)')       '  ramDir           = ',trim(ramDir)
          write(*,'(a,i5)')      '  timestep         = ',timestep
       endif
       write(*,'(a,a)')       '  output_path      = ',trim(output_path)
       write(*,'(a,L1)')      '  verbose          = ',verbose
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_Si_fractions_params


END PROGRAM Si_fractions
