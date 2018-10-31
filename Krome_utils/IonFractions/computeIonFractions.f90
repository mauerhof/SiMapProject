PROGRAM computeSiFractions
  use py_cell_utils
  use krome_main
  use krome_user
  use mpi
  use SED_module,only:getSEDcsn


  implicit none

  !Variables for mpi
  integer(kind=4)	 			:: ierr, rank, npsize, status0(MPI_STATUS_SIZE), chunksize
  integer,parameter			:: master=0
  real(kind=8)	 			:: t1=0d0, t2=0d0
  character*10     			:: ranktxt


  !Variables for Krome
  !Many variables of size nBins-3 because here Krome doesn't ionize HI, HeI and HeII so I have
  integer,parameter			:: nBinK=krome_nPhotoBins, nsp=krome_nmols, nIon=nBinK !(nIon is not strictly useful but helps understand)    !nBinK is the number of Bins in Krome, which needs not be greater than the number of metal ions that you ionize (for example SiI, SiII, SiIII, SiIV  -> nBinK=4). 
  real(kind=8) 				:: binsK(nBinK+1), energyGaps(nBinK), middleEnergies(nBinK), kromeSigma(nBinK), kromeJ(nBinK), dt, oldE, nH, nHe, Tgas, densities(nsp), kromeCoeffs(nIon,nBinK)


  !Variable for Ramses
  integer,parameter			:: nBinR = 7, nfields = 5+nBinR
  !nBinR is the number of bins used in Ramses, which could be in principle way bigger than nBinK
  integer(kind=4)				:: ts, hydroIndices(5), binIndices(nBinR), usefulIndices(nfields), levelmax, ncells
  real(kind=8),allocatable	:: ages(:), Zs(:), Ls(:), cells(:,:), cell_pos(:,:), SEDs(:,:,:)           								 ! SEDs f(lambda,age,met)
  integer(kind=4),allocatable	:: cell_level(:)
  real(kind=8)				:: center(3), radius, unit_length, unit_time, ratesRamses(nIon), binsR(nBinR+1), sigmaN(nIon,nBinR)   !First index : species,  second index : bin
  character(len=1000)			:: sed_folder, ramDir
  logical						:: readRT


  !General variables
  real(kind=8)				:: xH, yHe, ev_to_hz, PI, testRates(7)
  integer(kind=4)				:: i , j, k, icells, l, orderIonE(nIon)!, orderBins(nBins-3)



  !Commands related to mpi
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  t1=MPI_WTIME()										!Initialize time to display the total time of execution
  open(unit=14,status='scratch')
  rewind(14)
  write (14,*) rank
  rewind(14)											!5 lines to create a character of the rank of the core,  to put it in names of output files   (more straightforward way ?)
  read  (14,*) ranktxt
  close (14)


!!!Can obviously improve this part to import numbers from ramses output,  would be more elegant
  ev_to_hz = 2.41799d14 ; PI=3.14159266
  xH = 0.76 ; yHe = 0.24 ; binsR = (/ 8.15d0, 1.36d1, 1.635d1, 2.459d1, 3.349d1, 4.514d1, 5.442d1, 1d3 /)  !Maybe should use infinity instead of 1d3,  changes almost nothing
  unit_length = 3.09d21 ; unit_time = 3.1556926d13


  sed_folder = '/Users/mauerhof/Documents/seds/bpass100/'
  call readSeds(sed_folder,ages,Zs,Ls,SEDs)


  orderIonE = (/ 4, 5, 6, 7 /)  !Indices of (Si) ions in rt_parameters

  do j=1,nIon
     sigmaN(j,:) = (/ (getSEDcsn(Ls, SEDs(:,10,2), size(Ls), binsR(i), binsR(i+1), orderIonE(j)), i=1,nBinR) /)    !!!!SiI and SiII show a rise of sigmaN, weird.   sigmaN is the number weighted average of the cross section over each bin, for each species of Si.   (SEDs(:,10,2) is the SED for bpass100 with z=0.002 and age=6.3 Myr (to check)

!!!THIS PART IS THE ONLY ONE THAT REALLY NEED TO BE CHANGED, in order to use the same sed as Ramses, which is an interpolation over all stellar particles. The easiest way is to ask ramses to print Csn for Si  (and also C, O, etc).
  end do


  ramDir = '/Users/mauerhof/Documents/RamsesFiles/idealized/G8/lmax11Ngroups5/'
  ts = 12; levelmax = 11 ; center(:) = (/ 5d-1, 5d-1, 5d-1/) ; radius = 3d-2 ; readRT = .true.
  hydroIndices = (/ 1, 6, 7, 8, 9 /)   !Density, pressure, xHII, xHeII, xHeIII
  binIndices = (/ (10+i*4, i=0,nBinR-1) /) !nBinR photon Bins
  usefulIndices(1:size(hydroIndices)) = hydroIndices
  usefulIndices(size(hydroIndices)+1:size(hydroIndices)+size(binIndices)) = binIndices


  call count_cells(ramDir,ts,levelmax,ncells,center,radius)   !return the total number of cells 

  chunksize = (ncells+npsize-1)/npsize							!Size of the array for each core.  A little bit smaller for the last core.
  if (rank==npsize-1) chunksize = ncells - (npsize-1)*chunksize	



  if(rank==master) then		!The master reads all the data
     allocate(cells(ncells,nfields))
     allocate(cell_pos(ncells,3))
     allocate(cell_level(ncells))
     call read_cells_hydro(ramDir,ts,levelmax,cells,cell_pos,cell_level,ncells,usefulIndices,nfields,center,radius,readRT)


     do i=1,npsize-2			!Then the master sends each piece of data to the other threads
        call MPI_SEND(cells(i*chunksize+1:(i+1)*chunksize,:),chunksize*nfields,MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,ierr)
     end do
     call MPI_SEND(cells((npsize-1)*chunksize+1:ncells,:),(ncells-(npsize-1)*chunksize)*nfields,MPI_DOUBLE_PRECISION,npsize-1,11,MPI_COMM_WORLD,ierr)  !The last thread data has another size
  else
     allocate(cells(chunksize,nfields))		!The other cores receive their piece of data
     call MPI_RECV(cells,chunksize*nfields,MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
  end if



  call krome_init() !init krome (mandatory)
  binsK = (/ 8.15d0, 1.635d1, 3.349d1, 4.514d1, 6d1 /)  !TO AUTOMATIZE !!!
  call krome_set_photobinE_limits(binsK)
  energyGaps(:) = krome_get_photoBinE_delta() ; middleEnergies(:) = krome_get_photoBinE_mid()


  do i=1,nIon
     kromeSigma(:) = krome_get_xsec(orderIonE(i))  !Photoionization cross section used in Krome, for each Si ion. Useful here to impose a given photoionzation rate (given by Ramses output).
     kromeCoeffs(i,:) = (/ (4*PI*energyGaps(j)*ev_to_hz*kromeSigma(j)/middleEnergies(j), j=1,nBinK) /)  !Coefficients of the matrix of a linear system to be solved to find the right input to give to Krome in order to have the wanted photoionization rate
  end do


  open(unit=rank,file='./G8l11/out'//trim(ranktxt)//'.dat',status='replace',form='formatted',action='write')  !To write the outputs
  !open(unit=rank+20,file='rates'//trim(ranktxt)//'.dat',status='replace',form='formatted',action='write')   !To check if the rates are reproduced correctly

  dt = 5d12 !best I found on 23.03.18,  may be inappropriate,  still have to think about the strategy to converge to equilibrium


  !Loop over all the cells
  do icells=1,chunksize
     ratesRamses(:) = (/ (sum(cells(icells,size(hydroIndices)+1:size(hydroIndices)+nBinR)*sigmaN(i,:))*unit_length/unit_time, i=1,nIon) /)  !Photoionization rate of each Si ion, computed from ramses output. This is the inhomogeneous part of the linear system cited above

     call solveTriangularMatrix(nIon,ratesRamses(:),kromeCoeffs(:,:),kromeJ(:))   !Solves the linear system to find to correct Krome inputs
     !kromeJ(:) = (/ 0d0, 0d0, 0d0, 0d0 /)
     call krome_set_photobinJ(kromeJ(:)) !Sets the Krome inut to have the correct photoionization rates (i.e ramsesRates)

     nH = xH*cells(icells,1)		
     nHe = 0.25*yHe*cells(icells,1)		

     densities(:) = (/ nH*cells(icells,3) + nHe*(cells(icells,4)+2*cells(icells,5)), nH*(1-cells(icells,3)), nHe*(1-cells(icells,4)-cells(icells,5)), 3.467d-6, nH*cells(icells,3), nHe*cells(icells,4), 1d-20, nHe*cells(icells,5), 1d-20, 1d-20, 1d-20 /) !(E, HI, HeI, SiI, HII, HeII, SiII, HeIII, SiIII, SiIV, SiV), as needed for Krome

     Tgas = cells(icells,2)/cells(icells,1)*1.66d-24/1.38062d-16*(unit_length/unit_time)**2/(xH*(1+cells(icells,3)) + 0.25*yHe*(1+cells(icells,4)+2*cells(icells,5)))  !Trust me

     !Loop to find the equilibrium state of Krome,  still to be improved
     oldE = densities(1)
     do k=1,10
        call krome(densities(:), Tgas, dt)				
        if (ABS(oldE-densities(1))/oldE < 1d-9) exit    !To be improved,   cells that are outside the stromgren sphere never reach this condition...
        oldE = densities(1)
     end do

     !		testRates = krome_get_photoBin_rates()  !To check if the rates are reproduced correctly
     !		write(rank+20,fmt=*) ratesRamses
     !		write(rank+20,fmt=*) testRates(4:7)
     !		write(rank+20,fmt=*) ''

     write(rank,fmt=*) densities(4),densities(7),densities(9),densities(10),densities(11)   !Writes output for one cell
  end do


  close(rank) !; close(rank+20)
  deallocate(ages, Zs, Ls, SEDs)
  if (rank==master) then
     deallocate(cells,cell_pos,cell_level)
  else
     deallocate(cells)
  end if

  t2=MPI_WTIME()
  write(*,*) rank,t2-t1

  call MPI_FINALIZE(ierr)




CONTAINS


  SUBROUTINE readSeds(path,ages,Zs,Ls,SEDs)

    implicit none

    character(len=128),intent(in)::path
    integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
    real(kind=8),allocatable,intent(out)::ages(:), Zs(:), Ls(:)
    real(kind=8),allocatable,intent(out)::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
    !real(kind=8),allocatable::tbl(:,:,:), tbl2(:,:,:), reb_tbl(:,:,:)
    integer::i,ia,iz,ip,ii,dum
    character(len=128)::fZs, fAges, fSEDs                        ! Filenames
    logical::ok,okAge,okZ


    !if(path.eq.'')call get_environment_variable('RAMSES_SED_DIR',path)
    inquire(FILE=TRIM(path)//'/all_seds.dat', exist=ok)
    if(.not. ok)then
       !if(myid.eq.1) then 
       write(*,*)'Cannot access SED directory ',TRIM(path)
       write(*,*)'Directory '//TRIM(path)//' not found'
       write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
            ' to the correct path, or use the namelist.'
       !endif
       !call clean_stop
    end if
    write(fZs,'(a,a)')   trim(path),"/metallicity_bins.dat"
    write(fAges,'(a,a)') trim(path),"/age_bins.dat"
    write(fSEDs,'(a,a)') trim(path),"/all_seds.dat"
    inquire(file=fZs, exist=okZ)
    inquire(file=fAges, exist=okAge)
    inquire(file=fSEDs, exist=ok)
    if(.not. ok .or. .not. okAge .or. .not. okZ) then
       !if(myid.eq.1) then 
       write(*,*) 'Cannot read SED files...'
       write(*,*) 'Check if SED-directory contains the files ',  &
            'metallicity_bins.dat, age_bins.dat, all_seds.dat'
       !endif
       !call clean_stop
    end if

    ! READ METALLICITY BINS-------------------------------------------------
    open(unit=10,file=fZs,status='old',form='formatted')
    read(10,'(i8)') nzs
    allocate(zs(nzs))
    do i = 1, nzs
       read(10,'(e14.6)') zs(i)
    end do
    close(10)
    ! READ AGE BINS---------------------------------------------------------
    open(unit=10,file=fAges,status='old',form='formatted')
    read(10,'(i8)') nAges
    allocate(ages(nAges))
    do i = 1, nAges
       read(10,'(e14.6)') ages(i)
    end do
    close(10)
    ages = ages*1.e-9                       !         Convert from yr to Gyr
    if(ages(1) .ne. 0.) ages(1) = 0.
    ! READ SEDS-------------------------------------------------------------
    open(unit=10,file=fSEDs,status='old',form='unformatted')
    read(10) nLs, dum
    allocate(Ls(nLs))
    read(10) Ls(:)
    allocate(SEDs(nLs,nAges,nzs))
    do iz = 1, nzs
       do ia = 1, nAges
          read(10) SEDs(:,ia,iz)
       end do
    end do
    close(10)

    return

    deallocate(ages, Zs, Ls, SEDs)


  END SUBROUTINE readSeds

END PROGRAM computeSiFractions




SUBROUTINE solveTriangularMatrix(n, a, B, x)  !Solve the linear system Bx = a,  with B a triangular matrix, with the lower left half equal zero.

  implicit none

  integer(kind=4),intent(in)	:: n
  real(kind=8),intent(in)		:: a(n), B(n,n)
  real(kind=8),intent(out)	:: x(n)
  integer(kind=4)				:: i,j

  do i=n,1,-1
     x(i) = a(i)/B(i,i)
     do j=n,i+1,-1
        x(i) = x(i) - B(i,j)*x(j)/B(i,i)
     end do
  end do

  return

END SUBROUTINE solveTriangularMatrix



