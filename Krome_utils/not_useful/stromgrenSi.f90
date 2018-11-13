program computeSiFractions
	use py_cell_utils
	use krome_main
	use krome_user
	use mpi
	implicit none

	integer			 :: ierr,rank,npsize,status0(MPI_STATUS_SIZE)
	real(kind=8)	 :: t1=0d0,t2=0d0
	character*10     :: ranktxt
	integer(kind=4)  :: chunksize
	integer,parameter:: master=0
    character(1000)  :: ramDir
    integer,parameter:: nBins=krome_nPhotoBins
    integer(kind=4)  :: ts,levelmax,i,j,k,icells
    integer(kind=4)  :: binIndices(nBins), hydroIndices(5)
    integer,parameter:: nfields = nBins + 5
    integer(kind=4)  :: usefulIndices(nfields)
    real(kind=8)     :: center(3),radius,xH,yHe,nH,nHe
    integer(kind=4)  :: ncells
    real(kind=8),allocatable  :: cells(:,:)
    real(kind=8),allocatable  :: cell_pos(:,:)
    integer(kind=4),allocatable :: cell_level(:)
    logical			 :: readRT
    integer,parameter::nsp=krome_nmols !number of species (common)
    real(kind=8),allocatable :: densities(:,:),finaldensities(:,:)
    real(kind=8)	 :: Tgas,dt,myJ(nBins),ramsesEgy(nBins),energyGaps(nBins),phbinLimits(nBins+1),oldE
    character(1000)  :: species(nsp)
	
	
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
	
	t1=MPI_WTIME()
	
	open(unit=14,status='scratch')
    rewind(14)
    write (14,*) rank
    rewind(14)
    read  (14,*) ranktxt
    close (14)
	
	
	ramDir = '/Users/mauerhof/Documents/RamsesFiles/stromgren/'
	ts = 17
	levelmax = 6
	center(:) = (/ 5d-1, 5d-1, 5d-1/)
	radius = 1d0
	readRT = .true.
	hydroIndices = (/ 1, 6, 7, 8, 9 /)   !Density, pressure, xHII, xHeII, xHeIII
	binIndices = (/ (10+i*4, i=0,nBins-1) /) !nBins photon Bins
	
	usefulIndices(1:size(hydroIndices)) = hydroIndices
	usefulIndices(size(hydroIndices)+1:size(hydroIndices)+size(binIndices)) = binIndices
	
	ramsesEgy = (/ 10.65, 14.85, 19.77, 28.35, 37.98, 48.27, 68.25 /)   !Mean energies of Ramses' simulation of n photoBins
	xH = 0.76  !Hydrogen mass fraction
	yHe = 0.24 !Helium mass fraction
	
	call count_cells(ramDir,ts,levelmax,ncells,center,radius)   !return the total number of cells    
	chunksize = (ncells+npsize-1)/npsize						!chunksize is the size of the array that each thread will treat. Except the last thread, which has potentially a smaller array
	
	if(rank==master) then										!The master allocates all the data
		allocate(cells(ncells,nfields))
		allocate(cell_pos(ncells,3))
		allocate(cell_level(ncells))
		call read_cells_hydro(ramDir,ts,levelmax,cells,cell_pos,cell_level,ncells,usefulIndices,nfields,center,radius,readRT)

		do i=1,npsize-2											!Then the master sends each piece of data to other processus
			call MPI_SEND(cells(i*chunksize+1:(i+1)*chunksize,:),chunksize*nfields,MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,ierr)
		end do
		call MPI_SEND(cells((npsize-1)*chunksize+1:ncells,:),(ncells-(npsize-1)*chunksize)*nfields,MPI_DOUBLE_PRECISION,npsize-1,11,MPI_COMM_WORLD,ierr)  !The last thread has another size
	elseif(rank==npsize-1) then
		allocate(cells(ncells-(npsize-1)*chunksize,nfields))																							  !The last thread has another size
		call MPI_RECV(cells,(ncells-(npsize-1)*chunksize)*nfields,MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
	else
		allocate(cells(chunksize,nfields))
		call MPI_RECV(cells,chunksize*nfields,MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
	end if
	
	if (rank /= npsize-1) then
		allocate(densities(chunksize,nsp))
	else
		allocate(densities(ncells-(npsize-1)*chunksize,nsp))																							  !The last thread has another size
	end if
	
	
	call krome_init() !init krome (mandatory)
	
	
	phbinLimits = (/ 8.15, 13.60, 16.35, 24.59, 33.49, 45.14, 54.42, 1e2 /)    !Ionization energies of H, He and Si,  except for 1e2, which is arbitrary,  may have to rethink that
	call krome_set_photobinE_limits(phbinLimits) 
	energyGaps(:) = krome_get_photoBinE_delta()
	
	dt = 5d12 !best I found on 23.03.18,  first simulation


	do icells=1,size(densities,1)
		
		myJ(:) = (/ (3.22254d-8*cells(icells,5+j)*ramsesEgy(j)/energyGaps(j), j=1,nBins) /)   !3.22254d-8 is kpc_to_cm/Myr_to_s*h_cgs/4Pi/ev_to_erg,  to get myJ in ev/s/cm2/sr/Hz
		call krome_set_photobinJ(myJ(:)) !initialize custom intensities

		
		nH = xH*cells(icells,1)
		nHe = 0.25*yHe*cells(icells,1)
		
		densities(icells,:) = (/ nH*cells(icells,3) + nHe*(cells(icells,4)+2*cells(icells,5)), nH*(1-cells(icells,3)), nHe*(1-cells(icells,4)-cells(icells,5)), 3.467d-6, nH*cells(icells,3), nHe*cells(icells,4), 0d0, nHe*cells(icells,5), 0d0, 0d0, 0d0 /) !(E, HI, HeI, SiI, HII, HeII, SiII, HeIII, SiIII, SiIV, SiV)

		Tgas = cells(icells,2)/cells(icells,1)*1.66d-24/1.38062d-16*(3.09d21/3.1556926d13)**2/(xH*(1+cells(icells,3)) + 0.25*yHe*(1+cells(icells,4)+2*cells(icells,5)))
		
		oldE = densities(icells,1)
		do k=1,10
			call krome(densities(icells,:), Tgas, dt)
			if (ABS(oldE-densities(icells,1))/oldE < 1d-9) exit
			oldE = densities(icells,1)
		end do
	end do
	
	
	open(unit=rank,file='outstromgrenSi'//trim(ranktxt)//'.dat',status='replace',form='formatted',action='write')
	do i=1,size(densities,1)
		write(rank,fmt=*) densities(i,4),densities(i,7),densities(i,9),densities(i,10),densities(i,11)
	end do
	
	
	t2=MPI_WTIME()
	write(*,*) rank,t2-t1
	
	close(rank)
	if (rank==master) then
		deallocate(cells,cell_pos,cell_level,densities)
	else
		deallocate(cells,densities)
	end if
	call MPI_FINALIZE(ierr)
end program computeSiFractions
