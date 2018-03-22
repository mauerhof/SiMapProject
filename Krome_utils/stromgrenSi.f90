program ramses
	use py_cell_utils
	use krome_main
	use krome_user
	use mpi
	implicit none

	integer			 :: ierr,rank,npsize,status0(MPI_STATUS_SIZE)
	character*10     :: ranktxt
	integer(kind=4)  :: chunksize
	integer,parameter:: master=0, msgtag1=11, msgtag2=12
    character(1000)  :: ramDir
    integer,parameter:: nBins=krome_nPhotoBins
    integer(kind=4)  :: ts,levelmax,i,j,k,icells
    integer(kind=4)  :: binIndices(nBins), hydroIndices(5)
    integer,parameter:: nfields = nBins + 5
    integer(kind=4)  :: usefulIndices(nfields)
    real(kind=8)     :: center(3),radius,xH,yHe,nH,nHe,t1=0d0,t2=0d0
    integer(kind=4)  :: ncells
    real(kind=8),allocatable  :: cells(:,:)
    real(kind=8),allocatable  :: cell_pos(:,:)
    integer(kind=4),allocatable :: cell_level(:)
    logical			 :: readRT
    integer,parameter::nsp=krome_nmols !number of species (common)
    real(kind=8),allocatable :: densities(:,:),finaldensities(:,:)
    real(kind=8)	 :: Tgas,dt,myJ(nBins),ramsesEgy(nBins),energyGaps(nBins),phbinLimits(nBins+1),oldSi,oldSiII
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
	
	
	ramDir = '/Users/mauerhof/Documents/Ramses/stromgren/'
	ts = 17
	levelmax = 6
	center(:) = (/ 5d-1, 5d-1, 5d-1/)
	radius = 1d0
	readRT = .true.
	hydroIndices = (/ 1, 6, 7, 8, 9 /)   !Density, pressure, xHII, xHeII, xHeIII
	binIndices = (/ (10+i*4, i=0,nBins-1) /) !nBins photon Bins
	
	usefulIndices(1:size(hydroIndices)) = hydroIndices
	usefulIndices(size(hydroIndices)+1:size(hydroIndices)+size(binIndices)) = binIndices
	
	ramsesEgy = (/ 10.65, 14.85, 19.77, 28.35, 37.98, 48.27, 68.25 /)
	xH = 0.76  !Hydrogen mass fraction
	yHe = 1.-xH !Helium mass fraction
	call count_cells(ramDir,ts,levelmax,ncells,center,radius)
	!ncells=1000
	
	chunksize = (ncells+npsize-1)/npsize
	allocate(densities(MIN(ncells,(rank+1)*chunksize)-rank*chunksize,nsp))
	allocate(cells(ncells,nfields))
	allocate(cell_pos(ncells,3))
	allocate(cell_level(ncells))
	
	!if(rank==master) then
	call read_cells_hydro(ramDir,ts,levelmax,cells,cell_pos,cell_level,ncells,usefulIndices,nfields,center,radius,readRT)
	!end if

	!call MPI_BCAST(cells, ncells*nfields, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
	!call MPI_BCAST(densities, ncells*nsp, MPI_REAL, master, MPI_COMM_WORLD, ierr)
	
	
	call krome_init() !init krome (mandatory)
	
	
	phbinLimits = (/ 8.15, 13.60, 16.35, 24.59, 33.49, 45.14, 54.42, 1e2 /)
	call krome_set_photobinE_limits(phbinLimits) 
	energyGaps(:) = krome_get_photoBinE_delta()
	
	dt = 1d12 !may be inapropriate
	

	
	do i = 1, size(densities,1)
		icells = i+rank*chunksize
		myJ(:) = (/ (3.22254d-8*cells(icells,5+j)*ramsesEgy(j)/energyGaps(j), j=1,nBins) /)   !3.22254d-8 is kpc_to_cm/Myr_to_s*h_cgs/4Pi/ev_to_erg,  to get myJ in ev/s/cm2/sr/Hz
		call krome_set_photobinJ(myJ(:)) !initialize custom intensities
		
		nH = xH*cells(icells,1)
		nHe = 0.25*yHe*cells(icells,1)
	
		densities(i,:) = (/ nH*cells(icells,3) + nHe*(cells(icells,4)+2*cells(icells,5)), nH*(1-cells(icells,3)), nHe*(1-cells(icells,4)-cells(icells,5)), 2.512d-6, nH*cells(icells,3), nHe*cells(icells,4), 0d0, nHe*cells(icells,5), 0d0, 0d0, 0d0, 0d0, 0d0 /) !(E, HI, HeI, SiI, HII, HeII, SiII, HeIII, SiIII, SiIV, SiV, SiVI, SiVII)
		
		Tgas = cells(icells,2)/cells(icells,1)*1.66d-24/1.38062d-16*(3.09d21/3.1556926d13)**2/(xH*(1+cells(icells,3)) + 0.25*yHe*(1+cells(icells,4)+2*cells(icells,5)))
		
		do k=1,2
			call krome(densities(i,:), Tgas, dt)
		end do
	
	end do
	
	open(unit=rank,file='outstromgrenSi'//trim(ranktxt)//'.dat',status='replace',form='formatted',action='write')
	do i=1,size(densities,1)
		write(rank,fmt=*) rank, densities(i,:)
	end do
	
	t2=MPI_WTIME()
	write(*,*) rank,t2-t1
!	if (rank /= 0) then
!		call MPI_SEND(densities,size(densities,1)*nsp,MPI_DOUBLE_PRECISION,master,msgtag2,MPI_COMM_WORLD,ierr)
!	else
!		allocate(finaldensities(ncells,nsp))
!		finaldensities(1:size(densities,1),:) = densities(:,:)
!		do i=1,npsize-1
!			call MPI_RECV(densities,size(densities,1)*nsp,MPI_DOUBLE_PRECISION,i,msgtag2,MPI_COMM_WORLD,status0,ierr)
!			finaldensities(chunksize+1:chunksize+1+size(densities,1),:) = densities(:,:)
!			deallocate(finaldensities)
!		end do
!	end if
	
!	if(rank==master) then
!		open(unit=1,file="stromgrenSi.dat",status='replace',form='formatted',action='write')
!		do i=1,1000
!			write(1,fmt=*) densities(i,:)
!		end do
!	end if


	deallocate(cells,cell_pos,cell_level,densities)
	
	call MPI_FINALIZE(ierr)
end program ramses
