PROGRAM computeSiFractions
	use py_cell_utils
!~ 	use mpi
	!use SED_module,only:getSEDcsn
	
	
	implicit none
	
	!Variables for mpi
!	integer(kind=4)	 			:: ierr, rank, npsize, status0(MPI_STATUS_SIZE), chunksize
	integer,parameter			:: master=0
	real(kind=8)	 			:: t1=0d0, t2=0d0
	character*10     			:: ranktxt
	
	
	!Variables for Cloudy
	integer,parameter			:: nBinC=4, nIon=nBinC !(nIon is not strictly useful but helps understand)  
	real(kind=8) 				:: EC(nBinC), cloudyJ(nBinC), nH, nHe, Tgas, cloudyCoeffs(nIon,nBinC), laserWidth
	character(len=1000)			:: cloudyString
	
	
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
	real(kind=8)				:: xH, yHe, ev_to_hz, ev_to_ryd, PI, testRates(7)
	integer(kind=4)				:: i , j, k, icells, l, orderIonE(nIon)!, orderBins(nBins-3)
	character*10				:: icellstxt
	
	
	!Commands related to mpi
!~ 	call MPI_INIT(ierr)
!~ 	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
!~ 	call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
!~ 	t1=MPI_WTIME()										!Initialize time to display the total time of execution
!~ 	open(unit=14,status='scratch')
!~     rewind(14)
!~     write (14,*) rank
!~     rewind(14)											!5 lines to create a character of the rank of the core,  to put it in names of output files   (more straightforward way ?)
!~     read  (14,*) ranktxt
!~     close (14)
	

	!!!Can obviously improve this part to import numbers from ramses output,  would be more elegant
	ev_to_hz = 2.41799d14 ; ev_to_ryd = 7.34989d-2 ; PI=3.14159266
	xH = 0.76 ; yHe = 0.24 ; binsR = (/ 8.15d0, 1.36d1, 1.635d1, 2.459d1, 3.349d1, 4.514d1, 5.442d1, 1d3 /)  !Maybe should use infinity instead of 1d3,  changes almost nothing
	unit_length = 3.09d21 ; unit_time = 3.1556926d13
	

	!sed_folder = '/Users/mauerhof/Documents/seds/bpass100/'
	!call readSeds(sed_folder,ages,Zs,Ls,SEDs)
	
	
	orderIonE = (/ 4, 5, 6, 7 /)  !Indices of (Si) ions in rt_parameters

!~ 	do j=1,nIon
!~ 		sigmaN(j,:) = (/ (getSEDcsn(Ls, SEDs(:,10,2), size(Ls), binsR(i), binsR(i+1), orderIonE(j)), i=1,nBinR) /)    !!!!SiI and SiII show a rise of sigmaN, weird.   sigmaN is the number weighted average of the cross section over each bin, for each species of Si.   (SEDs(:,10,2) is the SED for bpass100 with z=0.002 and age=6.3 Myr (to check)
!~ 	end do
	open(unit = 2 , file = 'SiCSN.dat')
	do i = 1,4
		read(2,*) sigmaN(i,:)
	end do
	
	print*, sigmaN
	
!~ 	ramDir = '/Users/mauerhof/Documents/RamsesFiles/stromgren/'
	ramDir = '../'
	ts = 16; levelmax = 13 ; center(:) = (/ 5d-1, 5d-1, 5d-1/) ; radius = 1d-1 ; readRT = .true.
	hydroIndices = (/ 1, 6, 7, 8, 9 /)   !Density, pressure, xHII, xHeII, xHeIII
	binIndices = (/ (10+i*4, i=0,nBinR-1) /) !nBinR photon Bins
	usefulIndices(1:size(hydroIndices)) = hydroIndices
	usefulIndices(size(hydroIndices)+1:size(hydroIndices)+size(binIndices)) = binIndices
	
	
!~ 	call count_cells(ramDir,ts,levelmax,ncells,center,radius)   !return the total number of cells 
	
!~ 	chunksize = (ncells+npsize-1)/npsize							!Size of the array for each core.  A little bit smaller for the last core.
!~ 	if (rank==npsize-1) chunksize = ncells - (npsize-1)*chunksize	
	

	
!~ 	if(rank==master) then		!The master reads all the data
!~ 		allocate(cells(ncells,nfields))
!~ 		allocate(cell_pos(ncells,3))
!~ 		allocate(cell_level(ncells))
!~ 		call read_cells_hydro(ramDir,ts,levelmax,cells,cell_pos,cell_level,ncells,usefulIndices,nfields,center,radius,readRT)
		

!~ 		do i=1,npsize-2			!Then the master sends each piece of data to the other threads
!~ 			call MPI_SEND(cells(i*chunksize+1:(i+1)*chunksize,:),chunksize*nfields,MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,ierr)
!~ 		end do
!~ 		call MPI_SEND(cells((npsize-1)*chunksize+1:ncells,:),(ncells-(npsize-1)*chunksize)*nfields,MPI_DOUBLE_PRECISION,npsize-1,11,MPI_COMM_WORLD,ierr)  !The last thread data has another size
!~ 	else
!~ 		allocate(cells(chunksize,nfields))		!The other cores receive their piece of data
!~ 		call MPI_RECV(cells,chunksize*nfields,MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
!~ 	end if
	
!~         laserWidth = 5d-3
!~         EC(:) = (/ 9d0, 1.7d1, 3.4d1, 4.6d1 /)
!~         do i=1,nBinC
!~            cloudyCoeffs(i,:) = (/ (getSEDcsn(Ls, SEDs(:,10,2), size(Ls), EC(j)-laserWidth*EC(j), EC(j)+laserWidth*EC(j), orderIonE(i))/(EC(j)*ev_to_erg), j=1,nBinC) /)
!~         print*, cloudyCoeffs(i,:)
!~         end do

!~ 	!open(unit=rank+20,file='rates'//trim(ranktxt)//'.dat',status='replace',form='formatted',action='write')   !To check if the rates are reproduced correctly
	
	
!~ 	!Loop over all the cells
!~ 	do icells=1,5
!~ 		open(unit=icells,file='./cloudy.in',status='replace',form='formatted',action='write')  !To write the outputs
!~ 		open(unit=icells+10000000,status='scratch')
!~ 		rewind(14)
!~ 		write (14,*) icells
!~ 		rewind(14)											!5 lines to create a character of the rank of the core,  to put it in names of output files   (more straightforward way ?)
!~ 		read  (14,*) icellstxt
!~ 		close (14)
!~ 		!print*, icellstxt
!~ 		ratesRamses(:) = (/ (sum(cells(icells,size(hydroIndices)+1:size(hydroIndices)+nBinR)*sigmaN(i,:))*unit_length/unit_time, i=1,nIon) /)  !Photoionization rate of each Si ion, computed from ramses output. This is the inhomogeneous part of the linear system cited above

!~ 		call solveTriangularMatrix(nIon,ratesRamses(:),cloudyCoeffs(:,:),cloudyJ(:))   !Solves the linear system to find to correct Krome inputs
!~ 		print*,cloudyJ(:)

!~ 		nH = xH*cells(icells,1)		
!~ 		nHe = 0.25*yHe*cells(icells,1)		
		
!~ 	!	densities(:) = (/ nH*cells(icells,3) + nHe*(cells(icells,4)+2*cells(icells,5)), nH*(1-cells(icells,3)), nHe*(1-cells(icells,4)-cells(icells,5)), 3.467d-6, nH*cells(icells,3), nHe*cells(icells,4), 1d-20, nHe*cells(icells,5), 1d-20, 1d-20, 1d-20 /) !(E, HI, HeI, SiI, HII, HeII, SiII, HeIII, SiIII, SiIV, SiV), as needed for Krome

!~ 		Tgas = cells(icells,2)/cells(icells,1)*1.66d-24/1.38062d-16*(unit_length/unit_time)**2/(xH*(1+cells(icells,3)) + 0.25*yHe*(1+cells(icells,4)+2*cells(icells,5)))  !Trust me
!~ 		!print*,Tgas
		
!~ 		write(icells,"(a)") 'hden -1'
!~ 		write(icells,"(a)")'abundances he =-1.1 li =-40 be =-40 b =-40 c =-40 n =-40 o =-40'
!~ 		write(icells,"(a)") 'continue f =-40 ne =-40 na =-40 mg =-40 '
!~ 		write(icells,"(a)") 'continue al =-40 si =-4.46 p =-40 s =-40 cl=-40 ar=-40 k =-40 '
!~ 		write(icells,"(a)") 'continue ca =-40 sc =-40 ti =-40 v =-40 cr =-40 mn =-40 fe =-40 '
!~ 		write(icells,"(a)") 'continue co =-40 ni =-40 cu =-40 zn =-40 '
!~ 		write(icells,"(a,F10.3,F10.3)") 'element hydrogen ionization ', 1 - cells(icells,3), cells(icells,3)
!~ 		write(icells,"(a,F10.3,F10.3,F10.3)") 'element helium ionization ', 1 - cells(icells,4) - cells(icells,5), cells(icells,4), cells(icells,5)
!~ 		do i=1,nBinC
!~ 			write(icells,"(a,F10.3,a)") 'Laser, frequency = ', EC(i)*ev_to_ryd, ' ryd rel width 0.005'
!~ 			write(icells,"(a,ES10.3,a)") 'intensity total ', cloudyJ(i), ' linear'
!~ 		end do
!~ 		write(icells,"(a,F10.3)") 'const temp', log10(Tgas)
!~ 		write(icells,"(a)") 'sphere'
!~ 		write(icells,"(a)") 'stop zone 1'
!~ 		write(icells,"(a)") 'save averages, file="ion'//trim(icellstxt)//'.avr", print last iteration'
!~ 		write(icells,"(a)") 'ionization, silicone 1'//NEW_LINE('')//'ionization, silicone 2'//NEW_LINE('')//'ionization, silicone 3'//NEW_LINE('')//'ionization, silicone 4'//NEW_LINE('')//'ionization, silicone 5'//NEW_LINE('')//'end of averages'
!~ 		close(icells)
!~ 		call system('/Users/mauerhof/Documents/c17.00/source/cloudy.exe < cloudy.in > cloudy.out')
		!CloudyString = 'hden -1 '//NEW_LINE('')// 'abundances he =-1.1 li =-40 be =-40 b =-40 c =-40 n =-40 o =-40 '//NEW_LINE('')// 'continue f =-40 ne =-40 na =-40 mg =-40 '//NEW_LINE('')// 'continue al =-40 si =-4.46 p =-40 s =-40 cl=-40 ar=-40 k =-40 '//NEW_LINE('')// 'continue ca =-40 sc =-40 ti =-40 v =-40 cr =-40 mn =-40 fe =-40 '//NEW_LINE('')// 'continue co =-40 ni =-40 cu =-40 zn =-40 '//NEW_LINE('')// 'element hydrogen ionization ' !+ 1-cells(icells,3) + ' ' + (data_ramses[j,1]) + ''//NEW_LINE('')// 'element helium ionization ' + (1-data_ramses[j, 2]-data_ramses[j,3]) + ' ' + (data_ramses[j,2]) + ' ' + (data_ramses[j,3]) + '\n'
!		for i in range(numberofBins):
!		CloudyString = CloudyString + 'Laser, frequency = ' + (energiesRyd[i]) + ' ryd rel width 0.005 '//NEW_LINE('')// 'intensity total ' + (9.79183e7*data_ramses[j,4+i]*energiesErg[i]) + ' linear \n'
!		CloudyString = CloudyString + 'constant temp 4 '//NEW_LINE('')// 'sphere '//NEW_LINE('')// 'stop zone 1 '//NEW_LINE('')// 'iterate to convergence max=3, error=0.1 '//NEW_LINE('')// 'save averages, file=\'ion'+str(j+1)+'.avr\', print last iteration '//NEW_LINE('')// 'ionization, hydrogen 1 '//NEW_LINE('')// 'ionization, hydrogen 2 '//NEW_LINE('')// 'ionization, helium 1 '//NEW_LINE('')// 'ionization, helium 2 '//NEW_LINE('')// 'ionization, helium 3 '//NEW_LINE('')// 'ionization, silicone 1 '//NEW_LINE('')// 'ionization, silicone 2 '//NEW_LINE('')// 'ionization, silicone 3 '//NEW_LINE('')// 'ionization, silicone 4 '//NEW_LINE('')// 'ionization, silicone 5 '//NEW_LINE('')// 'ionization, silicone 6 '//NEW_LINE('')// 'end of averages'
		
!		testRates = krome_get_photoBin_rates()  !To check if the rates are reproduced correctly
!		write(rank+20,fmt=*) ratesRamses
!		write(rank+20,fmt=*) testRates(4:7)
!		write(rank+20,fmt=*) ''
		
!		write(rank,fmt=*) densities(4),densities(7),densities(9),densities(10),densities(11)   !Writes output for one cell
!~ 	end do
	!cloudyString = 'test a'//NEW_LINE('')//'test '
	!print*,cloudyString
!	close(rank) !; close(rank+20)
!~ 	deallocate(ages, Zs, Ls, SEDs)
!	if (rank==master) then
!~ 		deallocate(cells,cell_pos,cell_level)
!	else
!		deallocate(cells)
!	end if
	
!	t2=MPI_WTIME()
!	write(*,*) rank,t2-t1

!	call MPI_FINALIZE(ierr)
	
	
	

END PROGRAM

!~ SUBROUTINE readSeds(path,ages,Zs,Ls,SEDs)

!~ 	implicit none
	
!~ 	character(len=128),intent(in)::path
!~ 	integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
!~ 	real(kind=8),allocatable,intent(out)::ages(:), Zs(:), Ls(:)
!~ 	real(kind=8),allocatable,intent(out)::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
!~ 	!real(kind=8),allocatable::tbl(:,:,:), tbl2(:,:,:), reb_tbl(:,:,:)
!~ 	integer::i,ia,iz,ip,ii,dum
!~ 	character(len=128)::fZs, fAges, fSEDs                        ! Filenames
!~ 	logical::ok,okAge,okZ
	
	
!~ 	!if(path.eq.'')call get_environment_variable('RAMSES_SED_DIR',path)
!~ 	inquire(FILE=TRIM(path)//'/all_seds.dat', exist=ok)
!~ 	if(.not. ok)then
!~ 		!if(myid.eq.1) then 
!~ 			write(*,*)'Cannot access SED directory ',TRIM(path)
!~ 			write(*,*)'Directory '//TRIM(path)//' not found'
!~ 			write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
!~ 			' to the correct path, or use the namelist.'
!~ 		!endif
!~ 		!call clean_stop
!~ 	end if
!~ 	write(fZs,'(a,a)')   trim(path),"/metallicity_bins.dat"
!~ 	write(fAges,'(a,a)') trim(path),"/age_bins.dat"
!~ 	write(fSEDs,'(a,a)') trim(path),"/all_seds.dat"
!~ 	inquire(file=fZs, exist=okZ)
!~ 	inquire(file=fAges, exist=okAge)
!~ 	inquire(file=fSEDs, exist=ok)
!~ 	if(.not. ok .or. .not. okAge .or. .not. okZ) then
!~ 		!if(myid.eq.1) then 
!~ 			write(*,*) 'Cannot read SED files...'
!~ 			write(*,*) 'Check if SED-directory contains the files ',  &
!~ 				'metallicity_bins.dat, age_bins.dat, all_seds.dat'
!~ 		!endif
!~ 			!call clean_stop
!~ 	end if

!~ ! READ METALLICITY BINS-------------------------------------------------
!~ 	open(unit=10,file=fZs,status='old',form='formatted')
!~ 	read(10,'(i8)') nzs
!~ 	allocate(zs(nzs))
!~ 	do i = 1, nzs
!~ 		read(10,'(e14.6)') zs(i)
!~ 	end do
!~ 	close(10)
!~ ! READ AGE BINS---------------------------------------------------------
!~ 	open(unit=10,file=fAges,status='old',form='formatted')
!~ 	read(10,'(i8)') nAges
!~ 		allocate(ages(nAges))
!~ 	do i = 1, nAges
!~ 		read(10,'(e14.6)') ages(i)
!~ 	end do
!~ 	close(10)
!~ 	ages = ages*1.e-9                       !         Convert from yr to Gyr
!~ 	if(ages(1) .ne. 0.) ages(1) = 0.
!~ ! READ SEDS-------------------------------------------------------------
!~ 	open(unit=10,file=fSEDs,status='old',form='unformatted')
!~ 	read(10) nLs, dum
!~ 	allocate(Ls(nLs))
!~ 	read(10) Ls(:)
!~ 	allocate(SEDs(nLs,nAges,nzs))
!~ 	do iz = 1, nzs
!~ 		do ia = 1, nAges
!~ 			read(10) SEDs(:,ia,iz)
!~ 		end do
!~ 	end do
!~ 	close(10)
	
!~ 	return
	
!~ 	deallocate(ages, Zs, Ls, SEDs)
	

!~ END SUBROUTINE





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
	
END SUBROUTINE



