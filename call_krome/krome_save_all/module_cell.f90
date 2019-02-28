module module_cell

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi_mine
  use module_krome

  implicit none

  type cell
     integer(kind=4)                             :: ID
     real(kind=8)                                :: nHI
     real(kind=8)                                :: nHII
     real(kind=8)                                :: nHeII
     real(kind=8)                                :: nHeIII
     real(kind=8)                                :: T
     real(kind=8)                                :: Z
     real(kind=8),dimension(nPhotoRea)           :: rates      !nPhotoRea is the number of photoionization reactions in the chemical network of Krome
     real(kind=8),dimension(nmols-natoms-3)      :: den_ions   !Density of the ionization stages of interest
  end type cell

  public :: init_cells, compute_cells, save_cells

contains  


  subroutine compute_cells(nbundle, cellpacket)

    implicit none

    integer(kind=4),intent(in)                   :: nbundle
    type(cell),dimension(nbundle),intent(inout)  :: cellpacket
    integer(kind=4)                              :: i, j, k, l, ion_state
    real(kind=8)                                 :: densities(nmols)

    do i=1,nbundle
       if(cellpacket(i)%ID > 0) then

          if(cellpacket(i)%T > 1d6) then
             cellpacket(i)%den_ions(:) = 1d-18
          else

             call krome_set_photoBin_rates(cellpacket(i)%rates)

             densities(:) = 1d-18
             densities(krome_idx_H)     = max(cellpacket(i)%nHI,1d-18)
             densities(krome_idx_Hj)    = max(cellpacket(i)%nHII,1d-18)
             densities(krome_idx_He)    = max(7.895d-2*(cellpacket(i)%nHI + cellpacket(i)%nHII) - cellpacket(i)%nHeII - cellpacket(i)%nHeIII,1d-18)
             densities(krome_idx_Hej)   = max(cellpacket(i)%nHeII,1d-18)
             densities(krome_idx_Hejj)  = max(cellpacket(i)%nHeIII,1d-18)
             densities(krome_idx_E) = densities(krome_idx_Hj) + densities(krome_idx_Hej) + 2*densities(krome_idx_Hejj)

             do j=1,n_elements
                k = get_non_zero_index(j,cellpacket(i)%T,ion_state)
                densities(k) = max(cellpacket(i)%Z/0.0134*abundances(j)*(densities(krome_idx_H) +  densities(krome_idx_Hj)), 1d-18)
                densities(krome_idx_E) = densities(krome_idx_E) + (ion_state-1)*densities(k)
             end do
             densities(krome_idx_E) = max(densities(krome_idx_E),1d-18)

             call krome_equilibrium(densities, cellpacket(i)%T)

             l=0
             do j=1,n_elements
                do k=1,n_ions(j)
                   if(indices(j,k) < 1) then
                      print*, 'Problem ! Trying to update a ionization state that does not exist'
                      call stop_mpi
                   end if
                   cellpacket(i)%den_ions(k+l) = max(densities(indices(j,k)), 1d-18)
                end do
                l = l + n_ions(j)
             end do
          end if
       end if
    end do

  end subroutine compute_cells



  subroutine init_cells(repository, snapnum, reading_method, cell_data_file, cellgrid, fileout)

    use module_ramses
    use module_krome

    implicit none

    character(2000),intent(in)              :: repository, reading_method, cell_data_file, fileout
    integer(kind=4),intent(in)              :: snapnum
    type(cell),allocatable,intent(inout)    :: cellgrid(:)
    integer(kind=4),allocatable             :: cpu_list(:), cell_level(:)
    real(kind=8),allocatable,dimension(:,:) :: csn, cells_rt, fractions, ions, cell_pos, cells, velocities
    real(kind=8),allocatable,dimension(:)   :: nH, Tgas, mets, nHI
    integer(kind=4)                         :: nSEDgroups, ncells, nfields, ncpu_read, i, j, k


    call compute_csn_in_box(repository, snapnum, n_elements, elements, n_ions, csn)

    ! -------------------- read/compute the data ------------------------
    if(reading_method=='fromlist') then
       open(unit=10, file=trim(cell_data_file), form='unformatted', action='read')
       read(10) ncells
       read(10) nSEDgroups
       allocate(nH(ncells), Tgas(ncells), mets(ncells), cells_rt(nSEDgroups,ncells), fractions(3,ncells))
       read(10) nH
       read(10) Tgas
       read(10) mets
       read(10) cells_rt
       read(10) fractions
       close(10)
    else
       nSEDgroups = get_nSEDgroups(repository,snapnum)
       ncpu_read = get_ncpu(repository,snapnum)
       allocate(cpu_list(1:ncpu_read))
       cpu_list(:) = (/ (i, i=1,ncpu_read) /)
       call read_leaf_cells_omp(repository, snapnum, ncpu_read, cpu_list, ncells, nfields, cell_pos, cells, cell_level)
       !call read_leaf_cells(repository, snapnum, ncells, nfields, cell_pos, cells, cell_level)
       allocate(nH(ncells), nHI(ncells), Tgas(ncells), mets(ncells), cells_rt(nSEDgroups,ncells), fractions(3,ncells), velocities(3,ncells))
       call ramses_get_velocity_cgs(repository,snapnum,ncells,nfields,cells,velocities)
       call ramses_get_nh_cgs(repository,snapnum,ncells,nfields,cells,nH)
       call ramses_get_T_nhi_cgs(repository,snapnum,ncells,nfields,cells,Tgas,nHI) ; deallocate(nHI)
       call ramses_get_metallicity(ncells,nfields,cells,mets)
       call ramses_get_fractions(ncells,nfields,cells,fractions)
       call ramses_get_flux(repository,snapnum,ncells,nfields,nSEDgroups,cells,cells_rt)

       call dump_cells_v_pos_l(fileout, cell_pos, cell_level, velocities)
       deallocate(cell_pos, cell_level, velocities)
    endif
    ! ------------------------------------------------------------------------------------------------

    allocate(cellgrid(ncells))
    do i=1,ncells
       cellgrid(i)%ID = i
       cellgrid(i)%nHI = nH(i)*(1d0-fractions(1,i))
       cellgrid(i)%nHII = nH(i)*fractions(1,i)
       cellgrid(i)%nHeII = 7.895d-2*nH(i)*fractions(2,i)
       cellgrid(i)%nHeIII = 7.895d-2*nH(i)*fractions(3,i)
       cellgrid(i)%T = Tgas(i)
       cellgrid(i)%Z = mets(i)
       do j=1,nPhotoRea
          cellgrid(i)%rates(j) = sum(cells_rt(:,i)*csn(:,j))
       end do
       cellgrid(i)%den_ions(:) = 0d0
    end do

    deallocate(nH, fractions, Tgas, mets, cells_rt)

    call save_cells(fileout, cellgrid, .false., 0)
    

  end subroutine init_cells



  subroutine dump_cells_v_pos_l(file, cell_pos, cell_l, v)

    implicit none

    character(2000),intent(in)              :: file
    real(kind=8),dimension(:,:),intent(in)  :: cell_pos, v
    integer(kind=4),dimension(:),intent(in) :: cell_l
    integer(kind=4)                         :: i,j,ncells

    ncells = size(cell_l)

    open(unit=14, file=trim(file)//'.pos_l_v', status='unknown', form='unformatted', action='write')
    write(14) ncells
    write(14) (cell_pos(i,1), i=1,ncells)
    write(14) (cell_pos(i,2), i=1,ncells)
    write(14) (cell_pos(i,3), i=1,ncells)
    write(14) (cell_l(i),     i=1,ncells)
    write(14) (v(1,i),        i=1,ncells)
    write(14) (v(2,i),        i=1,ncells)
    write(14) (v(3,i),        i=1,ncells)
    close(14)
  end subroutine dump_cells_v_pos_l


  subroutine save_cells(file,cgrid,not_first_save,ncells_done)

    character(2000),intent(in)                   :: file
    type(cell),dimension(:),intent(in)           :: cgrid
    logical,intent(in)                           :: not_first_save
    integer(kind=4),intent(in)                   :: ncells_done
    logical                                      :: file_exists
    integer(kind=4)                              :: i,j,k,l,ncells
    character(2000)                              :: file_ion, filebak, command

    ncells = size(cgrid)

    if(not_first_save) then
       l=0
       do i=1,n_elements
          do j=1,n_ions(i)
             file_ion = trim(file)//'.'//trim(element_names(elements(i)))//trim(roman_num(j))
             INQUIRE(FILE=file_ion, EXIST=file_exists)
             if(file_exists) then
                filebak = trim(file_ion)//'.bak'
                command = 'mv '//trim(file_ion)//' '//trim(filebak)
                call system(command)
             end if
             open(unit=16, file=file_ion, status='unknown', form='unformatted', action='write')
             write(16) ncells_done
             write(16) (cgrid(k)%den_ions(l+j), k=1,ncells_done)
             close(16)
          end do
          l = l + n_ions(i)
       end do
    end if

    !To accelerate restarts, save all the properties of the cells
    if(.not. not_first_save) then
       open(unit=16, file=trim(file)//'.cells', status='unknown', form='unformatted', action='write')
       write(16) ncells
       write(16) (cgrid(i)%nHI,          i=1,ncells)
       write(16) (cgrid(i)%nHII,         i=1,ncells)
       write(16) (cgrid(i)%nHeII,        i=1,ncells)
       write(16) (cgrid(i)%nHeIII,       i=1,ncells)
       write(16) (cgrid(i)%T,            i=1,ncells)
       write(16) (cgrid(i)%Z,            i=1,ncells)
       write(16) (cgrid(i)%rates(:),     i=1,ncells)
       close(16)
    end if

  end subroutine save_cells
  

  subroutine restore_cells(file,cgrid)

    character(2000),intent(in)                                 :: file
    type(cell),dimension(:),allocatable, intent(out)           :: cgrid
    integer(kind=4)                                            :: i,j,k,l, ncells, ncells_done

  
    ! restore cells from saved file
    open(unit=16, file=trim(file)//'.cells', status='old', form='unformatted', action='read')
    read(16) ncells
    allocate(cgrid(ncells))

    do i=1,ncells
       cgrid(i)%ID = i
    end do
    read(16) (cgrid(i)%nHI,          i=1,ncells)
    read(16) (cgrid(i)%nHII,         i=1,ncells)
    read(16) (cgrid(i)%nHeII,        i=1,ncells)
    read(16) (cgrid(i)%nHeIII,       i=1,ncells)
    read(16) (cgrid(i)%T,            i=1,ncells)
    read(16) (cgrid(i)%Z,            i=1,ncells)
    read(16) (cgrid(i)%rates(:),     i=1,ncells)
    close(16)

    l=0
    do i=1,n_elements
       do j=1,n_ions(i)
          open(unit=16, file=trim(file)//'.'//trim(element_names(elements(i)))//trim(roman_num(j)), status='old', form='unformatted', action='read')
          read(16) ncells_done
          read(16) (cgrid(k)%den_ions(l+j), k=1,ncells_done)
          close(16)
       end do
       l = l + n_ions(i)
    end do
    do i=ncells_done+1,ncells
       cgrid(i)%den_ions(:) = 0d0
    end do

  end subroutine restore_cells





end module module_cell


