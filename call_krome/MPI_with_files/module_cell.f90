module module_cell

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi_mine
  use module_krome
  use module_ramses
  use module_random


  implicit none

  type cell
     real(kind=8)                                :: nHI
     real(kind=8)                                :: nHII
     real(kind=8)                                :: nHeII
     real(kind=8)                                :: nHeIII
     real(kind=8)                                :: T
     real(kind=8)                                :: Z
     real(kind=8),dimension(nPhotoRea)           :: rates      !nPhotoRea is the number of photoionization reactions in the chemical network of Krome
     real(kind=8),dimension(nmols-natoms-3)      :: den_ions   !Density of the ionization stages of interest
  end type cell

  real(kind=8),allocatable,dimension(:,:),public :: csn
  type(cell),allocatable                         :: cellgrid(:)

  public :: init_cells, compute_cells, write_ion_files

contains

  subroutine init_csn(repository, snapnum, csn_file, nGroups, number_ions)

    implicit none

    character(2000),intent(in)              :: repository, csn_file
    integer(kind=4),intent(in)              :: snapnum, nGroups, number_ions
    integer(kind=4)                         :: i,j, test_nGroups, test_number_ions
    logical                                 :: exist

    inquire(file=csn_file, exist=exist)
    if(exist) then
       open(unit=10, file=csn_file, form='unformatted', action='read')
       read(10) test_nGroups
       if(test_nGroups /= nGroups) then
          print*, 'Problem, the number of photon groups in the csn_file is not correct. Check your parameters or your csn_file'
          call stop_mpi
       end if
       read(10) test_number_ions
       if(test_number_ions /= number_ions) then
          print*, 'Problem, the number of ions in the csn_file is not correct. Check your parameters or your csn_file'
          call stop_mpi
       end if

       allocate(csn(nGroups, number_ions))
       
       do i=1,number_ions
          read(10) (csn(j,i), j=1,nGroups)
       end do
       close(10)

       print*, 'csn (nSEDgroups * nIons)'
       do j=1,nGroups
          print*, csn(j,:)
       end do

    else
       print*, 'Beginning computation of csn'
       call compute_csn_in_box(repository, snapnum, n_elements, elements, n_ions, csn)

       open(unit=10, file=csn_file, form='unformatted', action='write')
       write(10) nGroups
       write(10) number_ions
       do i=1,number_ions
          write(10) (csn(j,i), j=1,nGroups)
       end do
       close(10)
       
    end if

  end subroutine init_csn


  subroutine init_cells(repository, snapnum, nvar, ncell, ramses_var, restart)

    use module_ramses
    use module_krome
    use module_spectra

    implicit none

    character(2000),intent(in)                    :: repository
    integer(kind=4),intent(in)                    :: snapnum, nvar, ncell
    real(kind=8),intent(in),dimension(ncell,nvar) :: ramses_var
    logical,intent(in)                            :: restart
    real(kind=8),allocatable,dimension(:,:)       :: cells_rt
    real(kind=8),dimension(ncell)                 :: nH, Tgas, mets, nHI
    real(kind=8),dimension(3,ncell)               :: fractions
    integer(kind=4)                               :: nSEDgroups, i, j, k, seed
    real(kind=8)   :: ran

    nSEDgroups = get_nOptBins()
    allocate(cells_rt(nSEDgroups,ncell))

    call ramses_get_nh_cgs(repository,snapnum,ncell,nvar,ramses_var,nH)
    call ramses_get_T_nhi_cgs(repository,snapnum,ncell,nvar,ramses_var,Tgas,nHI)
    call ramses_get_metallicity(ncell,nvar,ramses_var,mets)
    call ramses_get_fractions(ncell,nvar,ramses_var,fractions)
    call ramses_get_flux(repository,snapnum,ncell,nvar,nSEDgroups,ramses_var,cells_rt)

    allocate(cellgrid(ncell))
    do i=1,ncell
       cellgrid(i)%nHI = nH(i)*(1d0-fractions(1,i))
       cellgrid(i)%nHII = nH(i)*fractions(1,i)
       cellgrid(i)%nHeII = 7.895d-2*nH(i)*fractions(2,i)
       cellgrid(i)%nHeIII = 7.895d-2*nH(i)*fractions(3,i)
       cellgrid(i)%T = Tgas(i)
       cellgrid(i)%Z = mets(i)
       do j=1,nPhotoRea
          cellgrid(i)%rates(j) = sum(cells_rt(:,i)*csn(:,j))
       end do
       if((.not. restart) .and. (elements(1)/=8)) cellgrid(i)%rates(1) = 1.65d-9
       
       !HACKKKK, HARD CODE
       !Put UVB instead of stellar radiation
       ! cellgrid(i)%rates(1) = 1.578d-9
       ! cellgrid(i)%rates(2) = 1.212d-13
       
       cellgrid(i)%den_ions(:) = 0d0
    end do

    deallocate(cells_rt)

  end subroutine init_cells


  subroutine compute_cells(icpu)

    implicit none

    integer(kind=4),intent(in)         :: icpu
    integer(kind=4)                    :: i, j, k, l, ion_state(n_elements), ncell, non_zero_index(n_elements)
    real(kind=8)                       :: densities(nmols), n_ion_save(n_elements)

    ncell = size(cellgrid)

    do i=1,ncell

       if(cellgrid(i)%T > 1d6) then
          cellgrid(i)%den_ions(:) = 1d-18
       else

          call krome_set_photoBin_rates(cellgrid(i)%rates)

          densities(:)               = 1d-18
          densities(krome_idx_H)     = max(cellgrid(i)%nHI,1d-18)
          densities(krome_idx_Hj)    = max(cellgrid(i)%nHII,1d-18)
          densities(krome_idx_He)    = max(7.895d-2*(cellgrid(i)%nHI + cellgrid(i)%nHII) - cellgrid(i)%nHeII - cellgrid(i)%nHeIII,1d-18)
          densities(krome_idx_Hej)   = max(cellgrid(i)%nHeII,1d-18)
          densities(krome_idx_Hejj)  = max(cellgrid(i)%nHeIII,1d-18)
          densities(krome_idx_E)     = densities(krome_idx_Hj) + densities(krome_idx_Hej) + 2*densities(krome_idx_Hejj)


          do j=1,n_elements
             n_ion_save(j) = cellgrid(i)%Z/0.0134*abundances(j)*(densities(krome_idx_H) + densities(krome_idx_Hj))
             non_zero_index(j) = get_non_zero_index(j,cellgrid(i)%T,ion_state(j))
             densities(non_zero_index(j)) = max(n_ion_save(j), 1d-18)
             densities(krome_idx_E) = densities(krome_idx_E) + (ion_state(j)-1)*densities(non_zero_index(j))
          end do
          densities(krome_idx_E) = max(densities(krome_idx_E),1d-18)

          call krome_equilibrium(densities, cellgrid(i)%T, icpu)

          l=0
          do j=1,n_elements
             !Rare cases of bugs :
             if(maxval(densities(indices(j,1:n_ions(j)))) > 1.0001*n_ion_save(j)) then
                print*, 'bug in file ', icpu
                print*, 'temperature, nHI, nHII, photorates, metallicity'
                print*, cellgrid(i)%T, cellgrid(i)%nHI, cellgrid(i)%nHII, cellgrid(i)%rates(:), cellgrid(i)%Z
                cellgrid(i)%den_ions(l+1:l+n_ions(j)) = 1d-18
                cellgrid(i)%den_ions(l+ion_state(j)) = max(n_ion_save(j),1d-18)
             !No bug
             else
                do k=1,n_ions(j)
                   if(indices(j,k) < 1) then
                      print*, 'Problem ! Trying to update a ionization state that does not exist'
                      call stop_mpi
                   end if
                   cellgrid(i)%den_ions(k+l) = max(densities(indices(j,k)), 1d-18)
                end do
             end if
             l = l + n_ions(j)
          end do
       end if
    end do
  end subroutine compute_cells


  subroutine write_ion_files(snapnum, icpu, output_path)
    
    implicit none

    character(2000),intent(in)              :: output_path
    integer(kind=4),intent(in)              :: snapnum, icpu
    character(1000)                         :: nomfich
    integer(kind=4)                         :: i,j,k,l, ncell

    ncell = size(cellgrid)

    l=0
    do j=1,n_elements
       do k=1,n_ions(j)
          write(nomfich,'(a,a,a,a,a,i5.5,a,i5.5)') trim(output_path),'/',trim(element_names(elements(j))),trim(roman_num(k)),'_',snapnum,'.out',icpu
          open(unit=10, file=nomfich, form='unformatted', action='write')
          write(10) ncell
          write(10) (cellgrid(i)%den_ions(k+l), i=1,ncell)
          close(10)
       end do
       l = l + n_ions(j)
    end do

    deallocate(cellgrid)

  end subroutine write_ion_files
    
end module module_cell


