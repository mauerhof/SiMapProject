module module_cell

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi_mine
  use module_krome
  use module_ramses


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

  subroutine init_csn(repository, snapnum)

    implicit none

    character(2000),intent(in)              :: repository
    integer(kind=4),intent(in)              :: snapnum

    !allocate(csn(3,2))
    !csn = 1d-18
    call compute_csn_in_box(repository, snapnum, n_elements, elements, n_ions, csn)

  end subroutine init_csn


  subroutine init_cells(repository, snapnum, nvar, ncell, ramses_var)

    use module_ramses
    use module_krome

    implicit none

    character(2000),intent(in)                    :: repository
    integer(kind=4),intent(in)                    :: snapnum, nvar, ncell
    real(kind=8),intent(in),dimension(ncell,nvar) :: ramses_var
    real(kind=8),allocatable,dimension(:,:)       :: cells_rt
    real(kind=8),dimension(ncell)                 :: nH, Tgas, mets, nHI
    real(kind=8),dimension(3,ncell)               :: fractions
    integer(kind=4)                               :: nSEDgroups, i, j, k

    nSEDgroups = get_nSEDgroups(repository,snapnum)
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
       !!!!!!!!! To remove after succesful restarts of simulations with Ramses, to add a "low-energy" photon bin'
       do j=1,n_elements
          if(elements(j) /= 8) then
             cellgrid(i)%rates(sum(n_ions(1:j-1))+1) = 1d-10
          end if
       end do
       !!!!!!!!!
       cellgrid(i)%den_ions(:) = 0d0
    end do

    deallocate(cells_rt)

  end subroutine init_cells


  subroutine compute_cells()

    implicit none

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

          call krome_equilibrium(densities, cellgrid(i)%T)

          l=0
          do j=1,n_elements
             !Rare cases of bugs :
             if(maxval(densities(indices(j,1:n_ions(j)))) > 1.0001*n_ion_save(j)) then
                print*, 'bug'
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
          write(nomfich,'(a,a,i5.5,a,a,a,a,i5.5,a,i5.5)') trim(output_path),'/output_',snapnum,'/',trim(element_names(elements(j))),trim(roman_num(k)),'_',snapnum,'.out',icpu
          open(unit=10, file=nomfich, form='unformatted', action='write')
          write(10) ncell
          ! write(10) (cellgrid(i)%den_ions(k+l), i=1,ncell)
          write(10) (cellgrid(i)%nHI+cellgrid(i)%nHII, i=1,ncell)
          close(10)
       end do
       l = l + n_ions(j)
    end do

    deallocate(cellgrid)

  end subroutine write_ion_files
    
end module module_cell

