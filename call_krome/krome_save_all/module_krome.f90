module module_krome

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi_mine

  private
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [worker] of the parameter file
  ! --------------------------------------------------------------------------
  integer(kind=4),public    :: n_elements
  character(2000)           :: krome_parameter_file
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------


  integer(kind=4),allocatable,public          :: elements(:), n_ions(:), indices(:,:)
  real(kind=8),allocatable,public             :: abundances(:)
  real(kind=8),parameter,dimension(14),public :: solar_abundances = (/ 1d0, 8.51d-02, 1.12d-11, 2.40d-11, 5.01d-10, 2.69d-04, 6.76d-05, 4.90d-04, 3.63d-08, 8.51d-05, 1.74d-06, 3.98d-05, 2.82d-06, 3.24d-05 /)
  integer(kind=4),parameter,public            :: nIons = nmols - 6 - (natoms-3)
  character(2),parameter,dimension(14),public :: element_names = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si' /)
  character(3),parameter,dimension(6),public  :: roman_num     = (/ 'I  ', 'II ', 'III', 'IV ', 'V  ', 'VI ' /)


  public :: read_krome_params, get_non_zero_index

contains

  !*************************************************************************
  subroutine init_krome

    implicit none

    integer(kind=4)    :: i


    call krome_init

    if(n_elements /= natoms-3) then
       if(rank==0) print*, 'Problem, the number elements in the parameter file is in conflict with the chemical network used in Krome.'
       call stop_mpi
    end if

    if(n_elements < 1) then
       if(rank==0) print*, 'n_elements must be at least 1.  Stopping the program'
       call stop_mpi
    end if

    allocate(elements(n_elements), n_ions(n_elements), abundances(n_elements), indices(n_elements,5))

    call read_krome_params_file

    !Set to solar abundances if the user wrote a negative abundances
    do i=1,n_elements
       if(n_ions(i) < 1) then
          print*, 'Please put at least 1 for n_ions, stopping the program.'
          call stop_mpi
       end if
       if(abundances(i) < 0d0) then
          if(elements(i) > 14) then
             if(rank==0) print*, 'I did not implement the solar abundances for elements heavier than Silicone (14), please but an abundance in the parameter file'
             call stop_mpi
          else
             abundances(i) = solar_abundances(elements(i))
          end if
       end if
    end do

    do i=1,n_elements
       indices(i,:) = get_indices(i)
    end do


    if(verbose .and. rank==0) then
       print*, 'Krome initialized'
       print*, ' '
    end if

  end subroutine init_krome



  subroutine read_krome_params_file
    implicit none
    integer(kind=4) :: i

    open(unit=10, file=krome_parameter_file, status='old',action='read',form='formatted')
    do i=1,n_elements
       read(10,*) elements(i), n_ions(i), abundances(i)
    end do
    close(10)
  end subroutine read_krome_params_file
  !*************************************************************************


  !*************************************************************************
  function get_indices(int)

    implicit none

    integer(kind=4), intent(in) :: int
    integer(kind=4)             :: get_indices(5), i, j

    get_indices = 0
    get_indices(1) = 3 + int
    get_indices(2) = 5 + n_elements + int
    j=0
    do i=3,n_ions(int)+1
       get_indices(i) = 6 + 2*n_elements + int + j - count(mask=(n_ions(1:int-1)<i-1))
       j = j + count(mask=(n_ions>i-2))
    end do

  end function get_indices
  !*************************************************************************

  !*************************************************************************
  function get_non_zero_index(int,T,ion_state)

    implicit none

    integer(kind=4),intent(in)    :: int
    real(kind=8),intent(in)       :: T
    integer(kind=4),intent(out)   :: ion_state
    integer(kind=4)               :: get_non_zero_index, state

    ion_state = get_ion_state(int,T)

    get_non_zero_index = indices(int,ion_state)
    

  end function get_non_zero_index
  !*************************************************************************

  !*************************************************************************
  !Find the probable dominating ionization state of an element, depending on temperature, to accelerate the convergence to equilibrium of Krome
  function get_ion_state(int, T)

    implicit none

    integer(kind=4),intent(in) :: int
    real(kind=8),intent(in)    :: T
    integer(kind=4)            :: get_ion_state

    select case(elements(int))
    case(6)
       if(T<4.3d4) then
          get_ion_state = 2
       else if(T<1.1d5) then
          get_ion_state = 3
       else
          get_ion_state = 5
       end if

    case(8)
       if(T<1.6d4) then
          get_ion_state = 1
       else if(T<5.5d4) then
          get_ion_state = 2
       else if(T<1.13d5) then
          get_ion_state = 3
       else if(T<2.1d5) then
          get_ion_state = 4
       else
          get_ion_state = 5
       end if

    case(12)
       if(T<1.75d4) then
          get_ion_state = 2
       else if(T<1.25d5) then
          get_ion_state = 3
       else if(T<2.2d5) then
          get_ion_state = 4
       else
          get_ion_state = 5
       end if

    case(14)
       if(T<2d4) then
          get_ion_state = 2
       else if(T<7.8d4) then
          get_ion_state = 3
       else
          get_ion_state = 5
       end if

    case default
       print*, 'Other elements than Carbon(6), Oxygen(8), Magnesium(12) and Silicone(14) have not been implemented yet'
       call stop_mpi

    end select

    get_ion_state = min(get_ion_state, n_ions(int)+1)

  end function get_ion_state
  !*************************************************************************



  !*************************************************************************
  subroutine read_krome_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:7) == '[krome]') then
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
          case ('n_elements')
             read(value,*) n_elements
          case ('krome_parameter_file')
             write(krome_parameter_file,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)

    call init_krome

  end subroutine read_krome_params
  !*************************************************************************

end module module_krome



