!Program that writes initial conditions with photons coming from a central point source, with random directions
program write_IC

  use module_constants
  use module_random
  use module_utils

  integer(kind=4)             :: nphoton, iseed, i, j, k, nquasars
  real(kind=8)                :: totalflux   !nb of real photons [# / s]
  integer(kind=4),allocatable :: ID(:), iran(:)
  real(kind=8),allocatable    :: nu_em(:), x_em(:,:), k_em(:,:)
  real(kind=8)                :: lambda_min, lambda_max, k_q(3), k_obs(3), halo_center(3), rvir, aperture
  character(1000)             :: repository, file

  !Parameters to choose
  nphoton = 500000
  totalflux = 1.2d52 !
  iseed = 476        !Random integer
  lambda_max = 1300 ; lambda_min = 1304
  nquasars = 5
  k_obs(:) = (/ 0.0d0, 0.0d0, -1.0d0 /)
  halo_center(:) = (/ 0.48030293d0, 0.49758723d0, 0.5158643d0 /)
  rvir = 0.00358d0
  repository = '/home/cral/mauerhofer/rascas/RascasFiles/simu_jeremy/00107/OI/quasar/k_obs_00-1/'
  aperture = 5d0  !in degrees

  allocate(ID(nphoton), iran(nphoton), nu_em(nphoton), x_em(3,nphoton), k_em(3,nphoton))

  ID(:) = (/ (i, i=1,nphoton) /)
  iran(:) = (/ (-i, i=1,nphoton) /)  !inspired from PhotonsFromStars.f90


  do j=6,6!1,nquasars

     nu_em(:) = (/ (clight*1d8/(lambda_max - (lambda_max-lambda_min)*ran3(iseed)), i=1,nphoton) /)
     
     do i=1,nphoton
        x_em(:,i) = halo_center(:) - 9d0/1.1d1*rvir*k_obs(:) + (j-1)*rvir/1d1*(/ 0d0, k_obs(3), -k_obs(2) /) 
        k_em(:,i) = k_obs(:)
     end do


     ! nu_em(:) = (/ (clight*1d8/(lambda_max - (lambda_max-lambda_min)*ran3(iseed)), i=1,nphoton) /)
     ! do
     !    call isotropic_direction(k_q, iseed)
     !    if(dot_product(k_q,-k_obs) > cos(aperture*3.14159d0/180d0)) exit
     ! end do

     ! print*, 'k_quasar  = ', k_q

     ! do i=1,nphoton
     !    x_em(:,i) = halo_center(:) + rvir*k_q(:)
     !    k_em(:,i) = k_obs(:)
     ! end do

     write(file, '(a,a,I2.2,a)') trim(repository),'/output_',j+18,'.out'

     open(unit=11, file=file, form='unformatted', action='write')
     write(11) nphoton
     write(11) totalflux
     write(11) iseed
     write(11) (ID(i), i=1,nphoton)
     write(11) (nu_em(i), i=1,nphoton)
     write(11) (x_em(:,i), i=1,nphoton)
     write(11) (k_em(:,i), i=1,nphoton)
     write(11) (iran(i), i=1,nphoton)
     close(11)

  end do

end program write_IC

  
