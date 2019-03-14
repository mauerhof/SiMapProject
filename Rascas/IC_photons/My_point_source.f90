!Program that writes initial conditions with photons coming from a central point source, with random directions
program write_IC

  use module_constants
  use module_random
  use module_utils

  integer(kind=4)             :: nphoton, iseed, i
  real(kind=8)                :: totalflux   !nb of real photons [# / s]
  integer(kind=4),allocatable :: ID(:), iran(:)
  real(kind=8),allocatable    :: nu_em(:), x_em(:,:), k_em(:,:)
  real(kind=8)                :: lambda_min, lambda_max

  !Parameters to choose
  nphoton = 1000000
  totalflux = 1.2d52 !
  iseed = 476        !Random integer
  lambda_max = 1300 ; lambda_min = 1308 
  
  
  allocate(ID(nphoton), iran(nphoton), nu_em(nphoton), x_em(3,nphoton), k_em(3,nphoton))

  ID(:) = (/ (i, i=1,nphoton) /)
  iran(:) = (/ (-i, i=1,nphoton) /)  !inspired from PhotonsFromStars.f90
  nu_em(:) = (/ (clight*1d8/(lambda_max - (lambda_max-lambda_min)*ran3(iseed)), i=1,nphoton) /)
  
  do i=1,nphoton
     x_em(:,i) = (/ 0.48030293d0, 0.49758723d0, 0.519d0 /)
     k_em(:,i) = (/ 0.0d0, 0.0d0, -1.0d0 /)
  end do

  open(unit=11, file='/Users/mauerhof/Documents/rascas/RascasFiles/simu_Jeremy/00107/OI/quasar/output_1.out', form='unformatted', action='write')
  write(11) nphoton
  write(11) totalflux
  write(11) iseed
  write(11) (ID(i), i=1,nphoton)
  write(11) (nu_em(i), i=1,nphoton)
  write(11) (x_em(:,i), i=1,nphoton)
  write(11) (k_em(:,i), i=1,nphoton)
  write(11) (iran(i), i=1,nphoton)
  close(11)
  

end program write_IC

  
