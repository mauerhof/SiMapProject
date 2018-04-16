module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(kind=8),allocatable,dimension(:)::msink,r2sink,v2sink,c2sink,oksink_new,oksink_all,tsink
  real(kind=8),allocatable,dimension(:)::msink_new,msink_all,r2k,v2sink_new,c2sink_new,tsink_new,tsink_all
  real(kind=8),allocatable,dimension(:)::v2sink_all,c2sink_all
  real(qdp)   ,allocatable,dimension(:)::xmsink
  real(kind=8),allocatable,dimension(:)::dMBHoverdt,dMEdoverdt,wdens,wvol,wc2,wdens_new,wvol_new,wc2_new,total_volume
  real(kind=8),allocatable,dimension(:,:)::wmom,wmom_new
  real(kind=8),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(kind=8),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(kind=8),allocatable,dimension(:,:)::weighted_density,weighted_volume,weighted_c2
  real(kind=8),allocatable,dimension(:,:)::jsink,jsink_new,jsink_all
  real(kind=8),allocatable,dimension(:)::dMBH_coarse,dMEd_coarse,dMsmbh,dMBH_coarse_new
  real(kind=8),allocatable,dimension(:)::dMEd_coarse_new,dMsmbh_new,dMBH_coarse_all,dMEd_coarse_all,dMsmbh_all
  real(kind=8),allocatable,dimension(:)::Esave,Esave_new,Esave_all
  real(kind=8),allocatable,dimension(:,:,:)::weighted_momentum
  real(kind=8),allocatable,dimension(:,:,:)::sink_stat,sink_stat_all

  integer ,allocatable,dimension(:)::idsink,idsink_new,idsink_all

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  real(dp),allocatable,dimension(:)  ::oxp      ! Birth O
  real(dp),allocatable,dimension(:)  ::fep      ! Birth Fe
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer ,allocatable,dimension(:)  ::idp      ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1
end module pm_commons
