program DMhalo_generator
  implicit none
  logical ::output_vel=.true.
  integer ::type_DM  =1      ! DM profile
  integer ::type_gas =1      ! gas profile
  integer ::type_star=0      ! star profile
  real(kind=8)::rmax=22.5d0  ! Maximum radius of the particles in r_s units
                             ! (characteristic scale of the profile)
!  real(kind=8)::rmax=1.733d0     ! Maximum radius of the particles in r_s units
!                             ! (characteristic scale of the profile)
  real(kind=8)::boxlen=60d0  ! Box length
  real(kind=8)::lsoft=0.8d0  ! Softening length for correction of particle
                             ! velocities
  real(kind=8)::v200=35d0    ! Virial velocity in km/s
  real(kind=8)::c=10d0       ! Concentration parameter
  real(kind=8)::h=0.7d0      ! Hubble constant
  real(kind=8)::fB=0.15d0    ! Baryon fraction
  real(kind=8),parameter::pi  = 3.14159265d0
  real(kind=8),parameter::rhoc= 1.88000000d-29
  integer,parameter::nbin=10000
  integer,parameter::npart=100000 ! Number of DM particles
  integer::ncum,i,ind
  real(kind=8),dimension(1:nbin)::rho_DM,rho_gas,rho_star,phig
  real(kind=8),dimension(1:npart,1:3)::xp,vp
  real(kind=8),dimension(1:npart)::mp
  real(kind=8)::r,x,RandNumber,P,q,a1,a2
  real(kind=8)::cosphi,sinphi,costeta,sinteta,phimin,phimax
  real(kind=8)::vel,vrms,ex,ey,ez,ux,uy,uz,utot,phir
  real(kind=8)::DM_density,gas_density,star_density,alpha
  real(kind=8)::r200,scale_d,scale_l,scale_t,scale_v,scale_m,Mtot,alphamax

  ! Initialize the seed of the random number generator
  !call init_random_seed()

  ! NFW scale radius
  r200 = v200*1d5/sqrt(2d0*2d0*pi*6.67d-8/3.*200.*rhoc*h**2)
  ! scale_d converts mass density from user units into g/cc
  scale_d = rhoc*h**2*200./3.0*c*c*c/(log(1d0+c)-c/(1d0+c))
  ! scale_t converts time from user units into seconds
  scale_t = 1.0/sqrt(6.67d-8*scale_d)
  ! scale_l converts distance from user units into cm
  scale_l = r200/c
  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t
  ! scale_m converts mass in user units into g
  scale_m = scale_d * scale_l**3d0

  if(type_DM==1)write(*,*)'Generating',npart,' positions for a NFW halo...'
  if(type_DM==2)write(*,*)'Generating',npart,' positions for an isothermal halo..'
  ! Generate DM random positions
  ncum=0
  
  alphamax=alpha(type_DM,nbin,rmax)
  do while(ncum .lt. npart)
     call random_number(RandNumber)
     r=RandNumber*rmax
     call random_number(RandNumber)
     x=RandNumber
     P=r**2*DM_density(r,type_DM,rmax)
     q=P/0.25d0
     !q=P/alphamax
     if(x<q) then
        call random_number(RandNumber)
        cosphi=cos(RandNumber*2d0*pi)
        sinphi=sin(RandNumber*2d0*pi)
        call random_number(RandNumber)
        costeta=2d0*RandNumber-1d0
        sinteta=sqrt(1d0-costeta**2)
        ncum=ncum+1
        xp(ncum,1)=r*sinteta*cosphi
        xp(ncum,2)=r*sinteta*sinphi
        xp(ncum,3)=r*costeta
     endif
  enddo
  
  if (output_vel)then

  write(*,*)'Computing gravitational potential...'
  call cmp_grav_potential(type_DM,type_gas,type_star,phig,nbin,rmax,c)

  write(*,*)'Computing particle velocities...'
  do i=1,npart
     
     ! Radius of the particle
     r=sqrt(xp(i,1)**2+xp(i,2)**2+xp(i,3)**2)
     ex=xp(i,1)/r
     ey=xp(i,2)/r
     ez=xp(i,3)/r
     
     ! Find the closest radial cell to that particle
     ind=nint(r*nbin/rmax)
     ! Compute the rms velocity by Vrms=sqrt(-Phi)
     vrms=v200*sqrt(-phig(ind))

     x=1d5
     do while(x .ge. 0.9d0*sqrt(2d0))
        call random_number(RandNumber)
        a1=RandNumber
        call random_number(RandNumber)
        a2=RandNumber
        ! Maxwellian distribution
        x=sqrt(-2d0*log(a1+1d-15))*cos(2d0*pi*a2)
        x=abs(x)
     enddo
     vel=x*vrms
     
     if (r .le. lsoft*rmax)then
        call random_number(RandNumber)
        ux=2d0*RandNumber-1d0
        call random_number(RandNumber)
        uy=2d0*RandNumber-1d0
        call random_number(RandNumber)
        uz=2d0*RandNumber-1d0
        utot=sqrt(ux**2+uy**2+uz**2)
        ! Renormalize such that the norm is equal to 1
        ux=ux/utot
        uy=uy/utot
        uz=uz/utot
     else
        phimin=pi/2d0*(r/rmax-lsoft)/(1d0-lsoft)
        phimax=pi-phimin
        do while(phir.lt.phimin .or. phir.gt.phimax)
           call random_number(RandNumber)
           ux=2d0*RandNumber-1d0
           call random_number(RandNumber)
           uy=2d0*RandNumber-1d0
           call random_number(RandNumber)
           uz=2d0*RandNumber-1d0
           utot=sqrt(ux**2+uy**2+uz**2)
           ! Renormalize such that the norm is equal to 1
           ux=ux/utot
           uy=uy/utot
           uz=uz/utot
           ! Compute the scalar product between velocity unit vector 
           ! and radial unit vector
           ! (simply gives the cosinus and the relative angle)
           phir=acos(ux*ex+uy*ey+uz*ez)
        enddo
     endif
     
     ! Finally assign its Maxwellian velocity
     vp(i,1)=ux*vel
     vp(i,2)=uy*vel
     vp(i,3)=uz*vel
  enddo

  else
     vp=0d0
  endif

  Mtot=4d0*pi*(log(1d0+rmax)-rmax/(1d0+rmax))*(1d0-fB)
  mp=Mtot/npart
  write(*,'(" Virial mass   =",1pe9.2," Msun")')4d0*pi*(log(1d0+c)-c/(1d0+c))*scale_d*scale_l**3d0/2d33
  write(*,'(" Total DM mass =",1pe9.2," Msun")')Mtot*scale_d*scale_l**3d0/2d33
  write(*,'(" mp            =",1pe9.2," Msun")')Mtot/npart*scale_d*scale_l**3d0/2d33
!  xp=xp+boxlen/2d0
!!$  do i=1,npart
!!$     xp(i,1)=xp(i,1)+boxlen/2d0
!!$     xp(i,2)=xp(i,2)+boxlen/2d0
!!$     xp(i,3)=xp(i,3)+boxlen/2d0
!!$  enddo
  vp=vp*1d5/scale_v
  open(unit=1,file='DM.dat')
  do i=1,npart
     write(1,'(7ES17.9)')xp(i,1),xp(i,2),xp(i,3),vp(i,1),vp(i,2),vp(i,3),mp(i)
  enddo
  close(1)
  write(*,*)'Particles are written in DM.dat'

end program DMhalo_generator
!################################################################
!################################################################
!################################################################
!################################################################
function alpha(type,nbin,rmax)
  implicit none
  integer::i,type,nbin
  real(kind=8)::alpha,rmax
  real(kind=8),dimension(1:nbin)::rho,r,P

  do i=1,nbin
     r(i)=i/nbin*rmax
     select case (type)
     case (1)
        ! NFW profile
        rho(i)=1/(r(i)*(1+r(i))**2)
     case (2)
        ! Isothermal profile
        rho(i)=1/r(i)**2
     end select
     P(i)=r(i)**2*rho(i)
  enddo
  alpha=maxval(P)

  return

end function alpha
!################################################################
!################################################################
!################################################################
!################################################################
function DM_density(r,type,rmax)
  implicit none
  integer::type
  real(kind=8)::r,rmax
  real(kind=8)::DM_density

  select case (type)
  case (1)
     ! NFW profile
     DM_density=1d0/(r*(1d0+r)**2)
  case (2)
     ! Isothermal profile
     DM_density=1d0/r**2
  end select
  return

end function DM_density
!################################################################
!################################################################
!################################################################
!################################################################
function gas_density(r,type,rmax)
  implicit none
  integer::type
  real(kind=8)::r,rmax
  real(kind=8)::gas_density

  select case (type)
  case (1)
     ! NFW profile
     gas_density=1d0/(r*(1d0+r)**2)
  case (2)
     ! Isothermal profile
     gas_density=1d0/r**2
  end select
  return

end function gas_density
!################################################################
!################################################################
!################################################################
!################################################################
function star_density(r,type,rmax)
  implicit none
  integer::type
  real(kind=8)::r,rmax
  real(kind=8)::star_density

  select case (type)
  case (0)
     star_density=0d0
  end select
  return

end function star_density
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cmp_grav_potential(type_DM,type_gas,type_star,phig,nbin,rmax,c)
  implicit none
  integer::type_DM,type_gas,type_star,i,nbin
  real(kind=8)::c,rmax
  real(kind=8)::f,r
  real(kind=8),dimension(1:nbin)::phig

  f=log(1d0+c)/c-1d0/(1d0+c)
  
  do i=1,nbin
     r=i*rmax/nbin
     phig(i)=-1d0/f * log(1d0+r)/r
  enddo

end subroutine cmp_grav_potential
!################################################################
!################################################################
!################################################################
!################################################################
