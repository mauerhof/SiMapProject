!================================================================
!================================================================
!================================================================
!================================================================
subroutine boundana(x,u,dx,ibound,nn)
  use amr_commons
  use pm_commons
  use random
  use amr_parameters
  use hydro_parameters
  use poisson_parameters, ONLY: gravity_params
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xx,yy,zz,r,rr,v,xc,yc,zc,M,rho,dzz,zint,HH,rdisk,dpdr,dmax
  real(dp)::rgal,sum,sum2,dmin,zmin,zmax,c,fgas,pi,tol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,M_b,az,eps
  real(dp)::rmin,rmax,Tiso
  real(dp)::RandNum,phi,Rrand,SS,CC,UU
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  real(dp)::xl,yl,zl,xr,yr,zr,xit,yit,zit
  real(dp)::Axrr,Axlr,Axrl,Axll,Ayrr,Aylr,Ayrl,Ayll,Azrr,Azlr,Azrl,Azll
  real(dp)::dens,Bfloor,dxmin,Bzero
  integer ::it,nticks

  do ivar=1,nvar+3
     do i=1,nn
        u(i,ivar)=boundary_var(ibound,ivar)
     end do
  end do

  write(*,*)'Entering boundana'
  
  ! Add here, if you wish, some user-defined initial conditions
  xc  =boxlen/2.0
  yc  =boxlen/2.0
  zc  =boxlen/2.0

  c   =gravity_params(2)
  fgas=gravity_params(3)
  Tiso=gravity_params(4)
  Bzero=gravity_params(5)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rgal=3.5*3.08d21/scale_l
  HH  =0.2*3.08d21/scale_l
  rdisk=15.*3.08d21/scale_l
  !rdisk=150000.*3.08d21/scale_l
  pi=acos(-1.)

  ! Initialise utils for B field
  dxmin=boxlen*0.5d0**nlevelmax
  nticks=max(ceiling(dx/dxmin),1)
  Bfloor=Bzero/sqrt(4d0*pi)/(sqrt(scale_d)*scale_l/scale_t)

  dmax=fgas*(log(1.+c)-c/(1.+c))/(HH*rgal*rgal) &
       & /(1.-(1.+rdisk/rgal)*exp(-rdisk/rgal))
  !dmax=fgas*(log(1.+c)-c/(1.+c))/(HH*rgal*rgal)
  dmin=1d-3*fgas
  tol=1d-4
  !rmin=1./2.**nlevelmax/10.
  rmin=0.

  ! If necessary, initialize random number generator
  !if(localseed(1)==-1)then
  !   call rans(ncpu,iseed,allseed)
  !   localseed=allseed(myid,1:IRandNumSize)
  !end if

  do i=1,nn
     xx=x(i,1)-xc
     yy=x(i,2)-yc
     zz=x(i,3)-zc

     r=sqrt(xx**2+yy**2)
     ! Compute density
     !q(i,1)=dmin+dmax*exp(-r/rgal)*(2d0/(exp(zz/HH)+exp(-zz/HH)))**2
     q(i,1)=dmin+dmax*exp(-r/rgal)*exp(-(r/rdisk)**5)*(2.0/(exp(zz/HH)+exp(-zz/HH)))**2

     ! Compute Pressure
     ! Pressure is initialized such that Temperature is equal to Tiso everywhere
     q(i,ndim+2)=Tiso/scale_T2*q(i,1)

     !zmin=log(abs(zz))
     !zmax=log(sqrt((2.*boxlen)**2-r**2))
     !zmax=log(boxlen/2d0*1.1)
     !sum =romberg(zmin,zmax,tol)
     !r=r*(1d0+1d-2)
     !sum2=romberg(zmin,zmax,tol)
     !q(i,ndim+2)=sum
     ! Compute dP/dr
     !dpdr=(sum2-sum)/(1d-2*r)


     ! Compute angular velocity
     r  =sqrt(xx*xx+yy*yy)
     rr =sqrt(xx*xx+yy*yy+zz*zz)

     M=4d0*pi*(log(1d0+rr)-rr/(1d0+rr))*(1d0-fgas)
     if(r.le.rdisk)then
        rmax=rr
        sum=romberg(rmin,rmax,tol)
        M_b= 4d0*pi*dmax*HH*sum     + dmin*4./3.*pi*rr**3
     else
        M_b= 4d0*pi*dmax*HH*rgal**2 + dmin*4./3.*pi*rr**3
     endif
     !az=abs(r)
     !if(r.le.rdisk)then
     !   M_b=4d0*pi*dmax*HH*rgal**2*(1d0-(1d0+r/rgal)*exp(-r/rgal)) &
     !        &   *(exp(az/HH)-exp(-az/HH))/(exp(az/HH)+exp(-az/HH)) &
     !        &   +dmin*4./3.*pi*rr**3
     !else
     !   M_b=4d0*pi*dmax*HH*rgal**2*(1d0-(1d0+rdisk/rgal)*exp(-rdisk/rgal)) &
     !        &   *(exp(rdisk/HH)-exp(-rdisk/HH))/(exp(rdisk/HH)+exp(-rdisk/HH)) &
     !        &   +dmin*4./3.*pi*rr**3
     !endif
     !M_b=4d0*pi*dmax*HH*rgal**2d0*(1d0-(1d0+r/rgal)*exp(-r/rgal)) &
     !     &   *(exp(az/HH)-exp(-az/HH))/(exp(az/HH)+exp(-az/HH)) &
     !     &   +dmin*4d0/3d0*pi*rr**3d0
     M=M+M_b
     v= M*r**2/rr**3! + r*dpdr/q(i,1)
     if (v .gt. 0d0) then
        v=sqrt(v)
     else
        v=0d0
     endif
     !if(r.le.rdisk.and.abs(zz).le.3d0*HH)then 
     q(i,ndim-1)=-v*yy/r
     q(i,ndim  )= v*xx/r
     q(i,ndim+1)=0d0
        ! Random velocities
        !call ranf(localseed,RandNum)
        !SS =(RandNum-0.5)*2.
        !call ranf(localseed,RandNum)
        !phi=(RandNum-0.5)*2.*pi
        !call ranf(localseed,RandNum)
        !UU =RandNum
        !Rrand=0.1*v*UU**(1./3.)
        !CC=Rrand*sqrt(1.-SS**2.)
        !q(i,2)=q(i,2)+CC*cos(phi)
        !q(i,3)=q(i,3)+CC*sin(phi)
        !q(i,4)=q(i,4)+Rrand*SS
     !else
     !   q(i,ndim-1)=0.
     !   q(i,ndim  )=0.
     !   q(i,ndim+1)=0d0
     !endif

     ! -------------------------------
     ! Initialise B field
     ! -------------------------------     
     ! B left
     q(i,6)     =Bfloor
     q(i,7)     =0d0
     q(i,8)     =0d0
     ! B right
     q(i,nvar+1)=Bfloor
     q(i,nvar+2)=0d0
     q(i,nvar+3)=0d0

     ! Passive scalar
     if(metal)q(i,9)=0.

  enddo

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  ! Magnetic field components
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

contains
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fffy(z)
    implicit none
    !      Computes the integrand
    real(dp)::fffy
    real(dp)::z,rprim,rprim_max,tangenth
    !real(dp)::z,rm,zzz,M_b,M_d

    rprim=z
    !rmax_loc=rmax
    !rprim_max=sqrt(rmax_loc**2.-rprim**2.)
    rprim_max=sqrt(rmax**2.-rprim**2.)
    tangenth=(exp(rprim_max/HH)-exp(-rprim_max/HH))/(exp(rprim_max/HH)+exp(-rprim_max/HH))
    fffy=rprim*exp(-rprim/rgal)*tangenth
    
    !zzz=exp(z)
    !rm=sqrt(r**2+zzz**2)
    !rho=dmin+dmax*exp(-r/rgal)*(2.0/(exp(zzz/HH)+exp(-zzz/HH)))**2
    !rho=dmin+dmax*exp(-r/rgal)*exp(-(r/rdisk)**5)*(2.0/(exp(zzz/HH)+exp(-zzz/HH)))**2
    !M_d=4d0*pi*(log(1d0+rm)-rm/(1d0+rm))*(1.-fgas)
    !if(r.le.rdisk)then
    !   M_b=4d0*pi*dmax*HH*rgal**2*(1d0-(1d0+r/rgal)*exp(-r/rgal)) &
    !        &   *(exp(zzz/HH)-exp(-zzz/HH))/(exp(zzz/HH)+exp(-zzz/HH)) &
    !        &   +dmin*4./3.*pi*rm**3
    !else
    !   M_b=4d0*pi*dmax*HH*rgal**2*(1d0-(1d0+rdisk/rgal)*exp(-rdisk/rgal)) &
    !        &   *(exp(zzz/HH)-exp(-zzz/HH))/(exp(zzz/HH)+exp(-zzz/HH))
    !endif
    !M=M_d+M_b
    !fffy=M*rho*zzz*zzz*rm**(-1.5)
    return
  end function fffy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccc
  function romberg(a,b,tol)
    implicit none
    real(dp)::romberg
    !
    !     Romberg returns the integral from a to b of f(x)dx using Romberg 
    !     integration. The method converges provided that f(x) is continuous 
    !     in (a,b). The function f must be double precision and must be 
    !     declared external in the calling routine.  
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fffy(a)+fffy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !	Calculate next trapezoidal rule approximation to integral.
       
       g0=0.0d0
       do k=1,nint
          g0=g0+fffy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0
       
       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    romberg=g0

    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Romberg failed to converge; integral, error=', &
         &    romberg,error

    

    return
  end function romberg

end subroutine boundana
