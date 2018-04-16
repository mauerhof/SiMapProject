! Additional stuff for outputting the SFR

module sfr_parameters
  use amr_parameters, ONLY: dp
  real(dp)::sfr_tot_old=0             	! Total SFR old
  real(dp)::sfr_tot_new=0		! Total SFR new

end module sfr_parameters

! Patched for outputting star formation rates
!################################################################
!################################################################
!################################################################
!################################################################
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  use random
  use sfr_parameters 		! see above!                                         
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled. 
  ! It modifies hydrodynamic variables according to mass conservation 
  ! and assumes an isothermal transformation... 
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::t0,d0,d00,e0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star,ndebris_tot
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  integer ::ntot,ntot_all,info,nstar_corrected,ideb,ndeb
  logical ::ok_free,ok_all
  real(dp)::d,x,y,z,u,v,w,e,zg,vdisp,dgas
  real(dp)::mstar,dstar,tstar,nISM,nCOM
  real(dp)::T2,nH,T_poly
  real(dp)::velc,uc,vc,wc,mass_load
  real(dp)::vxgauss,vygauss,vzgauss,birth_epoch
  real(kind=8)::mlost,mtot,mlost_all,mtot_all
  real(kind=8)::RandNum,GaussNum,PoissMean   
  real(dp)::vsn,costheta,sintheta,phi,cosphi,sinphi
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::mdebris,vdebris,zdebris,rdebris
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::list_debris,ind_debris
  logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  
  
  !!!!!! variables to add up SFRs in each cell !!!!!!
  real(dp)::sfr_local_old, sfr_local_new
    
  ! Parameters for the Kennicutt-Schidt relation
  real(dp)::gammaSH 		   ! ratio of specific heats
  real(dp)::f_gas	   	   ! mass fraction in gas f_g
  real(dp)::SchmidtLawCoefficient  ! =A, units implicit of M_dot pc^-2 yr^-1
  real(dp)::SchmidtLawExponent	   ! =n
  
  ! some constants
  real(dp)::Gcgs   	! in cm^3 g^{-1} s^{-2}
  real(dp)::MSUNcgs	! Solar mass in g
  real(dp)::PCcgs	! Parsec in cm
  real(dp)::YRcgs	! Year in seconds

  real(dp),dimension(1:nvector),save::gasDensity, P_tot, sfr
  real(dp),dimension(1:nvector),save::t_gas	! gas consumption time scale. previously t_star
  real(dp)::MSUNpPC2cgs, preFactor

  ! variables to calc the pressure of cells
  real(dp)::eDens, eKin
  integer:: i_dim

  sfr_local_new 	= 0.0
  sfr_local_old 	= 0.0
  !gammaSH		= 5.0/3.0
  gammaSH		= gamma
  f_gas 		= 1.0
  SchmidtLawCoefficient = 1.515e-4/1.e6		! 1.e6 fixes the messed up units from EAGLE 
  SchmidtLawExponent	= 1.4
  Gcgs   		= 6.6725985e-8 	   	! in cm^3 g^{-1} s^{-2}
  MSUNcgs		= 1.9891e33	   	! Solar mass in g
  PCcgs			= 3.08567758e18    	! Parsec in cm
  YRcgs			= 365.25*24.*3600. 	! Year in seconds


  !!!!!! things we only need to compute once !!!!!!
  ! convert SchmidtLawCoefficient to cgs
  SchmidtLawCoefficient = SchmidtLawCoefficient * MSUNcgs / (YRcgs * PCcgs**2) ! = A 
  
  ! convert (1 M_dot pc^-2) =   MSUNpPC2 to cgs
  MSUNpPC2cgs 	= MSUNcgs /( PCcgs**2)
  
  ! calculate prefactor (1 M_dot pc^-2)^n / A
  preFactor 	= ((MSUNpPC2cgs)**SchmidtLawExponent)/SchmidtLawCoefficient
  
  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering star_formation'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Star formation time scale from Gyr to code units
  ! SFR apply here for long lived stars only
  t0=t_star*(1d9*365.*24.*3600.)/scale_t

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
  nISM = MAX(nCOM,nISM)
  d0   = nISM/scale_nH
  d00  = n_star/scale_nH

  ! Initial star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  dstar=mstar/vol_loc

  ! Birth epoch
  birth_epoch=t

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

#if NDIM==3
  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)/d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e-0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do
     end do
  end do

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0
  ndebris_tot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Density criterion
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           uold(ind_cell(i),7) = 0.0
           if(d<=d0)ok(i)=.false. 
        end do
!        ! Temperature criterion                    <--- turned off for now, fk
!         do i=1,ngrid
!            T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0)
!            nH=uold(ind_cell(i),1)*scale_nH
!            T_poly=T2_star*(nH/nISM)**(g_star-1.0)
!            T2=T2-T_poly
!            if(T2>2e4)ok(i)=.false. 
!         end do
        ! Geometrical criterion
        if(ivar_refine>0)then
           do i=1,ngrid
              d=uold(ind_cell(i),ivar_refine)
              if(d<=var_cut_refine)ok(i)=.false.
           end do
        endif
        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i)=0
           if(ok(i))then
              ! Compute mean number of events
              d=uold(ind_cell(i),1)
              mcell=d*vol_loc
              tstar=t0*sqrt(d00/d)	
              
	      !! RAMSES-prescription
	      if((d/tstar)*dtnew(ilevel)  > 0.5*d)then
		sfr_local_old = sfr_local_old 		! SF density rate is zero
	      else
		sfr_local_old = sfr_local_old+d/tstar	! SF density rate (code units)              
              endif
              
              !! EAGLE-prescription -for now, this is done in cgs units
              gasDensity(i) = d*scale_d

	      ! pressure = energy density - e_kin density , see cooling_fine
	      eDens  	= 0.0
! 	      eKin	= 0.0
	      eDens	= uold(ind_cell(i),ndim+2)
! 	      not needed:	
! 	      do i_dim=1,ndim
! 		 eKin=eKin+0.5*uold(ind_cell(i),i_dim+1)**2/d		  
! 	      end do
! 	      write (*,'(A, 1pe12.5)') 'DEBUG: density = ', d
! 	      write (*,'(A, 1pe12.5)') 'DEBUG: eKin    = ', eKin
! 	      write (*,'(A, 1pe12.5)') 'DEBUG: eDens   = ', eDens
              ! FIX , CHECK in the morning
	      ! P_tot(i) = (eDens - eKin) * scale_d*scale_v**2

              P_tot(i) = (gamma-1.0)*(eDens) * scale_d*scale_v**2
              
	      P_tot(i) = max(P_tot(i),0d0)
	      
	      ! gas-consumption time-scale
              t_gas(i)= preFactor * ((gammaSH/Gcgs)*f_gas*P_tot(i))**(0.5*(1.-SchmidtLawExponent))
	          
              sfr(i)= gasDensity(i)/t_gas(i) 			! star formation density, conversion below!
              
              if(sfr(i)*dtnew(ilevel)*scale_t  > 0.5*d*scale_d)then
		sfr_local_new = sfr_local_new 		
              else
		sfr_local_new = sfr_local_new + sfr(i)
              endif
              
	      ! Here, we write the cell's SFR to passive scalar.
	      ! Passive scalars have to be multiplied by density d and
	      ! they should also be in code units, thus * (scale_d * scale_l^3 / scale_t)^-1
	      ! But since sfr(i) has units mass/(volume*time), we can skip the scale_l^3 term	       
	      uold(ind_cell(i),7) = sfr(i) * d * (scale_t/scale_d)	
 
	      ! Now, let us actually change the SFR recipe: we replace tstar with t_gas. 
	      ! Note that t_gas is given in seconds and tstar should be in code units
	      tstar = t_gas(i) / scale_t 
	      
	       
              PoissMean=dtnew(ilevel)/tstar*mcell/mstar
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas mass
              mgas=nstar(i)*mstar
              ! Security to prevent more than 50% of gas depletion
              if (mgas > 0.5*mcell) then
                 nstar_corrected=int(0.5*mcell/mstar)
                 mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
                 nstar(i)=nstar_corrected
              endif
              ! Compute new stars local statistics 
              mstar_tot=mstar_tot+nstar(i)*mstar
              if(nstar(i)>0)then
                 ntot=ntot+1
                 if(f_w>0)ndebris_tot=ndebris_tot+1
              endif
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i))=nstar(i)
        end do
     end do
  end do
  
  !!!!!!!!! Add up SFRs !!!!!!!!!

  sfr_tot_old = sfr_tot_old + sfr_local_old * vol_loc * dtnew(ilevel) 
  sfr_tot_new = sfr_tot_new + sfr_local_new * vol_loc * dtnew(ilevel) * scale_l**3 

  ! Now sfr_tot_* is the stellar mass formed during the coarse step in code units 
  if(ilevel .eq. levelmin) call output_sfr ! Once per coarse step
  
  !!!!!!!!! Done adding up SFRs !!!!!!!!!

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot-ndebris_tot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF4.1,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/mtot_all*100.
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do

        ! Gather new star arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Update linked list for stars
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Update linked list for debris
        if(f_w>0)then
           call remove_free(ind_debris,nnew)
           call add_list(ind_debris,ind_grid_new,ok_new,nnew)
        endif

        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star=index_star+1

           ! Get gas variables
           n=flag2(ind_cell_new(i))
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
           if(metal)zg=uold(ind_cell_new(i),imetal)

           ! Set star particle variables
           tp(ind_part(i))=birth_epoch  ! Birth epoch
           mp(ind_part(i))=n*mstar      ! Mass
           levelp(ind_part(i))=ilevel   ! Level
           idp(ind_part(i))=index_star  ! Star identity
           xp(ind_part(i),1)=x
           xp(ind_part(i),2)=y
           xp(ind_part(i),3)=z
           vp(ind_part(i),1)=u
           vp(ind_part(i),2)=v
           vp(ind_part(i),3)=w
           if(metal)zp(ind_part(i))=zg  ! Initial star metallicity

           ! Set GMC particle variables
           if(f_w>0)then
              ! Compute GMC mass without more than 50% of gas depletion
              mdebris=min(f_w*n*mstar,0.5*d*vol_loc-n*mstar)
              ! Add supernova ejecta
              mdebris=mdebris+eta_sn*n*mstar
              ! Remove ejecta from the long lived star mass
              mp(ind_part(i))=n*mstar-eta_sn*n*mstar
              ! Set GMC particle variables
              tp(ind_debris(i))=birth_epoch  ! Birth epoch
              mp(ind_debris(i))=mdebris      ! Mass
              levelp(ind_debris(i))=ilevel   ! Level
              idp(ind_debris(i))=-n          ! Number of individual stars
              xp(ind_debris(i),1)=x
              xp(ind_debris(i),2)=y
              xp(ind_debris(i),3)=z
              vp(ind_debris(i),1)=u
              vp(ind_debris(i),2)=v
              vp(ind_debris(i),3)=w
              ! GMC metallicity + yield from ejecta 
              if(metal)zp(ind_debris(i))=zg+eta_sn*yield*(1-zg)*n*mstar/mdebris
           endif

        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
           if(flag2(ind_cell(i))>0)then
              n=flag2(ind_cell(i))
              d=uold(ind_cell(i),1)
              uold(ind_cell(i),1)=max(d-n*dstar*(1.0+f_w),0.5*d)
           endif
        end do

     end do
     ! End loop over cells
  end do
  ! End loop over grids
  
  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=d*e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

#endif

end subroutine star_formation 
!################################################################
!################################################################
!################################################################
!################################################################





!*************************************************************************
subroutine output_sfr

! Print the average star formation rate in the simulation box during the
! coarse timestep [ Dubois+Teyssier02 (old) and Eagle (new) version]
!-------------------------------------------------------------------------
  use amr_commons
  use sfr_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::  info
  real(dp):: sfr_tot_old_all = 0.0 
  real(dp):: sfr_tot_new_all = 0.0 
  real(dp):: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),parameter:: yr2s = 3.1556926d7
  real(dp),parameter:: Msun = 1.98d33 ! Solar mass in grams
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(sfr_tot_old, sfr_tot_old_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info) 
  call MPI_ALLREDUCE(sfr_tot_new, sfr_tot_new_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)            
#endif

  ! sfr_tot_old_all is total mass formed during the coarse step in code units
  ! Convert to a star formation rate (still code units)

  sfr_tot_old_all = sfr_tot_old_all / dtnew(levelmin)
  sfr_tot_new_all = sfr_tot_new_all / dtnew(levelmin)

  ! Convert to g/s (sfr_tot_new_all is already in cgs)
  sfr_tot_old_all = sfr_tot_old_all * scale_d * scale_l**3 / scale_t
  
  ! Convert to Msun/yr
  sfr_tot_old_all = sfr_tot_old_all / Msun * yr2s
  sfr_tot_new_all = sfr_tot_new_all / Msun * yr2s
  
  if(myid .eq. 1)then 
       write(*, 113) t*scale_t / yr2s, sfr_tot_old_all
       write(*, 114) t*scale_t / yr2s, sfr_tot_new_all
  endif

  sfr_tot_old = 0d0
  sfr_tot_new = 0d0

113 format(1pe12.5, 1pe12.5, '  # SFR_RAMSES (t [yr] | SFR [Msun/yr])')
114 format(1pe12.5, 1pe12.5, '  # SFR_EAGLE  (t [yr] | SFR [Msun/yr])')
! 113 format('SFR_old:   t [yr] = ',1pe12.5, '   SFR [Msun/yr])= ',1pe12.5)
! 114 format('SFR_new:   t [yr] = ',1pe12.5, '   SFR [Msun/yr])=   ',1pe12.5)
end subroutine output_sfr

