!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine stellar_wind_kinetic
   use amr_commons
   use pm_commons
   use hydro_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
   integer::nCW_tot_mpi
   integer,dimension(1:ncpu)::nCW_icpu_mpi
   real(dp),dimension(:),allocatable::mCW_mpi,ZCW_mpi,OxCW_mpi,FeCW_mpi
   real(dp),dimension(:),allocatable::eCW_mpi,dtCW_mpi
   real(dp),dimension(:,:),allocatable::xCW_mpi,vCW_mpi
   integer,dimension(:),allocatable::icpuCW_mpi
#endif
   !-----------------------------------------------------------------
   ! This routine injects blast waves from stellar wind+SNIa 
   !   if the condition (r_Bcw > rmas) is met 
   ! Taysun Kimm 30/Apr/2011
   !             21/Aug/2011 (added Oxygen/Iron abundance)
   ! Abbreviations and descriptions of some variables
   !   CW: Continuous Winds, which includes stelalr winds and SNIa
   !   Ox: Oxygen
   !   Fe: Iron 
   !-----------------------------------------------------------------
   ! local scope
   integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part,ilevel
   integer::nCW,nCW_tot,info,iCW_ilevel,ivar,isort,ista,iend,nCW_myid,iCW,idpdum
   integer,dimension(1:ncpu)::nCW_icpu
   logical::ok_free
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0,mass_load
   integer,dimension(:),allocatable::ind_CW,iCW_myid
   real(dp),dimension(:),allocatable::mCW,ZCW,vol_gas,vol_rho,ekBlast,mloadCW
   real(dp),dimension(:),allocatable::ZloadCW,OxloadCW,FeloadCW
   real(dp),dimension(:),allocatable::dtCW,eCW,OxCW,FeCW
   real(dp),dimension(:,:),allocatable::xCW,vCW,dq,vloadCW,xCW_tot,vCW_tot
   integer,dimension(:),allocatable::ind_CWtot,itemp,idpCW,icpuCW,icpuCW_tot
   real(dp),dimension(:),allocatable::mCW_tot,ZCW_tot,dtCW_tot,eCW_tot,OxCW_tot,FeCW_tot
   real(dp),dimension(:),allocatable::mpinit,mCWind,eCWind,ZCWind,dtlast,OxCWind,FeCWind
   real(dp),dimension(:),allocatable::zpsub
   real(dp)::scale_msun,ttsta,ttend
   integer::nx_loc

   if(.not.hydro) return
   if(ndim.ne.3)  return
   if(myid.eq.1) write(*,*)'Entering stellar_wind_kinetic'

!!$   if (myid.eq.1) ttsta=MPI_WTIME(info)

   ! Conversion factor from user units to cgs units
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   t0=t_ctw*(365d0*24d0*3600d0)/scale_t
   scale_msun=(scale_d*scale_l**3d0)/2d33

   ! Count the number of star particles that can *potentially* undergo the blast phase
   nCW_tot=0
   do icpu=1,ncpu
   ! Loop over cpus
       igrid=headl(icpu,levelmin)
       ! Loop over grids
       do jgrid=1,numbl(icpu,levelmin)
          npart1=numbp(igrid)  ! Number of particles in the grid
          npart2=0
          ! Count star particles younger than t_ctw
          if(npart1>0)then
             ipart=headp(igrid)
             ! Loop over particles
             do jpart=1,npart1
                ! Save next particle   <--- Very important !!!
                next_part=nextp(ipart)
                if( (tp(ipart).ne.0).and.(tp(ipart)>=(t-t0)) ) then
                  npart2=npart2+1
                endif
                ipart=next_part  ! Go to next particle
             end do
          endif
  
          nCW_tot=nCW_tot+npart2   ! Add the number of particles to the total
          igrid=next(igrid)   ! Go to next grid
      end do
   enddo

   nCW_icpu=0
   nCW_icpu(myid)=nCW_tot
   
#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(nCW_icpu,nCW_icpu_mpi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
   nCW_icpu=nCW_icpu_mpi
#endif
   nCW_tot=sum(nCW_icpu(1:ncpu))

   if (nCW_tot.eq.0) return

   if (myid.eq.1) then
   write(*,*) ' -------------------------------------------'
   write(*,*) ' Multi-blast candidate (all cpus):', nCW_tot
   write(*,*) ' -------------------------------------------'
   endif
   ! Allocate arrays for the position and the mass of the CW
   allocate(xCW_tot(1:nCW_tot,1:3),vCW_tot(1:nCW_tot,1:3),mCW_tot(1:nCW_tot),itemp(1:nCW_tot))
   allocate(ZCW_tot(1:nCW_tot),OxCW_tot(1:nCW_tot),FeCW_tot(1:nCW_tot))
   allocate(dtCW_tot(1:nCW_tot),eCW_tot(1:nCW_tot),icpuCW_tot(1:nCW_tot))
   allocate(idpCW(1:nCW_icpu(myid)))


   xCW_tot=0d0;vCW_tot=0d0;mCW_tot=0d0;ZCW_tot=0d0;OxCW_tot=0d0;FeCW_tot=0d0
   dtCW_tot=0d0;eCW_tot=0d0;icpuCW_tot=0

   ! Synchronise position and mass of the star to CW array
   if (myid==1)then
     iCW=0
   else
     iCW=sum(nCW_icpu(1:myid-1))
   endif
   ista=iCW+1


   do icpu=1,ncpu
   ! Loop over cpus
       igrid=headl(icpu,levelmin)
       ! Loop over grids
       do jgrid=1,numbl(icpu,levelmin)
          npart1=numbp(igrid)  ! Number of particles in the grid
          npart2=0
          ! Count star particles younger than t_ctw
          if(npart1>0)then
             ipart=headp(igrid)
             ! Loop over particles
             do jpart=1,npart1
                ! Save next particle   <--- Very important !!!
                next_part=nextp(ipart)
                if( (tp(ipart).ne.0).and.(tp(ipart).ge.(t-t0)) ) then
                   iCW=iCW+1
                   xCW_tot(iCW,1)=xp(ipart,1) 
                   xCW_tot(iCW,2)=xp(ipart,2) 
                   xCW_tot(iCW,3)=xp(ipart,3) 
                   vCW_tot(iCW,1)=vp(ipart,1)
                   vCW_tot(iCW,2)=vp(ipart,2)
                   vCW_tot(iCW,3)=vp(ipart,3)
                   idpCW(iCW-ista+1)=ipart
                endif
                ipart=next_part  ! Go to next particle
             end do
          endif
          igrid=next(igrid)   ! Go to next grid
       end do
   enddo

   iend=iCW
   nCW_myid=nCW_icpu(myid)

   if (nCW_myid.gt.0) then

     allocate(mpinit(1:nCW_myid),mCWind(1:nCW_myid),eCWind(1:nCW_myid),dtlast(1:nCW_myid))
     allocate(ZCWind(1:nCW_myid),OxCWind(1:nCW_myid),FeCWind(1:nCW_myid))
     allocate(zpsub(1:nCW_myid))

     mpinit=0d0;mcWind=0d0;eCwind=0d0;dtlast=0d0
     ZCWind=0d0;OxCWind=0d0;FeCWind=0d0
     zpsub=0d0

     if(metal)then
        zpsub(1:nCW_myid)=zp(idpCW(1:nCW_myid))
     endif


     ! Compute the Continuous Wind properties (E,mdot,Z)
     call Mstar_Init(nCW_myid,mp(idpCW),tp(idpCW),zpsub,mpinit)

     call Wind_Prop(nCW_myid,mp(idpCW)/mpinit,tp(idpCW),zpsub,idp(idpCW),&
                   &mCWind,eCWind,ZCWind,OxCWind,FeCWind,dtlast,.false.)

     do iCW=ista,iend
        mCW_tot (iCW)=mCWind(iCW-ista+1)*mpinit(iCW-ista+1)  ! mCwind: fraction
        ZCW_tot (iCW)=ZCWind(iCW-ista+1)                     ! stellar yield (Z)
        OxCW_tot(iCW)=OxCWind(iCW-ista+1)                    ! stellar yield (O)
        FeCW_tot(iCW)=FeCWind(iCW-ista+1)                    ! stellar yield (Fe)
        dtCW_tot(iCW)=dtlast(iCW-ista+1)                     ! tstar-tlast
        eCW_tot (iCW)=eCwind(iCW-ista+1)*(mpinit(iCW-ista+1)*scale_msun) ! unit: [10^40 erg]
               ! Notice: this already takes into account mass loss. don't need to multiply dmloss
        icpuCW_tot(iCW)=myid
     enddo
     deallocate(mpinit,mCwind,eCwind,dtlast)
     deallocate(ZCWind,OxCWind,FeCWind)
     deallocate(zpsub)
   endif


#ifndef WITHOUTMPI
   allocate(xCW_mpi(1:nCW_tot,1:3),vCW_mpi(1:nCW_tot,1:3),mCW_mpi(1:nCW_tot))
   allocate(ZCW_mpi(1:nCW_tot),OxCW_mpi(1:nCW_tot),FeCW_mpi(1:nCW_tot))
   allocate(dtCW_mpi(1:nCW_tot),eCW_mpi(1:nCW_tot),icpuCW_mpi(1:nCW_tot))
   call MPI_ALLREDUCE(xCW_tot,xCW_mpi,nCW_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(vCW_tot,vCW_mpi,nCW_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(mCW_tot,mCW_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(ZCW_tot,ZCW_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(OxCW_tot,OxCW_mpi,nCW_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(FeCW_tot,FeCW_mpi,nCW_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(dtCW_tot,dtCW_mpi,nCW_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(eCW_tot,eCW_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(icpuCW_tot,icpuCW_mpi,nCW_tot,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
   xCW_tot=xCW_mpi
   vCW_tot=vCW_mpi
   mCW_tot=mCW_mpi
   ZCW_tot=ZCW_mpi
   OxCW_tot=OxCW_mpi
   FeCW_tot=FeCW_mpi
   dtCW_tot=dtCW_mpi
   eCW_tot=eCW_mpi
   icpuCW_tot=icpuCW_mpi
   deallocate(xCW_mpi,vCW_mpi,mCW_mpi)
   deallocate(ZCW_mpi,OxCW_mpi,FeCW_mpi)
   deallocate(dtCW_mpi,eCW_mpi,icpuCW_mpi)
#endif

   call getSNonmyid(itemp,nCW,xCW_tot,nCW_tot)
   !   nCW_icpu(myid)!=nCW if blast(s) from a neighbouring CPU affect the gas properties of cells

   ! Allocate the arrays for the position and the mass of the CW
   allocate(xCW(1:nCW,1:3),vCW(1:nCW,1:3),mCW(1:nCW),iCW_myid(1:nCW))
   allocate(ZCW(1:nCW),OxCW(1:nCW),FeCW(1:nCW))
   allocate(dtCW(1:nCW),eCW(1:nCW),icpuCW(1:nCW))

   xCW=0d0;  vCW=0d0;  mCW=0d0; iCW_myid=0; dtCW=0d0;eCW=0d0
   ZCW=0d0; OxCW=0d0; FeCW=0d0 

   do iCW=1,nCW
      isort=itemp(iCW)
      iCW_myid(iCW)=isort
      xCW(iCW,1)=xCW_tot(isort,1)
      xCW(iCW,2)=xCW_tot(isort,2)
      xCW(iCW,3)=xCW_tot(isort,3)
      vCW(iCW,1)=vCW_tot(isort,1)
      vCW(iCW,2)=vCW_tot(isort,2)
      vCW(iCW,3)=vCW_tot(isort,3)
      mCW(iCW)  =mCW_tot(isort)
      ZCW(iCW)  =ZCW_tot(isort)
      OxCW(iCW) =OxCW_tot(isort)
      FeCW(iCW) =FeCW_tot(isort)
      dtCW(iCW) =dtCW_tot(isort)
      eCW(iCW)  =eCW_tot(isort)
      icpuCW(iCW)=icpuCW_tot(isort)
   enddo

   deallocate(xCW_tot,vCW_tot,mCW_tot,itemp)
   deallocate(ZCW_tot,OxCW_tot,FeCW_tot)
   deallocate(dtCW_tot,eCW_tot,icpuCW_tot)

   allocate(vol_gas(1:nCW),vol_rho(1:nCW),dq(1:nCW,1:3),ekBlast(1:nCW),ind_CW(1:nCW))
   allocate(mloadCW(1:nCW),ZloadCW(1:nCW),OxloadCW(1:nCW),FeloadCW(1:nCW),vloadCW(1:nCW,1:3))


   ! Compute the grid discretization effects 
   !    AND check whether shock fronts by stellar winds (or SNIa) has propagated distant enough to launch a blast wave 
   call average_CW(xCW,vCW,mCW,ZCW,OxCW,FeCW,eCW,dtCW,vol_gas,vol_rho,dq,ekBlast,ind_CW&
                     &,nCW,nCW_tot,iCW_myid,mloadCW,ZloadCW,OxloadCW,FeloadCW,vloadCW)

!   if (debug) then 
!     do iCW=1,nCW
!       write(800+myid,*) sngl(vol_gas(iCW)),sngl(mCW(iCW)), sngl(mloadCW(iCW))
!     enddo
!     call MPI_BARRIER(MPI_COMM_WORLD,info)
!     call clean_stop
!   endif

   ! Modify hydro quantities to account for a Sedov blast wave
   call Sedov_blast_CW(xCW,mCW,eCW,ind_CW,vol_gas,vol_rho,dq,ekBlast,nCW,mloadCW,ZloadCW,OxloadCW,FeloadCW,vloadCW)

   ! Update stellar mass
   do iCW=1,nCW
      if (icpuCW(iCW).eq.myid.and.vol_gas(iCW)>=0d0) then  ! to avoid subtracting mass twice
        isort=iCW_myid(iCW)
        idpdum=idpCW(isort-ista+1)
        mp(idpdum)=mp(idpdum)-mCW(iCW)
      endif
   enddo 


   deallocate(xCW,vCW,mCW,ZCW,OxCW,FeCW,iCW_myid,dtCW,eCW)
   deallocate(vol_gas,vol_rho,dq,ekBlast,ind_CW,mloadCW,ZloadCW,OxloadCW,FeloadCW,vloadCW)
   deallocate(idpCW,icpuCW)

!!$   if (myid.eq.1) then
!!$      ttend=MPI_WTIME(info)
!!$      write(*,*) ' Time elapsed in stellar_wind [sec]:', sngl(ttend-ttsta)
!!$   endif



   ! Update hydro quantities for split cells
   do ilevel=nlevelmax,levelmin,-1
      call upload_fine(ilevel)
      do ivar=1,nvar
         call make_virtual_fine_dp(uold(1,ivar),ilevel)
      enddo 
   enddo


end subroutine stellar_wind_kinetic
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine average_CW(xCW,vCW,mCW,ZCW,OxCW,FeCW,eCW,dtCW,vol_gas,vol_rho,dq,ekBlast,ind_blast&
                     &,nCW,nCW_tot,iCW_myid,mloadCW,ZloadCW,OxloadCW,FeloadCW,vloadCW)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  ! and do the mass loading process.
  ! In addition, this estimates the size of r_Bcw (Blast of stellar wind) and
  ! return physical properties if r_Bcw > dx_loc 
  ! r_Bcw = (E0/rho0)^1/5 t^1/5 
  ! Taysun Kimm
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nCW,nCW_tot,j,iCW,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_CW,d,u,v,w,ek,u2,v2,w2,dr_cell,dr_cell_tmp
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_g
  real(dp)::eint,ekk,ekk1,ekk2,heat,mload,Zload,Oxload,Feload
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nCW)::ind_blast
  real(dp),dimension(1:nCW)::vol_gas,ekBlast,vol_loc_Blast,vol_rho
  real(dp),dimension(1:nCW)::mCW,ZCW,OxCW,FeCW,eCW,dtCW,m_gas
  real(dp),dimension(1:nCW)::mloadCW,ZloadCW,OxloadCW,FeloadCW
  real(dp),dimension(1:nCW,1:3)::xCW,vCW,dq,u2Blast,vloadCW
#ifndef WITHOUTMPI
  real(dp),dimension(1:nCW_tot)::vol_gas_mpi,m_gas_mpi,vol_rho_mpi
  real(dp),dimension(1:nCW_tot)::vol_gas_tot,m_gas_tot,vol_rho_tot
  real(dp),dimension(1:nCW_tot)::mloadCW_mpi,mloadCW_tot,ZloadCW_mpi,ZloadCW_tot
  real(dp),dimension(1:nCW_tot)::OxloadCW_mpi,OxloadCW_tot,FeloadCW_mpi,FeloadCW_tot
  real(dp),dimension(1:nCW_tot,1:3)::dq_mpi,u2Blast_mpi,vloadCW_mpi
  real(dp),dimension(1:nCW_tot,1:3)::dq_tot,u2Blast_tot,vloadCW_tot
#endif
  logical ,dimension(1:nvector),save::ok
  logical ,dimension(1:nCW)::okBlast
  integer ,dimension(1:nCW)::iCW_myid
  integer::ind_CW
  real(dp)::r_Bcw,rho0_CW,r_windmax,rlocmax,yr2s
  real(dp)::pxmi,pxma,pymi,pyma,pzmi,pzma
  if(verbose)write(*,*)'Entering average_CW'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  yr2s=365d0*24d0*3600d0
  scale_g=scale_d*scale_l**3d0

  ! Maximum radius of the ejecta
  rmax=MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax


  ! To reduce time spent to search for cells
  pxmi=minval(xCW(:,1))
  pxma=maxval(xCW(:,1))
  pymi=minval(xCW(:,2))
  pyma=maxval(xCW(:,2))
  pzmi=minval(xCW(:,3))
  pzma=maxval(xCW(:,3))


  ! Initialize the averaged variables
  vol_gas=0d0;dq=0d0;u2Blast=0d0;ekBlast=0d0;ind_blast=-1;m_gas=0d0;vol_loc_Blast=0d0
  mloadCW=0d0;vloadCW=0d0;ZloadCW=0d0;OxloadCW=0d0;FeloadCW=0d0

  vol_rho=0d0
  dr_cell_tmp=1d0

  !----------------------------
  ! Compute dq(:,1:3)
  !----------------------------
  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim

     ! compute the rough estimate of the region where stars are located 
     rlocmax=MAX(dx_loc/2d0,rmax)

     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rlocmax).or.(x.gt.pxma+rlocmax).or.&
                   &(y.lt.pymi-rlocmax).or.(y.gt.pyma+rlocmax).or.&
                   &(z.lt.pzmi-rlocmax).or.(z.gt.pzma+rlocmax)) then
                    ok(i)=.false.
                 endif
              endif
          enddo

          do i=1,ngrid
             if(ok(i))then
                ! Get gas cell position
                x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                do iCW=1,nCW
                   ! Check if the cell lies within the CW radius
                   dxx=x-xCW(iCW,1)
                   dyy=y-xCW(iCW,2)
                   dzz=z-xCW(iCW,3)
                   dr_CW=dxx**2+dyy**2+dzz**2
                   dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                   if(dr_CW.lt.rmax2)then
                      vol_gas(iCW)=vol_gas(iCW)+vol_loc
                      m_gas(iCW)=m_gas(iCW)+vol_loc*uold(ind_cell(i),1)

                      ! Take account for grid effects on the conservation of the
                      ! normalized linear momentum
                      u=dxx/rmax
                      v=dyy/rmax
                      w=dzz/rmax

                      ! Add the local normalized linear momentum to the total linear
                      ! momentum of the blast wave (should be zero with no grid effect)
                      if(pseudo_sedov)then
                         dr_cell_tmp=dsqrt(u*u + v*v + w*w)
                         vol_rho(iCW)=vol_rho(iCW)+vol_loc*dr_cell_tmp
                      endif

                      dq(iCW,1)=dq(iCW,1)+u*vol_loc*dr_cell_tmp
                      dq(iCW,2)=dq(iCW,2)+v*vol_loc*dr_cell_tmp
                      dq(iCW,3)=dq(iCW,3)+w*vol_loc*dr_cell_tmp
                   endif

                end do
             endif
          end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels


#ifndef WITHOUTMPI
  vol_gas_tot=0d0; dq_tot=0d0; m_gas_tot=0d0; vol_rho_tot=0d0
  ! Put the nCW size arrays into nCW_tot size arrays to synchronize processors
  do iCW=1,nCW
     ind_CW=iCW_myid(iCW)
     vol_gas_tot(ind_CW)=vol_gas(iCW)
     vol_rho_tot(ind_CW)=vol_rho(iCW)
     dq_tot     (ind_CW,1)=dq     (iCW,1)
     dq_tot     (ind_CW,2)=dq     (iCW,2)
     dq_tot     (ind_CW,3)=dq     (iCW,3)
     m_gas_tot(ind_CW)=m_gas(iCW)
  enddo
  call MPI_ALLREDUCE(vol_gas_tot,vol_gas_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vol_rho_tot,vol_rho_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq_tot     ,dq_mpi     ,nCW_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m_gas_tot,m_gas_mpi,nCW_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_tot=vol_gas_mpi
  vol_rho_tot=vol_rho_mpi
  dq_tot     =dq_mpi
  m_gas_tot  =m_gas_mpi

  ! Put the nCW_tot size arrays into nCW size arrays
  do iCW=1,nCW
     ind_CW=iCW_myid(iCW)
     vol_gas(iCW)  =vol_gas_tot(ind_CW)
     vol_rho(iCW)  =vol_rho_tot(ind_CW)
     dq     (iCW,1)=dq_tot     (ind_CW,1)
     dq     (iCW,2)=dq_tot     (ind_CW,2)
     dq     (iCW,3)=dq_tot     (ind_CW,3)
     m_gas(iCW)=m_gas_tot(ind_CW)
  enddo
#endif


  if (pseudo_sedov)then
     do iCW=1,nCW
     if(vol_gas(iCW)>0d0.and.mCW(iCW).gt.0d0)then
        dq(iCW,1)=dq(iCW,1)/vol_rho(iCW)
        dq(iCW,2)=dq(iCW,2)/vol_rho(iCW)
        dq(iCW,3)=dq(iCW,3)/vol_rho(iCW)
     endif
     enddo
  else
     do iCW=1,nCW
     if(vol_gas(iCW)>0d0.and.mCW(iCW).gt.0d0)then
        dq(iCW,1)=dq(iCW,1)/vol_gas(iCW)
        dq(iCW,2)=dq(iCW,2)/vol_gas(iCW)
        dq(iCW,3)=dq(iCW,3)/vol_gas(iCW)
     endif
     enddo
  endif



  !----------------------------
  ! Compute ekBlast(:,1:3)
  !----------------------------
  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim

     ! compute the rough estimate of the region where stars are located 
     rlocmax=MAX(dx_loc/2d0,rmax)

     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rlocmax).or.(x.gt.pxma+rlocmax).or.&
                   &(y.lt.pymi-rlocmax).or.(y.gt.pyma+rlocmax).or.&
                   &(z.lt.pzmi-rlocmax).or.(z.gt.pzma+rlocmax)) then
                    ok(i)=.false.
                 endif
              endif
          enddo

          do i=1,ngrid
             if(ok(i))then
                ! Get gas cell position
                x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                do iCW=1,nCW
                   ! Check if the cell lies within the CW radius
                   dxx=x-xCW(iCW,1)
                   dyy=y-xCW(iCW,2)
                   dzz=z-xCW(iCW,3)
                   dr_CW=dxx**2+dyy**2+dzz**2
                   dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                   if(dr_CW.lt.rmax2)then
                      ! Take account for grid effects on the conservation of the
                      ! normalized linear momentum (now use updated value!!)
                      u=dxx/rmax - dq(iCW,1)
                      v=dyy/rmax - dq(iCW,2)
                      w=dzz/rmax - dq(iCW,3)

                      ! Add the local normalized linear momentum to the total linear
                      ! momentum of the blast wave (should be zero with no grid effect)
                      if(pseudo_sedov)then
                         dr_cell_tmp=dsqrt(u*u + v*v + w*w)
                      endif

                      u2Blast(iCW,1)=u2Blast(iCW,1)+u*u*vol_loc*dr_cell_tmp
                      u2Blast(iCW,2)=u2Blast(iCW,2)+v*v*vol_loc*dr_cell_tmp
                      u2Blast(iCW,3)=u2Blast(iCW,3)+w*w*vol_loc*dr_cell_tmp
                   endif

                   if(dr_cell.le.dx_loc/2.0)then
                      ind_blast(iCW)=ind_cell(i)
                      ekBlast  (iCW)=vol_loc
                      vol_loc_Blast(iCW)=vol_loc
                   endif
                end do
             endif
          end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels


#ifndef WITHOUTMPI
  u2Blast_tot=0d0
  ! Put the nCW size arrays into nCW_tot size arrays to synchronize processors
  do iCW=1,nCW
     ind_CW=iCW_myid(iCW)
     u2Blast_tot(ind_CW,1)=u2Blast(iCW,1)
     u2Blast_tot(ind_CW,2)=u2Blast(iCW,2)
     u2Blast_tot(ind_CW,3)=u2Blast(iCW,3)
  enddo
  call MPI_ALLREDUCE(u2Blast_tot,u2Blast_mpi,nCW_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  u2Blast_tot=u2Blast_mpi

  ! Put the nCW_tot size arrays into nCW size arrays
  do iCW=1,nCW
     ind_CW=iCW_myid(iCW)
     u2Blast(iCW,1)=u2Blast_tot(ind_CW,1)
     u2Blast(iCW,2)=u2Blast_tot(ind_CW,2)
     u2Blast(iCW,3)=u2Blast_tot(ind_CW,3)
  enddo
#endif

  if (pseudo_sedov)then
     do iCW=1,nCW
     if(vol_gas(iCW)>0d0.and.mCW(iCW).gt.0d0)then
        u2Blast(iCW,1)=u2Blast(iCW,1)/vol_rho(iCW)
        u2Blast(iCW,2)=u2Blast(iCW,2)/vol_rho(iCW)
        u2Blast(iCW,3)=u2Blast(iCW,3)/vol_rho(iCW)
     endif
     enddo
  else
     do iCW=1,nCW
     if(vol_gas(iCW)>0d0.and.mCW(iCW).gt.0d0)then
        u2Blast(iCW,1)=u2Blast(iCW,1)/vol_gas(iCW)
        u2Blast(iCW,2)=u2Blast(iCW,2)/vol_gas(iCW)
        u2Blast(iCW,3)=u2Blast(iCW,3)/vol_gas(iCW)
     endif
     enddo
  endif



  okBlast=.false.
  do iCW=1,nCW
     rho0_CW=0d0
     if(vol_gas(iCW)>0d0.and.mCW(iCW).gt.0d0)then
        rho0_CW=dble(m_gas(iCW)/vol_gas(iCW))*dble(scale_d)
        r_Bcw=(dble(eCW(iCW))/rho0_CW)**0.2d0*(dtCW(iCW)*yr2s)**0.4d0*(1d8/scale_l)
                                     ! 1d8=1d40^0.2 (unit of eCW)
        r_windmax=dble(dtCW(iCW)*yr2s)*dsqrt(2d0*eCW(iCW)/(mCW(iCW)*scale_g))*(1d20/scale_l)
                                     ! 1d20=1d40^0.5
        r_Bcw=min(r_Bcw,r_windmax)

        if (r_Bcw>=rmax) then 
          u2=u2Blast(iCW,1)
          v2=u2Blast(iCW,2)
          w2=u2Blast(iCW,3)
          ekBlast(iCW)=(u2+v2+w2)/2d0
           
          if(ekBlast(iCW).gt.0.01) then  ! for numerical stability
             okBlast(iCW)=.true.
          else
             vol_gas(iCW)=-1
             mloadCW(iCW)=0
          endif
        else
          vol_gas(iCW)=-1
          mloadCW(iCW)=0
        endif
     else
        vol_gas(iCW)=-1
        mloadCW(iCW)=0
     endif

  enddo 

  do iCW=1, nCW
    if (okBlast(iCW).and.ind_blast(iCW)>0) then 
             ! Note that one star can have a blast across CPUs
       vol_loc=vol_loc_Blast(iCW)
       d=uold(ind_blast(iCW),1)
       u=uold(ind_blast(iCW),2)/d
       v=uold(ind_blast(iCW),3)/d
       w=uold(ind_blast(iCW),4)/d
       ekk=0.5d0*d*(u*u+v*v+w*w)
       eint=uold(ind_blast(iCW),5)-ekk
       ! Mass loading factor of the Sedov explosion
       ! Ensure that no more that 25% of the gas content is removed
       mload=min(f_w*mCW(iCW),0.25d0*d*vol_loc)
       mloadCW(iCW)=mCW(iCW)+mload
       ! Update gas mass and metal content in the cell
       if(metal)then
          Zload=uold(ind_blast(iCW),imetal)/d
          ZloadCW(iCW)=( mload*Zload + ZCW(iCW)*mCW(iCW)) / mloadCW(iCW)  ! metallicity, NOT yield
          uold(ind_blast(iCW),imetal)=uold(ind_blast(iCW),imetal)-Zload*mload/vol_loc
          if(chemo)then
             Oxload=uold(ind_blast(iCW),imetal+1)/d
             OxloadCW(iCW)=(mload*Oxload + OxCW(iCW)*mCW(iCW)) / mloadCW(iCW)
             uold(ind_blast(iCW),imetal+1)=uold(ind_blast(iCW),imetal+1)-Oxload*mload/vol_loc

             Feload=uold(ind_blast(iCW),imetal+2)/d
             FeloadCW(iCW)=(mload*Feload + FeCW(iCW)*mCW(iCW)) / mloadCW(iCW)
             uold(ind_blast(iCW),imetal+2)=uold(ind_blast(iCW),imetal+2)-Feload*mload/vol_loc
          endif
       endif
       
       d=uold(ind_blast(iCW),1)-mload/vol_loc
       uold(ind_blast(iCW),1)=d
       uold(ind_blast(iCW),2)=d*u
       uold(ind_blast(iCW),3)=d*v
       uold(ind_blast(iCW),4)=d*w
       uold(ind_blast(iCW),5)=eint+0.5d0*d*(u*u+v*v+w*w)
   
       vloadCW(iCW,1)=(mCW(iCW)*vCW(iCW,1)+mload*u)/mloadCW(iCW)
       vloadCW(iCW,2)=(mCW(iCW)*vCW(iCW,2)+mload*v)/mloadCW(iCW)
       vloadCW(iCW,3)=(mCW(iCW)*vCW(iCW,3)+mload*w)/mloadCW(iCW)
    endif
  enddo


#ifndef WITHOUTMPI
mloadCW_tot=0d0; ZloadCW_tot=0d0; OxloadCW_tot=0d0; FeloadCW_tot=0d0;vloadCW_tot=0d0

  do iCW=1,nCW
     ind_CW=iCW_myid(iCW)
     mloadCW_tot(ind_CW)=mloadCW(iCW)
     ZloadCW_tot(ind_CW)=ZloadCW(iCW)
     OxloadCW_tot(ind_CW)=OxloadCW(iCW)
     FeloadCW_tot(ind_CW)=FeloadCW(iCW)
     vloadCW_tot(ind_CW,1)=vloadCW(iCW,1)
     vloadCW_tot(ind_CW,2)=vloadCW(iCW,2)
     vloadCW_tot(ind_CW,3)=vloadCW(iCW,3)
  enddo
  call MPI_ALLREDUCE(mloadCW_tot,mloadCW_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZloadCW_tot,ZloadCW_mpi,nCW_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(OxloadCW_tot,OxloadCW_mpi,nCW_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(FeloadCW_tot,FeloadCW_mpi,nCW_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vloadCW_tot,vloadCW_mpi,nCW_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  mloadCW_tot=mloadCW_mpi
  ZloadCW_tot=ZloadCW_mpi
  OxloadCW_tot=OxloadCW_mpi
  FeloadCW_tot=FeloadCW_mpi
  vloadCW_tot=vloadCW_mpi

  ! Put the nCW_tot size arrays into nCW size arrays
  do iCW=1,nCW
     ind_CW=iCW_myid(iCW)
     mloadCW(iCW)=mloadCW_tot(ind_CW)
     ZloadCW(iCW)=ZloadCW_tot(ind_CW)
     OxloadCW(iCW)=OxloadCW_tot(ind_CW)
     FeloadCW(iCW)=FeloadCW_tot(ind_CW)
     vloadCW(iCW,1)=vloadCW_tot(ind_CW,1)
     vloadCW(iCW,2)=vloadCW_tot(ind_CW,2)
     vloadCW(iCW,3)=vloadCW_tot(ind_CW,3)
  enddo
#endif


  if(verbose)write(*,*)'Exiting average_CW'

end subroutine average_CW
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine Sedov_blast_CW(xCW,mCW,eCW,indCW,vol_gas,vol_rho,dq,ekBlast,nCW,mloadCW,ZloadCW,OxloadCW,FeloadCW,vloadCW)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,j,iCW,nCW,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_CW,d,u,v,w,ek,u_r,d_gas
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,eCW1,dr_cell_tmp
  real(dp)::scale_nh,scale_t2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nCW)::mCW,eCW,vol_gas,vol_rho,usedov,ekBlast
  real(dp),dimension(1:nCW)::mloadCW,ZloadCW,OxloadCW,FeloadCW,p_gas
  real(dp),dimension(1:nCW,1:3)::xCW,dq,vloadCW
  integer ,dimension(1:nCW)::indCW
  logical ,dimension(1:nvector),save::ok
  real(dp)::pxmi,pxma,pymi,pyma,pzmi,pzma
  real(dp)::scale_m,dum
  integer::its(1:1)
  if(verbose)write(*,*)'entering sedov_blast_CW'

  ! mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**nlevelmax

  ! conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)

  ! maximum radius of the ejecta
  rmax=max(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  scale_m=dble(scale_d)*dble(scale_l)**3d0

  ! To reduce time spent to search for cells
  pxmi=minval(xCW(:,1))
  pxma=maxval(xCW(:,1))
  pymi=minval(xCW(:,2))
  pyma=maxval(xCW(:,2))
  pzmi=minval(xCW(:,3))
  pzma=maxval(xCW(:,3))

  usedov=0
  do iCW=1,nCW
     ! ejecta specific energy (accounting for dilution)
     eCW1=(dble(eCW(iCW))/dble(mCW(iCW)*scale_m))*(1d40/dble(scale_v)**2d0)
     if(vol_gas(iCW)>0d0)then
        d_gas      =mCW(iCW)/vol_gas(iCW)
        p_gas (iCW)=(1d0-f_ek)*d_gas*eCW1
        usedov(iCW)=dsqrt(2*f_ek*mCW(iCW)*eCW1/ekBlast(iCW)/mloadCW(iCW))
!        dum=dsqrt(2*f_ek*eCW1)*scale_v*1d-5
!  if (myid.eq.1.and.usedov(iCW)*scale_v.gt.3d8)write(800,*)'USEDOV', sngl(usedov(iCW)*scale_v*1e-5), sngl(dum),&
!                    & sngl(mCW(iCW)/mloadCW(iCW)),sngl(ekBlast(iCW)),sngl(d_gas*scale_nh)
     endif
  end do


  
  dr_cell_tmp=1d0

  ! loop over levels
  do ilevel=levelmin,nlevelmax
     ! computing local volume (important for averaging hydro quantities) 
     dx=0.5d0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim

     ! cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx
        xc(ind,2)=(dble(iy)-0.5d0)*dx
        xc(ind,3)=(dble(iz)-0.5d0)*dx
     end do

     ! loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=min(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rmax).or.(x.gt.pxma+rmax).or.&
                   &(y.lt.pymi-rmax).or.(y.gt.pyma+rmax).or.&
                   &(z.lt.pzmi-rmax).or.(z.gt.pzma+rmax)) then
                    ok(i)=.false.
                 endif
              endif
          enddo

          do i=1,ngrid
             if(ok(i))then
                ! get gas cell position
                x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                do iCW=1,nCW
                   if (vol_gas(iCW)>0d0) then  !added by TS
                      ! check if the cell lies within the CW radius
                      dxx=x-xCW(iCW,1)
                      dyy=y-xCW(iCW,2)
                      dzz=z-xCW(iCW,3)
                      dr_CW=dxx**2+dyy**2+dzz**2
                      if(dr_CW.lt.rmax2)then
                         if (pseudo_sedov)then
                             dr_cell_tmp=dsqrt(dr_CW/rmax2)
                             d_gas=mloadCW(iCW)/vol_rho(iCW)*dr_cell_tmp
                         else
                             d_gas=mloadCW(iCW)/vol_gas(iCW)
                         endif
                         ! compute the density and the metal density of the cell
                         uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas
                         if(metal)then
                            uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_gas*ZloadCW(iCW)
                            if(chemo)then
                               uold(ind_cell(i),imetal+1)=uold(ind_cell(i),imetal+1)+d_gas*OxloadCW(iCW)
                               uold(ind_cell(i),imetal+2)=uold(ind_cell(i),imetal+2)+d_gas*FeloadCW(iCW)
                            endif
                         endif
                         ! velocity at a given dr_CW linearly interpolated between zero and usedov
                         u=usedov(iCW)*(dxx/rmax-dq(iCW,1))+vloadCW(iCW,1)
                         v=usedov(iCW)*(dyy/rmax-dq(iCW,2))+vloadCW(iCW,2)
                         w=usedov(iCW)*(dzz/rmax-dq(iCW,3))+vloadCW(iCW,3)
                         ! add each momentum component of the blast wave to the gas
                         uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u
                         uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v
                         uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w
                         ! finally update the total energy of the gas
                         uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)+p_gas(iCW)
                      endif
                   endif
                end do
             endif
          end do
           
        end do
        ! end loop over cells
     end do
     ! end loop over grids
  end do
  ! end loop over levels

  do iCW=1,nCW
     if(vol_gas(iCW)==0d0)then
        d_gas=mloadCW(iCW)/ekBlast(iCW) 
        u=vloadCW(iCW,1)
        v=vloadCW(iCW,2)
        w=vloadCW(iCW,3)
        if(indCW(iCW)>0)then
           uold(indCW(iCW),1)=uold(indCW(iCW),1)+d_gas
           uold(indCW(iCW),2)=uold(indCW(iCW),2)+d_gas*u
           uold(indCW(iCW),3)=uold(indCW(iCW),3)+d_gas*v
           uold(indCW(iCW),4)=uold(indCW(iCW),4)+d_gas*w
           uold(indCW(iCW),5)=uold(indCW(iCW),5)+d_gas*0.5*(u*u+v*v+w*w)+p_gas(iCW)
           if(metal)then
              uold(indCW(iCW),imetal)=uold(indCW(iCW),imetal)+d_gas*ZloadCW(iCW)
              if(chemo)then
                 uold(indCW(iCW),imetal+1)=uold(indCW(iCW),imetal+1)+d_gas*OxloadCW(iCW)
                 uold(indCW(iCW),imetal+2)=uold(indCW(iCW),imetal+2)+d_gas*FeloadCW(iCW)
              endif
           endif
        endif
     endif
  end do

  if(verbose)write(*,*)'exiting sedov_blast_CW'

end subroutine Sedov_blast_CW
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine stellar_wind_thermal(ilevel)
   use amr_commons
   use pm_commons
   use hydro_commons
   implicit none
   integer::ilevel
   integer::igrid,jgrid,ipart,jpart,next_part
   integer::i,ig,ip,npart1,npart2,icpu
   integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

   if (numbtot(1, ilevel)==0) return
   ! Gather star particles only (for each cpu?)

#if NDIM==3
   ! Loop over cpus
   do icpu=1,ncpu
      igrid=headl(icpu, ilevel) !TS: headl: head grid in the level
      ig=0
      ip=0
      ! Loop over grids
      do jgrid=1,numbl(icpu, ilevel)   !TS: numbl: number of grid in the level
         npart1=numbp(igrid) ! Number of particles in the grid
         npart2=0
         
         ! Count star particles
         if(npart1>0) then
            ipart=headp(igrid)        !TS: headp: head particle in the grid
            ! Loop over particles
            do jpart=1,npart1
               ! Save next particle <--- Very important !!!
               next_part=nextp(ipart)
               if (tp(ipart).ne.0d0)then
                  npart2=npart2+1
               end if
               ipart=next_part ! Go to next particle
            end do
         endif

         ! Gather nvector star particles and compute the feedback 
         !   until we iterate this process for all(=npart2) particles
         if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(tp(ipart).ne.0d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              end if
              if(ip==nvector) then
                 !main loop for physics
                 call feedbk_thermal(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to the next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)  ! Go to the next grid
     end do
     ! End loop over grids
     if(ip>0) then 
        call feedbk_thermal(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
     end if
  end do
  ! End loop over cpus
  
#endif                
111 format('   Entering Stellar_Wind for level ', I2)

end subroutine stellar_wind_thermal
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine feedbk_thermal(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
   use amr_commons
   use pm_commons
   use hydro_commons
   implicit none
   integer::ng,np,ilevel
   integer,dimension(1:nvector)::ind_grid
   integer,dimension(1:nvector)::ind_grid_part,ind_part
   !----------------------------------------------------------------
   ! This routine is called by subroutine stellar_wind. Each star 
   ! particle dumps mass, momentum, and energy in the nearest grid
   ! cell using array uold. This is basically the similar with 
   ! the subroutine feedbk
   ! Taysun Kimm
   !----------------------------------------------------------------
   integer::i,j,idim,nx_loc
   real(dp)::dx,dx_loc,scale,vol_loc,mejecta,zwind,zloss
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_e
   logical ::error
   ! Grid based arrays
   real(dp),dimension(1:nvector,1:ndim),save::x0
   integer ,dimension(1:nvector), save::ind_cell
   integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
   integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
   ! Particle based arrays
   integer ,dimension(1:nvector),save::igrid_son,ind_son
   integer ,dimension(1:nvector),save::list1
   logical ,dimension(1:nvector),save::ok
   real(dp),dimension(1:nvector),save::mloss,ethermal,ekinetic,dteff
   real(dp),dimension(1:nvector),save::mzloss,moxloss,mfeloss
   real(dp),dimension(1:nvector,1:ndim),save::x
   integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
   integer ,dimension(1:nvector),save::igrid,icell,indp,kg
   real(dp),dimension(1:3)::skip_loc
   real(dp),dimension(1:np)::mCWind_pum,eCWind_pum,mpinit,dt_yr
   real(dp),dimension(1:np)::ZCWind,OxCWind,FeCWind
   real(dp),dimension(1:np)::zpsub

   ! Conversion factor from user units to cgs units
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! Energy release from stellar winds+SNIa from erg to code units 
   scale_e=1d40/2d33/scale_v**2   ! the factor 1d40 from cwind table 

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


#if NDIM==3
   ! Lower left corner of 3x3x3 grid-cube
   do idim=1,ndim
      do i=1,ng
         x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
      end do
   end do

   ! Gather 27 neighbouring father cells (should be present anytime !) 
   do i=1,ng
      ind_cell(i)=father(ind_grid(i))
   end do
   call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)
  
   ! Rescale position at level ilevel
   do idim=1,ndim
      do j=1,np
         x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
      end do
   end do
   do idim=1,ndim
      do j=1,np
         x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
      end do
   end do
   do idim=1,ndim
      do j=1,np
         x(j,idim)=x(j,idim)/dx
      end do
   end do


   ! Check for illegal moves
   error=.false.
   do idim=1,ndim
      do j=1,np
         if(x(j,idim)<=0D0.or.x(j,idim)>=6D0)error=.true.
      end do
   end do

   if(error)then
      write(*,*)'problem in stellar_wind'
      write(*,*)ilevel,ng,np
      stop
   endif

   !NGP at level ilevel
   do idim=1,ndim
      do j=1,np
         id(j,idim)=x(j,idim)
      end do
   end do

   ! Compute parent grids
   do idim=1,ndim
      do j=1,np
         igd(j,idim)=id(j,idim)/2
      end do
   end do
   do j=1,np
      kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
   end do
   do j=1,np
      igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
   enddo

   ! Check if particles are entirely in level ilevel
   ok(1:np)=.true.
   do j=1,np
      ok(j)=ok(j).and.igrid(j)>0
   end do

   ! Compute parent cell position
   do idim=1,ndim
      do j=1,np 
         if(ok(j)) then
            icd(j,idim)=id(j,idim)-2*igd(j,idim)
         end if
      end do
   end do
   do j=1,np
      if(ok(j)) then
         icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
      end if
   end do

   ! Compute parent cell addresses
   do j=1,np
      if(ok(j))then
         indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
      end if
   end do 


   ! Compute individual time steps
   do j=1, np
      if(ok(j)) then
         if(levelp(ind_part(j))>=ilevel) then
             dteff(j)=dtnew(levelp(ind_part(j)))
         else
             dteff(j)=dtold(levelp(ind_part(j)))
         end if
      end if
   end do

   ! Reset ejected mass, metallicity, thermal energy
    do j=1, np
       if(ok(j)) then
          mloss(j)=0d0
          mzloss(j)=0d0
          moxloss(j)=0d0
          mfeloss(j)=0d0
          ethermal(j)=0d0
       end if
    end do

    ! Compute stellar mass loss and thermal feedback due to continuous feedback 
    ! (stellar winds + snIa)
    zpsub=0d0
    if(metal)then
       zpsub(1:np)=zp(ind_part(1:np))
    endif

    call Mstar_Init(np,mp(ind_part(1:np)),tp(ind_part(1:np)),zpsub,mpinit)

    call Wind_Prop(np,mp(ind_part(1:np))/mpinit,tp(ind_part(1:np)),zpsub,idp(ind_part(1:np)),&
                  &mCWind_pum,eCWind_pum,ZCWind,OxCWind,FeCWind,dt_yr,.true.) ! per unit mass

    do j=1,np
       if(ok(j)) then
          ! Stellar mass loss + snIa
          mejecta = mCWind_pum(j)*mpinit(j)!mp(ind_part(j))
          mloss(j)= mloss(j) + mejecta/vol_loc  !--> density

          ! Thermal energy
          ethermal(j)=ethermal(j)+mpinit(j)*eCWind_pum(j)*scale_e/vol_loc
          ! Metallicity
          if(metal) then
             zloss=ZCWind(j)
             mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc

             if(chemo)then
                moxloss(j)=moxloss(j)+mejecta*OxCWind(j)/vol_loc
                mfeloss(j)=mfeloss(j)+mejecta*FeCWind(j)/vol_loc
             endif
          endif
        
          ! Reduce star particle mass
          mp(ind_part(j))=mp(ind_part(j))-mejecta

       end if
    end do

    ! Update hydro variable due to feedback
    do j=1,np
       if(ok(j)) then
           ! Specific kinetic energy of the star
           ekinetic(j) = 0.5*(vp(ind_part(j),1)**2 &
                             +vp(ind_part(j),2)**2 &
                             +vp(ind_part(j),3)**2)
           ! update hydro variable in NGP call
           uold(indp(j),1)=uold(indp(j),1)+mloss(j)
           uold(indp(j),2)=uold(indp(j),2)+mloss(j)*vp(ind_part(j),1)
           uold(indp(j),3)=uold(indp(j),3)+mloss(j)*vp(ind_part(j),2)
           uold(indp(j),4)=uold(indp(j),4)+mloss(j)*vp(ind_part(j),3)
           uold(indp(j),5)=uold(indp(j),5)+mloss(j)*ekinetic(j)+ethermal(j)
       end if
    end do
    if(metal)then
       do j=1,np
          if(ok(j))then
             uold(indp(j),imetal)=uold(indp(j),imetal)+mzloss(j)
          end if
       end do
       if(chemo)then
          do j=1,np
             if(ok(j))then
                uold(indp(j),imetal+1)=uold(indp(j),imetal+1)+moxloss(j)
             end if
          end do
          do j=1,np
             if(ok(j))then
                uold(indp(j),imetal+2)=uold(indp(j),imetal+2)+mfeloss(j)
             end if
          end do
       endif
    end if
#endif

end subroutine feedbk_thermal
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine Wind_Prop(np,mpfrac,tpsub,zpsub,idpsub,mCWind,eCWind,&
                    &ZCWind,OxCWind,FeCWind,dtlast,isThermal) 
   use amr_commons
   use stellar_commons 
   implicit none
   real(kind=dp),dimension(1:np),intent(inout):: mCWind,eCWind,dtlast
   real(kind=dp),dimension(1:np),intent(inout):: ZCWind,OxCWind,FeCWind
   real(kind=dp),dimension(1:np),intent(in)   :: mpfrac,tpsub,zpsub
   integer,intent(in):: np,idpsub(1:np)
   logical,intent(in):: isThermal

   real(kind=dp):: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
   integer:: i,j,k,ii,izc(1:1)
   real(kind=8):: f_cml, ixtmp, age_star, t_last
   real(kind=8):: t_uni,t_H0,cml_1,cml_2,tp_lb,s2yr,cel_1,cel_2,cmz_1,cmz_2
   real(kind=8):: cmox_1,cmox_2,cmfe_1,cmfe_2

   ZCWind = 0

   s2yr=1d0/(3600d0*24d0*365d0)
   t_H0=1d0/h0*(3.08d13*1d6)/(3600d0*24d0*365d0)

   ! Compute a couple of useful conformal times and look-back times
   if (cosmo) then
      call cv_tau_tlb(dble(t), t_uni)
   else
      call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2) 
      t_uni = dble(t*(scale_t*s2yr))
   end if

   ! Main loop (age_star > t_last)
   do j=1, np
     ! Age of a star in yr unit
     if(cosmo) then
        call cv_tau_tlb(dble(tpsub(j)),tp_lb)
        age_star=(t_uni-tp_lb)*t_H0      ! yr
     else
        age_star=dble((t-tpsub(j))*(scale_t*s2yr)) ! not yet updated age
     end if

     ! Metallicity grid for a star (no interpolation)
     if (abs(sngl(zpsub(j))).le.1e-8) then
        izc=1
     else
        izc=minloc(dabs(dlog10(Z_wind)-dlog10(zpsub(j))))
     end if

     ! Compute mass loss, energy release, and yield at t=age_star
     if(age_star.le.t_wind(1)) then
        cml_1=0d0
        cel_1=0d0
        cmz_1=0d0
     else
        ! Binary search for age_star
        call binary_search(t_wind,age_star,nbint,i)
        ixtmp = (dlog10(t_wind(i+1))-dlog10(age_star))/&
               &(dlog10(t_wind(i+1))-dlog10(t_wind(i)))
        if (ixtmp.lt.0) ixtmp = 0    ! no extrapolation
        if (ixtmp.gt.1) ixtmp = 1    ! no extrapolation
        cml_1 = 10d0**(cMwind(i,izc(1))*ixtmp + cMwind(i+1,izc(1))*(1d0-ixtmp))
        cel_1 = cEwind(i,izc(1))*ixtmp + cEwind(i+1,izc(1))*(1d0-ixtmp)
        if (metal)then
           cmz_1 = 10d0**(cMZwind(i,izc(1))*ixtmp + cMZwind(i+1,izc(1))*(1d0-ixtmp))
           if(chemo)then
              cmox_1 = 10d0**(cMOxwind(i,izc(1))*ixtmp + cMOxwind(i+1,izc(1))*(1d0-ixtmp))
              cmfe_1 = 10d0**(cMFewind(i,izc(1))*ixtmp + cMFewind(i+1,izc(1))*(1d0-ixtmp))
           endif
        endif

     endif

     ! Find t_last
     if (dabs(mpfrac(j)-1d0).le.1d-6.or.age_star.le.t_wind(1)) then
         t_last = 0
     else 
        call binary_search(cMwind(:,izc(1)),dlog10(1d0-mpfrac(j)),nbint,i)
        ixtmp = (cMwind(i+1,izc(1))-dlog10(1d0-mpfrac(j)))/&
               &(cMwind(i+1,izc(1))-cMwind(i,izc(1)))
        if (ixtmp.lt.0) ixtmp = 0    ! no extrapolation
        if (ixtmp.gt.1) ixtmp = 1    ! no extrapolation
        t_last = 10d0**(dlog10(t_wind(i))*ixtmp+dlog10(t_wind(i+1))*(1d0-ixtmp))

     endif


     ! Compute mass loss, energy release, and yield at t=t_last
     if(t_last.le.t_wind(1).or.t_last.gt.age_star) then
        cml_2=0d0; cel_2=0d0; cmz_2=0d0; cmox_2=0d0; cmfe_2=0d0
     else 
        ixtmp     = (dlog10(t_wind(i+1))-dlog10(t_last))/&
                   &(dlog10(t_wind(i+1))-dlog10(t_wind(i)))
        if (ixtmp.lt.0) ixtmp=0    ! no extrapolation
        if (ixtmp.gt.1) ixtmp=1    ! no extrapolation
        cml_2  = 10**(cMwind(i,izc(1))*ixtmp + cMwind(i+1,izc(1))*(1d0-ixtmp))
        cel_2  = cEwind(i,izc(1))*ixtmp + cEwind(i+1,izc(1))*(1d0-ixtmp)
        if (metal)then
           cmz_2 = 10d0**(cMZwind(i,izc(1))*ixtmp + cMZwind(i+1,izc(1))*(1d0-ixtmp))
           if(chemo)then
              cmox_2 = 10d0**(cMOxwind(i,izc(1))*ixtmp + cMOxwind(i+1,izc(1))*(1d0-ixtmp))
              cmfe_2 = 10d0**(cMFewind(i,izc(1))*ixtmp + cMFewind(i+1,izc(1))*(1d0-ixtmp))
           endif
        endif
     end if
  
     dtlast(j)=age_star-t_last   ! [yr]

     if  ((isThermal.and.(age_star.le.t_ctw)).or.&
       & (dtlast(j)<0.1).or.&
       & (cml_1-cml_2).le.1d-40) then
                      ! 0.1 (instead of 0) is for numerical stability
        mCwind(j)=0d0
        eCwind(j)=0d0
        dtlast(j)=0d0
        if(metal)then
           ZCWind(j)=0d0
           if(chemo)then
              OxCWind(j)=0d0;FeCWind(j)=0d0
           endif
        endif
     else
        mCWind(j) = (cml_1 - cml_2)
        eCWind(j) = (10d0**(cel_1-40d0) - 10d0**(cel_2-40d0))*eff_sfbk
        if(metal)then 
           ZCWind(j) = (cmz_1-cmz_2)/(cml_1-cml_2)
           if(chemo)then
              OxCWind(j) = (cmox_1-cmox_2)/(cml_1-cml_2)
              FeCWind(j) = (cmfe_1-cmfe_2)/(cml_1-cml_2)
           endif
        endif
        if(mCWind(j)<0)then ! suppress any type of numerical error
           mCWind(j)=0d0
           eCWind(j)=0d0
           dtlast(j)=0d0
           if(metal)then
              ZCWind(j)=0d0
              if(chemo)then
                 OxCWind(j)=0d0;FeCWind(j)=0d0
              endif
           endif
        endif
     endif
 

   enddo
end subroutine Wind_Prop
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine binary_search(database,xtarget,ndata,i)
   implicit none
   integer::i,j,k
   integer,intent(in)::ndata
   real(kind=8),intent(in)::database(1:ndata),xtarget

   i=1
   j=ndata
   do
     k=(i+j)/2
     if (xtarget<database(k)) then
         j=k
     else
         i=k
     end if
     if (i+1>=j) exit
   end do

end subroutine binary_search
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine cv_tlb_tau(tinput, toutput)
   use amr_commons
   implicit none
   real(kind=8),intent(in)   :: tinput 
   real(kind=8),intent(inout):: toutput 
   real(kind=8)              :: t_H0, t_uni, ixtmp
   integer                    :: i,j,k

   t_H0 = 1d0/h0*(3.08d13*1d6)/(3600d0*24d0*365d0)
   i = 0
   j = 1000 
   do
     k = (i+j)/2
     if (dabs(tinput) < dabs(t_frw(k))) then
         j = k
     else
         i = k
     end if
     if (i+1 >= j) exit
   end do
   ixtmp   = dble((t_frw(i+1)-tinput)/(t_frw(i+1)-t_frw(i)))
   toutput = dble((tau_frw(i)*ixtmp + tau_frw(i+1)*(1d0-ixtmp)))
end subroutine cv_tlb_tau
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine cv_tau_tlb(tinput, toutput)
   use amr_commons
   implicit none
   real(kind=8),intent(in):: tinput 
   real(kind=8),intent(out):: toutput 
   real(kind=8):: t_H0, t_uni, ixtmp
   integer:: i,j,k

   i = 0
   j = 1000 
   do
     k = (i+j)/2
     if (dabs(tinput) < dabs(tau_frw(k))) then
         j = k
     else
         i = k
     end if
     if (i+1 >= j) exit
   end do
   ixtmp   = dble((tau_frw(i+1)-tinput)/(tau_frw(i+1)-tau_frw(i)))
   toutput = dble((t_frw(i)*ixtmp + t_frw(i+1)*(1d0-ixtmp)))  ! 13.7Gyr <-> -1
   if (toutput.lt.-1d0) toutput=-1d0 
end subroutine cv_tau_tlb
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine Mstar_Init(np,mpsub,tpsub,zpsub,mpinit)
  use amr_commons
  use stellar_commons
  implicit none
  integer::np
  real(kind=dp),dimension(1:np)::mpsub,tpsub,zpsub,mpinit
  real(kind=dp)::scale_nh,scale_t2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx_min,vol_min,scale,mstar1,mp_estimate
  integer::nx_loc

  integer:: i,j,ii,izc(1:1)
  real(kind=8):: f_cml, ixtmp, age_star, mpfrac
  real(kind=8):: t_uni,t_H0,tp_lb,s2yr

  ! conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2) 

  ! mesh spacing in that level
  nx_loc  = (icoarse_max-icoarse_min+1)
  scale   = boxlen/dble(nx_loc)
  dx_min  = (0.5d0**nlevelmax)*scale
  vol_min = dx_min**ndim
  mstar1  = n_star/(scale_nh*aexp**3d0)*vol_min

  s2yr=1d0/(3600d0*24d0*365d0)
  t_H0=1d0/h0*(3.08d13*1d6)/(3600d0*24d0*365d0)

  ! Compute a couple of useful conformal times and look-back times
  if (cosmo) then
     call cv_tau_tlb(dble(t), t_uni)
  else
    t_uni = dble(t*(scale_t*s2yr))
  end if


  ! Main loop 
  do j=1, np
    ! Age of a star in yr unit
    if(cosmo) then
       call cv_tau_tlb(dble(tpsub(j)),tp_lb)
       age_star=(t_uni-tp_lb)*t_H0      ! yr
    else
       age_star=dble((t-tpsub(j))*(scale_t*s2yr)) ! yr 
    end if

    ! Metallicity grid for a star (no interpolation)
    if (abs(sngl(zpsub(j))).le.1e-8) then
       izc=1
    else
       izc=minloc(dabs(dlog10(Z_wind)-dlog10(zpsub(j))))
    end if

    ! Compute mass loss, energy release, and yield at t=age_star
    if(age_star.le.t_wind(1)) then
       mpfrac=1d0
    else
       ! Binary search for age_star
       call binary_search(t_wind,age_star,nbint,i)
       ixtmp = (dlog10(t_wind(i+1))-dlog10(age_star))/&
              &(dlog10(t_wind(i+1))-dlog10(t_wind(i)))
       if (ixtmp.lt.0) ixtmp = 0    ! no extrapolation
       if (ixtmp.gt.1) ixtmp = 1    ! no extrapolation
       mpfrac = 1d0-10d0**(cMwind(i,izc(1))*ixtmp + cMwind(i+1,izc(1))*(1d0-ixtmp))
    endif

    mp_estimate = mstar1*nint(mpsub(j)/mpfrac/mstar1)
    mpinit(j)   = mp_estimate
  end do
end subroutine Mstar_Init
!####################################################################
!####################################################################
!####################################################################
!####################################################################
