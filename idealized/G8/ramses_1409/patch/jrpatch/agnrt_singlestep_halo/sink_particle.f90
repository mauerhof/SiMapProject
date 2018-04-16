!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !-------------------------------------------------------------------------------
  ! Description: This subroutine create sink particle in cells where some density 
  ! threshold has been crossed. It also removes from the gas the corresponding 
  ! particle mass. On exit, all fluid variables in the cell are modified.
  ! This routine is called only once per coarse step by routine amr_step.
  ! Romain Teyssier, October 7th, 2007
  !-------------------------------------------------------------------------------
  ! local constants
  integer::ilevel,ivar,info,icpu,igrid,npartbound,isink,nelvelmax_loc
  real(dp)::dx_min,vol_min,dx,dx_temp
  integer::nx_loc

  if(verbose)write(*,*)' Entering create_sink'

  ! Remove particles to finer levels
  do ilevel=levelmin,nlevelmax
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! Get the star density value in each cell
  do ilevel=levelmin,nlevelmax
     call get_rho_star(ilevel)
  enddo

  ! Create new sink particles
  ! and gather particle from the grid
  call make_sink(nlevelmax)
  do ilevel=nlevelmax-1,1,-1
     if(ilevel>=levelmin)call make_sink(ilevel)
     call merge_tree_fine(ilevel)
  end do

  ! Remove particle clouds around old sinks 
  call kill_cloud(1)

  ! update sink position before merging sinks and creating clouds
  call update_sink_position_velocity

  ! Merge sink using FOF
  call merge_sink(1)

  ! Create new particle clouds
  call create_cloud(1)

  ! Scatter particle to the grid
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! Update hydro quantities for split cells
  if(hydro)then
     do ilevel=nlevelmax,levelmin,-1
        call upload_fine(ilevel)
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        ! Update boundaries 
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end do
  end if
  
  jsink=0d0
  ! Compute Bondi parameters and gather particle
  do ilevel=nlevelmax,levelmin,-1
     if(bondi)call bondi_hoyle(ilevel)
     call merge_tree_fine(ilevel)
  end do

!!$  if(myid==1)then
!!$     if(nsink>0)then
!!$        write(*,*)'Sink properties'
!!$        do isink=1,nsink
!!$           if(msink(isink)>0)then
!!$              write(*,'(I8,4(1X,1PE10.3))')isink,msink(isink),xsink(isink,1:ndim)
!!$
!!$              if(sink_AGN)write(*,*)'jsink=',jsink(isink,1:ndim)
!!$           endif
!!$        end do
!!$     endif
!!$  endif

end subroutine create_sink
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_part_from_sink

  ! joki put this here for this patch to compile with 
  ! Ramses version 3.11

end subroutine create_part_from_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: H=>X, rhoc, mH, twopi 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine create sink particle in cells where some
  ! density threshold is crossed. It also removes from the gas the
  ! corresponding particle mass. On exit, all fluid variables in the cell
  ! are modified. This is done only in leaf cells.
  ! Romain Teyssier, October 7th, 2007
  !----------------------------------------------------------------------
  ! local constants
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ::ncache,nnew,ivar,ngrid,icpu,index_sink,index_sink_tot,icloud
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,isink,inew,nx_loc
  integer ::ii,jj,kk,ind_cloud,ncloud
  integer ::ntot,ntot_all,info,ntot_tmp,izero_myid,ninc,ntot_myid
  logical ::ok_free
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::d,x,y,z,u,v,w,e,temp,zg
  real(dp)::velc,uc,vc,wc,l_jeans,d_jeans,d_thres,d_sink
  real(dp)::birth_epoch,xx,yy,zz,rr,rmax_sink2,x_half,y_half,z_half,x_box,y_box,z_box
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min,dxx,dyy,dzz,dr_sink,rmax_sink,rmax
  real(dp)::bx1,bx2,by1,by2,bz1,bz2,factG,pi,nCOM,d_star,star_ratio

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::ind_part_cloud,ind_grid_cloud
  logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
  real(dp),dimension(:,:),allocatable::x_tmp,x_tmp_all
  real(dp),dimension(:)  ,allocatable::dens_tmp,dens_tmp_all
  integer ,dimension(:)  ,allocatable::flag_tmp,flag_tmp_all,point2flag2

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering make_sink for level ',ilevel
  
  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  ! Minimum radius to create a new sink from another
  ! (physical kpc -> code units) 
  rmax_sink=r_gal*3.08d21/scale_l

  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp
  
  ! Density threshold for sink particle creation
  d_sink=n_sink/scale_nH
  d_star=0d0
  if (star)d_star=n_star/scale_nH

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  rmax=4.0D0*dx_min

  x_half=scale*xbound(1)/2.0; y_half=scale*xbound(2)/2.0; z_half=scale*xbound(3)/2.0
  x_box =scale*xbound(1); y_box =scale*xbound(2); z_box =scale*xbound(3)
  rmax_sink2=rmax_sink**2

  ! Birth epoch
  if(use_proper_time)then
     birth_epoch=texp
  else
     birth_epoch=t
  endif

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  xx=0.0; yy=0.0;zz=0.0
  ncloud=0
#if NDIM==3
  do kk=-8,8
     zz=dble(kk)*dx_min/2.0
#endif
#if NDIM>1
     do jj=-8,8
	yy=dble(jj)*dx_min/2.0
#endif
	do ii=-8,8
	   xx=dble(ii)*dx_min/2.0
	   rr=sqrt(xx*xx+yy*yy+zz*zz)
	   if(rr<=rmax)ncloud=ncloud+1
	end do
#if NDIM>1
     end do
#endif
#if NDIM==3
  end do
#endif

  ! Set new sink variables to old ones
  msink_new=0d0; xsink_new=0d0; dMsmbh_new=0d0; dMsmbh_coarse_new=0d0; Esave_new=0d0; vsink_new=0d0; oksink_new=0d0; tsink_new=0d0; idsink_new=0

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
           ! No AGN formation site by default
           flag2(ind_cell(i))=1
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        enddo
     end do
  end do

  !------------------------------
  ! Determine AGN formation sites
  !------------------------------
  call quenching(ilevel)

  !----------------------------
  ! Compute number of new sinks
  !----------------------------
  ntot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Density threshold crossed ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Create new sink if the gas density exceed some threshold
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           star_ratio=rho_star(ind_cell(i))/(d+rho_star(ind_cell(i)))

           ! Jeans length related density threshold
           temp=max(uold(ind_cell(i),5)*(gamma-1.0),smallc**2)
           d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
           d_thres=d_jeans

           ! User defined density threshold
           !d_thres=d_sink
           
           ! Check if the stellar density is higher than a star formation threshold
           if(star.and.rho_star(ind_cell(i))<d_star)ok(i)=.false.
           ! Check if the star ratio is higher than a certain fraction
           !if(star.and.star_ratio<star_ratio_floor)ok(i)=.false.
           ! Check if the density is higher than star formation threshold
           if(d    <d_star )ok(i)=.false.
           ! Check if gas is Jeans unstable
           if(d    <d_thres)ok(i)=.false.
           ! Quenching criterion
           !if(flag2(ind_cell(i))==1)ok(i)=.false.
           
           if(ok(i))then
              x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
              y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
              z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
              do isink=1,nsink
                 dxx=x-xsink(isink,1)
                 if(dxx> x_half)then
                    dxx=dxx-x_box
                 endif
                 if(dxx<-x_half)then
                    dxx=dxx+x_box
                 endif
                 dr_sink=dxx*dxx
                 dyy=y-xsink(isink,2)
                 if(dyy> y_half)then
                    dyy=dyy-y_box
                 endif
                 if(dyy<-y_half)then
                    dyy=dyy+y_box
                 endif
                 dr_sink=dyy*dyy+dr_sink
                 dzz=z-xsink(isink,3)
                 if(dzz> z_half)then
                    dzz=dzz-z_box
                 endif
                 if(dzz<-z_half)then
                    dzz=dzz+z_box
                 endif
                 dr_sink=dzz*dzz+dr_sink
                 if(dr_sink .le. rmax_sink2)ok(i)=.false.
              enddo
           endif

        end do
        ! Calculate number of new sinks in each cell plus cloud particles
        do i=1,ngrid
           flag2(ind_cell(i))=0
           if(ok(i))then
              ntot=ntot+1
              flag2(ind_cell(i))=1
           endif
        enddo
     end do
  end do

  !--------------------------------------------------------------------------------------
  !------NEW: This part avoids multiple sink creation at same coarse time step ----------
  !--------------------------------------------------------------------------------------
  allocate(point2flag2(1:ntot))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_tmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_tmp=ntot
#endif
  allocate(x_tmp(1:ntot_tmp,1:3),dens_tmp(1:ntot_tmp),flag_tmp(1:ntot_tmp))
  x_tmp=0d0;dens_tmp=0d0;flag_tmp=0
#ifndef WITHOUTMPI
  ntot_sink_cpu=0; ntot_sink_all=0
  ntot_sink_cpu(myid)=ntot
  call MPI_ALLREDUCE(ntot_sink_cpu,ntot_sink_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_sink_cpu(1)=ntot_sink_all(1)
  do icpu=2,ncpu
     ntot_sink_cpu(icpu)=ntot_sink_cpu(icpu-1)+ntot_sink_all(icpu)
  end do
#endif
  if(myid.gt.1)then
     izero_myid=ntot_sink_cpu(myid-1)
  else
     izero_myid=0
  endif
  ninc=0
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
	! Gather cells with a new sink
        nnew=0
        do i=1,ngrid
           if (flag2(ind_cell(i))>0)then
              ninc=ninc+1
              x_tmp(izero_myid+ninc,1)=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
              x_tmp(izero_myid+ninc,2)=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
              x_tmp(izero_myid+ninc,3)=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
              dens_tmp(izero_myid+ninc)=uold(ind_cell(i),1)
              flag_tmp(izero_myid+ninc)=1
              point2flag2(ninc)=ind_cell(i) ! This is a local pointer that is not shared with other cpus
           end if
        end do
     enddo
  enddo
#ifndef WITHOUTMPI
  allocate(x_tmp_all(1:ntot_tmp,1:3),dens_tmp_all(1:ntot_tmp),flag_tmp_all(1:ntot_tmp))
  call MPI_ALLREDUCE(x_tmp   ,x_tmp_all   ,ntot_tmp*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dens_tmp,dens_tmp_all,ntot_tmp  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(flag_tmp,flag_tmp_all,ntot_tmp  ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  x_tmp   =x_tmp_all
  dens_tmp=dens_tmp_all
  flag_tmp=flag_tmp_all
  deallocate(x_tmp_all,dens_tmp_all,flag_tmp_all)
#endif
  do i=1,ntot_tmp
     if(flag_tmp(i).eq.1)then
        x=x_tmp(i,1);y=x_tmp(i,2);z=x_tmp(i,3)
        do j=i+1,ntot_tmp
           if(flag_tmp(j).eq.1)then
              dxx=x-x_tmp(j,1)
              if(dxx> x_half)then
                 dxx=dxx-x_box
              endif
              if(dxx<-x_half)then
                 dxx=dxx+x_box
              endif
              dr_sink=dxx*dxx
              dyy=y-x_tmp(j,2)
              if(dyy> y_half)then
                 dyy=dyy-y_box
              endif
              if(dyy<-y_half)then
                 dyy=dyy+y_box
              endif
              dr_sink=dyy*dyy+dr_sink
              dzz=z-x_tmp(j,3)
              if(dzz> z_half)then
                 dzz=dzz-z_box
              endif
              if(dzz<-z_half)then
                 dzz=dzz+z_box
              endif
              dr_sink=dzz*dzz+dr_sink
              if(dr_sink .le. rmax_sink2)then
                 ! Keep the largest gas density region
                 if(dens_tmp(i).ge.dens_tmp(j))then
                    flag_tmp(j)=0
                 else
                    flag_tmp(i)=0                    
                 endif
              endif
           endif
        enddo
     endif
  enddo
  ntot_myid=ntot
  ntot=0
  do i=1,ntot_myid
     if(flag_tmp(izero_myid+i).eq.1)then
        ntot=ntot+1
     else
        flag2(point2flag2(i))=0
     endif
  enddo
  deallocate(x_tmp,dens_tmp,flag_tmp,point2flag2)
  !--------------------------------------------------------------------------------------
  !------NEW: This part avoids multiple sink creation at same coarse time step ----------
  !--------------------------------------------------------------------------------------

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------

  ok_free=(numbp_free-ntot*ncloud)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)ncloud,ntot
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if
  
  !---------------------------------
  ! Compute global sink statistics
  !---------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
#endif
#ifndef WITHOUTMPI
  ntot_sink_cpu=0; ntot_sink_all=0
  ntot_sink_cpu(myid)=ntot
  call MPI_ALLREDUCE(ntot_sink_cpu,ntot_sink_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_sink_cpu(1)=ntot_sink_all(1)
  do icpu=2,ncpu
     ntot_sink_cpu(icpu)=ntot_sink_cpu(icpu-1)+ntot_sink_all(icpu)
  end do
#endif

  nsink=nsink+ntot_all
  nindsink=nindsink+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level = ",I6," New sink = ",I6," Tot =",I8)')ilevel,ntot_all,nsink
     endif
  end if

  !------------------------------
  ! Create new sink particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_sink=nsink-ntot_all
     index_sink_tot=nindsink-ntot_all
  else
     index_sink=nsink-ntot_all+ntot_sink_cpu(myid-1)
     index_sink_tot=nindsink-ntot_all+ntot_sink_cpu(myid-1)
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
        
        ! Gather cells with a new sink
        nnew=0
        do i=1,ngrid
           if (flag2(ind_cell(i))>0)then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do
        
        ! Create new sink particles
        do i=1,nnew
           index_sink=index_sink+1
           index_sink_tot=index_sink_tot+1

           ! Get gas variables
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           e=uold(ind_cell_new(i),5)

           ! Get gas cell position
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale

           ! Mass of the new sink

           ! Jeans length related density threshold
           temp=max(e*(gamma-1.0),smallc**2)
           d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
           d_thres=d_jeans
           !d_thres=0.25d0*d

           ! User defined density threshold
           !d_thres=d_sink

           msink_new (index_sink)=min((d-d_thres/4.0)*dx_loc**3,Mseed*2d33/scale_m)
           dMsmbh_new(index_sink)=0d0
           dMsmbh_coarse_new(index_sink)=0d0
           Esave_new (index_sink)=0d0
           oksink_new(index_sink)=1d0

           ! Update linked list
           ind_grid_cloud(1)=ind_grid_new(i)
           call remove_free(ind_part_cloud,1)
           call add_list(ind_part_cloud,ind_grid_cloud,ok_true,1)
           ind_cloud=ind_part_cloud(1)

           ! Set new sink particle variables
           tp(ind_cloud)=0d0             ! Birth epoch
           mp(ind_cloud)=msink_new(index_sink) ! Mass
           levelp(ind_cloud)=ilevel      ! Level
           idp(ind_cloud)=-index_sink    ! Identity
           xp(ind_cloud,1)=x
           xp(ind_cloud,2)=y
           xp(ind_cloud,3)=z
           vp(ind_cloud,1)=u
           vp(ind_cloud,2)=v
           vp(ind_cloud,3)=w
           idsink_new(index_sink) =index_sink_tot
           tsink_new(index_sink)  =birth_epoch
           xsink_new(index_sink,1)=x
           xsink_new(index_sink,2)=y
           xsink_new(index_sink,3)=z
           vsink_new(index_sink,1)=u
           vsink_new(index_sink,2)=v
           vsink_new(index_sink,3)=w
           
           uold(ind_cell_new(i),1)=d-msink_new(index_sink)
           
        end do
        ! End loop over new sink particle cells

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

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax     ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dMsmbh_new,dMsmbh_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dMsmbh_coarse_new,dMsmbh_coarse_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Esave_new ,Esave_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  oksink_all=oksink_new
  msink_all =msink_new
  xsink_all =xsink_new
  vsink_all =vsink_new
  tsink_all =tsink_new
  idsink_all=idsink_new
  dMsmbh_all=dMsmbh_new
  dMsmbh_coarse_all=dMsmbh_coarse_new
  Esave_all =Esave_new
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        tsink(isink) =tsink_all(isink)
        idsink(isink)=idsink_all(isink)
        msink(isink) =msink_all(isink)
        dMsmbh(isink)=dMsmbh_all(isink)
        dMsmbh_coarse(isink)=dMsmbh_coarse_all(isink)
        Esave (isink)=Esave_all (isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
     endif
  end do

#endif

end subroutine make_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_sink(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine merges sink usink the FOF algorithm.
  ! It keeps only the group centre of mass and remove other sinks.
  !------------------------------------------------------------------------
  integer::j,isink,ii,jj,kk,ind,idim,new_sink
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax2,rmax
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer::igrp,icomp,gndx,ifirst,ilast,indx
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  integer,dimension(:),allocatable::psink,gsink
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp)::dx,vol_min

  integer,dimension(:),allocatable::rank_old,idsink_old
  real(dp),dimension(:),allocatable::tsink_old
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::d0,mstar,nISM,nCOM,vr,cs

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(nsink==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  rmax=4.0D0*dx_min
  rmax2=rmax*rmax

  allocate(psink(1:nsink),gsink(1:nsink))

  allocate(rank_old(1:nsink),idsink_old(1:nsink),tsink_old(1:nsink))
  
  !-------------------------------
  ! Merge sinks using FOF
  !-------------------------------
  do isink=1,nsink
     psink(isink)=isink
     gsink(isink)=0
  end do
  
  igrp=0
  icomp=1
  ifirst=2
  do while(icomp.le.nsink)
     gndx=psink(icomp)
     if(gsink(gndx)==0)then
        igrp=igrp+1
        gsink(gndx)=igrp
     endif
     ilast=nsink
     do while((ilast-ifirst+1)>0)
        indx=psink(ifirst)
        xx=xsink(indx,1)-xsink(gndx,1)
        if(xx>scale*xbound(1)/2.0)then
           xx=xx-scale*xbound(1)
        endif
        if(xx<-scale*xbound(1)/2.0)then
           xx=xx+scale*xbound(1)
        endif
        rr=xx**2
#if NDIM>1
        yy=xsink(indx,2)-xsink(gndx,2)
        if(yy>scale*xbound(2)/2.0)then
           yy=yy-scale*xbound(2)
        endif
        if(yy<-scale*xbound(2)/2.0)then
           yy=yy+scale*xbound(2)
        endif
        rr=yy**2+rr
#endif
#if NDIM>2
        zz=xsink(indx,3)-xsink(gndx,3)
        if(zz>scale*xbound(3)/2.0)then
           zz=zz-scale*xbound(3)
        endif
        if(zz<-scale*xbound(3)/2.0)then
           zz=zz+scale*xbound(3)
        endif
        rr=zz**2+rr
#endif
!!$        vr=(vsink(indx,1)-vsink(gndx,1))**2+(vsink(indx,2)-vsink(gndx,2))**2+(vsink(indx,3)-vsink(gndx,3))**2
!!$        cs=0.5d0*(c2sink(indx)+c2sink(gndx))
!!$        if(rr.le.rmax2.and.(vr.le.cs .or. cs.eq.0.0))then
        if(rr.le.rmax2)then
           ifirst=ifirst+1
           gsink(indx)=igrp
        else
           psink(ifirst)=psink(ilast)
           psink(ilast)=indx
           ilast=ilast-1
        endif
     end do
     icomp=icomp+1
  end do
  new_sink=igrp
  if(myid==1)then
     write(*,*)'Found ',new_sink,' groups'
     !do isink=1,nsink
     !   write(*,'(3(I4,1x),3(1PE10.3))')isink,psink(isink),gsink(isink),xsink(isink,1:ndim)
     !end do
  endif
  
  !----------------------------------------------------
  ! Compute group centre of mass and average velocity
  !----------------------------------------------------
  xsink_new=0d0; vsink_new=0d0; msink_new=0d0; dMsmbh_new=0d0; dMsmbh_coarse_new=0d0
  Esave_new=0d0; idsink_new=0
  oksink_all=0d0; oksink_new=0d0; tsink_new=0d0
  rank_old=0d0; idsink_old=0d0; tsink_old=0d0
  do isink=1,nsink
     igrp=gsink(isink)
     
     !----------------------------------------------------
     ! This is done to keep track of the most massive sink
     ! after a merger with a companion
     !----------------------------------------------------
     if ( rank_old(igrp) .eq. 0)then
        rank_old(igrp)=isink
        idsink_old(igrp)=idsink(isink)
        tsink_old(igrp) =tsink(isink)
     endif
     if ( msink(isink) .gt. msink(rank_old(igrp)) )then
        rank_old(igrp)=isink
        idsink_new(igrp)=idsink(isink)
        idsink_old(igrp)=idsink(isink)
        tsink_new(igrp) =tsink(isink)
        tsink_old(igrp) =tsink(isink)
     else
        idsink_new(igrp)=idsink_old(igrp)
        tsink_new(igrp) =tsink_old(igrp)
     endif
     !----------------------------------------------------
     !----------------------------------------------------

     !idsink_new(igrp)=idsink(isink)
     if(oksink_new(igrp)==0d0)then
        oksink_all(isink)=igrp
        oksink_new(igrp)=isink
     endif
     msink_new (igrp)=msink_new (igrp)+msink (isink)
     dMsmbh_new(igrp)=dMsmbh_new(igrp)+dMsmbh(isink)
     dMsmbh_coarse_new(igrp)=dMsmbh_coarse_new(igrp)+dMsmbh_coarse(isink)
     Esave_new (igrp)=Esave_new (igrp)+Esave (isink)
     xx=xsink(isink,1)-xsink(int(oksink_new(igrp)),1)
     if(xx>scale*xbound(1)/2.0)then
        xx=xx-scale*xbound(1)
     endif
     if(xx<-scale*xbound(1)/2.0)then
        xx=xx+scale*xbound(1)
     endif
     xsink_new(igrp,1)=xsink_new(igrp,1)+msink(isink)*xx
     vsink_new(igrp,1)=vsink_new(igrp,1)+msink(isink)*vsink(isink,1)
#if NDIM>1
     yy=xsink(isink,2)-xsink(int(oksink_new(igrp)),2)
     if(yy>scale*xbound(2)/2.0)then
        yy=yy-scale*xbound(2)
     endif
     if(yy<-scale*xbound(2)/2.0)then
        yy=yy+scale*xbound(2)
     endif
     xsink_new(igrp,2)=xsink_new(igrp,2)+msink(isink)*yy
     vsink_new(igrp,2)=vsink_new(igrp,2)+msink(isink)*vsink(isink,2)
#endif
#if NDIM>2
     zz=xsink(isink,3)-xsink(int(oksink_new(igrp)),3)
     if(zz>scale*xbound(3)/2.0)then
        zz=zz-scale*xbound(3)
     endif
     if(zz<-scale*xbound(3)/2.0)then
        zz=zz+scale*xbound(3)
     endif
     xsink_new(igrp,3)=xsink_new(igrp,3)+msink(isink)*zz
     vsink_new(igrp,3)=vsink_new(igrp,3)+msink(isink)*vsink(isink,3)
#endif
  end do
  do isink=1,new_sink
     xsink_new(isink,1)=xsink_new(isink,1)/msink_new(isink)+xsink(int(oksink_new(isink)),1)
     vsink_new(isink,1)=vsink_new(isink,1)/msink_new(isink)
#if NDIM>1
     xsink_new(isink,2)=xsink_new(isink,2)/msink_new(isink)+xsink(int(oksink_new(isink)),2)
     vsink_new(isink,2)=vsink_new(isink,2)/msink_new(isink)
#endif
#if NDIM>2
     xsink_new(isink,3)=xsink_new(isink,3)/msink_new(isink)+xsink(int(oksink_new(isink)),3)
     vsink_new(isink,3)=vsink_new(isink,3)/msink_new(isink)
#endif
  end do
  nsink=new_sink
  msink (1:nsink)=msink_new (1:nsink)
  dMsmbh(1:nsink)=dMsmbh_new(1:nsink)
  dMsmbh_coarse(1:nsink)=dMsmbh_coarse_new(1:nsink)
  Esave (1:nsink)=Esave_new (1:nsink)
  idsink(1:nsink)=idsink_new(1:nsink)
  tsink (1:nsink)=tsink_new (1:nsink)
  xsink(1:nsink,1:ndim)=xsink_new(1:nsink,1:ndim)
  vsink(1:nsink,1:ndim)=vsink_new(1:nsink,1:ndim)
  ! Periodic boundary conditions
  do isink=1,nsink
     xx=xsink(isink,1)
     if(xx<-scale*skip_loc(1))then
        xx=xx+scale*(xbound(1)-skip_loc(1))
     endif
     if(xx>scale*(xbound(1)-skip_loc(1)))then
        xx=xx-scale*(xbound(1)-skip_loc(1))
     endif
     xsink(isink,1)=xx
#if NDIM>1
     yy=xsink(isink,2)
     if(yy<-scale*skip_loc(2))then
        yy=yy+scale*(xbound(2)-skip_loc(2))
     endif
     if(yy>scale*(xbound(2)-skip_loc(2)))then
        yy=yy-scale*(xbound(2)-skip_loc(2))
     endif
     xsink(isink,2)=yy
#endif
#if NDIM>2
     zz=xsink(isink,3)
     if(zz<-scale*skip_loc(3))then
        zz=zz+scale*(xbound(3)-skip_loc(3))
     endif
     if(zz>scale*(xbound(3)-skip_loc(3)))then
        zz=zz-scale*(xbound(3)-skip_loc(3))
     endif
     xsink(isink,3)=zz
#endif
  enddo

  deallocate(psink,gsink)

  deallocate(rank_old,idsink_old,tsink_old)
  
  !-----------------------------------------------------
  ! Remove sink particles that are part of a FOF group.
  !-----------------------------------------------------
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering merge_sink for level ',I2)

end subroutine merge_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_sink(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine merge_sink
  ! It removes sink particles that are part of a FOF group.
  !-----------------------------------------------------------------------
  integer::j,isink,ii,jj,kk,ind,idim,isink_new
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok

  do j=1,np
     isink=-idp(ind_part(j))
     ok(j)=(oksink_all(isink)==0)
     if(.not. ok(j))then
        isink_new=oksink_all(isink)
        idp(ind_part(j))=-isink_new
        mp(ind_part(j))=msink(isink_new)
        xp(ind_part(j),1)=xsink(isink_new,1)
        vp(ind_part(j),1)=vsink(isink_new,1)
#if NDIM>1
        xp(ind_part(j),2)=xsink(isink_new,2)
        vp(ind_part(j),2)=vsink(isink_new,2)
#endif
#if NDIM>2
        xp(ind_part(j),3)=xsink(isink_new,3)
        vp(ind_part(j),3)=vsink(isink_new,3)
#endif                 
     endif
  end do

  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_grid_part,ok,np)
  call add_free_cond(ind_part,ok,np)

end subroutine kill_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine creates a cloud of test particle around each sink particle.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather sink particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call mk_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call mk_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering create_cloud for level ',I2)

end subroutine create_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine mk_cloud(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine create_cloud.
  !-----------------------------------------------------------------------
  logical::error
  integer::j,isink,ii,jj,kk,ind,idim,nx_loc,ncloud
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax
  ! Particle-based arrays
  integer ,dimension(1:nvector),save::ind_cloud
  logical ,dimension(1:nvector),save::ok_true=.true.
  real(dp)::vol_min,dx
  integer::i

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  rmax=4.0D0*dx_min

  xx=0.0; yy=0.0;zz=0.0
  ncloud=0
#if NDIM==3
  do kk=-8,8
     zz=dble(kk)*dx_min/2.0
#endif
#if NDIM>1
     do jj=-8,8
        yy=dble(jj)*dx_min/2.0
#endif
        do ii=-8,8
           xx=dble(ii)*dx_min/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=rmax)ncloud=ncloud+1
        end do
#if NDIM>1
     end do
#endif
#if NDIM==3
  end do
#endif

#if NDIM==3
  do kk=-8,8
     zz=dble(kk)*dx_min/2.0
#endif
#if NDIM>1
     do jj=-8,8
        yy=dble(jj)*dx_min/2.0
#endif
        do ii=-8,8
           xx=dble(ii)*dx_min/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr>0.and.rr<=rmax)then
              call remove_free(ind_cloud,np)
              call add_list(ind_cloud,ind_grid_part,ok_true,np)
              do j=1,np
                 isink=-idp(ind_part(j))
                 idp(ind_cloud(j))=-isink
                 mp(ind_cloud(j))=msink(isink)/dble(ncloud)
                 xp(ind_cloud(j),1)=xp(ind_part(j),1)+xx
                 vp(ind_cloud(j),1)=vsink(isink,1)
#if NDIM>1
                 xp(ind_cloud(j),2)=xp(ind_part(j),2)+yy
                 vp(ind_cloud(j),2)=vsink(isink,2)
#endif
#if NDIM>2
                 xp(ind_cloud(j),3)=xp(ind_part(j),3)+zz
                 vp(ind_cloud(j),3)=vsink(isink,3)
#endif                 
              end do
           end if
        end do
#if NDIM>1
     end do
#endif
#if NDIM>2
  end do
#endif

  ! Reduce sink particle mass
  do j=1,np
     isink=-idp(ind_part(j))
     mp(ind_part(j))=msink(isink)/dble(ncloud)
     vp(ind_part(j),1)=vsink(isink,1)
     vp(ind_part(j),2)=vsink(isink,2)
     vp(ind_part(j),3)=vsink(isink,3)
  end do

end subroutine mk_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine removes from the list cloud particles and keeps only
  ! sink particles. 
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather sink and cloud particles.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call rm_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        
        igrid=next(igrid)   ! Go to next grid
     end do
     
     ! End loop over grids
     if(ip>0)call rm_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do

111 format('   Entering kill_cloud for level ',I2)

end subroutine kill_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rm_cloud(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine kill_cloud.
  !-----------------------------------------------------------------------
  logical::error
  integer::j,isink,ii,jj,kk,ind,idim,nx_loc
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,r2,r2_eps
  ! Particle-based arrays
  logical,dimension(1:nvector),save::ok
  integer::i
  real(dp)::dx

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  r2_eps=(1d-15*dx_min)**2

  do j=1,np
     isink=-idp(ind_part(j))
     r2=0d0
     do idim=1,ndim
        r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
     end do
     ok(j)=r2>r2_eps
  end do

  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_grid_part,ok,np)
  call add_free_cond(ind_part,ok,np)

end subroutine rm_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_hoyle(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH , twopi
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the parameters of Bondi-Hoyle
  ! accretion for sink particles.
  ! It calls routine bondi_veocity and average_density.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,idim,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::r2,dx_loc,dx_min,scale
  real(dp)::dx,factG,pi
  real(dp)::RandNum,phi,Rrand,SS,CC,UU,csound,turb
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**nlevelmax/aexp

  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Reset new sink variables
  v2sink_new=0d0; c2sink_new=0d0; oksink_new=0d0

  ! Gather sink particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count only sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    npart2=npart2+1
                 end if
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather only sink particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if
                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig   
                 endif
              endif
              if(ip==nvector)then
                 call bondi_velocity(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call bondi_velocity(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(c2sink_new,c2sink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(v2sink_new,v2sink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     oksink_all=oksink_new
     c2sink_all=c2sink_new
     v2sink_all=v2sink_new
#endif
  endif

  do isink=1,nsink
     if(oksink_all(isink)==1d0)then
        c2sink(isink)=c2sink_all(isink)
        v2sink(isink)=v2sink_all(isink)
        ! Compute sink radius
        r2sink(isink)=(factG*msink(isink)/(v2sink(isink)+c2sink(isink)))**2
        !r2sink(isink)=(factG*msink(isink)/(c2sink(isink)))**2

        ! If radius is far smaller than the resolution, spurious velocity 
        ! must not be taken into account
        !if (sqrt(r2sink(isink)) .lt. dx_min/4d0) then
        !   r2sink(isink)=(factG*msink(isink)/c2sink(isink))**2
        !endif

        r2k(isink)   =min(max(r2sink(isink),(dx_min/4.0)**2),(2.*dx_min)**2)
     endif
  end do

  ! Gather sink and cloud particles.
  wdens=0d0; wvol =0d0; wc2=0d0; wmom=0d0; jsink_new=0d0

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call average_density(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 if(sink_AGN .and.(.not.random_jet)) call jet_AGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call average_density(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
     if(ip>0 .and. sink_AGN .and.(.not.random_jet)) call jet_AGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(random_jet)then
     if(myid==1)then
        do isink=1,nsink
           ! Random directions
           call ranf(localseed,RandNum)
           SS =(RandNum-0.5)*2.
           call ranf(localseed,RandNum)
           phi=(RandNum-0.5)*2.*pi
           call ranf(localseed,RandNum)
           UU =RandNum
           Rrand=UU**(1./3.)
           CC=Rrand*sqrt(1.-SS**2.)
           jsink_new(isink,1)=CC*cos(phi)
           jsink_new(isink,2)=CC*sin(phi)
           jsink_new(isink,3)=Rrand*SS
        enddo
     endif
  endif

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(wdens,wdens_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wvol ,wvol_new ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wc2  ,wc2_new  ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wmom ,wmom_new ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(jsink_new,jsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     wdens_new=wdens
     wvol_new=wvol
     wc2_new =wc2
     wmom_new=wmom
     jsink_all=jsink_new
#endif
  endif

  do isink=1,nsink
     weighted_density(isink,ilevel)=wdens_new(isink)
     weighted_volume (isink,ilevel)=wvol_new (isink)
     weighted_momentum(isink,ilevel,1:ndim)=wmom_new(isink,1:ndim)
     weighted_c2     (isink,ilevel)=wc2_new  (isink)
     do i=1,ndim
        jsink(isink,i)=jsink(isink,i)+jsink_all(isink,i)
     enddo
  end do

111 format('   Entering bondi_hoyle for level ',I2)

end subroutine bondi_hoyle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_velocity(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  ! It computes the gas velocity and soud speed in the cell
  ! each sink particle sits in.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::xxx,mmm,r2,v2,c2,d,u,v,w,e,bx1,bx2,by1,by2,bz1,bz2
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::meff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx_min,vol_min

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
  dx_min=scale*0.5d0**nlevelmax/aexp

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
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
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in bondi_velocity'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
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
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
        
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Gather hydro variables
  do j=1,np
     if(ok(j))then
        d=uold(indp(j),1)
        u=uold(indp(j),2)/d
        v=uold(indp(j),3)/d
        w=uold(indp(j),4)/d
        e=uold(indp(j),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        isink=-idp(ind_part(j))
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2
        c2=MAX(gamma*(gamma-1.0)*e,smallc**2)
        ! Relative velocity of the gas in regards of the sink
        u=u-vsink(isink,1)
        v=v-vsink(isink,2)
        w=w-vsink(isink,3)
        v2=(u**2+v**2+w**2)
        v2sink_new(isink)=v2
        c2sink_new(isink)=c2
        oksink_new(isink)=1d0
     endif
  end do

#endif
  
end subroutine bondi_velocity
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_density(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle. Each cloud particle
  ! reads up the value of density, sound speed and velocity from its
  ! position in the grid.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,length,scale,weight,r2
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::dgas,ugas,vgas,wgas,c2gas
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::u,v,w,d,e,v2,c2
  real(dp)::bx1,bx2,by1,by2,bz1,bz2

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
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
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in average_density'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)   
     end if
  end do
#endif
        
  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Gather gas density
  dgas(1:np)=0.0D0
  ugas(1:np)=0.0D0
  vgas(1:np)=0.0D0
  wgas(1:np)=0.0D0
  c2gas(1:np)=0.0D0
  do ind=1,twotondim
     do j=1,np
        dgas(j)=dgas(j)+uold(indp(j,ind),1)*vol(j,ind)
        d=uold(indp(j,ind),1)
        u=uold(indp(j,ind),2)/d
        v=uold(indp(j,ind),3)/d
        w=uold(indp(j,ind),4)/d
        e=uold(indp(j,ind),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j,ind),6)
        by1=uold(indp(j,ind),7)
        bz1=uold(indp(j,ind),8)
        bx2=uold(indp(j,ind),nvar+1)
        by2=uold(indp(j,ind),nvar+2)
        bz2=uold(indp(j,ind),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2
        c2=MAX(gamma*(gamma-1.0)*e,smallc**2)
        ! --------------------------
        ! Volume-weighted quantities
        ! (if you change this, be sure to renormalize
        ! correctly these quantities in grow_bondi)
        ! --------------------------
        !ugas(j)=ugas(j)+u*vol(j,ind)
        !vgas(j)=vgas(j)+v*vol(j,ind)
        !wgas(j)=wgas(j)+w*vol(j,ind)
        !c2gas(j)=c2gas(j)+c2*vol(j,ind)
        ! --------------------------
        ! Mass-weighted quantities
        ! --------------------------
        ugas(j)=ugas(j)+d*u*vol(j,ind)
        vgas(j)=vgas(j)+d*v*vol(j,ind)
        wgas(j)=wgas(j)+d*w*vol(j,ind)
        c2gas(j)=c2gas(j)+d*c2*vol(j,ind)
     end do
  end do

  do j=1,np
     isink=-idp(ind_part(j))
     r2=0d0
     do idim=1,ndim
        r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
     end do
     weight=exp(-r2/r2k(isink))
     wdens(isink)=wdens(isink)+weight*dgas(j)
     wmom(isink,1)=wmom(isink,1)+weight*ugas(j)
     wmom(isink,2)=wmom(isink,2)+weight*vgas(j)
     wmom(isink,3)=wmom(isink,3)+weight*wgas(j)
     wc2  (isink)=wc2  (isink)+weight*c2gas(j)
     wvol (isink)=wvol (isink)+weight
  end do

end subroutine average_density
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_bondi(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH, twopi 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine performs Bondi-Hoyle accretion of the gas onto 
  ! sink particles. On exit, sink mass and velocity are modified.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,idim,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::r2,density,volume
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  integer::ind,ivar,iskip
  real(dp)::alpha,d_star,pi,factG,c2mean,nfloor,sigmav2,v2mean
  real(dp),dimension(1:3)::velocity

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nfloor = omega_b*rhoc/aexp**3*XH/mH/scale_nH

  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Reset new sink variables
  msink_new=0d0; vsink_new=0d0; dMBH_coarse_new=0d0; dMEd_coarse_new=0d0; dMsmbh_new=0d0

  ! Set unew to uold for myid cells
  ! Need unew to get initial density before Bondi accreting mass 
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        enddo
     enddo
  enddo

  d_star=1d100
  if (star)d_star=n_star/scale_nH
  ! From km/s to user units: assume a 10 km/s vel disp in the ISM
  sigmav2=(sigmav_max*1d5/scale_v)**2d0

  c_avgptr=0d0;v_avgptr=0d0;d_avgptr=0d0

  ! Compute Bondi-Hoyle accretion rate
  do isink=1,nsink
     density=0d0
     c2mean=0d0
     velocity=0d0
     volume=0d0
     do i=levelmin,nlevelmax
        density=density+weighted_density(isink,i)
        c2mean =c2mean +weighted_c2     (isink,i)
        velocity(1)=velocity(1)+weighted_momentum(isink,i,1)
        velocity(2)=velocity(2)+weighted_momentum(isink,i,2)
        velocity(3)=velocity(3)+weighted_momentum(isink,i,3)
        volume =volume +weighted_volume (isink,i)
     end do
     density=density/volume
     ! --------------------------
     ! If volume-weighted
     ! --------------------------
     !velocity(1:3)=velocity(1:3)/volume
     !c2mean =c2mean /volume
     ! --------------------------
     ! If mass-weighted
     ! --------------------------
     velocity(1:3)=velocity(1:3)/volume/density
     c2mean =c2mean /volume/density
     ! --------------------------
     v2mean =min(SUM((velocity(1:3)-vsink(isink,1:3))**2),sigmav2)
     total_volume(isink)=volume
     alpha=max((density/d_star)**boost,1d0)
     if(Esave(isink).eq.0d0)then
        dMBHoverdt(isink)=alpha * 4d0*3.1415926*density* (factG*msink(isink))**2 &
             & / (c2mean+v2mean)**1.5d0
     else
        ! Prevent the accretion of new material onto the BH if
        ! energy has not been released in the previous coarse time step
        dMBHoverdt(isink)=0d0
     endif
     dMEdoverdt(isink)=4.*3.1415926*6.67d-8*msink(isink)*scale_m*1.66d-24/(0.1*6.652d-25*3d10) &
          & /(scale_m/scale_t)
     if(dMBHoverdt(isink)/dMEdoverdt(isink) .lt. X_floor)dMBHoverdt(isink)=f_bondi*dMBHoverdt(isink)
     c_avgptr(isink)=dsqrt(c2mean)
     d_avgptr(isink)=density
  end do
  
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call accrete_bondi(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call accrete_bondi(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dMBH_coarse_new,dMBH_coarse_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dMEd_coarse_new,dMEd_coarse_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dMsmbh_new     ,dMsmbh_all     ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     msink_all=msink_new
     vsink_all=vsink_new
     dMBH_coarse_all=dMBH_coarse_new
     dMEd_coarse_all=dMEd_coarse_new
     dMsmbh_all     =dMBsmbh_new
#endif
  endif
  do isink=1,nsink
     vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)+vsink_all(isink,1:ndim)
     msink(isink)       =msink(isink)       +msink_all(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)
     dMBH_coarse(isink)=dMBH_coarse(isink)+dMBH_coarse_all(isink)
     dMEd_coarse(isink)=dMEd_coarse(isink)+dMEd_coarse_all(isink)
     dMsmbh     (isink)=dMsmbh     (isink)+dMsmbh_all     (isink)
  end do

111 format('   Entering grow_bondi for level ',I2)

end subroutine grow_bondi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_bondi(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink,ivar
  real(dp)::xxx,mmm,r2,v2,c2,d,u,v,w,e,bx1,bx2,by1,by2,bz1,bz2
  real(dp),dimension(1:nvar)::z
  real(dp)::dx,dx_loc,scale,vol_loc,weight,acc_mass,temp
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::meff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::floorB,d_ini,dmsink,X_radio,density,volume
  real(dp)::E_acc,d_x,d_y,d_z,uacc,vacc,wacc,ekin,dx_min,vol_min,rmax
  real(dp)::xx,yy,zz,rr,ddebris,drcyl,fekAGN,E_centre,psy,dx_temp
  integer ::ncloud,ii,jj,kk
  real(dp)::E_sp,j_sp,uc,vc,wc,r_close,dr,j_x,j_y,j_z,jtot,dzjet,drjet

  ! Floor to avoid empty cells
  floorB=0.75d0

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
  dx_min=scale*0.5d0**nlevelmax/aexp
  rmax=4.0d0*dx_min

  xx=0.0; yy=0.0;zz=0.0
  ncloud=0
#if NDIM==3
  do kk=-8,8
     zz=dble(kk)*dx_min/2.0
#endif
#if NDIM>1
     do jj=-8,8
        yy=dble(jj)*dx_min/2.0
#endif
        do ii=-8,8
           xx=dble(ii)*dx_min/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=rmax)ncloud=ncloud+1
        end do
#if NDIM>1
     end do
#endif
#if NDIM==3
  end do
#endif

 ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
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
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in accrete_bondi'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
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
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif
        
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Remove mass from hydro cells
  do j=1,np
     if(ok(j))then
        isink=-idp(ind_part(j))
        r2=0d0
        do idim=1,ndim
           r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
        end do
        weight=exp(-r2/r2k(isink))
           
        d=uold(indp(j),1)
        u=uold(indp(j),2)/d
#if NDIM>1
        v=uold(indp(j),3)/d
#endif
#if NDIM>2
        w=uold(indp(j),4)/d
#endif
        e=uold(indp(j),ndim+2)/d
        do ivar=imetal,nvar
           z(ivar)=uold(indp(j),ivar)/d
        end do
         
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        d_ini=unew(indp(j),1)
        acc_mass=min(dMBHoverdt(isink),dMEdoverdt(isink)) &
             & *weight/total_volume(isink)*dtnew(ilevel)
        dmsink=max( min(acc_mass, (d-floorB*d_ini)*vol_loc), 0d0)

        ! Add the accreted mass to the total accreted mass over
        ! a coarse time step
        if(sink_AGN)then
           dMBH_coarse_new(isink)=dMBH_coarse_new(isink) + &
                & dMBHoverdt(isink)*weight/total_volume(isink)*dtnew(ilevel)
           dMEd_coarse_new(isink)=dMEd_coarse_new(isink) + &
                & dMEdoverdt(isink)*weight/total_volume(isink)*dtnew(ilevel)
           dMsmbh_new     (isink)=dMsmbh_new     (isink) + dmsink
        endif

        msink_new(isink  )=msink_new(isink  )+dmsink
        vsink_new(isink,1)=vsink_new(isink,1)+dmsink*u
#if NDIM>1
        vsink_new(isink,2)=vsink_new(isink,2)+dmsink*v
#endif
#if NDIM>2
        vsink_new(isink,3)=vsink_new(isink,3)+dmsink*w
#endif

        vp(ind_part(j),1)=mp(ind_part(j))*vp(ind_part(j),1)+dmsink*u
        vp(ind_part(j),2)=mp(ind_part(j))*vp(ind_part(j),2)+dmsink*v
        vp(ind_part(j),3)=mp(ind_part(j))*vp(ind_part(j),3)+dmsink*w
        mp(ind_part(j))=mp(ind_part(j))+dmsink
        vp(ind_part(j),1)=vp(ind_part(j),1)/mp(ind_part(j))
        vp(ind_part(j),2)=vp(ind_part(j),2)/mp(ind_part(j))
        vp(ind_part(j),3)=vp(ind_part(j),3)/mp(ind_part(j))

        d=d-dmsink/vol_loc
#ifdef SOLVERmhd
        e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        uold(indp(j),1)=d  
        uold(indp(j),2)=d*u
#if NDIM>1
        uold(indp(j),3)=d*v
#endif
#if NDIM>2
        uold(indp(j),4)=d*w
#endif
        uold(indp(j),ndim+2)=d*e
        do ivar=imetal,nvar
           uold(indp(j),ivar)=d*z(ivar)
        end do
     endif

  end do

end subroutine accrete_bondi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_jeans(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine determines if a cell covered by an old sink particle
  ! cross the density threshold for sink particle formation.If so, the 
  ! density is reduced and the corresponding mass is given to the sink.
  ! On exit, sink mass and velocity are modified.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,idim,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::r2,density,volume
  integer::nlevelmax_loc
  real(dp)::dx_min,vol_min,dx,dx_temp
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale
  real(dp)::d0,mstar,nISM,nCOM

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Reset new sink variables
  msink_new=0d0; vsink_new=0d0

  ! Finest cell size
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim
  ! Typical ISM mass density from H/cc to code units
  nISM = n_star
  nCOM = del_star*omega_b*rhoc/aexp**3*XH/mH
  nISM = MAX(nCOM,nISM)
  d0   = nISM/scale_nH
  ! Star particle mass
  mstar=MAX(del_star*omega_b*rhoc*XH/mH,n_star)/(scale_nH*aexp**3)*vol_min
  do i=1,nlevelmax
     dx_temp=scale*0.5D0**i
     ! Test is designed so that nlevelmax is activated at aexp \simeq 0.8
     if(d0*(dx_temp/2.0)**ndim.ge.mstar/2d0)nlevelmax_loc=i+1
  enddo

  ! Ensure that sinks accrete according to jeans only from the most refined cells
  if(ilevel .lt. nlevelmax_loc)return

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call accrete_jeans(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call accrete_jeans(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     msink_all=msink_new
     vsink_all=vsink_new
#endif
  endif
  do isink=1,nsink
     vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)+vsink_all(isink,1:ndim)
     msink(isink)       =msink(isink)       +msink_all(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)
  end do

111 format('   Entering grow_jeans for level ',I2)

end subroutine grow_jeans
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_jeans(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH, twopi 
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink,ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::xxx,mmm,r2,v2,c2,d,u,v=0d0,w=0d0,e,bx1,bx2,by1,by2,bz1,bz2
  real(dp),dimension(1:nvar)::z
  real(dp)::dx,dx_loc,scale,vol_loc,temp,d_jeans,acc_mass,d_sink,d_thres,nCOM
  real(dp)::factG,pi
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

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

  ! Density threshold for sink particle formation
  d_sink=n_sink/scale_nH

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
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
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in accrete_jeans'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
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
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif
        
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Gather hydro variables
  do j=1,np
     if(ok(j))then

        ! Convert to primitive variables
        d=uold(indp(j),1)
        u=uold(indp(j),2)/d
#if NDIM>1
        v=uold(indp(j),3)/d
#endif
#if NDIM>2
        w=uold(indp(j),4)/d
#endif
        e=uold(indp(j),ndim+2)/d
        do ivar=imetal,nvar
           z(ivar)=uold(indp(j),ivar)/d
        end do
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2

        ! Jeans length related density threshold
        temp=max(e*(gamma-1.0),smallc**2)
        d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
        d_thres=d_jeans

        ! User defined density threshold 
        !d_thres=d_sink

        if(d .ge. d_thres)then
           isink=-idp(ind_part(j))
           acc_mass=min((d-d_thres/4.0)*vol_loc,1d5*2d33/scale_m)
           msink_new(isink  )=msink_new(isink  )+acc_mass
           vsink_new(isink,1)=vsink_new(isink,1)+acc_mass*u
#if NDIM>1
           vsink_new(isink,2)=vsink_new(isink,2)+acc_mass*v
#endif
#if NDIM>2
           vsink_new(isink,3)=vsink_new(isink,3)+acc_mass*w
#endif           
           !d=d_thres/4.0
           d=d-acc_mass

           mp(ind_part(j))=mp(ind_part(j))+acc_mass

           ! Convert back to conservative variable
#ifdef SOLVERmhd
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(indp(j),1)=d
           uold(indp(j),2)=d*u
#if NDIM>1
           uold(indp(j),3)=d*v
#endif
#if NDIM>2
           uold(indp(j),4)=d*w
#endif
           uold(indp(j),ndim+2)=d*e
           do ivar=imetal,nvar
              uold(indp(j),ivar)=d*z(ivar)
           end do
        endif
          
     endif
  end do

#endif
  
end subroutine accrete_jeans
!################################################################
!################################################################
!################################################################
!################################################################
subroutine jet_AGN(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine computes the angular momentum vector of the gas around
  ! the sink 
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::xxx,mmm,r2,v2,c2,d,u,v,w,e,bx1,bx2,by1,by2,bz1,bz2,z
  real(dp)::dx,dx_loc,scale,vol_loc,weight,acc_mass,temp
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::meff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::d_x,d_y,d_z,uacc,vacc,wacc,ekin,dx_min,vol_min,rmax
  real(dp)::uc,vc,wc,r_close,dr,j_sp,E_sp
  integer ::ii,jj,kk

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

 ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
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
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in jet_AGN'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
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
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif
        
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Remove mass from hydro cells
  do j=1,np
     if(ok(j))then
        isink=-idp(ind_part(j))
        r2=0d0
        do idim=1,ndim
           r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
        end do
        weight=exp(-r2/r2k(isink))
           
        d=uold(indp(j),1)
        u=uold(indp(j),2)/d
#if NDIM>1
        v=uold(indp(j),3)/d
#endif
#if NDIM>2
        w=uold(indp(j),4)/d
#endif
        
        d_x=xp(ind_part(j),1)-xsink(isink,1)
        d_y=xp(ind_part(j),2)-xsink(isink,2)
        d_z=xp(ind_part(j),3)-xsink(isink,3)
        dr =sqrt(d_x*d_x+d_y*d_y+d_z*d_z)
        uc =u-vsink(isink,1)
        vc =v-vsink(isink,2)
        wc =w-vsink(isink,3)

        j_sp= d_x*uc + d_y*vc + d_z*wc
        jsink_new(isink,1)=jsink_new(isink,1) + (d_y*wc-d_z*vc)*d*vol_loc
        jsink_new(isink,2)=jsink_new(isink,2) + (d_z*uc-d_x*wc)*d*vol_loc
        jsink_new(isink,3)=jsink_new(isink,3) + (d_x*vc-d_y*uc)*d*vol_loc

     endif
  enddo

  do j=1,np
     isink=-idp(ind_part(j))
  enddo

end subroutine jet_AGN
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine get_rho_star(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Particles that are not entirely in
  ! level ilevel contribute also to the level density field
  ! (boundary particles) using buffer grids.
  ! Array rho_star is stored with:
  ! - rho_star containing the poisson source term fo stars
  !------------------------------------------------------------------
  integer::iskip,icpu,ind,i,info,nx_loc,ibound,idim
  real(dp)::dx,d_scale,scale,dx_loc
  real(dp)::dx_min,vol_min
  real(kind=8)::total,total_all,total2,total2_all,tms
  real(kind=8),dimension(2)::totals_in,totals_out

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  !--------------------------
  ! Initialize rho_star to zero
  !--------------------------
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho_star(active(ilevel)%igrid(i)+iskip)=0.0D0
     end do
  end do
  
  !-------------------------------------------------------
  ! Initialize rho_star to zero in virtual boundaries
  !-------------------------------------------------------
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           rho_star(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
        end do
     end do
  end do

  !---------------------------------------------------------
  ! Compute star particle contribution to density field
  !---------------------------------------------------------
  ! Compute density due to current level particles
  call rhostar_from_current_level(ilevel)
  ! Update boudaries
  call make_virtual_reverse_dp(rho_star(1),ilevel)
  call make_virtual_fine_dp   (rho_star(1),ilevel)

  !----------------------------------------------------
  ! Reset rho_star in physical boundaries
  !----------------------------------------------------
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho_star(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
        end do
     end do
  end do

end subroutine get_rho_star
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rhostar_from_current_level(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu
  integer::i,ig,ip,npart1
  real(dp)::dx

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0
    
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0   
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           
           ! Loop over particles
           do jpart=1,npart1
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              if(ip==nvector)then
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ig
                       x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
                    end do
                 end do
                 do i=1,ig
                    ind_cell(i)=father(ind_grid(i))
                 end do
                 call cic_star(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=nextp(ipart)  ! Go to next particle
           end do
           ! End loop over particles
           
        end if

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ig
              x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
           end do
        end do
        do i=1,ig
           ind_cell(i)=father(ind_grid(i))
        end do
        call cic_star(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
     end if

  end do
  ! End loop over cpus

end subroutine rhostar_from_current_level
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_star(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: mass_sph
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::ttt=0d0
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

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

  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
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

  ! Gather particle mass
  do j=1,np
     mmm(j)=mp(ind_part(j))
  end do

  ! Gather particle birth epoch
  if(star)then
     do j=1,np
        ttt(j)=tp(ind_part(j))
     end do
  endif

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in cic_amr'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif
        
  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icg(j,idim)=ig(j,idim)-2*igg(j,idim)
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
     icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

  ! Compute parent cell adress
  do ind=1,twotondim
     do j=1,np
        indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
     end do
  end do

  ! Update mass density field
  do ind=1,twotondim
     do j=1,np
        ok(j)=igrid(j,ind)>0
     end do

     do j=1,np
        vol2(j)=mmm(j)*vol(j,ind)/vol_loc
     end do

     do j=1,np
        if(ok(j).and.ttt(j).ne.0d0)then
           rho_star(indp(j,ind))=rho_star(indp(j,ind))+vol2(j)
        end if
     end do
  end do

end subroutine cic_star
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine quenching(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine selects regions which are eligible for SMBH formation.
  ! It is based on a stellar density threshold and on a stellar velocity
  ! dispersion threshold.
  ! On exit, flag2 array is set to 0 for AGN sites and to 1 otherwise.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::str_d,tot_m,ave_u,ave_v,ave_w,sig_u,sig_v,sig_w
  integer::igrid,jgrid,ipart,jpart,next_part,ind_cell,iskip,ind
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  real(dp),dimension(1:3)::skip_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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

#if NDIM==3
  ! Gather star particles only.

  ! Loop over grids
  do i=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(i)
     ! Number of particles in the grid
     npart1=numbp(igrid)
     npart2=0
     
     ! Reset velocity moments
     str_d=0.0
     tot_m=0.0
     ave_u=0.0
     ave_v=0.0
     ave_w=0.0
     sig_u=0.0
     sig_v=0.0
     sig_w=0.0
     
     ! Count star particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
              npart2=npart2+1
              tot_m=tot_m+mp(ipart)
              ave_u=ave_u+mp(ipart)*vp(ipart,1)
              ave_v=ave_v+mp(ipart)*vp(ipart,2)
              ave_w=ave_w+mp(ipart)*vp(ipart,3)
              sig_u=sig_u+mp(ipart)*vp(ipart,1)**2
              sig_v=sig_v+mp(ipart)*vp(ipart,2)**2
              sig_w=sig_w+mp(ipart)*vp(ipart,3)**2
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif
     
     ! Normalize velocity moments
     if(npart2.gt.0)then
        ave_u=ave_u/tot_m
        ave_v=ave_v/tot_m
        ave_w=ave_w/tot_m
        sig_u=sqrt(sig_u/tot_m-ave_u**2)*scale_v/1d5
        sig_v=sqrt(sig_v/tot_m-ave_v**2)*scale_v/1d5
        sig_w=sqrt(sig_w/tot_m-ave_w**2)*scale_v/1d5
        str_d=tot_m/(2**ndim*vol_loc)*scale_nH
     endif
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        ! AGN formation sites
        ! if n_star>0.1 H/cc and v_disp>100 km/s
        if(str_d>0.1.and.MAX(sig_u,sig_v,sig_w)>100.)then
           flag2(ind_cell)=0
        else
           flag2(ind_cell)=1
        end if
     end do
  end do
  ! End loop over grids

#endif

111 format('   Entering quenching for level ',I2)

end subroutine quenching
!################################################################
!################################################################
!################################################################
!################################################################
subroutine update_sink_position_velocity
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine updates position and velocity of sink particles.
  !------------------------------------------------------------------------
  integer::isink,idim,size_mpi,info,nx_loc
  real(dp)::vdum,ncloud,scale
  real(dp),dimension(1:3)::xbound

  ! Mesh spacing in that level
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
 
  ! update the sink particle position based on total sink particles
#ifndef WITHOUTMPI
  size_mpi=nsinkmax*(nlevelmax-levelmin+1)*(ndim*2+1)
  call MPI_ALLREDUCE(sink_stat,sink_stat_all,size_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  sink_stat_all=sink_stat
#endif
  do isink=1,nsink
     ncloud = sum(sink_stat_all(isink,levelmin:nlevelmax,ndim*2+1))
     if(ncloud>1)then
        do idim=1,ndim
           vdum = sum(sink_stat_all(isink,levelmin:nlevelmax,idim))
           xsink(isink,idim)=vdum/ncloud
           if(xsink(isink,idim)>scale*xbound(idim))then
              xsink(isink,idim)=xsink(isink,idim)-scale*xbound(idim)
           endif
           if(xsink(isink,idim)<0.0)then
              xsink(isink,idim)=xsink(isink,idim)+scale*xbound(idim)
           endif
        enddo
        do idim=1,ndim
           vdum = sum(sink_stat_all(isink,levelmin:nlevelmax,idim+ndim))
           vsink(isink,idim)=vdum/ncloud
        enddo
     endif
  enddo

end subroutine update_sink_position_velocity
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine AGN_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
#ifdef RT
  use rt_parameters,only: nRTvar
  use rt_hydro_commons,only: rtUold
  use rt_cooling_module,only: rt_AGN
#endif
  implicit none
  !----------------------------------------------------------------------
  ! Description: This subroutine checks AGN events in cells where a
  ! sink particle lies.
  ! Yohan Dubois, December 15th, 2010
  !----------------------------------------------------------------------
  ! local constants
  integer::nAGN,iAGN,ilevel,ivar,info
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical ,dimension(:),allocatable::ok_blast_agn
  integer ,dimension(:),allocatable::ind_blast,iAGN_myid
  real(dp),dimension(:),allocatable::mAGN,ZAGN,vol_gas,mass_gas,vol_blast,mass_blast,psy_norm
  real(dp),dimension(:),allocatable::dMBH_AGN,dMEd_AGN,dMsmbh_AGN,EsaveAGN,Msmbh
  real(dp),dimension(:,:),allocatable::xAGN,vAGN,jAGN
  real(dp)::temp_blast
  integer,dimension(1:nsink)::itemp
  integer::isort,isink,idim,ilun
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,Mfrac
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::filename,filedir,filecmd
  real(dp),allocatable,dimension(:)::xdp

  if(.not. hydro)return
  if(ndim.ne.3)return
  if(nsink.eq.0)return

  if(verbose)write(*,*)'Entering make_AGN'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(.not.rt_AGN)then
     do isink=1,nsink
        dMsmbh_coarse(isink)=dMsmbh(isink)
     enddo
  endif

  if(myid==1.and.nsink>0.and.sinkprops)then
     if(.not.(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0))then
     call title(nstep_coarse,nchar)
     filename='sink_'//TRIM(nchar)//'.dat'
     ilun=ncpu*4+11
     open(unit=ilun,file=TRIM(filename),form='unformatted')
     write(ilun)nsink     ! Number of sink
     write(ilun)ndim      ! Number of dimensions
     write(ilun)aexp      ! expansion factor
     write(ilun)scale_l   ! length scale
     write(ilun)scale_d   ! density scale
     write(ilun)scale_t   ! time scale
     allocate(xdp(1:nsink))
     write(ilun)idsink(1:nsink) ! Identities
     write(ilun)msink (1:nsink) ! Masses
     do idim=1,ndim
        xdp(1:nsink)=xsink(1:nsink,idim)
        write(ilun)xdp    ! Positions
     enddo
     do idim=1,ndim
        xdp(1:nsink)=vsink(1:nsink,idim)
        write(ilun)xdp    ! Velocities
     enddo
     do idim=1,ndim
        xdp(1:nsink)=jsink(1:nsink,idim)
        write(ilun)xdp    ! gas AM
     enddo
     write(ilun)dMBHoverdt(1:nsink) ! Bondi accretion rate
     write(ilun)dMEdoverdt(1:nsink) ! Eddington accretion rate
     write(ilun)dMsmbh_coarse(1:nsink) ! Total accreted mass
     write(ilun)d_avgptr  (1:nsink) ! Mean gas density
     write(ilun)c_avgptr  (1:nsink) ! Mean sound speed
     write(ilun)v_avgptr  (1:nsink) ! Relative BH-gas velocity
     write(ilun)Esave     (1:nsink) ! Energy saved from last coarse step
     write(ilun)t         ! Simulation time
     close(ilun)
     deallocate(xdp)
     endif
  endif

  ! Get only AGN that are in my CPU or close to the border of my cpu (huge speed-up!)
  call getAGNonmyid(itemp,nAGN)

  ! Allocate the arrays for the position and the mass of the AGN in myid
  allocate(xAGN(1:nAGN,1:3),vAGN(1:nAGN,1:3),mAGN(1:nAGN),ZAGN(1:nAGN),iAGN_myid(1:nAGN) &
       & ,dMBH_AGN(1:nAGN),dMEd_AGN(1:nAGN),dMsmbh_AGN(1:nAGN),Msmbh(1:nAGN),EsaveAGN(1:nAGN),jAGN(1:nAGN,1:3))
  xAGN=0d0; vAGN=0d0; mAGN=0d0; ZAGN=0d0; iAGN_myid=0; dMBH_AGN=0d0; dMEd_AGN=0d0; EsaveAGN=0d0; jAGN=0d0
  
!!$  if(myid==1)then
!!$     do isink=1,nsink
!!$        write(*,*)'================Before AGN input========================'
!!$        write(*,*)'m_BH =',msink(isink)*scale_d*scale_l**3d0/2d33,' M_sun'
!!$        write(*,*)'dm_bh /dt=',dMBHoverdt(isink)*scale_d*scale_l**3d0/2d33/scale_t*3600d0*24d0*365d0,' M_sun/yr'
!!$        write(*,*)'dm_edd/dt=',dMEdoverdt(isink)*scale_d*scale_l**3d0/2d33/scale_t*3600d0*24d0*365d0,' M_sun/yr'
!!$        write(*,*)'Esave    =',Esave(isink)*scale_d*scale_l**5/scale_t**2,' erg'
!!$        write(*,*)'========================================================'
!!$     enddo
!!$  endif

  do iAGN=1,nAGN
     isort=itemp(iAGN)
     iAGN_myid(iAGN)=isort
     xAGN(iAGN,1)=xsink(isort,1)
     xAGN(iAGN,2)=xsink(isort,2)
     xAGN(iAGN,3)=xsink(isort,3)
     vAGN(iAGN,1)=vsink(isort,1)
     vAGN(iAGN,2)=vsink(isort,2)
     vAGN(iAGN,3)=vsink(isort,3)
     jAGN(iAGN,1)=jsink(isort,1)
     jAGN(iAGN,2)=jsink(isort,2)
     jAGN(iAGN,3)=jsink(isort,3)
     Msmbh     (iAGN)=msink(isort)
     dMBH_AGN  (iAGN)=dMBH_coarse(isort)
     dMEd_AGN  (iAGN)=dMEd_coarse(isort)
     dMsmbh_AGN(iAGN)=dMsmbh_coarse(isort)
     EsaveAGN  (iAGN)=Esave (isort)
  enddo
     
  ! Allocate arrays that are outputs of average_AGN (initialised in average_AGN)
  allocate(vol_gas(1:nAGN),mass_gas(1:nAGN),psy_norm(1:nAGN),vol_blast(1:nAGN),mass_blast(1:nAGN))
  allocate(ind_blast(1:nAGN),ok_blast_agn(1:nAGN))
  
  ! Check if AGN goes into jet mode
  ok_blast_agn(1:nAGN)=.false.
  do iAGN=1,nAGN
     if(dMBH_AGN(iAGN)/dMEd_AGN(iAGN).lt.X_floor .and. EsaveAGN(iAGN).eq.0d0)then
        Mfrac=dMsmbh_AGN(iAGN)/(Msmbh(iAGN)-dMsmbh_AGN(iAGN))
        if(Mfrac.ge.jetfrac)ok_blast_agn(iAGN)=.true.
     endif
  enddo

  ! Compute some averaged quantities before doing the AGN energy input
  call average_AGN(xAGN,dMBH_AGN,dMEd_AGN,mAGN,ZAGN,jAGN,vol_gas,mass_gas,psy_norm,vol_blast &
       & ,mass_blast,ind_blast,nAGN,iAGN_myid,ok_blast_agn,EsaveAGN)
     
  ! Check if AGN goes into thermal blast wave mode
  do iAGN=1,nAGN
     if(dMBH_AGN(iAGN)/dMEd_AGN(iAGN) .ge. X_floor .and. EsaveAGN(iAGN).eq.0d0)then
        ! Compute estimated average temperature in the blast
        temp_blast=0.0
        if(vol_gas(iAGN)>0.0)then
           temp_blast=eAGN_T*1d12*dMsmbh_AGN(iAGN)/mass_gas(iAGN)
        else
           if(ind_blast(iAGN)>0)then
              temp_blast=eAGN_T*1d12*dMsmbh_AGN(iAGN)/mass_blast(iAGN)
           endif
        endif
        if(temp_blast>=TAGN)then
           ok_blast_agn(iAGN)=.true.
        endif
     endif
  end do
  
  ! Modify hydro quantities to account for a Sedov blast wave
  call AGN_blast(xAGN,vAGN,dMsmbh_AGN,dMBH_AGN,dMEd_AGN,mAGN,ZAGN,jAGN,ind_blast,vol_gas &
       & ,psy_norm,vol_blast,nAGN,iAGN_myid,ok_blast_agn,EsaveAGN)
  
  ! Reset total accreted mass if AGN input has been done
  do iAGN=1,nAGN
     if(ok_blast_agn(iAGN))then
        isort=iAGN_myid(iAGN)
        dMsmbh_coarse(isort)=0d0
     endif
  end do
  ! Important: initialise coarse Bondi and Eddington mass for the next coarse step
#ifndef WITHOUTMPI
  dMsmbh_new=dMsmbh_coarse
  call MPI_ALLREDUCE(dMsmbh_new,dMsmbh_all,nsink,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  dMsmbh_coarse=dMsmbh_all  
#endif
  if(.not. rt_AGN)dMsmbh=dMsmbh_coarse
  dMBH_coarse=0d0; dMEd_coarse=0d0
  ! Important: save the energy for the next time step that has not been put in the current time step
#ifndef WITHOUTMPI
  Esave_new=0d0
  do iAGN=1,nAGN
     isink=iAGN_myid(iAGN)
     Esave_new(isink)=EsaveAGN(iAGN)
  enddo
  call MPI_ALLREDUCE(Esave_new,Esave_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  Esave= Esave_all
#else
  Esave=EsaveAGN
#endif

!!$  if(myid==1)then
!!$     do isink=1,nsink
!!$        write(*,*)'================After AGN input========================='
!!$        write(*,*)'Esave    =',Esave(isink)*scale_d*scale_l**5/scale_t**2,' erg'
!!$        write(*,*)'========================================================'
!!$     enddo
!!$  endif

  ! Deallocate everything
  deallocate(vol_gas,mass_gas,psy_norm,vol_blast,mass_blast)
  deallocate(ind_blast)
  deallocate(xAGN,vAGN,mAGN,ZAGN,jAGN)
  deallocate(iAGN_myid,ok_blast_agn,dMBH_AGN,dMEd_AGN,dMsmbh_AGN,Msmbh,EsaveAGN)
  
  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo
  
#ifdef RT
  do ilevel=nlevelmax,levelmin,-1
     call rt_upload_fine(ilevel)
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     enddo
  enddo
#endif

end subroutine AGN_feedback
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine average_AGN(xAGN,dMBH_AGN,dMEd_AGN,mAGN,ZAGN,jAGN,vol_gas,mass_gas,psy_norm,vol_blast &
     & ,mass_blast,ind_blast,nAGN,iAGN_myid,ok_blast_agn,EsaveAGN)
  use pm_commons
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the AGN bubble
  ! Jet case    : remove some gas from the BH's cell
  ! Thermal case: get the mass of gas in the AGN bubble
  ! Yohan Dubois, December 15th, 2010
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nAGN,j,iAGN,ind,ix,iy,iz,ngrid,iskip
  integer::i,isink,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_AGN,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::x_half,y_half,z_half,x_box,y_box,z_box
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nAGN)::ind_blast,iAGN_myid
  real(dp),dimension(1:nAGN)::vol_gas,mass_gas,vol_blast,mass_blast,mAGN
  real(dp),dimension(1:nAGN)::dMBH_AGN,dMEd_AGN,ZAGN,psy_norm,EsaveAGN,X_radio
  real(dp),dimension(1:nAGN)::OxAGN,FeAGN
  real(dp),dimension(1:nAGN,1:3)::xAGN,vAGN,jAGN
#ifndef WITHOUTMPI
  real(dp),dimension(1:nsink)::vol_gas_mpi,mass_gas_mpi,mAGN_mpi,psy_norm_mpi,ZAGN_mpi
  real(dp),dimension(1:nsink)::vol_gas_all,mass_gas_all,mAGN_all,psy_norm_all,ZAGN_all
#endif
  logical ,dimension(1:nAGN)::ok_blast_agn
  logical ,dimension(1:nvector),save::ok
  real(dp)::jtot,j_x,j_y,j_z,drjet,dzjet,psy
  real(dp)::eint,ekk

  if(verbose)write(*,*)'Entering average_AGN'

  ! Mesh spacing in that level
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
  x_half=scale*xbound(1)/2.0; y_half=scale*xbound(2)/2.0; z_half=scale*xbound(3)/2.0
  x_box =scale*xbound(1); y_box =scale*xbound(2); z_box =scale*xbound(3)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do iAGN=1,nAGN
     X_radio(iAGN)=dMBH_AGN(iAGN)/dMEd_AGN(iAGN)
  enddo

  ! Maximum radius of the ejecta
  rmax=MAX(1d0*dx_min*scale_l/aexp,rAGN*3.08d21)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0d0;mass_gas=0d0;vol_blast=0d0;mass_blast=0d0;ind_blast=-1;psy_norm=0d0;ZAGN=0d0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
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

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iAGN=1,nAGN
                    dxx=x-xAGN(iAGN,1)
                    dyy=y-xAGN(iAGN,2)
                    dzz=z-xAGN(iAGN,3)
                    dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz

                    ! ------------------------------------------
                    ! case 0: Some energy has not been released 
                    ! ------------------------------------------
                    if(EsaveAGN(iAGN).gt.0d0)then
                       if(dr_AGN.le.rmax2)then
                          vol_gas (iAGN)=vol_gas (iAGN)+vol_loc
                       endif
                       dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                       if(dr_cell.le.dx_loc/2.0)then
                          ind_blast (iAGN)=ind_cell(i)
                          vol_blast (iAGN)=vol_loc
                       endif

                    ! ------------------------------------------
                    ! case 1: All energy has been released in 
                    ! previous time step -> choose between kinetic
                    ! or thermal feedback
                    ! ------------------------------------------                    
                    else 

                       ! ------------------------------------------                    
                       ! Jet feedback
                       ! ------------------------------------------                    
                       if(X_radio(iAGN).lt.X_floor)then
                          
                          if(ok_blast_agn(iAGN))then
                             jtot=sqrt(jAGN(iAGN,1)**2+jAGN(iAGN,2)**2+jAGN(iAGN,3)**2)
                             if(jtot.gt.0d0)then
                                j_x=jAGN(iAGN,1)/jtot
                                j_y=jAGN(iAGN,2)/jtot
                                j_z=jAGN(iAGN,3)/jtot
                                dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                                dzjet= dxx*j_x + dyy*j_y + dzz*j_z
                                drjet=sqrt(dr_AGN-dzjet*dzjet)
                                ! Check if the cell lies within the AGN jet cylindre
                                if (drjet .le. rmax .and. abs(dzjet) .le. rmax)then
                                   vol_gas(iAGN)=vol_gas(iAGN)+vol_loc
                                   psy=exp(-drjet**2/2d0/rmax**2d0)
                                   psy_norm(iAGN)=psy_norm(iAGN)+psy
                                endif
                                if(dr_cell.le.dx_loc/2.0)then
                                   ind_blast(iAGN)=ind_cell(i)
                                   d=uold(ind_blast(iAGN),1)
                                   u=uold(ind_blast(iAGN),2)/d
                                   v=uold(ind_blast(iAGN),3)/d
                                   w=uold(ind_blast(iAGN),4)/d
                                   ekk=0.5d0*d*(u*u+v*v+w*w)
                                   eint=uold(ind_blast(iAGN),5)-ekk  
                                   vol_blast  (iAGN)=vol_loc
                                   mAGN(iAGN)=min(mloadAGN*dMsmbh_coarse(iAGN),0.25d0*d*vol_loc)
                                   if(metal)then
                                      ZAGN(iAGN)=uold(ind_blast(iAGN),imetal)/d
                                      uold(ind_blast(iAGN),imetal)=uold(ind_blast(iAGN),imetal) &
                                           & - ZAGN(iAGN)*mAGN(iAGN)/vol_loc
                                   endif
                                   if(chemo)then
                                      OxAGN(iAGN)=uold(ind_blast(iAGN),imetal+1)/d
                                      FeAGN(iAGN)=uold(ind_blast(iAGN),imetal+2)/d
                                      uold(ind_blast(iAGN),imetal+1)=uold(ind_blast(iAGN),imetal+1) &
                                           & - OxAGN(iAGN)*mAGN(iAGN)/vol_loc
                                      uold(ind_blast(iAGN),imetal+2)=uold(ind_blast(iAGN),imetal+2) &
                                           & - FeAGN(iAGN)*mAGN(iAGN)/vol_loc
                                   endif

                                   d=uold(ind_blast(iAGN),1)-mAGN(iAGN)/vol_loc
                                   uold(ind_blast(iAGN),1)=d
                                   uold(ind_blast(iAGN),2)=d*u
                                   uold(ind_blast(iAGN),3)=d*v
                                   uold(ind_blast(iAGN),4)=d*w
                                   uold(ind_blast(iAGN),5)=eint+0.5d0*d*(u*u+v*v+w*w)
                                endif
                                ! If no spin for the jet then put all thermal
                             else
                                if(dr_AGN.le.rmax2) then
                                   vol_gas(iAGN)=vol_gas(iAGN)+vol_loc
                                endif
                                dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                             endif
                          endif
                          
                       ! ------------------------------------------                    
                       ! Thermal feedback   
                       ! ------------------------------------------                    
                       else
                          if(dr_AGN.le.rmax2)then
                             vol_gas (iAGN)=vol_gas (iAGN)+vol_loc
                             mass_gas(iAGN)=mass_gas(iAGN)+vol_loc*uold(ind_cell(i),1)
                          endif
                          dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                          if(dr_cell.le.dx_loc/2.0)then
                             ind_blast (iAGN)=ind_cell(i)
                             vol_blast (iAGN)=vol_loc
                             mass_blast(iAGN)=vol_loc*uold(ind_cell(i),1)
                          endif
                       endif
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

  !################################################################
#ifndef WITHOUTMPI
  vol_gas_mpi=0d0; mass_gas_mpi=0d0; mAGN_mpi=0d0; ZAGN_mpi=0d0; psy_norm_mpi=0d0
  ! Put the nAGN size arrays into nsink size arrays to synchronize processors
  do iAGN=1,nAGN
     isink=iAGN_myid(iAGN)
     vol_gas_mpi (isink)=vol_gas (iAGN)
     mass_gas_mpi(isink)=mass_gas(iAGN)
     mAGN_mpi    (isink)=mAGN    (iAGN)
     ZAGN_mpi    (isink)=ZAGN    (iAGN)
     psy_norm_mpi(isink)=psy_norm(iAGN)
  enddo
  call MPI_ALLREDUCE(vol_gas_mpi ,vol_gas_all ,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mass_gas_mpi,mass_gas_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mAGN_mpi    ,mAGN_all    ,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZAGN_mpi    ,ZAGN_all    ,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(psy_norm_mpi,psy_norm_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_mpi =vol_gas_all
  mass_gas_mpi=mass_gas_all
  mAGN_mpi    =mAGN_all
  ZAGN_mpi    =ZAGN_all
  psy_norm_mpi=psy_norm_all
  ! Put the nsink size arrays into nAGN size arrays
  do iAGN=1,nAGN
     isink=iAGN_myid(iAGN)
     vol_gas (iAGN)=vol_gas_mpi (isink)
     mass_gas(iAGN)=mass_gas_mpi(isink)
     mAGN    (iAGN)=mAGN_mpi    (isink)
     ZAGN    (iAGN)=ZAGN_mpi    (isink)
     psy_norm(iAGN)=psy_norm_mpi(isink)
  enddo
#endif
  !################################################################

  if(verbose)write(*,*)'Exiting average_AGN'

end subroutine average_AGN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine AGN_blast(xAGN,vAGN,dMsmbh_AGN,dMBH_AGN,dMEd_AGN,mAGN,ZAGN,jAGN,ind_blast,vol_gas &
     & ,psy_norm,vol_blast,nAGN,iAGN_myid,ok_blast_agn,EsaveAGN)
  use pm_commons
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH
#ifdef RT
  use rt_hydro_commons,only: rtUold
  use rt_cooling_module,only: rt_AGN, group_agnegy_frac
  use rt_parameters,only: nGroups, iGroups, group_egy
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine do the AGN energy inupt
  ! Jet case: do a kinetic jet solution as in Dubois et al. (2010)
  ! by depositing mass, momentum and energy into a small cylindre around the BH
  ! Thermal case: do a thermal energy input as in Teyssier et al. (2010)
  ! by depositing internal energy into small bubble
  ! Yohan Dubois, December 15th, 2010
  !------------------------------------------------------------------------
  integer::ilevel,j,iAGN,nAGN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info,ncache,ip,mygroup
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_AGN,d,u,v,w,ek,u_r,d_gas,dT
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,x_half,y_half,z_half,x_box,y_box,z_box
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::scale_Np,scale_Fp,scale_vol,Np_inj,scale_evtocode
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nAGN)::mAGN,ZAGN,p_gas,vol_gas,uBlast,vol_blast,OxAGN,FeAGN
  real(dp),dimension(1:nAGN)::EAGN,X_radio,psy_norm,dMsmbh_AGN,dMBH_AGN,dMEd_AGN,EsaveAGN,Eini
  real(dp),dimension(1:nAGN,1:3)::xAGN,vAGN,jAGN
  integer ,dimension(1:nAGN)::ind_blast,iAGN_myid
  logical ,dimension(1:nAGN)::ok_blast_agn,ok_save
  logical ,dimension(1:nvector),save::ok
  real(dp)::jtot,j_x,j_y,j_z,drjet,dzjet,psy,nCOM,T2_1,T2_2,ekk,eint,vm_,dr_cell,T2,etot
  real(dp)::eint1,ekkold
  integer::idim,isink
#ifndef WITHOUTMPI
  real(dp),dimension(1:nsink)::EsaveAGN_mpi,EsaveAGN_all
#endif

  if(verbose)write(*,*)'Entering AGN_blast'

  ! Mesh spacing in that level
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
  x_half=scale*xbound(1)/2.0; y_half=scale*xbound(2)/2.0; z_half=scale*xbound(3)/2.0
  x_box =scale*xbound(1); y_box =scale*xbound(2); z_box =scale*xbound(3)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif
  scale_vol=scale_l**ndim
  scale_m=scale_d*scale_l**3d0
  nCOM = 1d-1*omega_b*rhoc/aexp**3*XH/mH/scale_nH

  ! Maximum radius of the ejecta
  rmax=MAX(1d0*dx_min*scale_l/aexp,rAGN*3.08d21)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  do iAGN=1,nAGN
     X_radio(iAGN)=dMBH_AGN(iAGN)/dMEd_AGN(iAGN)
  enddo

  ! Energy of a 2keV photon into code units
  !E2keV=2d0*1d3*1.60217646d-12/(scale_d*scale_l**5/scale_t**2)
  scale_evtocode=1.60217646d-12/(scale_d*scale_l**5/scale_t**2)

  uBlast=0d0
  ok_save=.false.
  do iAGN=1,nAGN
     if(EsaveAGN(iAGN).gt.0d0)then
        ok_save(iAGN)=.true.
        EAGN (iAGN)=EsaveAGN(iAGN)
        p_gas(iAGN)=EAGN    (iAGN) / vol_gas(iAGN)
     else if(ok_blast_agn(iAGN))then
        if(X_radio(iAGN).lt.X_floor)then
           EAGN  (iAGN)=eAGN_K*0.1d0*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
           p_gas (iAGN)=(1d0-f_ekAGN)*EAGN(iAGN) / vol_gas(iAGN)
           if(mAGN(iAGN).gt.0d0)uBlast(iAGN)=sqrt(2*f_ekAGN*EAGN(iAGN)/mAGN(iAGN))
        else
           Eini  (iAGN)=0.1d0*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
           EAGN  (iAGN)=eAGN_T*0.1d0*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
           p_gas (iAGN)=EAGN(iAGN) / vol_gas(iAGN)
        endif
     endif
  end do
  EsaveAGN=0d0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
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

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale

                 do iAGN=1,nAGN

                    ! ------------------------------------------
                    ! case 0: Some energy has not been released 
                    ! in previous time step -> do thermal input
                    ! ------------------------------------------                    
                    if(ok_save(iAGN))then
                       dxx=x-xAGN(iAGN,1)
                       dyy=y-xAGN(iAGN,2)
                       dzz=z-xAGN(iAGN,3)
                       dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz
                       if(dr_AGN.le.rmax2)then
                          ekk=0d0
                          d=uold(ind_cell(i),1)
                          do idim=1,ndim
                             ekk=ekk+0.5d0*uold(ind_cell(i),idim+1)**2/d
                          end do
                          etot=uold(ind_cell(i),5)
                          eint=etot - ekk
                          T2_1=(gamma-1d0)*eint/d*scale_T2
                          if(T2_1 .lt. T2maxAGN)then
                             eint=eint + p_gas(iAGN)
                             T2_2=(gamma-1d0)*eint/d*scale_T2
                             if(T2_2 .le. T2maxAGN)then
                                uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_gas(iAGN)
                             else
                                uold(ind_cell(i),5)=uold(ind_cell(i),5)+T2maxAGN/scale_T2/(gamma-1d0)*d
                                EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGN)/scale_T2/(gamma-1d0)*d*vol_loc
                             endif
                          else
                             EsaveAGN(iAGN)=EsaveAGN(iAGN)+p_gas(iAGN) *vol_loc
                          endif
                       endif

                    ! ------------------------------------------
                    ! case 1: All energy has been released in 
                    ! previous time step -> choose between kinetic
                    ! or thermal feedback
                    ! ------------------------------------------                    
                    else if(ok_blast_agn(iAGN))then
                       ! ------------------------------------------
                       ! Jet case
                       ! ------------------------------------------
                       if(X_radio(iAGN).lt.X_floor)then
                          
                          dxx=x-xAGN(iAGN,1)
                          dyy=y-xAGN(iAGN,2)
                          dzz=z-xAGN(iAGN,3)
                          dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz
                          jtot=sqrt(jAGN(iAGN,1)**2+jAGN(iAGN,2)**2+jAGN(iAGN,3)**2)
                          if(jtot.gt.0d0)then
                             j_x=jAGN(iAGN,1)/jtot
                             j_y=jAGN(iAGN,2)/jtot
                             j_z=jAGN(iAGN,3)/jtot
                             dzjet= dxx*j_x + dyy*j_y + dzz*j_z
                             drjet=sqrt(dr_AGN-dzjet*dzjet)
                             
                             ! Check if the cell lies within the AGN jet cylindre
                             if (drjet .le. rmax .and. abs(dzjet) .le. rmax)then
                                psy=exp(-drjet**2/2d0/rmax**2d0) / psy_norm(iAGN)
                                if (dzjet.lt.0d0)then
                                   u=-j_x*uBlast(iAGN)
                                   v=-j_y*uBlast(iAGN)
                                   w=-j_z*uBlast(iAGN)
                                endif
                                if (dzjet.gt.0d0)then
                                   u= j_x*uBlast(iAGN)
                                   v= j_y*uBlast(iAGN)
                                   w= j_z*uBlast(iAGN)
                                endif
                                
                                ekk=0d0
                                ! Compute kinetic energy before velocity input
                                do idim=1,ndim
                                   ekk=ekk+0.5d0*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
                                end do
                                ekkold=ekk
                                ! Compute total energy before energy input
                                etot=uold(ind_cell(i),5)
                                ! Compute internal energy before energy input
                                eint=etot - ekk
                                ! Compute temperature T2=T/mu in Kelvin before energy input
                                T2_1=(gamma-1d0)*eint/uold(ind_cell(i),1)*scale_T2

                                d_gas=mAGN(iAGN)/vol_gas(iAGN)
                                ! Compute the density and the metal density of the cell
                                uold(ind_cell(i),1)=uold(ind_cell(i),1) + d_gas * psy
                                d=uold(ind_cell(i),1)
                                if(metal)then
                                   uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+ZAGN(iAGN)*d_gas*psy
                                   if(chemo)then
                                       uold(ind_cell(i),imetal+1)=uold(ind_cell(i),imetal+1)+OxAGN(iAGN)*d_gas*psy
                                       uold(ind_cell(i),imetal+2)=uold(ind_cell(i),imetal+2)+FeAGN(iAGN)*d_gas*psy
                                   endif
                                endif
                                ! Velocity at a given dr_AGN linearly interpolated between zero and uBlast
                                u= u + vAGN(iAGN,1)
                                v= v + vAGN(iAGN,2)
                                w= w + vAGN(iAGN,3)
                                
                                ! Add each momentum component of the jet to the gas
                                uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u * psy
                                uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v * psy
                                uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w * psy
                                
                                ekk=0d0
                                ! Compute kinetic energy after velocity input
                                do idim=1,ndim
                                   ekk=ekk+0.5*uold(ind_cell(i),idim+1)**2/d
                                end do
                                if(T2_1 .ge. T2maxAGN)then
                                   etot=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)*psy &
                                        & + p_gas(iAGN)
                                   ! Update total energy with new kinetic energy and 
                                   ! old internal energy (Temperature does not increase!)
                                   uold(ind_cell(i),5)=ekk+eint
                                   T2_2=(gamma-1d0)*(etot-uold(ind_cell(i),5))/d*scale_T2
                                   EsaveAGN(iAGN)=EsaveAGN(iAGN)+T2_2/scale_T2/(gamma-1d0)*d*vol_loc
                                else
                                   ! Compute new total energy
                                   etot=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)*psy &
                                        & + p_gas(iAGN)
                                   ! Compute new internal energy
                                   eint=etot - ekk
                                   ! Compute T2=T/mu in Kelvin
                                   T2_2=(gamma-1d0)*eint/d*scale_T2
                                   if(T2_2 .le. T2maxAGN)then
                                      uold(ind_cell(i),5)=ekk+T2_2/scale_T2/(gamma-1d0)*d
                                   else
                                      uold(ind_cell(i),5)=ekk+T2maxAGN/scale_T2/(gamma-1d0)*d
                                      EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGN)/scale_T2/(gamma-1d0)*d*vol_loc
                                   endif
                                endif
                               
                             endif
                             ! Jet case with jsink=0
                          else
                             if(dr_AGN.le.rmax2)then
                                uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_gas(iAGN)
                             endif
                          endif

                          ! ------------------------------------------
                          ! Thermal case
                          ! ------------------------------------------
                       else
                          dxx=x-xAGN(iAGN,1)
                          dyy=y-xAGN(iAGN,2)
                          dzz=z-xAGN(iAGN,3)
                          dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz
                          if(dr_AGN.le.rmax2)then
                             ekk=0d0
                             d=uold(ind_cell(i),1)
                             do idim=1,ndim
                                ekk=ekk+0.5d0*uold(ind_cell(i),idim+1)**2/d
                             end do
                             etot=uold(ind_cell(i),5)
                             eint=etot - ekk
                             T2_1=(gamma-1d0)*eint/d*scale_T2
                             if(T2_1 .lt. T2maxAGN)then
                                eint=eint + p_gas(iAGN)
                                T2_2=(gamma-1d0)*eint/d*scale_T2
                                if(T2_2 .le. T2maxAGN)then
                                   uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_gas(iAGN)
                                else
                                   uold(ind_cell(i),5)=uold(ind_cell(i),5)+T2maxAGN/scale_T2/(gamma-1d0)*d
                                   EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGN)/scale_T2/(gamma-1d0)*d*vol_loc
                                endif
                             else
                                EsaveAGN(iAGN)=EsaveAGN(iAGN)+p_gas(iAGN) *vol_loc
                             endif

                          endif
                       endif
                       
                    endif
                    !End of ok_blast_agn

                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels
  
  do iAGN=1,nAGN

     if(ind_blast(iAGN)>0)then

     ! ------------------------------------------
     ! case 0: Some energy has not been released 
     ! in previous time step -> do thermal input
     ! ------------------------------------------                    
     if(ok_save(iAGN))then
        if(vol_gas(iAGN)==0d0)then
           ekk=0d0
           d=uold(ind_blast(iAGN),1)
           do idim=1,ndim
              ekk=ekk+0.5d0*uold(ind_blast(iAGN),idim+1)**2/d
           end do
           etot=uold(ind_blast(iAGN),5)
           eint=etot - ekk
           T2_1=(gamma-1d0)*eint/d*scale_T2
           if(T2_1 .lt. T2maxAGN)then
              eint=eint + EAGN(iAGN)/vol_blast(iAGN)
              T2_2=(gamma-1d0)*eint/d*scale_T2
              if(T2_2 .le. T2maxAGN)then
                 uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+EAGN(iAGN)/vol_blast(iAGN)
              else
                 uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+T2maxAGN/scale_T2/(gamma-1d0)*d
                 EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGN)/scale_T2/(gamma-1d0)*d*vol_loc
              endif
           else
              EsaveAGN(iAGN)=EsaveAGN(iAGN)+EAGN(iAGN)
           endif
        endif
     ! ------------------------------------------
     ! case 1: All energy has been released in 
     ! previous time step -> choose between kinetic
     ! or thermal feedback
     ! ------------------------------------------                    
     else
        if(vol_gas(iAGN)==0d0.and.ok_blast_agn(iAGN))then
           ! ------------------------------------------
           ! Jet case
           ! ------------------------------------------
           if(X_radio(iAGN).lt.X_floor)then
              ! Here vol_blast lies for the cell volume where the AGN sits in
              d_gas=mAGN(iAGN)/vol_blast(iAGN)
              u=vAGN(iAGN,1)
              v=vAGN(iAGN,2)
              w=vAGN(iAGN,3)
              uold(ind_blast(iAGN),1)=uold(ind_blast(iAGN),1)+d_gas
              if(metal)then
                 uold(ind_blast(iAGN),imetal)=uold(ind_blast(iAGN),imetal)+d_gas*ZAGN(iAGN)
                 if(chemo)then
                    uold(ind_blast(iAGN),imetal+1)=uold(ind_blast(iAGN),imetal+1)+d_gas*OxAGN(iAGN)
                    uold(ind_blast(iAGN),imetal+2)=uold(ind_blast(iAGN),imetal+2)+d_gas*FeAGN(iAGN)
                 endif
              endif
              ekk=0d0
              d=uold(ind_blast(iAGN),1)
              do idim=1,ndim
                 ekk=ekk+0.5d0*uold(ind_blast(iAGN),idim+1)**2/d
              end do
              etot=uold(ind_blast(iAGN),5)
              eint=etot - ekk
              T2_1=(gamma-1d0)*eint/d*scale_T2
              if(T2_1 .lt. T2maxAGN)then
                 eint=eint + EAGN(iAGN)/vol_blast(iAGN)
                 T2_2=(gamma-1d0)*eint/d*scale_T2
                 if(T2_2 .le. T2maxAGN)then
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+EAGN(iAGN)/vol_blast(iAGN)
                 else
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+T2maxAGN/scale_T2/(gamma-1d0)*d
                    EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGN)/scale_T2/(gamma-1d0)*d*vol_loc
                 endif
              else
                 EsaveAGN(iAGN)=EsaveAGN(iAGN)+EAGN(iAGN)
              endif

           ! ------------------------------------------
           ! Thermal case
           ! ------------------------------------------
           else
              ekk=0d0
              d=uold(ind_blast(iAGN),1)
              do idim=1,ndim
                 ekk=ekk+0.5d0*uold(ind_blast(iAGN),idim+1)**2/d
              end do
              etot=uold(ind_blast(iAGN),5)
              eint=etot - ekk
              T2_1=(gamma-1d0)*eint/d*scale_T2
              if(T2_1 .lt. T2maxAGN)then
                 eint=eint + EAGN(iAGN)/vol_blast(iAGN)
                 T2_2=(gamma-1d0)*eint/d*scale_T2
                 if(T2_2 .le. T2maxAGN)then
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+EAGN(iAGN)/vol_blast(iAGN)
                 else
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+T2maxAGN/scale_T2/(gamma-1d0)*d
                    EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGN)/scale_T2/(gamma-1d0)*d*vol_loc
                 endif
              else
                 EsaveAGN(iAGN)=EsaveAGN(iAGN)+EAGN(iAGN)
              endif
              
           endif
        endif
     endif

     endif
  end do ! end loop over iAGN

  if(verbose)write(*,*)'Exiting AGN_blast'

end subroutine AGN_blast
!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedbackAGN(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  real(dp)::scale,r2
  integer::igrid,jgrid,ipart,jpart,next_part,isink,idim,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  oksink_new=0d0

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0.and.tp(ipart).eq.0)then
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    npart2=npart2+1
                 endif
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0.and.tp(ipart).eq.0)then
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if
                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig
                 endif
              endif
              if(ip==nvector)then
                 call feedbkAGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call feedbkAGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     oksink_all=oksink_new
#endif
  endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        dMsmbh_coarse(isink)=dMsmbh_coarse(isink)+dMsmbh(isink)
        dMsmbh(isink)=0d0
     endif
  enddo

111 format('   Entering thermal_feedbackAGN for level ',I2)

end subroutine thermal_feedbackAGN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbkAGN(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons,only: rtUold
  use rt_cooling_module,only: rt_AGN, group_agnegy_frac
  use rt_parameters,only: nGroups, iGroups, group_egy
#endif
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine thermal_feedbackAGN.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink,mygroup
  real(dp)::xxx,dx,dx_loc,scale,vol_loc,birth_time,current_time,dert
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_Np,scale_Fp
  real(dp)::scale_vol,scale_evtocode,Np_inj
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

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

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif
  scale_vol=scale_l**ndim

  ! Energy of a 2keV photon into code units
  !E2keV=2d0*1d3*1.60217646d-12/(scale_d*scale_l**5/scale_t**2)
  scale_evtocode=1.60217646d-12/(scale_d*scale_l**5/scale_t**2)

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
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
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in thermalAGN'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
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
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Update hydro variables due to feedback
  do j=1,np
     if(ok(j))then
        isink=-idp(ind_part(j))
        dert=0.1d0*dMsmbh(isink)*(3d10/scale_v)**2
        do mygroup=1,nGroups
           Np_inj=dert * group_agnegy_frac(mygroup) / (scale_evtocode*group_egy(mygroup)) / (vol_loc*scale_vol) / scale_Np
           rtuold(indp(j),iGroups(mygroup))=rtuold(indp(j),iGroups(mygroup))+Np_inj
        enddo
        oksink_new(isink)=1d0
     endif
  end do

#endif
  
end subroutine feedbkAGN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getAGNonmyid(isink_myid,nsink_myid)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine check which BHs stand in the cpu myid
  ! Yohan Dubois, December 15th, 2010
  !------------------------------------------------------------------------
  integer,dimension(1:nsink)::isink_myid
  integer::nsink_myid,ii,info
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,isink,nx_loc,ilevel,lmax,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drsink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  
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

  ! Maximum radius of the ejecta
  drsink=2d0*MAX(1d0*dx_min*scale_l/aexp,rAGN*3.08d21)
  drsink=drsink/scale_l
  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  isink_myid=0
  ii=0
  do isink=1,nsink

     cpu_read=.false.
     ! Compute boundaries for the sink cube of influence
     xxmin=(xsink(isink,1)-drsink)/scale ; xxmax=(xsink(isink,1)+drsink)/scale
     yymin=(xsink(isink,2)-drsink)/scale ; yymax=(xsink(isink,2)+drsink)/scale
     zzmin=(xsink(isink,3)-drsink)/scale ; zzmax=(xsink(isink,3)+drsink)/scale

     if(TRIM(ordering).eq.'hilbert')then
        
        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif
        
        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax
        
        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min=0.0d0
           endif
           bounding_min(i)=(order_min)*dkey
           bounding_max(i)=(order_min+oneqdp)*dkey
        end do
        
        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do
        
        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if
     
     ! Create the index array for sinks in processor myid
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
           ii=ii+1
           isink_myid(ii)=isink
        endif
     enddo

  enddo

  ! Number of sinks in processor myid
  nsink_myid=ii

end subroutine getAGNonmyid
!################################################################
!################################################################
!################################################################
!################################################################
