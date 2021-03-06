! Ifront_writer patch (see !IF comments): 
!        Added subroutine write_ifront to find ifront position in a
!        stromgren sphere/shadow experiment and write it to std output.
! RT patch: Plenty of change here. To see it, do a diff.
!*************************************************************************
subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use rt_parameters
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,idim,ivar,info
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem,real_mem_tot
  logical::doCool !----------------------------------------------------!RT
                  ! This is to initialize the RT cooling-module,
                  ! which is used even if there is no actual cooling,
                  ! since we still do ionization in that module.
                  !----------------------------------------------------!RT

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  doCool=cooling ; cooling=.true.    ! --------------------------------!RT
  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables    
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid
  call set_table(dble(aexp))         ! --------------------------------!RT
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again
  call rt_init                       ! Initialize RT variables---------!RT
  cooling=doCool                     ! --------------------------------!RT

#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse
  if(rt_pp) nstep_coarse=nstep_coarse+1 ! To output the initial state!-!RT

  if(myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if((.not. rt_pp .or. rt_refine) .and. levelmin.lt.nlevelmax)then!-!RT
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     !call write_ifront_stromgren !-------------------------------------!IF
     call write_ifront_shadow !----------------------------------------!IF

     ! Call base level
     !Begin !RT-----------------------------------------------------------
     if ( rt .and. rt_pp ) then                    
        call rt_pp_step                   !             Post-processing RT
     else                                 
        call amr_step(levelmin,1,.true.)  ! Possibly coupled RT (or no RT)
     endif
     if( rt .and. .not. rt_pp .and. rt_coupling_pp )                     &
        call rt_coupled_pp_step(dtnew(levelmin))         ! Semi-coupled RT
     
     call output_rt_stats
     !End !RT-------------------------------------------------------------
     if((.not. rt_pp .or. rt_refine) .and. levelmin.lt.nlevelmax)then!-!RT
        ! Hydro book-keeping
        if(hydro)then
           do ilevel=levelmin-1,1,-1
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end do
        end if
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX      &
             ,MPI_COMM_WORLD,info)
        if(myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop

!************************************************************************
SUBROUTINE write_ifront_stromgren

! Find and write to standard output the position of an ionization front
! from the x=y=z=0 origin of the box.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ind,info
  integer::i,j,ilevel,ncache,ix,iy,iz
  integer::nx_loc,nbins,j0
  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer,dimension(:),allocatable::ind_grid,ind_cell
  real(kind=8)::xion,rr
  real(kind=8),dimension(:,:),allocatable::bins,bins_all  !  r,r_avg,x_avg
  integer,dimension(:),allocatable::bin_count,bin_count_all
  real(dp)::ifront_rad
!-------------------------------------------------------------------------

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  rr=0.0D0; xion=0.0D0
  nbins=100
  allocate(bin_count(nbins), bin_count_all(nbins))  ; bin_count=0
  allocate(bins(3,nbins)   , bins_all(3,nbins))     ; bins=0
  do i=1,nbins !initialize the bins
     bins(1,i) = (i-1) * boxlen/(nbins) ! Lower radial values for bins
  end do

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 ! Now need to bin the cell's ionization state according 
                 ! to it's radius
                 rr = sqrt(   &
                      ((xg(ind_grid(i),1)+xc(ind,1)-dble(icoarse_min))   &
                       *scale)**2 &
                      +((xg(ind_grid(i),2)+xc(ind,2)-dble(jcoarse_min))  &
                       *scale)**2 &
                      +((xg(ind_grid(i),3)+xc(ind,3)-dble(kcoarse_min))  &
                       *scale)**2 )
                 xion = uold(ind_cell(i),iIons)/uold(ind_cell(i),1)
                 do j=1,nbins
                    if(bins(1,j) .ge. rr) then
                       bins(2,j) = bins(2,j) + rr
                       bins(3,j) = bins(3,j) + xion
                       bin_count(j)=bin_count(j)+1
                       exit
                    endif
                 end do

              end if
           end do
        end do
        deallocate(ind_grid, ind_cell)
     end if
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(bins,      bins_all,      nbins*3, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(bin_count, bin_count_all, nbins, &
                     MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
#endif

  if(myid==1) then
     bins(2:3,:)=bins_all(2:3,:); bin_count=bin_count_all
     bins(2,:)=bins(2,:)/bin_count ! Average radius per radius bin
     bins(3,:)=bins(3,:)/bin_count ! Average xion per radius bin
     ifront_rad=0.
     j0=0
     do j=1,nbins
        if(bin_count(j) .gt. 0) then
           if(bins(3,j) .le. 0.5) then
              exit ! Found the front. j points at first bin w x<=0.5
           endif
           j0=j ! Adjacent nonempty bin to the left for interpolation
        endif
     end do
     ! Interpolate
     if(j0 .eq. 0) then
        ifront_rad = 0.
     else
        ifront_rad = bins(2,j0) + (0.5d0-bins(3,j0))                     &
                         / (bins(3,j)-bins(3,j0)) * (bins(2,j)-bins(2,j0))
     endif
     
     write(*, 111) t, ifront_rad

  end if

  ! Deallocate local arrays
  deallocate(bin_count,bins)
  deallocate(bin_count_all,bins_all)

 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

111 format('IFRONT     ', f21.6, f21.6)

END SUBROUTINE write_ifront_stromgren

!*************************************************************************
SUBROUTINE write_ifront_shadow
! Find and write to standard output the position of an ionization front
! along the x-axis at the center of the yz-plane, i.e. at the center of 
! the box.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ind,info
  integer::i,j,k,ilevel,ncache,ix,iy,iz
  integer::nx_loc,nbins,j0
  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer,dimension(:),allocatable::ind_grid,ind_cell
  real(kind=8)::xHII, rx, ry, rz, dr_max
  real(kind=8),dimension(:,:),allocatable::bins,bins_all!x,x_avg,xHII_avg
  integer,dimension(:),allocatable::bin_count,bin_count_all
  real(dp)::ifront_x
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  rx=0.0D0; xHII=0.0D0
  nbins=200
  allocate(bin_count(nbins), bin_count_all(nbins))  ; bin_count=0
  allocate(bins(3,nbins)   , bins_all(3,nbins))     ; bins=0
  do i=1,nbins !initialize the bins
     bins(1,i) = (i-1) * boxlen/(nbins) ! Lower x values for bins
  end do

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  ! width of a cell
  !dr_max = sqrt(2.)*scale/2**(levelmin)                   ! Use this in Iliev3
  dr_max = sqrt(2.)*scale/2**(nlevelmax) ! Careful here  ! Use this in Iliev7
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then ! only bin leaf cells
                 ! bin the cell's ionization state according to it's x-pos
                 rx= (xg(ind_grid(i),1)+xc(ind,1)-dble(icoarse_min))*scale
                 ry= (xg(ind_grid(i),2)+xc(ind,2)-dble(jcoarse_min))*scale
                 rz= (xg(ind_grid(i),3)+xc(ind,3)-dble(kcoarse_min))*scale
                 if(sqrt((ry-boxlen/2)**2+(rz-boxlen/2)**2) .gt. dr_max) &
                      cycle ! Only consider cells along x in box middle
                 xHII = uold(ind_cell(i),iIons)/uold(ind_cell(i),1)
                 do j=1,nbins
                    if(bins(1,j) .ge. rx) then
                       bins(2,j) = bins(2,j) + rx
                       bins(3,j) = bins(3,j) + xHII
                       bin_count(j)=bin_count(j)+1
                       exit
                    endif
                 end do
              end if ! End leaf cells
           end do ! End ncache loop
        end do ! End twotondim loop
        deallocate(ind_grid, ind_cell)
     end if ! End ncache>0?
  end do ! End ilevel loop

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(bins,      bins_all,      nbins*3,                  &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(bin_count, bin_count_all, nbins,                    &
                     MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
#endif

  if(myid==1) then
     bins(2:3,:)=bins_all(2:3,:); bin_count=bin_count_all
     bins(2,:)=bins(2,:)/bin_count ! Average x per x bin
     bins(3,:)=bins(3,:)/bin_count ! Average xHII per x bin
     ifront_x=0.
     j0=0
     do j=1,nbins
        if(bin_count(j) .gt. 0) then
           if(bins(3,j) .le. 0.5) then
              exit ! Found the front. j points at first bin w xHII<=0.5
           endif
           j0=j ! Adjacent nonempty bin to the left for interpolation
        endif
     end do
     ! Interpolate
     if(j0 .eq. 0) then
        ifront_x = 0.
     else
        ifront_x = bins(2,j0) + (0.5d0-bins(3,j0))                       &
                         / (bins(3,j)-bins(3,j0)) * (bins(2,j)-bins(2,j0))
     endif
    
     write(*, 111) t, ifront_x
  end if

  ! Deallocate local arrays
  deallocate(bin_count,bins)
  deallocate(bin_count_all,bins_all)

 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

111 format('IFRONT     ', f21.6, f21.6)

END SUBROUTINE write_ifront_shadow
