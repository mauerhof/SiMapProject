!RT patch: Plenty of change here. To see it, do a diff.
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
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
