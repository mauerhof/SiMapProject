!*************************************************************************
RECURSIVE SUBROUTINE rt_step(ilevel, dt)

! This routine is the adaptive-mesh main driver for radiative transfer.
! Updates RT variables and pressure through cooling on ilevel and all 
! levels below in the given timestep length.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  use pm_commons
  use SED_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount,ivar,i,idim
  real(dp):: dt
  logical,save::first_step=.true.
!-------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,999)ilevel

  if(rt_pp .and. rt_refine) then
     !-------------------------------------------
     ! Make new refinements and update boundaries
     !-------------------------------------------
     if(levelmin.lt.nlevelmax)then
        if(ilevel==levelmin.or.icount>1)then
           do i=ilevel,nlevelmax
              if(i>levelmin)then
                 
                 !--------------------------
                 ! Build communicators
                 !--------------------------
                 call build_comm(i)

                 !--------------------------
                 ! Update boundaries
                 !--------------------------
                 call make_virtual_fine_int(cpu_map(1),i)
                 do ivar=1,nvar
                    call make_virtual_fine_dp(uold(1,ivar),i)
                 end do
                 if(simple_boundary)call make_boundary_hydro(i)
              end if

              !--------------------------
              ! Refine grids
              !--------------------------
              call refine_fine(i)
           end do
        end if
     end if

     !--------------------------
     ! Load balance
     !--------------------------
     if(levelmin.lt.nlevelmax)then
        if(ilevel==levelmin)then
           if(nremap>0)then
              ! Skip first load balance because it has been 
              ! performed before file dump
              if(nrestart>0.and.first_step)then
                 first_step=.false.
              else
                 if(MOD(nstep_coarse,nremap)==0)then
                    call load_balance
                    call defrag
                 endif
              end if
           end if
        endif
     end if

  endif

  call rt_set_unew(ilevel)                        ! Set unew equal to uold

  !---------------------------
  ! Recursive call to rt_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        call rt_step(ilevel+1,dt/2.)
        call rt_step(ilevel+1,dt/2.)
     else if(rt_pp) then
        call update_time(ilevel)
     else if(rt_coupling_pp) then
        t=t+dt                ! We're only updating a temporary time here!
        call getProperTime(t, t_proper)
     end if
  else if(rt_pp) then
     call update_time(ilevel)
  else if(rt_coupling_pp) then
     t=t+dt                   ! We're only updating a temporary time here!
     call getProperTime(t, t_proper)
  endif

  !-----------
  ! rt step
  !-----------
  call rt_godunov_fine(ilevel, dt)             !         Hyperbolic solver
  if(rt_star) call star_RT_feedback(ilevel, dt)!              Updates unew
  call add_rt_sources(ilevel,dt)               !              Updates unew

  do ivar=1,nvar                               ! Reverse update boundaries
     call make_virtual_reverse_dp(unew(1,ivar),ilevel)
  end do

  call rt_set_uold(ilevel)                    !     Set uold equal to unew

  if(rt_pp .or. rt_coupling_pp) then          !      Post-processing: 
     call upload_fine(ilevel)                 !          Upload everything 
  else if (rt_coupling_full) then             !      Coupled: 
     call rt_upload_fine(ilevel)              !        Upload only rt-vars
  endif

  call cooling_fine(ilevel,dt)                !     Ionization and cooling

  call add_cont_sources(ilevel,dt)            !                     ->uold

  do ivar=1,nvar                              !          Update boundaries 
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
  if(simple_boundary)call make_boundary_hydro(ilevel)

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(rt_pp .and. rt_refine) call flag_fine(ilevel,1)

999 format(' Entering rt_step for level',i1)

END SUBROUTINE rt_step


!*************************************************************************
SUBROUTINE rt_pp_step

! Post-processing radiative transfer step, includes deciding the timestep
! length and determining whether or not to do output.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
  real(dp):: rt_dt
  integer::ilevel
  character(LEN=5)::nchar
!-------------------------------------------------------------------------
  if(mod(nstep_coarse,foutput)==0  .or. nstep_coarse .eq. 0 .or.         &
     aexp>=aout(iout) .or. t>=tout(iout) .or. output_iState)       then
     ! Output current state of affairs to files
     output_iState=.false.
     if(myid .eq. 1) then
        call title(ifout,nchar)
     endif
     call defrag ;  call dump_all
     if(myid .eq. 1) then
        write(*,*) 'Wrote current state to output_',nchar
     endif
  endif
  if(cosmo) call update_rt_c
  call get_rt_courant_coarse(rt_dt)
  do ilevel = levelmin, nlevelmax
     dtnew(ilevel)=rt_dt/(2.**(ilevel-levelmin))
     dtold(ilevel)=rt_dt/(2.**(ilevel-levelmin))
  enddo

  call rt_step(levelmin, rt_dt)
END SUBROUTINE rt_pp_step


!*************************************************************************
SUBROUTINE rt_coupled_pp_step(dt)

! Semi-coupled RT step. Completely separated from amr-step, and run only
! once a coarse hydro timestep is finished.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use pm_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp):: dt, rt_dt_courant, t_save
  real(dp),save:: dt_left=0.  ! timestep left over, can also be from 
                              ! previous rt-step
  integer::i,info, ilevel
!-------------------------------------------------------------------------
  call get_rt_courant_coarse(rt_dt_courant)
  dt_left = dt_left + dt
  ! we do a little trick with the time, to get evolution of star
  ! ages within the hydro timestep:
  t_save=t ; t=t-dt_left ! Shift the time backwards

  if(pic .and. rt_star) then ! Set up grid-to-star pointers for all levels
     do ilevel = levelmin, nlevelmax
        call kill_tree_fine(ilevel)
        call virtual_tree_fine(ilevel)
     enddo
  endif
  i=0
  do
     !if(myid==1) print*,'Time:',dt,dt_left,rt_dt_courant
     if(dt_left .lt. rt_dt_courant) then
        t=t_save                                        ! Restore the time
        exit                  ! Keep rest of timestep for the next rt_step
     endif
     call rt_step(levelmin, rt_dt_courant)
     dt_left = dt_left - rt_dt_courant
     i=i+1
  enddo

  if(pic .and. rt_star) then               ! Tidy up grid-to-star pointers
     do ilevel = levelmin, nlevelmax
        call merge_tree_fine(ilevel)
     enddo
  endif
  if(myid==1) write(*,999)i
999 format(' RT pp-coupling did ',I4,' substeps')

END SUBROUTINE rt_coupled_pp_step


!*************************************************************************
SUBROUTINE rt_coupled_step(dt)

! Coupled radiative transfer step for all levels, should
! take place within the usual AMR step at the finest level, and
! evolve the RT over all levels.
! dt: Timestep length, over which to evolve RT (should be the hydro
! levelmax - i.e. finest - timestep length).
! This is actually not implemented, and not immediately possible because
! of the way the AMR substepping is structured...the routine sits here
! ready though if we evetually come up with a way of doing RT substepping
! at the finest hydro level.
! Joki - Sept. 2011.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
  real(dp),intent(in):: dt
  real(dp):: rt_dt_courant
  integer:: i
  real(dp),save:: dt_left=0.  ! timestep left over, can also be from 
                              ! previous rt-step
!-------------------------------------------------------------------------
  call get_rt_courant_coarse(rt_dt_courant)
  dt_left = dt_left + dt
  do
     if(dt_left .lt. rt_dt_courant) then
        exit  !keep rest of timestep for the next rt_step
     endif
     call rt_step(levelmin, rt_dt_courant)
     dt_left = dt_left - rt_dt_courant
  enddo

  ! update the hydro timestep length at all levels, since the pressure
  ! and thereby the Courant condition has changed in all cells.
  !do i = levelmin, nlevelmax
  !   call newdt_fine(i)
  !   if(i>levelmin)then
  !      dtnew(i)=MIN(dtnew(i-1)/real(nsubcycle(i-1)),dtnew(i))
  !      if(nsubcycle(i-1)==1) dtnew(i-1) = dtnew(ilevel)
  !      if(icount==2) dtnew(i-1) = dtold(i) + dtnew(i)
  !end if

END SUBROUTINE rt_coupled_step


!*************************************************************************
SUBROUTINE get_rt_courant_coarse(dt)

! Determine the coarse RT timestep length set by the Courant condition
!-------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  implicit none
  integer:: nx_loc
  real(dp):: dt, scale, dx
!-------------------------------------------------------------------------
  ! Mesh spacing at coarse level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**levelmin*scale
  dt = rt_courant_factor*dx/3.d0/rt_c
END SUBROUTINE get_rt_courant_coarse

!*************************************************************************
SUBROUTINE output_rt_stats

! Output and reset rt statistics. These are cooling statistics and
! star rt feedback statistics
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  implicit none
  integer*8:: max_all, tot_all, cells_all,loopCodes_tot
  integer*8:: loopCodes_all(4)
  integer::info
  real(dp)::step_nPhot_all, step_nStar_all, step_mStar_all
  real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nh, scale_T2
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Cooling statistics:
  cells_all=0 ; tot_all=0 ; max_all=0 ; loopCodes_all=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(n_cool_cells,         cells_all,     1, &
                     MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(tot_cool_loopcnt,     tot_all,       1, &
                     MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(max_cool_loopcnt,     max_all,       1, &
                     MPI_INTEGER,          MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(loopCodes,            loopCodes_all, 4, &
                     MPI_INTEGER,          MPI_MAX, MPI_COMM_WORLD, info)
  n_cool_cells     = cells_all ; tot_cool_loopcnt = tot_all
  max_cool_loopcnt = max_all   ; loopCodes        = loopCodes_all
#endif
  if(myid .eq. 1) then
     if(n_cool_cells .eq. 0) n_cool_cells=1.
     write(*, 111) dble(tot_cool_loopcnt)/n_cool_cells,max_cool_loopcnt,rt
     loopCodes_tot = SUM(loopCodes)
     if(loopCodes_tot .gt. 0) then
        write(*, 112) dble(loopCodes)/dble(loopCodes_tot)
     else
        write(*, 112) dble(loopCodes)
     endif
  endif
  max_cool_loopcnt=0; tot_cool_loopcnt=0; n_cool_cells=0; loopCodes(:)=0
111 format(' Coolstats: Avg. # loops = ', f21.6, ', max. # loops = ', I10, ', rt=',L)
112 format(' Subcycling codes [Np, T, xH, xHe]% = ', 4(f7.3, ''))

  ! Stellar rt feedback statistics:
  if(showSEDstats .and. rt_star) then
     step_nPhot_all=0.d0 ; step_nStar_all=0.d0 ; ; step_mStar_all=0.d0
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(step_nPhot,           step_nPhot_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(step_nStar,           step_nStar_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(step_mStar,           step_mStar_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     step_nPhot  = step_nPhot_all
     step_nStar  = step_nStar_all
     step_mStar  = step_mStar_all
#endif
     tot_nPhot = tot_nPhot + step_nPhot
     if(myid .eq. 1)                                                     &
          write(*, 113) step_nPhot, tot_nPhot, step_nStar/dtnew(levelmin)&
          ,step_mStar/dtnew(levelmin), dtnew(levelmin)*scale_t/(3.15569d7)
     step_nPhot = 0.d0 ; step_nStar = 0.d0 ; step_mStar = 0.d0
  endif
113 format(' SED feedback(phot/step/1d50, phot/tot/1d50, *, */Msun , dt[yr])= '  &
                                                             ,10(1pe9.2))
END SUBROUTINE output_rt_stats
