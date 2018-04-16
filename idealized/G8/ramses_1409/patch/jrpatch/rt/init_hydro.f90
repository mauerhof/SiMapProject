!RT-- Added reading of RT variables from file at restart.
!     Also added reading of rt_infoXXXXX.txt to read the rt vars correctly.
!     Also added subroutine read_int to read the info parameter file.
!*************************************************************************

subroutine init_hydro
  use amr_commons
  use hydro_commons
  use rt_parameters !--------------------------------------------------!RT
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  integer::nhvar2=-1,nRTvar2=0,varCount
  logical::ok

  if(verbose)write(*,*)'Entering init_hydro'

  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0 ; unew=0.0d0
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     ! Try to read rt info file and retrieve nhvar2 and nrtvar2        !RT
     fileloc='output_'//TRIM(nchar)//'/info_rt_'//TRIM(nchar)//'.txt'  !RT
     inquire(file=fileloc, exist=ok)                                   !RT
     if (.not. ok) then                                                !RT
        nhvar2=-1 ; nRTvar2=0                                          !RT
     else                                                              !RT
        open(unit=ilun,file=fileloc,status='old',form='formatted')     !RT
        call read_int( ilun, 'nHvar', nhvar2)                          !RT
        call read_int( ilun, 'nRTvar', nRTvar2)                        !RT
        close(ilun)                                                    !RT
     endif                                                             !RT

     fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     if(nhvar2 .eq. -1) nhvar2=nvar2 !---------------------------------!RT
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(nhvar2.gt.nhvar .and. myid==1)then ! OK to drop variables -----!RT
        write(*,*)'File hydro.tmp is not compatible (1)'               !RT
        write(*,*)'Found nhvar  =',nhvar2                              !RT
        write(*,*)'Expected=',nhvar, ' or less'                        !RT
        write(*,*)'So reading the first nhvar and ignoring the rest'   !RT
     end if                                                            !RT
     if(nhvar2.lt.nhvar )then ! OK to add new vars --------------------!RT
        if(myid==1) write(*,*)'File hydro.tmp is not compatible (2)'   !RT
        if(myid==1) write(*,*)'Found nhvar2  =',nhvar2                 !RT
        if(myid==1) write(*,*)'Expected=',nhvar                        !RT
        if(myid==1) write(*,*)'..so only reading first ',nhvar2, &     !RT
                  'variables and setting the rest to zero'             !RT
        ! If we can't read ionization states -> initialize them from T !RT
        varCount=2+ndim                                                !RT
        if(metal) varCount=varCount+1                                  !RT
        if(delayed_cooling) varCount=varCount+1                        !RT
        if(aton) varCount=varCount+1                                   !RT
        if(nhvar2 .le. varCount) rt_is_init_xion = .true.              !RT
     end if                                                            !RT
     if(nrtvar2.lt.nrtvar .and. myid==1)then ! OK to add new pacs -----!RT
        write(*,*)'File hydro.tmp is not compatible (1)'               !RT
        write(*,*)'Found nrtvar  =',nrtvar2                            !RT
        write(*,*)'Expected=',nrtvar                                   !RT
        write(*,*)'..so only reading first ',nrtvar2, &                !RT
                  'rt variables and setting the rest to zero'          !RT
     end if                                                            !RT
     if(nrtvar2.gt.nrtvar)then ! NOT OK to remove pacs ----------------!RT
        if(myid==1) write(*,*)'File hydro.tmp is not compatible (1)'   !RT
        if(myid==1) write(*,*)'Found nrtvar  =',nrtvar2                !RT
        if(myid==1) write(*,*)'Expected=',nrtvar,' so QUITTING'        !RT
        call clean_stop
     end if
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 ! Loop over conservative variables
                 !do ivar=1,min(nvar,nvar2)  !-------------------------!RT
                 do ivar=1,min(nhvar,nhvar2)  !------------------------!RT
                 !do ivar=1,nhvar2 !-----------------------------------!RT
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then        !Momentum
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*uold(ind_grid(i)+iskip,1)
                       end do
                    else if(ivar==ndim+2)then             ! Energy density
                       do i=1,ncache
                          ! Initialized as pressure so conv. to E_u+E_kin:
                          xx(i)=xx(i)/(gamma-1d0)
                          if(uold(ind_grid(i)+iskip,1) .ne. 0 ) then   !RT 
                             !if added by joki to avoid floating point exp
                             xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,2)**2/uold(ind_grid(i)+iskip,1)
#if NDIM>1
                             xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,3)**2/uold(ind_grid(i)+iskip,1)
#endif
#if NDIM>2
                             xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,4)**2/uold(ind_grid(i)+iskip,1)
#endif
                          endif                                        !RT
                          uold(ind_grid(i)+iskip,ivar)=xx(i)
                       end do
                    else                                 ! Passive scalars
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*uold(ind_grid(i)+iskip,1)
                       end do
                    endif
                 end do

                 ! Loop over ignored hydro vars                        !RT
                 do ivar=min(nhvar,nhvar2)+1,nhvar2                    !RT
                    read(ilun)xx                                       !RT
                 end do                                                !RT
                 ! Loop over RT variables                              !RT
                 do ivar=1,nRTvar2                                     !RT
                    read(ilun)xx                                       !RT
                    do i=1,ncache                                      !RT
                       uold(ind_grid(i)+iskip,nhvar+ivar)=xx(i)        !RT
                    end do                                             !RT
                 end do                                                !RT
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)
#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  else                        ! Not a restart -------------------------!RT
     rt_is_init_xion = .true. !----------------------------------------!RT
  end if

end subroutine init_hydro

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!RT
SUBROUTINE read_int(lun, param_name, value)                            !RT
                                                                       !RT
! Try to read a parameter from lun                                     !RT
!----------------------------------------------------------------------!RT
  integer::lun                                                         !RT
  character(*)::param_name                                             !RT
  character(128)::line,tmp                                             !RT
  integer::value                                                       !RT
!----------------------------------------------------------------------!RT
  rewind(unit=lun)                                                     !RT
  do                                                                   !RT
     read(lun, '(A128)', end=223) line                                 !RT
     if(index(line,trim(param_name)) .eq. 1) then                      !RT
        read(line,'(A13,I30)') tmp, value                              !RT
        return                                                         !RT
     endif                                                             !RT
  end do                                                               !RT
223 return   ! eof reached, didn't find the parameter                  !RT
                                                                       !RT
END SUBROUTINE read_int                                                !RT






