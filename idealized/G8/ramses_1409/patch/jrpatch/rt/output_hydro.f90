!RT patch: RT variables are output as they are and not divided by gas  
!          density. Also there are checks on zero division to avoid 
!          floating point exceptions.
!          Also added call to output_rtInfo.
!************************************************************************
SUBROUTINE backup_hydro(filename)

!------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters !--------------------------------------------------!RT
  implicit none
  character(LEN=80)::filename,filedir,rt_filename !--------------------!RT

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  logical,allocatable :: not_zero(:) !--------------------------------!RT
  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10
     
  if(myid==1)then                                                      !RT
     call title(ifout-1,nchar)                                         !RT
     filedir='output_'//TRIM(nchar)//'/'                               !RT
     rt_filename=TRIM(filedir)//'info_rt_'//TRIM(nchar)//'.txt'        !RT
     call output_rtInfo(rt_filename)                                   !RT
  endif                                                                !RT

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           allocate(not_zero(ncache)) !--------------------------------!RT
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do ivar=1,nvar
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                       not_zero(i) = (xdp(i) /= 0.0) !-----------------!RT
                    end do
                 else if(ivar>=2.and.ivar<=ndim+1)then  ! Write vel. field
                    do i=1,ncache
                       if (not_zero(i)) then !-------------------------!RT
                          xdp(i) = uold(ind_grid(i)+iskip,ivar)          &
                                   / uold(ind_grid(i)+iskip,1)
                       end if !----------------------------------------!RT
                    end do
                 else if(ivar==ndim+2) then               ! Write pressure
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ndim+2)
                       if (not_zero(i)) then !-------------------------!RT
                          xdp(i)=xdp(i)                                  &
                                -0.5d0*uold(ind_grid(i)+iskip,2)**2      &
                                 /uold(ind_grid(i)+iskip,1)
#if NDIM>1
                          xdp(i)=xdp(i)                                  &
                                -0.5d0*uold(ind_grid(i)+iskip,3)**2      &
                                 /uold(ind_grid(i)+iskip,1)
#endif
#if NDIM>2
                          xdp(i)=xdp(i)                                  &
                                -0.5d0*uold(ind_grid(i)+iskip,4)**2      &
                                 /uold(ind_grid(i)+iskip,1)
#endif
                          xdp(i)=(gamma-1d0)*xdp(i)
                       end if !----------------------------------------!RT
                    end do
                 else if(ivar .le. nhvar) then !-----------------------!RT 
                    ! Write passive scalars if any
                    do i=1,ncache
                       if (not_zero(i)) then !-------------------------!RT
                          xdp(i)=uold(ind_grid(i)+iskip,ivar)            &
                                / uold(ind_grid(i)+iskip,1)
                       end if !----------------------------------------!RT
                    end do
                 else  ! ivar>hnvar -- RT variables -------------------!RT
                    do i=1,ncache                                      !RT
                        xdp(i)=uold(ind_grid(i)+iskip,ivar)            !RT
                    end do
                 endif
                 write(ilun)xdp
              end do
           end do
           deallocate(ind_grid, xdp)
           deallocate(not_zero) !--------------------------------------!RT
        end if
     end do
  end do
  close(ilun)
     
end subroutine backup_hydro

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! BEGIN RT SPECIFIC ROUTINES

!************************************************************************
SUBROUTINE output_rtInfo(filename)

! Output rt information into info_rt_XXXXX.txt
!------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  use cooling_module
  implicit none
  character(LEN=80)::filename
  integer::ilun
  real(dp)::scale_np,scale_pf
  character(LEN=80)::fileloc
!------------------------------------------------------------------------
  if(verbose)write(*,*)'Entering output_rtInfo'

  ilun=myid+10

  ! Conversion factor from user units to cgs units
  call rt_units(scale_np, scale_pf)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write run parameters
  write(ilun,'("nHvar       =",I11)')nHvar
  write(ilun,'("nRTvar      =",I11)')nRTvar
  write(ilun,'("nIons       =",I11)')nIons
  write(ilun,'("nPacs       =",I11)')nPacs
  write(ilun,'("iIons       =",I11)')iIons
  write(ilun,*)

  ! Write cooling parameters
  write(ilun,'("X_fraction  =",E23.15)')X
  write(ilun,'("Y_fraction  =",E23.15)')Y
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("unit_np     =",E23.15)')scale_np
  write(ilun,'("unit_pf     =",E23.15)')scale_pf
  write(ilun,'("rt_c_frac   =",E23.15)')rt_c_fraction
  write(ilun,*)

  ! Write polytropic parameters
  write(ilun,'("n_star      =",E23.15)')n_star
  write(ilun,'("T2_star     =",E23.15)')T2_star
  write(ilun,'("g_star      =",E23.15)')g_star
  write(ilun,*)
  call write_PacProps(.false.,ilun)

  close(ilun)

end subroutine output_rtInfo

!************************************************************************
SUBROUTINE write_PacProps(update,lun)

! Write photon package properties to file or std output.
! lun => File identifier (use 6 for std. output)
!------------------------------------------------------------------------
  use rt_parameters
  use amr_commons,only:myid
  implicit none
  logical::update
  integer::ip,lun
!------------------------------------------------------------------------
  if(myid .ne. 1) RETURN
  if(.not. update) then
     write(lun,*) 'Photon package properties------------------------------ '
  else
     write(lun,*) 'Photon properties have been changed to----------------- '
  endif
  write(lun,901) pacL0(:)
  write(lun,902) pacL1(:)
  write(lun,903) spec2pac(:)
  do ip=1,nPacs
     write(lun,906) ip
     write(lun,904) pac_csn(ip,:)
     write(lun,905) pac_egy(ip,:)
  enddo
  write (lun,*) '-------------------------------------------------------'

901 format ('  pacL0    [eV] =', 20f12.3)
902 format ('  pacL1    [eV] =', 20f12.3)
903 format ('  spec2pac      =', 20I12)
904 format ('  csn    [cm^2] =', 20(1pe12.3))
905 format ('  egy      [ev] =', 20f12.3)
906 format ('  ---Package',I2)

END SUBROUTINE write_PacProps

! END RT SPECIFIC ROUTINES
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX





