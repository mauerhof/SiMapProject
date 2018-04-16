!################################################################
!################################################################
!################################################################
!################################################################
subroutine rt_make_boundary_hydro(ilevel)
  use amr_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref,alt
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref
  
  real(dp)::switch,dx,dx_loc,scale
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nrtvar),save::uu

  integer::rtType !RT-- Type of rt variable =0,1,2,3
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2 &
       ,scale_Np,scale_Fp
  real(dp):: F_star, cE_up, Fy_limit, Fx_up, Fy_up, Fmag_up, fred_up, chi_up  &
       ,nbold, cPy_up, eps, Fx_down, Fy_down, cE_down, Fy_half, fred_down       &
       ,cPy_down, cPy_half, factor, Fmag
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Loop over boundaries
  do ibound=1,nboundary
     ! Compute direction of reference neighbors
     boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
     if(boundary_dir==1)inbor=2
     if(boundary_dir==2)inbor=1
     if(boundary_dir==3)inbor=4
     if(boundary_dir==4)inbor=3
     if(boundary_dir==5)inbor=6
     if(boundary_dir==6)inbor=5

     ! Compute index of reference cells
     ! Reflexive boundary
     if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     ! Free boundary
     if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! Imposed boundary (used only for flag1)
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)

     ! Velocity sign switch for reflexive boundary conditions
     gs=(/1,1,1/)
     if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
     if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
     if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector 
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do
        
        ! Gather neighboring reference grid
        do i=1,ngrid
           ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
              
           ! Gather neighboring reference cell
           iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
           do i=1,ngrid
              ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
           end do

           ! Wall and free boundary conditions
           if((boundary_type(ibound)/10).ne.2)then
              
              ! Gather reference hydro variables
              do ivar=1,nrtvar
                 do i=1,ngrid
                    uu(i,ivar)=rtuold(ind_cell_ref(i),ivar)
                 end do
              end do

              ! BEGIN GONZALEZ 1------------------------------------------
              if(boundary_dir .ne. 1) then
                 do i=1,ngrid ! zero-valued boundaries by default
                    uu(i,1) = 1d-50
                    uu(i,2) = 0d0
                    uu(i,3) = 0d0
                 end do
              endif

!#define GONZALEZ
#ifdef GONZALEZ
              if(boundary_dir .eq. 1) then ! Left (emitting) boundary
                 do i=1,ngrid
                    uu(i,1)   = 3.395d16 / scale_Fp /rt_c  ! Gonzalez 1
                    uu(i,2)   = 3.395d16 / scale_fp        ! Gonzalez 1
                    uu(i,3)   = 0d0
                 end do
              endif
#endif       

#define DAVIS
#ifdef DAVIS
              if(boundary_dir .eq. 1) then ! Left (emitting) boundary
                 do i=1,ngrid
                    !! Fx_down=F_star, and adjust cE_down accordingly
                    F_star = 3.395d16/scale_Fp
                    cE_up = uu(i,1)*rt_c
                    Fx_up = uu(i,2)
                    Fy_up = uu(i,3)
                    Fx_down = F_star
                    cE_down = Fx_down - Fx_up + cE_up
                    uu(i,1) = cE_down/rt_c
                    uu(i,2) = Fx_down
                    uu(i,3) = 0d0

                    !uu(i,1) = cE_down/rt_c
                    !uu(i,2) = Fx_down
                    !uu(i,3) = 0d0
                 end do
              endif       
#endif       
              ! END GONZALEZ 1--------------------------------------------

              ! Scatter to boundary region
              do ivar=1,nrtvar
                 switch=1
                 ! ------------------------------------------------------------
                 ! need to switch photon flux in RT...
                 ! rtVar= ivar = 1,2,3,4,.. = Np, Fx, Fy, Fz, Np, Fx,..
                 rtType = mod(ivar-1, ndim+1)
                 if(rtType .ne. 0) switch=gs(rtType) !0=Np, 1=Fx, 2=Fy, 3=Fz
                 ! ------------------------------------------------------------
                 do i=1,ngrid
                    rtuold(ind_cell(i),ivar)=uu(i,ivar)*switch
                 end do
                      
              end do
              
              ! Imposed boundary conditions
           else
              
              ! Compute cell center in code units
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
                 end do
              end do
              
              ! Rescale position from code units to user units
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                 end do
              end do
              
              call rt_boundana(xx,uu,dx_loc,ibound,ngrid)

              ! Scatter variables
              do ivar=1,nrtvar
                 do i=1,ngrid
                    rtuold(ind_cell(i),ivar)=uu(i,ivar)
                 end do
              end do
                 
           end if
              
        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering rt_make_boundary_hydro for level ',I2)

end subroutine rt_make_boundary_hydro


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE flux_root_func(Fy_down, F_star, cE_up, Fy_up, cPy_up, f, df)

! Returns equilibrium radiation/dust temperature, and T derivative, given
! initial radiation energy density Ei, temperature Ti, T-E conversion Cv
! light speed fraction fc (all in cgs).
!-------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::Fy_down, F_star, cE_up, Fy_up, cPy_up, f, df
  real(dp)::cE_down,cPy_down,f_red, chinum, chidiv, chi
!-------------------------------------------------------------------------
  cE_down = 2d0 * F_star + cE_up - Fy_up - Fy_down
  f_red = abs(Fy_down)/cE_down
  chinum = (3d0+4d0*f_red**2)
  chidiv = 5d0+2d0*sqrt(4d0-3d0*f_red**2)
  chi = chinum / chidiv
  cPy_down = cE_down * chi
  f = 2d0*F_star - cPy_down - cPy_up + Fy_up - Fy_down
!print*,'f=',f,F_star,cPy_down,cE_down,chi,cPy_up,Fy_up,Fy_down
  df = -chi + (1d0+Fy_down/cE_down) /chidiv**2 &
            * (8d0*f_red*chidiv + 6d0*f_red*chinum/sqrt(4d0-3d0*f_red**2))
END SUBROUTINE flux_root_func

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE root_fyghost(F_y_limit, eps, F_star, cE_up, Fx_up, Fy_up, cPy_up &
     ,F_y_down) 

! Find y-directed photon flux in the lower ghost zone.
! The flux is a polynomial root, derived from the following constraints:
! ==F_x_ghost = 0
! ==Fy_half = F_star (Gives E_ghost(F_y_ghost))
! ==Pyy_half = F_star/c (Gives F_y_ghost from the neighbour)
! Additionally there are the constraints:
! ==cE_ghost .ge. |F_y_ghost|
! ==cE_ghost .ge. 0d0
!-------------------------------------------------------------------------
  use amr_commons
  implicit none
  INTEGER,parameter::MAXIT=100
  REAL(dp)::F_y_limit,eps,F_star,cE_up,Fx_up,Fy_up,cPy_up
  INTEGER::j 
  REAL(dp)::df,dFy,dFyold,f,fh,fl,temp,Fyh,Fyl
  REAL(dp)::F_y_down,epssq
!-------------------------------------------------------------------------
  epssq=eps**2
  !Fyl=-F_y_limit ; Fyh=F_y_limit
  Fyl=F_star ; Fyh=F_y_limit  
  call flux_root_func(Fyl, F_star, cE_up, Fy_up, cPy_up, fl, df) 
  call flux_root_func(Fyh, F_star, cE_up, Fy_up, cPy_up, fh, df)

  if(abs(fl) .le. epssq) then 
     F_y_down=Fyl
     return 
  else if(abs(fh) .le. epssq) then 
     F_y_down=Fyh 
     return 
  endif
  if((fl>0..and.fh>0.).or.(fl<0..and.fh<0.)) then
     write(*,*) 'root must be bracketed in root_fyghost' 
     print*,'fl,fh: ',fl,fh,F_y_limit/F_star
     print*,'Up cond: ',cE_up/F_star, Fx_up/F_star, Fy_up/F_star, cPy_up/F_star
     !print*,'Dn cond: ',cE_down/F_star, 0d0, Fy_down/F_star, cPy_down/F_star
     !print*,'Half: ',Fy_half/F_star,cPy_half/F_star
     stop
  endif
  if(fl.ge.0.) then 
     temp=Fyh
     Fyh=Fyl
     Fyl=temp
  endif
  F_y_down=.5*(Fyl+Fyh) 
  dFyold=abs(Fyh-Fyl) 
  dFy=dFyold 
  call flux_root_func(F_y_down, F_star, cE_up, Fy_up, cPy_up, f, df)
  do j=1,MAXIT 
     !print*,Fyl,Fyh,f
     if(((F_y_down-Fyh)*df-f)*((F_y_down-Fyl)*df-f)>0.& 
          .or.abs(2.*f)>abs(dFyold*df) ) then 
        dFyold=dFy
        dFy=0.5*(Fyh-Fyl) 
        F_y_down=Fyl+dFy
        if(Fyl==F_y_down) return 
     else 
        dFyold=dFy
        dFy=f/df 
        temp=F_y_down 
        F_y_down=F_y_down-dFy
        if(temp==F_y_down) return 
     endif
     if(abs(dFy)<epssq) return 
     call flux_root_func(F_y_down, F_star, cE_up, Fy_up, cPy_up, f, df)
     if(f<0.) then 
        Fyl=F_y_down
     else 
        Fyh=F_y_down
     endif
  end do
  write(*,*) 'root_fyghost exceeding maximum iterations' 
  stop
END SUBROUTINE root_fyghost 

