!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,            only: WP
   use geometry,             only: cfg
   use hypre_str_class,      only: hypre_str
   use ddadi_class,          only: ddadi
   use tpns_class,           only: tpns
   use vfs_class,            only: vfs
   use timetracker_class,    only: timetracker
   use ensight_class,        only: ensight
   use surfmesh_class,       only: surfmesh
   use vtk_class,         only: vtk
   use event_class,          only: event
   use monitor_class,        only: monitor
   use temp_transport
   use irl_fortran_interface
   use conservative_st
   
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),      public :: ps
   type(ddadi),          public :: vs
   type(tpns),           public :: fs
   type(vfs),            public :: vf
   type(tads),           public :: ts
   type(timetracker),    public :: time
   type(conservative_st_type), public :: cst

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(vtk)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,bubblefile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resH
   real(WP), dimension(:,:,:), allocatable :: resHG,resHL
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rho
   
   !> Problem definition
   logical :: second_bubble
   real(WP), dimension(3) :: center,center2,gravity
   real(WP) :: volume,radius,radius2,Ycent,Vrise
   real(WP) :: Vin,Vin_old,Vrise_ref,Ycent_ref,G,ti
   real(WP) :: Tbot,Ttop,Ly,Lx,dtps
contains


   !> Function that defines a level set function for a rising bubble problem
   function levelset_rising_bubble(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the bubble
      G=radius-sqrt(sum((xyz-center)**2))
      if (second_bubble) G=min(G,-radius2+sqrt(sum((xyz-center2)**2)))
   end function levelset_rising_bubble
   
   function levelset_RT(xyz,t) result(G)
      use param, only: param_read
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), dimension(3) :: dr
      real(WP), intent(in) :: t
      real(WP), PARAMETER :: PI = 3.141592653589793
      real(WP) :: G
      
      ! write(*,'(A,3F35.10)') '   center: ', center
      G=0.05*COS(2*PI*xyz(1))+xyz(2)
   end function levelset_RT

   !> Function that defines a level set function for a Sphere
   function levelset_sphere(xyz,t) result(G)
      use param, only: param_read
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G

      G=radius-sqrt(sum((xyz-center)**2))
   end function levelset_sphere

   !> Function that defines a level set function for a cylinder
   function levelset_cylinder(xyz,t) result(G)
      use param, only: param_read
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), dimension(3) :: dr
      real(WP), intent(in) :: t
      real(WP) :: G
      dr = (xyz-center)**2
      ! write(*,'(A,3F35.10)') '   center: ', center
      G=radius-sqrt((xyz(1)-center(1))**2+(xyz(2)-center(2))**2)
   end function levelset_cylinder

   !> Function that localizes the y+ side of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator
   
   
   !> Function that localizes the y- side of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator
   
   
   !> Routine that computes rise velocity
   subroutine rise_vel()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr
      real(WP) :: myYcent,myVrise,myvol,bubble_vol
      myYcent=0.0_WP
      myVrise=0.0_WP
      myvol=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               myYcent=myYcent+vf%cfg%ym(j)*(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
               myVrise=myVrise+Vi(i,j,k)*(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
               myvol=myvol+(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myYcent,Ycent     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myVrise,Vrise     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myvol  ,bubble_vol,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      Ycent=Ycent/bubble_vol
      Vrise=Vrise/bubble_vol
   end subroutine
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resH(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resHG(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resHL(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rho (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: plicnet,r2p,VFhi,VFlo,remap,flux_storage,remap_storage
         use mathtools, only: Pi
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=flux_storage,name='VOF')
         !vf%cons_correct=.false.
         !vf%thin_thld_max=1.5_WP
         !vf%twoplane_thld2=0.8_WP
         ! Initialize a bubble
         call param_read('Bubble position',center,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Bubble Radius',radius);
         ! Add a second one if needed
         second_bubble=param_exists('Bubble 2 position')
         if (second_bubble) then
            call param_read('Bubble 2 position',center2)
            call param_read('Bubble 2 volume',radius2,default=4.0_WP/3.0_WP*Pi*radius**3)
            radius2=(radius2*3.0_WP/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
         end if
         ! Generate interface
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_rising_bubble,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: clipped_neumann,dirichlet
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',gravity); fs%gravity=gravity
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         !ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Zero inflow velocity
         Vin=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      ! Create Surface Tension Object
      create_surface_tension_solver : block
         call cst%init(fs,vf,time)
         call param_read('Marangoni',cst%MarangoniOption)
         call param_read('Pressure',cst%PressureOption)
         call param_read('Smoothing Option',cst%SmoothingOption)
         call param_read('Surface Tension Option',cst%SurfaceTensionOption)
         call param_read('Curvature Option',cst%CurvatureOption)
         call param_read('PU Spread',cst%PU_spread)
         call cst%temp()
      end block create_surface_tension_solver
      
      ! Temperature Solver
      create_and_initialize_temperature_solver : block 
         integer :: i,j,k
         call ts%init(fs,vf,time)
         call ts%temp()
         ts%rhoL = fs%rho_l
         ts%rhoG = fs%rho_g 
         ts%cpL = 1.0_WP
         ts%cpG = 1.0_WP 
         ts%kL = 0.1_WP
         ts%kG = 0.1_WP
         ! Parameters
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         
         call param_read('k1',ts%kL)
         call param_read('k2',ts%kG)

         call param_read('cp1',ts%cpL)
         call param_read('cp2',ts%cpG)

         call param_read('Top Temperature',Ttop)
         call param_read('Bottom Temperature',Tbot)
         ! Initial Condition
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Note that the exact linear profile is:
                  ! (Ttop-Tbot)*y/Ly + Tbot.
                  ! Since it is linear, the value in a cell is equal to that at the center
                  ts%T(i,j,k) = (Ttop - Tbot) * (fs%cfg%ym(j)+Ly/2.0_WP)/Ly + Tbot
               enddo
            enddo
         enddo
         call ts%populate_enthalpy()
         ts%TG = ts%T
         ts%TL = ts%T
      end block create_and_initialize_temperature_solver
     


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=vtk(cfg=cfg,name='Marangoni-rise')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         call ens_out%add_scalar('Enthalpy',ts%H)
         call ens_out%add_vector('sigma_3d_x',cst%sigma_3D(:,:,:,1,1),cst%sigma_3D(:,:,:,1,2),cst%sigma_3D(:,:,:,1,3))
         call ens_out%add_vector('sigma_3d_y',cst%sigma_3D(:,:,:,2,1),cst%sigma_3D(:,:,:,2,2),cst%sigma_3D(:,:,:,2,3))
         call ens_out%add_vector('sigma_3d_z',cst%sigma_3D(:,:,:,3,1),cst%sigma_3D(:,:,:,3,2),cst%sigma_3D(:,:,:,3,3))

         call ens_out%add_scalar('sigma_xx_NoP',cst%sigma_xx_NoP)
         call ens_out%add_scalar('sigma_xy_NoP',cst%sigma_xy_NoP)
         call ens_out%add_scalar('sigma_yx_NoP',cst%sigma_yx_NoP)
         call ens_out%add_scalar('sigma_yy_NoP',cst%sigma_yy_NoP)

         call ens_out%add_scalar('cp',ts%cp)
         call ens_out%add_scalar('One Fluid Temp',ts%T)
         call ens_out%add_scalar('Interface Temperature',ts%Tinterface)
         call ens_out%add_vector('Temperatures',ts%TPmix,ts%TG,ts%TL)
         call ens_out%add_vector('Extrapolated Temperatures',ts%TPmix,ts%TGExtrap,ts%TLExtrap)

         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call rise_vel()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create bubble monitor
         bubblefile=monitor(fs%cfg%amRoot,'bubble')
         call bubblefile%add_column(time%n,'Timestep number')
         call bubblefile%add_column(time%t,'Time')
         call bubblefile%add_column(Ycent,'Y centroid')
         call bubblefile%add_column(Vrise,'Rise velocity')
         call bubblefile%add_column(Vin,'Inflow velocity')
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: harmonic_visc
      implicit none
      integer :: i,j,k
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Remember old scalar
         ts%Hold=ts%H
         call ts%populate_temperature()
         ts%Told = ts%T
         ts%TGold = ts%TG
         ts%TLold = ts%TL

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            ! mid-time Scalar
            ts%H = 0.5_WP*(ts%H+ts%Hold)
            ts%TG = 0.5_WP*(ts%TG+ts%TGold)
            ts%TL = 0.5_WP*(ts%TL+ts%TLold)
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! ============= SCALAR SOLVER =======================
            rho = vf%VF*fs%rho_l + (1.0_WP-vf%VF)*fs%rho_g
            resU = fs%U; resV = fs%V; resW = fs%W;
            ! Explicit Calculation of dHdt
            ! call ts%get_dHdt_SL(resH,resU,resV,resW,vf%detailed_face_flux, time%dt)
            call ts%get_dHdt(resH,resU,resV,resW)
            ! Explicit Update
            ts%H = ts%Hold + resH * time%dt 
            call ts%populate_temperature()
            ! ! Explicit Resdiual
            ! resH = -2.0_WP *(ts%H - ts%Hold)+time%dt*resH
            ! ! Form implicit residual
            ! call ts%solve_implicit(time%dt,resH,resU,resV,resW)
            ! ! Apply residual
            ! ts%H = 2*ts%H - ts%Hold + resH
            ! ! Apply other boundary conditions
            call ts%apply_bcond(time%t,time%dt)
            ! Extrapolate Values
            dtps = 1e-2
            call ts%extrapolate_fields_palmore(ts%TL,vf%VF,ts%TLExtrap,dtps)
            call ts%extrapolate_fields_palmore(ts%TG,1.0_WP - vf%VF,ts%TGExtrap,dtps)
            ts%TL = ts%TLExtrap
            ts%TG = ts%TGExtrap
            ! Step
            call ts%step_temperature_palmore(resHG,resHL,resU,resV,resW,time%dt)
            call ts%mix_temperature_palmore()

            call ts%extrapolate_fields_palmore(ts%TL,vf%VF,ts%TLExtrap,dtps)
            call ts%extrapolate_fields_palmore(ts%TG,1.0_WP - vf%VF,ts%TGExtrap,dtps)
            ts%TL = ts%TLExtrap
            ts%TG = ts%TGExtrap

            ! Dirichlet Boundary Condition
            do k=vf%cfg%kmino_,vf%cfg%kmaxo_
               do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                  do i=vf%cfg%imino_,vf%cfg%imaxo_
                     ! Note that the exact linear profile is:
                     ! (Ttop-Tbot)*y/Ly + Tbot.
                     ! Since it is linear, the value in a cell is equal to that at the center
                     if(fs%cfg%ym(j) .lt. -Ly/2.0 + fs%cfg%dy(j)) then 
                        ts%T(:,j,:) = Tbot 
                     endif

                     if(fs%cfg%ym(j) .gt. Ly/2.0 - fs%cfg%dy(j)) then 
                        ts%T(:,j,:) = Ttop 
                     endif
                  enddo
               enddo
            enddo
            call ts%populate_enthalpy()
            ! ===================================================

            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call cst%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms - adjust gravity if accelerating frame of reference
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call cst%add_surface_tension_jump()
            !call fs%add_surface_tension_jump_thin(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call rise_vel()
         call mfile%write()
         call cflfile%write()
         call bubblefile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation