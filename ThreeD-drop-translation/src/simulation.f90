
!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,            only: WP
   use mathtools,            only: Pi
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
   use config_class, only: config
   ! use temp_transport
   use irl_fortran_interface
   use conservative_st
   
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),      public :: ps
   type(ddadi),          public :: vs
   type(tpns),           public :: fs
   type(vfs),            public :: vf
   ! type(tads),           public :: ts
   type(timetracker),    public :: time
   type(conservative_st_type), public :: cst

   !> Ensight postprocessing
   type(surfmesh) :: smesh,smesh_PU
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
   real(WP), dimension(3) :: center,center2,gravity,U1,U2
   real(WP) :: volume,radius,radius2,Ycent,Vrise
   real(WP) :: Vin,Vin_old,Vrise_ref,Ycent_ref,G,ti
   real(WP) :: Tbot,Ttop,Ly,Lx,dtps,RootMeanSquareVelocity,MaxVelocity
   real(WP) :: Fx_sum,Fy_sum,Fz_sum
   logical :: twoD

   real(WP) :: loc1,loc2,loc3,loc4,loc5
contains
   ! Locator Functions
   !> Function that localizes the left (x-) of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator
   
   !> Function that localizes the right (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   !> Function that localizes the top (y+) of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator

   !> Function that localizes the top (y-) of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator

   !> Function that localizes the end (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator

   !> Function that localizes the start (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator

   !> Enforce no-slip boundary conditions at all six domain boundaries
   subroutine set_wall_velocity()
      use tpns_class, only: bcond
      implicit none

      type(bcond), pointer :: mybc
      integer :: n,i,j,k,stag

      ! ============================================================
      ! X+ boundary
      ! ============================================================
      call fs%get_bcond('bc_xp',mybc)
      stag=min(mybc%dir,0)

      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n)
         j=mybc%itr%map(2,n)
         k=mybc%itr%map(3,n)

         ! Normal velocity: no penetration
         fs%U(i,j,k)=0.0_WP

         ! Tangential velocities: no slip
         fs%V(i+stag,j:j+1,k) = &
            -fs%V(i-mybc%dir+stag,j:j+1,k)

         fs%W(i+stag,j,k:k+1) = &
            -fs%W(i-mybc%dir+stag,j,k:k+1)
      end do

      ! ============================================================
      ! X- boundary
      ! ============================================================
      call fs%get_bcond('bc_xm',mybc)
      stag=min(mybc%dir,0)

      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n)
         j=mybc%itr%map(2,n)
         k=mybc%itr%map(3,n)

         ! Normal velocity
         fs%U(i,j,k)=0.0_WP

         ! Tangential velocities
         fs%V(i+stag,j:j+1,k) = &
            -fs%V(i-mybc%dir+stag,j:j+1,k)

         fs%W(i+stag,j,k:k+1) = &
            -fs%W(i-mybc%dir+stag,j,k:k+1)
      end do

      ! ============================================================
      ! Y+ boundary
      ! ============================================================
      call fs%get_bcond('bc_yp',mybc)
      stag=min(mybc%dir,0)

      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n)
         j=mybc%itr%map(2,n)
         k=mybc%itr%map(3,n)

         ! Normal velocity
         fs%V(i,j,k)=0.0_WP

         ! Tangential velocities
         fs%U(i:i+1,j+stag,k) = &
            -fs%U(i:i+1,j-mybc%dir+stag,k)

         fs%W(i,j+stag,k:k+1) = &
            -fs%W(i,j-mybc%dir+stag,k:k+1)
      end do

      ! ============================================================
      ! Y- boundary
      ! ============================================================
      call fs%get_bcond('bc_ym',mybc)
      stag=min(mybc%dir,0)

      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n)
         j=mybc%itr%map(2,n)
         k=mybc%itr%map(3,n)

         ! Normal velocity
         fs%V(i,j,k)=0.0_WP

         ! Tangential velocities
         fs%U(i:i+1,j+stag,k) = &
            -fs%U(i:i+1,j-mybc%dir+stag,k)

         fs%W(i,j+stag,k:k+1) = &
            -fs%W(i,j-mybc%dir+stag,k:k+1)
      end do

      ! ============================================================
      ! Z+ boundary
      ! ============================================================
      call fs%get_bcond('bc_zp',mybc)
      stag=min(mybc%dir,0)

      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n)
         j=mybc%itr%map(2,n)
         k=mybc%itr%map(3,n)

         ! Normal velocity
         fs%W(i,j,k)=0.0_WP

         ! Tangential velocities
         fs%U(i:i+1,j,k+stag) = &
            -fs%U(i:i+1,j,k-mybc%dir+stag)

         fs%V(i,j:j+1,k+stag) = &
            -fs%V(i,j:j+1,k-mybc%dir+stag)
      end do

      ! ============================================================
      ! Z- boundary
      ! ============================================================
      call fs%get_bcond('bc_zm',mybc)
      stag=min(mybc%dir,0)

      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n)
         j=mybc%itr%map(2,n)
         k=mybc%itr%map(3,n)

         ! Normal velocity
         fs%W(i,j,k)=0.0_WP

         ! Tangential velocities
         fs%U(i:i+1,j,k+stag) = &
            -fs%U(i:i+1,j,k-mybc%dir+stag)

         fs%V(i,j:j+1,k+stag) = &
            -fs%V(i,j:j+1,k-mybc%dir+stag)
      end do

      ! ============================================================
      ! Synchronize velocity fields
      ! ============================================================
      call fs%cfg%sync(fs%U)
      call fs%cfg%sync(fs%V)
      call fs%cfg%sync(fs%W)

   end subroutine set_wall_velocity

   !> Function that defines a level set function for a rising bubble problem
   function levelset_rising_bubble(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the bubble
      G=radius-sqrt(sum((xyz-center)**2))
      G=max(G,radius2-sqrt(sum((xyz-center2)**2)))
   end function levelset_rising_bubble

   function levelset_PU(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), dimension(3) :: dr
      real(WP), intent(in) :: t
      real(WP) :: G

      call cst%getValuePU(xyz,G)
      G = -G
   end function levelset_PU

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
               myYcent=myYcent+vf%cfg%ym(j)*(vf%VF(i,j,k))*cfg%vol(i,j,k)
               myVrise=myVrise+Vi(i,j,k)*(vf%VF(i,j,k))*cfg%vol(i,j,k)
               myvol=myvol+(vf%VF(i,j,k))*cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myYcent,Ycent     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myVrise,Vrise     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myvol  ,bubble_vol,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      Ycent=Ycent/bubble_vol
      Vrise=Vrise/bubble_vol
   end subroutine
   
   function Calc_RMS_Velocity() result(VrmsGlobal)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MAX
      implicit none

      integer :: i,j,k,ierr
      real(WP) :: Vrms,count,VrmsGlobal,countGlobal
      real(WP) :: Vmax,Vmag2,VmaxGlobal
      real(WP) :: Urel,Vrel,Wrel
      real(WP) :: dvol

      Vrms = 0.0_WP
      Vmax = 0.0_WP
      count = 0.0_WP

      do k = cfg%kmin_,cfg%kmax_
         do j = cfg%jmin_,cfg%jmax_
            do i = cfg%imin_,cfg%imax_

               ! Velocity relative to imposed droplet velocity
               Urel = fs%U(i,j,k) - U1(1)
               Vrel = fs%V(i,j,k) - U1(2)
               Wrel = fs%W(i,j,k) - U1(3)

               ! Velocity magnitude squared
               Vmag2 = Urel**2 + Vrel**2 + Wrel**2

               ! Cell volume
               dvol = cfg%dx(i)*cfg%dy(j)*cfg%dz(k)

               ! RMS velocity
               Vrms = Vrms + Vmag2*dvol
               count = count + dvol

               ! Maximum velocity
               Vmax = max(Vmax,sqrt(Vmag2))

            end do
         end do
      end do

      ! MPI reductions
      call MPI_ALLREDUCE(Vrms,VrmsGlobal,1,MPI_DOUBLE_PRECISION,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(count,countGlobal,1,MPI_DOUBLE_PRECISION,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(Vmax,VmaxGlobal,1,MPI_DOUBLE_PRECISION,MPI_MAX,cfg%comm,ierr)

      VrmsGlobal = sqrt(VrmsGlobal/countGlobal)
      MaxVelocity = VmaxGlobal

   end function Calc_RMS_Velocity

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
         use vfs_class,only: r2p,lvira,VFhi,VFlo,jibben,plicnet,jibbenPU,jibbenSq,PUppic,jibbenPLIC
         use mathtools, only: Pi
         integer :: i,j,k,n,si,sj,sk,InterfaceTypes
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         real(WP) :: B
         ! Create a VOF solver
         call param_read('Interface Type',InterfaceTypes,default=1)
         SELECT CASE (InterfaceTypes)
            case (1)
               call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
            case (2) 
               call vf%initialize(cfg=cfg,reconstruction_method=jibben,name='VOF')
            case (3) 
               call vf%initialize(cfg=cfg,reconstruction_method=plicnet,name='VOF')
            case (4) 
               call vf%initialize(cfg=cfg,reconstruction_method=jibbenPU,name='VOF')
            case (5) 
               call vf%initialize(cfg=cfg,reconstruction_method=jibbenSq,name='VOF')
            case (6) 
               call vf%initialize(cfg=cfg,reconstruction_method=jibbenPLIC,name='VOF')
            case default 
               call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
         END SELECT
      
         
         !vf%cons_correct=.false.
         !vf%thin_thld_max=1.5_WP
         !vf%twoplane_thld2=0.8_WP
         ! Initialize a bubble
         call param_read('Droplet Radius',radius)
         call param_read('Two Dimensional',twoD)
         ! Now that we have radius and Impact Parameter
         ! Assume the droplets move in the X directio only, and are separated in the Y direction. 
         center = (/0.0_WP,0.0_WP,0.0_WP/)
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
                  if(TwoD) then 
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cylinder,0.0_WP,amr_ref_lvl)
                  else
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
                  endif
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
         use tpns_class,      only: clipped_neumann,dirichlet,wall
         integer :: i,j,k
         real(WP),dimension(3) :: xyz
         real(WP) :: Re,La,U,delta,dtCap,minTc,minTv,nx,time_scale_factor,Tv,Tc

         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Density Ratio',fs%rho_g); fs%rho_g=fs%rho_l/fs%rho_g
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
         ! initial field
         call param_read('Droplet Velocity',U1)
         fs%U = U1(1)
         fs%V = U1(2)
         fs%W = U1(3)
         ! X boundaries
         ! call fs%add_bcond(name='bc_xp',type=dirichlet,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         ! call fs%add_bcond(name='bc_xm',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)

         ! ! Y boundaries
         ! call fs%add_bcond(name='bc_yp',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         ! call fs%add_bcond(name='bc_ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)

         ! ! Z boundaries
         ! call fs%add_bcond(name='bc_zp',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=zp_locator)
         ! call fs%add_bcond(name='bc_zm',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=zm_locator)

         ! call set_wall_velocity()
         ! Surface Tension
         call param_read('Surface tension coefficient',fs%sigma)
         call param_read('Laplace Number',La)
         fs%visc_l = sqrt(fs%rho_l*fs%sigma*2*radius/La)
         call param_read('Viscosity Ratio',fs%visc_g); fs%visc_g = fs%visc_l/fs%visc_g

         call param_read('Min Viscous Time',minTv)
         call param_read('Min Capillary Time', minTc)
         call param_read('Lx',Lx)
         call param_read('nx',nx)
         Tv = 4*radius*radius*fs%rho_l/fs%visc_l
         Tc = sqrt(fs%rho_l*8*radius*radius*radius/fs%sigma)
         delta = Lx/nx
         dtCap = sqrt(fs%rho_l*delta*delta*delta/(Pi*fs%sigma))

         call param_read('Timestep factor',time_scale_factor)
         time%dtmax = min(time_scale_factor*dtCap,time%dtmax)
         time%dt=time%dtmax
         call fs%apply_bcond(time%t,time%dt)
         time%tmax = max(minTv*Tv,minTc*Tc)
         if(cfg%amRoot) then 
            print *, "============== TIME INFO ================="
            print *, "Viscous Time = ",Tv
            print *, "Capillary Time = ",Tc
            print *, "dt = ",time%dt
            print *, "tmax = ",time%tmax

            print *, "============== FLUID INFO ================"
            print *, "Liquid Density = ",fs%rho_l
            print *, "Gas Density = ",fs%rho_g
            print *, "Density Ratio = ",fs%rho_l/fs%rho_g
            print *, "Liquid Viscosity = ",fs%visc_l
            print *, "Gas Viscosity = ",fs%visc_g
            print *, "Viscosity Ratio = ",fs%visc_l/fs%visc_g
            print *, "Surface Tension Coefficient = ",fs%sigma
         endif
         ! Zero inflow velocity
         Vin=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      ! Create Surface Tension Object
      create_surface_tension_solver : block
         call cst%init(fs,vf,time)
         call param_read('Marangoni',cst%MarangoniOption,default = [0.0_WP,0.0_WP,0.0_WP])
         call param_read('Pressure',cst%PressureOption,default = 0.0_WP)
         call param_read('Smoothing Option',cst%SmoothingOption,default = 1)
         call param_read('Surface Tension Option',cst%SurfaceTensionOption,default = 8)
         call param_read('Levelset Option',cst%LevelsetOption,default = 2)
         call param_read('Matching Option',cst%MatchingOption,default = 1)
         call param_read('Curvature Option',cst%CurvatureOption,default = 1)
         call param_read('PU Spread',cst%PU_spread,default = 2.0_WP)
         call param_read('Two Dimensional',cst%TwoD,default = .false.)
         if(cfg%amRoot) then 
            print *, "============== SURFACE TENSION INFO ================"
            print *, "Marangoni Option = ",cst%MarangoniOption
            print *, "Pressure Option = ",cst%PressureOption
            print *, "Smoothing Option = ",cst%SmoothingOption
            print *, "Surface Tension Option = ",cst%SurfaceTensionOption
            print *, "Levelset Option = ",cst%LevelsetOption
            print *, "Matching Option = ",cst%MatchingOption
            print *, "Curvature Option = ",cst%CurvatureOption
            print *, "PU Spread = ",cst%PU_spread
            print *, "Two Dimensional = ",cst%TwoD
            print *, "===================================================="
         endif
      end block create_surface_tension_solver
     


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=2,name='plic')
         smesh%varname(1)='curv'
         smesh%varname(2)='interface_type'

         ! smesh%varname(3)='STx'
         ! smesh%varname(4)='STy'
         ! smesh%varname(5)='STz'
         ! smesh%varname(6)='STmag'

         ! smesh%varname(7)='Sigma_xx'
         ! smesh%varname(8)='Sigma_xy'
         ! smesh%varname(9)='Sigma_xz'
         ! smesh%varname(10)='Sigma_yx'
         ! smesh%varname(11)='Sigma_yy'
         ! smesh%varname(12)='Sigma_yz'
         ! smesh%varname(13)='Sigma_zx'
         ! smesh%varname(14)='Sigma_zy'
         ! smesh%varname(15)='Sigma_zz'

         ! smesh%varname(16)='CSFx'
         ! smesh%varname(17)='CSFy'
         ! smesh%varname(18)='CSFz'
         ! smesh%varname(19)='CSFmag'


         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         call add_surfgrid_variable(smesh,1,vf%curv,vf)  
         call add_surfgrid_variable(smesh,2,vf%interface_type,vf)    

         ! call add_surfgrid_variable(smesh,3,fs%Pjx,vf)  
         ! call add_surfgrid_variable(smesh,4,fs%Pjy,vf)    
         ! call add_surfgrid_variable(smesh,5,fs%Pjz,vf)  
         ! call add_surfgrid_variable(smesh,6,cst%Fst_mag,vf)    

         ! call add_surfgrid_variable(smesh,7,cst%sigma_3D(:,:,:,1,1),vf)  
         ! call add_surfgrid_variable(smesh,8,cst%sigma_3D(:,:,:,1,2),vf)    
         ! call add_surfgrid_variable(smesh,9,cst%sigma_3D(:,:,:,1,3),vf)  
         ! call add_surfgrid_variable(smesh,10,cst%sigma_3D(:,:,:,2,1),vf)  
         ! call add_surfgrid_variable(smesh,11,cst%sigma_3D(:,:,:,2,2),vf)  
         ! call add_surfgrid_variable(smesh,12,cst%sigma_3D(:,:,:,2,3),vf)    
         ! call add_surfgrid_variable(smesh,13,cst%sigma_3D(:,:,:,3,1),vf)  
         ! call add_surfgrid_variable(smesh,14,cst%sigma_3D(:,:,:,3,2),vf)    
         ! call add_surfgrid_variable(smesh,15,cst%sigma_3D(:,:,:,3,3),vf)    

         ! call add_surfgrid_variable(smesh,16,cst%PjxD,vf)  
         ! call add_surfgrid_variable(smesh,17,cst%PjyD,vf)    
         ! call add_surfgrid_variable(smesh,18,cst%PjzD,vf)  
         ! call add_surfgrid_variable(smesh,19,cst%CSF_Mag,vf)    

      end block create_smesh
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=vtk(cfg=cfg,name='DropCollision')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output 
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         call ens_out%add_vector('sigma_3d_x',cst%sigma_3D(:,:,:,1,1),cst%sigma_3D(:,:,:,1,2),cst%sigma_3D(:,:,:,1,3))
         call ens_out%add_vector('sigma_3d_y',cst%sigma_3D(:,:,:,2,1),cst%sigma_3D(:,:,:,2,2),cst%sigma_3D(:,:,:,2,3))
         call ens_out%add_vector('sigma_3d_z',cst%sigma_3D(:,:,:,3,1),cst%sigma_3D(:,:,:,3,2),cst%sigma_3D(:,:,:,3,3))

         call ens_out%add_vector('OLD sigma_3d_x',cst%sigma_3D_OLD(:,:,:,1,1),cst%sigma_3D_OLD(:,:,:,1,2),cst%sigma_3D_OLD(:,:,:,1,3))
         call ens_out%add_vector('OLD sigma_3d_y',cst%sigma_3D_OLD(:,:,:,2,1),cst%sigma_3D_OLD(:,:,:,2,2),cst%sigma_3D_OLD(:,:,:,2,3))
         call ens_out%add_vector('OLD sigma_3d_z',cst%sigma_3D_OLD(:,:,:,3,1),cst%sigma_3D_OLD(:,:,:,3,2),cst%sigma_3D_OLD(:,:,:,3,3))

         call ens_out%add_vector('d sigma_3d_x',cst%dSigma(:,:,:,1,1),cst%dSigma(:,:,:,1,2),cst%dSigma(:,:,:,1,3))
         call ens_out%add_vector('d sigma_3d_y',cst%dSigma(:,:,:,2,1),cst%dSigma(:,:,:,2,2),cst%dSigma(:,:,:,2,3))
         call ens_out%add_vector('d sigma_3d_z',cst%dSigma(:,:,:,3,1),cst%dSigma(:,:,:,3,2),cst%dSigma(:,:,:,3,3))
         
         call ens_out%add_vector('diff ST Force',cst%Fst_x_3D_diff,cst%Fst_y_3D_diff,cst%Fst_z_3D_diff)

         call ens_out%add_vector('OLD ST Force',cst%Fst_x_3D_OLD,cst%Fst_y_3D_OLD,cst%Fst_z_3D_OLD)

         call ens_out%add_vector('ST Force',fs%Pjx,fs%Pjy,fs%Pjz)
         call ens_out%add_vector('CSF Force',cst%PjxD,cst%PjyD,cst%PjzD)

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
         call mfile%add_column(RootMeanSquareVelocity,'RMS Velocity')
         call mfile%add_column(MaxVelocity,'Max Velocity')
         call mfile%add_column(fs%divmax,'Max Div')
         call mfile%add_column(loc1,'VFint Loc 1')
         call mfile%add_column(loc2,'VFint Loc 2')
         call mfile%add_column(loc3,'VFint Loc 3')
         call mfile%add_column(loc4,'VFint Loc 4')
         call mfile%add_column(loc5,'VFint Loc 5')
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
      use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank
      use iso_fortran_env, only: error_unit
      implicit none

      integer :: i,j,k
      integer :: myrank, ierr

      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
      ! Perform time integration
      ! ============================================================
      ! Main time integration loop
      ! ============================================================  
      do while (.not.time%done())

         ! call debug_print('START TIMESTEP')

         ! ----------------------------------------------------------
         ! Time integration
         ! ----------------------------------------------------------
         ! call debug_print('START get_cfl')
         call fs%get_cfl(time%dt,time%cfl)
         ! call debug_print('END get_cfl')

         ! call debug_print('START adjust_dt')
         call time%adjust_dt()
         ! call debug_print('END adjust_dt')

         ! call debug_print('START increment')
         call time%increment()
         ! call debug_print('END increment')

         ! ----------------------------------------------------------
         ! Save previous timestep
         ! ----------------------------------------------------------
         ! call debug_print('START save old fields')

         vf%VFold=vf%VF

         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         call vf%get_max()
         loc1 = vf%VFint
         ! call debug_print('END save old fields')

         ! ----------------------------------------------------------
         ! Density
         ! ----------------------------------------------------------
         ! call debug_print('START get_olddensity')
         call fs%get_olddensity(vf=vf)
         ! call debug_print('END get_olddensity')

         ! ----------------------------------------------------------
         ! VOF advection
         ! ----------------------------------------------------------
         ! call debug_print('START vf advance')
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         ! call debug_print('END vf advance')
         call vf%get_max()
         loc2 = vf%VFint
         ! ----------------------------------------------------------
         ! Viscosity
         ! ----------------------------------------------------------
         ! call debug_print('START get_viscosity')
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         ! call debug_print('END get_viscosity')

         ! ==========================================================
         ! Sub-iterations
         ! ==========================================================
         do while (time%it.le.time%itmax)

            ! call debug_print('START SUBITERATION')

            ! Mid-time velocity
            ! call debug_print('START mid-time velocity')

            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! call debug_print('END mid-time velocity')

            ! Advection
            ! call debug_print('START prepare_advection_upwind')
            call fs%prepare_advection_upwind(dt=time%dt)
            ! call debug_print('END prepare_advection_upwind')

            ! Surface tension
            ! call debug_print('START get_dmomdt')
            call cst%get_dmomdt(resU,resV,resW)
            call vf%get_max()
            loc3 = vf%VFint
            ! call debug_print('END get_dmomdt')

            ! Gravity
            ! call debug_print('START addsrc_gravity')
            ! call fs%addsrc_gravity(resU,resV,resW)
            ! call debug_print('END addsrc_gravity')

            ! Explicit residual
            ! call debug_print('START explicit residual')

            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW

            ! call debug_print('END explicit residual')

            ! Implicit solve
            ! call debug_print('START solve_implicit')
            call fs%solve_implicit(time%dt,resU,resV,resW)
            ! call debug_print('END solve_implicit')

            ! Apply residuals
            ! call debug_print('START apply residuals')

            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! call debug_print('END apply residuals')

            ! Boundary conditions
            ! call debug_print('START apply_bcond')
            call fs%apply_bcond(time%t,time%dt)
            ! call debug_print('END apply_bcond')

            ! Laplacian
            ! call debug_print('START update_laplacian')
            call fs%update_laplacian()
            ! call debug_print('END update_laplacian')

            ! Mass flux
            ! call debug_print('START correct_mfr')
            call fs%correct_mfr()
            ! call debug_print('END correct_mfr')

            ! Divergence
            ! call debug_print('START get_div')
            call fs%get_div()
            ! call debug_print('END get_div')
            call vf%get_max()
            loc4 = vf%VFint
            ! Surface tension jump
            ! call debug_print('START add_surface_tension_jump')
            call cst%add_surface_tension_jump()
            call vf%get_max()
            loc5 = vf%VFint
            ! call debug_print('END add_surface_tension_jump')

            ! Pressure RHS
            ! call debug_print('START pressure RHS')

            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP

            ! call debug_print('END pressure RHS')

            ! Pressure solve
            ! call debug_print('START pressure solve')
            call fs%psolv%solve()
            ! call debug_print('END pressure solve')

            ! Pressure shift
            ! call debug_print('START shift_p')
            call fs%shift_p(fs%psolv%sol)
            ! call debug_print('END shift_p')

            ! Pressure gradient
            ! call debug_print('START get_pgrad')
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            ! call debug_print('END get_pgrad')

            ! Velocity correction
            ! call debug_print('START velocity correction')

            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W

            ! call debug_print('END velocity correction')

            time%it=time%it+1

            ! call debug_print('END SUBITERATION')

         end do

         ! ==========================================================
         ! Postprocessing
         ! ==========================================================

         ! call debug_print('START interp_vel')
         call fs%interp_vel(Ui,Vi,Wi)
         ! call debug_print('END interp_vel')

         ! call debug_print('START final get_div')
         call fs%get_div()
         ! call debug_print('END final get_div')

         ! ----------------------------------------------------------
         ! Ensight output
         ! ----------------------------------------------------------
         if (ens_evt%occurs()) then

            ! call debug_print('START update_surfmesh')
            call vf%update_surfmesh(smesh)
            ! call debug_print('END update_surfmesh')

            ! ! call debug_print('START add_surfgrid_variable 1')
            ! call add_surfgrid_variable(smesh,1,vf%curv,vf)
            ! ! call debug_print('END add_surfgrid_variable 1')

            ! ! call debug_print('START add_surfgrid_variable 2')
            ! call add_surfgrid_variable(smesh,2,vf%interface_type,vf)
            ! ! call debug_print('END add_surfgrid_variable 2')

            ! call debug_print('START ens_out write_data')
            call ens_out%write_data(time%t)
            ! call debug_print('END ens_out write_data') 

         endif

         ! ----------------------------------------------------------
         ! Monitoring        
         ! ----------------------------------------------------------
         ! call debug_print('START fs get_max')
         call fs%get_max()
         ! call debug_print('END fs get_max') 

         ! call debug_print('START vf get_max')
         call vf%get_max()
         ! call debug_print('END vf get_max')

         ! call debug_print('START RMS velocity')
         RootMeanSquareVelocity=Calc_RMS_Velocity()
         ! call debug_print('END RMS velocity')

         ! call debug_print('START monitor writes')

         call mfile%write()
         call cflfile%write()
         call bubblefile%write()

         ! call debug_print('END monitor writes')

         ! call debug_print('END TIMESTEP')

      end do
            
      contains

      subroutine debug_print(message)
         character(len=*), intent(in) :: message

         write(error_unit,'(A,I4,A,I8,A,I4,A,A)') &
            '[Rank ',myrank,'] Step ',time%n, &
            ' Iter ',time%it,': ',trim(message)

         flush(error_unit)
      end subroutine debug_print

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
   
   

   !> Make a surface scalar variable
   subroutine add_surfgrid_variable(smesh,var_index,A,vf_local)
      use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
      use vfs_class,only: VFhi,VFlo
      implicit none
      class(surfmesh), intent(inout) :: smesh
      integer, intent(in) :: var_index
      type(vfs) :: vf_local
      real(WP), dimension(vf_local%cfg%imino_:vf_local%cfg%imaxo_,vf_local%cfg%jmino_:vf_local%cfg%jmaxo_,vf_local%cfg%kmino_:vf_local%cfg%kmaxo_), intent(in) :: A 
      integer :: i,j,k,n,shape,np,nplane,nbt
      
      ! Fill out arrays
      ! print *,"Fill Arrays"
      if ((smesh%nPoly+smesh%nBezierTri).gt.0) then
         np=0; nbt = 0
         ! Start with quadratic surfaces
         do k=vf_local%cfg%kmin_,vf_local%cfg%kmax_
            do j=vf_local%cfg%jmin_,vf_local%cfg%jmax_
               do i=vf_local%cfg%imin_,vf_local%cfg%imax_
                  shape=getNumberOfTriangles(vf_local%interface_mixed_surface(i,j,k))
                  if (shape.gt.0) then
                     do n=1,shape
                        nbt=nbt+1
                        smesh%trivar(var_index,nbt)=A(i,j,k)
                     end do
                  end if
               end do
            end do
         end do
         ! Then do planes
         ! print *,"Fill Planes"
         do k=vf_local%cfg%kmin_,vf_local%cfg%kmax_
            do j=vf_local%cfg%jmin_,vf_local%cfg%jmax_
               do i=vf_local%cfg%imin_,vf_local%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf_local%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(vf_local%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        ! Set nplane variable
                        smesh%polyvar(var_index,np)=A(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      else
         smesh%polyvar(var_index,1)=1
      end if      
   ! print *,"Fill End"
   end subroutine add_surfgrid_variable


end module simulation