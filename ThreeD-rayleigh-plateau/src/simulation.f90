!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use mathtools,         only: Pi
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use vtk_class,         only: vtk
   use partmesh_class,    only: partmesh
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use irl_fortran_interface
   use conservative_st 
   implicit none
   private

   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs           !< DDADI solver for velocity
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   type(conservative_st_type), public :: cst

   !> Ensight postprocessing
   type(surfmesh) :: smesh,smesh_x,smesh_y,smesh_z
   type(vtk)      :: vtk_out
   type(event)    :: vtk_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,curvfile
   character(len = 1000) :: MonitorFileName

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Problem definition
   real(WP), dimension(3) :: center,vel
   real(WP) :: radius
   real(WP) :: pert_amp    !< perturbation amplitude
   integer  :: pert_mode
   integer  :: npart
   integer :: InterfaceTypes
   !> output monitor variables
   real(WP) :: r_max, r_min, time_cap, time_nd

contains

   !> computing maximum radius during growth of ligament for rayleigh pletau instability test case
   subroutine compute_rp_max_radius(vf)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_MAX, MPI_SUM, MPI_DOUBLE_PRECISION
      use mathtools, only: Pi
      implicit none

      class(vfs), intent(inout) :: vf
      integer :: i,j,k,ierr,nk
      real(WP) :: radius_k
      real(WP) :: max_radius_local, max_radius_global
      real(WP), allocatable, dimension(:) :: slice_volume_local, slice_volume_global

      allocate(slice_volume_local(vf%cfg%kmin:vf%cfg%kmax))
      allocate(slice_volume_global(vf%cfg%kmin:vf%cfg%kmax))

      slice_volume_local  = 0.0_WP
      slice_volume_global = 0.0_WP

      ! Accumulate local contribution to each global z-slice
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               slice_volume_local(k) = slice_volume_local(k) + &
                  vf%VF(i,j,k)*vf%cfg%vol(i,j,k)
            end do
         end do
      end do

      nk = vf%cfg%kmax - vf%cfg%kmin + 1

      ! Sum slice volumes entrywise across all ranks
      call MPI_ALLREDUCE(slice_volume_local, slice_volume_global, nk, &
         MPI_DOUBLE_PRECISION, MPI_SUM, vf%cfg%comm, ierr)

      max_radius_local = 0.0_WP

      do k=vf%cfg%kmin,vf%cfg%kmax
         if (slice_volume_global(k) .gt. 0.0_WP) then
            radius_k = sqrt(slice_volume_global(k)/(Pi*vf%cfg%dz(k)))
            if (radius_k .gt. max_radius_local) max_radius_local = radius_k
         end if
      end do

      call MPI_ALLREDUCE(max_radius_local, max_radius_global, 1, &
         MPI_DOUBLE_PRECISION, MPI_MAX, vf%cfg%comm, ierr)

      r_max = max_radius_global

      deallocate(slice_volume_local, slice_volume_global)
   end subroutine compute_rp_max_radius

   
   !> computing minimum radius during growth of ligament for rayleigh pletau instability test case
   subroutine compute_rp_min_radius(vf)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_MIN, MPI_SUM, MPI_DOUBLE_PRECISION
      use mathtools, only: Pi
      implicit none

      class(vfs), intent(inout) :: vf
      integer :: i,j,k,ierr,nk
      real(WP) :: dz, radius_k
      real(WP) :: min_radius_local, min_radius_global
      real(WP), allocatable :: slice_volume_local(:), slice_volume_global(:)

      allocate(slice_volume_local(vf%cfg%kmin:vf%cfg%kmax))
      allocate(slice_volume_global(vf%cfg%kmin:vf%cfg%kmax))

      slice_volume_local  = 0.0_WP
      slice_volume_global = 0.0_WP

      ! Accumulate local contribution to each global z-slice
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               slice_volume_local(k) = slice_volume_local(k) + &
                  vf%VF(i,j,k)*vf%cfg%vol(i,j,k)
            end do
         end do
      end do

      nk = vf%cfg%kmax - vf%cfg%kmin + 1

      ! Sum slice volumes entrywise across all ranks
      call MPI_ALLREDUCE(slice_volume_local, slice_volume_global, nk, &
         MPI_DOUBLE_PRECISION, MPI_SUM, vf%cfg%comm, ierr)

      min_radius_local = huge(1.0_WP)

      do k=vf%cfg%kmin,vf%cfg%kmax
         if (slice_volume_global(k) .gt. 0.0_WP) then
            radius_k = sqrt(slice_volume_global(k)/(Pi*vf%cfg%dz(k)))
            if (radius_k .lt. min_radius_local) min_radius_local = radius_k
         end if
      end do

      call MPI_ALLREDUCE(min_radius_local, min_radius_global, 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, vf%cfg%comm, ierr)

      r_min = min_radius_global

      deallocate(slice_volume_local, slice_volume_global)
   end subroutine compute_rp_min_radius


   subroutine compute_rp_stats(vf)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_MAX, MPI_SUM, MPI_DOUBLE_PRECISION,MPI_IN_PLACE,MPI_MIN
      use mathtools, only: Pi
      implicit none

      class(vfs), intent(inout) :: vf
      integer :: i,j,k,ierr
      real(WP) :: max_radius_local, max_radius_global
      real(WP) :: x,y,r
      logical :: found_interface
      real(WP), dimension(1:3) :: cen
      r_max = -huge(1.0_WP)
      r_min = huge(1.0_WP)
      found_interface=.false.

      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if(vf%VF(i,j,k) .gt. 1e-12 .and. vf%VF(i,j,k).lt.1.0_WP - 1e-12) then 
                  cen = calculateCentroid(vf%interface_polygon(1,i,j,k))
                  x = cen(1)-center(1)
                  y = cen(2)-center(2)
                  r = sqrt(x*x+y*y)

                  r_max=max(r_max,r)
                  r_min=min(r_min,r)
                  found_interface = .true. 
               endif
            end do
         end do
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE,r_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,r_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,vf%cfg%comm,ierr)
   end subroutine compute_rp_stats


   !> modified bessel functions of first kind
   function bessel_i0(x) result(val)
      implicit none
      real(WP), intent(in) :: x
      real(WP) :: val
      real(WP) :: ax, y

      ax = abs(x)

      if (ax < 3.75_WP) then
         y = (x/3.75_WP)**2
         val = 1.0_WP + y*(3.5156229_WP + y*(3.0899424_WP + y*(1.2067492_WP + &
               y*(0.2659732_WP + y*(0.0360768_WP + y*0.0045813_WP)))))
      else
         y = 3.75_WP/ax
         val = (exp(ax)/sqrt(ax)) * (0.39894228_WP + y*(0.01328592_WP + &
               y*(0.00225319_WP + y*(-0.00157565_WP + y*(0.00916281_WP + &
               y*(-0.02057706_WP + y*(0.02635537_WP + y*(-0.01647633_WP + &
               y*0.00392377_WP))))))))
      end if
   end function bessel_i0

   function bessel_i1(x) result(val)
      implicit none
      real(WP), intent(in) :: x
      real(WP) :: val
      real(WP) :: ax, y

      ax = abs(x)

      if (ax < 3.75_WP) then
         y = (x/3.75_WP)**2
         val = x*(0.5_WP + y*(0.87890594_WP + y*(0.51498869_WP + &
               y*(0.15084934_WP + y*(0.02658733_WP + y*(0.00301532_WP + &
               y*0.00032411_WP))))))
      else
         y = 3.75_WP/ax
         val = (exp(ax)/sqrt(ax)) * (0.39894228_WP + y*(-0.03988024_WP + &
               y*(-0.00362018_WP + y*(0.00163801_WP + y*(-0.01031555_WP + &
               y*(0.02282967_WP + y*(-0.02895312_WP + y*(0.01787654_WP - &
               y*0.00420059_WP))))))))
         if (x < 0.0_WP) val = -val
      end if
   end function bessel_i1

   !> Level set function for sinusoidally perturbed ligament along z-axis
   function levelset_ligament(xyz,t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP) :: xrel,yrel,zrel,r,rloc,kz

      xrel = xyz(1) - center(1)
      yrel = xyz(2) - center(2)
      zrel = xyz(3) - center(3)

      r = sqrt(xrel*xrel + yrel*yrel)

      ! Axial perturbation
      kz   = 2.0_WP*Pi*real(pert_mode,WP)/(cfg%zL)
      ! kz   = 2.0_WP*Pi
      rloc = radius*(1.0_WP + pert_amp*cos(kz*zrel))

      ! Positive inside liquid
      G = rloc - r
   end function levelset_ligament


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
         use mms_geom, only: cube_refine_vol
         use vfs_class,only: r2p,lvira,VFhi,VFlo,jibben,plicnet,jibbenPU,jibbenSq,PUppic,jibbenPLIC
         use mpi_f08,  only: MPI_WTIME
         use string,   only: str_medium,lowercase
         integer :: i,j,k,n,si,sj,sk,curvature_method,stencil_size,hf_backup_method
         character(len=str_medium) :: read_curvature_method
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=5
         real(WP) :: start, finish
         ! Create a VOF solver with r2p reconstruction
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
                  call vf%initialize(cfg=cfg,reconstruction_method=jibbenPU,name='VOF')
               case default 
                  call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
            END SELECT
         ! Reading ligament parameters
         call param_read('Ligament radius',radius)
         call param_read('Ligament center',center,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Perturbation amplitude',pert_amp)
         call param_read('Perturbation mode',pert_mode,default=1)      
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_ligament,0.0_WP,amr_ref_lvl)
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
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Perform PPIC reconstruction
         if (vf%ppic) call vf%build_quadratic_interface()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use mathtools, only: Pi
         integer :: i,j,k
         real(WP), dimension(3) :: xyz
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
         !Capilary time scale
         time_cap = sqrt(fs%rho_l*radius**3/fs%sigma)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initial velocity field
         call param_read('Initial velocity',vel,default=[0.0_WP,0.0_WP,0.0_WP])
         fs%U=vel(1)
         fs%V=vel(2)
         fs%W=vel(3)
         initialize_rp_velocity: block
            real(WP) :: kz, kr0, I0kr0, I1kr0, I0kr, I1kr
            real(WP), parameter :: tiny_r = 1.0e-14_WP
            real(WP) :: xrel, yrel, zrel, r
            real(WP) :: B, c
            real(WP), dimension(3) :: xyz

            ! wave number
            kz = 2.0_WP*Pi*real(pert_mode,WP)/cfg%zL
            kr0   = kz*radius
            I0kr0 = bessel_i0(kr0)
            I1kr0 = bessel_i1(kr0)

            ! Linear inviscid growth-rate parameter
            c = sqrt((fs%sigma/(fs%rho_l*radius**3)) * (I1kr0/I0kr0) * kr0 * (1.0_WP - kr0**2))
            !c  = sqrt(max(c,0.0_WP))

            ! amplitude of velocity potential
            B = pert_amp*radius*c/(kz*I1kr0)

            ! U velocity at x-faces
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     xyz = [fs%cfg%x(i), fs%cfg%ym(j), fs%cfg%zm(k)]

                     ! if (levelset_ligament(xyz,0.0_WP) >= 0.0_WP) then
                        xrel = xyz(1) - center(1)
                        yrel = xyz(2) - center(2)
                        zrel = xyz(3) - center(3)
                        r    = sqrt(xrel*xrel + yrel*yrel)

                        if (r > tiny_r) then
                           I1kr = bessel_i1(kz*r)
                           fs%U(i,j,k) = B*kz*I1kr*(xrel/r)*cos(kz*zrel)
                        else
                           fs%U(i,j,k) = 0.0_WP
                        ! end if
                     end if
                  end do
               end do
            end do

            ! V velocity at y-faces
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     xyz = [fs%cfg%xm(i), fs%cfg%y(j), fs%cfg%zm(k)]

                     ! if (levelset_ligament(xyz,0.0_WP) >= 0.0_WP) then
                        xrel = xyz(1) - center(1)
                        yrel = xyz(2) - center(2)
                        zrel = xyz(3) - center(3)
                        r    = sqrt(xrel*xrel + yrel*yrel)

                        if (r > tiny_r) then
                           I1kr = bessel_i1(kz*r)
                           fs%V(i,j,k) = B*kz*I1kr*(yrel/r)*cos(kz*zrel)
                        else
                           fs%V(i,j,k) = 0.0_WP
                        ! end if
                     end if
                  end do
               end do
            end do

            ! W velocity at z-faces
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     xyz = [fs%cfg%xm(i), fs%cfg%ym(j), fs%cfg%z(k)]

                     ! if (levelset_ligament(xyz,0.0_WP) >= 0.0_WP) then
                        xrel = xyz(1) - center(1)
                        yrel = xyz(2) - center(2)
                        zrel = xyz(3) - center(3)
                        r    = sqrt(xrel*xrel + yrel*yrel)

                        I0kr = bessel_i0(kz*r)
                        fs%W(i,j,k) = -B*kz*I0kr*sin(kz*zrel)
                     ! end if
                  end do
               end do
            end do

         end block initialize_rp_velocity
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
         call param_read('Curvature Option',cst%CurvatureOption,default = 1)
         call param_read('PU Spread',cst%PU_spread,default = 2.0_WP)
         call param_read('Two Dimensional',cst%TwoD,default = .false.)
      end block create_surface_tension_solver

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=2,name='plic')
         smesh%varname(1)='curv'
         smesh%varname(2)='interface_type'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         call add_surfgrid_variable(smesh,1,vf%curv,vf)  
         call add_surfgrid_variable(smesh,1,vf%interface_type,vf)    
      end block create_smesh

      create_smesh_x: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh_x=surfmesh(nvar=0,name='X_PLIC')
         ! smesh_x%varname(1)='curv'
         ! smesh_x%varname(1)='interface_type'
         ! Transfer polygons to smesh
         call cst%vf_x%update_surfmesh(smesh_x)
         ! Also populate nplane variable
         ! call add_surfgrid_variable(smesh_x,1,cst%vf_x%curv)  
         ! call add_surfgrid_variable(smesh_x,1,cst%vf_x%interface_type,cst%vf_x)    
      end block create_smesh_x
      print *,"smesh_x end"
      ! create_smesh_y: block
      !    use irl_fortran_interface
      !    integer :: i,j,k,nplane,np
      !    ! Include an extra variable for number of planes
      !    smesh_y=surfmesh(nvar=2,name='X_PLIC')
      !    smesh_y%varname(1)='curv'
      !    smesh_y%varname(2)='interface_type'
      !    ! Transfer polygons to smesh
      !    call cst%vf_y%update_surfmesh(smesh_y)
      !    ! Also populate nplane variable
      !    ! call add_surfgrid_variable(smesh_y,1,cst%vf_x%curv)  
      !    ! call add_surfgrid_variable(smesh_y,1,cst%vf_x%interface_type)    
      ! end block create_smesh_y

      ! create_smesh_z: block
      !    use irl_fortran_interface
      !    integer :: i,j,k,nplane,np
      !    ! Include an extra variable for number of planes
      !    smesh_z=surfmesh(nvar=2,name='X_PLIC')
      !    smesh_z%varname(1)='curv'
      !    smesh_z%varname(2)='interface_type'
      !    ! Transfer polygons to smesh
      !    call cst%vf_z%update_surfmesh(smesh_z)
      !    ! Also populate nplane variable
      !    ! call add_surfgrid_variable(smesh_z,1,cst%vf_x%curv)  
      !    ! call add_surfgrid_variable(smesh_z,1,cst%vf_x%interface_type)    
      ! end block create_smesh_z

      ! Add Ensight output
      create_vtk: block
         ! Create Ensight output from cfg
         vtk_out=vtk(cfg=cfg,name='RayleighPlateau')
         ! Create event for Ensight output
         vtk_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',vtk_evt%tper)
         ! Add variables to output
         call vtk_out%add_vector('velocity',Ui,Vi,Wi)
         call vtk_out%add_vector('sigma_3d_x',cst%sigma_3D(:,:,:,1,1),cst%sigma_3D(:,:,:,1,2),cst%sigma_3D(:,:,:,1,3))
         call vtk_out%add_vector('sigma_3d_y',cst%sigma_3D(:,:,:,2,1),cst%sigma_3D(:,:,:,2,2),cst%sigma_3D(:,:,:,2,3))
         call vtk_out%add_vector('sigma_3d_z',cst%sigma_3D(:,:,:,3,1),cst%sigma_3D(:,:,:,3,2),cst%sigma_3D(:,:,:,3,3))

         call vtk_out%add_vector('OLD sigma_3d_x',cst%sigma_3D_OLD(:,:,:,1,1),cst%sigma_3D_OLD(:,:,:,1,2),cst%sigma_3D_OLD(:,:,:,1,3))
         call vtk_out%add_vector('OLD sigma_3d_y',cst%sigma_3D_OLD(:,:,:,2,1),cst%sigma_3D_OLD(:,:,:,2,2),cst%sigma_3D_OLD(:,:,:,2,3))
         call vtk_out%add_vector('OLD sigma_3d_z',cst%sigma_3D_OLD(:,:,:,3,1),cst%sigma_3D_OLD(:,:,:,3,2),cst%sigma_3D_OLD(:,:,:,3,3))

         call vtk_out%add_vector('d sigma_3d_x',cst%dSigma(:,:,:,1,1),cst%dSigma(:,:,:,1,2),cst%dSigma(:,:,:,1,3))
         call vtk_out%add_vector('d sigma_3d_y',cst%dSigma(:,:,:,2,1),cst%dSigma(:,:,:,2,2),cst%dSigma(:,:,:,2,3))
         call vtk_out%add_vector('d sigma_3d_z',cst%dSigma(:,:,:,3,1),cst%dSigma(:,:,:,3,2),cst%dSigma(:,:,:,3,3))
         
         call vtk_out%add_vector('diff ST Force',cst%Fst_x_3D_diff,cst%Fst_y_3D_diff,cst%Fst_z_3D_diff)

         call vtk_out%add_vector('OLD ST Force',cst%Fst_x_3D_OLD,cst%Fst_y_3D_OLD,cst%Fst_z_3D_OLD)
         call vtk_out%add_vector('ST Force',fs%Pjx,fs%Pjy,fs%Pjz)
         call vtk_out%add_vector('CSF Force',cst%PjxD,cst%PjyD,cst%PjzD)
         call vtk_out%add_scalar('VOF',vf%VF)
         call vtk_out%add_scalar('VOFX',cst%vf_x%VF)
         call vtk_out%add_scalar('VOFY',cst%vf_y%VF)
         call vtk_out%add_scalar('VOFZ',cst%vf_z%VF)
         call vtk_out%add_scalar('curvature',vf%curv)
         call vtk_out%add_surface('plic',smesh) 
         call vtk_out%add_surface('plic_x',smesh_x) 
         ! call vtk_out%add_surface('plic_y',smesh_y) 
         ! call vtk_out%add_surface('plic_z',smesh_z) 
         ! Output to vtk
         if (vtk_evt%occurs()) call vtk_out%write_data(time%t)
      end block create_vtk

      create_monitor: block
         ! Create simulation monitor
         call param_read('Monitor File Name',MonitorFileName)
         ! Initial nondimensional monitor values
         time_nd = 0.0_WP
         r_max = 0.0_WP
         r_min = 0.0_WP
         mfile=monitor(fs%cfg%amRoot,MonitorFileName)
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time_nd,'t/tc')
         call mfile%add_column(r_max,'Rmax')
         call mfile%add_column(r_min,'Rmin')

         call mfile%write()         
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

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

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

         ! monitor variables for diagonistics
         ! call compute_rp_max_radius(vf)
         ! call compute_rp_min_radius(vf)

         call compute_rp_stats(vf)
         time_nd = time%t / time_cap
         call mfile%write()

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)

            ! Explicit calculation of drho*u/dt from NS
            call cst%get_dmomdt(resU,resV,resW)

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

         ! Output to vtk
         if (vtk_evt%occurs()) then
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               call add_surfgrid_variable(smesh,1,vf%curv,vf)   
               call add_surfgrid_variable(smesh,2,vf%interface_type,vf)   
               ! Transfer polygons to smesh
               call cst%vf_x%update_surfmesh(smesh_x)
               ! ! call add_surfgrid_variable(smesh_x,1,cst%vf_x%curv)   
               ! call add_surfgrid_variable(smesh_x,1,cst%vf_x%interface_type,cst%vf_x)   
               ! ! Transfer polygons to smesh
               ! call cst%vf_y%update_surfmesh(smesh_y)
               ! ! call add_surfgrid_variable(smesh_y,1,cst%vf_y%curv)   
               ! ! call add_surfgrid_variable(smesh_y,2,cst%vf_y%interface_type)   
               ! ! Transfer polygons to smesh
               ! call cst%vf_z%update_surfmesh(smesh_z)
               ! ! call add_surfgrid_variable(smesh_z,1,cst%vf_z%curv)   
               ! ! call add_surfgrid_variable(smesh_z,2,cst%vf_z%interface_type)   
            end block update_smesh  
            call vtk_out%write_data(time%t)
         end if

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! call print_curvature_error()
      ! Get rid of all objects - need destructors
      ! monitor
      ! vtk
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
                        smesh%polyvar(var_index,nbt+np)=A(i,j,k)
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
 
