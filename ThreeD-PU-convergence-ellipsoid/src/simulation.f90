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
   type(tpns),           public,allocatable :: fs
   type(vfs),            public,allocatable :: vf
   type(vfs),            public,allocatable :: vf_PU
   type(tads),           public :: ts
   type(timetracker),    public :: time
   type(conservative_st_type), public,allocatable :: cst

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(vtk)  :: vtk_out
   type(event)    :: vtk_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,bubblefile
   
   public :: simulation_init,simulation_run,simulation_final
   
   logical,parameter :: debug = .false.
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resH
   real(WP), dimension(:,:,:), allocatable :: resHG,resHL
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,errorMatrix_radius,errorMatrix_tangent,errorMatrix_curvature
   real(WP), dimension(:,:,:), allocatable :: errorMatrix_force_x,errorMatrix_force_y,errorMatrix_force_z
   real(WP), dimension(:,:,:,:,:), allocatable :: errorMatrix_stress
   logical :: twoD
   !> Problem definition
   logical :: second_bubble
   real(WP), dimension(3) :: center2,gravity
   real(WP) :: volume,radius,radius2,Ycent,Vrise
   real(WP) :: Vin,Vin_old,Vrise_ref,Ycent_ref,G,ti
   real(WP) :: Tbot,Ttop,Ly,Lx,dtps,Ncurr,NperDiameter,spreadCurr
   real(WP) :: errorNorm_radius, errorNorm_tangent, errorNorm_curvature,errorNorm_stress_total
   real(WP) :: errorNorm_force_x, errorNorm_force_y, errorNorm_force_z
   real(WP), dimension(3,3) :: errorNorm_stress
   ! Geometry
   real(WP), dimension(3) :: center
   real(WP) :: polarAngle, azimuthAngle,axis1size,axis2size,axis3size
   real(WP),dimension(3,3) :: A,D,Q 
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

   function levelset_PU(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), dimension(3) :: dr
      real(WP), intent(in) :: t
      real(WP) :: G

      call cst%getValuePU(xyz,G)
   end function levelset_PU
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
   
   function levelset_ellipsoid(xyz,t) result(G)
      use param, only: param_read
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), dimension(3) :: dx
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP), dimension(3) :: step1

      
      dx = xyz - center

      ! Get Value
      step1 = matmul(A,dx)
      G = 1-dot_product(dx,step1)
   end function levelset_ellipsoid

   
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
   
   subroutine random_stduniform(u)
      implicit none
      real(WP),intent(out) :: u
      real(WP) :: r
      call random_number(r)
      u = 1 - r
   end subroutine random_stduniform

   subroutine random_uniform(a,b,x)
      implicit none
      real(WP),intent(in) :: a,b
      real(WP),intent(out) :: x
      real(WP) :: u
      call random_stduniform(u)
      x = (b-a)*u + a
   end subroutine random_uniform

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none
      ! Create a monitor file
      create_monitor: block
         ! Create simulation monitor
         mfile=monitor(cfg%amRoot,'Error_Values')
         call mfile%add_column(Ncurr,'Grid Size')
         call mfile%add_column(NperDiameter,'Cell/Diam.')
         call mfile%add_column(spreadCurr,'PU Spread')
         call mfile%add_column(center(1),'Center X')
         call mfile%add_column(center(2),'Center Y')
         call mfile%add_column(center(3),'Center Z')
         call mfile%add_column(errorNorm_radius,'Radius Error')
         call mfile%add_column(errorNorm_tangent,'Tangent Error')
         call mfile%add_column(errorNorm_curvature,'Curvature Error')

         call mfile%add_column(errorNorm_stress(1,1),'Sigma_xx Error')
         call mfile%add_column(errorNorm_stress(1,2),'Sigma_xy Error')
         call mfile%add_column(errorNorm_stress(1,3),'Sigma_xz Error')

         call mfile%add_column(errorNorm_stress(2,1),'Sigma_yx Error')
         call mfile%add_column(errorNorm_stress(2,2),'Sigma_yy Error')
         call mfile%add_column(errorNorm_stress(2,3),'Sigma_yz Error')

         call mfile%add_column(errorNorm_stress(3,1),'Sigma_zx Error')
         call mfile%add_column(errorNorm_stress(3,2),'Sigma_zy Error')
         call mfile%add_column(errorNorm_stress(3,3),'Sigma_zz Error')

         call mfile%add_column(errorNorm_stress_total,'StressT Error')
         
         call mfile%add_column(errorNorm_force_x,'XForce Error')
         call mfile%add_column(errorNorm_force_Y,'YForce Error')
         call mfile%add_column(errorNorm_force_Z,'ZForce Error')

         call param_read('Two Dimensional',twoD)
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform the Convergence Test
   subroutine simulation_run
      use config_class, only: config
      use precision,    only: WP
      use param, only: param_read,param_exists
      use sgrid_class, only: sgrid
      implicit none
      
      !> Single config
      type(config) :: cfg
      real(WP) :: Lx,Ly,Lz
      integer :: nn,mm,Nmult,paramSteps
      real(WP), dimension(:), allocatable :: Nset
      real(WP), dimension(:), allocatable :: paramSet
      character(len=20) :: suffix
      real(WP) :: paramMin, paramMax
      real(WP) :: Nmin, ratio

      call param_read('Lx',Lx); 
      call param_read('Ly',Ly); 
      call param_read('Lz',Lz); 
      
      ! Grid Parameters
      call param_read('Nmin',Nmin)
      call param_read('Nmult', Nmult)
      call param_read('Ratio', ratio)

      allocate(Nset(Nmult+1))
      do nn = 1, Nmult+1
         Nset(nn) = Nmin * (ratio ** (nn-1))
      enddo
      if(cfg%amRoot) then 
         print *,"Nset = (",Nset,")"
      endif

      ! Variable Setup
      call param_read('Param Min',paramMin)
      call param_read('Param Max', paramMax)
      call param_read('Param Steps', paramSteps)

      allocate(paramSet(paramSteps+1))
      if(paramSteps .eq. 0) then 
         paramSet(1) = paramMin
      else
         do nn = 1, paramSteps+1
            paramSet(nn) = (paramMin *(paramSteps+1-nn) + paramMax *(nn-1))/paramSteps
         enddo
      endif
      if(cfg%amRoot) then 
         print *,"paramSet = (",paramSet,")"
      endif

      print*,""
      print*,""
      print*,""
      print*,""

      call param_read('Bubble position',center,default=[0.0_WP,0.0_WP,0.0_WP])
      ! Get Parameters
      call param_read('Polar Angle',polarAngle)
      call param_read('Azimuthal Angle',azimuthAngle)
      
      call param_read('Semiaxis 1',axis1size)
      call param_read('Semiaxis 2',axis2size)
      call param_read('Semiaxis 3',axis3size)

      do nn = 1, size(Nset)
      do mm = 1, size(paramSet)
         write(suffix,'(F0.0)') Nset(nn)

         ! Where we make the grid and appropriate Config File      
         if(cfg%amRoot .and. debug) then
            print *, "Making Geometry"
         endif
         geometry_init : block 
            type(sgrid) :: grid
            ! Create a grid from input params 
            create_grid: block
               use sgrid_class, only: cartesian
               integer :: i,j,k,nx,ny,nz
               
               real(WP), dimension(:), allocatable :: x,y,z
               
               ! Read in grid definition
               
               nx = Nset(nn)
               ny = Nset(nn)
               nz = Nset(nn)

               
               if(twoD) then 
                  Lz = Lx/nx 
                  nz = 1
               endif

               allocate(z(nz+1))
               allocate(y(ny+1))
               allocate(x(nx+1))
               ! Create simple rectilinear grid 
               do i=1,nx+1
                  x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
               end do
               do j=1,ny+1
                  y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
               end do
               do k=1,nz+1
                  z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
               end do
               
               ! General serial grid object - add boundary conditions
               grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='Convergence')
            end block create_grid
            
            
            ! Create a config from that grid on our entire group
            create_cfg: block
               use parallel, only: group
               integer, dimension(3) :: partition
               ! Read in partition
               call param_read('Partition',partition,short='p')
               ! Create partitioned grid
               cfg=config(grp=group,decomp=partition,grid=grid)
            end block create_cfg
            
            
            ! Create masks for this config
            create_walls: block
               cfg%VF=1.0_WP
            end block create_walls
         end block geometry_init
         


         ! Initialize time tracker with 2 subiterations  
         if(cfg%amRoot .and. debug) then
            print *, "Making Time Tracker" 
         endif
         initialize_timetracker: block
            time=timetracker(amRoot=cfg%amRoot)
            call param_read('Max timestep size',time%dtmax)
            call param_read('Max cfl number',time%cflmax)
            call param_read('Max time',time%tmax)
            time%dt=time%dtmax
            time%itmax=2
         end block initialize_timetracker
            
         
         ! Construct Rotation Matrix
         if(cfg%amRoot .and. debug) then
            print *, "Making Matrix"
         endif
         create_rotation_matrix : block 
            real(WP), dimension(3) :: axis1,axis2,axis3
            D = 0.0_WP 
            D(1,1) = 1/(axis1size*axis1size)
            D(2,2) = 1/(axis2size*axis2size)
            D(3,3) = 1/(axis3size*axis3size)
            if(twoD) then 
               D(3,3) = 0.0_WP
               azimuthAngle = 0.0_WP 
            endif
            axis1 = (/COS(polarAngle)*COS(azimuthAngle),SIN(polarAngle)*cos(azimuthAngle),-sin(azimuthAngle)/)
            axis2 = (/-SIN(polarAngle),COS(polarAngle),0.0_WP/)
            axis3 = (/COS(polarAngle)*SIN(azimuthAngle),SIN(polarAngle)*SIN(azimuthAngle),COS(azimuthAngle)/)
            Q(:,1)=axis1
            Q(:,2)=axis2
            Q(:,3)=axis3 

            ! Get A Matrix
            A = matmul(Q,matmul(D,transpose(Q)))
         end block create_rotation_matrix

         ! Initialize our VOF solver and field  
         allocate(vf)
         if(cfg%amRoot .and. debug) then
            print *, "Making VOF"
         endif
         create_and_initialize_vof: block
            use mms_geom,  only: cube_refine_vol
            use vfs_class, only: plicnet,r2p,VFhi,VFlo,remap,flux_storage,remap_storage,lvira,jibben
            use mathtools, only: Pi
            integer :: i,j,k,n,si,sj,sk
            real(WP), dimension(3,8) :: cube_vertex
            real(WP), dimension(3) :: v_cent,a_cent
            real(WP) :: vol,area
            integer, parameter :: amr_ref_lvl=4
            
            ! Create a VOF solver
            call vf%initialize(cfg=cfg,reconstruction_method=jibben,transport_method=flux_storage,name='VOF')
            !vf%cons_correct=.false.
            !vf%thin_thld_max=1.5_WP
            !vf%twoplane_thld2=0.8_WP 
            ! Initialize a bubble
            call param_read('Bubble Radius',radius);

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
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_ellipsoid,0.0_WP,amr_ref_lvl)
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
            ! Perform PPIC reconstruction
            if (vf%ppic) call vf%build_quadratic_interface()
            ! Reset moments to guarantee compatibility with interface reconstruction
            call vf%reset_volume_moments()
         end block create_and_initialize_vof

         
         allocate(fs)
         if(cfg%amRoot .and. debug) then
            print *, "Making FS"
         endif
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
            ! Calculate cell-centered velocities and divergence
            ! call fs%interp_vel(Ui,Vi,Wi)  
            ! call fs%get_div()
         end block create_and_initialize_flow_solver
         ! Create Surface Tension Object

         allocate(cst)
         if(cfg%amRoot .and. debug) then
            print *, "Making CST"
         endif
         create_surface_tension_solver : block
            call cst%init(fs,vf,time)
            call param_read('Marangoni',cst%MarangoniOption)
            call param_read('Pressure',cst%PressureOption)
            call param_read('Smoothing Option',cst%SmoothingOption)
            call param_read('Surface Tension Option',cst%SurfaceTensionOption)
            call param_read('Curvature Option',cst%CurvatureOption)
            call param_read('PU Spread',cst%PU_spread)
            cst%PU_spread = paramSet(mm)
         end block create_surface_tension_solver
         if(cfg%amRoot .and. debug) then
            print *, "Making Errors"
         endif
         allocate(errorMatrix_radius(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(errorMatrix_tangent(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(errorMatrix_curvature(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(errorMatrix_stress(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3,1:3))

         allocate(errorMatrix_force_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(errorMatrix_force_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(errorMatrix_force_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         errorMatrix_stress = 0.0_WP
         errorMatrix_radius = 0.0_WP;
         errorMatrix_tangent = 0.0_WP;
         errorMatrix_curvature =0.0_WP;
         errorMatrix_force_x = 0.0_WP
         errorMatrix_force_y = 0.0_WP
         errorMatrix_force_z = 0.0_WP
         if(cfg%amRoot .and. debug) then
            print *, "Updating Errors"
         endif
         update_distance_fields : block
            
            real(WP) :: xEval, yEval, zEval, dist, curv, curv_exact
            real(WP), dimension(3) :: pos,dPos,normal,posExact,normal_exact
            integer :: i,j,k,count 
            integer, dimension(3) :: initialIndex
            if(cfg%amRoot .and. debug) then
               print *,"Update PU"
            endif
            call cst%updatePartitionOfUnity()
            count = 0 
            
            if(cfg%amRoot .and. debug) then
               print *,"Update Stress"
            endif
            call cst%updateStresses(A,center) 
            if(cfg%amRoot .and. debug) then
               print *,"Update Force" 
            endif
            call cst%updateForces(A,center)
            ! print *, "full sum stress",sqrt(sum(cst%sigma_3D_Exact**2))
            ! Output True Distance 
            if(cfg%amRoot .and. debug) then
               print *, "Loop"
            endif
            do k=vf%cfg%kmin_,vf%cfg%kmax_
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  do i=vf%cfg%imin_,vf%cfg%imax_
                     if(vf%VF(i,j,k) .gt. 1e-12 .and. vf%VF(i,j,k) .lt. 1.0_WP - 1e-12) then
                        count = count +1

                        initialIndex = (/i,j,k/)
                        ! Get Projected Values    
                        call cst%getProjectedGeometry(initialIndex,pos,normal,curv)
                        ! Get Exact Values
                        call cst%getProjectedGeometryEllipsoid(initialIndex,A,center,posExact,normal_exact,curv_exact,pos)
                        ! Radius Error
                        dPos = pos - posExact
                        dist = sqrt(sum(dPos**2))
                        errorMatrix_radius(i,j,k) = dist

                        ! Tangent Error
                        ! Normalize Tangents
                        normal = normal / sqrt(sum(normal**2))
                        normal_exact = normal_exact / sqrt(sum(normal_exact**2))

                        errorMatrix_tangent(i,j,k) = 1.0_WP + sum(normal*normal_exact)

                        ! Curvature Error
                        errorMatrix_curvature(i,j,k) = (abs(curv)-abs(curv_exact))/abs(curv_exact)
                     endif

                     ! Stress Error
                     errorMatrix_stress(i,j,k,:,:) = cst%sigma_3D(i,j,k,:,:) - cst%sigma_3D_Exact(i,j,k,:,:)

                     ! Force Error
                     errorMatrix_force_x(i,j,k) = cst%Fst_x_3D(i,j,k) - cst%Fst_x_3D_Exact(i,j,k)
                     errorMatrix_force_y(i,j,k) = cst%Fst_y_3D(i,j,k) - cst%Fst_y_3D_Exact(i,j,k)
                     errorMatrix_force_z(i,j,k) = cst%Fst_z_3D(i,j,k) - cst%Fst_z_3D_Exact(i,j,k)
                  enddo
               enddo
            enddo 
         call vf%cfg%sync(errorMatrix_radius)
         call vf%cfg%sync(errorMatrix_tangent)
         call vf%cfg%sync(errorMatrix_curvature) 
         call vf%cfg%sync(errorMatrix_stress(:,:,:,1,1)) 
         call vf%cfg%sync(errorMatrix_stress(:,:,:,1,2))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,1,2))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,2,1))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,2,2))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,2,3))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,3,1))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,3,2))
         call vf%cfg%sync(errorMatrix_stress(:,:,:,3,3))
         call vf%cfg%sync(errorMatrix_force_x)
         call vf%cfg%sync(errorMatrix_force_y)
         call vf%cfg%sync(errorMatrix_force_z) 

         errorNorm_radius = sqrt(sum(errorMatrix_radius**2)/count)  
         errorNorm_tangent = sqrt(sum(errorMatrix_tangent**2)/count)
         errorNorm_curvature = sqrt(sum(errorMatrix_curvature**2)/count)

         errorNorm_stress = sqrt(sum(sum(sum(errorMatrix_stress**2, dim=1), dim=1), dim=1) / count)
         errorNorm_stress_total = sqrt(sum(errorNorm_stress**2)/9)


         errorNorm_force_x = sqrt(sum(errorMatrix_force_x**2)/count)
         errorNorm_force_y = sqrt(sum(errorMatrix_force_y**2)/count)
         errorNorm_force_z = sqrt(sum(errorMatrix_force_z**2)/count)

         end block update_distance_fields 

               
         ! Visual Output
         ! allocate(vf_PU)
         ! print *, "Allocated" 
         ! create_and_initialize_vof_for_PU: block
         !    use mms_geom,  only: cube_refine_vol 
         !    use vfs_class, only: plicnet,r2p,VFhi,VFlo,remap,flux_storage,remap_storage,lvira,jibben
         !    use mathtools, only: Pi
         !    integer :: i,j,k,n,si,sj,sk
         !    real(WP), dimension(3,8) :: cube_vertex
         !    real(WP), dimension(3) :: v_cent,a_cent
         !    real(WP) :: vol,area
         !    integer, parameter :: amr_ref_lvl=4
            
         !    ! Create a VOF solver
         !    call vf_PU%initialize(cfg=cfg,reconstruction_method=jibben,transport_method=flux_storage,name='VOF')
         !    ! Generate interface  
         !    i = -5
         !    j = -5
         !    k = -5
         !    do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         !       do j=vf%cfg%jmino_,vf%cfg%jmaxo_
         !          do i=vf%cfg%imino_,vf%cfg%imaxo_
         !             ! Set cube vertices
         !             n=0
         !             do sk=0,1
         !                do sj=0,1
         !                   do si=0,1
         !                      n=n+1; cube_vertex(:,n)=[vf_PU%cfg%x(i+si),vf_PU%cfg%y(j+sj),vf_PU%cfg%z(k+sk)]
         !                   end do
         !                end do
         !             end do
         !             ! Call adaptive refinement code to get volume and barycenters recursively
         !             vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
         !             ! print *, "AMRing",i,j,k-
         !             call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_PU,0.0_WP,amr_ref_lvl)
         !             ! print *, "Done"
         !             vf_PU%VF(i,j,k)=vol/vf_PU%cfg%vol(i,j,k)
         !             if (vf_PU%VF(i,j,k).ge.VFlo.and.vf_PU%VF(i,j,k).le.VFhi) then
         !                vf_PU%Lbary(:,i,j,k)=v_cent
         !                vf_PU%Gbary(:,i,j,k)=([vf_PU%cfg%xm(i),vf_PU%cfg%ym(j),vf_PU%cfg%zm(k)]-vf_PU%VF(i,j,k)*vf_PU%Lbary(:,i,j,k))/(1.0_WP-vf_PU%VF(i,j,k))
         !             else
         !                vf_PU%Lbary(:,i,j,k)=[vf_PU%cfg%xm(i),vf_PU%cfg%ym(j),vf_PU%cfg%zm(k)]
         !                vf_PU%Gbary(:,i,j,k)=[vf_PU%cfg%xm(i),vf_PU%cfg%ym(j),vf_PU%cfg%zm(k)]
         !             end if
         !          end do
         !       end do
         !    end do
         !    ! Update the band
         !    call vf_PU%update_band()
         !    ! Perform interface reconstruction from VOF field 
         !    call vf_PU%build_interface()
         !    ! Set interface planes at the boundaries
         !    call vf_PU%set_full_bcond()
         !    ! Create discontinuous polygon mesh from IRL interface
         !    call vf_PU%polygonalize_interface()
         !    ! Calculate distance from polygons
         !    call vf_PU%distance_from_polygon()
         !    ! Calculate subcell phasic volumes
         !    call vf_PU%subcell_vol()
         !    ! Calculate curvature
         !    call vf_PU%get_curvature()
         !    ! Perform PPIC reconstruction
         !    if (vf_PU%ppic) call vf_PU%build_quadratic_interface()
         !    ! Reset moments to guarantee compatibility with interface reconstruction  
         !    call vf_PU%reset_volume_moments()

         ! end block create_and_initialize_vof_for_PU

         ! create_smesh: block
         !    smesh=surfmesh(nvar=0,name='PU')
         !    call vf_PU%update_surfmesh(smesh)
         ! end block create_smesh

         ! Add Ensight output
         if(cfg%amRoot .and. debug) then
            print *, "Making VTK"
         endif
         create_vtk: block
            ! Create Ensight output from cfg
            vtk_out=vtk(cfg=cfg,name='ConvergenceVisualization' // trim(suffix))
            ! Create event for Ensight output
            vtk_evt=event(time=time,name='Ensight output')
            call param_read('Ensight output period',vtk_evt%tper)
            ! Add variables to output    
            call vtk_out%add_scalar('VOF',vf%VF)
            ! call vtk_out%add_scalar('VOF_PU',vf_PU%VF)            
            call vtk_out%add_scalar('curvature',vf%curv)
            call vtk_out%add_scalar('pressure',fs%P)
            call vtk_out%add_scalar('PU distance',cst%DistanceField)
            call vtk_out%add_scalar('PU Value',cst%PartitionOfUnityValue)
            call vtk_out%add_scalar('PU Weight',cst%PartitionOfUnityWeight)
            call vtk_out%add_scalar('PU Grad Mag',cst%PUTangent_mag)
            call vtk_out%add_scalar('Point Error',errorMatrix_radius)
            call vtk_out%add_scalar('Tangent Error',errorMatrix_tangent)
            call vtk_out%add_scalar('Curvature Error',errorMatrix_curvature)
            call vtk_out%add_vector('Sigma_x*',cst%sigma_3D(:,:,:,1,1),cst%sigma_3D(:,:,:,1,2),cst%sigma_3D(:,:,:,1,3))
            call vtk_out%add_vector('Sigma_y*',cst%sigma_3D(:,:,:,2,1),cst%sigma_3D(:,:,:,2,2),cst%sigma_3D(:,:,:,2,3))
            call vtk_out%add_vector('Sigma_z*',cst%sigma_3D(:,:,:,3,1),cst%sigma_3D(:,:,:,3,2),cst%sigma_3D(:,:,:,3,3))
            call vtk_out%add_vector('Sigma_x*_exact',cst%sigma_3D_Exact(:,:,:,1,1),cst%sigma_3D_Exact(:,:,:,1,2),cst%sigma_3D_Exact(:,:,:,1,3))
            call vtk_out%add_vector('Sigma_y*_exact',cst%sigma_3D_Exact(:,:,:,2,1),cst%sigma_3D_Exact(:,:,:,2,2),cst%sigma_3D_Exact(:,:,:,2,3))
            call vtk_out%add_vector('Sigma_z*_exact',cst%sigma_3D_Exact(:,:,:,3,1),cst%sigma_3D_Exact(:,:,:,3,2),cst%sigma_3D_Exact(:,:,:,3,3))

            call vtk_out%add_vector('errorX',errorMatrix_stress(:,:,:,1,1),errorMatrix_stress(:,:,:,1,2),errorMatrix_stress(:,:,:,1,3))
            call vtk_out%add_vector('errorY',errorMatrix_stress(:,:,:,2,1),errorMatrix_stress(:,:,:,2,2),errorMatrix_stress(:,:,:,2,3))
            call vtk_out%add_vector('errorZ',errorMatrix_stress(:,:,:,3,1),errorMatrix_stress(:,:,:,3,2),errorMatrix_stress(:,:,:,3,3))

            call vtk_out%add_vector('Force_Exact',cst%Fst_x_3D_Exact(:,:,:),cst%Fst_y_3D_Exact(:,:,:),cst%Fst_z_3D_Exact(:,:,:))
            call vtk_out%add_vector('Force',cst%Fst_x_3D(:,:,:),cst%Fst_y_3D(:,:,:),cst%Fst_z_3D(:,:,:))
            ! call vtk_out%add_surface('PU_SURF',smesh)   
            ! Output to vtk
            call vtk_out%write_data(time%t)

         ! Update for File
         Ncurr = Nset(nn)
         NperDiameter = Ncurr * 2 * axis1size / Lx
         spreadCurr = paramSet(mm)
         call mfile%write()
         end block create_vtk

         ! 
         ! All objects are made. Now we can compute
         ! Kill VF and FS
         if(cfg%amRoot .and. debug) then
            print *, "Deallocate"
         endif
         deallocate(vf)
         deallocate(fs)
         deallocate(cst)
         ! deallocate(vf_PU)

         deallocate(errorMatrix_radius)
         deallocate(errorMatrix_tangent)
         deallocate(errorMatrix_curvature)
         deallocate(errorMatrix_stress)
         deallocate(errorMatrix_force_x)
         deallocate(errorMatrix_force_y)
         deallocate(errorMatrix_force_z)
         if(cfg%amRoot) then
            print *, " COMPLETE for (N,T) = (",Nset(nn),paramSet(mm),")"
         endif
      enddo
      enddo
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
      ! deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation