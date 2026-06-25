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
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resH
   real(WP), dimension(:,:,:), allocatable :: resHG,resHL
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,errorMatrix
   
   !> Problem definition
   logical :: second_bubble
   real(WP), dimension(3) :: center,center2,gravity
   real(WP) :: volume,radius,radius2,Ycent,Vrise
   real(WP) :: Vin,Vin_old,Vrise_ref,Ycent_ref,G,ti
   real(WP) :: Tbot,Ttop,Ly,Lx,dtps,errorNorm,Ncurr,NperDiameter

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
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none
      
      
      ! Create a monitor file
      create_monitor: block
         ! Create simulation monitor
         mfile=monitor(cfg%amRoot,'Error Values')
         call mfile%add_column(Ncurr,'Grid Size')
         call mfile%add_column(NperDiameter,'Cell/Diam.')
         call mfile%add_column(errorNorm,'Distance Error')

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
      integer nn
      real(WP), dimension(3) :: Nset
      character(len=20) :: suffix
      
      Nset = (/32,64,128/)

      call param_read('Lx',Lx); 
      call param_read('Ly',Ly); 
      call param_read('Lz',Lz); 
      
      do nn = 1, size(Nset) 
         write(suffix,'(F0.0)') Nset(nn)
        
         ! Where we make the grid and appropriate Config File
         geometry_init : block 
            type(sgrid) :: grid
            logical :: twoD
            ! Create a grid from input params
            create_grid: block
               use sgrid_class, only: cartesian
               integer :: i,j,k,nx,ny,nz
               
               real(WP), dimension(:), allocatable :: x,y,z
               
               ! Read in grid definition
               
               nx = Nset(nn)
               ny = Nset(nn)
               nz = Nset(nn)

               call param_read('Two Dimensional',twoD)
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
         if(cfg%amRoot) then
            print *, " Geometry made for N = ",Nset(nn)
         endif


         ! Initialize time tracker with 2 subiterations
         initialize_timetracker: block
            time=timetracker(amRoot=cfg%amRoot)
            call param_read('Max timestep size',time%dtmax)
            call param_read('Max cfl number',time%cflmax)
            call param_read('Max time',time%tmax)
            time%dt=time%dtmax
            time%itmax=2
         end block initialize_timetracker
         if(cfg%amRoot) then
            print *, " Timetracker made for N = ",Nset(nn)
         endif

         ! Initialize our VOF solver and field
         allocate(vf)
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

         if(cfg%amRoot) then
            print *, " VF made for N = ",Nset(nn)
         endif

         
         allocate(fs)
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
         if(cfg%amRoot) then
            print *, " FS made for N = ",Nset(nn)
         endif
         ! ! Create Surface Tension Object

         allocate(cst)
         create_surface_tension_solver : block
            call cst%init(fs,vf,time)
            call param_read('Marangoni',cst%MarangoniOption)
            call param_read('Pressure',cst%PressureOption)
            call param_read('Smoothing Option',cst%SmoothingOption)
            call param_read('Surface Tension Option',cst%SurfaceTensionOption)
            call param_read('Curvature Option',cst%CurvatureOption)
            call param_read('PU Spread',cst%PU_spread)
         end block create_surface_tension_solver
         if(cfg%amRoot) then
            print *, " CST made for N = ",Nset(nn)
         endif

         allocate(errorMatrix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         update_distance_fields : block
            real(WP) :: xEval, yEval, zEval, dist 
            real(WP), dimension(3) :: pos,dPos
            integer :: i,j,k,count 
            call cst%updatePartitionOfUnity()
            count = 0
            ! Output True Distance
            do k=vf%cfg%kmin_,vf%cfg%kmax_
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  do i=vf%cfg%imin_,vf%cfg%imax_
                     ! write(*,'(A,3I10.5)') ' ==== location: ', i,j,k 
                     ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then  
                     ! Get Points
                     xEval = vf%cfg%xm(i)
                     yEval = vf%cfg%ym(j)
                     zEval = vf%cfg%zm(k)
                     pos = (/xEval,yEval,zEval/)

                     dPos = pos - center

                     dist = sqrt(sum(dPos**2))
                     if(abs(cst%DistanceField(i,j,k)).gt. 1e-8) then 
                        cst%trueDistanceField(i,j,k) = dist - radius
                        count = count +1
                     endif

                     errorMatrix(i,j,k) = cst%DistanceField(i,j,k) - cst%trueDistanceField(i,j,k)
                  enddo
               enddo
            enddo 
         errorNorm = sqrt(sum(errorMatrix**2)/count)

         end block update_distance_fields

               
         ! Visual Output
         ! Add Ensight output
         create_vtk: block
            ! Create Ensight output from cfg
            vtk_out=vtk(cfg=cfg,name='ConvergenceVisualization' // trim(suffix))
            ! Create event for Ensight output
            vtk_evt=event(time=time,name='Ensight output')
            call param_read('Ensight output period',vtk_evt%tper)
            ! Add variables to output
            call vtk_out%add_scalar('VOF',vf%VF)
            call vtk_out%add_scalar('curvature',vf%curv)
            call vtk_out%add_scalar('pressure',fs%P)
            call vtk_out%add_scalar('PU distance',cst%DistanceField)
            call vtk_out%add_scalar('PU Value',cst%PartitionOfUnityValue)
            call vtk_out%add_scalar('PU Weight',cst%PartitionOfUnityWeight)
            call vtk_out%add_scalar('PU Grad Mag',cst%PUTangent_mag)
            call vtk_out%add_scalar('True Distance Field',cst%trueDistanceField)
            call vtk_out%add_scalar('Error',errorMatrix)
            ! Output to vtk
            call vtk_out%write_data(time%t)

            ! Update for File
            Ncurr = Nset(nn)
            NperDiameter = Ncurr * 2 * radius / Lx
            call mfile%write()
         end block create_vtk

         ! 
         ! All objects are made. Now we can compute
         ! Kill VF and FS
         deallocate(vf)
         deallocate(fs)
         deallocate(cst)
         deallocate(errorMatrix)
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