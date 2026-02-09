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
   use tracer_class,      only: tracer
   use sgrid_class,       only: sgrid
   use config_class,      only: config
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
   type(tracer),      public :: pt
   type(conservative_st_type), public :: cst

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(partmesh) :: pmesh
   type(vtk)      :: vtk_out,GhostVtk_Out,PuVtk_Out
   type(event)    :: vtk_evt,Ghostvtk_evt,PuVtk_evt
   type(sgrid) :: GhostGrid,PuGrid
   type(config), public :: GhostCfg,PuCfg
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,curvfile
   
   public :: simulation_init,simulation_run,simulation_final
   integer :: dir
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: PartitionOfUnityValue,PartitionOfUnityWeight
   real(WP), dimension(:,:,:), allocatable :: PUTangent_x,PUTangent_y,PUTangent_z
   
   !> Problem definition
   real(WP), dimension(3) :: center,vel,COM,planeNormal,MarangoniOption
   real(WP) :: time_scale_factor,dtCap,Nx,Lx,delta,RootMeanSquareVelocity,amp,omega,omegaT,Nomega,lambda,a
   integer :: npart,SurfaceTensionOption,BoundaryConditionOption,CurvatureOption,ReconstructionOption,i,j,k
   integer :: SmoothingOption
   real(WP) :: PressureOption
   integer :: ierr2, myrank

   character(len = 1000) :: MonitorFileName
   real(WP), parameter :: VFlo=1.0e-10_WP                         !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo                        !< Maximum VF value considered
   
   ! Rank Data
   real(WP), dimension(:,:,:), allocatable :: Ranks
   ! Conservative Surface Tension Values
   ! These are the Surface Tension Stress values, which will be held in the pressure cell.
   real(WP) :: PU_spread
   ! real(WP), dimension(:,:,:), allocatable :: sigma_xx,sigma_yy,sigma_xy,sigma_yx 
   ! real(WP), dimension(:,:,:), allocatable :: sigma_xx_P,sigma_yy_P,sigma_xy_P,sigma_yx_P 
   ! real(WP), dimension(:,:,:), allocatable :: sigma_xx_NoP,sigma_yy_NoP,sigma_xy_NoP,sigma_yx_NoP 
   ! real(WP), dimension(:,:,:), allocatable :: force_potential_field,poisson_source
   ! These are the force values, which are held at the U and V cells, respectively.
   ! They are equivalent to Pjx and Pjy in the current code, but I have them here for clearer use.
   ! We will be designing a new 
   ! real(WP), dimension(:,:,:), allocatable :: Fst_x,Fst_y,Fst_z,PjxD,PjyD,PjzD,Pjx_ST,Pjy_ST,Pjz_ST
   ! real(WP), dimension(:,:,:), allocatable :: Pjx_Marangoni,Pjy_Marangoni,Pjz_Marangoni
   ! real(WP), dimension(:,:,:), allocatable :: grad_vf_x, grad_vf_y,grad_vf_z
   ! integer, dimension(:,:,:), allocatable :: x_smoothing_tracker,y_smoothing_tracker
   ! real(WP), dimension(:,:,:), allocatable :: SurfaceTensionDiv

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
      if (j.eq.pg%imin) isIn=.true.
   end function ym_locator


   function levelset_capwave(xyz,t) result(G)
      use param, only: param_read
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), dimension(3) :: dr
      real(WP), intent(in) :: t
      real(WP), PARAMETER :: PI = 3.141592653589793
      real(WP) :: G

      call param_read('Initial Amplitude',a)
      call param_read('Wavelength',lambda)
      
      ! write(*,'(A,3F35.10)') '   center: ', center
      G=a*COS(2*PI*xyz(1)/lambda)-xyz(2)
   end function levelset_capwave

   ! Function that finds the RMS Velocity Magnitude given U,V,W
   function Calc_RMS_Velocity() result(Vrms)
      implicit none
      integer :: i,j,k
      real(WP) :: Vrms,Umag,count
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      
      ! Get Values
      U = fs%U - vel(1)
      V = fs%V - vel(2)
      W = fs%W - vel(3)

      Vrms = 0.0_WP
      count = 0.0_WP
      k = 4
      ! Loop Over all points, sum magnitudes
      do i = 4,size(U,1)-3
         do j = 4,size(U,2)-3
            Vrms = Vrms + U(i,j,k)*U(i,j,k) + V(i,j,k)*V(i,j,k)

            count = count + 1
         end do
      end do

      Vrms = sqrt(Vrms/count)
   end function Calc_RMS_Velocity
   
   function Calc_Amplitude() result(height)
      implicit none
      integer :: i,j,k 
      real(WP), dimension(:,:,:), allocatable :: VOF
      real(WP) :: height, temp
      real(WP), PARAMETER :: PI = 3.141592653589793

      omega = fs%sigma * (2*PI/lambda)*(2*PI/lambda)*(2*PI/lambda)/(2*fs%rho_g)
      omega = sqrt(omega)
      omegaT = time%t*omega
      
      temp = 0.0
      VOF = vf%VF
      height = 0.0
      k = 1
      do i = 1,size(VOF,1)-3 ! Loop in x direction
         temp = 0.0
         do j = 1,size(VOF,2)-3 ! Loop in y direction, finding height in column
            temp = temp + VOF(i,j,k)
         enddo

         if(temp .gt. height) then ! Find maximum height
            height = temp 
         endif
      enddo
      ! write(*,'(A,F10.5)') "Height Found: ",height
      ! Now that we have maximum height, multiply by dy to get to dimensional terms
      height = height * cfg%dy(2)
      ! write(*,'(A,F10.5)') "Height Dimensional: ",height
      ! Now subtract the half domain size to go to amplitude from 0 (assuming unit domain size)
      height = height - 0.5
      ! write(*,'(A,F10.5)') "Height Ret: ",height
   end function Calc_Amplitude

   !> Function that updates the COM
   function Calc_Center_of_mass() result(CenterOfMass)
      implicit none
      integer :: i,j,k 
      real(WP), dimension(:,:,:), allocatable :: VOF
      real(WP), dimension(:,:,:,:), allocatable :: LiquidBarycenter
      real(WP) :: count,VFTotal
      real(WP), dimension(3) :: CenterOfMass,CellCenter
      

      VOF = vf%VF
      LiquidBarycenter = vf%Lbary
      VFTotal = 0
      CenterOfMass = (/0.0_WP,0.0_WP,0.0_WP/)
      count = 0
      k = 1
      ! Calculate Liquid COM
      do i = 1,size(VOF,1)-3
         do j = 1,size(VOF,2)-3
            VFTotal = VFTotal + VOF(i,j,k)
            CellCenter = [vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
            CenterOfMass = CenterOfMass + VOF(i,j,k) * CellCenter
            count = count + 1
         end do
      end do
      CenterOfMass = CenterOfMass/VFTotal
   end function Calc_Center_of_mass
   
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
         allocate(Ranks(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Conservative Surface Tension
         ! allocate(sigma_xx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_yy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_xy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_yx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(sigma_xx_P(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_yy_P(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_xy_P(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_yx_P(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(sigma_xx_NoP(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_yy_NoP(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_xy_NoP(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(sigma_yx_NoP(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(Fst_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Fst_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Fst_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(PjxD(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(PjyD(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(PjzD(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(Pjx_ST(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Pjy_ST(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Pjz_ST(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(Pjx_Marangoni(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Pjy_Marangoni(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Pjz_Marangoni(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(SurfaceTensionDiv(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         ! allocate(grad_vf_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(grad_vf_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(grad_vf_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(x_smoothing_tracker(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(y_smoothing_tracker(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(force_potential_field(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(poisson_source(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
         call param_read('Inner Time Iterations',time%itmax)
         ! write(*,'(A,I6.2,$)') "Inner Iterations: ", time%itmax
         ! STOP
      end block initialize_timetracker
      
      ! Initialize our VOF solver and field
      write(*,'(A)') '=========================== Initilize VOF'
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class,only: r2p,lvira,VFhi,VFlo,jibben,neumann,dirichlet,elvira,youngs
         !use mpi_f08,  only: MPI_WTIME
         use string,   only: str_medium,lowercase
         integer :: i,j,k,n,si,sj,sk,curvature_method,stencil_size,hf_backup_method
         character(len=str_medium) :: read_curvature_method
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=5
         real(WP) :: start, finish
         ! Create a VOF solver with lvira reconstruction
         call param_Read('Reconstruction Option',ReconstructionOption)
         SELECT CASE (ReconstructionOption)
            CASE(1)
               call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
            CASE(2)
               call vf%initialize(cfg=cfg,reconstruction_method=elvira,name='VOF')
            CASE(3)
               call vf%initialize(cfg=cfg,reconstruction_method=youngs,name='VOF')
         END SELECT
         
         ! Neumann Conditions if needed
         call param_Read('Boundary Condition Option',BoundaryConditionOption)
         call param_Read('Curvature Option',CurvatureOption)

         SELECT CASE (BoundaryConditionOption)
            CASE (2)
               call vf%add_bcond(name='left',type=neumann,locator=xp_locator,dir='+x')
               call vf%add_bcond(name='right',type=neumann,locator=xm_locator,dir='-x')
               call vf%add_bcond(name='top',type=neumann,locator=yp_locator,dir='+y')
               call vf%add_bcond(name='bottom',type=neumann,locator=ym_locator,dir='-y')
            CASE (3)
               call vf%add_bcond(name='top',type=neumann,locator=yp_locator,dir='+y')
               call vf%add_bcond(name='bottom',type=neumann,locator=ym_locator,dir='-y')
            case(4)
               call vf%add_bcond(name='left',type=neumann,locator=xp_locator,dir='+x')
               call vf%add_bcond(name='right',type=neumann,locator=xm_locator,dir='-x')
            CASE DEFAULT
               
         END SELECT
         ! Initialize cap wave
         
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_capwave,0.0_WP,amr_ref_lvl)
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
         ! Set interface planes at the boundaries?
         ! call vf%set_full_bcond()
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

         SELECT CASE (BoundaryConditionOption)
            CASE (2)
               VOF_Neumann : block
               use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
               &                                setNumberOfPlanes,setPlane,matchVolumeFraction
               real(WP), dimension(1:4) :: plane
               real(WP), dimension(1)  :: curv
               integer :: i,j,k
               type(RectCub_type) :: cell
               real(WP),dimension(5) :: hs

               call new(cell)
               ! Left Neumann Condition
               ! write(*,'(A)') '=================================== Left Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(i,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Right Neumann Condition
               write(*,'(A)') '=================================== Right Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                        call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Bottom Neumann Condition
               write(*,'(A)') '================================= Bot Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino_,vf%cfg%imaxo_
                        ! Match VOF
                        vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmin-j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Top Neumann Condition
               write(*,'(A)') '================================= Top Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino_,vf%cfg%imaxo_
                        
                        ! Match VOF
                        vf%VF(i,vf%cfg%jmax_-j+1,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        if(vf%VF(i,vf%cfg%jmin-j,k) .gt. 0.0_WP) then
                           ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmax+j,k),vf%VF(i,vf%cfg%jmax_-j+1,k)
                        endif
                        
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax_-j+1),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                        call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Create discontinuous polygon mesh from IRL interface
               call vf%polygonalize_interface()
               ! Calculate distance from polygons
               call vf%distance_from_polygon()
               ! Calculate subcell phasic volumes
               call vf%subcell_vol()
               ! Calculate curvature
               call vf%get_curvature()

               if(CurvatureOption .eq. 2) then
                  write(*,'(A)') 'Getting Curvature'
                  do k=vf%cfg%kmin_,vf%cfg%kmax_
                     do j=vf%cfg%jmin_,vf%cfg%jmax_
                        do i=vf%cfg%imin_,vf%cfg%imax_
                           if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                              ! Calculate normal component
                              planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                              if(planeNormal(1) .gt. planeNormal(2)) then
                                 dir = 1
                              else
                                 dir = 2
                              endif
                              ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                              ! Update Curvature using the new method, if the curvature is not zero
                              
                              ! write(*,'(A)')  '========= Eval Curvature'
                              ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                              call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                              vf%curv(i,j,k) = curv(1)
                              ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                           endif
                        end do
                     end do
                  end do
               endif
               ! Reset moments to guarantee compatibility with interface reconstruction
               call vf%reset_volume_moments()

            end block VOF_Neumann
            CASE (3)
               VOF_NeumannB : block
               use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
               &                                setNumberOfPlanes,setPlane,matchVolumeFraction
               real(WP), dimension(1:4) :: plane
               real(WP), dimension(1)  :: curv
               integer :: i,j,k
               type(RectCub_type) :: cell
               real(WP),dimension(5) :: hs

               call new(cell)
               !Bottom Neumann Condition
               ! write(*,'(A)') '================================= Bot Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino,vf%cfg%imaxo_
                        ! Match VOF
                        vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmin-j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Top Neumann Condition
               ! write(*,'(A)') '================================= Top Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino,vf%cfg%imaxo_
                        ! write(*,'(A,3I10.3)') 'Ghost Index: ',i,vf%cfg%jmax-j+1,k
                        ! write(*,'(A,3I10.3)') 'Inside Index: ',i,vf%cfg%jmax+j,k
                        ! Match VOF
                        vf%VF(i,vf%cfg%jmax-j+1,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax-j+1),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                        call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Left Periodic
               ! write(*,'(A)') '=================================== Left Periodic'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(i,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                     
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane 
                        ! Mirror Plane on left wall. This means nx->-nx
                        ! plane(1) = -plane(1)
                        ! write(*,'(A,4F25.5)') '   Changed Plane: ', plane
                        ! Construct the ghost cell
                        ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )
                        ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! write(*,'(A,4F25.5)') '  Cell : ',plane 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Right Periodic
               ! write(*,'(A)') '=================================== Right Periodic'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0) 
                        ! Make VF
                        vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        ! plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                        call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Create discontinuous polygon mesh from IRL interface
               call vf%polygonalize_interface()
               ! Calculate distance from polygons
               call vf%distance_from_polygon()
               ! Calculate subcell phasic volumes
               call vf%subcell_vol()
               ! Calculate curvature
               call vf%get_curvature()
               if(CurvatureOption .eq. 2) then
                  write(*,'(A)') 'Getting Curvature'
                  do k=vf%cfg%kmin_,vf%cfg%kmax_
                     do j=vf%cfg%jmin_,vf%cfg%jmax_
                        do i=vf%cfg%imin_,vf%cfg%imax_
                           if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                              ! Calculate normal component
                              planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                              if(planeNormal(1) .gt. planeNormal(2)) then
                                 dir = 1
                              else
                                 dir = 2
                              endif
                              ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                              ! Update Curvature using the new method, if the curvature is not zero
                              
                              ! write(*,'(A)')  '========= Eval Curvature'
                              ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                              call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                              vf%curv(i,j,k) = curv(1)
                              ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                           endif
                        end do
                     end do
                  end do
               endif
               ! Reset moments to guarantee compatibility with interface reconstruction
               call vf%reset_volume_moments()

            end block VOF_NeumannB

            CASE (4)
               VOF_NeumannC : block
               use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
               &                                setNumberOfPlanes,setPlane,matchVolumeFraction
               real(WP), dimension(1:4) :: plane
               real(WP), dimension(1)  :: curv
               integer :: i,j,k
               type(RectCub_type) :: cell
               real(WP),dimension(5) :: hs

               call new(cell)
               ! Left Neumann Condition
               ! write(*,'(A)') '=================================== Left Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(i,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Right Neumann Condition
               ! write(*,'(A)') '=================================== Right Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                        call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Bottom Periodic
                  ! write(*,'(A)') '================================= Bot Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  !Top Periodic
                  ! write(*,'(A)') '================================= Top Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,vf%cfg%jmax_-j+1,k) = vf%VF(i,vf%cfg%jmin-j,k)
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane n the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax_-j+1),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                           call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

               ! Create discontinuous polygon mesh from IRL interface
               call vf%polygonalize_interface()
               ! Calculate distance from polygons
               call vf%distance_from_polygon()
               ! Calculate subcell phasic volumes
               call vf%subcell_vol()
               ! Calculate curvature
               call vf%get_curvature()
               if(CurvatureOption .eq. 2) then
                  write(*,'(A)') 'Getting Curvature'
                  do k=vf%cfg%kmin_,vf%cfg%kmax_
                     do j=vf%cfg%jmin_,vf%cfg%jmax_
                        do i=vf%cfg%imin_,vf%cfg%imax_
                           if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                              ! Calculate normal component
                              planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                              if(planeNormal(1) .gt. planeNormal(2)) then
                                 dir = 1
                              else
                                 dir = 2
                              endif
                              ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                              ! Update Curvature using the new method, if the curvature is not zero
                              
                              ! write(*,'(A)')  '========= Eval Curvature'
                              ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                              call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                              vf%curv(i,j,k) = curv(1)
                              ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                           endif
                        end do
                     end do
                  end do
               endif
               ! Reset moments to guarantee compatibility with interface reconstruction
               call vf%reset_volume_moments()

            end block VOF_NeumannC
            CASE DEFAULT
               ! Periodic Boundaries
               Fully_Periodic : block
                  use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
                  &                                setNumberOfPlanes,setPlane,matchVolumeFraction
                  real(WP), dimension(1:4) :: plane
                  real(WP), dimension(1)  :: curv
                  integer :: i,j,k
                  type(RectCub_type) :: cell
                  real(WP),dimension(5) :: hs

                  call new(cell)
                  ! Left Periodic
                  ! write(*,'(A)') '=================================== Left Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                        do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                           ! Make VF
                           vf%VF(i,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                        
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane 
                           ! Mirror Plane on left wall. This means nx->-nx
                           ! plane(1) = -plane(1)
                           ! write(*,'(A,4F25.5)') '   Changed Plane: ', plane
                           ! Construct the ghost cell
                           ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )
                           ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! write(*,'(A,4F25.5)') '  Cell : ',plane 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  ! Right Periodic
                  ! write(*,'(A)') '=================================== Right Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                        do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0) 
                           ! Make VF
                           vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                           ! Mirror Plane on left wall. This means nx->-nx
                           ! plane(1) = -plane(1)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                           call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  !Bottom Periodic
                  ! write(*,'(A)') '================================= Bot Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  !Top Periodic
                  ! write(*,'(A)') '================================= Top Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,vf%cfg%jmax_-j+1,k) = vf%VF(i,vf%cfg%jmin-j,k)
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane n the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax_-j+1),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                           call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  ! Create discontinuous polygon mesh from IRL interface
                  call vf%polygonalize_interface()
                  ! Calculate distance from polygons
                  call vf%distance_from_polygon()
                  ! Calculate subcell phasic volumes
                  call vf%subcell_vol()
                  ! Calculate curvature
                  call vf%get_curvature()

                  if(CurvatureOption .eq. 2) then
                     write(*,'(A)') 'Getting Curvature'
                     do k=vf%cfg%kmin_,vf%cfg%kmax_
                        do j=vf%cfg%jmin_,vf%cfg%jmax_
                           do i=vf%cfg%imin_,vf%cfg%imax_
                              if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                                 ! Calculate normal component
                                 planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                                 if(planeNormal(1) .gt. planeNormal(2)) then
                                    dir = 1
                                 else
                                    dir = 2
                                 endif
                                 ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                                 ! Update Curvature using the new method, if the curvature is not zero
                                 
                                 ! write(*,'(A)')  '========= Eval Curvature'
                                 ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                                 call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                                 vf%curv(i,j,k) = curv(1)
                                 ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                              endif
                           end do
                        end do
                     end do
                  endif
                  ! Reset moments to guarantee compatibility with interface reconstruction
                  call vf%reset_volume_moments()

               end block Fully_Periodic
         END SELECT
      end block create_and_initialize_vof
      ! write(*,'(A)') '=========================== End Initilize VOF'


      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use mathtools, only: Pi
         use tpns_class, only : neumann,clipped_neumann,dirichlet,slip
         integer :: i,j,k
         real(WP), dimension(3) :: xyz
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Add Boundary Conditions
         SELECT CASE (BoundaryConditionOption)
            CASE (2)
               call fs%add_bcond(name='left',type=slip,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
               call fs%add_bcond(name='right',type=slip,face='x',dir=1,canCorrect=.false.,locator=xp_locator)
               call fs%add_bcond(name='top',type=slip,face='y',dir=1,canCorrect=.false.,locator=yp_locator)
               call fs%add_bcond(name='bottom',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
            CASE (3)
               call fs%add_bcond(name='top',type=slip,face='y',dir=1,canCorrect=.false.,locator=yp_locator)
               call fs%add_bcond(name='bottom',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
            CASE (4)
               call fs%add_bcond(name='left',type=slip,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
               call fs%add_bcond(name='right',type=slip,face='x',dir=1,canCorrect=.false.,locator=xp_locator)
            CASE DEFAULT
               
         END SELECT

         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         call param_read('Surface Tension Option',SurfaceTensionOption);
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         
         
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
         fs%U=0.0_WP
         fs%V=0.0_WP
         fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         call param_read('Frequency Count',Nomega)
         
         call fs%apply_bcond(time%t,time%dt)
         amp = Calc_Amplitude()
         time%tmax = Nomega/omega
      end block create_and_initialize_flow_solver

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

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='curv'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         call add_surfgrid_variable(smesh,1,vf%curv)      
      end block create_smesh
      
      ! Initialize our tracer particle tracker
      initialize_pt: block
         ! Get number of particles
         call param_read('Number of tracers',npart,default=500)
         ! Create solver
         pt=tracer(cfg=cfg,name='PT')
         ! Initialize with zero particles
         call pt%resize(0)
      end block initialize_pt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=0,nvec=2,name='tracers')
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='acceleration'
         call pt%update_partmesh(pmesh)
         do i=1,pt%np_
            pmesh%vec(:,1,i)=pt%p(i)%vel
            pmesh%vec(:,2,i)=pt%p(i)%acc
         end do
      end block create_pmesh

      
      ! Ghost Cell Output
      create_ghost_vtk : block 
         use sgrid_class, only: cartesian
         use parallel, only: group
         integer, dimension(3) :: partition
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         integer :: BoundaryConditionOption
         real(WP), dimension(:), allocatable :: x,y,z
         write(*,'(A)') 'Grid Size'
         ! Change nx,ny,nz to get Ghost Cells
         nx =cfg%nx + 6
         ny = cfg%ny + 6
         nz = cfg%nz + 6
         write(*,'(A)') 'Domain Size'
         ! Change Domain Size
         Lx = cfg%xL + cfg%dx(1)*6
         Ly = cfg%yL + cfg%dy(1)*6
         Lz = cfg%zL + cfg%dz(1)*6

         allocate(x(nx+1))
         allocate(y(ny+1))
         allocate(z(nz+1))
         write(*,'(A)') 'Fill Domain'
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz 
         end do
         write(*,'(A)') 'Make Grid'
         GhostGrid=sgrid(coord=cartesian,no=0,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='GhostDropTranslation')

         ! Now that we have a grid, create a new config
         ! Read in partition
         write(*,'(A)') 'Get Partition'
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         write(*,'(A)') 'Get Config'
         GhostCfg=config(grp=group,decomp=partition,grid=GhostGrid)

         ! Now, create the vtk output
         GhostVtk_Out = vtk(cfg=GhostCfg,name ='DropTranslation_Ghost')
         Ghostvtk_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',Ghostvtk_evt%tper)
         ! Add Variables
         call GhostVtk_Out%add_vector('velocity',Ui,Vi,Wi)
         call GhostVtk_Out%add_scalar('VOF',vf%VF)
         call GhostVtk_Out%add_scalar('curvature',vf%curv)
         call GhostVtk_Out%add_scalar('pressure',fs%P)
         call GhostVtk_Out%add_scalar('sigma_xx',cst%sigma_xx)
         call GhostVtk_Out%add_scalar('sigma_xy',cst%sigma_xy)
         call GhostVtk_Out%add_scalar('sigma_yx',cst%sigma_yx)
         call GhostVtk_Out%add_scalar('sigma_yy',cst%sigma_yy)
         call GhostVtk_Out%add_surface('plic',smesh) 
         if (Ghostvtk_evt%occurs()) call GhostVtk_Out%write_data(time%t)
      end block create_ghost_vtk

      create_PU_vtk : block 
         use sgrid_class, only: cartesian
         use parallel, only: group
         integer, dimension(3) :: partition
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,Scale
         integer :: BoundaryConditionOption
         real(WP), dimension(:), allocatable :: x,y,z
         write(*,'(A)') 'Grid Size'
         ! Change nx,ny,nz to get Ghost Cells
         call param_read('Partition Of Unity Scale',Scale)
         nx = Scale*cfg%nx
         ny = Scale*cfg%ny 
         nz = 1
         write(*,'(A)') 'Domain Size'
         ! Change Domain Size
         Lx = cfg%xL
         Ly = cfg%yL
         Lz = cfg%zL
         ! Allocate Space
         allocate(PartitionOfUnityValue  (1:nx,1:ny,1:nz))
         allocate(PartitionOfUnityWeight  (1:nx,1:ny,1:nz))
         allocate(PUTangent_x  (1:nx,1:ny,1:nz))
         allocate(PUTangent_y  (1:nx,1:ny,1:nz))
         allocate(PUTangent_z  (1:nx,1:ny,1:nz))
         PartitionOfUnityValue = 0.0_WP
         PartitionOfUnityWeight = 0.0_WP
         PUTangent_x = 0.0_WP
         PUTangent_y = 0.0_WP
         PUTangent_z = 0.0_WP
         allocate(x(nx+1))
         allocate(y(ny+1))
         allocate(z(nz+1))
         write(*,'(A)') 'Fill Domain'
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz 
         end do
         write(*,'(A)') 'Make Grid'
         PuGrid=sgrid(coord=cartesian,no=0,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='DropTranslation_PU')

         ! Now that we have a grid, create a new config
         ! Read in partition
         write(*,'(A)') 'Get Partition'
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         write(*,'(A)') 'Get Config'
         PuCfg=config(grp=group,decomp=partition,grid=PuGrid)

         ! Now, create the vtk output
         PuVtk_Out = vtk(cfg=PuCfg,name ='DropTranslation_PartitionOfUnity')
         PuVtk_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',PuVtk_evt%tper)
         call param_read('PU Spread',PU_spread)
         ! Add Variables
         call PuVtk_Out%add_scalar('Partition Of Unity Value',PartitionOfUnityValue)
         call PuVtk_Out%add_scalar('Partition Of Unity Weight',PartitionOfUnityWeight)
         call PuVTK_Out%add_vector('Partition Of Unity Tangent',PUTangent_x,PUTangent_y,PUTangent_z)
         if (PuVtk_evt%occurs()) call PuVtk_Out%write_data(time%t)
      end block create_PU_vtk
      
      ! Add Ensight output
      create_vtk: block
         ! Create Ensight output from cfg
         vtk_out=vtk(cfg=cfg,name='DropTranslation')
         ! Create event for Ensight output
         vtk_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',vtk_evt%tper)
         ! Add variables to output
         call vtk_out%add_vector('velocity',Ui,Vi,Wi)
         call vtk_out%add_scalar('VOF',vf%VF)
         call vtk_out%add_scalar('curvature',vf%curv)
         call vtk_out%add_scalar('pressure',fs%P)
         call vtk_out%add_scalar('sigma_xx',cst%sigma_xx)
         call vtk_out%add_scalar('sigma_xy',cst%sigma_xy)
         call vtk_out%add_scalar('sigma_yx',cst%sigma_yx)
         call vtk_out%add_scalar('sigma_yy',cst%sigma_yy)
         call vtk_out%add_scalar('Fst_x',cst%Fst_x)
         call vtk_out%add_scalar('Fst_y',cst%Fst_y)
         call vtk_out%add_scalar('Fst_z',cst%Fst_z)
         call vtk_out%add_scalar('PjxD',cst%PjxD)
         call vtk_out%add_scalar('PjyD',cst%PjyD)
         call vtk_out%add_scalar('PjzD',cst%PjzD)
         call vtk_out%add_scalar('Pjx_ST',cst%Pjx_ST)
         call vtk_out%add_scalar('Pjy_ST',cst%Pjy_ST)
         call vtk_out%add_scalar('Pjz_ST',cst%Pjz_ST)
         call vtk_out%add_scalar('sigma_xx_P',cst%sigma_xx_P)
         call vtk_out%add_scalar('sigma_xy_P',cst%sigma_xy_P)
         call vtk_out%add_scalar('sigma_yx_P',cst%sigma_yx_P)
         call vtk_out%add_scalar('sigma_yy_P',cst%sigma_yy_P)

         call vtk_out%add_scalar('sigma_xx_NoP',cst%sigma_xx_NoP)
         call vtk_out%add_scalar('sigma_xy_NoP',cst%sigma_xy_NoP)
         call vtk_out%add_scalar('sigma_yx_NoP',cst%sigma_yx_NoP)
         call vtk_out%add_scalar('sigma_yy_NoP',cst%sigma_yy_NoP)

         call vtk_out%add_scalar('Force Potential',cst%force_potential_field)
         call vtk_out%add_scalar('Force Potential Source',cst%poisson_source)

         call vtk_out%add_scalar('grad_vf_x',cst%grad_vf_x)
         call vtk_out%add_scalar('grad_vf_y',cst%grad_vf_y)
         call vtk_out%add_scalar('grad_vf_z',cst%grad_vf_z)

         call vtk_out%add_vector('Fst',cst%Fst_x,cst%Fst_y,cst%Fst_z)
         call vtk_out%add_surface('plic',smesh) 
         call vtk_out%add_particle('tracers',pmesh)
         ! Output to vtk
         if (vtk_evt%occurs()) call vtk_out%write_data(time%t)
      end block create_vtk
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! call print_curvature_error()
         call vf%get_max()
         ! Create simulation monitor
         call param_read('Monitor File Name',MonitorFileName)
         mfile=monitor(fs%cfg%amRoot,MonitorFileName)
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(amp,'Viscous Time Scale')
         call mfile%add_column(amp,'Capillary Time Scale')
         call mfile%add_column(RootMeanSquareVelocity,'Vrms')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(COM(1),'COM X')
         call mfile%add_column(COM(2),'COM Y')
         call mfile%add_column(COM(3),'COM Z')
         call mfile%add_column(amp,'Wave Amplitude')
         call mfile%add_column(omega,'Wave Frequency')
         call mfile%add_column(omegaT,'Wave Time')
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
         ! mixedYoungCentered: block 
         !    use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
         !       &                                setNumberOfPlanes,setPlane,matchVolumeFraction
         !    real(WP), dimension(1:4) :: plane
         !    real(WP), dimension(1:3) :: norm
         !    type(RectCub_type) :: cell
         !    ! Here we do the mixed-Young-Centered Method
         !    call new(cell)
         !    do k = vf%cfg%kmin_,vf%cfg%kmax_
         !       do j = vf%cfg%jmin_,vf%cfg%jmax_
         !          do i = vf%cfg%imin_,vf%cfg%imax_
         !             if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then
         !                ! Make Cell
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
         !                ! Get Normal
         !                call mixedYoungCenterNormal(vf,i,j,k,norm)
         !                ! Update Plane
         !                plane(1:3) = norm 
         !                ! Start at center
         !                plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                ! Set Number of Plnes and Plane
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
         !                ! Match Volume Fraction
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             endif
         !          enddo
         !       enddo
         !    enddo
         ! end block mixedYoungCentered
         
         ! Apply Boundary Conditions

         SELECT CASE (BoundaryConditionOption)
            CASE (2)
               VOF_Neumann : block
               use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
               &                                setNumberOfPlanes,setPlane,matchVolumeFraction
               real(WP), dimension(1:4) :: plane
               real(WP), dimension(1)  :: curv
               integer :: i,j,k
               type(RectCub_type) :: cell
               real(WP),dimension(5) :: hs

               call new(cell)
               ! Left Neumann Condition
               write(*,'(A)') '=================================== Left Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(i,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Right Neumann Condition
               write(*,'(A)') '=================================== Right Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                        call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Bottom Neumann Condition
               write(*,'(A)') '================================= Bot Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino_,vf%cfg%imaxo_
                        ! Match VOF
                        vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmin-j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Top Neumann Condition
               write(*,'(A)') '================================= Top Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino_,vf%cfg%imaxo_
                        
                        ! Match VOF
                        vf%VF(i,vf%cfg%jmax_-j+1,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        if(vf%VF(i,vf%cfg%jmin-j,k) .gt. 0.0_WP) then
                           ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmax+j,k),vf%VF(i,vf%cfg%jmax_-j+1,k)
                        endif
                        
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax_-j+1),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                        call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Create discontinuous polygon mesh from IRL interface
               call vf%polygonalize_interface()
               ! Calculate distance from polygons
               call vf%distance_from_polygon()
               ! Calculate subcell phasic volumes
               call vf%subcell_vol()
               ! Calculate curvature
               call vf%get_curvature()

               if(CurvatureOption .eq. 2) then
                  write(*,'(A)') 'Getting Curvature'
                  do k=vf%cfg%kmin_,vf%cfg%kmax_
                     do j=vf%cfg%jmin_,vf%cfg%jmax_
                        do i=vf%cfg%imin_,vf%cfg%imax_
                           if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                              ! Calculate normal component
                              planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                              if(planeNormal(1) .gt. planeNormal(2)) then
                                 dir = 1
                              else
                                 dir = 2
                              endif
                              ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                              ! Update Curvature using the new method, if the curvature is not zero
                              
                              ! write(*,'(A)')  '========= Eval Curvature'
                              ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                              call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                              vf%curv(i,j,k) = curv(1)
                              ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                           endif
                        end do
                     end do
                  end do
               endif
               ! Reset moments to guarantee compatibility with interface reconstruction
               call vf%reset_volume_moments()

            end block VOF_Neumann
            CASE (3)
               VOF_NeumannB : block
               use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
               &                                setNumberOfPlanes,setPlane,matchVolumeFraction
               real(WP), dimension(1:4) :: plane
               real(WP), dimension(1)  :: curv
               integer :: i,j,k
               type(RectCub_type) :: cell
               real(WP),dimension(5) :: hs

               call new(cell)
               !Bottom Neumann Condition
               ! write(*,'(A)') '================================= Bot Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino,vf%cfg%imaxo_
                        ! Match VOF
                        vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmin-j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Top Neumann Condition
               ! write(*,'(A)') '================================= Top Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmin-1
                     do i = vf%cfg%imino,vf%cfg%imaxo_
                        ! write(*,'(A,3I10.3)') 'Ghost Index: ',i,vf%cfg%jmax-j+1,k
                        ! write(*,'(A,3I10.3)') 'Inside Index: ',i,vf%cfg%jmax+j,k
                        ! Match VOF
                        vf%VF(i,vf%cfg%jmax-j+1,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                        ! Mirror Plane on left wall. This means ny->-ny
                        plane(2) = -plane(2)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax-j+1),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                        ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                        call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Left Periodic
               ! write(*,'(A)') '=================================== Left Periodic'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(i,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                     
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane 
                        ! Mirror Plane on left wall. This means nx->-nx
                        ! plane(1) = -plane(1)
                        ! write(*,'(A,4F25.5)') '   Changed Plane: ', plane
                        ! ! Construct the ghost cell
                        ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )
                        ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! write(*,'(A,4F25.5)') '  Cell : ',plane 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Right Periodic
               ! write(*,'(A)') '=================================== Right Periodic'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0) 
                        ! Make VF
                        vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        ! plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                        call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Create discontinuous polygon mesh from IRL interface
               call vf%polygonalize_interface()
               ! Calculate distance from polygons
               call vf%distance_from_polygon()
               ! Calculate subcell phasic volumes
               call vf%subcell_vol()
               ! Calculate curvature
               call vf%get_curvature()
               if(CurvatureOption .eq. 2) then
                  write(*,'(A)') 'Getting Curvature'
                  do k=vf%cfg%kmin_,vf%cfg%kmax_
                     do j=vf%cfg%jmin_,vf%cfg%jmax_
                        do i=vf%cfg%imin_,vf%cfg%imax_
                           if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                              ! Calculate normal component
                              planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                              if(planeNormal(1) .gt. planeNormal(2)) then
                                 dir = 1
                              else
                                 dir = 2
                              endif
                              ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                              ! Update Curvature using the new method, if the curvature is not zero
                              
                              ! write(*,'(A)')  '========= Eval Curvature'
                              ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                              call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                              vf%curv(i,j,k) = curv(1)
                              ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                           endif
                        end do
                     end do
                  end do
               endif
               ! Reset moments to guarantee compatibility with interface reconstruction
               call vf%reset_volume_moments()

            end block VOF_NeumannB

            CASE (4)
               VOF_NeumannC : block
               use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
               &                                setNumberOfPlanes,setPlane,matchVolumeFraction
               real(WP), dimension(1:4) :: plane
               real(WP), dimension(1)  :: curv
               integer :: i,j,k
               type(RectCub_type) :: cell
               real(WP),dimension(5) :: hs

               call new(cell)
               ! Left Neumann Condition
               ! write(*,'(A)') '=================================== Left Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(i,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               ! Right Neumann Condition
               ! write(*,'(A)') '=================================== Right Neumann'
               do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i = vf%cfg%imino,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                        ! Make VF
                        vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                        ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                        ! Get the Plane in the corresponding cell within the domain
                        plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                        ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                        ! Mirror Plane on left wall. This means nx->-nx
                        plane(1) = -plane(1)
                        ! Construct the ghost cell
                        call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                        ! Reset Plane Distance to match center of cell
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                        ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                        ! Set Plane
                        call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                        call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                        ! Match Volume Fraction
                        call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                        ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                     enddo
                  enddo
               enddo

               !Bottom Periodic
                  ! write(*,'(A)') '================================= Bot Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  !Top Periodic
                  ! write(*,'(A)') '================================= Top Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,vf%cfg%jmax_-j+1,k) = vf%VF(i,vf%cfg%jmin-j,k)
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane n the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax_-j+1),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                           call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

               ! Create discontinuous polygon mesh from IRL interface
               call vf%polygonalize_interface()
               ! Calculate distance from polygons
               call vf%distance_from_polygon()
               ! Calculate subcell phasic volumes
               call vf%subcell_vol()
               ! Calculate curvature
               call vf%get_curvature()
               if(CurvatureOption .eq. 2) then
                  write(*,'(A)') 'Getting Curvature'
                  do k=vf%cfg%kmin_,vf%cfg%kmax_
                     do j=vf%cfg%jmin_,vf%cfg%jmax_
                        do i=vf%cfg%imin_,vf%cfg%imax_
                           if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                              ! Calculate normal component
                              planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                              if(planeNormal(1) .gt. planeNormal(2)) then
                                 dir = 1
                              else
                                 dir = 2
                              endif
                              ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                              ! Update Curvature using the new method, if the curvature is not zero
                              
                              ! write(*,'(A)')  '========= Eval Curvature'
                              ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                              call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                              vf%curv(i,j,k) = curv(1)
                              ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                           endif
                        end do
                     end do
                  end do
               endif
               ! Reset moments to guarantee compatibility with interface reconstruction
               call vf%reset_volume_moments()

            end block VOF_NeumannC
            CASE DEFAULT
               ! Periodic Boundaries
               Fully_Periodic : block
                  use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
                  &                                setNumberOfPlanes,setPlane,matchVolumeFraction
                  real(WP), dimension(1:4) :: plane
                  real(WP), dimension(1)  :: curv
                  integer :: i,j,k
                  type(RectCub_type) :: cell
                  real(WP),dimension(5) :: hs

                  call new(cell)
                  ! Left Periodic
                  ! write(*,'(A)') '=================================== Left Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                        do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0)
                           ! Make VF
                           vf%VF(i,j,k) = vf%VF(vf%cfg%imax+i,j,k)
                        
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(vf%cfg%imax+i,j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane 
                           ! Mirror Plane on left wall. This means nx->-nx
                           ! plane(1) = -plane(1)
                           ! write(*,'(A,4F25.5)') '   Changed Plane: ', plane
                           ! Construct the ghost cell
                           ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )
                           ! write(*,'(A,3F25.5)') '   x,y,z: ', vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! write(*,'(A,4F25.5)') '  Cell : ',plane 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  ! Right Periodic
                  ! write(*,'(A)') '=================================== Right Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmaxo_
                        do i = vf%cfg%imino_,vf%cfg%imin-1 ! This makes it loop through only the left ghost cells (Indicies -2,-1,0) 
                           ! Make VF
                           vf%VF(vf%cfg%imax-i+1,j,k) = vf%VF(vf%cfg%imin-i,j,k)
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(vf%cfg%imin-i,j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane
                           ! Mirror Plane on left wall. This means nx->-nx
                           ! plane(1) = -plane(1)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(vf%cfg%imax-i+1),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(vf%cfg%imax-i+2),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(vf%cfg%imax-i+1),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),1)
                           call setPlane(vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(vf%cfg%imax-i+1,j,k),vf%liquid_gas_interface(vf%cfg%imax-i+1,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  !Bottom Periodic
                  ! write(*,'(A)') '================================= Bot Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,j,k) = vf%VF(i,vf%cfg%jmax+j,k)
                        
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! write(*,'(A,2F25.5)') '   VFinside, VFghost: ', vf%VF(i,vf%cfg%jmin-j,k),vf%VF(i,j,k)
                           ! Get the Plane in the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmax+j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  !Top Periodic
                  ! write(*,'(A)') '================================= Top Periodic'
                  do k = vf%cfg%kmino_,vf%cfg%kmaxo_
                     do j = vf%cfg%jmino_,vf%cfg%jmin-1
                        do i = vf%cfg%imino_,vf%cfg%imaxo_
                           ! Match VOF
                           vf%VF(i,vf%cfg%jmax_-j+1,k) = vf%VF(i,vf%cfg%jmin-j,k)
                           ! write(*,'(A,4I25.5)') '   i,j,k,imin: ', i,j,k,vf%cfg%imin
                           ! Get the Plane n the corresponding cell within the domain
                           plane = getPlane(vf%liquid_gas_interface(i,vf%cfg%jmin-j,k),0) 
                           ! write(*,'(A,4F25.5)') '   Init Plane: ', plane(1),plane(2),plane(3),plane(4)
                           ! Mirror Plane on left wall. This means ny->-ny
                           ! plane(2) = -plane(2)
                           ! Construct the ghost cell
                           call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(vf%cfg%jmax_-j+1),vf%cfg%z(k  )],&
                           &                       [vf%cfg%x(i+1),vf%cfg%y(vf%cfg%jmax-j+2),vf%cfg%z(k+1)]) 
                           ! Reset Plane Distance to match center of cell
                           plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(vf%cfg%jmax-j+1),vf%cfg%zm(k)])
                           ! write(*,'(A,F25.5)') '   New  Plane4: ', plane(4)
                           ! write(*,'(A,4F25.5)') '   New  Plane: ', plane
                           ! Set Plane
                           call setNumberOfPlanes(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),1)
                           call setPlane(vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k),0,plane(1:3),plane(4))
                           ! Match Volume Fraction
                           call matchVolumeFraction(cell,vf%VF(i,vf%cfg%jmax-j+1,k),vf%liquid_gas_interface(i,vf%cfg%jmax-j+1,k))
                           ! write(*,'(A,4F25.5)') '   Fin  Plane: ', plane
                        enddo
                     enddo
                  enddo

                  ! Create discontinuous polygon mesh from IRL interface
                  call vf%polygonalize_interface()
                  ! Calculate distance from polygons
                  call vf%distance_from_polygon()
                  ! Calculate subcell phasic volumes
                  call vf%subcell_vol()
                  ! Calculate curvature
                  call vf%get_curvature()

                  if(CurvatureOption .eq. 2) then
                     write(*,'(A)') 'Getting Curvature'
                     do k=vf%cfg%kmin_,vf%cfg%kmax_
                        do j=vf%cfg%jmin_,vf%cfg%jmax_
                           do i=vf%cfg%imin_,vf%cfg%imax_
                              if(vf%VF(i,j,k) .gt. VFlo .and. vf%VF(i,j,k) .lt. VFhi) then
                                 ! Calculate normal component
                                 planeNormal = calculateNormal(vf%interface_polygon(1,i,j,k))
                                 if(planeNormal(1) .gt. planeNormal(2)) then
                                    dir = 1
                                 else
                                    dir = 2
                                 endif
                                 ! write(*,'(A,I25.10)')  '==================================== Direction: ',dir
                                 ! Update Curvature using the new method, if the curvature is not zero
                                 
                                 ! write(*,'(A)')  '========= Eval Curvature'
                                 ! write(*,'(A,F25.10)')  '==================================== Regular Curvature: ',vf%curv(i,j,k)
                                 call heightFunctionCurvature(i,j,k,vf,curv,dir,hs)
                                 vf%curv(i,j,k) = curv(1)
                                 ! write(*,'(A,F25.10)')  '==================================== Height Function Curvature: ',vf%curv(i,j,k)
                              endif
                           end do
                        end do
                     end do
                  endif
                  ! Reset moments to guarantee compatibility with interface reconstruction
                  call vf%reset_volume_moments()

               end block Fully_Periodic
         END SELECT

         ! Advance and project tracer particles
         ! call pt%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W)
         ! call project_tracers(vf=vf,pt=pt)

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
            ! call fs%get_dmomdt(resU,resV,resW)
            ! write(*,'(A)') 'CASE SELECTION ================================================'
            
            ! SELECT CASE (SurfaceTensionOption)
            !    CASE (2)
            !       ! write(*,'(A)') 'CASE 2 ================================================'
            !       call add_surface_tension_jump_no_ST(fs,dt=time%dt,div=fs%div,vf=vf)
            !       call get_dmomdt_full_ST(fs,resU,resV,resW)
            !    CASE (3)
            !       ! write(*,'(A)') 'CASE 3 ================================================'
            !       call get_dmomdt_no_ST(fs,resU,resV,resW)
            !    CASE (4)
            !       ! Conservative Surface Tension Approach for get_dmomdt
            !       call fs%get_dmomdt(resU,resV,resW)
            !    CASE DEFAULT
            !       ! write(*,'(A)') 'CASE DEFAULT ================================================'
            !       call fs%get_dmomdt(resU,resV,resW)
            ! END SELECT
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
            
            ! write(*,'(A)') '============================================ FLOW SOLVER BC'
            call fs%apply_bcond(time%t,time%dt)
            ! STOP
            ! Solve Poisson equation
            call fs%update_laplacian(pinpoint=[fs%cfg%imin,fs%cfg%jmin,fs%cfg%kmin])

            call fs%correct_mfr()
            call fs%get_div()
            ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            ! SELECT CASE (SurfaceTensionOption)
            !    CASE (2)
            !       call add_surface_tension_jump_no_ST(fs,dt=time%dt,div=fs%div,vf=vf)
            !    CASE (3)
            !       call add_surface_tension_jump_full_ST(fs,dt=time%dt,div=fs%div,vf=vf)
            !    CASE (4) 
            !       call add_conservative_surface_tension_jump(fs,dt=time%dt,div=fs%div,vf=vf)
            !    CASE DEFAULT
            !       call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            ! END SELECT
            call cst%add_surface_tension_jump()

            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            if (cfg%amRoot) fs%psolv%rhs(cfg%imin,cfg%jmin,cfg%kmin)=0.0_WP
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
               call add_surfgrid_variable(smesh,1,vf%curv)      
            end block update_smesh  
            ! Update partmesh object
            update_pmesh: block
              integer :: i
              call pt%update_partmesh(pmesh)
              do i=1,pt%np_
                  pmesh%vec(:,1,i)=pt%p(i)%vel
                  pmesh%vec(:,2,i)=pt%p(i)%acc
              end do
            end block update_pmesh
            call updatePartitionOfUnity(vf,fs)
            call vtk_out%write_data(time%t)
            call GhostVtk_Out%write_data(time%t)
            call PuVtk_Out%write_data(time%t)
         end if
         
         RootMeanSquareVelocity = Calc_RMS_Velocity()
         amp = Calc_Amplitude()
         COM = Calc_Center_of_mass()
         call mfile%write()
         ! write (*,'(A,F10.5)') '======== Surface Tension Coefficient: ', fs%sigma
         ! STOP
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
      ! ,sigma_xx,sigma_yy,sigma_xy,sigma_yx,Fst_x,Fst_y,Fst_z,PjxD,PjyD,PjzD,SurfaceTensionDiv
      deallocate(resU,resV,resW,Ui,Vi,Wi) 
   end subroutine simulation_final



   ! ################################################################## HELPER FUNCTIONS ###########################################################################
   !> Make a surface scalar variable
   subroutine add_surfgrid_variable(smesh,var_index,A)
      use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
      use vfs_class,only: VFhi,VFlo
      implicit none
      class(surfmesh), intent(inout) :: smesh
      integer, intent(in) :: var_index
      real(WP), dimension(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_), intent(in) :: A 
      integer :: i,j,k,n,shape,np,nplane,nbt
      
      ! Fill out arrays
      if ((smesh%nPoly+smesh%nBezierTri).gt.0) then
         np=0; nbt = 0
         ! Start with quadratic surfaces 
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  shape=getNumberOfTriangles(vf%interface_mixed_surface(i,j,k))
                  if (shape.gt.0) then
                     do n=1,shape
                        nbt=nbt+1
                        smesh%var(var_index,nbt)=A(i,j,k)
                     end do
                  end if
               end do
            end do
         end do
         ! Then do planes
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(vf%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        ! Set nplane variable
                        smesh%var(var_index,nbt+np)=A(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      else
         smesh%var(var_index,1)=1
      end if      
      
   end subroutine add_surfgrid_variable

   ! Project tracers on the interface
   subroutine project_tracers(vf,pt)
      use mathtools, only: normalize
      implicit none
      class(vfs), intent(inout) :: vf
      class(tracer), intent(inout) :: pt
      integer :: i, j, maxit
      real(WP) :: distx, disty, distz, dist, maxdist
      real(WP), dimension(3) :: oldvel, oldpos, normal
      ! Maximum projection iterations
      maxit=10
      ! Maximum distance from the interface wanted
      maxdist=1.0E-9_WP*vf%cfg%meshsize(0,0,0)
      ! Get gradient of distance function to the interface 
      call fs%get_pgrad(vf%G,resU,resV,resW)
      ! Loop over local tracers
      do j=1,maxit
         do i=1,pt%np_
            ! Avoid particles with id=0
            if (pt%p(i)%id.eq.0) cycle
            oldpos=pt%p(i)%pos
            ! Interpolate distance function
            dist=cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=vf%G,bc='n')
            if (abs(dist).gt.maxdist) then
               ! Interpolate interface normal
               normal=normalize([cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resU,bc='n'),&
               &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resV,bc='n'),&
               &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resW,bc='n')])
               ! Move tracer along distance function gradient
               pt%p(i)%pos=pt%p(i)%pos-dist*normal
               ! Correct the position to take into account periodicity
               pt%p(i)%pos(1)=cfg%x(cfg%imin)+modulo(pt%p(i)%pos(1)-cfg%x(cfg%imin),cfg%xL)
               pt%p(i)%pos(2)=cfg%y(cfg%jmin)+modulo(pt%p(i)%pos(2)-cfg%y(cfg%jmin),cfg%yL)
               pt%p(i)%pos(3)=cfg%z(cfg%kmin)+modulo(pt%p(i)%pos(3)-cfg%z(cfg%kmin),cfg%zL)
               ! Relocalize the tracer
               pt%p(i)%ind=cfg%get_ijk_global(pt%p(i)%pos,pt%p(i)%ind)
            end if
         end do
         ! Communicate tracers across procs
         call pt%sync()
      end do
      do i=1,pt%np_
         ! Interpolate fluid quantities to particle location
         pt%p(i)%vel=cfg%get_velocity(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
      end do
   end subroutine project_tracers
   
   ! ################################################################## MOMENTUM EQUATION REPLACEMENT SUBROUTINES ##################################################################

   ! subroutine get_dmomdt_full_ST(this,drhoUdt,drhoVdt,drhoWdt)
   !    implicit none
   !    class(tpns), intent(inout) :: this
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    integer :: i,j,k,ii,jj,kk

   !    ! Zero out drhoUVW/dt arrays
   !    drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
   !    ! Flux of rhoU
   !    do kk=this%cfg%kmin_,this%cfg%kmax_+1
   !       do jj=this%cfg%jmin_,this%cfg%jmax_+1
   !          do ii=this%cfg%imin_,this%cfg%imax_+1
   !             ! Fluxes on x-face
   !             i=ii-1; j=jj-1; k=kk-1
   !             this%FUX(i,j,k)=this%FUX(i,j,k)-this%P(i,j,k) &
   !             &              -sum(this%hybu_x(:,i,j,k)*this%rho_Uold(i:i+1,j,k))*sum(this%hybu_x(:,i,j,k)*this%U(i:i+1,j,k))*sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k)) &
   !             &              +this%visc   (i,j,k)*(sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))+sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k)))
   !             ! Fluxes on y-face
   !             i=ii; j=jj; k=kk
   !             this%FUY(i,j,k)=this%FUY(i,j,k) &
   !             &              -sum(this%hybu_y(:,i,j,k)*this%rho_Uold(i,j-1:j,k))*sum(this%hybu_y(:,i,j,k)*this%U(i,j-1:j,k))*sum(this%itpv_x(:,i,j,k)*this%V(i-1:i,j,k)) &
   !             &              +this%visc_xy(i,j,k)*(sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))+sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k)))
   !             ! Fluxes on z-face
   !             i=ii; j=jj; k=kk
   !             this%FUZ(i,j,k)=this%FUZ(i,j,k) &
   !             &              -sum(this%hybu_z(:,i,j,k)*this%rho_Uold(i,j,k-1:k))*sum(this%hybu_z(:,i,j,k)*this%U(i,j,k-1:k))*sum(this%itpw_x(:,i,j,k)*this%W(i-1:i,j,k)) &
   !             &              +this%visc_zx(i,j,k)*(sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))+sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k)))
   !          end do
   !       end do
   !    end do
   !    ! Time derivative of rhoU
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             drhoUdt(i,j,k)=sum(this%divu_x(:,i,j,k)*this%FUX(i-1:i,j,k))+&
   !             &              sum(this%divu_y(:,i,j,k)*this%FUY(i,j:j+1,k))+&
   !             &              sum(this%divu_z(:,i,j,k)*this%FUZ(i,j,k:k+1))+this%Pjx(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! Sync it
   !    call this%cfg%sync(drhoUdt)
      
   !    ! Flux of rhoV
   !    do kk=this%cfg%kmin_,this%cfg%kmax_+1
   !       do jj=this%cfg%jmin_,this%cfg%jmax_+1
   !          do ii=this%cfg%imin_,this%cfg%imax_+1
   !             ! Fluxes on x-face
   !             i=ii; j=jj; k=kk
   !             this%FVX(i,j,k)=this%FVX(i,j,k) &
   !             &              -sum(this%hybv_x(:,i,j,k)*this%rho_Vold(i-1:i,j,k))*sum(this%hybv_x(:,i,j,k)*this%V(i-1:i,j,k))*sum(this%itpu_y(:,i,j,k)*this%U(i,j-1:j,k)) &
   !             &              +this%visc_xy(i,j,k)*(sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))+sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k)))
   !             ! Fluxes on y-face
   !             i=ii-1; j=jj-1; k=kk-1
   !             this%FVY(i,j,k)=this%FVY(i,j,k)-this%P(i,j,k) &
   !             &              -sum(this%hybv_y(:,i,j,k)*this%rho_Vold(i,j:j+1,k))*sum(this%hybv_y(:,i,j,k)*this%V(i,j:j+1,k))*sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k)) &
   !             &              +this%visc   (i,j,k)*(sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))+sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k)))
   !             ! Fluxes on z-face
   !             i=ii; j=jj; k=kk
   !             this%FVZ(i,j,k)=this%FVZ(i,j,k) &
   !             &              -sum(this%hybv_z(:,i,j,k)*this%rho_Vold(i,j,k-1:k))*sum(this%hybv_z(:,i,j,k)*this%V(i,j,k-1:k))*sum(this%itpw_y(:,i,j,k)*this%W(i,j-1:j,k)) &
   !             &              +this%visc_yz(i,j,k)*(sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))+sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k)))
   !          end do
   !       end do
   !    end do
   !    ! Time derivative of rhoV
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             drhoVdt(i,j,k)=sum(this%divv_x(:,i,j,k)*this%FVX(i:i+1,j,k))+&
   !             &              sum(this%divv_y(:,i,j,k)*this%FVY(i,j-1:j,k))+&
   !             &              sum(this%divv_z(:,i,j,k)*this%FVZ(i,j,k:k+1))+this%Pjy(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! Sync it
   !    call this%cfg%sync(drhoVdt)
      
   !    ! Flux of rhoW
   !    do kk=this%cfg%kmin_,this%cfg%kmax_+1
   !       do jj=this%cfg%jmin_,this%cfg%jmax_+1
   !          do ii=this%cfg%imin_,this%cfg%imax_+1
   !             ! Fluxes on x-face
   !             i=ii; j=jj; k=kk
   !             this%FWX(i,j,k)=this%FWX(i,j,k) &
   !             &              -sum(this%hybw_x(:,i,j,k)*this%rho_Wold(i-1:i,j,k))*sum(this%hybw_x(:,i,j,k)*this%W(i-1:i,j,k))*sum(this%itpu_z(:,i,j,k)*this%U(i,j,k-1:k)) &
   !             &              +this%visc_zx(i,j,k)*(sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))+sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k)))
   !             ! Fluxes on y-face
   !             i=ii; j=jj; k=kk
   !             this%FWY(i,j,k)=this%FWY(i,j,k) &
   !             &              -sum(this%hybw_y(:,i,j,k)*this%rho_Wold(i,j-1:j,k))*sum(this%hybw_y(:,i,j,k)*this%W(i,j-1:j,k))*sum(this%itpv_z(:,i,j,k)*this%V(i,j,k-1:k)) &
   !             &              +this%visc_yz(i,j,k)*(sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))+sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k)))
   !             ! Fluxes on z-face
   !             i=ii-1; j=jj-1; k=kk-1
   !             this%FWZ(i,j,k)=this%FWZ(i,j,k)-this%P(i,j,k) &
   !             &              -sum(this%hybw_z(:,i,j,k)*this%rho_Wold(i,j,k:k+1))*sum(this%hybw_z(:,i,j,k)*this%W(i,j,k:k+1))*sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1)) &
   !             &              +this%visc   (i,j,k)*(sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))+sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1)))
   !          end do
   !       end do
   !    end do
   !    ! Time derivative of rhoW
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             drhoWdt(i,j,k)=sum(this%divw_x(:,i,j,k)*this%FWX(i:i+1,j,k))+&
   !             &              sum(this%divw_y(:,i,j,k)*this%FWY(i,j:j+1,k))+&
   !             &              sum(this%divw_z(:,i,j,k)*this%FWZ(i,j,k-1:k))+this%Pjz(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! Sync it
   !    call this%cfg%sync(drhoWdt)
   ! end subroutine get_dmomdt_full_ST

   ! subroutine get_dmomdt_no_ST(this,drhoUdt,drhoVdt,drhoWdt)
   !    implicit none
   !    class(tpns), intent(inout) :: this
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    integer :: i,j,k,ii,jj,kk

      
   !    ! Zero out drhoUVW/dt arrays
   !    drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
   !    ! Flux of rhoU
   !    do kk=this%cfg%kmin_,this%cfg%kmax_+1
   !       do jj=this%cfg%jmin_,this%cfg%jmax_+1
   !          do ii=this%cfg%imin_,this%cfg%imax_+1
   !             ! Fluxes on x-face
   !             i=ii-1; j=jj-1; k=kk-1
   !             this%FUX(i,j,k)=this%FUX(i,j,k)-this%P(i,j,k) &
   !             &              -sum(this%hybu_x(:,i,j,k)*this%rho_Uold(i:i+1,j,k))*sum(this%hybu_x(:,i,j,k)*this%U(i:i+1,j,k))*sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k)) &
   !             &              +this%visc   (i,j,k)*(sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))+sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k)))
   !             ! Fluxes on y-face
   !             i=ii; j=jj; k=kk
   !             this%FUY(i,j,k)=this%FUY(i,j,k) &
   !             &              -sum(this%hybu_y(:,i,j,k)*this%rho_Uold(i,j-1:j,k))*sum(this%hybu_y(:,i,j,k)*this%U(i,j-1:j,k))*sum(this%itpv_x(:,i,j,k)*this%V(i-1:i,j,k)) &
   !             &              +this%visc_xy(i,j,k)*(sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))+sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k)))
   !             ! Fluxes on z-face
   !             i=ii; j=jj; k=kk
   !             this%FUZ(i,j,k)=this%FUZ(i,j,k) &
   !             &              -sum(this%hybu_z(:,i,j,k)*this%rho_Uold(i,j,k-1:k))*sum(this%hybu_z(:,i,j,k)*this%U(i,j,k-1:k))*sum(this%itpw_x(:,i,j,k)*this%W(i-1:i,j,k)) &
   !             &              +this%visc_zx(i,j,k)*(sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))+sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k)))
   !          end do
   !       end do
   !    end do
   !    ! Time derivative of rhoU
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             drhoUdt(i,j,k)=sum(this%divu_x(:,i,j,k)*this%FUX(i-1:i,j,k))+&
   !             &              sum(this%divu_y(:,i,j,k)*this%FUY(i,j:j+1,k))+&
   !             &              sum(this%divu_z(:,i,j,k)*this%FUZ(i,j,k:k+1))+this%Pjx(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! Sync it
   !    call this%cfg%sync(drhoUdt)
      
   !    ! Flux of rhoV
   !    do kk=this%cfg%kmin_,this%cfg%kmax_+1
   !       do jj=this%cfg%jmin_,this%cfg%jmax_+1
   !          do ii=this%cfg%imin_,this%cfg%imax_+1
   !             ! Fluxes on x-face
   !             i=ii; j=jj; k=kk
   !             this%FVX(i,j,k)=this%FVX(i,j,k) &
   !             &              -sum(this%hybv_x(:,i,j,k)*this%rho_Vold(i-1:i,j,k))*sum(this%hybv_x(:,i,j,k)*this%V(i-1:i,j,k))*sum(this%itpu_y(:,i,j,k)*this%U(i,j-1:j,k)) &
   !             &              +this%visc_xy(i,j,k)*(sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))+sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k)))
   !             ! Fluxes on y-face
   !             i=ii-1; j=jj-1; k=kk-1
   !             this%FVY(i,j,k)=this%FVY(i,j,k)-this%P(i,j,k) &
   !             &              -sum(this%hybv_y(:,i,j,k)*this%rho_Vold(i,j:j+1,k))*sum(this%hybv_y(:,i,j,k)*this%V(i,j:j+1,k))*sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k)) &
   !             &              +this%visc   (i,j,k)*(sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))+sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k)))
   !             ! Fluxes on z-face
   !             i=ii; j=jj; k=kk
   !             this%FVZ(i,j,k)=this%FVZ(i,j,k) &
   !             &              -sum(this%hybv_z(:,i,j,k)*this%rho_Vold(i,j,k-1:k))*sum(this%hybv_z(:,i,j,k)*this%V(i,j,k-1:k))*sum(this%itpw_y(:,i,j,k)*this%W(i,j-1:j,k)) &
   !             &              +this%visc_yz(i,j,k)*(sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))+sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k)))
   !          end do
   !       end do
   !    end do
   !    ! Time derivative of rhoV
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             drhoVdt(i,j,k)=sum(this%divv_x(:,i,j,k)*this%FVX(i:i+1,j,k))+&
   !             &              sum(this%divv_y(:,i,j,k)*this%FVY(i,j-1:j,k))+&
   !             &              sum(this%divv_z(:,i,j,k)*this%FVZ(i,j,k:k+1))+this%Pjy(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! Sync it
   !    call this%cfg%sync(drhoVdt)
      
   !    ! Flux of rhoW
   !    do kk=this%cfg%kmin_,this%cfg%kmax_+1
   !       do jj=this%cfg%jmin_,this%cfg%jmax_+1
   !          do ii=this%cfg%imin_,this%cfg%imax_+1
   !             ! Fluxes on x-face
   !             i=ii; j=jj; k=kk
   !             this%FWX(i,j,k)=this%FWX(i,j,k) &
   !             &              -sum(this%hybw_x(:,i,j,k)*this%rho_Wold(i-1:i,j,k))*sum(this%hybw_x(:,i,j,k)*this%W(i-1:i,j,k))*sum(this%itpu_z(:,i,j,k)*this%U(i,j,k-1:k)) &
   !             &              +this%visc_zx(i,j,k)*(sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))+sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k)))
   !             ! Fluxes on y-face
   !             i=ii; j=jj; k=kk
   !             this%FWY(i,j,k)=this%FWY(i,j,k) &
   !             &              -sum(this%hybw_y(:,i,j,k)*this%rho_Wold(i,j-1:j,k))*sum(this%hybw_y(:,i,j,k)*this%W(i,j-1:j,k))*sum(this%itpv_z(:,i,j,k)*this%V(i,j,k-1:k)) &
   !             &              +this%visc_yz(i,j,k)*(sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))+sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k)))
   !             ! Fluxes on z-face
   !             i=ii-1; j=jj-1; k=kk-1
   !             this%FWZ(i,j,k)=this%FWZ(i,j,k)-this%P(i,j,k) &
   !             &              -sum(this%hybw_z(:,i,j,k)*this%rho_Wold(i,j,k:k+1))*sum(this%hybw_z(:,i,j,k)*this%W(i,j,k:k+1))*sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1)) &
   !             &              +this%visc   (i,j,k)*(sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))+sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1)))
   !          end do
   !       end do
   !    end do
   !    ! Time derivative of rhoW
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             drhoWdt(i,j,k)=sum(this%divw_x(:,i,j,k)*this%FWX(i:i+1,j,k))+&
   !             &              sum(this%divw_y(:,i,j,k)*this%FWY(i,j:j+1,k))+&
   !             &              sum(this%divw_z(:,i,j,k)*this%FWZ(i,j,k-1:k)) ! +this%Pjz(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! Sync it
   !    call this%cfg%sync(drhoWdt)
      
   ! end subroutine get_dmomdt_no_ST
   
   ! ################################################################## SURFACE TENSION REPLACEMENT SUBROUTINES ##################################################################
   ! Full ST in Poisson
   ! subroutine add_surface_tension_jump_full_ST(this,dt,div,vf,contact_model)
   !    use messager,  only: die
   !    use vfs_class, only: vfs
   !    implicit none
   !    class(tpns), intent(inout) :: this
   !    real(WP), intent(inout) :: dt     !< Timestep size over which to advance
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    class(vfs), intent(inout) :: vf
   !    integer, intent(in), optional :: contact_model
   !    integer :: i,j,k
   !    real(WP) :: mycurv,mysurf
      
   !    ! Store old jump
   !    this%DPjx=this%Pjx
   !    this%DPjy=this%Pjy
   !    this%DPjz=this%Pjz
      
   !    ! Calculate pressure jump
   !    do k=this%cfg%kmin_,this%cfg%kmax_+1
   !       do j=this%cfg%jmin_,this%cfg%jmax_+1
   !          do i=this%cfg%imin_,this%cfg%imax_+1
   !             ! X face
   !             mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k))
   !             ! Y face
   !             mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*vf%VF(i,j-1:j,k))
   !             ! Z face
   !             mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k))
   !          end do
   !       end do
   !    end do
      
   !    ! Add wall contact force to pressure jump
   !    if (present(contact_model)) then
   !       select case (contact_model)
   !       case (1)
   !          call this%add_static_contact(vf=vf)
   !       case default
   !          call die('[tpns: add_surface_tension_jump] Unknown contact model!')
   !       end select
   !    end if
      
   !    ! Compute jump of DP
   !    this%DPjx=this%Pjx-this%DPjx
   !    this%DPjy=this%Pjy-this%DPjy
   !    this%DPjz=this%Pjz-this%DPjz
      
   !    ! Add div(Pjump) to RP
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             div(i,j,k)=div(i,j,k)+dt*(sum(this%divp_x(:,i,j,k)*this%Pjx(i:i+1,j,k)/this%rho_U(i:i+1,j,k))&
   !             &                        +sum(this%divp_y(:,i,j,k)*this%Pjy(i,j:j+1,k)/this%rho_V(i,j:j+1,k))&
   !             &                        +sum(this%divp_z(:,i,j,k)*this%Pjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
   !          end do
   !       end do
   !    end do
   ! end subroutine add_surface_tension_jump_full_ST

   ! ! No ST in Poisson, just recalc Pj
   ! subroutine add_surface_tension_jump_no_ST(this,dt,div,vf,contact_model)
   !    use messager,  only: die
   !    use vfs_class, only: vfs
   !    implicit none
      
   !    class(tpns), intent(inout) :: this
   !    real(WP), intent(inout) :: dt     !< Timestep size over which to advance
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_,this%cfg%jmino_,this%cfg%kmino_) :: addToDiv
   !    class(vfs), intent(inout) :: vf
   !    integer, intent(in), optional :: contact_model
   !    integer :: i,j,k
   !    real(WP) :: mycurv,mysurf, magx,magy,magz,magdiv
   !    ! Store old jump
   !    this%DPjx=this%Pjx
   !    this%DPjy=this%Pjy
   !    this%DPjz=this%Pjz
      
   !    ! Calculate pressure jump
   !    do k=this%cfg%kmin_,this%cfg%kmax_+1
   !       do j=this%cfg%jmin_,this%cfg%jmax_+1
   !          do i=this%cfg%imin_,this%cfg%imax_+1
   !             ! X face
   !             mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k))
   !             ! Y face
   !             mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*vf%VF(i,j-1:j,k))
               ! Z face
   !             mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k))
   !          end do
   !       end do
   !    end do
      
   !    ! Add wall contact force to pressure jump
   !    if (present(contact_model)) then
   !       select case (contact_model)
   !       case (1)
   !          call this%add_static_contact(vf=vf)
   !       case default
   !          call die('[tpns: add_surface_tension_jump] Unknown contact model!')
   !       end select
   !    end if
      
   !    ! Compute jump of DP
   !    this%DPjx=this%Pjx-this%DPjx
   !    this%DPjy=this%Pjy-this%DPjy
   !    this%DPjz=this%Pjz-this%DPjz
      
   !    ! Calculate Magnitude of DPs
   !    magx = NORM2(this%DPjx)
   !    magy = NORM2(this%DPjy)
   !    magz = NORM2(this%DPjz)
   !    magdiv = NORM2(this%div)
   !    addToDiv = dt*(sum(this%divp_x(:,i,j,k)*this%DPjx(i:i+1,j,k)/this%rho_U(i:i+1,j,k))&
   !             &                        +sum(this%divp_y(:,i,j,k)*this%DPjy(i,j:j+1,k)/this%rho_V(i,j:j+1,k))&
   !             &                        +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
      
   !    magz = NORM2(addToDiv)
   !    write(*,'(A,4F35.10)') '   magx,magy,magz,magdiv: ', magx,magy,magz,magdiv
   !    if(magz.gt.0) then
   !         write(*,'(A,F10.5)')  '====================================MagToDiv: ',magz
   !    endif
   ! end subroutine add_surface_tension_jump_no_ST
   
   subroutine computeInterfaceHeight(vf,i0,j0,k0,dir,ret)
      use messager,  only: die
      use vfs_class, only: vfs
      implicit none
      class(vfs), intent(inout) :: vf
      integer :: i0,j0,k0 !i0,j0,k0 tell us the original cell location
      integer :: dir  ! 1 = x oriented, 2 = y oriented, 3 = z oriented. Sign gives direction
      integer :: iplus,jplus,kplus ! Tells us how to add
      integer :: i,j,k,count ! Current Cell Location
      real(WP),dimension(4) :: ret
      real(WP) :: H,ct,cb
      LOGICAL :: bound,whileTest

      ! Initial Values
      ! write(*,'(A)') 'Start Function'
      i = i0
      j = j0
      k = k0
      ! write(*,'(A)') 'Get Volume Fractions'
      ct = vf%VF(i,j,k)
      H = vf%VF(i0,j0,k0)
      ! write(*,'(A)') 'End Volume Fractions'
      ! write(*,'(A,3I6.2,F6.2,I6.2)') 'Initial i,j,k,ct,dir = ',i,j,k,ct,dir
      count = 0
      ! Convert dir to index additions
      SELECT CASE (abs(dir))
         CASE (1) ! x oriented
            iplus = sign(1,dir)
            jplus = 0
            kplus = 0
         CASE (2) ! y oriented
            iplus = 0
            jplus = sign(1,dir)
            kplus = 0
         CASE (3) ! z oriented
            iplus = 0
            jplus = 0
            kplus = sign(1,dir)
      END SELECT
      ! write(*,'(A,3I6.2)') 'iplus,jplus,kplus = ', iplus,jplus,kplus
      ! Set initial Bound
      if(ct .lt. VFhi) then 
         bound = .TRUE.
      else
         bound = .FALSE.
      endif

      ! Find Top Out
      ! write(*,'(A)') 'Start Top Loop'
      ! write(*,'(A,F6.2)') 'Pre Loop H = ',H
      do while((.not. bound) .or. (ct .gt. VFlo) .and. count .lt. 3) ! If not at interface, or at interface exactly, then go to top neighbor
         ! write(*,'(A)')  '========= while 1'
         ! Go To Top N eighbor
         i = i + iplus
         j = j + jplus
         k = k + kplus
         ! Update ct
         ct = vf%VF(i,j,k)
         ! write(*,'(A,3I6.2,F6.2)') 'i,j,k,ct = ',i,j,k,ct
         ! Update H
         H = H + ct
         ! write(*,'(A,F6.2)') 'H = ',H
         if((ct .gt. VFlo) .and. (ct .le. VFhi)) then
            bound = .true.
         endif
         count = count + 1
      enddo
      ! write(*,'(A,F6.2)') 'End Top Loop COu = ',H
      ! Catch Inconsistent result
      if((ct .gt. VFlo)) then
         H = -1.0_WP
         ret = [-1.0_WP,-1.0_WP,-1.0_WP,-1.0_WP]
         RETURN
      endif

      ! Reset to center
      i = i0
      j = j0
      k = k0
      cb = vf%VF(i,j,k)
      ! write(*,'(A,3I6.2,F6.2)') 'Initial i,j,k,cb = ',i,j,k,cb
      count = 0
      ! Set initial Bound
      if(cb .gt. VFlo) then 
         bound = .TRUE.
      else
         bound = .FALSE.
      endif
      whileTest = (.not. bound) .or. (((cb .lt. VFhi) .and. (cb .gt. VFlo)) .and. count .lt. 3)
      
      ! Find Bottom in
      ! write(*,'(A)') 'Start Bot Loop'
      do while((.not. bound) .or. (((cb .lt. VFhi) .and. (cb .gt. VFlo)) .and. count .lt. 3)) ! If not at interface, or at interface exactly, then go to Bottom neighbor
         ! write(*,'(A)')  '========= while 2' 
         
         ! Go To Bottom Neighbor
         i = i - iplus
         j = j - jplus
         k = k - kplus
         ! write(*,'(A,3I10.5)')  '========= i,j,k = ',i,j,k 
         ! Update cb
         cb = vf%VF(i,j,k)
         ! write(*,'(A,3I6.2,F6.2)') 'i,j,k,cb = ',i,j,k,cb
         ! Update H
         H = H + cb
         if((cb .gt. VFlo) .and. (cb .le. VFhi)) then
            bound = .true.
         endif
         count = count + 1
      enddo
      ! write(*,'(A,F6.2)') 'End Bot Loop H = ',H
      ! Catch Inconsistent result
      if(cb .lt. VFhi) then
         ! print *, "Inconsistent 2"
         H = -2.0_WP
         ret = [-2.0_WP,-2.0_WP,-2.0_WP,-2.0_WP]
         RETURN
      endif
      !return N and H
      ret = [0.0_WP,0.0_WP,0.0_WP,0.0_WP]
      ret(1) = i
      ret(2) = j
      ret(3) = k
      ret(4) = H
   end subroutine computeInterfaceHeight

   subroutine heightFunctionCurvature(i,j,k,vf,mycurv,dir,hs)
      use messager,  only: die
      use vfs_class, only: vfs
      implicit none
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,n,iplus,jplus,kplus ! Current Cell Location
      integer :: dir  ! 1 = x oriented, 2 = y oriented, 3 = z oriented. Sign gives direction
      real(WP), dimension(1) :: mycurv,mysurf
      real(WP), dimension(1,3) :: mynorm
      real(WP) :: h0,hm1,hp1,dx,hP,hPP,shift ! h0 is the center height, hm1 is the cell to the "left" and hp1 is the cell to the right.
      real(WP), dimension(4) :: cellAndHeight
      integer, dimension(3) :: cell0,cellm1,cellp1
      real(WP),dimension(5),optional :: hs

      ! Convert dir to index additions in perpendicular direction
      ! write (*,'(A,I6.2)') '============================ START HEIGHT FUNCTION CURVATURE= ',dir
      ! write (*,'(A,I6.2)') 'dir = ',dir
      SELECT CASE (abs(dir))
         CASE (1) ! x oriented
            iplus = 0
            jplus = -(sign(int(1),dir))
            kplus = -(sign(int(1),dir))
         CASE (2) ! y oriented
            iplus = (sign(int(1),dir))
            jplus = 0
            kplus = (sign(int(1),dir))
         CASE (3) ! z oriented
            iplus = (sign(int(1),dir))
            jplus = (sign(int(1),dir))
            kplus = 0
      END SELECT

      ! Calculate H0,N0
      ! write (*,'(A)') 'COMPUTE H0 ================================================='
      call computeInterfaceHeight(vf,i,j,k,dir,cellAndHeight)
      cell0(1) = cellAndHeight(1)
      cell0(2) = cellAndHeight(2)
      cell0(3) = cellAndHeight(3)
      h0 = cellAndHeight(4)
      
      ! Right now, focus on 2D, so there are only 2
      ! Calculate hp1, Np1
      ! write (*,'(A)') 'COMPUTE Hp1 ================================================='
      ! write (*,'(A,3I6.2)') 'PLUS i,j,dir = ',iplus,jplus,dir
      call computeInterfaceHeight(vf,i+iplus,j+jplus,k,dir,cellAndHeight)
      cellp1(1) = cellAndHeight(1)
      cellp1(2) = cellAndHeight(2)
      cellp1(3) = cellAndHeight(3)
      hp1 = cellAndHeight(4)
      
      ! Calculate hm1, Nm1
      ! write (*,'(A)') 'COMPUTE Hm1 ================================================='
      call computeInterfaceHeight(vf,i-iplus,j-jplus,k,dir,cellAndHeight)
      cellm1(1) = cellAndHeight(1)
      cellm1(2) = cellAndHeight(2)
      cellm1(3) = cellAndHeight(3)
      hm1 = cellAndHeight(4)
      ! write (*,'(A)') 'COMPUTE COMPLETE ================================================='
      ! Check for Consistency
      hs = [hm1,h0,hp1,hP,hPP]
      if((sign(1.0_WP,hp1) .ne. -1.0_WP) .and. (sign(1.0_WP,hm1) .ne. -1.0_WP)) then 
         ! consistent
         ! write (*,'(A)') 'dx calc'
         dx = (vf%cfg%xm(cell0(1)) - vf%cfg%xm(cell0(1)+iplus)) + (vf%cfg%ym(cell0(2)) - vf%cfg%ym(cell0(2)+jplus))
         dx = -dx
         ! write (*,'(A,F10.2)') 'dx = ',dx
         ! write (*,'(A,F10.2)') 'dx = ',dx
         ! write (*,'(A,F10.2)') 'dx = ',dx
         ! write (*,'(A)') 'START CONSISTENCY CHECK'
         h0 = h0*abs(dx)
         ! Move Heights to a common origin
         ! write(*,'(A)')''
         ! write(*,'(A)')''
         ! write(*,'(A)')''
         ! write(*,'(A,F10.5)') 'hp=m1 Val Before mult',hm1
         ! write(*,'(A,F10.5)') 'dx',dx
         hm1 = hm1*abs(dx)
         shift = vf%cfg%xm(cellm1(1)) * (2-abs(dir))+ vf%cfg%ym(cellm1(2)) * (abs(dir)-1)- vf%cfg%xm(cell0(1)) * (2-abs(dir)) - vf%cfg%ym(cell0(2)) * (abs(dir)-1)
         ! write(*,'(A,F10.5)') 'hm1 Val Before Shift',hm1
         ! write(*,'(A,3I10.1)') 'cellm1: ',cellm1(1),cellm1(2),cellm1(3)
         ! write(*,'(A,3I10.1)') 'cell0: ',cell0
         ! write(*,'(A,I10.1)') 'dir: ',dir
         hm1 = hm1 + shift*sign(1.0_WP,real(dir,WP))
         ! write(*,'(A,F10.5)') 'hm1 Val After Shift2',hp1
         ! write(*,'(A,F10.5)') 'term1: ',vf%cfg%xm(cellm1(1)) * (2-abs(dir))
         ! write(*,'(A,F10.5)') 'term2: ',vf%cfg%ym(cellm1(2)) * (abs(dir)-1)
         ! write(*,'(A,F10.5)') 'term3: ',vf%cfg%xm(cell0(1)) * (2-abs(dir))
         ! write(*,'(A,F10.5)') 'term4: ',vf%cfg%ym(cell0(2)) * (abs(dir)-1)
         
         ! write(*,'(A,F10.5)') 'Net: ',shift
         ! hm1 = hm1 + vf%cfg%xm(cellm1(1)) * (2-abs(dir)) + vf%cfg%ym(cellm1(2)) * (abs(dir)-1)
         ! hm1 = hm1 - vf%cfg%xm(cell0(1)) * (2-abs(dir)) - vf%cfg%ym(cell0(2)) * (abs(dir)-1)
         ! write (*,'(A)') 'hm shift'
         ! write(*,'(A)')''
         ! write(*,'(A)')''
         ! write(*,'(A)')''
         ! write(*,'(A,F10.5)') 'hp1 Val Before mult',hp1
         ! write(*,'(A,F10.5)') 'dx',dx
         hp1 = hp1*abs(dx)
         ! write(*,'(A,F10.5)') 'hp1 Val Before Shift',hp1
         ! write(*,'(A,3I10.1)') 'cellp1: ',cellp1(1),cellp1(2),cellp1(3)
         ! write(*,'(A,3I10.1)') 'cell0: ',cell0
         ! write(*,'(A,I10.1)') 'dir: ',dir
         shift = vf%cfg%xm(cellp1(1)) * (2-abs(dir))+ vf%cfg%ym(cellp1(2)) * (abs(dir)-1)- vf%cfg%xm(cell0(1)) * (2-abs(dir)) - vf%cfg%ym(cell0(2)) * (abs(dir)-1)
         hp1 = hp1 + shift*sign(1.0_WP,real(dir,WP)) 
         ! hp1 = hp1 + vf%cfg%xm(cellp1(1)) * (2-abs(dir))+ vf%cfg%ym(cellp1(2)) * (abs(dir)-1)
         ! write(*,'(A,F20.8)') 'hp1 Val After Shift1',hp1
         ! hp1 = hp1 - vf%cfg%xm(cell0(1)) * (2-abs(dir)) - vf%cfg%ym(cell0(2)) * (abs(dir)-1)
         ! write(*,'(A,F10.5)') 'hp1 Val After Shift2',hp1
         ! write(*,'(A,F10.5)') 'term1: ',vf%cfg%xm(cellp1(1)) * (2-abs(dir))
         ! write(*,'(A,F10.5)') 'term2: ',vf%cfg%ym(cellp1(2)) * (abs(dir)-1)
         ! write(*,'(A,F10.5)') 'term3: ',vf%cfg%xm(cell0(1)) * (2-abs(dir))
         ! write(*,'(A,F10.5)') 'term4: ',vf%cfg%ym(cell0(2)) * (abs(dir)-1)
         
         ! write(*,'(A,F10.5)') 'Net: ',shift
         ! write (*,'(A)') 'hp shift'
         ! Now that everything is in the commmon frame, finite differences
         
         ! write(*,'(A)')''
         ! write(*,'(A)')''
         ! write(*,'(A)')''
         ! write (*,'(A,F10.5)') 'h0 = ',h0
         ! write (*,'(A,F10.5)') 'hp1 = ',hp1
         ! write (*,'(A,F10.5)') 'hm1 = ',hm1

         hP = (hp1-hm1)/(2*abs(dx))
         ! write (*,'(A,F10.5)') 'hP calc',hP
         
         hPP = (hp1 - 2* h0 + hm1)/(dx*dx)
         ! write (*,'(A,F10.5)') 'hPP calc',hPP

         hs = [hm1/dx,h0/dx,hp1/dx,hP,hPP]
         ! Calculate Curvature
         mycurv = -hPP/((1.0_WP+hP*hP)**(1.5_WP))
         ! I have included a factor of -1 here. This is because all the curvatures are exactly negative what they should be.
         ! The current theory is that I should not multiply by abs(dx) for h0,hm1,hp1 but simply dx. This could fix that problem.
         ! However, doing this would require me to then reevaluate how the shift is done, and i think that would hurt me right now.
         ! As such, I will simply not be doing that right, and instead adding in the -1. I think in an ideal world I would fix this, 
         ! but we do not live in an ideal world. 
         ! write(*,'(A,F6.2)') 'Final Curvature: ',mycurv
      else 
         mycurv = vf%curv(i,j,k)
      endif
      
   end subroutine heightFunctionCurvature

   ! ============================== Conservative Surface Tension Subroutines
   ! subroutine updateSurfaceTensionStresses(fs,vf)
   !    use irl_fortran_interface
   !    use f_PUNeigh_RectCub_class
   !    use f_SeparatorVariant_class
   !    use f_PUSolve_RectCub_class

   !    implicit none
   !    class(vfs), intent(inout) :: vf ! Volume Fraction Solver
   !    class(tpns), intent(inout) :: fs ! Two Phase Flow Solver
   !    integer :: i,j,k,j_in,i_in ! Current Cell Location
   !    type(PUNeigh_RectCub_type) :: neighborhood
   !    type(PU_RectCub_type) :: solver
      
   !    ! Temp Items
   !    type(SeparatorVariant_type) :: plane
   !    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
   !    integer, dimension(1:3) :: ind,indguess
   !    real(WP) :: dx,dy,C,xEval,yEval,zEval

   !    dx = vf%cfg%dx(1)
   !    dy = vf%cfg%dy(1)
   !    C = 4.0_WP
   !    ! Create Neighborhood and solver
   !    call new(neighborhood) 
   !    call new(solver)
   !    ! Loop over real domain 
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then

   !                ! Empty Neighborhood
   !                call emptyNeighborhood(neighborhood) 
   !                ! Add 7x7 (3 on each side) stencil, in plane
   !                ! write(*,'(A)') '=============== New: '
   !                do j_in = -3,3
   !                   do i_in = -3,3
   !                      if(vf%VF(i+i_in,j+j_in,k) .gt. vf%VFmin .and. vf%VF(i+i_in,j+j_in,k) .lt. vf%VFmax) then
   !                         ! For now, compute centroid as cell center
   !                         ! Gcen = vf%Gbary(:,i+i_in,j+j_in,k)
   !                         ! Lcen = vf%Lbary(:,i+i_in,j+j_in,k)
   !                         ! ! cen = vf%VF(i+i_in,j+j_in,k) * Lcen + (1-vf%VF(i+i_in,j+j_in,k)) * Gcen
   !                         ! ! cen = (Lcen+Gcen)/2
   !                         ! cen = (/vf%cfg%xm(i+i_in),vf%cfg%ym(j+j_in),vf%cfg%zm(k)/)
   !                         cen = calculateCentroid(vf%interface_polygon(1,i+i_in,j+j_in,k))
   !                         ! Get Plane
   !                         plane = vf%liquid_gas_interface(i+i_in,j+j_in,k)
   !                         ! if((vf%VF(i+i_in,j+j_in,k) .gt. VFlo .and. vf%VF(i+i_in,j+j_in,k) .lt. VFhi) .and. (i .eq. 16) .and. (j .eq. 24)) then
   !                         !    write(*,'(A,3I10.5)') 'Mid Cell: ', i,j,k
   !                         !    write(*,'(A,3I10.5)') 'Curr Cell: ', i+i_in,j+j_in,k 
   !                         !    write(*,'(A,3F10.5)') 'Cen: ', cen 
   !                         !    call printToScreen(plane)
   !                         call addMember(neighborhood,cen,1.0_WP,plane)
   !                      endif
   !                   end do
   !                end do
   !                ! Set Neighborhood in solver
   !                call setNeighborhood(solver,neighborhood)
   !                call setThreshold(solver,0.1875_WP) ! This is the of the wendland function at 0.5 (is radius is 1).  
   !                ! ====== Get Stresses

   !                ! sigma_xx, we go on right of u cell
   !                ! The right of the u cell goes y(j) to y(j+1) at xm(i)
   !                startPoint = (/cfg%xm(i),cfg%y(j),cfg%zm(k)/)
   !                endPoint = (/cfg%xm(i),cfg%y(j+1),cfg%zm(k)/)
   !                force = (/0.0_WP,0.0_WP,0.0_WP/)
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,PressureOption,MarangoniOption,force)
                  
   !                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
   !                ! We only want the x component of this force, so take it and store it as sigma_xx
   !                sigma_xx(i,j,k) = force(1)

   !                ! ====== For Debugging
   !                ! No Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,0.0_WP,MarangoniOption,force)
   !                sigma_xx_NoP(i,j,k) = force(1)
   !                ! With Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,1.0_WP,MarangoniOption,force)
   !                sigma_xx_P(i,j,k) = force(1)
   !                sigma_xx_P(i,j,k) = sigma_xx_P(i,j,k) - sigma_xx_NoP(i,j,k)


   !                ! sigma_xy, we go on bottom of u cell
   !                ! The bottom of the u cell goes xm(i-1) to xm(i) at y(j)
   !                startPoint = (/cfg%xm(i-1),cfg%y(j),cfg%zm(k)/)
   !                endPoint = (/cfg%xm(i),cfg%y(j),cfg%zm(k)/)
   !                force = (/0.0_WP,0.0_WP,0.0_WP/)
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,PressureOption,MarangoniOption,force)
   !                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
   !                ! We only want the x component of this force, so take it and store it as sigma_xy 
   !                sigma_xy(i,j,k) = force(1)

   !                ! ====== For Debugging
   !                ! No Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,0.0_WP,MarangoniOption,force)
   !                sigma_xy_NoP(i,j,k) = force(1)
   !                ! With Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,1.0_WP,MarangoniOption,force)
   !                sigma_xy_P(i,j,k) = force(1)
   !                sigma_xy_P(i,j,k) = sigma_xy_P(i,j,k) - sigma_xy_NoP(i,j,k)

   !                ! sigma_yx, we go on left of v cell
   !                ! The bottom of the u cell goes ym(j-1) to ym(j) at x(i)
   !                startPoint = (/cfg%x(i),cfg%ym(j),cfg%zm(k)/)
   !                endPoint = (/cfg%x(i),cfg%ym(j-1),cfg%zm(k)/)
   !                force = (/0.0_WP,0.0_WP,0.0_WP/)
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,PressureOption,MarangoniOption,force)
   !                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
   !                ! We only want the y component of this force, so take it and store it as sigma_yx 
   !                sigma_yx(i,j,k) = force(2)
                  
   !                ! ====== For Debugging
   !                ! No Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,0.0_WP,MarangoniOption,force)
   !                sigma_yx_NoP(i,j,k) = force(2)
   !                ! With Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,1.0_WP,MarangoniOption,force)
   !                sigma_yx_P(i,j,k) = force(2)
   !                sigma_yx_P(i,j,k) = sigma_yx_P(i,j,k) - sigma_yx_NoP(i,j,k)

   !                ! sigma_yy, we go on top of v cell
   !                ! The bottom of the u cell goes x(i) to x(i+1) at ym(j) 
   !                startPoint = (/cfg%x(i+1),cfg%ym(j),cfg%zm(k)/)
   !                endPoint = (/cfg%x(i),cfg%ym(j),cfg%zm(k)/)
   !                force = (/0.0_WP,0.0_WP,0.0_WP/)
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,PressureOption,MarangoniOption,force)
   !                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
   !                ! We only want the y component of this force, so take it and store it as sigma_yy
   !                sigma_yy(i,j,k) = force(2)

   !                ! ====== For Debugging  
   !                ! No Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,0.0_WP,MarangoniOption,force)
   !                sigma_yy_NoP(i,j,k) = force(2)
   !                ! With Pressure
   !                call solveEdge(solver,fs%sigma,startPoint,endPoint,PU_spread*dx,1.0_WP,MarangoniOption,force)
   !                sigma_yy_P(i,j,k) = force(2)
   !                sigma_yy_P(i,j,k) = sigma_yy_P(i,j,k) - sigma_yy_NoP(i,j,k)
   !             ! endif
   !          end do
   !       end do
   !    end do
   !    ! Boundary Conditions
   !    call vf%cfg%sync(sigma_xx)
   !    call vf%cfg%sync(sigma_xy)
   !    call vf%cfg%sync(sigma_yx)
   !    call vf%cfg%sync(sigma_yy)

   !    call vf%cfg%sync(sigma_xx_NoP)
   !    call vf%cfg%sync(sigma_xy_NoP)
   !    call vf%cfg%sync(sigma_yx_NoP)
   !    call vf%cfg%sync(sigma_yy_NoP)

   !    call vf%cfg%sync(sigma_xx_P)
   !    call vf%cfg%sync(sigma_xy_P)
   !    call vf%cfg%sync(sigma_yx_P)
   !    call vf%cfg%sync(sigma_yy_P)
   ! end subroutine updateSurfaceTensionStresses

   ! subroutine updateSurfaceTensionForces(fs)
   !    use irl_fortran_interface
   !    use f_PUNeigh_RectCub_class
   !    use f_SeparatorVariant_class
   !    use f_PUSolve_RectCub_class

   !    class(tpns), intent(inout) :: fs ! Two Phase Flow Solver
   !    integer :: i,j,k 
      
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             ! write(*,'(A,F10.5)') 'DX,i,j,k: ', cfg%dxmi(i)
   !             Fst_x(i,j,k) = (sigma_xx(i,j,k)-sigma_xx(i-1,j,k)) - (sigma_xy(i,j+1,k)-sigma_xy(i,j,k))
   !             Fst_x(i,j,k) = Fst_x(i,j,k) * cfg%dxi(i)

   !             Fst_y(i,j,k) = (sigma_yy(i,j,k)-sigma_yy(i,j-1,k)) - (sigma_yx(i+1,j,k)-sigma_yx(i,j,k))
   !             Fst_y(i,j,k) = Fst_y(i,j,k) * cfg%dyi(j)

   !             Fst_z(i,j,k) = 0.0_WP
   !          end do
   !       end do
   !    end do
   !    fs%Pjx = Fst_x 
   !    fs%Pjy = Fst_y
   !    fs%Pjz = Fst_z

   !    Pjx_ST = Fst_x
   !    Pjy_ST = Fst_y
   !    Pjz_ST = Fst_z
   ! end subroutine updateSurfaceTensionForces


   ! subroutine applyLaplacianSmoothing(fs)
   !    use irl_fortran_interface
   !    use f_PUNeigh_RectCub_class
   !    use f_SeparatorVariant_class
   !    use f_PUSolve_RectCub_class

   !    class(tpns), intent(inout) :: fs ! Two Phase Flow Solver
   !    integer :: i,j,k 
      
   !    Fst_x = 0.0_WP;
   !    Fst_y = 0.0_WP;
   !    Fst_z = 0.0_WP;

   !    ! Smoothing
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             Fst_x(i,j,k) = 0.125*(fs%Pjx(i+1,j,k)+fs%Pjx(i-1,j,k) + fs%Pjx(i,j+1,k) + fs%Pjx(i,j-1,k)) + 0.5*fs%Pjx(i,j,k)
   !             Fst_y(i,j,k) = 0.125*(fs%Pjy(i+1,j,k)+fs%Pjy(i-1,j,k) + fs%Pjy(i,j+1,k) + fs%Pjy(i,j-1,k)) + 0.5*fs%Pjy(i,j,k)
   !          end do
   !       end do
   !    end do

   !    fs%Pjx = Fst_x 
   !    fs%Pjy = Fst_y
   !    fs%Pjz = Fst_z
      
      
   ! end subroutine applyLaplacianSmoothing

   ! subroutine applyGradientSmoothing(fs,vf)
   !    use irl_fortran_interface
   !    use f_PUNeigh_RectCub_class
   !    use f_SeparatorVariant_class
   !    use f_PUSolve_RectCub_class

   !    class(tpns), intent(inout) :: fs ! Two Phase Flow Solver
   !    class(vfs), intent(inout) :: vf ! Volume Fraction Solver
   !    integer :: i,j,k,ii,jj,kk,negativeIndex,positiveIndex
   !    real(WP) :: tot_grad_x,tot_grad_y,tot_force

   !    Fst_x = 0.0_WP;
   !    Fst_y = 0.0_WP;
   !    Fst_z = 0.0_WP;
      
   !    negativeIndex = -1
   !    positiveIndex = 1
   !    ! Calculate Gradients
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             grad_vf_x(i,j,k) = abs((vf%VF(i  ,j  ,k) - vf%VF(i-1,j  ,k))  * vf%cfg%dxi(i))
   !             grad_vf_y(i,j,k) = abs((vf%VF(i  ,j  ,k) - vf%VF(i  ,j-1,k)) * vf%cfg%dyi(j))
   !             ! While we are looping everything, also check to see if everything has been processed
   !             if(abs(fs%Pjx(i,j,k)) .gt. 1e-10) then
   !                x_smoothing_tracker(i,j,k) = 1 ! Mark as needing smoothing
   !             endif

   !             if(abs(fs%Pjy(i,j,k)) .gt. 1e-10) then
   !                y_smoothing_tracker(i,j,k) = 1 ! Mark as needing smoothing
   !             endif

   !          end do
   !       end do
   !    end do
   !    ! Smoothing
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             tot_grad_x = 0.0_WP
   !             tot_grad_y = 0.0_WP
   !             ! X Direction
   !             if(abs(fs%Pjx(i,j,k)) .gt. 1e-12 .and. x_smoothing_tracker(i,j,k) == 1) then  ! This checks to see if we are in a mixed cell that has not yet been smoothed
   !                ! Current Cell
   !                tot_grad_x = grad_vf_x(i,j,k)
   !                tot_force = fs%Pjx(i,j,k)
   !                ! Going Left
   !                do while(abs(fs%Pjx(i+negativeIndex,j,k)) .gt. 1e-12 .and. (i+negativeIndex) .ge. cfg%imin_)
   !                   tot_grad_x = tot_grad_x + grad_vf_x(i+negativeIndex,j,k)
   !                   tot_force = tot_force + fs%Pjx(i+negativeIndex,j,k)
   !                   negativeIndex = negativeIndex -1
   !                end do
   !                ! Going Right
   !                do while(abs(fs%Pjx(i+positiveIndex,j,k)) .gt. 1e-12 .and. (i+positiveIndex) .le. cfg%imax_)
   !                   tot_grad_x = tot_grad_x + grad_vf_x(i+positiveIndex,j,k)
   !                   tot_force = tot_force + fs%Pjx(i+positiveIndex,j,k)
   !                   positiveIndex = positiveIndex +1
   !                end do
   !                ! Now Reassign Force Values in between 
   !                do ii = negativeIndex,positiveIndex
   !                   Fst_x(i+ii,j,k) = grad_vf_x(i+ii,j,k)*tot_force/(tot_grad_x + 1e-12_WP) ! Summing Components
                     
   !                   ! Mark as Smoothed
   !                   x_smoothing_tracker(i+ii,j,k) = 0
   !                end do
   !                ! Reset Indicies
   !                negativeIndex = -1
   !                positiveIndex = 1
   !                tot_force = 0.0_WP
   !             endif
               
   !             ! Y Direction
   !             if(abs(fs%Pjy(i,j,k)) .gt. 1e-12 .and. y_smoothing_tracker(i,j,k) == 1) then  
   !                ! Current Cell
   !                tot_grad_y = grad_vf_y(i,j,k)
   !                tot_force = fs%Pjy(i,j,k)
   !                ! Going Down
   !                do while(abs(fs%Pjy(i,j+negativeIndex,k)) .gt. 1e-12 .and. (j+negativeIndex) .ge. cfg%jmin_)
   !                   tot_grad_y = tot_grad_y + grad_vf_y(i,j+negativeIndex,k)
   !                   tot_force = tot_force + fs%Pjy(i,j+negativeIndex,k)
   !                   negativeIndex = negativeIndex -1
   !                end do
   !                ! Going Up
   !                do while(abs(fs%Pjy(i,j+positiveIndex,k)) .gt. 1e-12 .and. (j+positiveIndex) .le. cfg%jmax_)
   !                   tot_grad_y = tot_grad_y + grad_vf_y(i,j+positiveIndex,k)
   !                   tot_force = tot_force + fs%Pjy(i,j+positiveIndex,k)
   !                   positiveIndex = positiveIndex +1
   !                end do
   !                ! Now Reassign Force Values in between 
   !                do jj = negativeIndex,positiveIndex
   !                   Fst_y(i,j+jj,k) = grad_vf_y(i,j+jj,k)*tot_force/(tot_grad_y + 1e-12_WP) ! Summing Components
   !                   ! Mark as Smoothed
   !                   y_smoothing_tracker(i,j+jj,k) = 0
   !                end do
   !                ! Reset Indicies
   !                negativeIndex = -1
   !                positiveIndex = 1
   !                tot_force = 0.0_WP
   !             endif 

   !          end do 
   !       end do
   !    end do
   !    ! write(*,'(A)') 'Gradient Smoothing Complete'
   !    fs%Pjx = Fst_x 
   !    fs%Pjy = Fst_y
   !    fs%Pjz = Fst_z
      
   ! end subroutine applyGradientSmoothing

   ! subroutine applyPoissonSmoothing(fs,vf)
   !    use irl_fortran_interface
   !    use f_PUNeigh_RectCub_class
   !    use f_SeparatorVariant_class
   !    use f_PUSolve_RectCub_class
   !    use hypre_str_class, only: pcg_smg,gmres

   !    class(tpns), intent(inout) :: fs ! Two Phase Flow Solver
   !    class(vfs), intent(inout) :: vf ! Volume Fraction Solver
   !    type(hypre_str) :: poisson_solver ! Poisson Solver Object

   !    integer :: i,j,k,ii,jj,kk
   !    real(WP) :: tot_grad_x,tot_grad_y


   !    ! Initialize Poisson Solver
   !    poisson_solver=hypre_str(cfg=cfg,name='Smoothing',method=gmres,nst=7)
   !    poisson_solver%maxit=fs%psolv%maxit
   !    poisson_solver%rcvg=fs%psolv%rcvg
   !    ! Setup Solver
   !    ! Set 7-pt stencil map for the pressure solver
   !    poisson_solver%stc(1,:)=[ 0, 0, 0]
   !    poisson_solver%stc(2,:)=[+1, 0, 0]
   !    poisson_solver%stc(3,:)=[-1, 0, 0]
   !    poisson_solver%stc(4,:)=[ 0,+1, 0]
   !    poisson_solver%stc(5,:)=[ 0,-1, 0]
   !    poisson_solver%stc(6,:)=[ 0, 0,+1]
   !    poisson_solver%stc(7,:)=[ 0, 0,-1]
   !    ! Set Up Laplacian Operator
   !    do k=fs%cfg%kmin_,fs%cfg%kmax_
   !       do j=fs%cfg%jmin_,fs%cfg%jmax_
   !          do i=fs%cfg%imin_,fs%cfg%imax_
   !             ! Set Laplacian
   !             poisson_solver%opr(1,i,j,k)=fs%divp_x(1,i,j,k)*fs%divu_x(-1,i+1,j,k)+&
   !             &                       fs%divp_x(0,i,j,k)*fs%divu_x( 0,i  ,j,k)+&
   !             &                       fs%divp_y(1,i,j,k)*fs%divv_y(-1,i,j+1,k)+&
   !             &                       fs%divp_y(0,i,j,k)*fs%divv_y( 0,i,j  ,k)+&
   !             &                       fs%divp_z(1,i,j,k)*fs%divw_z(-1,i,j,k+1)+&
   !             &                       fs%divp_z(0,i,j,k)*fs%divw_z( 0,i,j,k  )
   !             poisson_solver%opr(2,i,j,k)=fs%divp_x(1,i,j,k)*fs%divu_x( 0,i+1,j,k)
   !             poisson_solver%opr(3,i,j,k)=fs%divp_x(0,i,j,k)*fs%divu_x(-1,i  ,j,k)
   !             poisson_solver%opr(4,i,j,k)=fs%divp_y(1,i,j,k)*fs%divv_y( 0,i,j+1,k)
   !             poisson_solver%opr(5,i,j,k)=fs%divp_y(0,i,j,k)*fs%divv_y(-1,i,j  ,k)
   !             poisson_solver%opr(6,i,j,k)=fs%divp_z(1,i,j,k)*fs%divw_z( 0,i,j,k+1)
   !             poisson_solver%opr(7,i,j,k)=fs%divp_z(0,i,j,k)*fs%divw_z(-1,i,j,k  )
   !             ! Scale it by the cell volume
   !             ! poisson_solver%opr(:,i,j,k)=-poisson_solver%opr(:,i,j,k)*this%cfg%vol(i,j,k) I don't think I need this right now.
   !          end do
   !       end do
   !    end do
   !    ! Initialize the Poisson solver
   !    call poisson_solver%init()
   !    call poisson_solver%setup()
   !    poisson_solver%sol=fs%Pjx
   !    ! We now need to calculate the right hand side, which will be the divergence of the force field
   !    do k=fs%cfg%kmin_,fs%cfg%kmax_
   !       do j=fs%cfg%jmin_,fs%cfg%jmax_
   !          do i=fs%cfg%imin_,fs%cfg%imax_
   !             SurfaceTensionDiv(i,j,k) = (sum(fs%divp_x(:,i,j,k)*fs%Pjx(i:i+1,j,k))&
   !                                        +sum(fs%divp_y(:,i,j,k)*fs%Pjy(i,j:j+1,k))) ! No z direction right now
   !          end do
   !       end do
   !    end do
      
   !    poisson_solver%rhs = SurfaceTensionDiv
   !    poisson_source = SurfaceTensionDiv
   !    ! Solve the Poisson Equation
   !    call poisson_solver%solve()
   !    force_potential_field = poisson_solver%sol
   !    ! Now we need to calculate the gradient of the smoothed pressure field to get the new forces
   !    ! Calculate Forces
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             Fst_x(i,j,k) = (poisson_solver%sol(i,j,k) - poisson_solver%sol(i-1,j,k)) * cfg%dxi(i)
   !             Fst_y(i,j,k) = (poisson_solver%sol(i,j,k) - poisson_solver%sol(i,j-1,k)) * cfg%dyi(j)
   !             Fst_z(i,j,k) = 0.0_WP
   !          end do
   !       end do
   !    end do   
   !    ! Calculate Magnitude of Fst_x
   !    ! write(*,'(A,F10.5)') 'Max Fst_x before assign: ', maxval(abs(Fst_x))
   !    ! Assign
   !    fs%Pjx = Fst_x 
   !    fs%Pjy = Fst_y
   !    fs%Pjz = Fst_z
   !    call poisson_solver%destroy()
   !    ! write(*,'(A)') 'Poisson Smoothing Complete'
   ! end subroutine applyPoissonSmoothing
   
   subroutine updatePartitionOfUnity(vf,fs)
      use irl_fortran_interface
      use f_PUNeigh_RectCub_class
      use f_SeparatorVariant_class
      use f_PUSolve_RectCub_class

      implicit none
      class(vfs), intent(inout) :: vf ! Volume Fraction Solver
      class(tpns), intent(inout) :: fs ! Two Phase Flow Solver
      integer :: i,j,k,j_in,i_in,count ! Current Cell Location 
      type(PUNeigh_RectCub_type) :: neighborhood
      type(PU_RectCub_type) :: solver
      
      ! Temp Items
      type(SeparatorVariant_type) :: plane
      real(WP), dimension(1:3) :: cen,Gcen,Lcen,force,pos,tangentHolder
      integer, dimension(1:3) :: ind,indguess
      real(WP) :: dx,dy,xEval,yEval,zEval,Scale
      Scale = PuCfg%nx/cfg%nx
      ! write(*,'(A,3F35.10)') ' ============================================================================== Scale value: ', Scale 
      dx = PuCfg%dx(1)
      dy = PuCfg%dy(1)
      count = 0
      ! Create Neighborhood and solver
      call new(neighborhood) 
      call new(solver)
      ! Loop over real domain 
      do k=PuCfg%kmin_,PuCfg%kmax_
         do j=PuCfg%jmin_,PuCfg%jmax_
            do i=PuCfg%imin_,PuCfg%imax_
               ! write(*,'(A,3I10.5)') ' ==== location: ', i,j,k 
               ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then  
                  ! Get Points
                  xEval = PuCfg%xm(i)
                  yEval = PuCfg%ym(j)
                  zEval = PuCfg%zm(k)
                  pos = (/xEval,yEval,zEval/)
                  ! write(*,'(A,3F10.5)') ' ==== pos: ', pos
                  ! Get Global Location
                  indGuess = (/int(i/Scale),int(j/Scale),int(k/Scale)/)
                  ind = cfg%get_ijk_global(pos,indGuess)
                  ! write(*,'(A,3I10.5)') ' ==== ind value: ', ind
                  
                  ! Now we have the global index in ind, so we can move around based on that.  
                  ! Empty Neighborhood
                  call emptyNeighborhood(neighborhood)
                  ! Add 7x7 (3 on each side) stencil, in plane
                  ! cen = (/vf%cfg%xm(ind(1)),vf%cfg%ym(ind(2)),vf%cfg%zm(ind(3))/)
                  ! call addMember(neighborhood,cen,vf%liquid_gas_interface(ind(1),ind(2),ind(3)))
                  do j_in = -3,3
                     do i_in = -3,3
                        if(vf%VF(ind(1)+i_in,ind(2)+j_in,ind(3)) .gt. vf%VFmin .or. vf%VF(ind(1)+i_in,ind(2)+j_in,ind(3)) .lt. vf%VFmax) then
                           ! write(*,'(A,2I10.5)') ' ====i in, j in: ', i_in,j_in
                           ! For now, compute centroid as cell center
                           ! cen = (/vf%cfg%xm(ind(1)+i_in),vf%cfg%ym(ind(2)+j_in),vf%cfg%zm(ind(3))/)
                           ! ! write(*,'(A,2I10.5)') ' ====i in, j in: ', i_in,j_in
                           ! Gcen = vf%Gbary(:,ind(1)+i_in,ind(2)+j_in,ind(3))
                           ! Lcen = vf%Lbary(:,ind(1)+i_in,ind(2)+j_in,ind(3))
                           ! ! cen = vf%VF(ind(1)+i_in,ind(2)+j_in,ind(3)) * Lcen + (1-vf%VF(ind(1)+i_in,ind(2)+j_in,ind(3))) * Gcen 
                           ! cen = (Lcen+Gcen)/2
                           cen = calculateCentroid(vf%interface_polygon(1,ind(1)+i_in,ind(2)+j_in,ind(3)))
                           ! Get Plane
                           plane = vf%liquid_gas_interface(ind(1)+i_in,ind(2)+j_in,ind(3))
                           ! Add Plane to Neighborhood 
                           call addMember(neighborhood,cen,1.0_WP,plane)
                        endif
                     end do
                  end do
                  ! Set Neighborhood in solver 
                  call setNeighborhood(solver,neighborhood)
                  ! if(ind(1) .eq. 11 .and. ind(2) .eq. 22 .and. count .lt. 1) then 
                  !    write(*,'(A,F10.8)') '> Time  =  ', time%t
                  !    ! call printSolver(solver)
                  !    count = count + 1
                  ! endif
                  ! ====== Get Value
                  call getValue(solver,xEval,yEval,zEval,PU_spread*vf%cfg%dx(1),PartitionOfUnityValue(i,j,k))   
                  ! write(*,'(A,F10.8)') '> Radius  =  ', radius
                  ! call getValue(solver,xEval,yEval,zEval,0.25_WP,(/0.0_WP,0.0_WP,0.0_WP/),PartitionOfUnityValue(i,j,k))
                  call getTangent(solver,xEval,yEval,zEval,PU_spread*vf%cfg%dx(1),tangentHolder) 
                  ! call getTangent(solver,xEval,yEval,zEval,radius,center,tangentHolder)
                  PUTangent_x(i,j,k) = tangentHolder(1)
                  PUTangent_y(i,j,k) = tangentHolder(2)
                  PUTangent_z(i,j,k) = tangentHolder(3)
                  call getWeight(solver,xEval,yEval,zEval,PU_spread*vf%cfg%dx(1),PartitionOfUnityWeight(i,j,k)) 
                  ! call getWeight(solver,xEval,yEval,zEval,radius,center,PartitionOfUnityWeight(i,j,k))
                  ! write(*,'(A,F35.10)') ' ================= Partition of unity value: ', PartitionOfUnityValue(i,j,k)   
               ! endif 
            end do
         end do
      end do
   end subroutine updatePartitionOfUnity

   ! subroutine add_conservative_surface_tension_jump(this,dt,div,vf,contact_model)
   !    use messager,  only: die
   !    use vfs_class, only: vfs
   !    use irl_fortran_interface
   !    use f_PUSTNeigh_RectCub_class
   !    use f_SeparatorVariant_class 
   !    use f_PUSolve_RectCub_class

   !    implicit none
   !    class(tpns), intent(inout) :: this
   !    class(vfs), intent(inout) :: vf
   !    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    integer, intent(in), optional :: contact_model
   !    integer :: i,j,k
   !    real(WP) :: mycurv,mysurf
      
   !    ! Store old jump
   !    this%DPjx=this%Pjx
   !    this%DPjy=this%Pjy
   !    this%DPjz=this%Pjz
   !    ! Calculate pressure jump
   !    do k=this%cfg%kmin_,this%cfg%kmax_+1
   !       do j=this%cfg%jmin_,this%cfg%jmax_+1
   !          do i=this%cfg%imin_,this%cfg%imax_+1
   !             ! X face
   !             mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k))
   !             ! Y face
   !             mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*vf%VF(i,j-1:j,k))
   !             ! Z face
   !             mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k))
   !          end do
   !       end do
   !    end do

   !    PjxD = this%Pjx
   !    PjyD = this%Pjy 
   !    PjzD = this%Pjz

   !    ! Update Values
   !    call updateSurfaceTensionStresses(this,vf)
   !    call updateSurfaceTensionForces(this)
      
   !    SELECT CASE (SmoothingOption)
   !       CASE(1)
            
   !       CASE(2)
   !          call applyLaplacianSmoothing(this)
   !       CASE(3)
   !          call applyGradientSmoothing(this,vf)
   !       CASE(4)
   !          call applyPoissonSmoothing(this,vf)
   !       CASE DEFAULT

   !    END SELECT

   !    ! Add wall contact force to pressure jump 
   !    if (present(contact_model)) then
   !       select case (contact_model)
   !       case (1)
   !          call this%add_static_contact(vf=vf)
   !       case default
   !          call die('[tpns: add_surface_tension_jump] Unknown contact model!')
   !       end select
   !    end if
      
   !    ! Compute jump of DP
   !    this%DPjx=this%Pjx-this%DPjx
   !    this%DPjy=this%Pjy-this%DPjy
   !    this%DPjz=this%Pjz-this%DPjz
      
   !    ! Add div(Pjump) to RP
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             SurfaceTensionDiv(i,j,k) = dt*(sum(this%divp_x(:,i,j,k)*this%DPjx(i:i+1,j,k)/this%rho_U(i:i+1,j,k))&
   !             &                             +sum(this%divp_y(:,i,j,k)*this%DPjy(i,j:j+1,k)/this%rho_V(i,j:j+1,k)))
   !             ! &                           +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
   !             ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
   !             !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
   !             ! endif
   !             this%div(i,j,k)=this%div(i,j,k)-SurfaceTensionDiv(i,j,k)
   !          end do
   !       end do
   !    end do
   ! end subroutine add_conservative_surface_tension_jump
   
   ! subroutine add_CSF_Shift_surface_tension_jump(this,dt,div,vf,contact_model)
   !    use messager,  only: die
   !    use vfs_class, only: vfs
   !    use irl_fortran_interface
   !    use f_PUSTNeigh_RectCub_class
   !    use f_SeparatorVariant_class 
   !    use f_PUSolve_RectCub_class

   !    implicit none
   !    class(tpns), intent(inout) :: this
   !    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    class(vfs), intent(inout) :: vf
   !    integer, intent(in), optional :: contact_model
   !    integer :: i,j,k
   !    real(WP) :: mycurv,mysurf
   !    real(WP), dimension(3) :: OldMarangoniOption
   !    !! THIS IS WHERE WE GET THE CSF VALUES !!
   !    ! Store old jump
   !    this%DPjx=this%Pjx
   !    this%DPjy=this%Pjy
   !    this%DPjz=this%Pjz
   !    ! Calculate pressure jump
   !    do k=this%cfg%kmin_,this%cfg%kmax_+1
   !       do j=this%cfg%jmin_,this%cfg%jmax_+1
   !          do i=this%cfg%imin_,this%cfg%imax_+1
   !             ! X face
   !             mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k))
   !             ! Y face
   !             mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*vf%VF(i,j-1:j,k))
   !             ! Z face
   !             mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
   !             if (mysurf.gt.0.0_WP) then
   !                mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
   !             else
   !                mycurv=0.0_WP
   !             end if
   !             this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k))
   !          end do
   !       end do
   !    end do

   !    PjxD = this%Pjx
   !    PjyD = this%Pjy 
   !    PjzD = this%Pjz

   !    !! THIS IS WHERE WE GET THE PCST VALUES, WITH THE CURRENT MARANGONI OPTION
   !    ! Update Values
   !    call updateSurfaceTensionStresses(this,vf)
   !    call updateSurfaceTensionForces(this)
   !    ! Store Force Values
   !    Pjx_Marangoni = this%Pjx
   !    Pjy_Marangoni = this%Pjy
   !    Pjz_Marangoni = this%Pjz

   !    !! THIS IS WHERE WE GET THE PCST VALUES WITH CONSTANT SURFACE TENSION COEFFICIENT
   !    OldMarangoniOption = MarangoniOption
   !    MarangoniOption = (/0.0_WP,0.0_WP,0.0_WP/)
   !    call updateSurfaceTensionStresses(this,vf)
   !    call updateSurfaceTensionForces(this)
   !    MarangoniOption = OldMarangoniOption   

   !    ! Find the Difference caused by Marangoni
   !    this%Pjx = Pjx_Marangoni-this%Pjx
   !    this%Pjy = Pjy_Marangoni-this%Pjy
   !    this%Pjz = Pjz_Marangoni-this%Pjz
   !    ! We might need to smooth this so that the shape of the force field is the same as the CSF model
   !    SELECT CASE (SmoothingOption)
   !       CASE(1)
            
   !       CASE(2)
   !          call applyLaplacianSmoothing(this)
   !       CASE(3)
   !          call applyGradientSmoothing(this,vf)
   !       CASE(4)
   !          call applyPoissonSmoothing(this,vf)
   !       CASE DEFAULT

   !    END SELECT

   !    ! Now we have (F_PCST - F_PCST_CONST_ST). We now add this to the CSF values and use that
   !    this%Pjx = this%Pjx + PjxD
   !    this%Pjy = this%Pjy + PjyD
   !    this%Pjz = this%Pjz + PjzD

   !    ! Add wall contact force to pressure jump 
   !    if (present(contact_model)) then
   !       select case (contact_model)
   !       case (1)
   !          call this%add_static_contact(vf=vf)
   !       case default
   !          call die('[tpns: add_surface_tension_jump] Unknown contact model!')
   !       end select
   !    end if
      
   !    ! Compute jump of DP
   !    this%DPjx=this%Pjx-this%DPjx
   !    this%DPjy=this%Pjy-this%DPjy
   !    this%DPjz=this%Pjz-this%DPjz
      
   !    ! Add div(Pjump) to RP
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             SurfaceTensionDiv(i,j,k) = dt*(sum(this%divp_x(:,i,j,k)*this%DPjx(i:i+1,j,k)/this%rho_U(i:i+1,j,k))&
   !             &                             +sum(this%divp_y(:,i,j,k)*this%DPjy(i,j:j+1,k)/this%rho_V(i,j:j+1,k)))
   !             ! &                           +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
   !             ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
   !             !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
   !             ! endif
   !             this%div(i,j,k)=this%div(i,j,k)-SurfaceTensionDiv(i,j,k)
   !          end do
   !       end do
   !    end do
   ! end subroutine add_CSF_Shift_surface_tension_jump

   subroutine YoungsNormalMethod(vf,i,j,k,YoungNormal)
      use messager,  only: die
      use vfs_class, only: vfs

      implicit none
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,n
      ! Young Method Variables
      real(WP) :: Cbar_i,Cbar_ip1
      real(WP), dimension(1:3,1:4) :: normals
      real(WP), dimension(1:3) :: YoungNormal

      ! Top Right Corner
      Cbar_i = (vf%VF(i,j,k) + vf%VF(i,j+1,k))/2
      Cbar_ip1 = (vf%VF(i+1,j,k) + vf%VF(i+1,j+1,k))/2 
      normals(1,1) = vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

      Cbar_i = (vf%VF(i,j,k) + vf%VF(i+1,j,k))/2
      Cbar_ip1 = (vf%VF(i,j+1,k) + vf%VF(i+1,j+1,k))/2 
      normals(2,1) = vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
      
      normals(3,1) = 0.0_WP

      ! Bottom Right Corner
      Cbar_i = (vf%VF(i,j,k) + vf%VF(i,j-1,k))/2
      Cbar_ip1 = (vf%VF(i+1,j,k) + vf%VF(i+1,j-1,k))/2 
      normals(1,2) = vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

      Cbar_i = (vf%VF(i,j-1,k) + vf%VF(i+1,j-1,k))/2 
      Cbar_ip1 = (vf%VF(i,j,k) + vf%VF(i+1,j,k))/2
      normals(2,2) = vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
      
      normals(3,2) = 0.0_WP

      ! Top Left Corner
      Cbar_i = (vf%VF(i-1,j,k) + vf%VF(i-1,j+1,k))/2
      Cbar_ip1 = (vf%VF(i,j,k) + vf%VF(i,j+1,k))/2 
      normals(1,3) = vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

      Cbar_i = (vf%VF(i,j,k) + vf%VF(i-1,j,k))/2 
      Cbar_ip1 = (vf%VF(i,j+1,k) + vf%VF(i-1,j+1,k))/2
      normals(2,3) = vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
      
      normals(3,3) = 0.0_WP

      ! Bottom Left Corner
      Cbar_i = (vf%VF(i-1,j,k) + vf%VF(i-1,j-1,k))/2
      Cbar_ip1 = (vf%VF(i,j,k) + vf%VF(i,j-1,k))/2 
      normals(1,4) = vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

      Cbar_i = (vf%VF(i,j-1,k) + vf%VF(i-1,j-1,k))/2 
      Cbar_ip1 = (vf%VF(i,j,k) + vf%VF(i-1,j,k))/2
      normals(2,4) = vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
      
      normals(3,4) = 0.0_WP

      ! Combine into a single normal
      YoungNormal(:) = (normals(:,1)+normals(:,2)+normals(:,3)+normals(:,4))/4

      ! Normalize
      YoungNormal = YoungNormal/(sqrt(YoungNormal(1)*YoungNormal(1) + YoungNormal(2)*YoungNormal(2)))
   end subroutine YoungsNormalMethod

   subroutine CenteredColumnsNormalMethod(vf,i,j,k,normal)
      use messager,  only: die
      use vfs_class, only: vfs

      implicit none
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,n
      ! Young Method Variables
      real(WP) :: Xmx,Xmy,Ymx,Ymy,norm
      real(WP), dimension(1:3,1:2) :: normals
      real(WP), dimension(1:3) :: normal

      ! X Major Orientation
      ! First get Sign of mx
      Xmx = sign(1.0_WP,(vf%VF(i+1,j,k)-vf%VF(i-1,j,k)))
      ! Now get value of my
      Xmy = 0.5 * (sum(vf%VF(i-1:i+1,j+1,k))-sum(vf%VF(i-1:i+1,j-1,k)))
      ! Normalize
      norm = abs(Xmx)+abs(Xmy)
      Xmx = Xmx/norm
      Xmy = Xmy/norm 
      ! Store
      normals(1,1) = Xmx
      normals(2,1) = Xmy
      normals(3,1) = 0.0_WP

      ! Y Major Orientation
      ! First get Sign of my
      Ymy = sign(1.0_WP,(vf%VF(i,j+1,k)-vf%VF(i,j-1,k)))
      ! Now get value of mx
      Ymx = 0.5 * (sum(vf%VF(i+1,j-1:j+1,k))-sum(vf%VF(i-1,j-1:j+1,k)))

      ! Normalize
      norm = abs(Ymx)+abs(Ymy)
      Ymx = Ymx/norm
      Ymy = Ymy/norm 
      ! Store
      normals(1,2) = Ymx
      normals(2,2) = Ymy
      normals(3,2) = 1.0_WP

      ! Pick the Better Normal
      if(abs(Ymy) .gt. abs(Xmx)) then
         normal = normals(:,2)
      else
         normal = normals(:,1)
      endif
      ! Normalize
      normal = -normal
      ! normal = -normal/(sqrt(normal(1)*normal(1) + normal(2)*normal(2)))
      ! ! STOP
   end subroutine CenteredColumnsNormalMethod
      
   subroutine mixedYoungCenterNormal(vf,i,j,k,normal_return)
      use messager,  only: die
      use vfs_class, only: vfs

      implicit none
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,n
      ! Variables
      real(WP), dimension(1:3) :: YoungNormal,CCNormal
      real(WP) :: mYoung,mCC
      ! Return Value
      real(WP), dimension(1:3) :: normal_return

      ! ================================== Young's Method
      call YoungsNormalMethod(vf,i,j,k,YoungNormal)
      ! if(abs(YoungNormal(1)) .gt. abs(YoungNormal(2))) then ! X Oriented
      !    ! Normalize to have Unit x
      !    YoungNormal = YoungNormal/abs(YoungNormal(1))
      !    mYoung = YoungNormal(2)
      ! else
      !    YoungNormal = YoungNormal/abs(YoungNormal(2))
      !    mYoung = YoungNormal(1)
      ! endif
      ! ================================== CC Method
      call CenteredColumnsNormalMethod(vf,i,j,k,CCNormal)
      ! if(abs(CCNormal(1)) .gt. abs(CCNormal(2))) then ! X Oriented
      !    ! Normalize to have Unit x
      !    CCNormal = CCNormal/abs(CCNormal(1))
      !    mCC = CCNormal(2)
      ! else
      !    CCNormal = CCNormal/abs(CCNormal(2))
      !    mCC = CCNormal(1)
      ! endif

      ! ================================== Pick the one with smaller angular offset from aligned
      normal_return = YoungNormal ! Default
      if(CCNormal(3) .eq. 0.0_WP) then
         if(abs(CCNormal(1)) .le. abs(YoungNormal(1))) then
            normal_return = CCNormal
         endif
      else 
         if(abs(CCNormal(2)).le. abs(YoungNormal(1))) then
            normal_return = CCNormal
            normal_return(3) = 0.0_WP
         endif
      endif

      normal_return = normal_return/sqrt(normal_return(1)*normal_return(1)+normal_return(2)*normal_return(2))
   end subroutine mixedYoungCenterNormal

end module simulation
