!> AMR Rayleigh-Plateau instability

module simulation

   use precision,         only: WP
   use mathtools,         only: Pi
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrmpinc_class,    only: amrmpinc
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use amrio_class,       only: amrio
   use string,            only: str_medium
   use irl_fortran_interface
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Time integration
   type(timetracker) :: time

   ! Solver data
   type(amrmpinc), target :: fs
   type(amrdata) :: dQdt,Umag

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regridding
   type(event) :: regrid_evt

   ! Checkpoint/restart
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile,rpfile

   ! Rayleigh-Plateau parameters
   real(WP), dimension(3) :: center
   real(WP), dimension(3) :: initial_vel
   real(WP) :: radius
   real(WP) :: pert_amp
   integer  :: pert_mode

   ! Material properties
   real(WP) :: viscL_mol,viscG_mol

   ! Diagnostics
   real(WP) :: r_max,r_min,time_cap,time_nd

   real(WP), parameter, public :: VFlo=1.0e-12_WP
   real(WP), parameter, public :: VFhi=1.0_WP-VFlo

contains

   !====================================================================
   ! Modified Bessel functions of the first kind
   !====================================================================

   function bessel_i0(x) result(val)
      implicit none
      real(WP), intent(in) :: x
      real(WP) :: val
      real(WP) :: ax,y

      ax=abs(x)

      if (ax.lt.3.75_WP) then
         y=(x/3.75_WP)**2
         val=1.0_WP+y*(3.5156229_WP+y*(3.0899424_WP+y*(1.2067492_WP+ &
             y*(0.2659732_WP+y*(0.0360768_WP+y*0.0045813_WP)))))
      else
         y=3.75_WP/ax
         val=(exp(ax)/sqrt(ax))*(0.39894228_WP+y*(0.01328592_WP+ &
             y*(0.00225319_WP+y*(-0.00157565_WP+y*(0.00916281_WP+ &
             y*(-0.02057706_WP+y*(0.02635537_WP+y*(-0.01647633_WP+ &
             y*0.00392377_WP))))))))
      end if
   end function bessel_i0


   function bessel_i1(x) result(val)
      implicit none
      real(WP), intent(in) :: x
      real(WP) :: val
      real(WP) :: ax,y

      ax=abs(x)

      if (ax.lt.3.75_WP) then
         y=(x/3.75_WP)**2
         val=x*(0.5_WP+y*(0.87890594_WP+y*(0.51498869_WP+ &
             y*(0.15084934_WP+y*(0.02658733_WP+y*(0.00301532_WP+ &
             y*0.00032411_WP))))))
      else
         y=3.75_WP/ax
         val=(exp(ax)/sqrt(ax))*(0.39894228_WP+y*(-0.03988024_WP+ &
             y*(-0.00362018_WP+y*(0.00163801_WP+y*(-0.01031555_WP+ &
             y*(0.02282967_WP+y*(-0.02895312_WP+y*(0.01787654_WP- &
             y*0.00420059_WP))))))))
         if (x.lt.0.0_WP) val=-val
      end if
   end function bessel_i1


   !====================================================================
   ! Sinusoidally perturbed ligament
   !====================================================================

   function levelset_ligament(xyz,t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP) :: xrel,yrel,zrel,r,rloc,kz,zL

      zL=amr%zhi-amr%zlo

      xrel=xyz(1)-center(1)
      yrel=xyz(2)-center(2)
      zrel=xyz(3)-center(3)

      r=sqrt(xrel*xrel+yrel*yrel)

      kz=2.0_WP*Pi*real(pert_mode,WP)/zL
      rloc=radius*(1.0_WP+pert_amp*cos(kz*zrel))

      ! Positive inside liquid
      G=rloc-r
   end function levelset_ligament


   !====================================================================
   ! Linear inviscid Rayleigh-Plateau perturbation velocity
   !====================================================================

   subroutine rp_velocity(xyz,vel)
      implicit none
      real(WP), dimension(3), intent(in)  :: xyz
      real(WP), dimension(3), intent(out) :: vel

      real(WP), parameter :: tiny_r=1.0e-14_WP
      real(WP) :: xrel,yrel,zrel,r
      real(WP) :: kz,zL,kr0,I0kr0,I1kr0,I0kr,I1kr
      real(WP) :: growth,B

      zL=amr%zhi-amr%zlo
      kz=2.0_WP*Pi*real(pert_mode,WP)/zL

      kr0=kz*radius
      I0kr0=bessel_i0(kr0)
      I1kr0=bessel_i1(kr0)

      ! Linear inviscid Rayleigh-Plateau growth-rate parameter.
      growth=sqrt(max(0.0_WP, &
           (fs%sigma/(fs%rhoL*radius**3))*(I1kr0/I0kr0)*kr0*(1.0_WP-kr0**2)))

      ! Velocity-potential amplitude.
      if (abs(kz*I1kr0).gt.tiny_r) then
         B=pert_amp*radius*growth/(kz*I1kr0)
      else
         B=0.0_WP
      end if

      xrel=xyz(1)-center(1)
      yrel=xyz(2)-center(2)
      zrel=xyz(3)-center(3)
      r=sqrt(xrel*xrel+yrel*yrel)

      if (r.gt.tiny_r) then
         I1kr=bessel_i1(kz*r)
         vel(1)= B*kz*I1kr*(xrel/r)*cos(kz*zrel)
         vel(2)= B*kz*I1kr*(yrel/r)*cos(kz*zrel)
      else
         vel(1)=0.0_WP
         vel(2)=0.0_WP
      end if

      I0kr=bessel_i0(kz*r)
      vel(3)=-B*kz*I0kr*sin(kz*zrel)
   end subroutine rp_velocity


   !====================================================================
   ! Harmonic phase viscosity
   !====================================================================

   subroutine get_viscosity()
      use amrex_amr_module, only: amrex_mfiter,amrex_box

      implicit none

      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc
      real(WP), parameter :: myeps=1.0e-15_WP

      do lvl=0,amr%clvl()

         call amr%mfiter_build(lvl,mfi)

         do while (mfi%next())

            pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)

            bx=mfi%growntilebox(fs%nover)

            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)

                     pVisc(i,j,k,1)=1.0_WP/( &
                          pVF(i,j,k,1)/max(viscL_mol,myeps) + &
                          (1.0_WP-pVF(i,j,k,1))/max(viscG_mol,myeps))

                  end do
               end do
            end do

         end do

         call amr%mfiter_destroy(mfi)

      end do
   end subroutine get_viscosity


   !====================================================================
   ! AMR initialization callback
   !====================================================================

   subroutine ligament_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap, &
           amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpinc_class, only: VFlo

      implicit none

      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm

      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx

      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF=>null(), &
           pCL=>null(),pCG=>null(),pU=>null(),pV=>null(),pW=>null(),pQ=>null()

      real(WP), dimension(3) :: BL,BG
      real(WP), dimension(3) :: xyz,vtmp
      real(WP) :: dx,dy,dz,VF
      integer :: i,j,k
      integer, parameter :: nref=5

      dx=solver%amr%dx(lvl)
      dy=solver%amr%dy(lvl)
      dz=solver%amr%dz(lvl)

      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)

      do while (mfi%next())

         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         pU =>solver%U%mf(lvl)%dataptr(mfi)
         pV =>solver%V%mf(lvl)%dataptr(mfi)
         pW =>solver%W%mf(lvl)%dataptr(mfi)
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)

         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if

         bx=mfi%growntilebox(solver%nover)

         do k=bx%lo(3),bx%hi(3)
            do j=bx%lo(2),bx%hi(2)
               do i=bx%lo(1),bx%hi(1)

                  !-----------------------------------------------------
                  ! Cell-centered velocity Q
                  !-----------------------------------------------------
                  xyz=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx, &
                       solver%amr%ylo+(real(j,WP)+0.5_WP)*dy, &
                       solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]

                  call rp_velocity(xyz,vtmp)
                  pQ(i,j,k,1)=vtmp(1)
                  pQ(i,j,k,2)=vtmp(2)
                  pQ(i,j,k,3)=vtmp(3)

                  !-----------------------------------------------------
                  ! Staggered x-face velocity
                  !-----------------------------------------------------
                  xyz=[solver%amr%xlo+real(i,WP)*dx, &
                       solver%amr%ylo+(real(j,WP)+0.5_WP)*dy, &
                       solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]

                  call rp_velocity(xyz,vtmp)
                  pU(i,j,k,1)=vtmp(1)

                  !-----------------------------------------------------
                  ! Staggered y-face velocity
                  !-----------------------------------------------------
                  xyz=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx, &
                       solver%amr%ylo+real(j,WP)*dy, &
                       solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]

                  call rp_velocity(xyz,vtmp)
                  pV(i,j,k,1)=vtmp(2)

                  !-----------------------------------------------------
                  ! Staggered z-face velocity
                  !-----------------------------------------------------
                  xyz=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx, &
                       solver%amr%ylo+(real(j,WP)+0.5_WP)*dy, &
                       solver%amr%zlo+real(k,WP)*dz]

                  call rp_velocity(xyz,vtmp)
                  pW(i,j,k,1)=vtmp(3)

                  !-----------------------------------------------------
                  ! Volume fraction and phase barycenters
                  !-----------------------------------------------------
                  call initialize_volume_moments( &
                       lo=[solver%amr%xlo+real(i  ,WP)*dx, &
                           solver%amr%ylo+real(j  ,WP)*dy, &
                           solver%amr%zlo+real(k  ,WP)*dz], &
                       hi=[solver%amr%xlo+real(i+1,WP)*dx, &
                           solver%amr%ylo+real(j+1,WP)*dy, &
                           solver%amr%zlo+real(k+1,WP)*dz], &
                       levelset=levelset_ligament,time=time,level=nref, &
                       VFlo=VFlo,VF=VF,BL=BL,BG=BG)

                  pVF(i,j,k,1)=VF

                  if (lvl.eq.solver%amr%maxlvl) then
                     pCL(i,j,k,:)=BL
                     pCG(i,j,k,:)=BG
                  end if

               end do
            end do
         end do

      end do

      call amrex_mfiter_destroy(mfi)

   end subroutine ligament_init


   !====================================================================
   ! Rayleigh-Plateau radius diagnostics
   !====================================================================

   subroutine compute_rp_stats()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX,MPI_MIN

      implicit none

      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
      real(WP) :: dx,dy,dz,x,y,r
      logical :: found_interface
      real(WP), dimension(3) :: lo,hi
      type(RectCub_type) :: cell
      type(PlanarSep_type) :: planar_sep
      type(Poly_type) :: polygon
      real(WP), dimension(1:3) :: cen
      
      lvl=amr%maxlvl
      dx=amr%dx(lvl)
      dy=amr%dy(lvl)
      dz=amr%dz(lvl)

      r_max=-huge(1.0_WP)
      r_min= huge(1.0_WP)
      found_interface=.false.

      call amr%mfiter_build(lvl,mfi)

      do while (mfi%next())

         pVF=>fs%VF%mf(lvl)%dataptr(mfi)
         pPLIC=>fs%PLIC%dataptr(mfi)
         bx=mfi%tilebox()

         do k=bx%lo(3),bx%hi(3)
            do j=bx%lo(2),bx%hi(2)
               do i=bx%lo(1),bx%hi(1)

                  if (pVF(i,j,k,1).le.VFlo .or. pVF(i,j,k,1).ge.VFhi) cycle

                  call setNumberOfPlanes(planar_sep,1)
                  call setPlane(planar_sep,0,pPLIC(i,j,k,1:3),pPLIC(i,j,k,4))
                  lo=[amr%xlo+real(i  ,WP)*dx,amr%ylo+real(j  ,WP)*dy,amr%zlo+real(k  ,WP)*dz]
                  hi=[amr%xlo+real(i+1,WP)*dx,amr%ylo+real(j+1,WP)*dy,amr%zlo+real(k+1,WP)*dz]
                  call construct_2pt(cell,lo,hi)
                  call getPoly(cell,planar_sep,0,polygon)
                  cen = calculateCentroid(polygon)

                  x=amr%xlo+(real(i,WP)+0.5_WP)*dx-center(1)
                  y=amr%ylo+(real(j,WP)+0.5_WP)*dy-center(2)
                  r=sqrt(x*x+y*y)

                  r_max=max(r_max,r)
                  r_min=min(r_min,r)
                  found_interface=.true.

               end do
            end do
         end do

      end do

      call amr%mfiter_destroy(mfi)

      call MPI_ALLREDUCE(MPI_IN_PLACE,r_max,1,MPI_REAL_WP,MPI_MAX,amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,r_min,1,MPI_REAL_WP,MPI_MIN,amr%comm,ierr)

      if (r_max.lt.0.0_WP .or. r_min.gt.0.5_WP*huge(1.0_WP)) then
         r_max=0.0_WP
         r_min=0.0_WP
      end if

      time_nd=time%t/time_cap

   end subroutine compute_rp_stats


   !====================================================================
   ! Initialization
   !====================================================================

   subroutine simulation_init()
      use param, only: param_read
      use amrmg_class, only: amrmg_outer_pcg_mlmg

      implicit none

      !---------------------------------------------------------------
      ! Create AMR grid
      !---------------------------------------------------------------
      create_amrgrid: block

         amr%name='rayleigh_plateau'

         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)

         ! Old Rayleigh-Plateau domain: Lx = Ly = Lz = 1.
         ! Grid resolution itself is read using the new AMReX input format.
         amr%xlo=-0.5_WP; amr%xhi=+0.5_WP
         amr%ylo=-0.5_WP; amr%yhi=+0.5_WP
         amr%zlo=-0.5_WP; amr%zhi=+0.5_WP

         amr%xper=.true.
         amr%yper=.true.
         amr%zper=.true.

         call param_read('Max level',amr%maxlvl)
         call param_read('Blocking factor',amr%nbloc)
         call param_read('Max grid size',amr%nmax)

         call amr%initialize()

      end block create_amrgrid


      !---------------------------------------------------------------
      ! Restart/checkpoint setup
      !---------------------------------------------------------------
      handle_restart: block

         integer :: restart_step

         call io%initialize(amr=amr,nfiles=1)

         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)

         if (restarted) then
            call io%read_header(dirname=trim(restart_dir), &
                 time=restart_time,step=restart_step)
         end if

      end block handle_restart


      !---------------------------------------------------------------
      ! Time integration
      !---------------------------------------------------------------
      initialize_time: block

         time=timetracker(amRoot=amr%amRoot,name='Rayleigh Plateau')

         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)

         time%dt=time%dtmax

         call param_read('Subiterations',time%itmax,default=2)

         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if

      end block initialize_time


      !---------------------------------------------------------------
      ! Create flow solver
      !---------------------------------------------------------------
      create_flow_solver: block

         call fs%initialize(amr,name='rayleigh_plateau')

         ! Rayleigh-Plateau geometry
         call param_read('Ligament radius',radius)
         call param_read('Ligament center',center,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Perturbation amplitude',pert_amp)
         call param_read('Perturbation mode',pert_mode,default=1)
         call param_read('Initial velocity',initial_vel, &
              default=[0.0_WP,0.0_WP,0.0_WP])

         ! Material properties: preserve the old Rayleigh-Plateau inputs.
         call param_read('Liquid density',fs%rhoL)
         call param_read('Gas density',fs%rhoG)
         call param_read('Liquid dynamic viscosity',viscL_mol)
         call param_read('Gas dynamic viscosity',viscG_mol)
         call param_read('Surface tension coefficient',fs%sigma)

         ! Old capillary time definition.
         time_cap=sqrt(fs%rhoL*radius**3/fs%sigma)

         ! User initialization callback.
         fs%user_init=>ligament_init

         ! Pressure solver setup following the new AMReX example.
         fs%psolver%outer_solver=amrmg_outer_pcg_mlmg
         fs%psolver%tol_rel=1.0e-5_WP

      end block create_flow_solver


      !---------------------------------------------------------------
      ! Workspace
      !---------------------------------------------------------------
      create_workspace: block
         use amrdata_class, only: interp_none

         call dQdt%initialize(amr,name='dQdt',ncomp=3,ng=0,interp=interp_none)
         call dQdt%register()

         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none)
         call Umag%register()

      end block create_workspace


      !---------------------------------------------------------------
      ! Regridding and initial hierarchy
      !---------------------------------------------------------------
      init_regridding: block

         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)

         if (restarted) then

            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)

         else

            call amr%init_from_scratch(time=time%t)

#ifdef USE_IRL
            call fs%build_ppic(time%t)
#else
            call fs%build_plic(time%t)
#endif

            call fs%build_subVF()

            ! ligament_init already initializes face velocities and Q.
            ! Recompute/fill face velocity consistently with the AMR hierarchy.
            call fs%get_face_velocity()
            call fs%average_down_velocity()
            call fs%fill_velocity(time=time%t)

         end if

         call get_viscosity()

         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q, &
              compX=1,compY=2,compZ=3)

      end block init_regridding


      !---------------------------------------------------------------
      ! Checkpointing
      !---------------------------------------------------------------
      init_checkpoint: block

         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)

         call fs%register_checkpoint(io)
         call io%add_scalar(name='dt',value=time%dt)

      end block init_checkpoint


      !---------------------------------------------------------------
      ! Visualization
      !---------------------------------------------------------------
      create_visualization: block

         call viz%initialize(amr,'rayleigh_plateau',use_hdf5=.false.)

         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%Q,1,'U')
         call viz%add_scalar(fs%Q,2,'V')
         call viz%add_scalar(fs%Q,3,'W')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%VF,1,'VF')

         call fs%smesh%write_as_vtu()
         call viz%add_surfmesh(fs%smesh,'interface')

         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)

         if (viz_evt%occurs()) call viz%write(time=time%t)

      end block create_visualization


      !---------------------------------------------------------------
      ! Monitors
      !---------------------------------------------------------------
      create_monitor: block

         call fs%get_info()
         call fs%get_cfl(time%dt,time%cfl)

         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(fs%CFL,'CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%psolver%res,'Pressure residual')
         call mfile%add_column(fs%psolver%niter,'Pressure iterations')
         call mfile%add_column(fs%divmax,'Divergence')
         call mfile%write()

         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLst,'CFLst')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()

         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()

         rpfile=monitor(amRoot=amr%amRoot,name='rayleigh_plateau')
         call rpfile%add_column(time%n,'Timestep')
         call rpfile%add_column(time%t,'Time')
         call rpfile%add_column(time_nd,'t/tc')
         call rpfile%add_column(r_max,'Rmax')
         call rpfile%add_column(r_min,'Rmin')

         call compute_rp_stats()
         call rpfile%write()

      end block create_monitor

   end subroutine simulation_init


   !====================================================================
   ! Run simulation
   !====================================================================

   subroutine simulation_run()

      implicit none

      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old interface and velocity fields
         call fs%store_old()

         ! Prepare surface output
         if (viz_evt%occurs()) fs%update_smesh=.true.

         !------------------------------------------------------------
         ! Subiterations
         !------------------------------------------------------------
         do while (time%it.le.time%itmax)

            ! Mid-time velocity
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Advection + viscous terms
            call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)

            call fs%Q%average_down()
            call fs%Q%fill(time%t)

            ! Rebuild interface
#ifdef USE_IRL
            call fs%build_ppic(time%t)
#else
            call fs%build_plic(time%t)
#endif

            call fs%build_subVF()

            ! Cell-centered -> face velocity
            call fs%get_face_velocity()

            ! Current pressure term
            call fs%add_pressure(scale=time%dt,phi=fs%P)

            ! Surface tension
            call fs%add_surface_tension(scale=time%dt)

            ! Fill hierarchy
            call fs%Q%average_down()
            call fs%Q%fill(time=time%t)

            call fs%average_down_velocity()
            call fs%fill_velocity(time=time%t)

            ! Pressure projection
            call fs%get_div()
            call fs%div%mult(val=1.0_WP/time%dt)

            call fs%prepare_psolver()
            call fs%psolver%solve(rhs=fs%div)

            call fs%add_pressure(scale=time%dt)
            call fs%P%add(src=fs%psolver%sol)

            call fs%Q%average_down()
            call fs%Q%fill(time=time%t)

            call fs%average_down_velocity()
            call fs%fill_velocity(time=time%t)

            time%it=time%it+1

         end do

         !------------------------------------------------------------
         ! Regrid
         !------------------------------------------------------------
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Update viscosity after interface motion/regridding
         call get_viscosity()

         ! Velocity magnitude
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q, &
              compX=1,compY=2,compZ=3)

         ! Monitors
         call fs%get_info()
         call mfile%write()
         call cflfile%write()

         call compute_rp_stats()
         call rpfile%write()

         ! Visualization
         if (viz_evt%occurs()) then
            call viz%write(time=time%t)
            fs%update_smesh=.false.
         end if

         ! Checkpoint
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write( &
                    dirname='restart/rayleigh_plateau_'// &
                    trim(adjustl(rtoa(time%t))), &
                    time=time%t,step=time%n)
            end block save_checkpoint
         end if

      end do

   end subroutine simulation_run


   !====================================================================
   ! Finalization
   !====================================================================

   subroutine simulation_final()

      implicit none

      call time%finalize()

      call amr%finalize()
      call regrid_evt%finalize()

      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()

      call viz%finalize()
      call viz_evt%finalize()

      call save_evt%finalize()
      call io%finalize()

      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
      call rpfile%finalize()

   end subroutine simulation_final

end module simulation