module conservative_st
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
    implicit none
    private

    public :: temp,conservative_st_type

    ! Parameters
    real(WP), parameter :: VFlo=1.0e-10_WP                         !< Minimum VF value considered
    real(WP), parameter :: VFhi=1.0_WP-VFlo                        !< Maximum VF value considered
    !> Conservative Surface Tension Type
    type :: conservative_st_type
        ! Things we store here
        ! Options
        real(WP), dimension(3) :: MarangoniOption
        integer :: SurfaceTensionOption,CurvatureOption,SmoothingOption
        real(WP) :: PressureOption
        ! Solvers
        type(tpns),pointer,public :: fs => null()
        type(vfs),pointer,public :: vf => null()
        type(timetracker),pointer,public :: time => null()
        ! Working and Display Arrays
        real(WP) :: PU_spread ! How far the PU spreads
        real(WP), dimension(:,:,:), allocatable :: sigma_xx,sigma_yy,sigma_xy,sigma_yx ! Used Stress Tensor components
        real(WP), dimension(:,:,:), allocatable :: sigma_xx_P,sigma_yy_P,sigma_xy_P,sigma_yx_P ! Pressure Component
        real(WP), dimension(:,:,:), allocatable :: sigma_xx_NoP,sigma_yy_NoP,sigma_xy_NoP,sigma_yx_NoP ! Tangent Component
        real(WP), dimension(:,:,:), allocatable :: force_potential_field,poisson_source ! Potential Field Tools
        real(WP), dimension(:,:,:), allocatable :: Fst_x,Fst_y,Fst_z ! Forces sent to momentum equation
        real(WP), dimension(:,:,:), allocatable :: PjxD,PjyD,PjzD ! CSF Forces
        real(WP), dimension(:,:,:), allocatable :: Pjx_ST,Pjy_ST,Pjz_ST ! Forces from UpdateSTForces
        real(WP), dimension(:,:,:), allocatable :: Pjx_Marangoni,Pjy_Marangoni,Pjz_Marangoni ! Storage
        real(WP), dimension(:,:,:), allocatable :: grad_vf_x, grad_vf_y,grad_vf_z ! Gradient of volume fraction, but absolute value
        integer, dimension(:,:,:), allocatable :: x_smoothing_tracker,y_smoothing_tracker ! Track where smoothing is applied
        real(WP), dimension(:,:,:), allocatable :: SurfaceTensionDiv ! Divergence of surface tension force
        ! Initailize Boolean
        logical :: initialized = .false.
    contains
        procedure :: temp 
        procedure :: init
        procedure :: get_dmomdt
        procedure :: add_surface_tension_jump
        ! procedure :: end ! this is to deallocate variables, clear pointers, etc.

        !Private Functions
        procedure, private :: get_dmomdt_full_ST
        procedure, private :: get_dmomdt_no_ST

        procedure, private :: add_surface_tension_jump_full_ST
        procedure, private :: add_surface_tension_jump_no_ST
        procedure, private :: add_conservative_surface_tension_jump
        procedure, private :: add_CSF_Shift_surface_tension_jump

        procedure, private :: updateSurfaceTensionStresses
        procedure, private :: updateSurfaceTensionForces

        procedure, private :: applyLaplacianSmoothing
        procedure, private :: applyGradientSmoothing
        procedure, private :: applyPoissonSmoothing

    end type conservative_st_type
contains 

subroutine temp(this)
    class(conservative_st_type), intent(inout) :: this
    print *, "This is a temporary subroutine in the conservative_st module."
end subroutine temp

subroutine init(this,fs_in,vf_in,time_in)
    class(conservative_st_type), intent(inout) :: this
    class(tpns),target, intent(in) :: fs_in
    class(vfs),target, intent(in) :: vf_in
    class(timetracker),target, intent(in) :: time_in
    ! Assign fs_in to fs
    this%fs => fs_in
    this%vf => vf_in
    this%time => time_in
    ! Allocate arrays
    allocate(this%sigma_xx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_xy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%sigma_xx_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yy_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_xy_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yx_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%sigma_xx_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yy_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_xy_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yx_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Fst_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%PjxD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PjyD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PjzD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Pjx_ST(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjy_ST(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjz_ST(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Pjx_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjy_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjz_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%SurfaceTensionDiv(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%grad_vf_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%grad_vf_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%grad_vf_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%x_smoothing_tracker(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%y_smoothing_tracker(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%force_potential_field(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%poisson_source(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
end subroutine init 

subroutine get_dmomdt(this,resU,resV,resW)
    class(conservative_st_type), intent(inout) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)

    SELECT CASE (this%SurfaceTensionOption)
        CASE (2)
            ! write(*,'(A)') 'CASE 2 ================================================'
            call this%add_surface_tension_jump_no_ST(dt=this%time%dt,div=this%fs%div)
            call this%get_dmomdt_full_ST(resU,resV,resW)
        CASE (3)
            ! write(*,'(A)') 'CASE 3 ================================================'
            call this%get_dmomdt_no_ST(resU,resV,resW)
        CASE (4)
            ! Conservative Surface Tension Approach for get_dmomdt
            call this%fs%get_dmomdt(resU,resV,resW)
        CASE DEFAULT
            ! write(*,'(A)') 'CASE DEFAULT ================================================'
            call this%fs%get_dmomdt(resU,resV,resW)
    END SELECT
end subroutine get_dmomdt 


subroutine add_surface_tension_jump(this)
    class(conservative_st_type), intent(inout) :: this

    SELECT CASE (this%SurfaceTensionOption)
        CASE (2)
            call this%add_surface_tension_jump_no_ST(dt=this%time%dt,div=this%fs%div)
        CASE (3)
            call this%add_surface_tension_jump_full_ST(dt=this%time%dt,div=this%fs%div)
        CASE (4) 
            call this%add_conservative_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        CASE DEFAULT
            call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
    END SELECT
end subroutine add_surface_tension_jump


! PRIVATE FUNCTIONS
! ########################################################### MOMENTUM FUNCTIONS #############################################
subroutine get_dmomdt_full_ST(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(conservative_st_type), intent(inout) :: this
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk

      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
      ! Flux of rhoU
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FUX(i,j,k)=this%fs%FUX(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybu_x(:,i,j,k)*this%fs%rho_Uold(i:i+1,j,k))*sum(this%fs%hybu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))*sum(this%fs%itpu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))+sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FUY(i,j,k)=this%fs%FUY(i,j,k) &
               &              -sum(this%fs%hybu_y(:,i,j,k)*this%fs%rho_Uold(i,j-1:j,k))*sum(this%fs%hybu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))*sum(this%fs%itpv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))+sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FUZ(i,j,k)=this%fs%FUZ(i,j,k) &
               &              -sum(this%fs%hybu_z(:,i,j,k)*this%fs%rho_Uold(i,j,k-1:k))*sum(this%fs%hybu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))*sum(this%fs%itpw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))+sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoUdt(i,j,k)=sum(this%fs%divu_x(:,i,j,k)*this%fs%FUX(i-1:i,j,k))+&
               &              sum(this%fs%divu_y(:,i,j,k)*this%fs%FUY(i,j:j+1,k))+&
               &              sum(this%fs%divu_z(:,i,j,k)*this%fs%FUZ(i,j,k:k+1))+this%fs%Pjx(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FVX(i,j,k)=this%fs%FVX(i,j,k) &
               &              -sum(this%fs%hybv_x(:,i,j,k)*this%fs%rho_Vold(i-1:i,j,k))*sum(this%fs%hybv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))*sum(this%fs%itpu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))+sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FVY(i,j,k)=this%fs%FVY(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybv_y(:,i,j,k)*this%fs%rho_Vold(i,j:j+1,k))*sum(this%fs%hybv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))*sum(this%fs%itpv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))+sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FVZ(i,j,k)=this%fs%FVZ(i,j,k) &
               &              -sum(this%fs%hybv_z(:,i,j,k)*this%fs%rho_Vold(i,j,k-1:k))*sum(this%fs%hybv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))*sum(this%fs%itpw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))+sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoVdt(i,j,k)=sum(this%fs%divv_x(:,i,j,k)*this%fs%FVX(i:i+1,j,k))+&
               &              sum(this%fs%divv_y(:,i,j,k)*this%fs%FVY(i,j-1:j,k))+&
               &              sum(this%fs%divv_z(:,i,j,k)*this%fs%FVZ(i,j,k:k+1))+this%fs%Pjy(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FWX(i,j,k)=this%fs%FWX(i,j,k) &
               &              -sum(this%fs%hybw_x(:,i,j,k)*this%fs%rho_Wold(i-1:i,j,k))*sum(this%fs%hybw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))*sum(this%fs%itpu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))+sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FWY(i,j,k)=this%fs%FWY(i,j,k) &
               &              -sum(this%fs%hybw_y(:,i,j,k)*this%fs%rho_Wold(i,j-1:j,k))*sum(this%fs%hybw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))*sum(this%fs%itpv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))+sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FWZ(i,j,k)=this%fs%FWZ(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybw_z(:,i,j,k)*this%fs%rho_Wold(i,j,k:k+1))*sum(this%fs%hybw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))*sum(this%fs%itpw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))+sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)))
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoWdt(i,j,k)=sum(this%fs%divw_x(:,i,j,k)*this%fs%FWX(i:i+1,j,k))+&
               &              sum(this%fs%divw_y(:,i,j,k)*this%fs%FWY(i,j:j+1,k))+&
               &              sum(this%fs%divw_z(:,i,j,k)*this%fs%FWZ(i,j,k-1:k))+this%fs%Pjz(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoWdt)
   end subroutine get_dmomdt_full_ST

   subroutine get_dmomdt_no_ST(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(conservative_st_type), intent(inout) :: this
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk

      
      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
      ! Flux of rhoU
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FUX(i,j,k)=this%fs%FUX(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybu_x(:,i,j,k)*this%fs%rho_Uold(i:i+1,j,k))*sum(this%fs%hybu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))*sum(this%fs%itpu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))+sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FUY(i,j,k)=this%fs%FUY(i,j,k) &
               &              -sum(this%fs%hybu_y(:,i,j,k)*this%fs%rho_Uold(i,j-1:j,k))*sum(this%fs%hybu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))*sum(this%fs%itpv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))+sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FUZ(i,j,k)=this%fs%FUZ(i,j,k) &
               &              -sum(this%fs%hybu_z(:,i,j,k)*this%fs%rho_Uold(i,j,k-1:k))*sum(this%fs%hybu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))*sum(this%fs%itpw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))+sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoUdt(i,j,k)=sum(this%fs%divu_x(:,i,j,k)*this%fs%FUX(i-1:i,j,k))+&
               &              sum(this%fs%divu_y(:,i,j,k)*this%fs%FUY(i,j:j+1,k))+&
               &              sum(this%fs%divu_z(:,i,j,k)*this%fs%FUZ(i,j,k:k+1))+this%fs%Pjx(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FVX(i,j,k)=this%fs%FVX(i,j,k) &
               &              -sum(this%fs%hybv_x(:,i,j,k)*this%fs%rho_Vold(i-1:i,j,k))*sum(this%fs%hybv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))*sum(this%fs%itpu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))+sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FVY(i,j,k)=this%fs%FVY(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybv_y(:,i,j,k)*this%fs%rho_Vold(i,j:j+1,k))*sum(this%fs%hybv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))*sum(this%fs%itpv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))+sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FVZ(i,j,k)=this%fs%FVZ(i,j,k) &
               &              -sum(this%fs%hybv_z(:,i,j,k)*this%fs%rho_Vold(i,j,k-1:k))*sum(this%fs%hybv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))*sum(this%fs%itpw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))+sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoVdt(i,j,k)=sum(this%fs%divv_x(:,i,j,k)*this%fs%FVX(i:i+1,j,k))+&
               &              sum(this%fs%divv_y(:,i,j,k)*this%fs%FVY(i,j-1:j,k))+&
               &              sum(this%fs%divv_z(:,i,j,k)*this%fs%FVZ(i,j,k:k+1))+this%fs%Pjy(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FWX(i,j,k)=this%fs%FWX(i,j,k) &
               &              -sum(this%fs%hybw_x(:,i,j,k)*this%fs%rho_Wold(i-1:i,j,k))*sum(this%fs%hybw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))*sum(this%fs%itpu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))+sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FWY(i,j,k)=this%fs%FWY(i,j,k) &
               &              -sum(this%fs%hybw_y(:,i,j,k)*this%fs%rho_Wold(i,j-1:j,k))*sum(this%fs%hybw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))*sum(this%fs%itpv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))+sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FWZ(i,j,k)=this%fs%FWZ(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybw_z(:,i,j,k)*this%fs%rho_Wold(i,j,k:k+1))*sum(this%fs%hybw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))*sum(this%fs%itpw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))+sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)))
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoWdt(i,j,k)=sum(this%fs%divw_x(:,i,j,k)*this%fs%FWX(i:i+1,j,k))+&
               &              sum(this%fs%divw_y(:,i,j,k)*this%fs%FWY(i,j:j+1,k))+&
               &              sum(this%fs%divw_z(:,i,j,k)*this%fs%FWZ(i,j,k-1:k)) ! +this%fs%Pjz(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoWdt)
      
   end subroutine get_dmomdt_no_ST

! ############################################################### SURFACE TENSION JUMP SUBROUTINES ##########################################################
subroutine add_surface_tension_jump_full_ST(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do
    
    ! Add wall contact force to pressure jump
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%Pjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                        +sum(this%fs%divp_y(:,i,j,k)*this%fs%Pjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                        +sum(this%fs%divp_z(:,i,j,k)*this%fs%Pjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
        end do
        end do
    end do
end subroutine add_surface_tension_jump_full_ST

subroutine add_surface_tension_jump_no_ST(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    implicit none
    
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_,this%fs%cfg%jmino_,this%fs%cfg%kmino_) :: addToDiv
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf, magx,magy,magz,magdiv
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do
    
    ! Add wall contact force to pressure jump
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Calculate Magnitude of DPs
    magx = NORM2(this%fs%DPjx)
    magy = NORM2(this%fs%DPjy)
    magz = NORM2(this%fs%DPjz)
    magdiv = NORM2(this%fs%div)
    addToDiv = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                        +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                        +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
    
    magz = NORM2(addToDiv)
    write(*,'(A,4F35.10)') '   magx,magy,magz,magdiv: ', magx,magy,magz,magdiv
    if(magz.gt.0) then
        write(*,'(A,F10.5)')  '====================================MagToDiv: ',magz
    endif
end subroutine add_surface_tension_jump_no_ST

subroutine add_conservative_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUSTNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUSolve_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    ! Update Values
    call this%updateSurfaceTensionStresses()
    call this%updateSurfaceTensionForces()
    
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call this%applyLaplacianSmoothing()
        CASE(3)
        call this%applyGradientSmoothing()
        CASE(4)
        call this%applyPoissonSmoothing()
        CASE DEFAULT

    END SELECT

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                                  +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k)))
            ! &                           +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)-this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_conservative_surface_tension_jump
   
subroutine add_CSF_Shift_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUSTNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUSolve_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    real(WP), dimension(3) :: OldMarangoniOption
    !! THIS IS WHERE WE GET THE CSF VALUES !!
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    !! THIS IS WHERE WE GET THE PCST VALUES, WITH THE CURRENT MARANGONI OPTION
    ! Update Values
    call updateSurfaceTensionStresses(this)
    call updateSurfaceTensionForces(this)
    ! Store Force Values
    this%Pjx_Marangoni = this%fs%Pjx
    this%Pjy_Marangoni = this%fs%Pjy
    this%Pjz_Marangoni = this%fs%Pjz

    !! THIS IS WHERE WE GET THE PCST VALUES WITH CONSTANT SURFACE TENSION COEFFICIENT
    OldMarangoniOption = this%MarangoniOption
    this%MarangoniOption = (/0.0_WP,0.0_WP,0.0_WP/)
    call updateSurfaceTensionStresses(this)
    call updateSurfaceTensionForces(this)
    this%MarangoniOption = OldMarangoniOption   

    ! Find the Difference caused by Marangoni
    this%fs%Pjx = this%Pjx_Marangoni-this%fs%Pjx
    this%fs%Pjy = this%Pjy_Marangoni-this%fs%Pjy
    this%fs%Pjz = this%Pjz_Marangoni-this%fs%Pjz
    ! We might need to smooth this so that the shape of the force field is the same as the CSF model
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call applyLaplacianSmoothing(this)
        CASE(3)
        call applyGradientSmoothing(this)
        CASE(4)
        call applyPoissonSmoothing(this)
        CASE DEFAULT

    END SELECT

    ! Now we have (F_PCST - F_PCST_CONST_ST). We now add this to the CSF values and use that
    this%fs%Pjx = this%fs%Pjx + this%PjxD
    this%fs%Pjy = this%fs%Pjy + this%PjyD
    this%fs%Pjz = this%fs%Pjz + this%PjzD

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                             +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k)))
            ! &                           +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%fs%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)-this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_CSF_Shift_surface_tension_jump

! ############################################# PartitionOfUnity Functions
subroutine updateSurfaceTensionStresses(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUSolve_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PU_RectCub_type) :: solver
    
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,C,xEval,yEval,zEval

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    C = 4.0_WP
    ! Create Neighborhood and solver
    call new(neighborhood) 
    call new(solver)
    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then

                ! Empty Neighborhood
                call emptyNeighborhood(neighborhood) 
                ! Add 7x7 (3 on each side) stencil, in plane
                ! write(*,'(A)') '=============== New: '
                do j_in = -3,3
                    do i_in = -3,3
                    if(this%vf%VF(i+i_in,j+j_in,k) .gt. this%vf%VFmin .and. this%vf%VF(i+i_in,j+j_in,k) .lt. this%vf%VFmax) then
                        ! For now, compute centroid as cell center
                        ! Gcen = vf%Gbary(:,i+i_in,j+j_in,k)
                        ! Lcen = vf%Lbary(:,i+i_in,j+j_in,k)
                        ! ! cen = vf%VF(i+i_in,j+j_in,k) * Lcen + (1-vf%VF(i+i_in,j+j_in,k)) * Gcen
                        ! ! cen = (Lcen+Gcen)/2
                        ! cen = (/vf%cfg%xm(i+i_in),vf%cfg%ym(j+j_in),vf%cfg%zm(k)/)
                        cen = calculateCentroid(this%vf%interface_polygon(1,i+i_in,j+j_in,k))
                        ! Get Plane
                        plane = this%vf%liquid_gas_interface(i+i_in,j+j_in,k)
                        ! if((vf%VF(i+i_in,j+j_in,k) .gt. VFlo .and. vf%VF(i+i_in,j+j_in,k) .lt. VFhi) .and. (i .eq. 16) .and. (j .eq. 24)) then
                        !    write(*,'(A,3I10.5)') 'Mid Cell: ', i,j,k
                        !    write(*,'(A,3I10.5)') 'Curr Cell: ', i+i_in,j+j_in,k 
                        !    write(*,'(A,3F10.5)') 'Cen: ', cen 
                        !    call printToScreen(plane)
                        call addMember(neighborhood,cen,1.0_WP,plane)
                    endif
                    end do
                end do
                ! Set Neighborhood in solver
                call setNeighborhood(solver,neighborhood)
                call setThreshold(solver,0.1875_WP) ! This is the of the wendland function at 0.5 (is radius is 1).  
                ! ====== Get Stresses

                ! sigma_xx, we go on right of u cell
                ! The right of the u cell goes y(j) to y(j+1) at xm(i)
                startPoint = (/this%fs%cfg%xm(i),this%fs%cfg%y(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%xm(i),this%fs%cfg%y(j+1),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the x component of this force, so take it and store it as sigma_xx
                this%sigma_xx(i,j,k) = force(1)
                ! ====== For Debugging
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_xx_NoP(i,j,k) = force(1)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_xx_P(i,j,k) = force(1)
                this%sigma_xx_P(i,j,k) = this%sigma_xx_P(i,j,k) - this%sigma_xx_NoP(i,j,k)


                ! sigma_xy, we go on bottom of u cell
                ! The bottom of the u cell goes xm(i-1) to xm(i) at y(j)
                startPoint = (/this%fs%cfg%xm(i-1),this%fs%cfg%y(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%xm(i),this%fs%cfg%y(j),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the x component of this force, so take it and store it as sigma_xy 
                this%sigma_xy(i,j,k) = force(1)

                ! ====== For Debugging
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_xy_NoP(i,j,k) = force(1)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_xy_P(i,j,k) = force(1)
                this%sigma_xy_P(i,j,k) = this%sigma_xy_P(i,j,k) - this%sigma_xy_NoP(i,j,k)

                ! sigma_yx, we go on left of v cell
                ! The bottom of the u cell goes ym(j-1) to ym(j) at x(i)
                startPoint = (/this%fs%cfg%x(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%x(i),this%fs%cfg%ym(j-1),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the y component of this force, so take it and store it as sigma_yx 
                this%sigma_yx(i,j,k) = force(2)
                
                ! ====== For Debugging
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_yx_NoP(i,j,k) = force(2)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_yx_P(i,j,k) = force(2)
                this%sigma_yx_P(i,j,k) = this%sigma_yx_P(i,j,k) - this%sigma_yx_NoP(i,j,k)

                ! sigma_yy, we go on top of v cell
                ! The bottom of the u cell goes x(i) to x(i+1) at ym(j) 
                startPoint = (/this%fs%cfg%x(i+1),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%x(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the y component of this force, so take it and store it as sigma_yy
                this%sigma_yy(i,j,k) = force(2)

                ! ====== For Debugging  
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_yy_NoP(i,j,k) = force(2)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_yy_P(i,j,k) = force(2)
                this%sigma_yy_P(i,j,k) = this%sigma_yy_P(i,j,k) - this%sigma_yy_NoP(i,j,k)
            ! endif
        end do
        end do
    end do
    ! Boundary Conditions
    call this%vf%cfg%sync(this%sigma_xx)
    call this%vf%cfg%sync(this%sigma_xy)
    call this%vf%cfg%sync(this%sigma_yx)
    call this%vf%cfg%sync(this%sigma_yy)

    call this%vf%cfg%sync(this%sigma_xx_NoP)
    call this%vf%cfg%sync(this%sigma_xy_NoP)
    call this%vf%cfg%sync(this%sigma_yx_NoP)
    call this%vf%cfg%sync(this%sigma_yy_NoP)

    call this%vf%cfg%sync(this%sigma_xx_P)
    call this%vf%cfg%sync(this%sigma_xy_P)
    call this%vf%cfg%sync(this%sigma_yx_P)
    call this%vf%cfg%sync(this%sigma_yy_P)
end subroutine updateSurfaceTensionStresses

subroutine updateSurfaceTensionForces(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUSolve_RectCub_class
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k 

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! write(*,'(A,F10.5)') 'DX,i,j,k: ', cfg%dxmi(i)
            this%Fst_x(i,j,k) = (this%sigma_xx(i,j,k)-this%sigma_xx(i-1,j,k)) - (this%sigma_xy(i,j+1,k)-this%sigma_xy(i,j,k))
            this%Fst_x(i,j,k) = this%Fst_x(i,j,k) * this%fs%cfg%dxi(i)

            this%Fst_y(i,j,k) = (this%sigma_yy(i,j,k)-this%sigma_yy(i,j-1,k)) - (this%sigma_yx(i+1,j,k)-this%sigma_yx(i,j,k))
            this%Fst_y(i,j,k) = this%Fst_y(i,j,k) * this%fs%cfg%dyi(j)
            this%Fst_z(i,j,k) = 0.0_WP
        end do
        end do
    end do
    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z

    this%Pjx_ST = this%Fst_x
    this%Pjy_ST = this%Fst_y
    this%Pjz_ST = this%Fst_z
end subroutine updateSurfaceTensionForces

! ####################################################### SMOOTHING FUNCTIONS ################################################
subroutine applyLaplacianSmoothing(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUSolve_RectCub_class

    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k 
    
    this%Fst_x = 0.0_WP;
    this%Fst_y = 0.0_WP;
    this%Fst_z = 0.0_WP;

    ! Smoothing
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%Fst_x(i,j,k) = 0.125*(this%fs%Pjx(i+1,j,k)+this%fs%Pjx(i-1,j,k) + this%fs%Pjx(i,j+1,k) + this%fs%Pjx(i,j-1,k)) + 0.5*this%fs%Pjx(i,j,k)
            this%Fst_y(i,j,k) = 0.125*(this%fs%Pjy(i+1,j,k)+this%fs%Pjy(i-1,j,k) + this%fs%Pjy(i,j+1,k) + this%fs%Pjy(i,j-1,k)) + 0.5*this%fs%Pjy(i,j,k)
        end do
        end do
    end do

    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z
end subroutine applyLaplacianSmoothing

subroutine applyGradientSmoothing(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUSolve_RectCub_class

    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,ii,jj,kk,negativeIndex,positiveIndex
    real(WP) :: tot_grad_x,tot_grad_y,tot_force

    this%Fst_x = 0.0_WP;
    this%Fst_y = 0.0_WP;
    this%Fst_z = 0.0_WP;
    
    negativeIndex = -1
    positiveIndex = 1
    ! Calculate Gradients
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%grad_vf_x(i,j,k) = abs((this%vf%VF(i  ,j  ,k) - this%vf%VF(i-1,j  ,k))  * this%vf%cfg%dxi(i))
            this%grad_vf_y(i,j,k) = abs((this%vf%VF(i  ,j  ,k) - this%vf%VF(i  ,j-1,k)) * this%vf%cfg%dyi(j))
            ! While we are looping everything, also check to see if everything has been processed
            if(abs(this%fs%Pjx(i,j,k)) .gt. 1e-10) then
                this%x_smoothing_tracker(i,j,k) = 1 ! Mark as needing smoothing
            endif

            if(abs(this%fs%Pjy(i,j,k)) .gt. 1e-10) then
                this%y_smoothing_tracker(i,j,k) = 1 ! Mark as needing smoothing
            endif

        end do
        end do
    end do
    ! Smoothing
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            tot_grad_x = 0.0_WP
            tot_grad_y = 0.0_WP
            ! X Direction
            if(abs(this%fs%Pjx(i,j,k)) .gt. 1e-12 .and. this%x_smoothing_tracker(i,j,k) == 1) then  ! This checks to see if we are in a mixed cell that has not yet been smoothed
                ! Current Cell
                tot_grad_x = this%grad_vf_x(i,j,k)
                tot_force = this%fs%Pjx(i,j,k)
                ! Going Left
                do while(abs(this%fs%Pjx(i+negativeIndex,j,k)) .gt. 1e-12 .and. (i+negativeIndex) .ge. this%fs%cfg%imin_)
                    tot_grad_x = tot_grad_x + this%grad_vf_x(i+negativeIndex,j,k)
                    tot_force = tot_force + this%fs%Pjx(i+negativeIndex,j,k)
                    negativeIndex = negativeIndex -1
                end do
                ! Going Right
                do while(abs(this%fs%Pjx(i+positiveIndex,j,k)) .gt. 1e-12 .and. (i+positiveIndex) .le. this%fs%cfg%imax_)
                    tot_grad_x = tot_grad_x + this%grad_vf_x(i+positiveIndex,j,k)
                    tot_force = tot_force + this%fs%Pjx(i+positiveIndex,j,k)
                    positiveIndex = positiveIndex +1
                end do
                ! Now Reassign Force Values in between 
                do ii = negativeIndex,positiveIndex
                    this%Fst_x(i+ii,j,k) = this%grad_vf_x(i+ii,j,k)*tot_force/(tot_grad_x + 1e-12_WP) ! Summing Components
                    
                    ! Mark as Smoothed
                    this%x_smoothing_tracker(i+ii,j,k) = 0
                end do
                ! Reset Indicies
                negativeIndex = -1
                positiveIndex = 1
                tot_force = 0.0_WP
            endif
            
            ! Y Direction
            if(abs(this%fs%Pjy(i,j,k)) .gt. 1e-12 .and. this%y_smoothing_tracker(i,j,k) == 1) then  
                ! Current Cell
                tot_grad_y = this%grad_vf_y(i,j,k)
                tot_force = this%fs%Pjy(i,j,k)
                ! Going Down
                do while(abs(this%fs%Pjy(i,j+negativeIndex,k)) .gt. 1e-12 .and. (j+negativeIndex) .ge. this%fs%cfg%jmin_)
                    tot_grad_y = tot_grad_y + this%grad_vf_y(i,j+negativeIndex,k)
                    tot_force = tot_force + this%fs%Pjy(i,j+negativeIndex,k)
                    negativeIndex = negativeIndex -1
                end do
                ! Going Up
                do while(abs(this%fs%Pjy(i,j+positiveIndex,k)) .gt. 1e-12 .and. (j+positiveIndex) .le. this%fs%cfg%jmax_)
                    tot_grad_y = tot_grad_y + this%grad_vf_y(i,j+positiveIndex,k)
                    tot_force = tot_force + this%fs%Pjy(i,j+positiveIndex,k)
                    positiveIndex = positiveIndex +1
                end do
                ! Now Reassign Force Values in between 
                do jj = negativeIndex,positiveIndex
                    this%Fst_y(i,j+jj,k) = this%grad_vf_y(i,j+jj,k)*tot_force/(tot_grad_y + 1e-12_WP) ! Summing Components
                    ! Mark as Smoothed
                    this%y_smoothing_tracker(i,j+jj,k) = 0
                end do
                ! Reset Indicies
                negativeIndex = -1
                positiveIndex = 1
                tot_force = 0.0_WP
            endif 

        end do 
        end do
    end do
    ! write(*,'(A)') 'Gradient Smoothing Complete'
    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z
    
end subroutine applyGradientSmoothing

subroutine applyPoissonSmoothing(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUSolve_RectCub_class
    use hypre_str_class, only: pcg_smg,gmres

    class(conservative_st_type), intent(inout) :: this
    type(hypre_str) :: poisson_solver ! Poisson Solver Object

    integer :: i,j,k,ii,jj,kk
    real(WP) :: tot_grad_x,tot_grad_y


    ! Initialize Poisson Solver
    poisson_solver=hypre_str(cfg=cfg,name='Smoothing',method=gmres,nst=7)
    poisson_solver%maxit=this%fs%psolv%maxit
    poisson_solver%rcvg=this%fs%psolv%rcvg
    ! Setup Solver
    ! Set 7-pt stencil map for the pressure solver
    poisson_solver%stc(1,:)=[ 0, 0, 0]
    poisson_solver%stc(2,:)=[+1, 0, 0]
    poisson_solver%stc(3,:)=[-1, 0, 0]
    poisson_solver%stc(4,:)=[ 0,+1, 0]
    poisson_solver%stc(5,:)=[ 0,-1, 0]
    poisson_solver%stc(6,:)=[ 0, 0,+1]
    poisson_solver%stc(7,:)=[ 0, 0,-1]
    ! Set Up Laplacian Operator
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! Set Laplacian
            poisson_solver%opr(1,i,j,k)=this%fs%divp_x(1,i,j,k)*this%fs%divu_x(-1,i+1,j,k)+&
            &                       this%fs%divp_x(0,i,j,k)*this%fs%divu_x( 0,i  ,j,k)+&
            &                       this%fs%divp_y(1,i,j,k)*this%fs%divv_y(-1,i,j+1,k)+&
            &                       this%fs%divp_y(0,i,j,k)*this%fs%divv_y( 0,i,j  ,k)+&
            &                       this%fs%divp_z(1,i,j,k)*this%fs%divw_z(-1,i,j,k+1)+&
            &                       this%fs%divp_z(0,i,j,k)*this%fs%divw_z( 0,i,j,k  )
            poisson_solver%opr(2,i,j,k)=this%fs%divp_x(1,i,j,k)*this%fs%divu_x( 0,i+1,j,k)
            poisson_solver%opr(3,i,j,k)=this%fs%divp_x(0,i,j,k)*this%fs%divu_x(-1,i  ,j,k)
            poisson_solver%opr(4,i,j,k)=this%fs%divp_y(1,i,j,k)*this%fs%divv_y( 0,i,j+1,k)
            poisson_solver%opr(5,i,j,k)=this%fs%divp_y(0,i,j,k)*this%fs%divv_y(-1,i,j  ,k)
            poisson_solver%opr(6,i,j,k)=this%fs%divp_z(1,i,j,k)*this%fs%divw_z( 0,i,j,k+1)
            poisson_solver%opr(7,i,j,k)=this%fs%divp_z(0,i,j,k)*this%fs%divw_z(-1,i,j,k  )
            ! Scale it by the cell volume
            ! poisson_solver%opr(:,i,j,k)=-poisson_solver%opr(:,i,j,k)*this%cfg%vol(i,j,k) I don't think I need this right now.
        end do
        end do
    end do
    ! Initialize the Poisson solver
    call poisson_solver%init()
    call poisson_solver%setup()
    poisson_solver%sol=this%fs%Pjx
    ! We now need to calculate the right hand side, which will be the divergence of the force field
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = (sum(this%fs%divp_x(:,i,j,k)*this%fs%Pjx(i:i+1,j,k))&
                                        +sum(this%fs%divp_y(:,i,j,k)*this%fs%Pjy(i,j:j+1,k))) ! No z direction right now
        end do
        end do
    end do
    
    poisson_solver%rhs = this%SurfaceTensionDiv
    this%poisson_source = this%SurfaceTensionDiv
    ! Solve the Poisson Equation
    call poisson_solver%solve()
    this%force_potential_field = poisson_solver%sol
    ! Now we need to calculate the gradient of the smoothed pressure field to get the new forces
    ! Calculate Forces
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%Fst_x(i,j,k) = (poisson_solver%sol(i,j,k) - poisson_solver%sol(i-1,j,k)) * this%fs%cfg%dxi(i)
            this%Fst_y(i,j,k) = (poisson_solver%sol(i,j,k) - poisson_solver%sol(i,j-1,k)) * this%fs%cfg%dyi(j)
            this%Fst_z(i,j,k) = 0.0_WP
        end do
        end do
    end do   
    ! Calculate Magnitude of Fst_x
    ! write(*,'(A,F10.5)') 'Max Fst_x before assign: ', maxval(abs(Fst_x))
    ! Assign
    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z
    call poisson_solver%destroy()
    ! write(*,'(A)') 'Poisson Smoothing Complete'
end subroutine applyPoissonSmoothing
end module conservative_st