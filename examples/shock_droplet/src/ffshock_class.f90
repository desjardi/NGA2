!> Definition of a farfield shock problem
module ffshock_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use spcomp_class,      only: spcomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: ffshock
   
   !> FFshock object
   type :: ffshock
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(spcomp) :: fs        !< Single-phase compressible solver
      type(timetracker) :: time !< Time info
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out
      
      !> Simulation monitor file
      type(monitor) :: mfile,cflfile,consfile
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
      real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc
      
      !> Constant phasic kinematic viscosity
      real(WP) :: cst_visc
      
   contains
      procedure :: initialize                      !< Initialize farfield shock simulation
      procedure :: finalize                        !< Finalize farfield shock simulation
      procedure :: step                            !< Advance farfield shock simulation by one time step
      procedure :: output_monitor                  !< Monitoring for farfield shock case
      procedure :: output_ensight                  !< Ensight output for farfield shock case
      procedure, private :: prepare_viscosities    !< Prepare viscosities
      procedure, private :: apply_bconds           !< Apply boundary conditions
   end type ffshock
   
contains
   
   
   !> Initialization of a far-field shock problem
   subroutine initialize(this,dx,meshsize,startloc,group,partition)
      use mpi_f08, only: MPI_Group
      implicit none
      class(ffshock), intent(inout) :: this
      real(WP), intent(in) :: dx
      integer , dimension(3), intent(in) :: meshsize
      real(WP), dimension(3), intent(in) :: startloc
      type(MPI_Group)       , intent(in) :: group
      integer , dimension(3), intent(in) :: partition
      
      ! Initialize config object
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         type(sgrid) :: grid
         integer  :: i,j,k
         real(WP) :: Lx,Ly,Lz
         logical  :: xper,yper,zper
         real(WP), dimension(:), allocatable :: x,y,z
         ! By default, the case is fully non-periodic
         xper=.false.; yper=.false.; zper=.false.
         ! Handle 2D cases
         if (meshsize(1).eq.1) xper=.true.
         if (meshsize(2).eq.1) yper=.true.
         if (meshsize(3).eq.1) zper=.true.
         ! Create simple rectilinear grid
         allocate(x(meshsize(1)+1)); do i=1,meshsize(1)+1; x(i)=real(i-1,WP)*dx+startloc(1); end do
         allocate(y(meshsize(2)+1)); do j=1,meshsize(2)+1; y(j)=real(j-1,WP)*dx+startloc(2); end do
         allocate(z(meshsize(3)+1)); do k=1,meshsize(3)+1; z(k)=real(k-1,WP)*dx+startloc(3); end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=xper,yper=yper,zper=zper,name='FarField')
         ! Deallocate x/y/z
         deallocate(x,y,z)
         ! Create config
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      ! Initialize time tracker
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
      end block initialize_timetracker
      
      ! Create multiphase compressible flow solver
      create_velocity_solver: block
         call this%fs%initialize(cfg=this%cfg,name='Compressible NS')
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%dQdt(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%fs%nQ,1:4))
         allocate(this%beta(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ma(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='FarField')
         ! No need to output FVF file
         this%ens_out%write_fvf=.false.
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('RHO',this%fs%Q(:,:,:,1))
         call this%ens_out%add_scalar('I',this%fs%I)
         call this%ens_out%add_scalar('P',this%fs%P)
         call this%ens_out%add_scalar('Mach',this%Ma)
         call this%ens_out%add_scalar('beta',this%beta)
         call this%ens_out%add_scalar('visc',this%visc)
      end block create_ensight
      
      ! Create monitor files
      create_monitor: block
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'farfield_sim')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%RHOmax,'max(RHO)')
         call this%mfile%add_column(this%fs%RHOmin,'min(RHO)')
         call this%mfile%add_column(this%fs%Imax  ,'max(I)'  )
         call this%mfile%add_column(this%fs%Imin  ,'min(I)'  )
         call this%mfile%add_column(this%fs%Pmax  ,'max(P)'  )
         call this%mfile%add_column(this%fs%Pmin  ,'min(P)'  )
         call this%mfile%add_column(this%fs%Tmax  ,'max(T)'  )
         call this%mfile%add_column(this%fs%Tmin  ,'min(T)'  )
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'farfield_cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLa_x,'Acoustic xCFL')
         call this%cflfile%add_column(this%fs%CFLa_y,'Acoustic yCFL')
         call this%cflfile%add_column(this%fs%CFLa_z,'Acoustic zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         ! Create conservation monitor
         this%consfile=monitor(this%fs%cfg%amRoot,'farfield_cons')
         call this%consfile%add_column(this%time%n,'Timestep number')
         call this%consfile%add_column(this%time%t,'Time')
         call this%consfile%add_column(this%fs%Qint(1),'Mass')
         call this%consfile%add_column(this%fs%Qint(2),'Energy')
         call this%consfile%add_column(this%fs%Qint(3),'U Momentum')
         call this%consfile%add_column(this%fs%Qint(4),'V Momentum')
         call this%consfile%add_column(this%fs%Qint(5),'W Momentum')
         call this%consfile%add_column(this%fs%RHOKint,'Kinetic Energy')
         call this%consfile%add_column(this%fs%RHOSint,'Entropy')
      end block create_monitor
      
   end subroutine initialize
   
   
   !> Take one time step
   subroutine step(this,dt)
      implicit none
      class(ffshock), intent(inout) :: this
      real(WP) :: dt
      
      ! Increment time
      this%time%dt=dt
      call this%fs%get_cfl(dt=this%time%dt,cfl=this%time%cfl)
      call this%time%increment()
      
      ! Remember conserved variables
      this%fs%Qold=this%fs%Q
      
      ! Prepare SGS viscosity models
      call this%prepare_viscosities()
      
      ! First RK step ====================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,1))
      this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,1)
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Second RK step ===================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,2))
      this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,2)
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Third RK step ====================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,3))
      this%fs%Q=this%fs%Qold+1.0_WP*this%time%dt*this%dQdt(:,:,:,:,3)
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Fourth RK step ===================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,4))
      this%fs%Q=this%fs%Qold+this%time%dt/6.0_WP*(this%dQdt(:,:,:,:,1)+2.0_WP*this%dQdt(:,:,:,:,2)+2.0_WP*this%dQdt(:,:,:,:,3)+this%dQdt(:,:,:,:,4))
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Apply boundary conditions
      call this%apply_bconds()
      
      ! Interpolate velocity
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      
      ! Compute local Mach number
      this%Ma=sqrt(this%Ui**2+this%Vi**2+this%Wi**2)/this%fs%C
      
   end subroutine step
   
   
   !> Perform and output monitoring for the far-field shock problem
   subroutine output_monitor(this)
      implicit none
      class(ffshock), intent(inout) :: this
      call this%fs%get_info()
      call this%mfile%write()
      call this%cflfile%write()
      call this%consfile%write()
   end subroutine output_monitor
   
   
   !> Output ensight files for the far-field shock problem
   subroutine output_ensight(this,t)
      implicit none
      class(ffshock), intent(inout) :: this
      real(WP), intent(in), optional :: t
      if (present(t)) then
         call this%ens_out%write_data(t)
      else
         call this%ens_out%write_data(this%time%t)
      end if
   end subroutine output_ensight
   
   
   !> Calculate viscosities
   subroutine prepare_viscosities(this)
      implicit none
      class(ffshock), intent(inout) :: this
      ! Get LAD
      call this%fs%get_viscartif(dt=this%time%dt,beta=this%beta); this%fs%BETA=this%fs%Q(:,:,:,1)*(this%beta              )
      ! Get eddy viscosity
      call this%fs%get_vreman   (dt=this%time%dt,visc=this%visc); this%fs%VISC=this%fs%Q(:,:,:,1)*(this%visc+this%cst_visc)
   end subroutine prepare_viscosities
   
   
   !> Apply boundary conditions
   subroutine apply_bconds(this)
      implicit none
      class(ffshock), intent(inout) :: this
      integer :: i,j,k
      
      ! Apply clipped Neumann on primitive variables in x+
      if (.not.this%fs%cfg%xper.and.this%fs%cfg%iproc.eq.this%fs%cfg%npx) then
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
            ! Copy over from imax to imax+1 and above
            do i=this%fs%cfg%imax+1,this%fs%cfg%imaxo
               ! Copy primitive variables
               this%fs%Q(i,j,k,1)=this%fs%Q(this%fs%cfg%imax,j,k,1)
               this%fs%P(i,j,k)=this%fs%P(this%fs%cfg%imax,j,k)
               this%fs%I(i,j,k)=this%fs%I(this%fs%cfg%imax,j,k)
               this%fs%U(i,j,k)=max(this%fs%U(this%fs%cfg%imax,j,k),0.0_WP)
               this%fs%V(i,j,k)=this%fs%V(this%fs%cfg%imax,j,k)
               this%fs%W(i,j,k)=this%fs%W(this%fs%cfg%imax,j,k)
            end do
         end do; end do
      end if
      
      ! Apply clipped Neumann on primitive variables in y+
      if (.not.this%fs%cfg%yper.and.this%fs%cfg%jproc.eq.this%fs%cfg%npy) then
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! Copy over from jmax to jmax+1 and above
            do j=this%fs%cfg%jmax+1,this%fs%cfg%jmaxo
               ! Copy primitive variables
               this%fs%Q(i,j,k,1)=this%fs%Q(i,this%fs%cfg%jmax,k,1)
               this%fs%P(i,j,k)=this%fs%P(i,this%fs%cfg%jmax,k)
               this%fs%I(i,j,k)=this%fs%I(i,this%fs%cfg%jmax,k)
               this%fs%U(i,j,k)=this%fs%U(i,this%fs%cfg%jmax,k)
               this%fs%V(i,j,k)=max(this%fs%V(i,this%fs%cfg%jmax,k),0.0_WP)
               this%fs%W(i,j,k)=this%fs%W(i,this%fs%cfg%jmax,k)
            end do
         end do; end do
      end if
      
      ! Apply clipped Neumann on primitive variables in y-
      if (.not.this%fs%cfg%yper.and.this%fs%cfg%jproc.eq.1) then
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! First copy over V from jmin+1 to jmin
            this%fs%V(i,this%fs%cfg%jmin,k)=min(this%fs%V(i,this%fs%cfg%jmin+1,k),0.0_WP)
            ! Then copy over from jmin to jmin-1 and below
            do j=this%fs%cfg%jmino,this%fs%cfg%jmin-1
               ! Copy primitive variables
               this%fs%Q(i,j,k,1)=this%fs%Q(i,this%fs%cfg%jmin,k,1)
               this%fs%P(i,j,k)=this%fs%P(i,this%fs%cfg%jmin,k)
               this%fs%I(i,j,k)=this%fs%I(i,this%fs%cfg%jmin,k)
               this%fs%U(i,j,k)=this%fs%U(i,this%fs%cfg%jmin,k)
               this%fs%V(i,j,k)=min(this%fs%V(i,this%fs%cfg%jmin,k),0.0_WP)
               this%fs%W(i,j,k)=this%fs%W(i,this%fs%cfg%jmin,k)
            end do
         end do; end do
      end if
      
      ! Apply clipped Neumann on primitive variables in z+
      if (.not.this%fs%cfg%zper.and.this%fs%cfg%kproc.eq.this%fs%cfg%npz) then
         do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! Copy over from kmax to kmax+1 and above
            do k=this%fs%cfg%kmax+1,this%fs%cfg%kmaxo
               ! Copy primitive variables
               this%fs%Q(i,j,k,1)=this%fs%Q(i,j,this%fs%cfg%kmax,1)
               this%fs%P(i,j,k)=this%fs%P(i,j,this%fs%cfg%kmax)
               this%fs%I(i,j,k)=this%fs%I(i,j,this%fs%cfg%kmax)
               this%fs%U(i,j,k)=this%fs%U(i,j,this%fs%cfg%kmax)
               this%fs%V(i,j,k)=this%fs%V(i,j,this%fs%cfg%kmax)
               this%fs%W(i,j,k)=max(this%fs%W(i,j,this%fs%cfg%kmax),0.0_WP)
            end do
         end do; end do
      end if
      
      ! Apply clipped Neumann on primitive variables in z-
      if (.not.this%fs%cfg%zper.and.this%fs%cfg%kproc.eq.1) then
         do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! First copy over W from kmin+1 to kmin
            this%fs%W(i,j,this%fs%cfg%kmin)=min(this%fs%W(i,j,this%fs%cfg%kmin+1),0.0_WP)
            ! Then copy over from kmin to kmin-1 and below
            do k=this%fs%cfg%kmino,this%fs%cfg%kmin-1
               ! Copy primitive variables
               this%fs%Q(i,j,k,1)=this%fs%Q(i,j,this%fs%cfg%kmin,1)
               this%fs%P(i,j,k)=this%fs%P(i,j,this%fs%cfg%kmin)
               this%fs%I(i,j,k)=this%fs%I(i,j,this%fs%cfg%kmin)
               this%fs%U(i,j,k)=this%fs%U(i,j,this%fs%cfg%kmin)
               this%fs%V(i,j,k)=this%fs%V(i,j,this%fs%cfg%kmin)
               this%fs%W(i,j,k)=min(this%fs%W(i,j,this%fs%cfg%kmin),0.0_WP)
            end do
         end do; end do
      end if
      
      ! Rebuild conserved quantities
      this%fs%Q(:,:,:,2)=this%fs%Q(:,:,:,1)*this%fs%I
      call this%fs%get_momentum()
      
   end subroutine apply_bconds
   
   
   !> Finalize the shockdrop problem
   subroutine finalize(this)
      implicit none
      class(ffshock), intent(inout) :: this
      ! Deallocate work arrays
      deallocate(this%dQdt,this%Ui,this%Vi,this%Wi,this%Ma,this%beta,this%visc)
      ! Finalize all objects here
      call this%cfg%finalize()
      call this%fs%finalize()
      call this%time%finalize()
      call this%ens_out%finalize()
      call this%mfile%finalize()
      call this%cflfile%finalize()
      call this%consfile%finalize()
   end subroutine finalize
   
   
end module ffshock_class
