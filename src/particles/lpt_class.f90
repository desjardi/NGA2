!> Basic Lagrangian particle solver class:
!> Provides support for Lagrangian-transported objects
module lpt_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: lpt
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   
   !> Particle object definition
   type :: part
      !> Particle center coordinates
      real(WP) :: x,y,z
      !> Index of cell containing particle center
      integer :: i,j,k
      ! Velocity of particle
      real(WP) :: u,v,w
      ! Time step size
      real(WP) :: dt
      ! Control parameter
      integer :: stop
   end type part
   
   
   !> Lagrangian particle tracking solver object definition
   type :: lpt
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_LPT'     !< Solver name (default=UNNAMED_LPT)
      
      ! Particle data
      integer :: np                                       !< Global number of particles
      integer :: np_                                      !< Local number of particles
      integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
      type(part), dimension(:), allocatable :: p          !< Array of particles of type part
      
      ! Particle density
      real(WP) :: rho                                    !< Density of particle
      
      ! Solver parameters
      real(WP) :: nstep=1                                 !< Number of substeps (default=1)
      logical  :: twoway=.false.                          !< Two-way coupling   (default=no)
      
      ! Monitoring info
      real(WP) :: Umin,Umax,Umean                         !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean                         !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean                         !< W velocity info
      integer  :: np_new,np_out                           !< Number of new and removed particles
      
   contains
      procedure :: print=>lpt_print                       !< Output solver info to the screen
   end type lpt
   
   
   !> Declare lpt solver constructor
   interface lpt
      procedure constructor
   end interface lpt
   
contains
   
   
   !> Default constructor for lpt solver
   function constructor(cfg,name) result(self)
      use messager, only: die
      implicit none
      type(lpt) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Allocate variables
      allocate(self%np_proc(self%cfg%nproc))
      
      ! Initialize MPI derived datatype for a particle
      call self%resize(2)
      prep_mpi_part: block
         integer, dimension(22) :: types,lengths
         integer(kind=MPI_ADDRESS_KIND) :: lower_bound, extent
         integer(kind=MPI_ADDRESS_KIND), dimension(22) :: displacement
         integer :: base,ierr
         
         ! Create the MPI structure to send particles
         types( 1: 1) = MPI_INTEGER8
         types( 2:18) = MPI_REAL_WP
         types(19:22) = MPI_INTEGER
         lengths(:)   = 1
         
         ! Hard-code displacement here
         displacement( 1)=  0
         displacement( 2)=  8
         displacement( 3)= 16
         displacement( 4)= 24
         displacement( 5)= 32
         displacement( 6)= 40
         displacement( 7)= 48
         displacement( 8)= 56
         displacement( 9)= 64
         displacement(10)= 72
         displacement(11)= 80
         displacement(12)= 88
         displacement(13)= 96
         displacement(14)=104
         displacement(15)=112
         displacement(16)=120
         displacement(17)=128
         displacement(18)=136
         displacement(19)=144
         displacement(20)=148
         displacement(21)=152
         displacement(22)=156
         
         ! Finalize by creating and commiting the new type
         call MPI_Type_create_struct(22,lengths,displacement,types,MPI_PARTICLE_TMP,ierr)
         call MPI_Type_get_extent(MPI_PARTICLE_TMP, lower_bound, extent, ierr)
         call MPI_Type_create_resized(MPI_PARTICLE_TMP,lower_bound,extent,MPI_PARTICLE,ierr)
         call MPI_Type_commit(MPI_PARTICLE,ierr)
         
         ! If problem, say it
         if (ierr.ne.0) call die("Problem with MPI_PARTICLE")
         
         ! Get the size of this type
         call MPI_type_size(MPI_PARTICLE,SIZE_MPI_PARTICLE,ierr)
         
      end block prep_mpi_part
      call self%resize(0)
      
      !
      
   end function constructor
   
   
   !> Adaptation of particle array size
   subroutine resize(this,size)
      implicit none
      class(lpt), intent(inout) :: this
      integer, intent(in) :: size
      type(part), dimension(:), allocatable :: tmp
      integer :: n_now,n_new,i
      real(WP) :: up,down
      ! Resize part array to size n
      if (.not.allocated(this%p)) then
         ! part is of size 0
         if (size.gt.0) then
            ! Allocate directly to size n
            allocate(this%p(size))
            this%p(1:n)%stop=1
         end if
      else if (size.eq.0) then
         ! part is associated, but we want to empty it
         deallocate(part)
         nullify(part)
      else
         ! Update non zero size to another non zero size
         n_now = size(part,dim=1)
         up  =real(n_now,WP)*coeff_up
         down=real(n_now,WP)*coeff_down
         if (n.gt.n_now) then
            ! Increase from n_now to n_new
            n_new = max(n,int(up))
            allocate(part_temp(n_new))
            do i=1,n_now
               part_temp(i) = part(i)
            end do
            deallocate(part)
            nullify(part)
            part => part_temp
            part(n_now+1:n_new)%stop=1
         else if (n.lt.int(down)) then
            ! Decrease from n_now to n_new
            allocate(part_temp(n))
            do i=1,n
               part_temp(i) = part(i)
            end do
            deallocate(part)
            nullify(part)
            part => part_temp
         end if
      end if
   end subroutine resize
   
   
end module lpt_class
