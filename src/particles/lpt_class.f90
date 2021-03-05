!> Basic Lagrangian particle solver class:
!> Provides support for Lagrangian-transported objects
module lpt_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype
   use parallel,       only: MPI_REAL_WP
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: lpt
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   
   
   !> Basic particle object definition
   type :: part
      integer(kind=8) :: id     !< Particle ID
      real(WP) :: d             !< Particle diameter
      real(WP) :: x,y,z         !< Particle center coordinates
      real(WP) :: u,v,w         !< Velocity of particle
      real(WP) :: dt            !< Time step size for the particle
      integer  :: i,j,k         !< Index of cell containing particle center
      integer  :: flag          !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=3
   integer           , dimension(part_nblock) :: part_lblock=[1,8,4]
   type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_INTEGER8,MPI_REAL_WP,MPI_INTEGER]
   !> MPI_PART derived datatype and size
   type(MPI_Datatype) :: MPI_PART
   integer :: MPI_PART_SIZE
   
   !> Lagrangian particle tracking solver object definition
   type :: lpt
      
      ! This is our underlying config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_LPT'     !< Solver name (default=UNNAMED_LPT)
      
      ! Particle data
      integer :: np                                       !< Global number of particles
      integer :: np_                                      !< Local number of particles
      integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
      type(part), dimension(:), allocatable :: p          !< Array of particles of type part
      
      ! Particle density
      real(WP) :: rho                                     !< Density of particle
      
      ! Solver parameters
      real(WP) :: nstep=1                                 !< Number of substeps (default=1)
      logical  :: twoway=.false.                          !< Two-way coupling   (default=no)
      
      ! Monitoring info
      real(WP) :: dmin,dmax,dmean                         !< Diameter info
      real(WP) :: Umin,Umax,Umean                         !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean                         !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean                         !< W velocity info
      integer  :: np_new,np_out                           !< Number of new and removed particles
      
   contains
      procedure :: print=>lpt_print                       !< Output solver info to the screen
      procedure :: resize                                 !< Resize particle array to given size
      procedure :: recycle                                !< Recycle particle array by removing flagged particles
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
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Allocate variables
      allocate(self%np_proc(self%cfg%nproc))
      
      ! Initialize MPI derived datatype for a particle
      call prepare_mpi_part()
      
      !
      
   end function constructor
   
   
   !> Creation of the MPI datatype for particle
   subroutine prepare_mpi_part()
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      use messager, only: die
      implicit none
      integer(MPI_ADDRESS_KIND), dimension(part_nblock) :: disp
      integer(MPI_ADDRESS_KIND) :: mydisp,lb,extent
      type(MPI_Datatype) :: MPI_PART_TMP
      integer :: i,size,ierr
      ! Prepare the displacement array
      disp(1)=0
      do i=2,part_nblock
         call MPI_type_size(part_tblock(i-1),mydisp,ierr)
         disp(i)=disp(i-1)+mydisp*int(part_lblock,MPI_ADDRESS_KIND)
      end do
      ! Create and commit the new type
      call MPI_Type_create_struct(part_nblock,part_lblock,disp,part_tblock,MPI_PART_TMP,ierr)
      call MPI_Type_get_extent(MPI_PART_TMP,lb,extent,ierr)
      call MPI_Type_create_resized(MPI_PART_TMP,lb,extent,MPI_PART,ierr)
      call MPI_Type_commit(MPI_PART,ierr)
      ! If a problem was encountered, say it
      if (ierr.ne.0) call die('[lpt prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Adaptation of particle array size
   subroutine resize(this,size)
      implicit none
      class(lpt), intent(inout) :: this
      integer, intent(in) :: size
      type(part), dimension(:), allocatable :: tmp
      integer :: size_now,size_new
      ! Resize part array to size n
      if (.not.allocated(this%p)) then
         ! part is of size 0
         if (size.gt.0) then
            ! Allocate directly to size n
            allocate(this%p(size))
            this%p(1:n)%flag=1
         end if
      else if (size.eq.0) then
         ! part is associated, but we want to empty it
         deallocate(this%p)
      else
         ! Update from a non-zero size to another non-zero size
         size_now=size(this%p,dim=1)
         if (size.gt.size_now) then
            size_new=max(size,int(real(size_now,WP)*coeff_up))
            allocate(tmp(size_new))
            tmp(1:size_now)=this%p
            tmp(size_now+1:)%flag=1
            call move_alloc(tmp,this%p)
         else if (size.lt.int(real(size_now,WP)*coeff_dn)) then
            allocate(tmp(size))
            tmp(1:size)=this%p(1:size)
            call move_alloc(tmp,this%p)
         end if
      end if
   end subroutine resize
   
   
   !> Clean-up of particle array by removing flag=1 particles
   subroutine recycle(this)
      use mpi_f08
      implicit none
      class(lpt), intent(inout) :: this
      integer :: new_size,ierr,i
      ! Compact all active particles at the beginning of the array
      new_size=0
      if (allocated(this%p)) then
         do i=1,size(this%p,dim=1)
            if (this%p(i)%flag.ne.1) then
               new_size=new_size+1
               if (i.ne.new_size) then
                  this%p(new_size)=this%p(i)
                  this%p(i)%flag=1
               end if
            end if
         end do
      end if
      ! Resize to new size
      call this%resize(new_size)
      ! Update number of particles
      this%np_=size(this%p,dim=1)
      call MPI_ALLGATHER(this%np_,1,MPI_INTEGER,this%np_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
      this%np=sum(this%np_proc)
   end subroutine recycle
   
   
end module lpt_class