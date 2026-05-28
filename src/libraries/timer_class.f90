module timer_class
   use precision, only: WP
   use string,    only: str_medium
   use mpi_f08,   only: MPI_Comm,MPI_BARRIER,MPI_Wtime
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: timer
   
   !> Timer object
   type :: timer
      character(len=str_medium) :: name='UNNAMED_TIMER'   !< Name for timer
      type(MPI_Comm) :: comm                              !< Communicator for timer
      real(WP) :: time=0.0_WP                             !< Elapsed time (this rank)
      real(WP), private :: stored_time=0.0_WP             !< Stored time
      logical :: is_started=.false.                       !< Is the timer started?
      logical :: sync=.true.                              !< Barrier-synchronize start/stop?
      integer :: ncall=0                                  !< Number of completed start/stop cycles
      ! Statistics (populated by get_stats)
      real(WP) :: tmin=0.0_WP                             !< Minimum time across ranks
      real(WP) :: tmax=0.0_WP                             !< Maximum time across ranks
      real(WP) :: tmean=0.0_WP                            !< Mean time across ranks
      real(WP) :: imbalance_ratio=0.0_WP                  !< Imbalance ratio (tmax/tmin)
      real(WP) :: efficiency=0.0_WP                       !< Load efficiency (tmean/tmax)
      real(WP) :: timepercall=0.0_WP                      !< Per-call average (time/ncall, this rank)
   contains
      procedure :: reset                                  !< Reset timer to zero
      procedure :: start                                  !< Start timer
      procedure :: stop                                   !< Stop timer
      procedure :: get_stats                              !< Compute and store statistics across ranks
      procedure :: finalize                               !< Finalize timer
   end type timer
   
   !> Declare timer constructor
   interface timer
      procedure constructor
   end interface timer   
   
contains
   
   !> Constructor for timer object
   function constructor(comm,name,sync) result(self)
      implicit none
      type(timer) :: self
      type(MPI_Comm), intent(in) :: comm
      character(len=*), optional :: name
      logical, intent(in), optional :: sync
      ! Store communicator
      self%comm=comm
      ! Set timer name
      if (present(name)) self%name=trim(adjustl(name))
      ! Set sync mode (default: .true. for backward compatibility)
      self%sync=.true.
      if (present(sync)) self%sync=sync
      ! Initialize timing
      self%time=0.0_WP
      self%stored_time=0.0_WP
      self%is_started=.false.
      self%ncall=0
      ! Initialize statistics
      self%tmin=0.0_WP; self%tmax=0.0_WP; self%tmean=0.0_WP
      self%imbalance_ratio=0.0_WP; self%efficiency=0.0_WP; self%timepercall=0.0_WP
   end function constructor
   
   !> Reset timer
   subroutine reset(this)
      implicit none
      class(timer), intent(inout) :: this
      this%time=0.0_WP
      this%stored_time=0.0_WP
      this%is_started=.false.
      this%ncall=0
   end subroutine reset
   
   !> Start timer
   subroutine start(this)
      use messager, only: die
      implicit none
      class(timer), intent(inout) :: this
      ! Can only start a stopped timer
      if (this%is_started) call die('[timer start] Timer '//trim(this%name)//' is already started')
      ! Synchronize all processes if requested
      if (this%sync) call MPI_BARRIER(this%comm)
      ! Store current time
      this%stored_time=MPI_Wtime()
      ! Set to started
      this%is_started=.true.
   end subroutine start
   
   !> Stop timer
   subroutine stop(this)
      use messager, only: die
      implicit none
      class(timer), intent(inout) :: this
      ! Can only stop a started timer
      if (.not.this%is_started) call die('[timer stop] Timer '//trim(this%name)//' is already stopped')
      ! Synchronize all processes if requested
      if (this%sync) call MPI_BARRIER(this%comm)
      ! Increment elapsed time by current time minus stored_time
      this%time=this%time+(MPI_Wtime()-this%stored_time)
      ! Set to stopped
      this%is_started=.false.
      ! Increment call counter
      this%ncall=this%ncall+1
   end subroutine stop
   
   !> Compute and store timing statistics across all ranks
   !> For synchronized timers, tmin≈tmax≈tmean and imbalance≈1
   !> For unsynchronized timers, these statistics reveal load imbalance
   subroutine get_stats(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(timer), intent(inout) :: this
      integer :: nproc,ierr
      ! Get communicator size
      call MPI_Comm_size(this%comm,nproc,ierr)
      ! Reduce min
      this%tmin=this%time
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%tmin,1,MPI_REAL_WP,MPI_MIN,this%comm,ierr)
      ! Reduce max
      this%tmax=this%time
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%tmax,1,MPI_REAL_WP,MPI_MAX,this%comm,ierr)
      ! Reduce sum for mean
      this%tmean=this%time
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%tmean,1,MPI_REAL_WP,MPI_SUM,this%comm,ierr)
      this%tmean=this%tmean/real(nproc,WP)
      ! Imbalance ratio
      if (this%tmin.gt.0.0_WP) then
         this%imbalance_ratio=this%tmax/this%tmin
      else
         this%imbalance_ratio=0.0_WP
      end if
      ! Load efficiency
      if (this%tmax.gt.0.0_WP) then
         this%efficiency=this%tmean/this%tmax
      else
         this%efficiency=0.0_WP
      end if
      ! Per-call average (local)
      if (this%ncall.gt.0) then
         this%timepercall=this%time/real(this%ncall,WP)
      else
         this%timepercall=0.0_WP
      end if
   end subroutine get_stats
   
   !> Finalize timer
   subroutine finalize(this)
      use mpi_f08, only: MPI_COMM_NULL
      implicit none
      class(timer), intent(inout) :: this
      call this%reset()
      this%name='UNNAMED_TIMER'
      this%comm=MPI_COMM_NULL
   end subroutine finalize
   
end module timer_class
