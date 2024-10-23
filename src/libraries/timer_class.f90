module timer_class
   use precision,         only: WP
   use string,            only: str_medium
   use mpi_f08,           only: MPI_Comm,MPI_BARRIER,MPI_Wtime
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: timer
   
   !> Timer object
   type :: timer
      character(len=str_medium) :: name='UNNAMED_TIMER'   !< Name for timer
      type(MPI_Comm) :: comm                              !< Communicator for timer
      real(WP) :: time                                    !< Elapsed time
      real(WP), private :: stored_time                    !< Stored time
      logical :: is_started                               !< Is the timer started?
   contains
      procedure :: reset                                  !< Reset timer to zero
      procedure :: start                                  !< Start timer
      procedure :: stop                                   !< Stop timer
   end type timer
   
   !> Declare timer constructor
   interface timer
      procedure constructor
   end interface timer   
   
contains
   
   !> Constructor for timer object
   function constructor(comm,name) result(self)
      implicit none
      type(timer) :: self
      type(MPI_Comm), intent(in) :: comm
      character(len=*), optional :: name
      ! Store communicator
      self%comm=comm
      ! Set timer name
      if (present(name)) self%name=trim(adjustl(name))
      ! Initialize
      self%time=0.0_WP
      self%stored_time=0.0_WP
      self%is_started=.false.
   end function constructor
   
   !> Reset timer
   subroutine reset(this)
      implicit none
      class(timer), intent(inout) :: this
      this%time=0.0_WP
      this%stored_time=0.0_WP
      this%is_started=.false.
   end subroutine reset
   
   !> Start timer
   subroutine start(this)
      use messager, only: die
      implicit none
      class(timer), intent(inout) :: this
      ! Can only start a stopped timer
      if (this%is_started) call die('[timer start] Timer '//trim(this%name)//' is already started')
      ! Synchronize all processes
      call MPI_BARRIER(this%comm)
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
      ! Synchronize all processes
      call MPI_BARRIER(this%comm)
      ! Increment elapsed time by current time minus stored_time
      this%time=this%time+(MPI_Wtime()-this%stored_time)
      ! Set to stopped
      this%is_started=.false.
   end subroutine stop
   
end module timer_class
