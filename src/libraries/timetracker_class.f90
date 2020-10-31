module timetracker_class
   use precision, only: WP
   use string,    only: str_medium,str_long
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: timetracker
   
   
   !> Smoothing parameter for dt selection
   real(WP), parameter :: alpha = 0.7_WP
   !> Give ourselves 10 minutes to clean up
   real(WP), parameter :: wtsafe= 600.0_WP
   
   
   !> Timetracker object
   type :: timetracker
      character(len=str_medium) :: name='UNNAMED_TIME' !< Name for timetracker
      integer  ::   n,  nmax                           !< Current and max timestep
      real(WP) ::   t,  tmax                           !< Current and max time
      real(WP) ::  dt, dtmax                           !< Current and max timestep size
      real(WP) :: cfl,cflmax                           !< Current and max CFL
      real(WP) ::  wt, wtmax                           !< Current and max wallclock time
      real(WP) :: told,dtold                           !< Old time and timestep size
   contains
      procedure :: increment                           !< Default method for incrementing time
      procedure :: adjust_dt                           !< Default method for adjusting timestep size
      procedure :: done                                !< Default termination check for time integration
   end type timetracker
   
   
   !> Declare timetracker constructor
   interface timetracker
      procedure constructor
   end interface timetracker
   
   
contains
   
   
   !> Constructor for timetracker object
   function constructor(name) result(self)
      implicit none
      type(timetracker) :: self
      character(len=*), optional :: name
      if (present(name)) self%name=trim(adjustl(name))
      self%n  =0;        self%nmax=huge(  self%nmax)
      self%t  =0.0_WP;   self%tmax=huge(  self%tmax)
      self%dt =0.0_WP;  self%dtmax=huge( self%dtmax)
      self%cfl=0.0_WP; self%cflmax=huge(self%cflmax)
      self%wt =0.0_WP;  self%wtmax=huge( self%wtmax)
      self%told=0.0_WP
      self%dtold=0.0_WP
   end function constructor
   
   
   !> Time increment
   subroutine increment(this)
      use parallel, only: wtinit,parallel_time
      implicit none
      class(timetracker), intent(inout) :: this
      this%told=this%t
      this%t=this%t+this%dt
      this%n=this%n+1
      this%wt=parallel_time()-wtinit
   end subroutine increment
   
   
   !> Adjust dt based on CFL
   subroutine adjust_dt(this)
      implicit none
      class(timetracker), intent(inout) :: this
      ! Store old timestep size
      this%dtold=this%dt
      ! Estimate new timestep size based on cfl info
      if (this%cfl.gt.0.0_WP) this%dt=min(this%dtold*this%cflmax/this%cfl,this%dtmax)
      ! Prevent dt from increasing too fast
      if (this%dt.gt.this%dtold) this%dt=alpha*this%dt+(1.0_WP-alpha)*this%dtold
   end subroutine adjust_dt
   
   
   !> Termination check
   logical function done(this)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(timetracker), intent(inout) :: this
      character(len=str_long) :: message
      ! Assume not done
      done=.false.
      ! Go through standard ending tests
      if (this%n.ge.this%nmax) then
         done=.true.
         write(message,'(" Timetracker [",a,"] reached maximum number of iterations ")') trim(this%name); call log(message)
      end if
      if (this%t.ge.this%tmax) then
         done=.true.
         write(message,'(" Timetracker [",a,"] reached maximum integration time ")') trim(this%name); call log(message)
      end if
      if ((this%wt+wtsafe).ge.this%wtmax) then
         done=.true.
         write(message,'(" Timetracker [",a,"] reached maximum wallclock time allowed ")') trim(this%name); call log(message)
      end if
   end function done
   
   
end module timetracker_class
