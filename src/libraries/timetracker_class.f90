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
   !> Give ourselves a safety margin to avoid missing end of simulation due to round-off
   real(WP), parameter :: dtsafe= 1.0e-6_WP
   
   
   !> Timetracker object
   type :: timetracker
      logical :: amRoot                                !< Timetracker needs to know who's the boss
      character(len=str_medium) :: name='UNNAMED_TIME' !< Name for timetracker
      integer  ::  it, itmax                           !< Current and max sub-iteration
      integer  ::   n,  nmax                           !< Current and max timestep
      real(WP) ::   t,  tmax                           !< Current and max time
      real(WP) ::  dt, dtmax                           !< Current and max timestep size
      real(WP) :: cfl,cflmax                           !< Current and max CFL
      real(WP) ::  wt, wtmax                           !< Current and max wallclock time
      real(WP) :: told,dtold,tmid,dtmid                !< Old/mid time and timestep size
      logical  :: print_info=.true.                    !< Should I print time information?
   contains
      procedure :: increment                           !< Default method for incrementing time
      procedure :: adjust_dt                           !< Default method for adjusting timestep size
      procedure :: done                                !< Default termination check for time integration
      procedure :: print=>timetracker_print            !< Output timetracker info to screen
      procedure :: log  =>timetracker_log              !< Output timetracker info to log
      procedure :: reset                               !< Reset timetracker to zero
   end type timetracker
   
   
   !> Declare timetracker constructor
   interface timetracker
      procedure constructor
   end interface timetracker
   
   
contains
   
   
   !> Constructor for timetracker object
   function constructor(amRoot,name,print_info) result(self)
      implicit none
      type(timetracker) :: self
      logical, intent(in) :: amRoot
      character(len=*), optional :: name
      logical, optional :: print_info
      self%amRoot=amRoot
      if (present(name))       self%name=trim(adjustl(name))
      if (present(print_info)) self%print_info=print_info
      self%n    =0;            self%nmax  =huge(self%nmax)
      self%t    =0.0_WP;       self%tmax  =huge(self%tmax)
      self%dt   =0.0_WP;       self%dtmax =huge(self%dtmax)
      self%cfl  =0.0_WP;       self%cflmax=huge(self%cflmax)
      self%wt   =0.0_WP;       self%wtmax =huge(self%wtmax)
      self%told =0.0_WP
      self%dtold=0.0_WP
      self%tmid =0.0_WP
      self%dtmid=0.0_WP
      self%it   =1;            self%itmax=1
   end function constructor
   
   
   !> Time increment
   subroutine increment(this)
      use param,    only: verbose
      use parallel, only: wtinit,parallel_time
      implicit none
      class(timetracker), intent(inout) :: this
      ! Move to next time step
      this%it=1
      this%told=this%t
      this%t=this%t+this%dt
      this%tmid=0.5_WP*(this%t+this%told)
      this%dtmid=0.5_WP*(this%dt+this%dtold)
      this%n=this%n+1
      this%wt=parallel_time()-wtinit
      ! If verbose run, log and or print info
      if (this%print_info) then
         if (verbose.gt.0) call this%log
         if (verbose.gt.1) call this%print
      end if
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
      if (this%t+dtsafe*this%dt.ge.this%tmax) then
         done=.true.
         write(message,'(" Timetracker [",a,"] reached maximum integration time ")') trim(this%name); call log(message)
      end if
      if ((this%wt+wtsafe).ge.this%wtmax) then
         done=.true.
         write(message,'(" Timetracker [",a,"] reached maximum wallclock time allowed ")') trim(this%name); call log(message)
      end if
   end function done
   
   
   !> Reset time to zero
   subroutine reset(this)
      implicit none
      class(timetracker), intent(inout) :: this
      this%n    =0
      this%t    =0.0_WP
      this%wt   =0.0_WP
      this%told =0.0_WP
      this%dtold=0.0_WP
      this%tmid =0.0_WP
      this%it=1
   end subroutine reset


   !> Log timetracker info
   subroutine timetracker_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(timetracker), intent(in) :: this
      character(len=str_long) :: message
      if (this%amRoot) then
         write(message,'("Timetracker [",a,"] status")') trim(this%name); call log(message)
         write(message,'(" >  it/ itmax = ",i0,"/",i0)')         this%it,this%itmax; call log(message)
         write(message,'(" >   n/  nmax = ",i0,"/",i0)')         this%n,this%nmax; call log(message)
         write(message,'(" >   t/  tmax = ",es12.5,"/",es12.5)') this%t,this%tmax; call log(message)
         write(message,'(" >  dt/ dtmax = ",es12.5,"/",es12.5)') this%dt,this%dtmax; call log(message)
         write(message,'(" > cfl/cflmax = ",es12.5,"/",es12.5)') this%cfl,this%cflmax; call log(message)
         write(message,'(" >  wt/ wtmax = ",es12.5,"/",es12.5)') this%wt,this%wtmax; call log(message)
      end if
   end subroutine timetracker_log
   
   
   !> Print out timetracker info
   subroutine timetracker_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(timetracker), intent(in) :: this
      if (this%amRoot) then
         write(output_unit,'("Timetracker [",a,"] status")') trim(this%name)
         write(output_unit,'(" >  it/ itmax = ",i0,"/",i0)')         this%it,this%itmax
         write(output_unit,'(" >   n/  nmax = ",i0,"/",i0)')         this%n,this%nmax
         write(output_unit,'(" >   t/  tmax = ",es12.5,"/",es12.5)') this%t,this%tmax
         write(output_unit,'(" >  dt/ dtmax = ",es12.5,"/",es12.5)') this%dt,this%dtmax
         write(output_unit,'(" > cfl/cflmax = ",es12.5,"/",es12.5)') this%cfl,this%cflmax
         write(output_unit,'(" >  wt/ wtmax = ",es12.5,"/",es12.5)') this%wt,this%wtmax
      end if
   end subroutine timetracker_print
   
   
end module timetracker_class
