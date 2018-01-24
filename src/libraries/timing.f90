!> NGA timing tool
module timing
  use precision
  use string
  implicit none
  
  ! Number of timer
  integer :: ntimers
  ! Definition of a timer
  type :: timer
     character(len=str_medium) :: name
     real(WP) :: time
     real(WP) :: time_in
     logical  :: started
  end type timer
  ! Array of timers
  type(timer), dimension(:), allocatable :: timers
  ! Timing for full timestep
  real(WP) :: full_time_in
  ! Array of values to transfer to monitor
  real(WP), dimension(:), allocatable :: mval
  
contains
  
  
  !> Find a timer of a given name
  subroutine find_timer(name,itimer)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(out) :: itimer
    loop: do itimer=1,ntimers
       if (trim(timers(itimer)%name).eq.trim(name)) exit loop
    end do loop
    if (itimer.gt.ntimers) then
       call die('[find_timer] Unknown timer:'//trim(name))
    end if
    return
  end subroutine find_timer
  
  
end module timing


!> Initialize the timing module
subroutine timing_init
  use timing
  use parallel
  implicit none
  ! Start with no timers
  ntimers=0
  ! Get time for full time step
  call parallel_time(full_time_in)
  return
end subroutine timing_init


!> Finalizes the timing module
subroutine timing_final
  use timing
  implicit none
  ! Remove all timers
  ntimers=0
  ! Deallocate
  deallocate(mval);   nullify(mval)
  deallocate(timers); nullify(timers)
  return
end subroutine timing_final


!> Get all the timers and prepare for monitoring
subroutine timing_monitor_init
  use timing
  implicit none
  integer :: itimer
  ! Allocate the array of data to transfer to monitor
  allocate(mval(2*ntimers+3))
  ! Create the file to monitor
  call monitor_create_file('timing',2*ntimers+3,1)
  call monitor_set_header(1,'Total [s]','r')
  do itimer=1,ntimers
     call monitor_set_header(itimer*2+0,trim(timers(itimer)%name)//' [s]','r')
     call monitor_set_header(itimer*2+1,trim(timers(itimer)%name)//' [%]','r')
  end do
  call monitor_set_header(2*ntimers+2,'Rest [s]','r')
  call monitor_set_header(2*ntimers+3,'Rest [%]','r')
  return
end subroutine timing_monitor_init


!> Create a new timer object
subroutine timing_create(name)
  use timing
  implicit none
  character(len=*), intent(in) :: name
  type(timer), dimension(:), allocatable :: temp
  ! Add a timer
  ntimers=ntimers+1
  allocate(temp(ntimers)); temp(1:size(timers))=timers; call move_alloc(temp,timers)
  ! Preset the values
  timers(ntimers)%name=trim(name)
  timers(ntimers)%time=0.0_WP
  timers(ntimers)%time_in=0.0_WP
  timers(ntimers)%started=.false.
  return
end subroutine timing_create


!> Start the timer
subroutine timing_start(name)
  use timing
  use parallel
  implicit none
  character(len=*), intent(in) :: name
  integer :: itimer
  ! Find the timer
  call timing_find_timer(name,itimer)
  ! Check if started already
  if (timers(itimer)%started) call warn('[timing_start] Timer already started: '//trim(name))
  ! Get the time and mark timer as started
  call parallel_time(timers(itimer)%time_in)
  timers(itimer)%started=.true.
  return
end subroutine timing_start


!> Stop the timer
subroutine timing_stop(name)
  use timing
  use parallel
  implicit none
  character(len=*), intent(in) :: name
  real(WP) :: newtime
  integer :: itimer
  ! Find the timer
  call timing_find_timer(name,itimer)
  ! Check timer was previously started
  if (.not.timers(itimer)%started) call warn('[timing_stop] Timer not started: '//trim(name))
  ! Get time and stop timer
  call parallel_time(newtime)
  timers(itimer)%time=timers(itimer)%time+newtime-timers(itimer)%time_in
  timers(itimer)%time_in=0.0_WP
  timers(itimer)%started=.false.
  return
end subroutine timing_stop


!> Monitor the timers
subroutine timing_monitor
  use timing
  use parallel
  implicit none
  integer :: itimer
  real(WP) :: full_time_out,rest_time
  ! Get time for full time step
  call parallel_time(full_time_out)
  ! Create the array of values
  rest_time=0.0_WP
  mval(1)=full_time_out-full_time_in
  do itimer=1,ntimers
     rest_time=rest_time+timers(itimer)%time
     mval(2*itimer+0)=timers(itimer)%time
     mval(2*itimer+1)=100.0_WP*timers(itimer)%time / mval(1)
  end do
  mval(2*ntimers+2)=mval(1)-rest_time
  mval(2*ntimers+3)=100.0_WP*(mval(1)-rest_time)/mval(1)
  ! Transfer data to monitor
  call monitor_select_file('timing')
  call monitor_set_array_values(mval)
  ! Reset all timers
  full_time_in=full_time_out
  do itimer=1,ntimers
     timers(itimer)%time=0.0_WP
     timers(itimer)%time_in=0.0_WP
     timers(itimer)%started=.false.
  end do
  return
end subroutine timing_monitor

