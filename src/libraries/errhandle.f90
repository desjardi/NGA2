!> This is an error handling module for NGA.
module errhandle
  implicit none
  private
  public :: die,warn
  
contains
  
  !> Subroutine that violently interrupts code execution 
  subroutine die(error_text)
    !use monitor,  only: monitor_log
    use parallel, only: parallel_kill
    implicit none
    character(len=*), intent(in) :: error_text
    !call monitor_log("[KILLED] "//trim(error_text))
    call parallel_kill(error_text)
  end subroutine die
  
  !> Subroutine that warns of issue
  subroutine warn(warn_text)
    !use monitor,  only: monitor_log
    implicit none
    character(len=*), intent(in) :: warn_text
    print*,'WARN WARN WARN'//trim(warn_text)
    !call monitor_log("[WARNING] "//trim(warn_text))
  end subroutine warn
  
end module errhandle
