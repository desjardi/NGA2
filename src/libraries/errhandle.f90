!> This is an error handling module for NGA.
module errhandle
  implicit none
  ! Dummy
end module errhandle


!> Subroutine that violently interrupts code execution 
subroutine die(error_text)
  use errhandle
  implicit none
  character(len=*), intent(in) :: error_text
  call monitor_log("[KILLED] "//trim(error_text))
  call parallel_kill(error_text)
  return
end subroutine die


!> Subroutine that warns of issue
subroutine warn(warn_text)
  use errhandle
  implicit none
  character(len=*), intent(in) :: warn_text
  call monitor_log("[WARNING] "//trim(warn_text))
  return
end subroutine warn
