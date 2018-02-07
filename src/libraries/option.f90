!> Module that handles the reading of command line options
!> Expected format is: nga.exe [options] filename.
module option
  use string, only: str_medium,str_long
  implicit none
  private
  
  public :: option_init,option_final,option_read,option_print,option_isinvoked
  
  !> Type for storing command line options
  type clopt
     character(len=str_medium) :: lgn !< Long name
     character(len=1)          :: shn !< Short name
     character(len=str_medium) :: val !< Argument value
  end type clopt
  
  !> Array of command line options
  type(clopt), dimension(:), allocatable :: opt
  integer :: nopts
  
  !> Obtain information about specific option
  interface option_read
     module procedure option_readlogical
     module procedure option_readint
     module procedure option_readfloat
     module procedure option_readchar
  end interface option_read
  
contains

  
  !> Add a new option in the option array:
  !>   - if the option has been found already, replace it
  !>   - if the option has not been found yet, add it
  !>   .
  subroutine option_addoption(shn,lgn,val)
    use errhandle, only: warn,die
    implicit none
    character(len=*), intent(in), optional :: shn
    character(len=*), intent(in), optional :: lgn
    character(len=*), intent(in), optional :: val
    integer :: pos
    logical :: is_def
    type(clopt), dimension(:), allocatable :: temp
    if (present(shn)) then ! Short option version
       ! Check if option exists already
       call option_isinvoked(is_def,shn=shn,pos=pos)
       if (is_def) then
          ! Warn the user
          call warn('[option_addoption] The following command line option appeared more than once: -'//trim(shn))
          ! Replace the value if one is provided - leave the old one if not is provided
          if (present(val)) opt(pos)%val=adjustl(trim(val))
       else
          ! Resize the opt array
          nopts=nopts+1
          allocate(temp(nopts))
          if (nopts.gt.1) temp(1:nopts-1)=opt
          call move_alloc(temp,opt)
          ! Add option
          opt(nopts)%shn=adjustl(trim(shn))
          opt(nopts)%lgn=''
          if (present(val)) then
             opt(nopts)%val=adjustl(trim(val))
          else
             opt(nopts)%val=''
          end if
       end if
    else if (present(lgn)) then ! Long option version
       ! Check if option exists already
       call option_isinvoked(is_def,lgn=shn,pos=pos)
       if (is_def) then
          ! Warn the user
          call warn('[option_addoption] The following command line option appeared more than once: --'//trim(lgn))
          ! Replace the value if one is provided - leave the old one if not is provided
          if (present(val)) opt(pos)%val=adjustl(trim(val))
       else
          ! Resize the opt array
          nopts=nopts+1
          allocate(temp(nopts))
          if (nopts.gt.1) temp(1:nopts-1)=opt
          call move_alloc(temp,opt)
          ! Add option
          opt(nopts)%lgn=adjustl(trim(lgn))
          opt(nopts)%shn=''
          if (present(val)) then
             opt(nopts)%val=adjustl(trim(val))
          else
             opt(nopts)%val=''
          end if
       end if
    else
       call die('[option_addoption] Subroutine was called without short nor long option name.')
    end if
  end subroutine option_addoption
  
  
  !> Initialize and parse all command line options and store them: for X=any character
  !>   - '-X' is a valid short option
  !>   - '-X val' is a valid short option with value 'val'
  !>   - '--XXXX' is a valid long option
  !>   - '--XXXX=val' is a valid long option with value 'val'
  !>   .
  subroutine option_init
    use errhandle, only: die
    implicit none
    character(len=str_medium) :: arg,val
    integer :: pos,eqpos
    
    ! Ensure options array is empty
    nopts=0
    if (allocated(opt)) deallocate(opt)
    
    ! Traverse all arguments and analyze them
    pos=1
    do while (pos.le.command_argument_count())
       
       ! Read the next argument
       call get_command_argument(pos,arg)
       
       ! Sort the arguments into short format and long format
       if (len_trim(arg).eq.2.and.arg(1:1).eq.'-'.and.arg(2:2).ne.'-') then ! This is a valid short option
          
          ! Check if option has an associated value
          if (pos.lt.command_argument_count()) then
             ! Move forward and read one more argument
             pos=pos+1; call get_command_argument(pos,val)
             ! If it looks this is another option, rewind
             if (val(1:1).eq.'-') then
                ! Rewind
                pos=pos-1
                ! Add short option without associated value
                call option_addoption(shn=arg(2:2))
             else
                ! Add short option with associated value
                call option_addoption(shn=arg(2:2),val=trim(val))
             end if
          else
             ! This was already the last option, so not associated value
             call option_addoption(shn=arg(2:2),val=trim(val))
          end if
          
       else if (len_trim(arg).gt.2.and.arg(1:1).eq.'-'.and.arg(2:2).eq.'-') then ! This is a valid long option
          
          ! Check if option has an associated value
          eqpos=scan(arg(3:),'=')
          if (eqpos.eq.0) then
             ! Add long option without associated value
             call option_addoption(lgn=arg(3:))
          else
             ! Add long option with associated value
             call option_addoption(lgn=arg(3:eqpos+1),val=arg(eqpos+3:))
          end if
          
       else
          ! Option does not follow expected format...
          call die('[option_init] The following argument does not follow the expected format: '//trim(arg))
       end if
       
       ! Increment position
       pos=pos+1
       
    end do
    
  end subroutine option_init
  
  
  !> Check if option has been invoked
  subroutine option_isinvoked(is_def,shn,lgn,pos)
    use errhandle, only: die
    implicit none
    logical,          intent(out)           :: is_def
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    integer,          intent(out), optional :: pos
    integer :: i
    
    ! Assume option is not defined
    is_def=.false.
    if (present(pos)) pos=0
    
    ! First localize the option with its short name
    if (present(shn)) then
       do i=1,nopts
          if (opt(i)%shn.eq.shn) then
             is_def=.true.
             if (present(pos)) pos=i
          end if
       end do
    end if
    
    ! Then localize the option with its long name
    if (present(lgn)) then
       do i=1,nopts
          if (opt(i)%lgn.eq.lgn) then
             if (is_def) call die('[option_isinvoked] Found both short and long version of option --'//trim(lgn))
             is_def=.true.
             if (present(pos)) pos=i
          end if
       end do
    end if
    
  end subroutine option_isinvoked


  !> Return the option's associated value
  subroutine option_getvalue(val,shn,lgn)
    implicit none
    character(len=*), intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    logical :: is_def
    integer :: pos
    ! Find the option
    call option_isinvoked(is_def,shn=shn,lgn=lgn,pos=pos)
    ! Extract its associated value
    if (is_def) then
       val=opt(pos)%val
    else
       val=''
    end if
  end subroutine option_getvalue
  
  
  !> Option_read for integer
  subroutine option_readint(val,shn,lgn,default)
    use errhandle, only: die
    implicit none
    integer,          intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    integer,          intent(in ), optional :: default
    character(len=str_medium) :: myval
    ! Find the option
    call option_getvalue(myval,shn=shn,lgn=lgn)
    ! Return appropriate information
    if (len_trim(myval).gt.0) then
       read(myval,*) val
    else
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('[option_read] Option -' //shn//' value not found')
          if (present(lgn)) call die('[option_read] Option --'//lgn//' value not found')
       end if
    end if
  end subroutine option_readint
  
  
  !> Option_read for logical
  subroutine option_readlogical(val,shn,lgn,default)
    use errhandle, only: die
    implicit none
    logical,          intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    logical,          intent(in ), optional :: default
    character(len=str_medium) :: myval
    ! Find the option
    call option_getvalue(myval,shn=shn,lgn=lgn)
    ! Return appropriate information
    if (len_trim(myval).gt.0) then
       read(myval,*) val
    else
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('[option_read] Option -' //shn//' value not found')
          if (present(lgn)) call die('[option_read] Option --'//lgn//' value not found')
       end if
    end if
  end subroutine option_readlogical

  
  !> Option_read for real
  subroutine option_readfloat(val,shn,lgn,default)
    use errhandle, only: die
    use precision, only: WP
    implicit none
    real(WP),         intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    real(WP),         intent(in ), optional :: default
    character(len=str_medium) :: myval
    ! Find the option
    call option_getvalue(myval,shn=shn,lgn=lgn)
    ! Return appropriate information
    if (len_trim(myval).gt.0) then
       read(myval,*) val
    else
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('[option_read] Option -' //shn//' value not found')
          if (present(lgn)) call die('[option_read] Option --'//lgn//' value not found')
       end if
    end if
  end subroutine option_readfloat

  
  !> Option_read for character
  subroutine option_readchar(val,shn,lgn,default)
    use errhandle, only: die
    implicit none
    character(len=*), intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    character(len=*), intent(in ), optional :: default
    ! Find the option
    call option_getvalue(val,shn=shn,lgn=lgn)
    ! Return appropriate information
    if (len_trim(val).eq.0) then
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('[option_read] Option -' //shn//' value not found')
          if (present(lgn)) call die('[option_read] Option --'//lgn//' value not found')
       end if
    end if
  end subroutine option_readchar
  
  
  !> Print all detected options
  subroutine option_print
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    integer :: i
    character(len=str_medium) :: message
    ! Loop over all options and print then out
    if (nopts.eq.0) then
       write(output_unit,'(a)') '> NGA was called without any options.'
    else
       write(output_unit,'(a,1x,i0,1x,a)') '> NGA was called with',nopts,'options:'
       do i=1,nopts
          if (opt(i)%shn.eq.'') then
             message='   - Option: --'//trim(opt(i)%lgn)
          else
             message='   - Option:  -'//trim(opt(i)%shn)
          end if
          if (len_trim(opt(i)%val).gt.0) then
             message(30:)='   |   Value = '//trim(opt(i)%val)
          else
             message(30:)='   |   No assigned value'
          end if
          write(output_unit,'(4x,a)') adjustl(trim(message))
       end do
    end if
  end subroutine option_print
  
  
  !> Finalize and clean up options module
  subroutine option_final
    implicit none
    deallocate(opt); nopts=0
  end subroutine option_final
  
end module option
