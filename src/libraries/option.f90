!> Module that handles the reading of command line options
!> Expected format is: nga.exe [options] filename.
module option
  use string, only: str_medium,str_long
  implicit none
  private
  
  public :: option_init,option_final
  public :: option_read,option_is_invoked,option_print
  
  !> Type for command line options
  type clopt
     character(len=str_medium) :: lgn     !< Long name
     character(len=1)          :: shn     !< Short name
     logical                   :: has_val !< Whether the option has an argument
     character(len=str_medium) :: val     !< Argument value
  end type clopt
  
  !> Array of command line options
  type(clopt), dimension(:), allocatable :: opt
  integer :: nopts
  
  !> Array of remaining command line options
  character(len=str_medium), dimension(:), allocatable :: ropt
  integer :: nropts
  
  !> Input file name
  logical :: given_inputfile
  character(len=str_medium) :: inputfile
  
  !> Obtain information about specific option
  interface option_read
     module procedure option_read_logical
     module procedure option_read_integer
     module procedure option_read_real
     module procedure option_read_character
  end interface option_read
  
contains
  
  !> Initialize and parse all command line options and store them: for X=any character
  !>   - '-X' is a valid short option
  !>   - '-X val' is a valid short option with value 'val'
  !>   - '--XXXX' is a valid long option
  !>   - '--XXXX=val' is a valid long option with value 'val'
  !>   .
  subroutine option_init
    use errhandle, only: die
    implicit none
    character(len=str_medium), save :: arg
    integer :: pos,eqpos
    logical :: is_there
    
    ! Allocate option arrays to sufficient size
    allocate( opt(command_argument_count()))
    
    ! Allocate remaining option array to sufficient size
    allocate(ropt(command_argument_count()))
    
    ! Traverse all arguments and analyze them
    pos=1
    nopts=0
    nropts=0
    given_inputfile=.false.
    do while (pos.le.command_argument_count())
       
       ! Read the next argument
       call get_command_argument(pos,arg)
       
       ! Sort the arguments
       if (len_trim(arg).eq.2.and.arg(1:1).eq.'-'.and.arg(2:2).ne.'-') then ! This is a valid short option

          ! Check if option appears more than once
          call option_is_invoked(is_there,shn=arg(2:2))
          if (is_there) call die('option_init: found multiple instance of option -'//arg(2:2))
          
          ! Add it to our array
          nopts=nopts+1
          opt(nopts)%shn=arg(2:2)
          opt(nopts)%lgn=''
          opt(nopts)%val=''
          opt(nopts)%has_val=.false.
          
          ! Check if option has an associated value
          ! Careful here about the input file, hence the -1
          if (pos.lt.command_argument_count()-1) then
             
             ! Move forward
             pos=pos+1
             
             ! Read argument
             call get_command_argument(pos,arg)
             
             ! If it looks like another option, rewind
             if (arg(1:1).eq.'-') then
                pos=pos-1
             else
                opt(nopts)%val=trim(arg(:))
                opt(nopts)%has_val=.true.
             end if
             
          end if
          
       else if (len_trim(arg).gt.2.and.arg(1:1).eq.'-'.and.arg(2:2).eq.'-') then ! This is a valid long option

          ! Check if option has a value
          eqpos=scan(arg(3:),'=')
          
          ! Check if option appears more than once
          if (eqpos.eq.0) then
             call option_is_invoked(is_there,lgn=trim(arg(3:)))
             if (is_there) call die('option_init: found multiple instance of option --'//trim(arg(3:)))
          else
             call option_is_invoked(is_there,lgn=trim(arg(3:eqpos+1)))
             if (is_there) call die('option_init: found multiple instance of option --'//trim(arg(3:eqpos+1)))
          end if
          
          ! Add it to our array
          nopts=nopts+1
          opt(nopts)%shn=''
          opt(nopts)%val=''
          opt(nopts)%has_val=.false.
          
          ! Careful about value
          if (eqpos.eq.0) then
             opt(nopts)%lgn=trim(arg(3:))
          else
             opt(nopts)%lgn=trim(arg(3:eqpos+1))
             opt(nopts)%val=trim(arg(eqpos+3:))
             opt(nopts)%has_val=.true.
          end if
          
       else if (pos.eq.command_argument_count().and.arg(1:1).ne.'-') then ! This is the input file
          
          ! This argument should be the input file
          given_inputfile=.true.
          inputfile=trim(arg)
          
       else ! Everything else
          
          ! Add argument to remaining options
          nropts=nropts+1
          ropt(nropts)=trim(arg)
          
       end if
       
       ! Increment position
       pos=pos+1
       
    end do
    
  end subroutine option_init
  
  
  !> Check if option has been invoked
  subroutine option_is_invoked(is_def,shn,lgn,pos)
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
             if (is_def) call die('option_is_invoked: found both short and long version of option --'//trim(lgn))
             is_def=.true.
             if (present(pos)) pos=i
          end if
       end do
    end if
    
  end subroutine option_is_invoked
  
  
  !> Option_read for integer
  subroutine option_read_integer(val,shn,lgn,has_val,default)
    use errhandle, only: die
    implicit none
    integer,          intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    logical,          intent(out), optional :: has_val
    integer,          intent(in ), optional :: default
    logical :: is_def
    integer :: pos
    
    ! Find the option
    call option_is_invoked(is_def,shn=shn,lgn=lgn,pos=pos)
    
    ! Return appropriate information
    if (is_def) then
       if (present(has_val)) has_val=opt(pos)%has_val
       if (opt(pos)%has_val) then
          read(opt(pos)%val,*) val
       else
          if (present(default)) then
             val=default
          else
             if (opt(pos)%shn.eq.'') call die('option_read: required value not found for option --'//trim(opt(pos)%lgn))
             if (opt(pos)%lgn.eq.'') call die('option_read: required value not found for option -' //trim(opt(pos)%shn))
          end if
       end if
    else
       if (present(has_val)) has_val=.false.
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('option_read: option -' //shn//' not found, not default provided')
          if (present(lgn)) call die('option_read: option --'//lgn//' not found, not default provided')
       end if
    end if
    
  end subroutine option_read_integer

  
  !> Option_read for logical
  subroutine option_read_logical(val,shn,lgn,has_val,default)
    use errhandle, only: die
    implicit none
    logical,          intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    logical,          intent(out), optional :: has_val
    logical,          intent(in ), optional :: default
    logical :: is_def
    integer :: pos
    
    ! Find the option
    call option_is_invoked(is_def,shn=shn,lgn=lgn,pos=pos)
    
    ! Return appropriate information
    if (is_def) then
       if (present(has_val)) has_val=opt(pos)%has_val
       if (opt(pos)%has_val) then
          read(opt(pos)%val,*) val
       else
          if (present(default)) then
             val=default
          else
             if (opt(pos)%shn.eq.'') call die('option_read: required value not found for option --'//trim(opt(pos)%lgn))
             if (opt(pos)%lgn.eq.'') call die('option_read: required value not found for option -' //trim(opt(pos)%shn))
          end if
       end if
    else
       if (present(has_val)) has_val=.false.
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('option_read: option -' //shn//' not found, not default provided')
          if (present(lgn)) call die('option_read: option --'//lgn//' not found, not default provided')
       end if
    end if
    
  end subroutine option_read_logical

  
  !> Option_read for real
  subroutine option_read_real(val,shn,lgn,has_val,default)
    use errhandle, only: die
    use precision, only: WP
    implicit none
    real(WP),         intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    logical,          intent(out), optional :: has_val
    real(WP),         intent(in ), optional :: default
    logical :: is_def
    integer :: pos
    
    ! Find the option
    call option_is_invoked(is_def,shn=shn,lgn=lgn,pos=pos)
    
    ! Return appropriate information
    if (is_def) then
       if (present(has_val)) has_val=opt(pos)%has_val
       if (opt(pos)%has_val) then
          read(opt(pos)%val,*) val
       else
          if (present(default)) then
             val=default
          else
             if (opt(pos)%shn.eq.'') call die('option_read: required value not found for option --'//trim(opt(pos)%lgn))
             if (opt(pos)%lgn.eq.'') call die('option_read: required value not found for option -' //trim(opt(pos)%shn))
          end if
       end if
    else
       if (present(has_val)) has_val=.false.
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('option_read: option -' //shn//' not found, not default provided')
          if (present(lgn)) call die('option_read: option --'//lgn//' not found, not default provided')
       end if
    end if
    
  end subroutine option_read_real

  
  !> Option_read for character
  subroutine option_read_character(val,shn,lgn,has_val,default)
    use errhandle, only: die
    implicit none
    character(len=*), intent(out)           :: val
    character(len=*), intent(in ), optional :: shn
    character(len=*), intent(in ), optional :: lgn
    logical,          intent(out), optional :: has_val
    character(len=*), intent(in ), optional :: default
    logical :: is_def
    integer :: pos
    
    ! Find the option
    call option_is_invoked(is_def,shn=shn,lgn=lgn,pos=pos)
    
    ! Return appropriate information
    if (is_def) then
       if (present(has_val)) has_val=opt(pos)%has_val
       if (opt(pos)%has_val) then
          read(opt(pos)%val,'(a)') val
       else
          if (present(default)) then
             val=default
          else
             if (opt(pos)%shn.eq.'') call die('option_read: required value not found for option --'//trim(opt(pos)%lgn))
             if (opt(pos)%lgn.eq.'') call die('option_read: required value not found for option -' //trim(opt(pos)%shn))
          end if
       end if
    else
       if (present(has_val)) has_val=.false.
       if (present(default)) then
          val=default
       else
          if (present(shn)) call die('option_read: option -' //shn//' not found, not default provided')
          if (present(lgn)) call die('option_read: option --'//lgn//' not found, not default provided')
       end if
    end if
    
  end subroutine option_read_character
  
  
  !> Print all detected options
  subroutine option_print
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    integer :: i
    
    ! Loop over all options and print then out
    if (nopts.eq.0) then
       write(output_unit,'(a)') 'NGA was called without any options.'
    else
       write(output_unit,'(a,1x,i0,1x,a)') 'NGA was called with',nopts,'options:'
       do i=1,nopts
          if (opt(i)%shn.eq.'') then
             write(output_unit,'(3x,a,1x,"--",a)') 'Option: ',trim(opt(i)%lgn)
          else
             write(output_unit,'(3x,a,1x," -",a)') 'Option: ',trim(opt(i)%shn)
          end if
          if (opt(i)%has_val) then
             write(output_unit,'(10x,a,1x,a)') 'Value set to:',trim(opt(i)%val)
          else
             write(output_unit,'(10x,a)')      'No corresponding value'
          end if
       end do
    end if
    
    ! Info on input file detected
    if (given_inputfile) then
       write(output_unit,'(a,1x,a)') 'The following input file was provided:',trim(inputfile)
    else
       write(output_unit,'(a,1x,a)') 'No input file was provided.'
    end if
    
    ! Info on remaining arguments
    if (nropts.gt.0) then
       write(output_unit,'(a)') 'The following arguments were not recognized, and ignored:'
       do i=1,nropts
          write(output_unit,'(3x,a)') trim(ropt(i))
       end do
    end if
    
  end subroutine option_print
  
  
  !> Finalize and clean up options module
  subroutine option_final
    implicit none
    deallocate( opt);  nopts=0
    deallocate(ropt); nropts=0
    inputfile=''; given_inputfile=.false.
  end subroutine option_final
  
end module option
