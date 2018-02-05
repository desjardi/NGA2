!> Module that handles the reading of command line options
!> Expected format is: nga.exe [options] filename.
!> @todo Deallocation/finalization
!> @todo Capability to seek a value from stored command line options
module options
  use string, only: str_medium,str_long
  implicit none
  private
  
  public :: options_read,options_print
  
  !> Type for command line options
  type clopt
     logical                   :: is_long !< Whether it is a long option
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
  
contains
  
  !> Parse all command line options and store them: for X=any character
  !>   - '-X' is a valid short option
  !>   - '-X val' is a valid short option with value 'val'
  !>   - '--XXXX' is a valid long option
  !>   - '--XXXX=val' is a valid long option with value 'val'
  !>   .
  subroutine options_read
    use errhandle, only: warn,die
    implicit none
    character(len=str_medium), save :: arg
    integer :: pos,stat,eqpos
    
    ! Check number of options
    nopts=0
    nropts=0
    do pos=1,command_argument_count()
       ! Read in the argument
       call get_command_argument(pos,arg,status=stat)
       ! Check it fits in our string
       if (stat.ne.0) call die('Command line option '//trim(arg)//' is too long.')
       ! Check if the argument is a valid option
       if (len_trim(arg).eq.2.and.arg(1:1).eq.'-') then
          nopts=nopts+1
       else if (len_trim(arg).gt.2.and.arg(1:1).eq.'-'.and.arg(2:2).eq.'-') then
          nopts=nopts+1
       else if (pos.eq.command_argument_count().and.arg(1:1).ne.'-') then
          ! Input file
       else
          ! Increment remaining argument
          nropts=nropts+1
       end if
    end do
    
    ! Allocate option array to appropriate size
    allocate(opt(nopts))
    
    ! Allocate remaining option array to appropriate size
    allocate(ropt(nropts))
    
    ! Traverse all arguments and analyze them
    pos=1
    nopts=0
    nropts=0
    given_inputfile=.false.
    do while (pos.le.command_argument_count())
       
       ! Read the next argument
       call get_command_argument(pos,arg)
       
       ! Sort the argument
       if (len_trim(arg).eq.2.and.arg(1:1).eq.'-') then
          ! This is a valid short option
          
          ! Add it to our array
          nopts=nopts+1
          opt(nopts)%is_long=.false.
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
          
       else if (len_trim(arg).gt.2.and.arg(1:1).eq.'-'.and.arg(2:2).eq.'-') then
          ! This is a valid long option
          
          ! Check if option has a value
          eqpos=scan(arg(3:),'=')
          
          ! Add it to our array
          nopts=nopts+1
          opt(nopts)%is_long=.true.
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
          
       else if (pos.eq.command_argument_count().and.arg(1:1).ne.'-') then
          ! This is the input file
          
          ! This argument should be the input file
          given_inputfile=.true.
          inputfile=trim(arg)
          
       else
          ! Everything else
          
          ! Add argument to remaining options
          nropts=nropts+1
          ropt(nropts)=trim(arg)
          
       end if
       
       ! Increment position
       pos=pos+1
       
    end do
    
  end subroutine options_read
  
  
  !> Print all detected options
  subroutine options_print
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    integer :: i
    ! Loop over all options and print then out
    if (nopts.eq.0) then
       write(output_unit,'(a)') 'NGA was called without any options.'
    else
       write(output_unit,'(a,x,i0,x,a)') 'NGA was called with',nopts,'options:'
       do i=1,nopts
          if (opt(i)%is_long) then
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
  end subroutine options_print
  
  
end module options
