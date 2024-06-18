!> This module handles user-defined parameters
!> from command line and input files.
module param
   use string, only: str_medium,str_long
   implicit none
   private
   
   public :: param_init,param_final,param_read
   public :: param_exists,param_getsize
   
   !> Verbosity of this run
   integer, public, protected :: verbose
   
   !> Type for storing parameters from command line and parser
   type :: param_type
      character(str_medium) :: tag !< Tag name for field
      character(str_long)   :: val !< Associated value
      character(str_medium) :: src !< Source of the parameter
   end type param_type
   
   !> This array stores all user-defined parameters
   type(param_type), dimension(:), allocatable :: params
   integer :: nparams
   
   !> Routine that reads a required value from the user
   interface param_read
      module procedure param_readlogical
      module procedure param_readint
      module procedure param_readintarray
      module procedure param_readfloat
      module procedure param_readfloatarray
      module procedure param_readfloatarray2D
      module procedure param_readchar
      module procedure param_readchararray
   end interface param_read
   
contains
   
   
   !> Initialize the module handling user parameters
   subroutine param_init
      use messager, only: die
      implicit none
      character(len=str_long) :: arg,val
      character(len=str_medium), dimension(:), allocatable :: files
      integer :: pos,eqpos,i
      
      ! First ensure empty storage
      call param_final
      
      ! Start by parsing the command line
      pos=1
      do while (pos.le.command_argument_count())
         ! Read the next argument
         call get_command_argument(pos,arg)
         ! Test whether argument conforms to one of the two accepted formats
         if (len_trim(arg).eq.2.and.arg(1:1).eq.'-'.and.arg(2:2).ne.'-') then ! This is a valid short-format parameter
            ! Check if parameter has an associated value
            if (pos.lt.command_argument_count()) then
               ! Move forward and read one more argument
               pos=pos+1; call get_command_argument(pos,val)
               ! If it looks like another option, rewind
               if (val(1:1).eq.'-') then
                  ! Rewind
                  pos=pos-1
                  ! Add parameter without associated value
                  call param_add(arg(2:2),'','command line: short')
               else
                  ! Add parameter with associated value
                  call param_add(arg(2:2),adjustl(trim(val)),'command line: short')
               end if
            else
               ! This was already the last option, so add parameter without associated value
               call param_add(arg(2:2),'','command line: short')
            end if
         else if (len_trim(arg).gt.2.and.arg(1:1).eq.'-'.and.arg(2:2).eq.'-') then ! This is a valid long-format parameter
            ! Check if parameter has an associated value
            eqpos=scan(arg,'=')
            if (eqpos.eq.0) then
               ! Add parameter without associated value
               call param_add(adjustl(trim(arg(3:))),'','command line: long')
            else
               ! Add parameter with associated value
               call param_add(adjustl(trim(arg(3:eqpos-1))),adjustl(trim(arg(eqpos+1:))),'command line: long')
            end if
         else
            ! Parameter does not follow expected format...
            call die('[param_init] The following argument does not follow the expected format: '//trim(arg))
         end if
         ! Increment position
         pos=pos+1
      end do
      
      ! Before anything else, check if help has been invoked
      if (param_exists('help',short='h')) call param_help
      
      ! Now continue with user-defined parameters coming from input files
      if (param_exists('input',short='i')) then
         ! Retrieve file name
         allocate(files(param_getsize('input','i')))
         call param_read('input',files,short='i')
         do i=1,size(files)
            call param_parsefile(files(i))
         end do
      end if
      
      ! Check verbosity level
      call param_read('verbose',verbose,short='v',default=1)
      if (verbose.gt.0) call param_print
      
   end subroutine param_init
   
   
   !> Read & parse parameter files
   subroutine param_parsefile(input)
      use messager, only: die
      implicit none
      character(len=*) :: input !< This is the name of the file to be read and parsed.
      integer :: iunit,ierr,limiter,nlines,i,j,ntags,comment
      integer, dimension(:), allocatable :: limit,line
      character(len=str_long) :: buffer
      character(len=str_long), dimension(:), allocatable :: file
      character(len=str_long) :: val
      character(len=str_medium) :: tag
      ! Open the file
      open(newunit=iunit,file=input,form='formatted',status='old',iostat=ierr)
      if (ierr.ne.0) call die('[param_parsefile] Could not open file: '//trim(input))
      ! Count the number of lines in the file
      ierr=0
      nlines=0
      do while (ierr.eq.0)
         read(iunit,'(a)',iostat=ierr) buffer
         nlines=nlines+1
      end do
      rewind(iunit)
      ! Allocate to the right size
      allocate(file(nlines+1),limit(nlines+1),line(nlines+1))
      ! Read everything in the buffer
      ierr=0
      nlines=0
      loop: do while (ierr.eq.0)
         read(iunit,'(a)',iostat=ierr) buffer
         if (ierr.ne.0) exit loop
         ! Find comments
         comment=scan(buffer,'!#%')
         ! Remove them
         if (comment.ne.0) buffer(comment:)=''
         ! Trim
         buffer=adjustl(buffer)
         ! Add line
         if (len_trim(buffer).ne.0) then
            nlines=nlines+1
            file(nlines)=buffer
         end if
      end do loop
      ! Close the file
      close(iunit)
      ! Get the tags
      ntags=0
      do i=1,nlines
         limiter=scan(file(i),':=')
         if (limiter.ne.0) then
            ntags=ntags+1
            line(ntags)=i
            line(ntags+1)=nlines+1
            limit(ntags)=limiter
         end if
      end do
      ! Read everything now
      do i=1,ntags
         buffer=''
         do j=line(i),line(i+1)-1
            if (j.eq.line(i)) then
               buffer=trim(buffer)//trim(file(j))
            else
               buffer=trim(buffer)//' '//trim(file(j))
            end if
         end do
         read(buffer(1:limit(i)-1),'(a)') tag
         read(buffer(limit(i)+1:),'(a)') val
         ! Add parameter
         call param_add(adjustl(trim(tag)),adjustl(trim(val)),'file: '//adjustl(trim(input)))
      end do
      ! Deallocate
      deallocate(file,limit,line)
   end subroutine param_parsefile
   
   
   !> Add a new parameter in the param array:
   !>   - if param has been found already, replace it
   !>   - if param has not been found yet, add it
   !>   .
   subroutine param_add(tag,val,src)
      use messager, only: warn
      implicit none
      character(len=*), intent(in) :: tag
      character(len=*), intent(in) :: val
      character(len=*), intent(in) :: src
      type(param_type), dimension(:), allocatable :: temp
      integer :: i
      ! Check if parameter exists already
      if (param_exists(tag)) then
         ! Warn the user
         call warn('[param_add] The following parameter appeared more than once: '//trim(tag))
         call warn('[param_add] Only the first occurence is considered!')
      else
         ! Resize the parameter array
         allocate(temp(nparams+1))
         if (nparams.gt.0) temp(1:nparams)=params
         call move_alloc(temp,params)
         ! Add parameter
         nparams=nparams+1
         params(nparams)%tag=adjustl(trim(tag))
         params(nparams)%val=adjustl(trim(val))
         params(nparams)%src=adjustl(trim(src))
         ! Remove unwanted characters from val string
         do i=1,len_trim(params(nparams)%val)
            if (ichar(params(nparams)%val(i:i)).eq.9.or.params(nparams)%val(i:i).eq.','.or.params(nparams)%val(i:i).eq.';') params(nparams)%val(i:i)=' '
         end do
      end if
   end subroutine param_add
   
   
   !> Check if a parameter already exists in params array
   function param_exists(tag,short) result(found)
      implicit none
      character(len=*), intent(in) :: tag
      character(len=*), intent(in), optional :: short
      logical :: found
      integer :: i
      ! Initializer
      found=.false.
      ! Traverse params array and look for tag
      do i=1,nparams
         if (params(i)%tag.eq.tag) found=.true.
      end do
      ! If short name was also provided, try to find it too
      if (present(short)) then
         do i=1,nparams
            if (params(i)%tag.eq.short) found=.true.
         end do
      end if
   end function param_exists
   
   
   !> Return parameter index for tag
   function param_index(tag) result(ind)
      implicit none
      character(len=*), intent(in) :: tag
      integer :: i,ind
      ! Initializer
      ind=0
      ! Traverse params array and look for tag
      do i=1,nparams
         if (params(i)%tag.eq.tag) ind=i
      end do
   end function param_index
   
   
   !> Get size of the parameter value for reading arrays
   function param_getsize(tag,short) result(count)
      use messager, only: die
      implicit none
      character(len=*), intent(in)           :: tag
      character(len=*), intent(in), optional :: short
      integer :: i,ind,count
      ! Initialize counter
      count=0
      ! Get the parameter index
      ind=param_index(tag)
      if (ind.gt.0) then
         ! Count the number of values provided
         if (params(ind)%val(1:1).ne.' ') count=1
         do i=2,len_trim(params(ind)%val)
            if (params(ind)%val(i:i).eq.' '.and.params(ind)%val(i-1:i-1).ne.' ') count=count+1
         end do
         return
      end if
      ! If we are still here, try the short name
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         if (ind.gt.0) then
            ! Count the number of values provided
            if (params(ind)%val(1:1).ne.' ') count=1
            do i=2,len_trim(params(ind)%val)
               if (params(ind)%val(i:i).eq.' '.and.params(ind)%val(i-1:i-1).ne.' ') count=count+1
            end do
            return
         end if
      end if
      ! Still here: we failed to find a size
      call die('[param_getsize] Did not find parameter: '//trim(tag))
   end function param_getsize
   
   
   !> Read integer value associated with parameter tag
   subroutine param_readint(tag,val,short,default)
      use messager, only: die
      implicit none
      character(len=*), intent(in)            :: tag
      integer         , intent(out)           :: val
      character(len=*), intent(in) , optional :: short
      integer         , intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            ! Parameter is defined
            if (len_trim(params(ind)%val).gt.0) then
               ! Value exists, so read it and return
               read(params(ind)%val,*) val
               return
            end if
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[param_read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine param_readint
   
   
   !> Read logical value associated with parameter tag
   subroutine param_readlogical(tag,val,short,default)
      use messager, only: die
      implicit none
      character(len=*), intent(in)            :: tag
      logical         , intent(out)           :: val
      character(len=*), intent(in) , optional :: short
      logical         , intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            ! Parameter is defined
            if (len_trim(params(ind)%val).gt.0) then
               ! Value exists, so read it and return
               read(params(ind)%val,*) val
               return
            end if
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[param_read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine param_readlogical
   
   
   !> Read real value associated with parameter tag
   subroutine param_readfloat(tag,val,short,default)
      use messager,  only: die
      use precision, only: WP
      implicit none
      character(len=*), intent(in)            :: tag
      real(WP)        , intent(out)           :: val
      character(len=*), intent(in) , optional :: short
      real(WP)        , intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            ! Parameter is defined
            if (len_trim(params(ind)%val).gt.0) then
               ! Value exists, so read it and return
               read(params(ind)%val,*) val
               return
            end if
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[param_read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine param_readfloat
   
   
   !> Read character value associated with parameter tag
   subroutine param_readchar(tag,val,short,default)
      use messager, only: die
      implicit none
      character(len=*), intent(in)            :: tag
      character(len=*), intent(out)           :: val
      character(len=*), intent(in) , optional :: short
      character(len=*), intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(params(ind)%val,'(a)') val
            return
         end if
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            ! Parameter is defined
            if (len_trim(params(ind)%val).gt.0) then
               ! Value exists, so read it and return
               read(params(ind)%val,'(a)') val
               return
            end if
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[param_read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine param_readchar
   
   
   !> Read integer array value associated with parameter tag
   subroutine param_readintarray(tag,val,short)
      use messager, only: die
      implicit none
      character(len=*), intent(in)            :: tag
      integer, dimension(:), intent(out)      :: val
      character(len=*), intent(in) , optional :: short
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(params(ind)%val,*) val
         return
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      call die('[param_read] Did not find required parameter: '//trim(tag))
   end subroutine param_readintarray
   
   
   !> Read float array value associated with parameter tag
   subroutine param_readfloatarray(tag,val,short,default)
      use messager,  only: die
      use precision, only: WP
      implicit none
      character(len=*), intent(in)            :: tag
      real(WP), dimension(:), intent(out)     :: val
      character(len=*), intent(in) , optional :: short
      real(WP), dimension(:), intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(params(ind)%val,*) val
         return
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[param_read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine param_readfloatarray
   
   
   !> Read float 2D array value associated with parameter tag
   subroutine param_readfloatarray2D(tag,val,short)
      use messager,  only: die
      use precision, only: WP
      implicit none
      character(len=*), intent(in)            :: tag
      real(WP), dimension(:,:), intent(out)   :: val
      character(len=*), intent(in) , optional :: short
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(params(ind)%val,*) val
         return
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      call die('[param_read] Did not find required parameter: '//trim(tag))
   end subroutine param_readfloatarray2D
   
   
   !> Read character array value associated with parameter tag
   subroutine param_readchararray(tag,val,short)
      use messager, only: die
      implicit none
      character(len=*), intent(in)                :: tag
      character(len=*), dimension(:), intent(out) :: val
      character(len=*), intent(in) , optional     :: short
      integer :: ind
      ! Get the parameter index
      ind=param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(params(ind)%val,*) val
         return
      end if
      ! If short name has been provided, check that too
      if (present(short)) then
         ! Get the parameter index
         ind=param_index(short)
         ! Read as appropriate
         if (ind.gt.0) then
            read(params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      call die('[param_read] Did not find required parameter: '//trim(tag))
   end subroutine param_readchararray
   
   
   !> Print out help message
   subroutine param_help
      use parallel, only: amroot,parallel_final
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      if (amroot) then
         write(output_unit,'(a)') ' ____________________________________________________________ '
         write(output_unit,'(a)') '|                    _   _  ____    _                        |'
         write(output_unit,'(a)') '|                   | \ | |/ ___|  / \                       |'
         write(output_unit,'(a)') '|                   |  \| | |  _  / _ \                      |'
         write(output_unit,'(a)') '|                   | |\  | |_| |/ ___ \                     |'
         write(output_unit,'(a)') '|                   |_| \_|\____/_/   \_\                    |'
         write(output_unit,'(a)') '|                                                            |'
         write(output_unit,'(a)') '| Multiphase Reactive Turbulent Flow Solver by O. Desjardins |'
         write(output_unit,'(a)') '|                   ---------------------                    |'
         write(output_unit,'(a)') '| The following format is expected when calling NGA:         |'
         write(output_unit,'(a)') '| ]> nga.().exe  [-X] [-X value] [--XXX] [--XXX=value]       |'
         write(output_unit,'(a)') '|                                                            |'
         write(output_unit,'(a)') '| Examples:                                                  |'
         write(output_unit,'(a)') '| ]> nga.().exe -i myinput or nga.().exe --input=myinput     |'
         write(output_unit,'(a)') '|        --> Read additional parameters from myinput         |'
         write(output_unit,'(a)') '| ]> nga.().exe -v or nga.().exe --verbose                   |'
         write(output_unit,'(a)') '|        --> Set run to higher verbosity level 1             |'
         write(output_unit,'(a)') '| ]> nga.().exe -i input --"Processors along X"=4            |'
         write(output_unit,'(a)') '|        --> Read input but set "Processors along X" to 4    |'
         write(output_unit,'(a)') '|        --> Command line takes precedence over input files  |'
         write(output_unit,'(a)') '|                                                            |'
         write(output_unit,'(a)') '| For a comprehensive list of parameters, run make param     |'
         write(output_unit,'(a)') '|____________________________________________________________|'
      end if
      call parallel_final; stop
   end subroutine param_help
   
   
   !> Print all user-defined parameters if root
   subroutine param_print
      use messager, only: log
      use parallel, only: amRoot
      implicit none
      integer :: i
      integer, parameter :: colwidth=30
      character(len=colwidth) :: tag,val,src,sep
      character(len=str_long) :: message
      ! Loop over all options and print then out
      if (amRoot) then
         if (nparams.eq.0) then
            call log('NGA was called without any user-defined parameter.')
         else
            write(message,'(a,1x,i0,1x,a)') 'NGA was called with',nparams,'user-defined parameters:'; call log(message)
            do i=1,colwidth
               sep(i:i)='_'
            end do
            write(message,'(1x,"_",a,"___",a,"___",a,"_",1x)') sep,sep,sep; call log(message)
            tag='Parameter tag name'; val='Assigned value'; src='Source'
            write(message,'("| ",a," | ",a," | ",a," |")') tag,val,src; call log(message)
            do i=1,nparams
               ! Check for truncation of strings
               if (len_trim(params(i)%tag).gt.colwidth) then
                  tag=params(i)%tag(1:colwidth-3)//'...'
               else
                  tag=''; tag=params(i)%tag(1:len_trim(params(i)%tag))
               end if
               if (len_trim(params(i)%val).gt.colwidth) then
                  val=params(i)%val(1:colwidth-3)//'...'
               else
                  val=''; val=params(i)%val(1:len_trim(params(i)%val))
               end if
               if (len_trim(params(i)%src).gt.colwidth) then
                  src=params(i)%src(1:colwidth-3)//'...'
               else
                  src=''; src=params(i)%src(1:len_trim(params(i)%src))
               end if
               write(message,'("| ",a," | ",a," | ",a," |")') tag,val,src; call log(message)
            end do
            write(message,'("|_",a,"_|_",a,"_|_",a,"_|")') sep,sep,sep; call log(message)
         end if
      end if
   end subroutine param_print
   
   
   !> Clean up the param handling module
   subroutine param_final
      implicit none
      nparams=0
      if (allocated(params)) deallocate(params)
   end subroutine param_final
   
   
end module param
