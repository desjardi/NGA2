!> This class parses an input file
module inputfile_class
   use string, only: str_medium,str_long
   implicit none
   private
   
   !> Type for storing input file parameters
   type :: param_type
      character(str_medium) :: tag !< Tag name for field
      character(str_long)   :: val !< Associated value
   end type param_type
   
   ! Expose type/constructor
   public :: inputfile
   
   !> inputfile object definition
   type :: inputfile
      
      !< We need to know who's the boss
      logical :: amRoot
      
      !> Name of file to parse
      character(len=str_medium) :: filename
      
      !> Array of parameters
      type(param_type), dimension(:), allocatable :: params
      integer :: nparams
      
   contains
      procedure, private :: add_param        !< Add a parameter
      procedure, private :: param_exists     !< Check if parameter exists
      procedure, private :: param_index      !< Return index of parameter
      procedure, private :: readlogical      !< Read in a logical
      procedure, private :: readint          !< Read in an integer
      procedure, private :: readintarray     !< Read in an integer array
      procedure, private :: readfloat        !< Read in a float
      procedure, private :: readfloatarray   !< Read in a float array
      procedure, private :: readfloatarray2D !< Read in a float 2D array
      procedure, private :: readchar         !< Read in a string
      procedure, private :: readchararray    !< Read in a string array
      generic :: read=>readlogical,readint,readintarray,readfloat,readfloatarray,readfloatarray2D,readchar,readchararray !< Generic routine to read a parameter
      procedure :: print                     !< Print out all parameters
      procedure :: log=>inputlog             !< Log all parameters
   end type inputfile
   
   
   !> Declare inputfile constructor
   interface inputfile
      procedure constructor
   end interface inputfile
   

contains
   
   
   !> Constructor of inputfile object
   function constructor(amRoot,filename) result(this)
      use messager, only: die
      implicit none
      type(inputfile) :: this
      logical, intent(in) :: amRoot
      character(len=*), intent(in) :: filename
      integer :: iunit,ierr,limiter,nlines,i,j,ntags,comment
      integer, dimension(:), allocatable :: limit,line
      character(len=str_long) :: buffer
      character(len=str_long), dimension(:), allocatable :: file
      character(len=str_long) :: val
      character(len=str_medium) :: tag
      
      ! Store the root and file name
      this%amRoot=amRoot
      this%filename=trim(adjustl(filename))
      
      ! Empty object
      this%nparams=0

      ! Open the file
      open(newunit=iunit,file=this%filename,form='formatted',status='old',iostat=ierr)
      if (ierr.ne.0) call die('[inputfile constructor] Could not open file: '//trim(this%filename))

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
         call this%add_param(adjustl(trim(tag)),adjustl(trim(val)))
      end do

      ! Deallocate
      deallocate(file,limit,line)
      
      ! Log this info
      call this%log()

   end function constructor
   
   
   !> Add a new parameter:
   !>   - if param has been found already, replace it
   !>   - if param has not been found yet, add it
   subroutine add_param(this,tag,val)
      use messager, only: warn
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in) :: tag
      character(len=*), intent(in) :: val
      type(param_type), dimension(:), allocatable :: temp
      integer :: i
      ! Check if parameter exists already
      if (this%param_exists(tag)) then
         ! Warn the user
         call warn('[inputfile] The following parameter appeared more than once: '//trim(tag))
         call warn('[inputfile] Only the first occurence is considered!')
      else
         ! Resize the parameter array
         allocate(temp(this%nparams+1))
         if (this%nparams.gt.0) temp(1:this%nparams)=this%params
         call move_alloc(temp,this%params)
         ! Add parameter
         this%nparams=this%nparams+1
         this%params(this%nparams)%tag=adjustl(trim(tag))
         this%params(this%nparams)%val=adjustl(trim(val))
         ! Remove unwanted characters from val string
         do i=1,len_trim(this%params(this%nparams)%val)
            if (ichar(this%params(this%nparams)%val(i:i)).eq.9.or.this%params(this%nparams)%val(i:i).eq.','.or.this%params(this%nparams)%val(i:i).eq.';') this%params(this%nparams)%val(i:i)=' '
         end do
      end if
   end subroutine add_param
   
   
   !> Check if a parameter already exists in params array
   function param_exists(this,tag) result(found)
      implicit none
      class(inputfile), intent(in) :: this
      character(len=*), intent(in) :: tag
      logical :: found
      integer :: i
      ! Initializer
      found=.false.
      ! Traverse params array and look for tag
      do i=1,this%nparams
         if (this%params(i)%tag.eq.tag) found=.true.
      end do
   end function param_exists
   
   
   !> Return parameter index for tag
   function param_index(this,tag) result(ind)
      implicit none
      class(inputfile), intent(in) :: this
      character(len=*), intent(in) :: tag
      integer :: i,ind
      ! Initializer
      ind=0
      ! Traverse params array and look for tag
      do i=1,this%nparams
         if (this%params(i)%tag.eq.tag) ind=i
      end do
   end function param_index
   
   
   !> Get size of the parameter value for reading arrays
   function param_getsize(this,tag) result(count)
      use messager, only: die
      implicit none
      class(inputfile), intent(in) :: this
      character(len=*), intent(in) :: tag
      integer :: i,ind,count
      ! Initialize counter
      count=0
      ! Get the parameter index
      ind=this%param_index(tag)
      if (ind.gt.0) then
         ! Count the number of values provided
         if (this%params(ind)%val(1:1).ne.' ') count=1
         do i=2,len_trim(this%params(ind)%val)
            if (this%params(ind)%val(i:i).eq.' '.and.this%params(ind)%val(i-1:i-1).ne.' ') count=count+1
         end do
         return
      end if
      ! Still here: we failed to find a size
      call die('[inputfile param_getsize] Did not find parameter: '//trim(tag))
   end function param_getsize
   
   
   !> Read integer value associated with parameter tag
   subroutine readint(this,tag,val,default)
      use messager, only: die
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)            :: tag
      integer         , intent(out)           :: val
      integer         , intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(this%params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(this%params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[inputfile read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine readint
   
   
   !> Read logical value associated with parameter tag
   subroutine readlogical(this,tag,val,default)
      use messager, only: die
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)            :: tag
      logical         , intent(out)           :: val
      logical         , intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(this%params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(this%params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[inputfile read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine readlogical
   
   
   !> Read real value associated with parameter tag
   subroutine readfloat(this,tag,val,default)
      use messager,  only: die
      use precision, only: WP
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)            :: tag
      real(WP)        , intent(out)           :: val
      real(WP)        , intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(this%params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(this%params(ind)%val,*) val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[inputfile read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine readfloat
   
   
   !> Read character value associated with parameter tag
   subroutine readchar(this,tag,val,default)
      use messager, only: die
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)            :: tag
      character(len=*), intent(out)           :: val
      character(len=*), intent(in) , optional :: default
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         ! Parameter is defined
         if (len_trim(this%params(ind)%val).gt.0) then
            ! Value exists, so read it and return
            read(this%params(ind)%val,'(a)') val
            return
         end if
      end if
      ! If still here, we have not found a value to read
      if (present(default)) then
         val=default
      else
         call die('[inputfile read] Did not find required parameter: '//trim(tag))
      end if
   end subroutine readchar
   
   
   !> Read integer array value associated with parameter tag
   subroutine readintarray(this,tag,val)
      use messager, only: die
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)            :: tag
      integer, dimension(:), intent(out)      :: val
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(this%params(ind)%val,*) val
         return
      end if
      ! If still here, we have not found a value to read
      call die('[inputfile read] Did not find required parameter: '//trim(tag))
   end subroutine readintarray
   
   
   !> Read float array value associated with parameter tag
   subroutine readfloatarray(this,tag,val)
      use messager,  only: die
      use precision, only: WP
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)        :: tag
      real(WP), dimension(:), intent(out) :: val
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(this%params(ind)%val,*) val
         return
      end if
      ! If still here, we have not found a value to read
      call die('[inputfile read] Did not find required parameter: '//trim(tag))
   end subroutine readfloatarray
   
   
   !> Read float 2D array value associated with parameter tag
   subroutine readfloatarray2D(this,tag,val)
      use messager,  only: die
      use precision, only: WP
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)          :: tag
      real(WP), dimension(:,:), intent(out) :: val
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(this%params(ind)%val,*) val
         return
      end if
      ! If still here, we have not found a value to read
      call die('[inputfile read] Did not find required parameter: '//trim(tag))
   end subroutine readfloatarray2D
   
   
   !> Read character array value associated with parameter tag
   subroutine readchararray(this,tag,val)
      use messager, only: die
      implicit none
      class(inputfile), intent(inout) :: this
      character(len=*), intent(in)                :: tag
      character(len=*), dimension(:), intent(out) :: val
      integer :: ind
      ! Get the parameter index
      ind=this%param_index(tag)
      ! Read as appropriate
      if (ind.gt.0) then
         read(this%params(ind)%val,*) val
         return
      end if
      ! If still here, we have not found a value to read
      call die('[inputfile read] Did not find required parameter: '//trim(tag))
   end subroutine readchararray


   !> Print out parameters
   subroutine print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(inputfile), intent(in) :: this
      integer :: i
      integer, parameter :: colwidth=30
      character(len=colwidth) :: tag,val,sep
      if (this%amRoot) then
         write(output_unit,'("Input file [",a,"] content:")') trim(this%filename)
         do i=1,colwidth
            sep(i:i)='_'
         end do
         write(output_unit,'(1x,"_",a,"___",a,"_",1x)') sep,sep
         tag='Parameter tag name'; val='Assigned value'
         write(output_unit,'("| ",a," | ",a," |")') tag,val
         do i=1,this%nparams
            ! Check for truncation of strings
            if (len_trim(this%params(i)%tag).gt.colwidth) then
               tag=this%params(i)%tag(1:colwidth-3)//'...'
            else
               tag=''; tag=this%params(i)%tag(1:len_trim(this%params(i)%tag))
            end if
            if (len_trim(this%params(i)%val).gt.colwidth) then
               val=this%params(i)%val(1:colwidth-3)//'...'
            else
               val=''; val=this%params(i)%val(1:len_trim(this%params(i)%val))
            end if
            write(output_unit,'("| ",a," | ",a," |")') tag,val
         end do
         write(output_unit,'("|_",a,"_|_",a,"_|")') sep,sep
      end if
   end subroutine print


   !> Log parameters
   subroutine inputlog(this)
      use messager, only: log
      implicit none
      class(inputfile), intent(in) :: this
      integer :: i
      integer, parameter :: colwidth=30
      character(len=colwidth) :: tag,val,sep
      character(len=str_long) :: message
      if (this%amRoot) then
         write(message,'("Input file [",a,"] content:")') trim(this%filename); call log(message)
         do i=1,colwidth
            sep(i:i)='_'
         end do
         write(message,'(1x,"_",a,"___",a,"_",1x)') sep,sep; call log(message)
         tag='Parameter tag name'; val='Assigned value'
         write(message,'("| ",a," | ",a," |")') tag,val; call log(message)
         do i=1,this%nparams
            ! Check for truncation of strings
            if (len_trim(this%params(i)%tag).gt.colwidth) then
               tag=this%params(i)%tag(1:colwidth-3)//'...'
            else
               tag=''; tag=this%params(i)%tag(1:len_trim(this%params(i)%tag))
            end if
            if (len_trim(this%params(i)%val).gt.colwidth) then
               val=this%params(i)%val(1:colwidth-3)//'...'
            else
               val=''; val=this%params(i)%val(1:len_trim(this%params(i)%val))
            end if
            write(message,'("| ",a," | ",a," |")') tag,val; call log(message)
         end do
         write(message,'("|_",a,"_|_",a,"_|")') sep,sep; call log(message)
      end if
   end subroutine inputlog
   
   
end module inputfile_class
