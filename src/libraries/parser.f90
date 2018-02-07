!> This module provides an input file parser capable
!> of handling various data type as well as arrays.
module parser
  use string, only: str_medium,str_long
  implicit none
  private

  public :: parser_init,parser_final
  public :: parser_read!,parser_print
  public :: parser_is_defined
  
  !> Type for storing file parser entries
  type entry_type
     character(str_medium) :: tag !< Tag name for field
     character(str_long)   :: val !< Associated value
  end type entry_type
  
  !> This array stores the field/value entries from the file
  type(entry_type), dimension(:), allocatable :: entries
  integer :: nfields
  
  !> Obtain information about a specific tag
  interface parser_read
     module procedure parser_readlogical
     module procedure parser_readint
     module procedure parser_readintarray
     module procedure parser_readfloat
     module procedure parser_readfloatarray
     module procedure parser_readfloatarray2D
     module procedure parser_readchar
     module procedure parser_readchararray
  end interface parser_read
  
contains
  
  !> Add a new entry in the structure
  subroutine parser_newentry(mytag,myval)
    implicit none
    character(len=*), intent(in) :: mytag
    character(len=*), intent(in) :: myval
    integer :: ifield
    logical :: isdef
    type(entry_type), dimension(:), allocatable :: temp
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) then
       ! Add a field
       nfields=nfields+1
       allocate(temp(nfields))
       if (nfields.gt.1) temp(1:nfields-1)=entries
       call move_alloc(temp,entries)
       ! Set value
       entries(nfields)%tag=mytag
       entries(nfields)%val=myval
    else
       entries(ifield)%val=myval
    end if
  end subroutine parser_newentry
  
  
  !> Get field number from tag
  subroutine parser_fieldfortag(mytag,myfield,isdef)
    implicit none
    character(*), intent(in) :: mytag
    integer, intent(out) :: myfield
    logical, optional, intent(out) :: isdef
    integer :: ifield
    isdef=.false.
    do ifield=1,nfields
       if (entries(ifield)%tag.eq.mytag) then
          myfield=ifield
          isdef=.true.
          return
       end if
    end do
  end subroutine parser_fieldfortag
  
  
  !> Check whether the field is defined
  subroutine parser_is_defined(mytag,isdef)
    implicit none
    character(*), intent(in) :: mytag
    logical, intent(out) :: isdef
    integer :: ifield
    call parser_fieldfortag(mytag,ifield,isdef)
  end subroutine parser_is_defined
  
  
  !> Read logicals
  subroutine parser_readlogical(mytag,val,default)
    use errhandle, only: die
    implicit none
    character(*), intent(in) :: mytag
    logical, intent(out) :: val
    logical, optional, intent(in) :: default
    integer :: ifield,conv
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef.and.present(default)) then   
       val=default
    else if (.not.isdef.and..not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       if (len_trim(entries(ifield)%val).eq.1) then          
          read(entries(ifield)%val,*) conv
          val=(conv.eq.1)
       else
          read(entries(ifield)%val,*) val
       end if
    end if
  end subroutine parser_readlogical
  
  
  !> Read integers
  subroutine parser_readint(mytag,val,default)
    use errhandle, only: die
    implicit none
    character(*), intent(in) :: mytag
    integer, intent(out) :: val
    integer, optional, intent(in) :: default
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef.and.present(default)) then   
       val=default
    else if (.not.isdef.and..not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       read(entries(ifield)%val,*) val
    end if
  end subroutine parser_readint
  
  
  !> Read floats
  subroutine parser_readfloat(mytag,val,default)
    use errhandle, only: die
    use precision, only: WP
    implicit none
    character(*), intent(in) :: mytag
    real(WP), intent(out) :: val
    real(WP), optional, intent(in) :: default
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef.and.present(default)) then   
       val=default
    else if (.not.isdef.and..not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       read(entries(ifield)%val,*) val
    end if
  end subroutine parser_readfloat
  
  
  !> Read characters
  subroutine parser_readchar(mytag,val,default)
    use errhandle, only: die
    implicit none
    character(*), intent(in) :: mytag
    character(len=*), intent(out) :: val
    character(len=*), optional, intent(in) :: default
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef.and.present(default)) then   
       val=default
    else if (.not.isdef.and..not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       read(entries(ifield)%val,'(a)') val
    end if
  end subroutine parser_readchar
  
  
  !> Count size of the arrays
  subroutine parser_getsize(mytag,numb)
    use errhandle, only: die
    implicit none
    character(*), intent(in) :: mytag
    integer, intent(out) :: numb
    integer :: ifield,i
    logical :: isdef
    integer, dimension(str_long) :: counter
    ! Read the field
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    ! Count the number of entries                                   
    counter=0
    do i=1,len_trim(entries(ifield)%val)
       if (entries(ifield)%val(i:i).eq.' ') counter(i)=1
    end do
    do i=1+1,len_trim(entries(ifield)%val)
       if (counter(i).eq.1.and.counter(i-1).eq.1) counter(i-1)=0
    end do
    numb=sum(counter)+1
  end subroutine parser_getsize
  
  
  !> Read integer arrays
  subroutine parser_readintarray(mytag,val)
    use errhandle, only: die
    implicit none
    character(*), intent(in) :: mytag
    integer, dimension(:), intent(out) :: val
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%val,*) val
  end subroutine parser_readintarray
  
  
  !> Read float arrays
  subroutine parser_readfloatarray(mytag,val)
    use errhandle, only: die
    use precision, only: WP
    implicit none
    character(*), intent(in) :: mytag
    real(WP), dimension(:), intent(out) :: val
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%val,*) val
  end subroutine parser_readfloatarray
  
  
  !> Read float arrays
  subroutine parser_readfloatarray2D(mytag,val)
    use errhandle, only: die
    use precision, only: WP
    implicit none
    character(*), intent(in) :: mytag
    real(WP), dimension(:,:), intent(out) :: val
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%val,*) val
  end subroutine parser_readfloatarray2D
  
  
  !> Read float arrays
  subroutine parser_readchararray(mytag,val)
    use errhandle, only: die
    implicit none
    character(*), intent(in) :: mytag
    character(*), dimension(:), intent(out) :: val
    integer :: ifield
    logical :: isdef
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%val,*) val
  end subroutine parser_readchararray
  
  
  !> This subroutine initializes the parser by resetting the
  !> number of fields to zero and emptying the entries array
  subroutine parser_init
    implicit none
    nfields=0
    if (allocated(entries)) deallocate(entries)
  end subroutine parser_init
  
  
  !> This subroutine empties the parser database
  subroutine parser_final
    implicit none
    nfields=0
    if (allocated(entries)) deallocate(entries)
  end subroutine parser_final
  
  
  !> Read & parse the input file
  subroutine parser_parsefile(input)
    use errhandle, only: die
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
    if (ierr.ne.0) call die('[parser_parsefile] Could not open file: '//trim(input))
    
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
       ! Remove the tabs
       do j=1,str_long
          if (ichar(buffer(j:j)).eq.9) buffer(j:j)=' '
       end do
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
       limiter=index(file(i),':')
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
       if (len_trim(val).ne.0) then
          val=adjustl(val)
          call parser_newentry(tag,val)
       end if
    end do
    
    ! Deallocate
    deallocate(file,limit,line)
    
  end subroutine parser_parsefile

  
end module parser
