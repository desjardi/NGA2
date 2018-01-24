!> This module provides an input file parser capable
!! of handling various data type as well as arrays.
module parser
  use precision
  use string
  implicit none
  
  ! Input values storage
  integer :: nfields !< This is the number of fields in the input file
  type entry_type
     character(str_medium) :: tag
     character(str_long)   :: value
  end type entry_type
  type(entry_type), dimension(:), pointer :: entries !< This array stores the field/value entries from the file
  
  ! Define interface for reading
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
  subroutine parser_newentry(mytag,myvalue)
    implicit none
    character(*), intent(in) :: mytag
    character(*), intent(in) :: myvalue
    integer :: ifield
    logical :: isdef 
    type(entry_type), dimension(:), pointer :: temp
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) then
       ! Add a field
       nfields=nfields+1
       allocate(temp(nfields))
       temp(1:size(entries))=entries
       call move_alloc(temp,entries)
       ! Set value
       entries(nfields)%tag=mytag
       entries(nfields)%value=myvalue
    else
       entries(ifield)%value=myvalue
    end if
    
    return
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
    
    return
  end subroutine parser_fieldfortag
  
  
  !> Check whether the field is defined
  subroutine parser_is_defined(mytag,isdef)
    implicit none
    
    character(*), intent(in) :: mytag
    logical, intent(out) :: isdef
    integer :: ifield
    
    call parser_fieldfortag(mytag,ifield,isdef)
    
    return
  end subroutine parser_is_defined
  
  
  !> Read logicals
  subroutine parser_readlogical(mytag,value,default)
    implicit none
    
    character(*), intent(in) :: mytag
    logical, intent(out) :: value
    logical, optional, intent(in) :: default
    integer :: ifield
    logical :: isdef
    integer :: conv
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       if (len_trim(entries(ifield)%value).eq.1) then          
          read(entries(ifield)%value,*) conv
          value = (conv.eq.1)
       else
          read(entries(ifield)%value,*) value
       end if
    end if
    
    return
  end subroutine parser_readlogical
  
  
  !> Read integers
  subroutine parser_readint(mytag,value,default)
    implicit none
    
    character(*), intent(in) :: mytag
    integer, intent(out) :: value
    integer, optional, intent(in) :: default
    integer :: ifield
    logical :: isdef

    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       read(entries(ifield)%value,*) value
    end if
    
    return
  end subroutine parser_readint
  
  
  !> Read floats
  subroutine parser_readfloat(mytag,value,default)
    implicit none
    
    character(*), intent(in) :: mytag
    real(WP), intent(out) :: value
    real(WP), optional, intent(in) :: default
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       read(entries(ifield)%value,*) value
    end if
    
    return
  end subroutine parser_readfloat
  
  
  !> Read characters
  subroutine parser_readchar(mytag,value,default)
    implicit none
    
    character(*), intent(in) :: mytag
    character(len=*), intent(out) :: value
    character(len=*), optional, intent(in) :: default
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       call die('[parser] Did not find tag: '//trim(mytag))
    else
       read(entries(ifield)%value,'(a)') value
    end if
    
    return
  end subroutine parser_readchar
  
  
  !> Count size of the arrays
  subroutine parser_getsize(mytag,numb)
    implicit none
    
    character(*), intent(in) :: mytag
    integer, intent(out) :: numb
    integer :: ifield
    logical :: isdef
    integer :: i
    integer, dimension(str_long) :: counter
    
    ! Read the field
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    
    ! Count the number of entries                                   
    counter=0
    do i=1,len_trim(entries(ifield)%value)
       if (entries(ifield)%value(i:i).eq.' ') counter(i)=1
    end do
    do i=1+1,len_trim(entries(ifield)%value)
       if (counter(i).eq.1.and.counter(i-1).eq.1) counter(i-1)=0
    end do
    numb=sum(counter)+1
    
    return
  end subroutine parser_getsize
  
  
  !> Read integer arrays
  subroutine parser_readintarray(mytag,value)
    implicit none
    
    character(*), intent(in) :: mytag
    integer, dimension(:), intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%value,*) value
    
    return
  end subroutine parser_readintarray
  
  
  !> Read float arrays
  subroutine parser_readfloatarray(mytag,value)
    implicit none

    character(*), intent(in) :: mytag
    real(WP), dimension(:), intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%value,*) value

    return
  end subroutine parser_readfloatarray
  
  
  !> Read float arrays
  subroutine parser_readfloatarray2D(mytag,value)
    implicit none
    
    character(*), intent(in) :: mytag
    real(WP), dimension(:,:), intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%value,*) value
    
    return
  end subroutine parser_readfloatarray2D
  
  
  !> Read float arrays
  subroutine parser_readchararray(mytag,value)
    implicit none

    character(*), intent(in) :: mytag
    character(*), dimension(:), intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef) call die('[parser] Did not find tag: '//trim(mytag))
    read(entries(ifield)%value,*) value

    return
  end subroutine parser_readchararray
  
  
end module parser


!> This subroutine initializes the parser by resetting the
!! number of fields to zero and emptying the entries array
subroutine parser_init
  use parser
  implicit none
  nfields=0
  if (associated(entries)) then
     deallocate(entries)
     nullify(entries)
  end if
  return
end subroutine parser_init


!> This subroutine empties the parser database
subroutine parser_final
  use parser
  implicit none
  nfields=0
  if (associated(entries)) then
     deallocate(entries)
     nullify(entries)
  end if
  return
end subroutine parser_final


!> Read & parse the input file
subroutine parser_parsefile(input)
  use parser
  implicit none
  integer :: iunit,ierr,limiter,nlines,i,j,ntags,comment
  integer, dimension(:), allocatable :: limit,line
  character(len=str_long) :: buffer
  character(len=str_long), dimension(:), allocatable :: file
  character(len=str_long) :: value
  character(len=str_medium) :: tag
  character(len=*) :: input !< This is the name of the file to be read and parsed.
  
  ! Open the file
  open(newunit=iunit,file=input,form='formatted',status='old',iostat=ierr)
  if (ierr.ne.0) call warn('[parser_parsefile] Could not open input file: '//trim(input))
  
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
     do j=1,line_length
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
     read(buffer(limit(i)+1:),'(a)') value
     if (len_trim(value).ne.0) then
        value=adjustl(value)
        call parser_newentry(tag,value)
     end if
  end do
  
  ! Deallocate
  deallocate(file,limit,line)
  
  return
end subroutine parser_parsefile
