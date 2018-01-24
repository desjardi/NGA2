!> Module handling standard output to the screen
!! or to text files for monitoring purposes.
module monitor
  use precision
  use string
  implicit none
  
  ! Log file
  integer :: logfile
  
  ! Preset some length and formats for the columns
  integer,          parameter :: col_len=14
  character(len=*), parameter :: f1='(a12)'
  character(len=*), parameter :: f2='(i12)'
  character(len=*), parameter :: f3='(ES12.5)'
  
  ! File numbers
  integer :: nfiles  !< Total number of files
  integer :: ifile   !< Current file
  integer :: iscreen !< Latest screen datafile
  
  ! Definition of type 'mfile'
  type :: mfile
     ! Data to dump
     integer :: ncol   !< Number of columns to output
     character(len=str_medium), dimension(:), allocatable :: col_head
     character(len=1),          dimension(:), allocatable :: col_type
     real(WP),                  dimension(:), allocatable :: val
     ! File information
     integer :: iunit
     character(len=str_medium) :: filename
     logical :: isnew  !< True: header needs to be written, False: no need to write header
     ! Output frequency
     integer :: type   !< type=0: monitoring to the screen (only latest definition is used)
                       !! type=1: file dumped once per timestep
                       !! type=2: file dumped once per subiteration     
  end type mfile
  
  ! Array of mfiles
  type(mfile), dimension(:), allocatable :: mfiles
  
contains
  
  
  !> Write the header
  subroutine monitor_first_dump(myfile)
    implicit none
    integer, intent(in) :: myfile
    integer :: ierr,icol,offset,index1
    character(len=col_len)    :: col
    character(len=str_long)   :: header
    character(len=str_medium) :: filename,buffer
    logical :: twolines
    
    ! Create the header
    select case(mfiles(myfile)%type)
    case (0)
       write(header(1+0*col_len:),f1) 'Step'
       write(header(1+1*col_len:),f1) 'Time'
       offset=2*col_len
    case (1)
       write(header(1+0*col_len:),f1) 'Step'
       write(header(1+1*col_len:),f1) 'Time'
       offset=2*col_len
    case(2)
       write(header(1+0*col_len:),f1) 'Step'
       write(header(1+1*col_len:),f1) 'Time'
       write(header(1+2*col_len:),f1) 'Niter'
       offset=3*col_len
    end select
    
    ! Extract the first line, detect if a second line is needed
    twolines=.false.
    do icol=1,mfiles(myfile)%ncol
       read(mfiles(myfile)%header(icol),'(a)') buffer
       index1=index(trim(buffer),' ')
       if (index1.ne.0 .and. index1.lt.col_len-1) then
          twolines=.true.
          read(buffer(1:index1),f1) col
       else
          read(buffer,f1) col
       end if
       write(header(offset+1+(icol-1)*col_len:),f1) trim(col)
    end do
    
    ! Dump the first line of header
    write(mfiles(myfile)%iunit,'(a)') trim(header)
    
    ! Dump second line if needed
    if (twolines) then
       header = ''
       do icol=1,mfiles(myfile)%ncol
          read (mfiles(myfile)%header(icol),'(a)') buffer
          index1=index(trim(buffer),' ')
          if (index1.ne.0.and. index1.lt.col_len-1) then
             read(buffer(index1:),f1) col
          else
             col=''
          end if
          write(header(offset+1+(icol-1)*col_len:),f1) trim(col)
       end do
       write(mfiles(myfile)%iunit,'(a)') trim(header)
    end if
    
    ! File is not new anymore
    mfiles(myfile)%isnew=.false.
    
    return
  end subroutine monitor_first_dump
  
  
  !> Dump the values to the files at each timestep
  subroutine monitor_dump_timestep
    use parallel
    use time_info
    implicit none
    integer :: icol,offset
    character(len=str_long) :: line
    ! Only the root process does something
    if (irank.ne.iroot) return
    ! Loop over all files
    do ifile=1,nfiles
       ! Test if we need to dump the values
       if (mfiles(ifile)%type.eq.1) then
          ! Check if we need to open the file and write the header
          if (mfiles(ifile)%isnew) call monitor_first_dump(ifile)
          ! Start the line to dump with time step info
          write(line(1+0*col_len:),f2) ntime
          write(line(1+1*col_len:),f3)  time
          offset=2*col_len
          ! Add all variables passed to monitor
          do icol=1,mfiles(ifile)%ncol
             select case(mfiles(ifile)%col_type(icol))
             case('i')
                write(line(offset+1+(icol-1)*col_len:),f2) int(mfiles(ifile)%val(icol))
             case('r')
                write(line(offset+1+(icol-1)*col_len:),f3) real(mfiles(ifile)%val(icol),WP)
             end select
          end do
          ! Dump the line
          write(mfiles(ifile)%iunit,'(a)') trim(line)
          ! Flush file
          call flush(mfiles(ifile)%iunit)
       end if
    end do
    return
  end subroutine monitor_dump_timestep
  
  
  !> Dump the values to the files at each subiteration
  subroutine monitor_dump_iteration
    use parallel
    use time_info
    implicit none
    integer :: icol,offset
    character(len=str_long) :: line
    ! Only the root process does something
    if (irank.ne.iroot) return
    ! Loop over all files
    do ifile=1,nfiles
       ! Test if we need to dump the values
       if (mfiles(ifile)%type.eq.2) then
          ! Check if we need to open the file and write the header
          if (mfiles(ifile)%isnew) call monitor_first_dump(ifile)
          ! Start the line to dump with time step and iteration info
          write(line(1+0*col_len:),f2) ntime
          write(line(1+1*col_len:),f3) time
          write(line(1+2*col_len:),f2) niter
          offset=3*col_len
          ! Add all variables passed to monitor
          do icol=1,mfiles(ifile)%ncol
             select case(mfiles(ifile)%col_type(icol))
             case('i')
                write(line(offset+1+(icol-1)*col_len:),f2) int(mfiles(ifile)%val(icol))
             case('r')
                write(line(offset+1+(icol-1)*col_len:),f3) real(mfiles(ifile)%val(icol),WP)
             end select
          end do
          ! Dump the line
          write(mfiles(ifile)%iunit,'(a)') trim(line)
          ! Flush file
          call flush(mfiles(ifile)%iunit)
       end if
    end do
    return
  end subroutine monitor_dump_iteration
  
  
  !> Print some values to the screen
  subroutine monitor_dump_screen
    use parallel
    use time_info
    implicit none
    ! Only the root process does something
    if (irank.ne.iroot) return
    ! Work only on last screen output defined
    if (mfiles(iscreen)%isnew) call monitor_first_dump(iscreen)
    ! Start the line to dump with time step and iteration info
    write(line(1+0*col_len:),f2) ntime
    write(line(1+1*col_len:),f3) time
    offset=2*col_len
    ! Add all variables passed to monitor
    do icol=1,mfiles(iscreen)%ncol
       select case(mfiles(iscreen)%col_type(icol))
       case('i')
          write(line(offset+1+(icol-1)*col_len:),f2) int(mfiles(iscreen)%val(icol))
       case('r')
          write(line(offset+1+(icol-1)*col_len:),f3) real(mfiles(iscreen)%val(icol),WP)
       end select
    end do
    ! Dump the line
    write(mfiles(iscreen)%iunit,'(a)') trim(line)
    ! Flush file
    call flush(mfiles(iscreen)%iunit)
    return
  end subroutine monitor_dump_screen
  
  
end module monitor


!> Initialize the monitor module
subroutine monitor_init
  use monitor
  use parallel
  implicit none
  ! Set number of files to zero
  nfiles=0
  ifiles=0
  iscreen=0
  ! Only the root process does the following
  if (irank.eq.iroot) then
     ! Create the monitor directory
     call CREATE_FOLDER("monitor")
     ! Open log file
     open(newunit=logfile,file="monitor/log",form="formatted",iostat=ierr,status="REPLACE")
  end if
  return
end subroutine monitor_init


!> Create a new file to monitor
subroutine monitor_create_file(filename,ncol,type)
  use monitor
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: ncol,type
  type(mfile), dimension(:), allocatable :: temp
  integer :: ierr
  ! Add a file
  nfiles=nfiles+1
  allocate(temp(nfiles)); temp(1:size(mfiles))=mfiles; call move_alloc(temp,mfiles)
  ! Preset the values
  mfiles(nfiles)%filename=trim(filename)
  mfiles(nfiles)%type    =type
  mfiles(nfiles)%ncol    =ncol
  ! Deal with screen output
  if (type.eq.0) then
     if (iscreen.ne.0) call warn('[monitor_create_file] Screen output is being redefined!')
     iscreen=nfiles
  end if
  ! Allocate the arrays
  allocate(mfiles(nfiles)%header(ncol))
  allocate(mfiles(nfiles)%col_type(ncol))
  allocate(mfiles(nfiles)%val(ncol))
  ! Set file as new
  mfiles(nfiles)%isnew=.true.
  ! Open the file
  if (irank.eq.iroot) then
     if (mfiles(nfiles)%type.ne.0) then
        open(newunit=mfiles(nfiles)%iunit,file='monitor/'//trim(mfiles(nfiles)%filename),form="formatted",iostat=ierr,status="REPLACE")
     else
        mfiles(nfiles)%iunit=6
     end if
  end if
  ! Automatically select the new file
  ifile=nfiles
  return
end subroutine monitor_create_file


! Search for the file and set the allocatable
subroutine monitor_select_file(filename)
  use monitor
  implicit none
  character(len=*), intent(in) :: filename
  ! Find the file
  loop: do ifile=1,nfiles
     if (trim(mfiles(ifile)%filename).eq.trim(filename)) exit loop
  end do loop
  if (ifile.gt.nfiles) call die('[monitor_select_file] Could not file following file: '//trim(filename))
  return
end subroutine monitor_select_file


!> Set the header of the file
subroutine monitor_set_header(icol,header,col_type)
  use monitor
  implicit none
  
  integer, intent(in) :: icol
  character(len=*), intent(in) :: header
  character(len=1), intent(in) :: col_type
  
  ! Test if enough columns
  if (icol.gt.mfiles(ifile)%ncol) then
     print*,'filename: ', trim(mfiles(ifile)%filename)
     print*,'column index:', icol
     call die('[monitor_set_header] Column index too large for file '//trim(mfiles(ifile)%filename))
  else
     mfiles(ifile)%header(icol)=trim(header)
  end if
  
  ! Test if column type is acceptable
  select case (lowercase(col_type))
  case ('i')
     mfiles(ifile)%col_type(icol)='i'
  case ('r')
     mfiles(ifile)%col_type(icol)='r'
  case default
     call die('[monitor_set_header] Column type unknown for file '//trim(mfiles(ifile)%filename))
  end select
  
  return
end subroutine monitor_set_header


!> Set the values to dump in the file
subroutine monitor_set_value(int,val)
  use monitor
  implicit none
  integer,  intent(in) :: int !< Column number for the new value
  real(WP), intent(in) :: val !< Value to dump in file
  if (int.gt.mfiles(ifile)%ncol) then
     call die('[monitor_set_single_value] Column index too large for file '//trim(mfiles(ifile)%filename))
  else
     mfiles(ifile)%val(int) = val
  end if
  return
end subroutine monitor_set_value


!> Set all the values to dump in the file simultaneously
subroutine monitor_set_array_values(val)
  use monitor
  implicit none
  real(WP), dimension(mfiles(ifile)%ncol) :: val
  mfiles(ifile)%val=val
  return
end subroutine monitor_set_array_values


!> Monitor basic NGA actions in a log file
subroutine monitor_log(text)
  use monitor
  use parallel
  implicit none
  character(len=*), intent(in) :: text !< Text to be printed to the log file
  integer, dimension(length=41) :: tstamp
  if (irank.eq.iroot) then
     call timestamp(tstamp)
     write(logfile,'(a41,1x,a2,1x,a60)') tstamp,'=>',text
     call flush(logfile)
  end if
  return
end subroutine monitor_log


!> Finalize monitoring by closing all files
subroutine monitor_final
  use monitor
  use parallel
  implicit none
  do ifile=1,nfiles
     ! Deallocate all arrays
     deallocate(mfiles(ifile)%header);   nullify(mfiles(ifile)%header)
     deallocate(mfiles(ifile)%col_type); nullify(mfiles(ifile)%col_type)
     deallocate(mfiles(ifile)%val);      nullify(mfiles(ifile)%val)
     ! Close files
     if (irank.eq.iroot.and.mfiles(ifile)%type.ne.0) close(mfiles(ifile)%iunit)
  end do
  ! Close logfiles
  if (irank.eq.iroot) close(logfile)
  ! Deallocate mfiles
  deallocate(mfiles); nullify(mfiles)
  return
end subroutine monitor_final


!> This routine outputs the current YMDHS date in a readeable fashion.
subroutine timestamp(tstamp)
  implicit none
  character(len=*), intent(out) :: tstamp !< This is the timestamp returned by the routine - should be 41 long
  character(len=8) :: ampm
  integer :: d,h,m,mm,n,s,y
  integer, dimension(8) :: values
  character(len=9), parameter, dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)
  
  ! Get date and time and make output more readable
  call date_and_time (values)
  y=values(1); m=values(2); d=values(3)
  h=values(5); n=values(6); s=values(7)
  mm=values(8)
  
  ! Assign AM/PM to time
  if (h.lt.12) then
     ampm='AM'
  elseif (h.eq.12) then
     if (n.eq.0.and.s.eq.0) then
        ampm='Noon'
     else
        ampm='PM'
     end if
  else
     h=h-12
     if (h.lt.12) then
        ampm='PM'
     elseif (h.eq.12) then
        if (n.eq.0.and.s.eq.0) then
           ampm='Midnight'
        else
           ampm='AM'
        end if
     end if
  end if
  
  ! Generate output
  write(tstamp,'(i2,1x,a9,1x,i4,1x,a1,1x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a8)') &
       d,month(m),y,,'@',h,':',n,':',s,'.',mm,ampm
  
  return
end subroutine timestamp
