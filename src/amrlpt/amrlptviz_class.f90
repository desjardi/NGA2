!> AMR particle visualization handler
!>
!> Usage:
!>   type(amrlptviz) :: pviz
!>   call pviz%initialize(lpt, name='particles')
!>   call pviz%select_comp('vel',   on=.true. )   ! write velocity
!>   call pviz%select_comp('Acol',  on=.false.)   ! skip collision accels
!>   ! In time loop:
!>   call pviz%write(time)
module amrlptviz_class
   use precision,   only: WP
   use string,      only: str_medium, str_long
   use iso_c_binding
   use amrlpt_class, only: amrlpt, AMRLPT_NREAL, AMRLPT_NINT
   implicit none
   private

   public :: amrlptviz

   ! C interface -- WritePlotFile wrapper in amrlpt_wrapper.cpp
   interface
      subroutine amrlpt_write_plotfile(pc, basedir, pname, write_real, write_int, time) bind(c)
         import
         type(c_ptr), value :: pc
         character(kind=c_char) :: basedir(*), pname(*)
         integer(c_int), intent(in) :: write_real(*), write_int(*)
         real(c_double), value :: time
      end subroutine
      function amrlpt_read_plotfile_time(basedir) result(t) bind(c)
         import
         character(kind=c_char) :: basedir(*)
         real(c_double) :: t
      end function
   end interface

   ! Component index map: name -> which rdata or idata indices to toggle.
   ! Order matches the part struct rdata layout (0-based C indices).
   ! rdata[0]    = d
   ! rdata[1..3] = vel (vx,vy,vz)
   ! rdata[4..6] = angVel (wx,wy,wz)
   ! rdata[7..9] = Acol (ax,ay,az)
   ! rdata[10..12]= Tcol (tx,ty,tz)
   ! rdata[13]   = dt
   ! idata[0]    = flag

   type :: amrlptviz

      !> Associated particle solver (non-owning pointer)
      class(amrlpt), pointer :: lpt => null()

      !> Output subdirectory name (inside amrviz/)
      character(len=str_medium) :: name = 'UNNAMED_AMRLPTVIZ'

      !> Time-series tracking (matches amrviz pattern)
      integer  :: ntime = 0
      real(WP), allocatable :: time(:)

      !> Component bitmasks: 1=write, 0=skip
      integer(c_int) :: write_real(AMRLPT_NREAL) = 1   !< all reals on by default
      integer(c_int) :: write_int (AMRLPT_NINT)  = 1   !< all ints  on by default

   contains
      procedure :: initialize    !< Set up output dir, restore time series on restart
      procedure :: select_comp   !< Toggle a named field group on/off
      procedure :: write         !< Write one plotfile snapshot
      procedure :: finalize      !< Clean up
   end type amrlptviz

contains

   !> Initialize: create output directory, restore time series on restart.
   !> Reads all existing plotfile times from disk (sequentially numbered).
   !> Rewind/truncation based on the restart time happens in write().
   subroutine initialize(this, lpt, name)
      use filesys,  only: makedir, isdir
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_BCAST, MPI_INTEGER
      implicit none
      class(amrlptviz), intent(inout) :: this
      class(amrlpt), target, intent(in) :: lpt
      character(len=*), intent(in) :: name

      character(len=str_long) :: pltdir
      integer :: ierr, n
      real(c_double) :: file_time

      this%lpt  => lpt
      this%name = trim(adjustl(name))
      this%ntime = 0
      this%write_real = 0          ! default: position (free) + diameter + velocity
      this%write_real(1) = 1       ! d
      this%write_real(2:4) = 1     ! vx, vy, vz
      this%write_int  = 0

      ! Create output directory
      if (lpt%amr%amRoot) then
         if (.not.isdir('amrviz')) call makedir('amrviz')
         if (.not.isdir('amrviz/'//trim(this%name))) &
            call makedir('amrviz/'//trim(this%name))
      end if

      ! Root probes for existing plotfiles via C++ function (mirrors amrviz pattern).
      ! amrlpt_read_plotfile_time returns -1.0 when directory/time file is absent.
      if (lpt%amr%amRoot) then
         n = 0
         find_files: do
            n = n + 1
            write(pltdir,'("amrviz/",a,"/plt.part.",i6.6)') trim(this%name), n
            file_time = amrlpt_read_plotfile_time(trim(pltdir)//c_null_char)
            if (file_time.lt.0.0_c_double) exit find_files
         end do find_files
         this%ntime = n - 1
         if (this%ntime.gt.0) then
            allocate(this%time(this%ntime))
            do n = 1, this%ntime
                write(pltdir,'("amrviz/",a,"/plt.part.",i6.6)') trim(this%name), n
               this%time(n) = real(amrlpt_read_plotfile_time(trim(pltdir)//c_null_char), WP)
            end do
         end if
      end if

      ! Broadcast ntime and time array to all ranks
      call MPI_BCAST(this%ntime, 1, MPI_INTEGER, 0, lpt%amr%comm, ierr)
      if (this%ntime.gt.0) then
         if (.not.lpt%amr%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time, this%ntime, MPI_REAL_WP, 0, lpt%amr%comm, ierr)
      end if

   end subroutine initialize

   !> Toggle a named field group on or off.
   !> Names: 'd', 'vel', 'angVel', 'Acol', 'Tcol', 'dt', 'flag'
   !> You can also pass individual component names: 'vx','vy','vz' etc.
   subroutine select_comp(this, name, on)
      implicit none
      class(amrlptviz), intent(inout) :: this
      character(len=*), intent(in) :: name
      logical, intent(in) :: on
      integer :: val

      val = merge(1, 0, on)

      select case(trim(name))
      ! Real components (1-based Fortran indexing into write_real)
      case('d')
         this%write_real(1) = val
      case('vel')
         this%write_real(2:4) = val
      case('vx')
         this%write_real(2) = val
      case('vy')
         this%write_real(3) = val
      case('vz')
         this%write_real(4) = val
      case('angVel')
         this%write_real(5:7) = val
      case('wx')
         this%write_real(5) = val
      case('wy')
         this%write_real(6) = val
      case('wz')
         this%write_real(7) = val
      case('Acol')
         this%write_real(8:10) = val
      case('ax')
         this%write_real(8) = val
      case('ay')
         this%write_real(9) = val
      case('az')
         this%write_real(10) = val
      case('Tcol')
         this%write_real(11:13) = val
      case('tx')
         this%write_real(11) = val
      case('ty')
         this%write_real(12) = val
      case('tz')
         this%write_real(13) = val
      case('dt')
         this%write_real(14) = val
      ! Integer components
      case('flag')
         this%write_int(1) = val
      ! Convenience groups
      case('all')
         this%write_real = val
         this%write_int  = val
      case('position_only')
         this%write_real = 0
         this%write_int  = 0
      end select

   end subroutine select_comp

   !> Write one particle plotfile snapshot.
   !> Directory: amrviz/<name>/plt<ntime> (8-digit zero-padded)
   subroutine write(this, time)
      use iso_c_binding, only: c_null_char, c_double
      implicit none
      class(amrlptviz), intent(inout) :: this
      real(WP), intent(in) :: time

      character(len=str_long) :: pltdir
      real(WP), allocatable :: tmp(:)
      integer :: i, n

      ! Update time array (same rewind logic as amrviz)
      if (this%ntime.eq.0) then
         this%ntime = 1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(1))
         this%time(1) = time
      else
         n = 1
         rewind: do i = this%ntime, 1, -1
            if (this%time(i).lt.time - 1.0e-6_WP) then
               n = i + 1; exit rewind
            end if
         end do rewind
         this%ntime = n
         allocate(tmp(this%ntime))
         tmp = [this%time(1:this%ntime-1), time]
         call move_alloc(tmp, this%time)
      end if

      ! Construct output directory
      write(pltdir,'("amrviz/",a,"/plt.part.",i6.6)') trim(this%name), this%ntime

      ! Write via C++ wrapper
      call amrlpt_write_plotfile(this%lpt%pc, trim(pltdir)//c_null_char, 'particles'//c_null_char, this%write_real, this%write_int, real(time, c_double))

      ! Write/rewrite JSON .series file for ParaView time association
      if (this%lpt%amr%amRoot) then
         open(newunit=n, file='amrviz/'//trim(this%name)//'/plt.part.series', status='replace', action='write')
         write(n,'(a)') '{ "file-series-version": "1.0",'
         write(n,'(a)') '  "files": ['
         do i = 1, this%ntime
            write(pltdir,'("plt.part.",i6.6)') i
            if (i.lt.this%ntime) then
               write(n,'(4x,a,a,a,es17.10,a)') '{ "name": "', trim(pltdir), '", "time": ', this%time(i), ' },'
            else
               write(n,'(4x,a,a,a,es17.10,a)') '{ "name": "', trim(pltdir), '", "time": ', this%time(i), ' }'
            end if
         end do
         write(n,'(a)') '  ]'
         write(n,'(a)') '}'
         close(n)
      end if

   end subroutine write

   !> Finalize: clean up allocations
   subroutine finalize(this)
      implicit none
      class(amrlptviz), intent(inout) :: this
      if (allocated(this%time)) deallocate(this%time)
      nullify(this%lpt)
      this%ntime = 0
   end subroutine finalize

end module amrlptviz_class
