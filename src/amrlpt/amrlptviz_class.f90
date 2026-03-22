!> AMR particle visualization handler
!> Mirrors amrviz_class structure: registration-based, time-series aware,
!> restart-capable. Outputs AMReX native particle plotfiles (VisIt/ParaView).
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

   ! -----------------------------------------------------------------------
   ! C interface -- WritePlotFile wrapper in amrlpt_wrapper.cpp
   ! -----------------------------------------------------------------------
   interface
      subroutine amrlpt_write_plotfile(pc, basedir, pname, write_real, write_int) bind(c)
         import
         type(c_ptr), value :: pc
         character(kind=c_char) :: basedir(*), pname(*)
         integer(c_int), intent(in) :: write_real(*), write_int(*)
      end subroutine
   end interface

   ! -----------------------------------------------------------------------
   ! Component index map: name → which rdata or idata indices to toggle.
   ! Order matches the part struct rdata layout (0-based C indices).
   ! -----------------------------------------------------------------------
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

   ! -----------------------------------------------------------------------
   !> Initialize: create output directory, scan for existing files on restart
   ! -----------------------------------------------------------------------
   subroutine initialize(this, lpt, name)
      use filesys,  only: makedir, isdir, isfile
      use parallel, only: MPI_REAL_WP, amRoot
      use mpi_f08,  only: MPI_BCAST, MPI_INTEGER
      implicit none
      class(amrlptviz), intent(inout) :: this
      class(amrlpt), target, intent(in) :: lpt
      character(len=*), intent(in) :: name

      character(len=str_long) :: timefile
      integer :: iunit, ierr, n
      real(WP) :: t
      real(WP), allocatable :: tmp(:)

      this%lpt  => lpt
      this%name = trim(adjustl(name))
      this%ntime = 0
      this%write_real = 1
      this%write_int  = 1

      ! Create output directory
      if (lpt%amRoot) then
         if (.not.isdir('amrviz')) call makedir('amrviz')
         if (.not.isdir('amrviz/'//trim(this%name))) &
            call makedir('amrviz/'//trim(this%name))
      end if

      ! Look for existing time index file (written by us on previous runs)
      timefile = 'amrviz/'//trim(this%name)//'/particle_times.txt'

      if (lpt%amRoot .and. isfile(trim(timefile))) then
         open(newunit=iunit, file=trim(timefile), status='old', action='read', iostat=ierr)
         if (ierr == 0) then
            ! Count lines
            n = 0
            do
               read(iunit, *, iostat=ierr)
               if (ierr /= 0) exit
               n = n + 1
            end do
            rewind(iunit)
            if (n > 0) then
               allocate(tmp(n))
               do n = 1, size(tmp)
                  read(iunit,*) tmp(n)
               end do
               this%ntime = size(tmp)
               call move_alloc(tmp, this%time)
            end if
            close(iunit)
         end if
      end if

      ! Broadcast ntime and time array
      call MPI_BCAST(this%ntime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (this%ntime > 0) then
         if (.not.lpt%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time, this%ntime, MPI_REAL_WP, 0, MPI_COMM_WORLD, ierr)
      end if

   end subroutine initialize

   ! -----------------------------------------------------------------------
   !> Toggle a named field group on or off.
   !> Names: 'd', 'vel', 'angVel', 'Acol', 'Tcol', 'dt', 'flag'
   !> You can also pass individual component names: 'vx','vy','vz' etc.
   ! -----------------------------------------------------------------------
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

   ! -----------------------------------------------------------------------
   !> Write one particle plotfile snapshot.
   !> Directory: amrviz/<name>/plt<ntime> (6-digit zero-padded)
   !> Sub-directory inside: 'particles' (AMReX particle name)
   ! -----------------------------------------------------------------------
   subroutine write(this, time)
      use iso_c_binding, only: c_null_char
      use parallel,      only: MPI_REAL_WP
      use mpi_f08,       only: MPI_BCAST
      implicit none
      class(amrlptviz), intent(inout) :: this
      real(WP), intent(in) :: time

      character(len=str_long) :: pltdir, timefile
      real(WP), allocatable :: tmp(:)
      integer :: iunit, ierr, i, n
      logical :: rewind_flag

      ! --------------- Update time array (same rewind logic as amrviz) -----
      if (this%ntime == 0) then
         this%ntime = 1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(1))
         this%time(1) = time
      else
         n = 1
         do i = this%ntime, 1, -1
            if (this%time(i) < time - 1.0e-6_WP) then
               n = i + 1; exit
            end if
         end do
         this%ntime = n
         allocate(tmp(this%ntime))
         tmp = [this%time(1:this%ntime-1), time]
         call move_alloc(tmp, this%time)
      end if

      ! --------------- Construct output directory ---------------------------
      write(pltdir,'(a,"/amrviz/",a,"/plt",i6.6)') &
           '', trim(this%name), this%ntime
      pltdir = adjustl(pltdir)
      ! Strip leading blank from write format
      pltdir = 'amrviz/'//trim(this%name)//'/plt'
      write(pltdir(len_trim(pltdir)+1:len_trim(pltdir)+6),'(i6.6)') this%ntime

      ! --------------- Write via C++ wrapper --------------------------------
      call amrlpt_write_plotfile(this%lpt%pc, &
           trim(pltdir)//c_null_char, 'particles'//c_null_char, &
           this%write_real, this%write_int)

      ! --------------- Update time index file (root only) ------------------
      if (this%lpt%amRoot) then
         timefile = 'amrviz/'//trim(this%name)//'/particle_times.txt'
         open(newunit=iunit, file=trim(timefile), status='replace', action='write')
         do i = 1, this%ntime
            write(iunit,'(ES23.16)') this%time(i)
         end do
         close(iunit)
      end if

   end subroutine write

   ! -----------------------------------------------------------------------
   !> Finalize: clean up allocations
   ! -----------------------------------------------------------------------
   subroutine finalize(this)
      implicit none
      class(amrlptviz), intent(inout) :: this
      if (allocated(this%time)) deallocate(this%time)
      nullify(this%lpt)
      this%ntime = 0
   end subroutine finalize

end module amrlptviz_class
