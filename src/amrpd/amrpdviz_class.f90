!> AMR peridynamics particle visualization handler
!>
!> Writes the particle container to AMReX plotfiles for ParaView/VisIt viewing.
!> Bond visualization is intentionally NOT supported yet — we expect particle-
!> only views to be sufficient for a long while; bonds are a heavy and rarely
!> needed view, and can be added as a parallel writer later.
!>
!> Usage:
!>   type(amrpdviz) :: pviz
!>   call pviz%initialize(pd, name='particles')
!>   call pviz%select_comp('vel',    on=.true. )   ! write velocity
!>   call pviz%select_comp('F_bond', on=.true. )   ! write internal force
!>   ! In time loop:
!>   call pviz%write(time)
module amrpdviz_class
   use precision,    only: WP
   use string,       only: str_medium,str_long
   use iso_c_binding
   use amrpd_class,  only: amrpd,AMRPD_NREAL_PART,AMRPD_NINT_PART
   implicit none
   private

   public :: amrpdviz

   ! C interface -- WritePlotFile wrapper in amrpd_wrapper.cpp
   interface
      subroutine amrpd_write_plotfile(pc,basedir,pname,write_real,write_int,time) bind(c)
         import
         type(c_ptr), value :: pc
         character(kind=c_char) :: basedir(*),pname(*)
         integer(c_int), intent(in) :: write_real(*),write_int(*)
         real(c_double), value :: time
      end subroutine
      function amrpd_read_plotfile_time(basedir) result(t) bind(c)
         import
         character(kind=c_char) :: basedir(*)
         real(c_double) :: t
      end function
   end interface

   ! Component layout (0-based C indices, matching part struct in amrpd_class):
   !   rdata[0..2]   vel       -> vx, vy, vz
   !   rdata[3..5]   F_bond    -> fbx, fby, fbz
   !   rdata[6..8]   F_fluid   -> ffx, ffy, ffz
   !   rdata[9]      mw
   !   rdata[10]     dil
   !   rdata[11]     damage
   !   rdata[12]     nb0
   !   idata[0]      flag

   type :: amrpdviz

      !> Associated peridynamics solver (non-owning pointer)
      class(amrpd), pointer :: pd => null()

      !> Output subdirectory name (inside amrviz/)
      character(len=str_medium) :: name='UNNAMED_AMRPDVIZ'

      !> Time-series tracking (matches amrviz / amrlptviz pattern)
      integer :: ntime=0
      real(WP), allocatable :: time(:)

      !> Component bitmasks: 1=write, 0=skip (Fortran 1-based; index n corresponds to rdata[n-1] / idata[n-1])
      integer(c_int) :: write_real(AMRPD_NREAL_PART)=0
      integer(c_int) :: write_int (AMRPD_NINT_PART) =0

   contains
      procedure :: initialize     !< Set up output dir, restore time series on restart
      procedure :: select_comp    !< Toggle a named field group on/off
      procedure :: write          !< Write one plotfile snapshot
      procedure :: finalize       !< Clean up
   end type amrpdviz

contains

   !> Initialize: create output directory, restore time series on restart.
   !> Reads all existing plotfile times from disk (sequentially numbered).
   !> Rewind/truncation based on the restart time happens in write().
   subroutine initialize(this,pd,name)
      use filesys,  only: makedir,isdir
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      implicit none
      class(amrpdviz), intent(inout) :: this
      class(amrpd), target, intent(in) :: pd
      character(len=*), intent(in) :: name

      character(len=str_long) :: pltdir
      integer :: ierr,n
      real(c_double) :: file_time

      this%pd => pd
      this%name = trim(adjustl(name))
      this%ntime = 0

      ! Default selection: write velocity (3 reals). Skip forces, mw, dil, flag.
      ! Position is always written by AMReX. Caller can override via select_comp.
      this%write_real = 0
      this%write_real(1:3) = 1     ! vel
      this%write_int  = 0

      ! Create output directory
      if (pd%amr%amRoot) then
         if (.not.isdir('amrviz')) call makedir('amrviz')
         if (.not.isdir('amrviz/'//trim(this%name))) &
            call makedir('amrviz/'//trim(this%name))
      end if

      ! Root probes for existing plotfiles via C++ function (mirrors amrviz pattern).
      ! amrpd_read_plotfile_time returns -1.0 when directory/time file is absent.
      if (pd%amr%amRoot) then
         n = 0
         find_files: do
            n = n + 1
            write(pltdir,'("amrviz/",a,"/plt.part.",i6.6)') trim(this%name),n
            file_time = amrpd_read_plotfile_time(trim(pltdir)//c_null_char)
            if (file_time.lt.0.0_c_double) exit find_files
         end do find_files
         this%ntime = n - 1
         if (this%ntime.gt.0) then
            allocate(this%time(this%ntime))
            do n = 1, this%ntime
               write(pltdir,'("amrviz/",a,"/plt.part.",i6.6)') trim(this%name),n
               this%time(n) = real(amrpd_read_plotfile_time(trim(pltdir)//c_null_char),WP)
            end do
         end if
      end if

      ! Broadcast ntime and time array to all ranks
      call MPI_BCAST(this%ntime,1,MPI_INTEGER,0,pd%amr%comm,ierr)
      if (this%ntime.gt.0) then
         if (.not.pd%amr%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time,this%ntime,MPI_REAL_WP,0,pd%amr%comm,ierr)
      end if

   end subroutine initialize


   !> Toggle a named field group on or off.
   !> Group names: 'vel', 'F_bond', 'F_fluid', 'mw', 'dil', 'damage', 'nb0', 'td2', 'td2a', 'flag'
   !> Individual components: 'vx','vy','vz','fbx','fby','fbz','ffx','ffy','ffz'
   !> Convenience: 'all', 'position_only'
   subroutine select_comp(this,name,on)
      implicit none
      class(amrpdviz), intent(inout) :: this
      character(len=*), intent(in) :: name
      logical, intent(in) :: on
      integer :: val

      val = merge(1, 0, on)

      select case(trim(name))
      ! Velocity (rdata[0..2])
      case('vel');     this%write_real(1:3) = val
      case('vx');      this%write_real(1)   = val
      case('vy');      this%write_real(2)   = val
      case('vz');      this%write_real(3)   = val
      ! Bond force (rdata[3..5])
      case('F_bond');  this%write_real(4:6) = val
      case('fbx');     this%write_real(4)   = val
      case('fby');     this%write_real(5)   = val
      case('fbz');     this%write_real(6)   = val
      ! Fluid force (rdata[6..8])
      case('F_fluid'); this%write_real(7:9) = val
      case('ffx');     this%write_real(7)   = val
      case('ffy');     this%write_real(8)   = val
      case('ffz');     this%write_real(9)   = val
      ! Scalar state (rdata[9..14])
      case('mw');      this%write_real(10)  = val
      case('dil');     this%write_real(11)  = val
      case('damage');  this%write_real(12)  = val
      case('nb0');     this%write_real(13)  = val
      case('td2');     this%write_real(14)  = val
      case('td2a');    this%write_real(15)  = val
      ! Integer (idata[0])
      case('flag');    this%write_int(1)    = val
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
   !> Directory: amrviz/<name>/plt.part.NNNNNN (6-digit zero-padded)
   subroutine write(this,time)
      use iso_c_binding, only: c_null_char,c_double
      implicit none
      class(amrpdviz), intent(inout) :: this
      real(WP), intent(in) :: time

      character(len=str_long) :: pltdir
      real(WP), allocatable :: tmp(:)
      integer :: i,n

      ! Update time array (same rewind logic as amrviz/amrlptviz: if new time is
      ! less than existing ones, truncate the series at the appropriate point)
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
      call amrpd_write_plotfile(this%pd%pcp, &
      &                         trim(pltdir)//c_null_char, &
      &                         'particles'//c_null_char,  &
      &                         this%write_real, this%write_int, &
      &                         real(time, c_double))

      ! Write/rewrite JSON .series file for ParaView time association
      if (this%pd%amr%amRoot) then
         open(newunit=n,file='amrviz/'//trim(this%name)//'/plt.part.series', &
         &    status='replace',action='write')
         write(n,'(a)') '{ "file-series-version": "1.0",'
         write(n,'(a)') '  "files": ['
         do i = 1, this%ntime
            write(pltdir,'("plt.part.",i6.6)') i
            if (i.lt.this%ntime) then
               write(n,'(4x,a,a,a,es24.17,a)') '{ "name": "', trim(pltdir), '", "time": ', this%time(i), ' },'
            else
               write(n,'(4x,a,a,a,es24.17,a)') '{ "name": "', trim(pltdir), '", "time": ', this%time(i), ' }'
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
      class(amrpdviz), intent(inout) :: this
      if (allocated(this%time)) deallocate(this%time)
      nullify(this%pd)
      this%ntime = 0
   end subroutine finalize

end module amrpdviz_class
