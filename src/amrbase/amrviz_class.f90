!> AMR Visualization class: plotfile output with field registration
!> Outputs HDF5 (Chombo format) or native AMReX plotfiles
!> readable by VisIt and ParaView
module amrviz_class
   use precision,      only: WP
   use string,         only: str_medium,str_long
   use iso_c_binding,  only: c_ptr,c_char,c_int,c_double,c_null_char,c_loc
   use amrgrid_class,  only: amrgrid
   use amrdata_class,  only: amrdata
   use surfmesh_class, only: surfmesh
   use mpi_f08,        only: MPI_Comm,MPI_BARRIER,MPI_BCAST,MPI_INTEGER,MPI_COMM_SIZE,MPI_COMM_RANK
   use parallel,       only: MPI_REAL_WP
   use amrex_interface, only: amrplotfile_write_hdf5,amrplotfile_write_native,amrplotfile_read_time
   use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box, &
   &                           amrex_multifab_build,amrex_multifab_destroy, &
   &                           amrex_mfiter_build,amrex_mfiter_destroy
   implicit none
   private

   public :: amrviz

   !> Scalar field registration node
   type :: scl
      type(scl), pointer :: next => null()
      character(len=str_medium) :: name
      type(amrdata), pointer :: ptr => null()
      integer :: comp
      logical :: recenter = .true.  !< If true, interpolate to cell centers
   end type scl

   !> Surface mesh registration node
   type :: srf
      type(srf), pointer :: next => null()
      character(len=str_medium) :: name
      type(surfmesh), pointer :: ptr => null()
   end type srf
   
   ! Base64 encoding table for VTP output
   character(len=64), parameter :: b64_table = &
      'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'

   !> AMR visualization handler
   type :: amrviz
      class(amrgrid), pointer :: amr => null()   !< Pointer to AMR grid
      character(len=str_medium) :: name          !< Subdirectory name for output
      type(scl), pointer :: first_scl => null()  !< Registered scalars
      type(srf), pointer :: first_srf => null()  !< Registered surface meshes
      integer :: ntime = 0                       !< File counter for output
      real(WP), allocatable :: time(:)           !< Time values for each file
      logical :: use_hdf5 = .true.               !< If true, HDF5/Chombo; if false, native AMReX
   contains
      procedure :: initialize                    !< Initialize with grid
      procedure :: add_scalar                    !< Register a scalar field
      procedure :: add_surfmesh                  !< Register a surface mesh
      procedure :: write                         !< Write all registered fields to single HDF5
      procedure :: finalize                      !< Clean up registered field lists
      procedure, private :: write_vtp            !< Write surface mesh to VTP
      procedure, private :: write_pvd            !< Write PVD collection file
   end type amrviz

contains

   !> Initialize visualization handler with AMR grid
   !> If output directory already exists with files, reads their times to continue series
   subroutine initialize(this, amr, name, use_hdf5)
      use filesys, only: makedir, isdir, isfile
      use mpi_f08, only: MPI_BCAST, MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrviz), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: use_hdf5

      character(len=str_medium) :: filename, dirname
      integer :: ierr, n
      real(WP) :: file_time

      this%amr => amr
      this%name = trim(adjustl(name))
      this%ntime = 0
      this%first_scl => null()
      if (present(use_hdf5)) this%use_hdf5 = use_hdf5

      ! Create output directory: amrviz/<name>/
      dirname = 'amrviz/'//trim(this%name)
      if (this%amr%amRoot) then
         if (.not.isdir('amrviz')) call makedir('amrviz')
         if (.not.isdir(trim(dirname))) call makedir(trim(dirname))
      end if

      ! Check for existing files and read their times (root only, then broadcast)
      ! Scan for any centering pattern - try each until we find files
      if (this%amr%amRoot) then
         centering_names: block
            character(len=8), parameter :: centerings(8) = &
               ['cell   ', 'xface  ', 'yface  ', 'zface  ', 'node   ', 'xyedge ', 'xzedge ', 'yzedge ']
            character(len=8) :: found_centering
            integer :: ic

            ! Find first centering type that has files or directories
            found_centering = ''
            find_centering: do ic = 1, 8
               write(filename,'(a,"/plt.nga2.",a,".",i6.6)') trim(dirname), trim(centerings(ic)), 1
               if (isfile(trim(filename)).or.isdir(trim(filename))) then
                  found_centering = centerings(ic)
                  exit find_centering
               end if
            end do find_centering

            if (len_trim(found_centering) > 0) then
               ! Count existing files using found centering
               n = 0
               do
                  n = n + 1
                  write(filename,'(a,"/plt.nga2.",a,".",i6.6)') trim(dirname), trim(found_centering), n
                  if (.not.isfile(trim(filename)).and..not.isdir(trim(filename))) exit
                  file_time = amrplotfile_read_time(trim(filename)//c_null_char)
                  if (file_time .lt. 0.0_WP) exit  ! File exists but unreadable
               end do
               this%ntime = n - 1

               ! Allocate and fill time array
               if (this%ntime .gt. 0) then
                  allocate(this%time(this%ntime))
                  do n = 1, this%ntime
                     write(filename,'(a,"/plt.nga2.",a,".",i6.6)') trim(dirname), trim(found_centering), n
                     this%time(n) = amrplotfile_read_time(trim(filename)//c_null_char)
                  end do
               end if
            end if
         end block centering_names
      end if

      ! Broadcast ntime and time array to all ranks
      call MPI_BCAST(this%ntime, 1, MPI_INTEGER, 0, this%amr%comm, ierr)
      if (this%ntime .gt. 0) then
         if (.not.this%amr%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time, this%ntime, MPI_REAL_WP, 0, this%amr%comm, ierr)
      end if
   end subroutine initialize


   !> Register a scalar field for output
   !> recenter: if .true. (default), interpolate to cell centers; if .false., use native centering
   subroutine add_scalar(this, data, comp, name, recenter)
      use messager, only: die
      implicit none
      class(amrviz), intent(inout) :: this
      type(amrdata), target, intent(in) :: data
      integer, intent(in) :: comp
      character(len=*), intent(in), optional :: name
      logical, intent(in), optional :: recenter
      type(scl), pointer :: new_scl

      ! Check that the component is meaningful
      if (comp.le.0) call die('[amrviz add_scalar] comp must be at least one')
      if (comp.gt.data%ncomp) call die('[amrviz add_scalar] comp is too large for provided amrdata')

      ! Create new scalar node
      allocate(new_scl)
      if (present(name)) then
         new_scl%name = trim(adjustl(name))
      else
         new_scl%name = trim(adjustl(data%name))
      end if
      new_scl%ptr => data
      new_scl%comp = comp
      if (present(recenter)) then
         new_scl%recenter = recenter
      else
         new_scl%recenter = .true.
      end if

      ! Insert at front of list
      new_scl%next => this%first_scl
      this%first_scl => new_scl
   end subroutine add_scalar


   !> Register a surface mesh for output (VTP format)
   subroutine add_surfmesh(this, smesh, name)
      implicit none
      class(amrviz), intent(inout) :: this
      type(surfmesh), target, intent(in) :: smesh
      character(len=*), intent(in) :: name
      type(srf), pointer :: new_srf

      ! Create new surface node
      allocate(new_srf)
      new_srf%name = trim(adjustl(name))
      new_srf%ptr => smesh

      ! Insert at front of list
      new_srf%next => this%first_srf
      this%first_srf => new_srf
   end subroutine add_surfmesh

   !> Write all registered fields to HDF5 plotfiles
   !> Fields are grouped by centering type - one file per centering
   !> File pattern: amrviz/<name>/plt.nga2.<centering>.<ntime>.h5
   !> Centerings: cell, xface, yface, zface, xyedge, xzedge, yzedge, node
   subroutine write(this, time)
      implicit none
      class(amrviz), intent(inout) :: this
      real(WP), intent(in) :: time

      integer :: nlev, lev, nscl_grp, ncomp, icomp, i, n, ngroups
      type(scl), pointer :: my_scl
      real(WP), dimension(:), allocatable :: temp_time
      character(len=str_medium) :: filename
      character(len=8) :: suffix

      ! Centering groups: store unique nodal patterns found
      ! Max 8 possible centerings (2^3), but typically only 1-4 used
      integer, parameter :: MAX_GROUPS = 8
      logical :: group_nodal(3, MAX_GROUPS)
      integer :: group_count(MAX_GROUPS)
      logical :: found, eff_nodal(3), needs_interp, native_nodal(3)
      integer :: ii, jj, kk, di, dj, dk

      ! Working arrays for current group
      type(amrex_multifab), allocatable :: combined(:)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: dst, src
      type(c_ptr), allocatable :: mf_ptrs(:), geom_ptrs(:)
      integer(c_int), allocatable :: level_steps(:), ref_ratios(:)
      character(len=64), target, allocatable :: varnames(:)
      type(c_ptr), allocatable :: varname_ptrs(:)

      nlev = this%amr%clvl() + 1
      ngroups = 0
      group_count = 0

      ! --- Phase 1: Collect unique centering types and count fields per group ---
      ! If recenter=.true., effective centering is cell-centered
      my_scl => this%first_scl
      do while (associated(my_scl))
         ! Determine effective centering for output
         if (my_scl%recenter) then
            eff_nodal = [.false., .false., .false.]  ! Recenter to cell
         else
            eff_nodal = my_scl%ptr%nodal             ! Native centering
         end if

         ! Check if this centering already exists in groups
         found = .false.
         do i = 1, ngroups
            if (all(eff_nodal .eqv. group_nodal(:,i))) then
               group_count(i) = group_count(i) + 1
               found = .true.
               exit
            end if
         end do
         ! New centering type
         if (.not.found) then
            ngroups = ngroups + 1
            group_nodal(:,ngroups) = eff_nodal
            group_count(ngroups) = 1
         end if
         my_scl => my_scl%next
      end do

      if (ngroups == 0) return  ! Nothing to write

      ! --- Time tracking (shared across all groups) ---
      if (this%ntime.eq.0) then
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(this%ntime))
         this%time(1)=time
      else
         n=1
         rewind: do i=this%ntime,1,-1
            if (this%time(i).lt.time-1.0e-6_WP) then
               n=i+1; exit rewind
            end if
         end do rewind
         this%ntime=n; allocate(temp_time(1:this%ntime))
         temp_time=[this%time(1:this%ntime-1),time]
         call move_alloc(temp_time,this%time)
      end if

      ! --- Phase 2: Write one file per centering group ---
      do i = 1, ngroups
         nscl_grp = group_count(i)
         ncomp = nscl_grp  ! For now, scalars only per group

         ! Get centering suffix
         suffix = centering_suffix(group_nodal(:,i))

         ! Allocate varnames for this group
         allocate(varnames(ncomp))
         allocate(varname_ptrs(ncomp))
         icomp = 0

         ! Collect varnames for matching scalars (using effective centering)
         my_scl => this%first_scl
         do while (associated(my_scl))
            if (my_scl%recenter) then
               eff_nodal = [.false., .false., .false.]
            else
               eff_nodal = my_scl%ptr%nodal
            end if
            if (all(eff_nodal .eqv. group_nodal(:,i))) then
               icomp = icomp + 1
               varnames(icomp) = trim(my_scl%name)//c_null_char
               varname_ptrs(icomp) = c_loc(varnames(icomp))
            end if
            my_scl => my_scl%next
         end do

         ! Allocate combined MultiFab for each level
         allocate(combined(0:this%amr%clvl()))

         ! Build combined MF and copy matching scalar data
         do lev = 0, this%amr%clvl()
            call this%amr%mfab_build(lvl=lev, mfab=combined(lev), ncomp=ncomp, nover=0, atface=group_nodal(:,i))

            icomp = 0
            my_scl => this%first_scl
            do while (associated(my_scl))
               ! Compute effective nodal for matching
               if (my_scl%recenter) then
                  eff_nodal = [.false., .false., .false.]
               else
                  eff_nodal = my_scl%ptr%nodal
               end if

               if (all(eff_nodal .eqv. group_nodal(:,i))) then
                  icomp = icomp + 1
                  ! Check if interpolation is needed (native != effective)
                  needs_interp = my_scl%recenter .and. any(my_scl%ptr%nodal)

                  call amrex_mfiter_build(mfi, combined(lev), tiling=this%amr%default_tiling)
                  do while (mfi%next())
                     bx = mfi%tilebox()
                     src => my_scl%ptr%mf(lev)%dataptr(mfi)
                     dst => combined(lev)%dataptr(mfi)

                     if (needs_interp) then
                        ! Interpolate nodal/face/edge data to cell centers
                        ! Average 2^n points where n = number of nodal directions
                        native_nodal = my_scl%ptr%nodal
                        do kk = bx%lo(3), bx%hi(3)
                           do jj = bx%lo(2), bx%hi(2)
                              do ii = bx%lo(1), bx%hi(1)
                                 ! Sum all corner contributions based on nodal directions
                                 dst(ii,jj,kk,icomp) = 0.0_WP
                                 do dk = 0, merge(1,0,native_nodal(3))
                                    do dj = 0, merge(1,0,native_nodal(2))
                                       do di = 0, merge(1,0,native_nodal(1))
                                          dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) + src(ii+di, jj+dj, kk+dk, my_scl%comp)
                                       end do
                                    end do
                                 end do
                                 ! Divide by number of points (2^nodal_count)
                                 dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) / real(2**count(native_nodal), WP)
                              end do
                           end do
                        end do
                     else
                        ! Direct copy (same centering)
                        dst(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), icomp) = &
                           src(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), my_scl%comp)
                     end if
                  end do
                  call amrex_mfiter_destroy(mfi)
               end if
               my_scl => my_scl%next
            end do
         end do

         ! Generate filename with centering type
         filename = 'amrviz/'//trim(this%name)//'/plt.nga2.'//trim(suffix)//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime

         ! Prepare pointers for HDF5 writer
         allocate(mf_ptrs(0:this%amr%clvl()))
         allocate(geom_ptrs(0:this%amr%clvl()))
         allocate(level_steps(0:this%amr%clvl()))
         allocate(ref_ratios(0:max(3*this%amr%clvl()-1,0)))

         do lev = 0, this%amr%clvl()
            mf_ptrs(lev) = combined(lev)%p
            geom_ptrs(lev) = this%amr%geom(lev)%p
            level_steps(lev) = this%ntime
         end do
         do lev = 0, this%amr%clvl()-1
            ref_ratios(3*lev+0) = this%amr%rrefx(lev)
            ref_ratios(3*lev+1) = this%amr%rrefy(lev)
            ref_ratios(3*lev+2) = this%amr%rrefz(lev)
         end do

         ! Write plotfile for this centering group
         if (this%use_hdf5) then
            call amrplotfile_write_hdf5(trim(filename)//c_null_char, nlev, mf_ptrs, &
               varname_ptrs, ncomp, geom_ptrs, real(time,c_double), level_steps, &
               ref_ratios, c_null_char)
         else
            call amrplotfile_write_native(trim(filename)//c_null_char, nlev, mf_ptrs, &
               varname_ptrs, ncomp, geom_ptrs, real(time,c_double), level_steps, &
               ref_ratios)
         end if

         ! Cleanup for this group
         do lev = 0, this%amr%clvl()
            call amrex_multifab_destroy(combined(lev))
         end do
         deallocate(combined)
         deallocate(mf_ptrs, geom_ptrs, level_steps, ref_ratios)
         deallocate(varnames, varname_ptrs)
      end do

      ! --- Write registered surface meshes to VTP ---
      write_surfaces: block
         type(srf), pointer :: my_srf
         my_srf => this%first_srf
         do while (associated(my_srf))
            call this%write_vtp(my_srf%name, my_srf%ptr)
            my_srf => my_srf%next
         end do
      end block write_surfaces

   contains

      !> Return filename suffix based on centering
      function centering_suffix(nodal) result(suf)
         logical, intent(in) :: nodal(3)
         character(len=8) :: suf
         if (.not.any(nodal)) then
            suf = 'cell'     ! Cell-centered
         else if (all(nodal)) then
            suf = 'node'     ! Node-centered
         else if (nodal(1) .and. .not.nodal(2) .and. .not.nodal(3)) then
            suf = 'xface'    ! X-face
         else if (.not.nodal(1) .and. nodal(2) .and. .not.nodal(3)) then
            suf = 'yface'    ! Y-face
         else if (.not.nodal(1) .and. .not.nodal(2) .and. nodal(3)) then
            suf = 'zface'    ! Z-face
         else if (nodal(1) .and. nodal(2) .and. .not.nodal(3)) then
            suf = 'xyedge'   ! XY-edge
         else if (nodal(1) .and. .not.nodal(2) .and. nodal(3)) then
            suf = 'xzedge'   ! XZ-edge
         else
            suf = 'yzedge'   ! YZ-edge
         end if
      end function centering_suffix

   end subroutine write


   !> Finalize: clean up registered field lists
   subroutine finalize(this)
      implicit none
      class(amrviz), intent(inout) :: this
      type(scl), pointer :: cur_scl, next_scl
      type(srf), pointer :: cur_srf, next_srf

      ! Clean up scalar list
      cur_scl => this%first_scl
      do while (associated(cur_scl))
         next_scl => cur_scl%next
         deallocate(cur_scl)
         cur_scl => next_scl
      end do
      nullify(this%first_scl)

      ! Clean up surface mesh list
      cur_srf => this%first_srf
      do while (associated(cur_srf))
         next_srf => cur_srf%next
         deallocate(cur_srf)
         cur_srf => next_srf
      end do
      nullify(this%first_srf)

      ! Reset state
      nullify(this%amr)
      this%ntime = 0
      if (allocated(this%time)) deallocate(this%time)
   end subroutine finalize


   !> Encode byte array to base64 string for VTP binary output
   function encode_base64(bytes, nbytes) result(encoded)
      implicit none
      integer(1), dimension(:), intent(in) :: bytes
      integer, intent(in) :: nbytes
      character(len=:), allocatable :: encoded
      integer :: i, j, nout, b1, b2, b3, idx
      integer :: npad
      
      ! Calculate output length (4 chars per 3 bytes, rounded up)
      nout = ((nbytes + 2) / 3) * 4
      allocate(character(len=nout) :: encoded)
      
      j = 1
      do i = 1, nbytes, 3
         ! Get up to 3 bytes (use 0 for padding)
         b1 = iand(int(bytes(i), kind=4), 255)
         if (i+1 <= nbytes) then
            b2 = iand(int(bytes(i+1), kind=4), 255)
         else
            b2 = 0
         end if
         if (i+2 <= nbytes) then
            b3 = iand(int(bytes(i+2), kind=4), 255)
         else
            b3 = 0
         end if
         
         ! Encode to 4 base64 characters
         idx = ishft(b1, -2) + 1
         encoded(j:j) = b64_table(idx:idx)
         
         idx = ior(ishft(iand(b1, 3), 4), ishft(b2, -4)) + 1
         encoded(j+1:j+1) = b64_table(idx:idx)
         
         idx = ior(ishft(iand(b2, 15), 2), ishft(b3, -6)) + 1
         encoded(j+2:j+2) = b64_table(idx:idx)
         
         idx = iand(b3, 63) + 1
         encoded(j+3:j+3) = b64_table(idx:idx)
         
         j = j + 4
      end do
      
      ! Add padding
      npad = mod(3 - mod(nbytes, 3), 3)
      if (npad >= 1) encoded(nout:nout) = '='
      if (npad >= 2) encoded(nout-1:nout-1) = '='
      
   end function encode_base64


   !> Write surface mesh to VTP file (binary base64 format)
   subroutine write_vtp(this, srf_name, smesh)
      implicit none
      class(amrviz), intent(in) :: this
      character(len=*), intent(in) :: srf_name
      type(surfmesh), intent(in) :: smesh
      
      character(len=str_long) :: filename, dirname
      character(len=str_medium) :: basename
      character(len=:), allocatable :: b64_data
      integer :: iunit, ierr, irank, n, v, conn_idx, vert_offset
      integer :: npts_bytes, nconn_bytes, noff_bytes, nvar_bytes
      integer(1), dimension(:), allocatable :: buffer
      integer(4) :: header_size
      real(WP), dimension(:), allocatable :: pts_data
      integer(4), dimension(:), allocatable :: conn_data, off_data
      integer :: nproc, rank
      
      ! Get MPI info
      call MPI_COMM_SIZE(this%amr%comm, nproc, ierr)
      call MPI_COMM_RANK(this%amr%comm, rank, ierr)
      
      ! Construct filename with timestep
      dirname = 'amrviz/'//trim(this%name)
      write(basename,'(A,"_",I6.6,".vtp")') trim(srf_name), this%ntime
      filename = trim(dirname)//'/'//trim(basename)
      
      ! Rank 0 creates header
      if (rank.eq.0) then
         open(newunit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr)
         write(iunit,'(a)') '<?xml version="1.0"?>'
         write(iunit,'(a)') '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt32">'
         write(iunit,'(a)') '  <PolyData>'
         close(iunit)
      end if
      call MPI_BARRIER(this%amr%comm, ierr)
      
      ! Each rank writes its data sequentially
      do irank = 0, nproc - 1
         if (irank.eq.rank) then
            open(newunit=iunit, file=trim(filename), status='old', position='append', action='write', iostat=ierr)
            
            ! Write this rank's piece
            if (smesh%nPoly.gt.0 .and. smesh%nVert.gt.0) then
               write(iunit,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', smesh%nVert, &
                  '" NumberOfPolys="', smesh%nPoly, '">'
               
               ! === Points (base64 binary) ===
               write(iunit,'(a)') '      <Points>'
               write(iunit,'(a)', advance='no') '        <DataArray type="Float64" NumberOfComponents="3" format="binary">'
               
               ! Pack points into byte buffer: [header(4 bytes)][data]
               npts_bytes = smesh%nVert * 3 * 8
               allocate(pts_data(smesh%nVert * 3))
               do v = 1, smesh%nVert
                  pts_data((v-1)*3 + 1) = smesh%xVert(v)
                  pts_data((v-1)*3 + 2) = smesh%yVert(v)
                  pts_data((v-1)*3 + 3) = smesh%zVert(v)
               end do
               allocate(buffer(4 + npts_bytes))
               header_size = int(npts_bytes, 4)
               buffer(1:4) = transfer(header_size, buffer(1:4))
               buffer(5:) = transfer(pts_data, buffer(5:))
               b64_data = encode_base64(buffer, size(buffer))
               write(iunit,'(a)', advance='no') b64_data
               deallocate(buffer, pts_data, b64_data)
               
               write(iunit,'(a)') '</DataArray>'
               write(iunit,'(a)') '      </Points>'
               
               ! === Polys (connectivity + offsets) ===
               write(iunit,'(a)') '      <Polys>'
               
               ! Connectivity
               write(iunit,'(a)', advance='no') '        <DataArray type="Int32" Name="connectivity" format="binary">'
               nconn_bytes = sum(smesh%polySize(1:smesh%nPoly)) * 4
               allocate(conn_data(sum(smesh%polySize(1:smesh%nPoly))))
               conn_idx = 1
               do n = 1, smesh%nPoly
                  do v = 1, smesh%polySize(n)
                     conn_data(conn_idx) = int(smesh%polyConn(conn_idx) - 1, 4)  ! 0-based for VTK
                     conn_idx = conn_idx + 1
                  end do
               end do
               allocate(buffer(4 + nconn_bytes))
               header_size = int(nconn_bytes, 4)
               buffer(1:4) = transfer(header_size, buffer(1:4))
               buffer(5:) = transfer(conn_data, buffer(5:))
               b64_data = encode_base64(buffer, size(buffer))
               write(iunit,'(a)', advance='no') b64_data
               deallocate(buffer, conn_data, b64_data)
               write(iunit,'(a)') '</DataArray>'
               
               ! Offsets
               write(iunit,'(a)', advance='no') '        <DataArray type="Int32" Name="offsets" format="binary">'
               noff_bytes = smesh%nPoly * 4
               allocate(off_data(smesh%nPoly))
               vert_offset = 0
               do n = 1, smesh%nPoly
                  vert_offset = vert_offset + smesh%polySize(n)
                  off_data(n) = int(vert_offset, 4)
               end do
               allocate(buffer(4 + noff_bytes))
               header_size = int(noff_bytes, 4)
               buffer(1:4) = transfer(header_size, buffer(1:4))
               buffer(5:) = transfer(off_data, buffer(5:))
               b64_data = encode_base64(buffer, size(buffer))
               write(iunit,'(a)', advance='no') b64_data
               deallocate(buffer, off_data, b64_data)
               write(iunit,'(a)') '</DataArray>'
               
               write(iunit,'(a)') '      </Polys>'
               
               ! === Cell data (per-polygon variables) ===
               if (smesh%nvar.gt.0 .and. allocated(smesh%var)) then
                  write(iunit,'(a)') '      <CellData>'
                  do v = 1, smesh%nvar
                     write(iunit,'(a,a,a)', advance='no') '        <DataArray type="Float64" Name="', &
                        trim(smesh%varname(v)), '" format="binary">'
                     nvar_bytes = smesh%nPoly * 8
                     allocate(buffer(4 + nvar_bytes))
                     header_size = int(nvar_bytes, 4)
                     buffer(1:4) = transfer(header_size, buffer(1:4))
                     buffer(5:) = transfer(smesh%var(v,1:smesh%nPoly), buffer(5:))
                     b64_data = encode_base64(buffer, size(buffer))
                     write(iunit,'(a)', advance='no') b64_data
                     deallocate(buffer, b64_data)
                     write(iunit,'(a)') '</DataArray>'
                  end do
                  write(iunit,'(a)') '      </CellData>'
               end if
               
               write(iunit,'(a)') '    </Piece>'
            end if
            
            close(iunit)
         end if
         call MPI_BARRIER(this%amr%comm, ierr)
      end do
      
      ! Rank 0 writes footer
      if (rank.eq.0) then
         open(newunit=iunit, file=trim(filename), status='old', position='append', action='write', iostat=ierr)
         write(iunit,'(a)') '  </PolyData>'
         write(iunit,'(a)') '</VTKFile>'
         close(iunit)
      end if
      call MPI_BARRIER(this%amr%comm, ierr)
      
      ! Write PVD collection file
      call this%write_pvd(srf_name, rank)
      
   end subroutine write_vtp


   !> Write PVD collection file for time series
   subroutine write_pvd(this, srf_name, rank)
      implicit none
      class(amrviz), intent(in) :: this
      character(len=*), intent(in) :: srf_name
      integer, intent(in) :: rank
      character(len=str_long) :: filename, dirname
      character(len=str_medium) :: basename
      integer :: iunit, ierr, n
      
      ! Only rank 0 writes PVD
      if (rank.ne.0) return
      
      dirname = 'amrviz/'//trim(this%name)
      filename = trim(dirname)//'/'//trim(srf_name)//'.pvd'
      
      open(newunit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr)
      
      write(iunit,'(a)') '<?xml version="1.0"?>'
      write(iunit,'(a)') '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">'
      write(iunit,'(a)') '  <Collection>'
      
      do n = 1, this%ntime
         write(basename,'(A,"_",I6.6,".vtp")') trim(srf_name), n
         write(iunit,'(a,es18.10,a,a,a)') '    <DataSet timestep="', this%time(n), &
            '" file="', trim(basename), '"/>'
      end do
      
      write(iunit,'(a)') '  </Collection>'
      write(iunit,'(a)') '</VTKFile>'
      
      close(iunit)
      
   end subroutine write_pvd


end module amrviz_class
