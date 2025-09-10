!> Ensight class concept is defined here: given a config object,
!> it provides parallel I/O access to an vtk file
module vtk_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype
   use surfmesh_class, only: surfmesh
   use partmesh_class, only: partmesh
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: vtk
   
   ! List types
   type :: scl !< Scalar field
      type(scl), pointer :: next
      character(len=str_medium) :: name
      real(WP), dimension(:,:,:), pointer :: rptr=>NULL()  !< real(WP) data
      integer , dimension(:,:,:), pointer :: iptr=>NULL()  !< integer  data
   end type scl
   type :: vct !< Vector field
      type(vct), pointer :: next
      character(len=str_medium) :: name
      real(WP), dimension(:,:,:), pointer :: ptrx
      real(WP), dimension(:,:,:), pointer :: ptry
      real(WP), dimension(:,:,:), pointer :: ptrz
   end type vct
   type :: srf !< Surface mesh
      type(srf), pointer :: next
      character(len=str_medium) :: name
      type(surfmesh), pointer :: ptr
   end type srf
   type :: prt !< Particle mesh
      type(prt), pointer :: next
      character(len=str_medium) :: name
      type(partmesh), pointer :: ptr
   end type prt
   
   !> Ensight object definition as list of pointers to arrays
   type :: vtk
      ! An vtk object has a name
      character(len=str_medium) :: name                               !< Name of vtk directory to read/write
      ! An vtk object stores time values
      integer :: ntime                                                !< Number of scalar values
      real(WP), dimension(:), allocatable :: time                     !< Time values
      ! An vtk object stores geometry data
      type(config), pointer :: cfg                                    !< Config for vtk geometry and parallel I/O
      ! An vtk object stores lists of pointers to data
      type(scl), pointer :: first_scl                                 !< Scalar list
      type(vct), pointer :: first_vct                                 !< Vector list
      type(srf), pointer :: first_srf                                 !< Surface list
      type(prt), pointer :: first_prt                                 !< Particle list
   contains
      procedure :: write_data                                         !< Write out data
      procedure :: write_case                                         !< Write out case file
      procedure :: write_surf                                         !< Write out surface mesh file
      procedure :: write_part                                         !< Write out particle mesh file
      generic :: add_scalar=>add_rscalar,add_iscalar                  !< Add a new scalar field
      procedure, private :: add_rscalar                               !< Add a new real(WP) scalar field
      procedure, private :: add_iscalar                               !< Add a new integer  scalar field
      procedure :: add_vector                                         !< Add a new vector field
      procedure :: add_surface                                        !< Add a new surface mesh
      procedure :: add_particle                                       !< Add a new particle mesh
   end type vtk
   
   
   !> Declare vtk constructor
   interface vtk
      procedure construct_vtk
   end interface vtk
   
   
contains
   
   !> Constructor for an empty vtk object
   function construct_vtk(cfg,name) result(self)
      use messager, only: die
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      use filesys,  only: makedir,isdir
      implicit none
      type(vtk) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      character(len=str_medium) :: line
      integer :: iunit,ierr,stat,n
      logical :: file_is_there,found
      character(len=2), parameter :: idt = '  '
      
      ! Link to config
      self%cfg=>cfg
      
      ! Store casename
      self%name=trim(adjustl(name))
      
      ! Start with no time stamps
      self%ntime=0
      
      ! Create directory
      if (self%cfg%amRoot) then
         if (.not.isdir('vtk')) &
         & call makedir('vtk')
         if (.not.isdir('vtk/'//trim(self%name))) &
         & call makedir('vtk/'//trim(self%name))
         if (.not.isdir('vtk/'//trim(self%name)//'/data')) &
         & call makedir('vtk/'//trim(self%name)//'/data')
      end if
      
      ! Empty pointer to lists for now
      self%first_scl=>NULL()
      self%first_vct=>NULL()
      self%first_srf=>NULL()
      self%first_prt=>NULL()
      
      ! Check if a case file exists already - root only
      if (self%cfg%amRoot) then
         inquire(file='vtk/'//trim(self%name)//'/nga.pvd',exist=file_is_there)
         if (file_is_there) then
            ! Open the case file
            open(newunit=iunit,file='vtk/'//trim(self%name)//'/nga.pvd',form='formatted',status='old',access='stream',iostat=ierr)
            ! Read lines until we find time values section
            stat=0; found=.false.
            do while (.not.found.and..not.is_iostat_end(stat))
               read(iunit,'(a)',iostat=stat) line
               if (line(3:20).eq.'<Collection nsteps=') found=.true.
               if (found) read(line(22:27),'(i6)') self%ntime
            end do
            allocate(self%time(self%ntime))
            do n=1,self%ntime
               read(iunit,'(a)') line
               read(line(22:27),'(es12.5)') self%time(n)
            end do
            ! Close the case file
            close(iunit)
         else
            ! Output a barebone case file with only the geometry
            open(newunit=iunit,file='vtk/'//trim(self%name)//'/nga.pvd',form='formatted',status='replace',access='stream',iostat=ierr)
            write(iunit,'(a)')   '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
            write(iunit,'(a)')   idt//'<Collection nsteps="000000">'
            write(iunit,'(a)')   idt//'</Collection>'
            write(iunit,'(a)')   '</VTKFile>'
            close(iunit)
         end if
      end if
      
      ! Communicate to all processors
      call MPI_BCAST(self%ntime,1,MPI_INTEGER,0,self%cfg%comm,ierr)
      if (self%ntime.gt.0) then
         if (.not.self%cfg%amRoot) allocate(self%time(self%ntime))
         call MPI_BCAST(self%time,self%ntime,MPI_REAL_WP,0,self%cfg%comm,ierr)
      end if
      
   end function construct_vtk
   
   
   !> Add a real scalar field for output
   subroutine add_rscalar(this,name,scalar)
      implicit none
      class(vtk), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: scalar
      type(scl), pointer :: new_scl
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%rptr=>scalar
      new_scl%iptr=>NULL()
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
   end subroutine add_rscalar
   
   
   !> Add an integer scalar field for output
   subroutine add_iscalar(this,name,scalar)
      implicit none
      class(vtk), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: scalar
      type(scl), pointer :: new_scl
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%rptr=>NULL()
      new_scl%iptr=>scalar
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
   end subroutine add_iscalar
   
   
   !> Add a vector field for output
   subroutine add_vector(this,name,vectx,vecty,vectz)
      implicit none
      class(vtk), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: vectx
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: vecty
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: vectz
      type(vct), pointer :: new_vct
      ! Prepare new vector
      allocate(new_vct)
      new_vct%name=trim(adjustl(name))
      new_vct%ptrx=>vectx
      new_vct%ptry=>vecty
      new_vct%ptrz=>vectz
      ! Insert it up front
      new_vct%next=>this%first_vct
      ! Point list to new object
      this%first_vct=>new_vct
   end subroutine add_vector
   
   
   !> Add a surface mesh for output
   subroutine add_surface(this,name,surface)
      use filesys,  only: makedir,isdir
      implicit none
      class(vtk), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(surfmesh), target, intent(in) :: surface
      type(srf), pointer :: new_srf
      ! Prepare new surface
      allocate(new_srf)
      new_srf%name=trim(adjustl(name))
      new_srf%ptr =>surface
      ! Insert it up front
      new_srf%next=>this%first_srf
      ! Point list to new object
      this%first_srf=>new_srf
      ! Also create the corresponding directory
      if (this%cfg%amRoot) then
         if (.not.isdir('vtk/'//trim(this%name)//'/'//trim(new_srf%name))) &
         & call makedir('vtk/'//trim(this%name)//'/'//trim(new_srf%name))
      end if
   end subroutine add_surface
   
   
   !> Add a particle mesh for output
   subroutine add_particle(this,name,particle)
      use filesys,  only: makedir,isdir
      implicit none
      class(vtk), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(partmesh), target, intent(in) :: particle
      type(prt), pointer :: new_prt
      ! Prepare new particle
      allocate(new_prt)
      new_prt%name=trim(adjustl(name))
      new_prt%ptr =>particle
      ! Insert it up front
      new_prt%next=>this%first_prt
      ! Point list to new object
      this%first_prt=>new_prt
      ! Also create the corresponding directory
      if (this%cfg%amRoot) then
         if (.not.isdir('vtk/'//trim(this%name)//'/'//trim(new_prt%name))) &
         & call makedir('vtk/'//trim(this%name)//'/'//trim(new_prt%name))
      end if
   end subroutine add_particle
   
   !> Output all data in the object
   subroutine write_data(this,time)
      use precision, only: SP,I4
      use messager,  only: die
      use parallel,  only: info_mpiio,MPI_REAL_SP
      use mpi_f08
      implicit none
      class(vtk), intent(inout) :: this
      real(WP), intent(in) :: time
      character(len=str_medium) :: filename
      integer :: iunit,ierr,n,i,offset
      type(MPI_File) :: ifile
      integer(kind=MPI_OFFSET_KIND) :: disp
      type(MPI_Status):: status
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      type(srf), pointer :: my_srf
      type(prt), pointer :: my_prt
      real(SP), dimension(:,:,:),   allocatable :: spbuff
      real(SP), dimension(:,:,:,:), allocatable :: sp3buff
      real(WP), dimension(:), allocatable :: temp_time
      character(len=str_medium) :: ctime
      real(WP) :: rtime
      character(len=1), parameter :: eol = char(10)
      character(len=2), parameter :: idt = '  '
      integer(I4) :: data_size

      ! Check provided time stamp and decide what to do
      if (this%ntime.eq.0) then
         ! First time stamp
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(this%ntime))
         this%time(1)=time
      else
         ! There are time stamps already, check where to insert
         n=1
         rewind: do i=this%ntime,1,-1
            ! Convert time to appropriate accuracy before comparing
            ctime=''; write(ctime,'(es12.5)') time; read(ctime,'(es12.5)') rtime
            if (this%time(i).lt.rtime) then
               n=i+1; exit rewind
            end if
         end do rewind
         this%ntime=n; allocate(temp_time(1:this%ntime))
         temp_time=[this%time(1:this%ntime-1),time]
         call move_alloc(temp_time,this%time)
      end if
      
      ! Prepare the SP and SP3 buffer
      if (associated(this%first_scl))allocate(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))
      if (associated(this%first_vct)) allocate(sp3buff(1:3,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))

      filename='vtk/'//trim(this%name)//'/data/data.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      filename=trim(filename)//'.vtr'
      
      ! Start by writing ASCII header. Keep track of raw binary data offset
      offset=0
      if (this%cfg%amRoot) then
         open(newunit=iunit,file=filename,status='replace',form='formatted',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
         write(iunit,'(a)')             '<?xml version="1.0"?>'
         write(iunit,'(a)')             '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'
         write(iunit,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)')  idt//'<RectilinearGrid WholeExtent="',this%cfg%imin-1,' ',this%cfg%imax,' ',this%cfg%jmin-1,' ',this%cfg%jmax,' ',this%cfg%kmin-1,' ',this%cfg%kmax,'">'
         write(iunit,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)')  idt//'<Piece Extent="',this%cfg%imin-1,' ',this%cfg%imax,' ',this%cfg%jmin-1,' ',this%cfg%jmax,' ',this%cfg%kmin-1,' ',this%cfg%kmax,'">'
         write(iunit,'(a)')               idt//idt//'<Coordinates>'
         write(iunit,'(a,i0,a,i0,a)')     idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="x" format="appended" offset="',offset,'">'
         write(iunit,'(a)')               idt//idt//idt//'</DataArray>'
         offset=offset+I4+(this%cfg%nx+1)*SP
         write(iunit,'(a,i0,a,i0,a)')     idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="y" format="appended" offset="',offset,'">'
         write(iunit,'(a)')               idt//idt//idt//'</DataArray>'
         offset=offset+I4+(this%cfg%ny+1)*SP
         write(iunit,'(a,i0,a,i0,a)')     idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="z" format="appended" offset="',offset,'">'
         write(iunit,'(a)')               idt//idt//idt//'</DataArray>'
         offset=offset+I4+(this%cfg%nz+1)*SP
         write(iunit,'(a)')               idt//idt//'</Coordinates>'
         write(iunit,'(a)')               idt//idt//'<CellData>'
         ! Traverse all datasets and print them all out - scalars first
         my_scl=>this%first_scl
         do while (associated(my_scl))
            ! Write data array header and store offset
            write(iunit,'(a,i0,a,i0,a)')  idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="'//trim(my_scl%name)//'" format="appended" offset="',offset,'">'
            write(iunit,'(a)')            idt//idt//idt//'</DataArray>'
            offset=offset+I4+(this%cfg%nx*this%cfg%ny*this%cfg%nz)*SP
            ! Continue on to the next scalar object
            my_scl=>my_scl%next
         end do
         ! Traverse all datasets and print them all out - vectors second
         my_vct=>this%first_vct
         do while (associated(my_vct))
            ! Write vector header and store offset
            write(iunit,'(a,i0,a,i0,a)')  idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="'//trim(my_vct%name)//'" NumberOfComponents="3" format="appended" offset="',offset,'">'
            write(iunit,'(a)')            idt//idt//idt//'</DataArray>'
            offset=offset+I4+3*(this%cfg%nx*this%cfg%ny*this%cfg%nz)*SP
            ! Continue on to the next vector object
            my_vct=>my_vct%next
         end do
         write(iunit,'(a)')               idt//idt//'</CellData>'
         write(iunit,'(a)')               idt//'</Piece>'
         write(iunit,'(a)')               idt//'</RectilinearGrid>'
         ! Close the file
         close(iunit)

         ! Reopen file and dump raw binary grid coordinates
         open(newunit=iunit,file=filename,status='old',form='unformatted',access='stream',position='append',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
         write(iunit) '<AppendedData encoding="raw">'//eol
         write(iunit) '_'
         data_size=(this%cfg%nx+1)*SP
         write(iunit) data_size
         write(iunit) real(this%cfg%x(this%cfg%imin:this%cfg%imax+1),SP)
         data_size=(this%cfg%ny+1)*SP
         write(iunit) data_size
         write(iunit) real(this%cfg%y(this%cfg%jmin:this%cfg%jmax+1),SP)
         data_size=(this%cfg%nz+1)*SP
         write(iunit) data_size
         write(iunit) real(this%cfg%z(this%cfg%kmin:this%cfg%kmax+1),SP)
         close(iunit)
      end if

      ! Traverse all datasets and print them all out - scalars first
      my_scl=>this%first_scl
      do while (associated(my_scl))
         ! Write size of array
         if (this%cfg%amRoot) then
            open(newunit=iunit,file=filename,status='old',form='unformatted',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
            data_size=(this%cfg%nx*this%cfg%ny*this%cfg%nz)*SP
            write(iunit) data_size
            close(iunit)
         end if
         
         ! Now parallel-write the actual data (note that we allow both real and integer fields!)
         call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[vtk write data] Problem encountered while parallel writing data file '//trim(filename))
         call MPI_FILE_GET_POSITION(ifile,disp,ierr)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
         if (associated(my_scl%rptr)) spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_scl%rptr(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         if (associated(my_scl%iptr)) spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_scl%iptr(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
         call MPI_FILE_CLOSE(ifile,ierr)

         ! Continue on to the next scalar object
         my_scl=>my_scl%next
      end do

      ! Traverse all datasets and print them all out - vectors second
      my_vct=>this%first_vct
      do while (associated(my_vct))
         ! Write size of array
         if (this%cfg%amRoot) then
            open(newunit=iunit,file=filename,status='old',form='unformatted',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
            data_size=3*(this%cfg%nx*this%cfg%ny*this%cfg%nz)*SP
            write(iunit) data_size
            close(iunit)
         end if
         
         ! Now parallel-write the actual data
         call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[vtk write data] Problem encountered while parallel writing data file '//trim(filename))
         call MPI_FILE_GET_POSITION(ifile,disp,ierr)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SP3view,'native',info_mpiio,ierr)
         sp3buff(1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptrx(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         sp3buff(2,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptry(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         sp3buff(3,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptrz(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         call MPI_FILE_WRITE_ALL(ifile,sp3buff,3*this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
         call MPI_FILE_CLOSE(ifile,ierr)

         ! Continue on to the next scalar object
         my_vct=>my_vct%next
      end do

      if (this%cfg%amRoot) then
         ! Write general ASCII header for timestep VTU file
         open(newunit=iunit,file=filename,status='old',form='unformatted',access='stream',position='append',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
         write(iunit) eol//'</AppendedData>'//eol
         write(iunit) '</VTKFile>'//eol
         close(iunit)
      end if


      ! Get rid of the SP buffer
      if (associated(this%first_scl)) deallocate(spbuff)
      if (associated(this%first_vct)) deallocate(sp3buff)
      
      ! Finally, re-write the case file
      call this%write_case()
      
      ! Now output all surface meshes
      my_srf=>this%first_srf
      do while (associated(my_srf))
         ! Output the surface mesh as a distinct directory
         call this%write_surf(my_srf)
         ! Continue on to the next surface mesh object
         my_srf=>my_srf%next
      end do
      
      ! Now output all particle meshes
      my_prt=>this%first_prt
      do while (associated(my_prt))
         ! Output the particle mesh as a distinct directory
         call this%write_part(my_prt)
         ! Continue on to the next particle mesh object
         my_prt=>my_prt%next
      end do
      
   end subroutine write_data

   
   !> Case description serial output to a text file
   subroutine write_case(this)
      use messager,  only: die
      implicit none
      class(vtk), intent(in) :: this
      integer :: iunit,ierr,n
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      type(srf), pointer :: my_srf
      type(prt), pointer :: my_prt
      character(len=str_medium) :: filename
      character(len=1), parameter :: eol = char(10)
      character(len=2), parameter :: idt = '  '

      ! Only the root does this work
      if (.not.this%cfg%amRoot) return
      
      ! Open the case file
      open(newunit=iunit,file='vtk/'//trim(this%name)//'/nga.pvd',form='formatted',status='replace',access='stream',iostat=ierr)
      if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//'vtk/'//trim(this%name)//'/nga.pvd')
      write(iunit,'(a)')   '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write(iunit,'(a,i6.6,a)')   idt//'<Collection nsteps="',this%ntime,'">'
      do n=1,this%ntime
         ! Add timestep data file
         filename='data/data.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') n
         filename=trim(filename)//'.vtr'
         write(iunit,'(a,es12.5,a)') idt//idt//'<DataSet timestep="',this%time(n),'" file="'//trim(filename)//'"/>'
      end do
      write(iunit,'(a)')   idt//'</Collection>'
      write(iunit,'(a)')   '</VTKFile>'
      close(iunit)

      ! Add surfaces
      my_srf=>this%first_srf
      do while (associated(my_srf))
         ! Open the case file
         open(newunit=iunit,file='vtk/'//trim(this%name)//'/'//trim(my_srf%name)//'.pvd',form='formatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//'vtk/'//trim(this%name)//'/nga.pvd')
         write(iunit,'(a)')   '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
         write(iunit,'(a,i6.6,a)')   idt//'<Collection nsteps="',this%ntime,'">'
         do n=1,this%ntime
            ! Write timestep surface mesh
            filename=trim(my_srf%name)//'/'//trim(my_srf%name)//'.'
            write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') n
            filename=trim(filename)//'.vtu'       
            write(iunit,'(a,es12.5,a)') idt//idt//'<DataSet timestep="',this%time(n),'" file="'//trim(filename)//'"/>'
         end do
         write(iunit,'(a)')   idt//'</Collection>'
         write(iunit,'(a)')   '</VTKFile>'
         close(iunit)
         ! Continue on to the next surface mesh object
         my_srf=>my_srf%next
      end do

      ! Add particles
      my_prt=>this%first_prt
      do while (associated(my_prt))
         ! Open the case file
         open(newunit=iunit,file='vtk/'//trim(this%name)//'/'//trim(my_prt%name)//'.pvd',form='formatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//'vtk/'//trim(this%name)//'/nga.pvd')
         write(iunit,'(a)')   '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
         write(iunit,'(a,i6.6,a)')   idt//'<Collection nsteps="',this%ntime,'">'
         do n=1,this%ntime
            ! Write timestep particle file
            filename=trim(my_prt%name)//'/'//trim(my_prt%name)//'.'
            write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') n
            filename=trim(filename)//'.vtp'       
            write(iunit,'(a,es12.5,a)') idt//idt//'<DataSet timestep="',this%time(n),'" file="'//trim(filename)//'"/>'
         end do
         write(iunit,'(a)')   idt//'</Collection>'
         write(iunit,'(a)')   '</VTKFile>'
         close(iunit)
         ! Continue on to the next particle mesh object
         my_prt=>my_prt%next
      end do

   end subroutine write_case
   
   !> Procedure that writes out a surface mesh in Ensight format
   subroutine write_surf(this,surf)
      use messager,  only: die
      use mpi_f08,   only: MPI_BARRIER,MPI_BCAST,MPI_INTEGER4
      use precision, only: SP,DP,I4,I8

      implicit none
      class(vtk), intent(in) :: this
      type(srf), pointer, intent(in) :: surf
      integer     :: iunit,ierr,rank,n
      integer     :: nvert,ntri,npoly,nconn_tri,nconn_poly,offset
      integer(I8) :: buffer_I8
      integer(I4) :: data_size
      integer(I4) :: VTK_POLYGON = 7
      integer(I4) :: VTK_BEZIER_TRIANGLE = 76
      character(len=1), parameter :: eol = char(10)
      character(len=2), parameter :: idt = '  '
      character(len=str_medium) :: filename
      character(len=50) :: time
      nvert      = surf%ptr%nVert
      ntri       = surf%ptr%nBezierTri
      npoly      = surf%ptr%nPoly
      nconn_tri  = size(surf%ptr%bezierTriConn)
      nconn_poly = size(surf%ptr%polyConn)
      offset     = 0
      
      filename='vtk/'//trim(this%name)//'/'//trim(surf%name)//'/'//trim(surf%name)//'.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      filename=trim(filename)//'.vtu'

      if (this%cfg%amRoot) then
         ! Write general ASCII header for timestep VTU file
         open(newunit=iunit,file=filename,status='replace',form='formatted',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
         write(iunit,'(a)')             '<?xml version="1.0"?>'
         write(iunit,'(a)')             '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
         write(iunit,'(a)')             idt//'<UnstructuredGrid>'
         ! Close the file
         close(iunit)
      end if
      
      ! Write ASCII header for local piece
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            open(newunit=iunit,file=filename,status='old',form='formatted',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//'<Piece NumberOfPoints="',nvert,'" NumberOfCells ="',ntri+npoly,'">'
            write(iunit,'(a)')             idt//idt//idt//'<Points>'
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Float',SP*8,'" NumberOfComponents="3" format="appended" offset="',offset,'">'
            write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
            offset=offset+I4+nvert*3*SP
            write(iunit,'(a)')             idt//idt//idt//'</Points>'
            write(iunit,'(a)')             idt//idt//idt//'<PointData RationalWeights="RationalWeights">'
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="RationalWeights" format="appended" offset="',offset,'">'
            write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
            offset=offset+I4+nvert*SP
            write(iunit,'(a)')             idt//idt//idt//'</PointData>'
            write(iunit,'(a)')             idt//idt//idt//'<Cells>'
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Int',I8*8,'" Name="connectivity" format="appended" offset="',offset,'">'
            write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
            offset=offset+I4+nconn_tri*I8+nconn_poly*I8
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Int',I8*8,'" Name="offsets" format="appended" offset="',offset,'">'
            write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
            offset=offset+I4+ntri*I8+npoly*I8
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Int',I4*8,'" Name="types" format="appended" offset="',offset,'">'
            write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
            offset=offset+I4+ntri*I4+npoly*I4
            write(iunit,'(a)')             idt//idt//idt//'</Cells>'
            if (surf%ptr%nvar.gt.0) then
               write(iunit,'(a)')          idt//idt//idt//'<CellData>'
               do n=1,surf%ptr%nvar
                  write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="'//trim(surf%ptr%varname(n))//'" format="appended" offset="',offset,'">'
                  write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
                  offset=offset+I4+ntri*SP+npoly*SP
               end do
               write(iunit,'(a)')          idt//idt//idt//'</CellData>'
            end if
            write(iunit,'(a)')             idt//idt//'</Piece>'
            if (rank.eq.this%cfg%nproc-1) then
               write(iunit,'(a)')          idt//'</UnstructuredGrid>'
            end if
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BCAST(offset,1,MPI_INTEGER4,rank,this%cfg%comm,ierr)
      end do

      ! Write binary data
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            open(newunit=iunit,file=filename,status='old',form='unformatted',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
            if (rank.eq.0) then
               write(iunit) '<AppendedData encoding="raw">'//eol
               write(iunit) '_'
            end if
            data_size=nvert*3*SP
            write(iunit) data_size
            do n=1,nvert
               write(iunit) real(surf%ptr%xVert(n),SP),real(surf%ptr%yVert(n),SP),real(surf%ptr%zVert(n),SP)
            end do
            data_size=nvert*SP
            write(iunit) data_size
            write(iunit) real(surf%ptr%wVert,SP)
            data_size=nconn_tri*I8+nconn_poly*I8
            write(iunit) data_size
            write(iunit) int(surf%ptr%bezierTriConn,I8)
            write(iunit) int(surf%ptr%polyConn,I8)
            data_size=ntri*I8+npoly*I8
            write(iunit) data_size
            do n=1,ntri
               buffer_I8=n*6
               write(iunit) buffer_I8
            end do
            buffer_I8=6*ntri
            do n=1,npoly
               buffer_I8=buffer_I8+surf%ptr%polySize(n)
               write(iunit) buffer_I8
            end do
            data_size=ntri*I4+npoly*I4
            write(iunit) data_size
            do n=1,ntri
               write(iunit) VTK_BEZIER_TRIANGLE
            end do
            do n=1,npoly
               write(iunit) VTK_POLYGON
            end do
            if (surf%ptr%nvar.gt.0) then
                do n=1,surf%ptr%nvar
                  data_size=ntri*SP+npoly*SP
                  write(iunit) data_size
                  write(iunit) real(surf%ptr%var(n,:),SP)
               end do
            end if
            if (rank.eq.this%cfg%nproc-1) then
               write(iunit) eol//'</AppendedData>'//eol
               write(iunit) '</VTKFile>'//eol
            end if
            close(iunit)
            end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do

   end subroutine write_surf
   
   
   !> Procedure that writes out a particle mesh in Ensight format
   subroutine write_part(this,part)
      use precision, only: SP,I4
      use messager,  only: die
      use mpi_f08,   only: MPI_BARRIER,MPI_INTEGER4
      implicit none
      class(vtk), intent(in) :: this
      type(prt), pointer, intent(in) :: part
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,n,offset
      integer(I4) :: data_size
      character(len=1), parameter :: eol = char(10)
      character(len=2), parameter :: idt = '  '

      filename='vtk/'//trim(this%name)//'/'//trim(part%name)//'/'//trim(part%name)//'.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      filename=trim(filename)//'.vtp'

      if (this%cfg%amRoot) then
         ! Write general ASCII header for timestep VTU file
         open(newunit=iunit,file=filename,status='replace',form='formatted',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
         write(iunit,'(a)')             '<?xml version="1.0"?>'
         write(iunit,'(a)')             '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
         write(iunit,'(a)')             idt//'<PolyData>'
         ! Close the file
         close(iunit)
      end if

      ! Write ASCII header for local piece
      offset=0
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            open(newunit=iunit,file=filename,status='old',form='formatted',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
            write(iunit,'(a,i0,a)')        idt//idt//'<Piece NumberOfPoints="',part%ptr%n,'">'
            write(iunit,'(a)')             idt//idt//idt//'<Points>'
            write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Float',SP*8,'" NumberOfComponents="3" format="appended" offset="',offset,'">'
            write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
            offset=offset+I4+part%ptr%n*3*SP
            write(iunit,'(a)')             idt//idt//idt//'</Points>'
            write(iunit,'(a)')             idt//idt//idt//'<PointData>'
            do n=1,part%ptr%nvar
               write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="'//trim(part%ptr%varname(n))//'" format="appended" offset="',offset,'">'
               write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
               offset=offset+I4+part%ptr%n*SP
            end do
            do n=1,part%ptr%nvec
               write(iunit,'(a,i0,a,i0,a)')   idt//idt//idt//idt//'<DataArray type="Float',SP*8,'" Name="'//trim(part%ptr%vecname(n))//'" NumberOfComponents="3" format="appended" offset="',offset,'">'
               write(iunit,'(a)')             idt//idt//idt//idt//'</DataArray>'
               offset=offset+I4+3*part%ptr%n*SP
            end do
            write(iunit,'(a)')             idt//idt//idt//'</PointData>'
            write(iunit,'(a)')             idt//idt//'</Piece>'
            if (rank.eq.this%cfg%nproc-1) then
               write(iunit,'(a)')          idt//'</PolyData>'
            end if
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BCAST(offset,1,MPI_INTEGER4,rank,this%cfg%comm,ierr)
      end do

      ! Write binary data
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            open(newunit=iunit,file=filename,status='old',form='unformatted',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[vtk write surf] Could not open file: '//trim(filename))
            if (rank.eq.0) then
               write(iunit) '<AppendedData encoding="raw">'//eol
               write(iunit) '_'
            end if
            data_size=part%ptr%n*3*SP
            write(iunit) data_size
            if (part%ptr%n.gt.0) write(iunit) real(part%ptr%pos,SP)
            do n=1,part%ptr%nvar
               data_size=part%ptr%n*SP
               write(iunit) data_size
               if (part%ptr%n.gt.0) write(iunit) real(part%ptr%var(n,:),SP)
            end do
            do n=1,part%ptr%nvec
               data_size=3*part%ptr%n*SP
               write(iunit) data_size
               if (part%ptr%n.gt.0) write(iunit) real(part%ptr%vec(:,n,:),SP)
            end do
            if (rank.eq.this%cfg%nproc-1) then
               write(iunit) eol//'</AppendedData>'//eol
               write(iunit) '</VTKFile>'//eol
            end if
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do

   end subroutine write_part
   
end module vtk_class
