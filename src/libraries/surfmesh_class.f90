module surfmesh_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: surfmesh
   
   !> Surface mesh object
   type :: surfmesh
      character(len=str_medium) :: name='UNNAMED_SURFMESH' !< Name for the surface mesh
      integer :: nVert                                     !< Number of vertices
      real(WP), dimension(:), allocatable :: xVert         !< X position of the vertices - size=nVert
      real(WP), dimension(:), allocatable :: yVert         !< Y position of the vertices - size=nVert
      real(WP), dimension(:), allocatable :: zVert         !< Z position of the vertices - size=nVert
      integer :: nPoly                                     !< Number of polygons
      integer,  dimension(:), allocatable :: polySize      !< Size of polygons - size=nPoly
      integer,  dimension(:), allocatable :: polyConn      !< Connectivity - size=sum(polySize)
      integer :: nvar=0                                                 !< Number of surface variables stored (default=0)
      real(WP), dimension(:,:), allocatable :: var                      !< Surface variable storage
      character(len=str_medium), dimension(:), allocatable :: varname   !< Name of surface variable fields
   contains
      procedure :: reset                                   !< Reset surface mesh to zero size
      procedure :: set_size                                !< Set surface mesh to provided size
      procedure :: finalize                                !< Finalize surfmesh object
      procedure :: add_polygon                             !< Append a polygon (grows arrays as needed)
   end type surfmesh
   
   
   !> Declare surface mesh constructor
   interface surfmesh
      module procedure construct_empty
      module procedure construct_from_ply_serial_read
      module procedure construct_from_ply_parallel_read
   end interface surfmesh
   
   
contains
   
   
   !> Constructor for surface mesh object from a .ply file - parallel read version
   function construct_from_ply_parallel_read(comm,plyfile,nvar,name) result(self)
      use parallel, only: info_mpiio,MPI_REAL_SP
      use messager, only: die
      use mpi_f08
      implicit none
      type(surfmesh) :: self
      type(MPI_Comm), intent(in) :: comm
      character(len=*), intent(in) :: plyfile
      integer, intent(in) :: nvar
      character(len=*), optional :: name
      type(MPI_File) :: ifile
      type(MPI_Status) :: status
      integer :: ierr
      
      ! Set the name of the surface mesh
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Default to 0 size
      self%nVert=0
      self%nPoly=0
      
      ! Initialize additional variables
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname='' !< Users will set the name themselves
      
      ! Open the ply file
      call MPI_FILE_OPEN(comm,trim(adjustl(plyfile)),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[surfmesh constructor from file] Could not open file: '//trim(plyfile))
      
      ! Read the ply header
      read_header: block
         use string, only: str_medium
         character(len=str_medium) :: line
         character :: cbuf
         integer :: i,sp1,sp2
         ! Read one line at a time
         line=''
         do while (trim(line).ne.'end_header')
            ! Prepare to read a new line
            i=0; line=''
            do
               call MPI_FILE_READ_ALL(ifile,cbuf,1,MPI_CHARACTER,status,ierr)
               if (cbuf.eq.new_line(cbuf)) then
                  exit
               else
                  i=i+1; line(i:i)=cbuf
               end if
            end do
            ! Understand the line
            sp1=scan(line,' ')
            if (line(1:sp1).eq.'element') then
               sp2=scan(line(sp1+1:),' '); sp2=sp1+sp2
               if (line(sp1+1:sp2).eq.'vertex') then
                  read(line(sp2+1:),*) self%nVert
               else if (line(sp1+1:sp2).eq.'face') then
                  read(line(sp2+1:),*) self%nPoly
               end if
            end if
         end do
      end block read_header
      
      ! Resize my surfmesh
      call self%set_size(nvert=self%nVert,npoly=self%nPoly)
      
      ! Read the ply vertices
      read_vertices: block
         use precision, only: SP
         !real(SP), dimension(4,self%nVert) :: myVert
         real(SP), dimension(3,self%nVert) :: myVert
         integer :: n
         call MPI_FILE_READ_ALL(ifile,myVert,3*self%nVert,MPI_REAL_SP,status,ierr)
         do n=1,self%nVert
            self%xVert(n)=myVert(1,n)
            self%yVert(n)=myVert(2,n)
            self%zVert(n)=myVert(3,n)
         end do
      end block read_vertices
      
      ! Read the ply faces
      read_faces: block
         !use precision, only: SP
         !real(SP) :: buffer
         character(1) :: ibuf
         integer :: np,current_size
         do np=1,self%nPoly
            ! Read polygon size
            call MPI_FILE_READ_ALL(ifile,ibuf,1,MPI_UNSIGNED_CHAR,status,ierr)
            self%polySize(np)=ichar(ibuf)
            ! Extend allocation of connectivity
            current_size=sum(self%polySize(1:np))
            call resize_polyConn(current_size)
            ! Read connectivity
            call MPI_FILE_READ_ALL(ifile,self%polyConn(current_size-self%polySize(np)+1:current_size),self%polySize(np),MPI_INTEGER,status,ierr)
         end do
         ! First vertex is number 1, not 0
         self%polyConn=self%polyConn+1
      end block read_faces
      
      ! Close the plyfile
      call MPI_FILE_CLOSE(ifile,ierr)
      
   contains
      
      !> Adaptation of polyConn array size
      subroutine resize_polyConn(n)
         implicit none
         integer, intent(in) :: n
         real(WP), parameter :: coeff_up=1.5_WP      !< Connectivity array size increase factor
         real(WP), parameter :: coeff_dn=0.7_WP      !< Connectivity array size decrease factor
         integer, dimension(:), allocatable :: tmp
         integer :: size_now,size_new
         ! Resize particle array to size n
         if (.not.allocated(self%polyConn)) then
            ! Allocate directly to size n
            allocate(self%polyConn(n))
         else
            ! Update from a non-zero size to another non-zero size
            size_now=size(self%polyConn,dim=1)
            if (n.gt.size_now) then
               size_new=max(n,int(real(size_now,WP)*coeff_up))
               allocate(tmp(size_new))
               tmp(1:size_now)=self%polyConn
               call move_alloc(tmp,self%polyConn)
            else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
               allocate(tmp(n))
               tmp(1:n)=self%polyConn(1:n)
               call move_alloc(tmp,self%polyConn)
            end if
         end if
      end subroutine resize_polyConn
      
   end function construct_from_ply_parallel_read
   

   !> Constructor for surface mesh object from a .ply file - serial read version
   function construct_from_ply_serial_read(plyfile,nvar,name) result(self)
      use messager, only: die
      implicit none
      type(surfmesh) :: self
      character(len=*), intent(in) :: plyfile
      integer, intent(in) :: nvar
      character(len=*), optional :: name
      integer :: iunit,ierr

      ! Set the name of the surface mesh
      if (present(name)) self%name=trim(adjustl(name))

      ! Default to 0 size
      self%nVert=0
      self%nPoly=0

      ! Initialize additional variables
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname='' !< Users will set the name themselves
      
      ! Open the ply file
      open(newunit=iunit,file=trim(adjustl(plyfile)),form='unformatted',status='old',access='stream',iostat=ierr)
      if (ierr.ne.0) call die('[surfmesh constructor from file] Could not open file: '//trim(plyfile))
      
      ! Read the ply header
      read_header: block
         use string, only: str_medium
         character(len=str_medium) :: line
         character :: cbuf
         integer :: i,sp1,sp2
         ! Read one line at a time
         line=''
         do while (trim(line).ne.'end_header')
            ! Prepare to read a new line
            i=0; line=''
            do
               read(iunit) cbuf
               if (cbuf.eq.new_line(cbuf)) then
                  exit
               else
                  i=i+1; line(i:i)=cbuf
               end if
            end do
            ! Understand the line
            sp1=scan(line,' ')
            if (line(1:sp1).eq.'element') then
               sp2=scan(line(sp1+1:),' '); sp2=sp1+sp2
               if (line(sp1+1:sp2).eq.'vertex') then
                  read(line(sp2+1:),*) self%nVert
               else if (line(sp1+1:sp2).eq.'face') then
                  read(line(sp2+1:),*) self%nPoly
               end if
            end if
         end do
      end block read_header
      
      ! Resize my surfmesh
      call self%set_size(nvert=self%nVert,npoly=self%nPoly)
      
      ! Read the ply vertices
      read_vertices: block
         use precision, only: SP
         !real(SP), dimension(4,self%nVert) :: myVert
         real(SP), dimension(3,self%nVert) :: myVert
         integer :: n
         read(iunit) myVert
         do n=1,self%nVert
            self%xVert(n)=myVert(1,n)
            self%yVert(n)=myVert(2,n)
            self%zVert(n)=myVert(3,n)
         end do
      end block read_vertices
      
      ! Read the ply faces
      read_faces: block
         !use precision, only: SP
         !real(SP) :: buffer
         character(1) :: ibuf
         integer :: np,current_size
         do np=1,self%nPoly
            ! Read polygon size
            read(iunit) ibuf
            self%polySize(np)=ichar(ibuf)
            ! Extend allocation of connectivity
            current_size=sum(self%polySize(1:np))
            call resize_polyConn(current_size)
            ! Read connectivity
            read(iunit) self%polyConn(current_size-self%polySize(np)+1:current_size)!,buffer
         end do
         ! First vertex is number 1, not 0
         self%polyConn=self%polyConn+1
      end block read_faces
      
      ! Close the plyfile
      close(iunit)
      
   contains
   
      !> Adaptation of polyConn array size
      subroutine resize_polyConn(n)
         implicit none
         integer, intent(in) :: n
         real(WP), parameter :: coeff_up=1.5_WP      !< Connectivity array size increase factor
         real(WP), parameter :: coeff_dn=0.7_WP      !< Connectivity array size decrease factor
         integer, dimension(:), allocatable :: tmp
         integer :: size_now,size_new
         ! Resize particle array to size n
         if (.not.allocated(self%polyConn)) then
            ! Allocate directly to size n
            allocate(self%polyConn(n))
         else
            ! Update from a non-zero size to another non-zero size
            size_now=size(self%polyConn,dim=1)
            if (n.gt.size_now) then
               size_new=max(n,int(real(size_now,WP)*coeff_up))
               allocate(tmp(size_new))
               tmp(1:size_now)=self%polyConn
               call move_alloc(tmp,self%polyConn)
            else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
               allocate(tmp(n))
               tmp(1:n)=self%polyConn(1:n)
               call move_alloc(tmp,self%polyConn)
            end if
         end if
      end subroutine resize_polyConn
      
   end function construct_from_ply_serial_read
   

   !> Constructor for an empty surface mesh object
   function construct_empty(nvar,name) result(self)
      implicit none
      type(surfmesh) :: self
      integer, intent(in) :: nvar
      character(len=*), optional :: name
      ! Set the name of the surface mesh
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 size
      self%nVert=0
      self%nPoly=0
      ! Initialize additional variables
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname='' !< Users will set the name themselves
   end function construct_empty
   
   
   !> Reset mesh storage
   subroutine reset(this)
      implicit none
      class(surfmesh), intent(inout) :: this
      this%nPoly=0; this%nVert=0
      if (allocated(this%xVert))    deallocate(this%xvert)
      if (allocated(this%yVert))    deallocate(this%yvert)
      if (allocated(this%zVert))    deallocate(this%zvert)
      if (allocated(this%polySize)) deallocate(this%polySize)
      if (allocated(this%polyConn)) deallocate(this%polyConn)
      if (allocated(this%var))      deallocate(this%var)
   end subroutine reset
   
   
   ! Set mesh storage size - leave connectivity alone
   subroutine set_size(this,nvert,npoly)
      implicit none
      class(surfmesh), intent(inout) :: this
      integer, intent(in) :: nvert,npoly
      this%nPoly=npoly; this%nVert=nvert
      allocate(this%xVert   (this%nVert))
      allocate(this%yVert   (this%nVert))
      allocate(this%zVert   (this%nVert))
      allocate(this%polySize(this%nPoly))
      allocate(this%var     (this%nvar,this%nPoly))
   end subroutine set_size
   
   
   !> Finalize surfmesh
   subroutine finalize(this)
      implicit none
      class(surfmesh), intent(inout) :: this
      call this%reset()
      this%nvar=0
      if (allocated(this%varname)) deallocate(this%varname)
      this%name='UNNAMED_SURFMESH'
   end subroutine finalize
   
   
   !> Append a single polygon to the surface mesh
   !> Grows internal arrays as needed using doubling strategy
   !> Optional vardata provides nvar surface variable values for the polygon
   subroutine add_polygon(this, nv, verts, vardata)
      implicit none
      class(surfmesh), intent(inout) :: this
      integer, intent(in) :: nv                       !< Number of vertices
      real(WP), dimension(3,nv), intent(in) :: verts  !< Polygon vertices (3 x nv)
      real(WP), dimension(this%nvar), intent(in), optional :: vardata !< Per-polygon surface variables
      real(WP), dimension(:), allocatable :: new_x, new_y, new_z
      real(WP), dimension(:,:), allocatable :: new_var
      integer, dimension(:), allocatable :: new_size, new_conn
      integer :: new_cap, i

      ! Initialize arrays if not allocated
      if (.not.allocated(this%xVert)) then
         allocate(this%xVert(1000), this%yVert(1000), this%zVert(1000))
         allocate(this%polySize(1000), this%polyConn(1000))
         this%nVert = 0
         this%nPoly = 0
      end if
      if (this%nvar .gt. 0 .and. .not.allocated(this%var)) allocate(this%var(this%nvar, 1000))

      ! Grow vertex arrays if needed
      if (this%nVert + nv .gt. size(this%xVert)) then
         new_cap = max(2*size(this%xVert), this%nVert + nv)
         allocate(new_x(new_cap), new_y(new_cap), new_z(new_cap))
         new_x(1:this%nVert) = this%xVert(1:this%nVert)
         new_y(1:this%nVert) = this%yVert(1:this%nVert)
         new_z(1:this%nVert) = this%zVert(1:this%nVert)
         call move_alloc(new_x, this%xVert)
         call move_alloc(new_y, this%yVert)
         call move_alloc(new_z, this%zVert)
      end if

      ! Grow polyConn array if needed
      if (this%nVert + nv .gt. size(this%polyConn)) then
         new_cap = max(2*size(this%polyConn), this%nVert + nv)
         allocate(new_conn(new_cap))
         new_conn(1:this%nVert) = this%polyConn(1:this%nVert)
         call move_alloc(new_conn, this%polyConn)
      end if

      ! Grow polySize and var arrays if needed
      if (this%nPoly + 1 .gt. size(this%polySize)) then
         new_cap = 2*size(this%polySize)
         allocate(new_size(new_cap))
         new_size(1:this%nPoly) = this%polySize(1:this%nPoly)
         call move_alloc(new_size, this%polySize)
      end if
      if (this%nvar .gt. 0 .and. this%nPoly + 1 .gt. size(this%var, 2)) then
         new_cap = 2*size(this%var, 2)
         allocate(new_var(this%nvar, new_cap))
         new_var(:, 1:this%nPoly) = this%var(:, 1:this%nPoly)
         call move_alloc(new_var, this%var)
      end if

      ! Append polygon
      this%nPoly = this%nPoly + 1
      this%polySize(this%nPoly) = nv
      do i = 1, nv
         this%nVert = this%nVert + 1
         this%xVert(this%nVert) = verts(1, i)
         this%yVert(this%nVert) = verts(2, i)
         this%zVert(this%nVert) = verts(3, i)
         this%polyConn(this%nVert) = this%nVert
      end do
      ! Store per-polygon surface variables
      if (this%nvar .gt. 0) then
         if (present(vardata)) then
            this%var(:, this%nPoly) = vardata
         else
            this%var(:, this%nPoly) = 0.0_WP
         end if
      end if

   end subroutine add_polygon

end module surfmesh_class
