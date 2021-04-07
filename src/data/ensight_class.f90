!> Ensight class concept is defined here: given a config object,
!> it provides parallel I/O access to an ensight file
module ensight_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype
   use surfmesh_class, only: surfmesh
   use partmesh_class, only: partmesh
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ensight
   
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
   type :: ensight
      ! An ensight object has a name
      character(len=str_medium) :: name                               !< Name of ensight directory to read/write
      ! An ensight object stores time values
      integer :: ntime                                                !< Number of scalar values
      real(WP), dimension(:), allocatable :: time                     !< Time values
      ! An ensight object stores geometry data
      type(config), pointer :: cfg                                    !< Config for ensight geometry and parallel I/O
      ! An ensight object stores lists of pointers to data
      type(scl), pointer :: first_scl                                 !< Scalar list
      type(vct), pointer :: first_vct                                 !< Vector list
      type(srf), pointer :: first_srf                                 !< Surface list
      type(prt), pointer :: first_prt                                 !< Particle list
   contains
      procedure :: write_geom                                         !< Write out geometry
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
   end type ensight
   
   
   !> Declare ensight constructor
   interface ensight
      procedure construct_ensight
   end interface ensight
   
   
contains
   
   !> Constructor for an empty ensight object
   function construct_ensight(cfg,name) result(self)
      use messager, only: die
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      type(ensight) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      character(len=str_medium) :: line
      integer :: iunit,ierr,stat
      logical :: file_is_there,found
      
      ! Link to config
      self%cfg=>cfg
      
      ! Store casename
      self%name=trim(adjustl(name))
      
      ! Start with no time stamps
      self%ntime=0
      
      ! Create directory
      if (self%cfg%amRoot) then
         call execute_command_line('mkdir -p ensight')
         call execute_command_line('mkdir -p ensight/'//trim(self%name))
      end if
      
      ! Write out the geometry
      call self%write_geom(cfg=self%cfg,name='geometry')
      
      ! Empty pointer to lists for now
      self%first_scl=>NULL()
      self%first_vct=>NULL()
      self%first_srf=>NULL()
      self%first_prt=>NULL()
      
      ! Check if a case file exists already - root only
      if (self%cfg%amRoot) then
         inquire(file='ensight/'//trim(self%name)//'/nga.case',exist=file_is_there)
         if (file_is_there) then
            ! Open the case file
            open(newunit=iunit,file='ensight/'//trim(self%name)//'/nga.case',form='formatted',status='old',access='stream',iostat=ierr)
            ! Read lines until we find time values section
            stat=0; found=.false.
            do while (.not.found.and..not.is_iostat_end(stat))
               read(iunit,'(a)',iostat=stat) line
               if (line(1:16).eq.'number of steps:') found=.true.
            end do
            if (found) read(line(17:),'(i6)') self%ntime
            ! Now read the time values
            if (self%ntime.gt.0) then
               allocate(self%time(self%ntime))
               found=.false.
               do while (.not.found)
                  read(iunit,'(a)') line
                  if (line(1:12).eq.'time values:') found=.true.
               end do
               if (found) read(iunit,'(999999(es12.5,/))') self%time
            end if
            ! Close the case file
            close(iunit)
         else
            ! Output a barebone case file with only the geometry
            open(newunit=iunit,file='ensight/'//trim(self%name)//'/nga.case',form='formatted',status='replace',access='stream',iostat=ierr)
            write(iunit,'(4(a,/))') 'FORMAT','type: ensight gold','GEOMETRY','model: geometry'
            close(iunit)
         end if
      end if
      
      ! Communicate to all processors
      call MPI_BCAST(self%ntime,1,MPI_INTEGER,0,self%cfg%comm,ierr)
      if (self%ntime.gt.0) then
         if (.not.self%cfg%amRoot) allocate(self%time(self%ntime))
         call MPI_BCAST(self%time,self%ntime,MPI_REAL_WP,0,self%cfg%comm,ierr)
      end if
      
   end function construct_ensight
   
   
   !> Add a real scalar field for output
   subroutine add_rscalar(this,name,scalar)
      implicit none
      class(ensight), intent(inout) :: this
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
      ! Also create the corresponding directory
      if (this%cfg%amRoot) call execute_command_line('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_scl%name))
   end subroutine add_rscalar
   
   
   !> Add an integer scalar field for output
   subroutine add_iscalar(this,name,scalar)
      implicit none
      class(ensight), intent(inout) :: this
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
      ! Also create the corresponding directory
      if (this%cfg%amRoot) call execute_command_line('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_scl%name))
   end subroutine add_iscalar
   
   
   !> Add a vector field for output
   subroutine add_vector(this,name,vectx,vecty,vectz)
      implicit none
      class(ensight), intent(inout) :: this
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
      ! Also create the corresponding directory
      if (this%cfg%amRoot) call execute_command_line('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_vct%name))
   end subroutine add_vector
   
   
   !> Add a surface mesh for output
   subroutine add_surface(this,name,surface)
      implicit none
      class(ensight), intent(inout) :: this
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
      if (this%cfg%amRoot) call execute_command_line('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_srf%name))
   end subroutine add_surface
   
   
   !> Add a particle mesh for output
   subroutine add_particle(this,name,particle)
      implicit none
      class(ensight), intent(inout) :: this
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
      if (this%cfg%amRoot) call execute_command_line('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_prt%name))
   end subroutine add_particle
   
   
   !> Output all data in the object
   subroutine write_data(this,time)
      use precision, only: SP
      use messager,  only: die
      use parallel,  only: info_mpiio,MPI_REAL_SP
      use mpi_f08
      implicit none
      class(ensight), intent(inout) :: this
      real(WP), intent(in) :: time
      character(len=str_medium) :: filename
      integer :: iunit,ierr,n,i
      integer :: ibuff
      character(len=80) :: cbuff
      type(MPI_File) :: ifile
      integer(kind=MPI_OFFSET_KIND) :: disp
      type(MPI_Status):: status
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      type(srf), pointer :: my_srf
      type(prt), pointer :: my_prt
      real(SP), dimension(:,:,:), allocatable :: spbuff
      real(WP), dimension(:), allocatable :: temp_time
      
      ! Check provided time stamp and decide what to do
      if (this%ntime.eq.0) then
         ! First time stamp
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(this%ntime))
         this%time(1)=time
      else
         ! There are time stamps already, check where to insert
         n=this%ntime+1
         do i=this%ntime,1,-1
            if (time.le.this%time(i)) n=n-1
         end do
         this%ntime=n; allocate(temp_time(1:this%ntime))
         temp_time=[this%time(1:this%ntime-1),time]
         call move_alloc(temp_time,this%time)
      end if
      
      ! Prepare the SP buffer
      allocate(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))
      
      ! Traverse all datasets and print them all out - scalars first
      my_scl=>this%first_scl
      do while (associated(my_scl))
         
         ! Create filename
         filename='ensight/'//trim(this%name)//'/'//trim(my_scl%name)//'/'//trim(my_scl%name)//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         
         ! Root process starts writing the file header
         if (this%cfg%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write data] Could not open file '//trim(filename))
            ! Write the header
            cbuff=trim(my_scl%name); write(iunit) cbuff
            cbuff='part'           ; write(iunit) cbuff
            ibuff=1                ; write(iunit) ibuff
            cbuff='block'          ; write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         
         ! Now parallel-write the actual data (note that we allow both real and integer fields!)
         call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[ensight write data] Problem encountered while parallel writing data file '//trim(filename))
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
         
         ! Create filename
         filename='ensight/'//trim(this%name)//'/'//trim(my_vct%name)//'/'//trim(my_vct%name)//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         
         ! Root process starts writing the file header
         if (this%cfg%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write data] Could not open file '//trim(filename))
            ! Write the header
            cbuff=trim(my_vct%name); write(iunit) cbuff
            cbuff='part'           ; write(iunit) cbuff
            ibuff=1                ; write(iunit) ibuff
            cbuff='block'          ; write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         
         ! Now parallel-write the actual data
         call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[ensight write data] Problem encountered while parallel writing data file '//trim(filename))
         call MPI_FILE_GET_POSITION(ifile,disp,ierr)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
         spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptrx(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
         disp=disp+int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)*int(SP,MPI_OFFSET_KIND)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
         spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptry(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
         disp=disp+int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)*int(SP,MPI_OFFSET_KIND)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
         spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptrz(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
         call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
         call MPI_FILE_CLOSE(ifile,ierr)
         
         ! Continue on to the next vector object
         my_vct=>my_vct%next
         
      end do
      
      ! Get rid of the SP buffer
      deallocate(spbuff)
      
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
         ! Continue on to the next surface mesh object
         my_prt=>my_prt%next
      end do
      
   end subroutine write_data
   
   
   !> Case description serial output to a text file
   subroutine write_case(this)
      implicit none
      class(ensight), intent(in) :: this
      integer :: iunit,ierr
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      
      ! Only the root does this work
      if (.not.this%cfg%amRoot) return
      
      ! Open the case file
      open(newunit=iunit,file='ensight/'//trim(this%name)//'/nga.case',form='formatted',status='replace',access='stream',iostat=ierr)
      
      ! Write all the geometry information
      write(iunit,'(a,/,a,/,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: geometry'
      
      ! Write all the variable information
      write(iunit,'(a)') 'VARIABLE'
      write(iunit,'(a)') 'scalar per element: wall geometry.wall'
      my_scl=>this%first_scl
      do while (associated(my_scl))
         write(iunit,'(a)') 'scalar per element: 1 '//trim(my_scl%name)//' '//trim(my_scl%name)//'/'//trim(my_scl%name)//'.******'
         my_scl=>my_scl%next
      end do
      my_vct=>this%first_vct
      do while (associated(my_vct))
         write(iunit,'(a)') 'vector per element: 1 '//trim(my_vct%name)//' '//trim(my_vct%name)//'/'//trim(my_vct%name)//'.******'
         my_vct=>my_vct%next
      end do
      
      ! Write the time information
      write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
      write(iunit,'(999999(es12.5,/))') this%time
      
      ! Close the case file
      close(iunit)
      
   end subroutine write_case
   
   
   !> Geometry output to a file in parallel
   subroutine write_geom(this,cfg,name)
      use precision, only: SP
      use messager,  only: die
      use parallel,  only: info_mpiio,MPI_REAL_SP
      use mpi_f08
      implicit none
      class(ensight),   intent(in) :: this
      class(config),    intent(in) :: cfg
      character(len=*), intent(in) :: name
      integer :: iunit,ierr
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: disp
      real(SP), dimension(:,:,:), allocatable :: spbuff
      
      ! Only cfg root does geometry I/O
      if (cfg%amRoot) then
         
         ! First create a new file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(name),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write geom] Could not open file: ensight/'//trim(this%name)//'/'//trim(name))
         
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'        ; write(iunit) cbuff
         cbuff=trim(adjustl(cfg%name))             ; write(iunit) cbuff
         cbuff='node id off'                       ; write(iunit) cbuff
         cbuff='element id off'                    ; write(iunit) cbuff
         
         ! Extents
         cbuff='extents'                           ; write(iunit) cbuff
         rbuff=real(cfg%x(cfg%imin  ),SP)          ; write(iunit) rbuff
         rbuff=real(cfg%x(cfg%imax+1),SP)          ; write(iunit) rbuff
         rbuff=real(cfg%y(cfg%jmin  ),SP)          ; write(iunit) rbuff
         rbuff=real(cfg%y(cfg%jmax+1),SP)          ; write(iunit) rbuff
         rbuff=real(cfg%z(cfg%kmin  ),SP)          ; write(iunit) rbuff
         rbuff=real(cfg%z(cfg%kmax+1),SP)          ; write(iunit) rbuff
         
         ! Part header
         cbuff='part'                              ; write(iunit) cbuff
         ibuff=1                                   ; write(iunit) ibuff
         cbuff='Complete geometry'                 ; write(iunit) cbuff  ! We have a single grid-cfg here
         cbuff='block rectilinear'                 ; write(iunit) cbuff  ! Not blanked
         
         ! Number of cells
         ibuff=cfg%nx+1                            ; write(iunit) ibuff
         ibuff=cfg%ny+1                            ; write(iunit) ibuff
         ibuff=cfg%nz+1                            ; write(iunit) ibuff
         
         ! Mesh
         write(iunit) real(cfg%x(cfg%imin:cfg%imax+1),SP)
         write(iunit) real(cfg%y(cfg%jmin:cfg%jmax+1),SP)
         write(iunit) real(cfg%z(cfg%kmin:cfg%kmax+1),SP)
         
         ! Close the file
         close(iunit)
         
      end if
      
      ! Root process starts writing the file header for wall data
      if (cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(name)//'.wall',form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write data] Could not open file: '//'ensight/'//trim(this%name)//'/'//trim(name)//'.wall')
         ! Write the header
         cbuff='wall'           ; write(iunit) cbuff
         cbuff='part'           ; write(iunit) cbuff
         ibuff=1                ; write(iunit) ibuff
         cbuff='block'          ; write(iunit) cbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Prepare the SP buffer
      allocate(spbuff(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_))
      
      ! Now parallel-write the wall data
      call MPI_FILE_OPEN(cfg%comm,'ensight/'//trim(this%name)//'/'//trim(name)//'.wall',IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[ensight write geom] Problem encountered while parallel writing wall data file: '//'ensight/'//trim(this%name)//'/'//trim(name)//'.wall')
      call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,cfg%SPview,'native',info_mpiio,ierr)
      spbuff(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)=real(cfg%VF(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_),SP)
      call MPI_FILE_WRITE_ALL(ifile,spbuff,cfg%nx_*cfg%ny_*cfg%nz_,MPI_REAL_SP,status,ierr)
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Deallocate SP buffer
      deallocate(spbuff)
      
   end subroutine write_geom
   
   
   !> Procedure that writes out a surface mesh in Ensight format
   subroutine write_surf(this,surf)
      use precision, only: SP
      use messager,  only: die
      use mpi_f08,   only: mpi_barrier
      implicit none
      class(ensight), intent(in) :: this
      type(srf), pointer, intent(in) :: surf
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,n
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff
      
      ! Write the case file from scratch in ASCII format
      if (this%cfg%amRoot) then
         ! Open the case file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(surf%name)//'.case',form='formatted',status='replace',access='stream',iostat=ierr)
         ! Write all the geometry information
         write(iunit,'(a,/,a,/,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: 1 '//trim(surf%name)//'/'//trim(surf%name)//'.******'
         ! Write the variables
         write(iunit,'(a)') 'VARIABLE'
         do n=1,surf%ptr%nvar
            write(iunit,'(a)') 'scalar per element: 1 '//trim(surf%ptr%varname(n))//' '//trim(surf%name)//'/'//trim(surf%ptr%varname(n))//'.******'
         end do
         ! Write the time information
         write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
         write(iunit,'(999999(es12.5,/))') this%time
         ! Close the case file
         close(iunit)
      end if
      
      ! Generate the surface geometry filename
      filename='ensight/'//trim(this%name)//'/'//trim(surf%name)//'/'//trim(surf%name)//'.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      
      ! Write the file header for Ensight Gold unstructured geometry
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'        ; write(iunit) cbuff
         cbuff=trim(adjustl(surf%ptr%name))        ; write(iunit) cbuff
         cbuff='node id off'                       ; write(iunit) cbuff
         cbuff='element id off'                    ; write(iunit) cbuff
         ! Extents
         cbuff='extents'                           ; write(iunit) cbuff
         rbuff=real(this%cfg%x(this%cfg%imin  ),SP); write(iunit) rbuff
         rbuff=real(this%cfg%x(this%cfg%imax+1),SP); write(iunit) rbuff
         rbuff=real(this%cfg%y(this%cfg%jmin  ),SP); write(iunit) rbuff
         rbuff=real(this%cfg%y(this%cfg%jmax+1),SP); write(iunit) rbuff
         rbuff=real(this%cfg%z(this%cfg%kmin  ),SP); write(iunit) rbuff
         rbuff=real(this%cfg%z(this%cfg%kmax+1),SP); write(iunit) rbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Write polygonal mesh in Ensight Gold 'nsided' format
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
            ! Part header
            cbuff='part'                              ; write(iunit) cbuff
            ibuff=rank+1                              ; write(iunit) ibuff
            cbuff='Surface geometry per processor #'
            write(cbuff(len_trim(cbuff)+1:len_trim(cbuff)+6),'(i6.6)') this%cfg%rank
            write(iunit) cbuff
            ! Write part info if it exists on the processor
            if (surf%ptr%nPoly.gt.0) then
               ! Write out vertices
               cbuff='coordinates'                       ; write(iunit) cbuff
               ibuff=surf%ptr%nVert                      ; write(iunit) ibuff
               write(iunit) real(surf%ptr%xVert,SP)
               write(iunit) real(surf%ptr%yVert,SP)
               write(iunit) real(surf%ptr%zVert,SP)
               ! Write out polygons
               cbuff='nsided'                            ; write(iunit) cbuff
               ibuff=surf%ptr%nPoly                      ; write(iunit) ibuff
               write(iunit) surf%ptr%polySize
               write(iunit) surf%ptr%polyConn
            end if
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do
      
      ! Generate the additional variable files
      do n=1,surf%ptr%nvar
         filename='ensight/'//trim(this%name)//'/'//trim(surf%name)//'/'//trim(surf%ptr%varname(n))//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         ! Root write the header
         if (this%cfg%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
            ! Write the header
            cbuff=trim(surf%name); write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         ! Write the surface variables
         do rank=0,this%cfg%nproc-1
            if (rank.eq.this%cfg%rank) then
               ! Open the file
               open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
               ! Part header
               cbuff='part'         ; write(iunit) cbuff
               ibuff=rank+1         ; write(iunit) ibuff
               ! Write surf info if it exists on the processor
               if (surf%ptr%nPoly.gt.0) then
                  cbuff='nsided'       ; write(iunit) cbuff
                  write(iunit) real(surf%ptr%var(n,:),SP)
               end if
               ! Close the file
               close(iunit)
            end if
            ! Force synchronization
            call MPI_BARRIER(this%cfg%comm,ierr)
         end do
      end do
      
   end subroutine write_surf
   
   
   !> Procedure that writes out a particle mesh in Ensight format
   subroutine write_part(this,part)
      use precision, only: SP
      use messager,  only: die
      use mpi_f08,   only: MPI_BARRIER,MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      implicit none
      class(ensight), intent(in) :: this
      type(prt), pointer, intent(in) :: part
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,n
      character(len=80) :: cbuff
      integer :: ibuff,npart
      
      ! Write the case file from scratch in ASCII format
      if (this%cfg%amRoot) then
         ! Open the case file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(part%name)//'.case',form='formatted',status='replace',access='stream',iostat=ierr)
         ! Write all the geometry information
         write(iunit,'(a,/,a,/,/,a,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: geometry','measured: 1 '//trim(part%name)//'/particle.******'
         ! Write the variables
         write(iunit,'(a)') 'VARIABLE'
         write(iunit,'(a)') 'scalar per element: wall geometry.wall'
         do n=1,part%ptr%nvar
            write(iunit,'(a)') 'scalar per measured node: 1 '//trim(part%ptr%varname(n))//' '//trim(part%name)//'/'//trim(part%ptr%varname(n))//'.******'
         end do
         ! Write the time information
         write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
         write(iunit,'(999999(es12.5,/))') this%time
         ! Close the case file
         close(iunit)
      end if
      
      ! Generate the particle geometry file
      filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/particle.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      ! Count the number of particles
      call MPI_ALLREDUCE(part%ptr%n,npart,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
      ! Root write the header
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff=trim(adjustl(part%ptr%name))        ; write(iunit) cbuff
         cbuff='particle coordinates'              ; write(iunit) cbuff
         ibuff=npart                               ; write(iunit) ibuff
         write(iunit) (ibuff,ibuff=1,npart)
         ! Close the file
         close(iunit)
      end if
      ! Write the particle coordinates
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
            ! Write part info if it exists on the processor
            if (part%ptr%n.gt.0) write(iunit) real(part%ptr%pos,SP)
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do
      
      ! Generate the additional variable files
      do n=1,part%ptr%nvar
         filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/'//trim(part%ptr%varname(n))//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         ! Root write the header
         if (this%cfg%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
            ! General header
            cbuff='particle '//trim(part%ptr%varname(n)); write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         ! Write the particle variable
         do rank=0,this%cfg%nproc-1
            if (rank.eq.this%cfg%rank) then
               ! Open the file
               open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
               ! Write part info if it exists on the processor
               if (part%ptr%n.gt.0) write(iunit) real(part%ptr%var(n,:),SP)
               ! Close the file
               close(iunit)
            end if
            ! Force synchronization
            call MPI_BARRIER(this%cfg%comm,ierr)
         end do
      end do
      
   end subroutine write_part
   
   
end module ensight_class
