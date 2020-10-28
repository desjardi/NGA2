!> Ensight class concept is defined here: given a config object,
!> it provides parallel I/O access to an ensight file
module ensight_class
   use precision,    only: WP
   use string,       only: str_short,str_medium
   use config_class, only: config
   use mpi_f08,      only: MPI_Datatype
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ensight
   
   ! List types
   type :: scl
      type(scl), pointer :: next
      character(len=str_short) :: name
      real(WP), dimension(:,:,:), pointer :: ptr
   end type scl
   type :: vct
      type(vct), pointer :: next
      character(len=str_short) :: name
      real(WP), dimension(:,:,:), pointer :: ptrx
      real(WP), dimension(:,:,:), pointer :: ptry
      real(WP), dimension(:,:,:), pointer :: ptrz
   end type vct
   
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
   contains
      procedure :: write_geom                                         !< Write out geometry
      procedure :: write_data                                         !< Write out data
      procedure :: write_case                                         !< Write out case file
      procedure :: add_scalar                                         !< Add a new scalar field
      procedure :: add_vector                                         !< Add a new vector field
   end type ensight
   
   
   !> Declare ensight constructor
   interface ensight
      procedure construct_ensight
   end interface ensight
   
   
contains
   
   !> Constructor for an empty ensight object
   function construct_ensight(cfg,name) result(self)
      use monitor,  only: die
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
      call self%write_geom()
      
      ! Empty pointer to lists for now
      self%first_scl=>NULL()
      self%first_vct=>NULL()
      
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
   
   
   !> Add a scalar field for output
   subroutine add_scalar(this,name,scalar)
      implicit none
      class(ensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target :: scalar
      type(scl), pointer :: new_scl
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%ptr =>scalar
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
      ! Also create the corresponding directory
      if (this%cfg%amRoot) call execute_command_line('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_scl%name))
   end subroutine add_scalar
   
   
   !> Add a vector field for output
   subroutine add_vector(this,name,vectx,vecty,vectz)
      implicit none
      class(ensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target :: vectx
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target :: vecty
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target :: vectz
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
   
   
   !> Output all data in the object
   subroutine write_data(this,time)
      use precision, only: SP
      use monitor,   only: die
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
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(1:this%ntime)); this%time=temp_time
         deallocate(temp_time)
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
         
         ! Now parallel-write the actual data
         call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[ensight write data] Problem encountered while parallel writing data file '//trim(filename))
         call MPI_FILE_GET_POSITION(ifile,disp,ierr)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
         spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_scl%ptr(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
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
      
   end subroutine
   
   
   !> Geometry output to a file in parallel
   subroutine write_geom(this)
      use precision, only: SP
      use monitor,   only: die
      use parallel,  only: info_mpiio
      use mpi_f08
      implicit none
      class(ensight), intent(in) :: this
      integer :: iunit,ierr
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff
      type(MPI_Datatype) :: view
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: disp
      integer, dimension(:,:,:), allocatable :: iblank
      integer :: i,i1,i2,j,j1,j2,k,k1,k2,datasize
      integer, dimension(3) :: gsizes,lsizes,lstart
      
      ! Root does most of the I/O at first
      if (this%cfg%amRoot) then
         
         ! First create a new file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/geometry',form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write geom] Could not open file: ensight/'//trim(this%name)//'/geometry')
         
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'        ; write(iunit) cbuff
         cbuff=trim(adjustl(this%cfg%name))        ; write(iunit) cbuff
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
         
         ! Part header
         cbuff='part'                              ; write(iunit) cbuff
         ibuff=1                                   ; write(iunit) ibuff
         cbuff='Complete geometry'                 ; write(iunit) cbuff  ! We have a single grid-cfg here
         cbuff='block rectilinear iblanked'        ; write(iunit) cbuff
         
         ! Number of cells
         ibuff=this%cfg%nx+1                       ; write(iunit) ibuff
         ibuff=this%cfg%ny+1                       ; write(iunit) ibuff
         ibuff=this%cfg%nz+1                       ; write(iunit) ibuff
         
         ! Mesh
         write(iunit) real(this%cfg%x(this%cfg%imin:this%cfg%imax+1),SP)
         write(iunit) real(this%cfg%y(this%cfg%jmin:this%cfg%jmax+1),SP)
         write(iunit) real(this%cfg%z(this%cfg%kmin:this%cfg%kmax+1),SP)
         
         ! Close the file
         close(iunit)
         
      end if
      
      ! Set array size
      i1=this%cfg%imin_; i2=this%cfg%imax_; if (this%cfg%iproc.eq.this%cfg%npx) i2=i2+1
      j1=this%cfg%jmin_; j2=this%cfg%jmax_; if (this%cfg%jproc.eq.this%cfg%npy) j2=j2+1
      k1=this%cfg%kmin_; k2=this%cfg%kmax_; if (this%cfg%kproc.eq.this%cfg%npz) k2=k2+1
      datasize=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      
      ! Allocate iblank array
      allocate(iblank(i1:i2,j1:j2,k1:k2))
      
      ! Build the blanking info from VF
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               iblank(i,j,k)=nint(maxval(this%cfg%VF(i-1:i,j-1:j,k-1:k)))
            end do
         end do
      end do
      
      ! We need to define a proper MPI-I/O view
      gsizes=[this%cfg%nx+1,this%cfg%ny+1,this%cfg%nz+1]
      lsizes=[i2-i1+1,j2-j1+1,k2-k1+1]
      lstart=[this%cfg%imin_-this%cfg%imin,this%cfg%jmin_-this%cfg%jmin,this%cfg%kmin_-this%cfg%kmin]
      call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_INTEGER,view,ierr)
      call MPI_TYPE_COMMIT(view,ierr)
      
      ! Only need to parallel write masks now
      call MPI_FILE_OPEN(this%cfg%comm,'ensight/'//trim(this%name)//'/geometry',IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[ensight write geom] Problem encountered while parallel writing geometry file')
      call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_INTEGER,view,'native',info_mpiio,ierr)
      call MPI_FILE_WRITE_ALL(ifile,iblank,datasize,MPI_INTEGER,status,ierr)
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Deallocate iblank array
      deallocate(iblank)
      
   end subroutine write_geom
   
   
end module ensight_class
