!> Basic Lagrangian particle solver class:
!> Provides support for Lagrangian-transported objects
module lpt_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
   use partmesh_class, only: partmesh
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: lpt
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   
   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing
   
   
   !> Basic particle object definition
   type :: part
      integer(kind=8) :: id                !< Particle ID
      real(WP) :: d                        !< Particle diameter
      real(WP), dimension(3) :: pos        !< Particle center coordinates
      real(WP), dimension(3) :: vel        !< Velocity of particle
      real(WP) :: dt                       !< Time step size for the particle
      integer , dimension(3) :: ind        !< Index of cell containing particle center
      integer  :: flag                     !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=3
   integer           , dimension(part_nblock) :: part_lblock=[1,8,4]
   type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_INTEGER8,MPI_DOUBLE_PRECISION,MPI_INTEGER]
   !> MPI_PART derived datatype and size
   type(MPI_Datatype) :: MPI_PART
   integer :: MPI_PART_SIZE
   
   !> Lagrangian particle tracking solver object definition
   type :: lpt
      
      ! This is our underlying config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_LPT'     !< Solver name (default=UNNAMED_LPT)
      
      ! Particle data
      integer :: np                                       !< Global number of particles
      integer :: np_                                      !< Local number of particles
      integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
      type(part), dimension(:), allocatable :: p          !< Array of particles of type part
      
      ! Particle density
      real(WP) :: rho                                     !< Density of particle
      
      ! Solver parameters
      real(WP) :: nstep=1                                 !< Number of substeps (default=1)
      logical  :: twoway=.false.                          !< Two-way coupling   (default=no)
      
      ! Basic representation of the particle mesh
      type(partmesh) :: pmesh
      
      ! Monitoring info
      real(WP) :: dmin,dmax,dmean,dvar                    !< Diameter info
      real(WP) :: Umin,Umax,Umean,Uvar                    !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean,Vvar                    !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean,Wvar                    !< W velocity info
      integer  :: np_new,np_out                           !< Number of new and removed particles
      
   contains
      procedure :: update_partmesh                        !< Create a simple particle mesh for visualization
      procedure :: advance                                !< Step forward the particle ODEs
      procedure :: resize                                 !< Resize particle array to given size
      procedure :: recycle                                !< Recycle particle array by removing flagged particles
      procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
      procedure :: read                                   !< Parallel read particles from file
      procedure :: write                                  !< Parallel write particles to file
      procedure :: get_max                                !< Extract various monitoring data
   end type lpt
   
   
   !> Declare lpt solver constructor
   interface lpt
      procedure constructor
   end interface lpt
   
contains
   
   
   !> Default constructor for lpt solver
   function constructor(cfg,name) result(self)
      implicit none
      type(lpt) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Allocate variables
      allocate(self%np_proc(1:self%cfg%nproc))
      
      ! Initialize MPI derived datatype for a particle
      call prepare_mpi_part()
      
      ! Initialize a particle mesh for output purposes
      self%pmesh=partmesh(name='lpt')
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%cfg%amRoot) then
            write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function constructor
   
   
   !> Advance the particle equations by a specified time step dt
   subroutine advance(this,dt,U,V,W)
      implicit none
      class(lpt), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      
      
      
   end subroutine advance
   
   
   !> Extract various monitoring data from particle field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(lpt), intent(inout) :: this
      real(WP) :: buf
      integer :: i,ierr
      
      ! Diameter and velocity min/max/mean
      this%dmin=huge(1.0_WP); this%dmax=-huge(1.0_WP); this%dmean=0.0_WP
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); this%Umean=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); this%Vmean=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); this%Wmean=0.0_WP
      do i=1,this%np_
         this%dmin=min(this%dmin,this%p(i)%d     ); this%dmax=max(this%dmax,this%p(i)%d     ); this%dmean=this%dmean+this%p(i)%d
         this%Umin=min(this%Umin,this%p(i)%vel(1)); this%Umax=max(this%Umax,this%p(i)%vel(1)); this%Umean=this%Umean+this%p(i)%vel(1)
         this%Vmin=min(this%Vmin,this%p(i)%vel(2)); this%Vmax=max(this%Vmax,this%p(i)%vel(2)); this%Vmean=this%Vmean+this%p(i)%vel(2)
         this%Wmin=min(this%Wmin,this%p(i)%vel(3)); this%Wmax=max(this%Wmax,this%p(i)%vel(3)); this%Wmean=this%Wmean+this%p(i)%vel(3)
      end do
      call MPI_ALLREDUCE(this%dmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%dmin =buf
      call MPI_ALLREDUCE(this%dmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%dmax =buf
      call MPI_ALLREDUCE(this%dmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dmean=buf/this%np
      call MPI_ALLREDUCE(this%Umin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Umin =buf
      call MPI_ALLREDUCE(this%Umax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Umax =buf
      call MPI_ALLREDUCE(this%Umean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Umean=buf/this%np
      call MPI_ALLREDUCE(this%Vmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vmin =buf
      call MPI_ALLREDUCE(this%Vmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vmax =buf
      call MPI_ALLREDUCE(this%Vmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vmean=buf/this%np
      call MPI_ALLREDUCE(this%Wmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wmin =buf
      call MPI_ALLREDUCE(this%Wmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wmax =buf
      call MPI_ALLREDUCE(this%Wmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wmean=buf/this%np
      
      ! Diameter and velocity variance
      this%dvar=0.0_WP
      this%Uvar=0.0_WP
      this%Vvar=0.0_WP
      this%Wvar=0.0_WP
      do i=1,this%np_
         this%dvar=this%dvar+(this%p(i)%d     -this%dmean)**2.0_WP
         this%Uvar=this%Uvar+(this%p(i)%vel(1)-this%Umean)**2.0_WP
         this%Vvar=this%Vvar+(this%p(i)%vel(2)-this%Vmean)**2.0_WP
         this%Wvar=this%Wvar+(this%p(i)%vel(3)-this%Wmean)**2.0_WP
      end do
      call MPI_ALLREDUCE(this%dvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dvar=buf/this%np
      call MPI_ALLREDUCE(this%Uvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Uvar=buf/this%np
      call MPI_ALLREDUCE(this%Vvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vvar=buf/this%np
      call MPI_ALLREDUCE(this%Wvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wvar=buf/this%np
      
   end subroutine get_max
   
   
   !> Update particle mesh for visualization
   subroutine update_partmesh(this)
      implicit none
      class(lpt), intent(inout) :: this
      integer :: i
      ! Reset particle mesh storage
      call this%pmesh%reset()
      ! Handle zero particle case by adding fictitious particle
      if (this%np_.eq.0) then
         call this%pmesh%set_size(1)
         this%pmesh%pos(:,1)=0.0_WP
         this%pmesh%d(1)    =0.0_WP
         return
      end if
      ! Copy particle info
      call this%pmesh%set_size(this%np_)
      do i=1,this%np_
         this%pmesh%pos(:,i)=this%p(i)%pos
         this%pmesh%d(i)    =this%p(i)%d
      end do
   end subroutine update_partmesh
   
   
   !> Creation of the MPI datatype for particle
   subroutine prepare_mpi_part()
      use mpi_f08
      use messager, only: die
      implicit none
      integer(MPI_ADDRESS_KIND), dimension(part_nblock) :: disp
      integer(MPI_ADDRESS_KIND) :: lb,extent
      type(MPI_Datatype) :: MPI_PART_TMP
      integer :: i,mysize,ierr
      ! Prepare the displacement array
      disp(1)=0
      do i=2,part_nblock
         call MPI_Type_size(part_tblock(i-1),mysize,ierr)
         disp(i)=disp(i-1)+int(mysize,MPI_ADDRESS_KIND)*int(part_lblock(i-1),MPI_ADDRESS_KIND)
      end do
      ! Create and commit the new type
      call MPI_Type_create_struct(part_nblock,part_lblock,disp,part_tblock,MPI_PART_TMP,ierr)
      call MPI_Type_get_extent(MPI_PART_TMP,lb,extent,ierr)
      call MPI_Type_create_resized(MPI_PART_TMP,lb,extent,MPI_PART,ierr)
      call MPI_Type_commit(MPI_PART,ierr)
      ! If a problem was encountered, say it
      if (ierr.ne.0) call die('[lpt prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Synchronize particle arrays across processors
   subroutine sync(this)
      use mpi_f08
      implicit none
      class(lpt), intent(inout) :: this
      integer, dimension(this%cfg%nproc) :: who_send,who_recv,counter
      integer :: i,nb_recv,nb_send,np_old
      integer :: rank_send,rank_recv,prank,ierr
      type(MPI_status) :: status
      type(part), dimension(:,:), allocatable :: buf_send
      ! Recycle first to minimize communication load
      call this%recycle()
      ! Prepare information about who sends what to whom
      who_send=0
      do i=1,this%np_
         prank=this%cfg%get_rank(this%p(i)%ind)+1
         who_send(prank)=who_send(prank)+1
      end do
      ! Remove the diagonal since we won't self-send
      who_send(this%cfg%rank+1)=0
      ! Prepare information about who receives what from whom
      do rank_recv=1,this%cfg%nproc
         call MPI_gather(who_send(rank_recv),1,MPI_INTEGER,who_recv,1,MPI_INTEGER,rank_recv-1,this%cfg%comm,ierr)
      end do
      ! Prepare the buffers
      nb_send=maxval(who_send)
      nb_recv=sum(who_recv)
      ! Allocate buffers to send particles
      allocate(buf_send(this%cfg%nproc,nb_send))
      ! Find and pack the particles to be sent
      counter=0
      do i=1,this%np_
         ! Get the proc
         prank=this%cfg%get_rank(this%p(i)%ind)+1
         if (prank.ne.this%cfg%rank+1) then
            ! Prepare for sending
            counter(prank)=counter(prank)+1
            buf_send(prank,counter(prank))=this%p(i)
            ! Need to remove the particle
            this%p(i)%flag=1
         end if
      end do
      ! Everybody resizes
      np_old=this%np_
      call this%resize(this%np_+nb_recv)
      ! Loop through the procs, pack particles in buf_send, send, unpack
      do rank_send=1,this%cfg%nproc
         if (this%cfg%rank+1.eq.rank_send) then
            ! I'm the sender of particles
            do rank_recv=1,this%cfg%nproc
               if (who_send(rank_recv).gt.0) then
                  call MPI_send(buf_send(rank_recv,:),who_send(rank_recv),MPI_PART,rank_recv-1,0,this%cfg%comm,ierr)
               end if
            end do
         else
            ! I'm not the sender, I receive
            if (who_recv(rank_send).gt.0) then
               call MPI_recv(this%p(np_old+1:np_old+who_recv(rank_send)),who_recv(rank_send),MPI_PART,rank_send-1,0,this%cfg%comm,status,ierr)
               np_old=np_old+who_recv(rank_send)
            end if
         end if
      end do
      ! Done, deallocate
      deallocate(buf_send)
      ! Recycle to remove duplicate particles
      call this%recycle()
   end subroutine sync
   
   
   !> Adaptation of particle array size
   subroutine resize(this,n)
      implicit none
      class(lpt), intent(inout) :: this
      integer, intent(in) :: n
      type(part), dimension(:), allocatable :: tmp
      integer :: size_now,size_new
      ! Resize part array to size n
      if (.not.allocated(this%p)) then
         ! part is of size 0
         if (n.gt.0) then
            ! Allocate directly to size n
            allocate(this%p(n))
            this%p(1:n)%flag=1
         end if
      else if (n.eq.0) then
         ! part is associated, but we want to empty it
         deallocate(this%p)
      else
         ! Update from a non-zero size to another non-zero size
         size_now=size(this%p,dim=1)
         if (n.gt.size_now) then
            size_new=max(n,int(real(size_now,WP)*coeff_up))
            allocate(tmp(size_new))
            tmp(1:size_now)=this%p
            tmp(size_now+1:)%flag=1
            call move_alloc(tmp,this%p)
         else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
            allocate(tmp(n))
            tmp(1:n)=this%p(1:n)
            call move_alloc(tmp,this%p)
         end if
      end if
   end subroutine resize
   
   
   !> Clean-up of particle array by removing flag=1 particles
   subroutine recycle(this)
      use mpi_f08
      implicit none
      class(lpt), intent(inout) :: this
      integer :: new_size,ierr,i
      ! Compact all active particles at the beginning of the array
      new_size=0
      if (allocated(this%p)) then
         do i=1,size(this%p,dim=1)
            if (this%p(i)%flag.ne.1) then
               new_size=new_size+1
               if (i.ne.new_size) then
                  this%p(new_size)=this%p(i)
                  this%p(i)%flag=1
               end if
            end if
         end do
      end if
      ! Resize to new size
      call this%resize(new_size)
      ! Update number of particles
      this%np_=size(this%p,dim=1)
      call MPI_ALLGATHER(this%np_,1,MPI_INTEGER,this%np_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
      this%np=sum(this%np_proc)
   end subroutine recycle
   
   
   !> Parallel write particles to file
   subroutine write(this,filename)
      use mpi_f08
      use messager, only: die
      use parallel, only: info_mpiio
      implicit none
      class(lpt), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: i,ierr,iunit
      
      ! Root serial-writes the file header
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[lpt write] Problem encountered while serial-opening data file: '//trim(filename))
         ! Number of particles and particle object size
         write(iunit) this%np,MPI_PART_SIZE
         ! Done with the header
         close(iunit)
      end if
      
      ! The rest is done in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[lpt write] Problem encountered while parallel-opening data file: '//trim(filename))
      
      ! Get current position
      call MPI_FILE_GET_POSITION(ifile,offset,ierr)
      
      ! Compute the offset and write
      do i=1,this%cfg%rank
         offset=offset+int(this%np_proc(i),MPI_OFFSET_KIND)*int(MPI_PART_SIZE,MPI_OFFSET_KIND)
      end do
      if (this%np_.gt.0) call MPI_FILE_WRITE_AT(ifile,offset,this%p,this%np_,MPI_PART,status,ierr)
      
      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (this%cfg%amRoot) then
            write(message,'("Wrote ",i0," particles to file [",a,"] on partitioned grid [",a,"]")') this%np,trim(filename),trim(this%cfg%name)
            if (verbose.gt.2) write(output_unit,'(a)') trim(message)
            if (verbose.gt.1) call log(message)
         end if
      end block logging
      
   end subroutine write
   
   
   !> Parallel read particles to file
   subroutine read(this,filename)
      use mpi_f08
      use messager, only: die
      use parallel, only: info_mpiio
      implicit none
      class(lpt), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
      integer :: i,j,ierr,npadd,psize,nchunk,count
      integer, dimension(:,:), allocatable :: ppp
      
      ! First open the file in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[lpt read] Problem encountered while reading data file: '//trim(filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
      
      ! Remember current position
      call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)
      
      ! Check compatibility of particle type
      if (psize.ne.MPI_PART_SIZE) call die('[lpt read] Particle type unreadable')
      
      ! Naively share reading task among all processors
      nchunk=int(npadd/(this%cfg%nproc*part_chunk_size))+1
      allocate(ppp(this%cfg%nproc,nchunk))
      ppp=int(npadd/(this%cfg%nproc*nchunk))
      count=0
      out:do j=1,nchunk
         do i=1,this%cfg%nproc
            count=count+1
            if (count.gt.mod(npadd,this%cfg%nproc*nchunk)) exit out
            ppp(i,j)=ppp(i,j)+1
         end do
      end do out
      
      ! Read by chunk
      do j=1,nchunk
         ! Find offset
         offset=header_offset+int(MPI_PART_SIZE,MPI_OFFSET_KIND)*int(sum(ppp(1:this%cfg%rank,:))+sum(ppp(this%cfg%rank+1,1:j-1)),MPI_OFFSET_KIND)
         ! Resize particle array
         call this%resize(this%np_+ppp(this%cfg%rank+1,j))
         ! Read this file
         call MPI_FILE_READ_AT(ifile,offset,this%p(this%np_+1:this%np_+ppp(this%cfg%rank+1,j)),ppp(this%cfg%rank+1,j),MPI_PART,status,ierr)
         ! Most general case: relocate every droplet
         do i=this%np_+1,this%np_+ppp(this%cfg%rank+1,j)
            this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
         end do
         ! Exchange all that
         call this%sync()
      end do
      
      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (this%cfg%amRoot) then
            write(message,'("Read ",i0," particles from file [",a,"] on partitioned grid [",a,"]")') npadd,trim(filename),trim(this%cfg%name)
            if (verbose.gt.2) write(output_unit,'(a)') trim(message)
            if (verbose.gt.1) call log(message)
         end if
      end block logging
      
   end subroutine read
   
   
end module lpt_class
