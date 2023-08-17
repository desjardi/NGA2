!> Basic tracer particle class:
!> Provides support for Lagrangian-transported objects
module tracer_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: tracer
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   
   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing
   
   !> Basic particle object definition
   type :: part
      !> MPI_INTEGER8 data
      integer(kind=8) :: id                !< Particle ID
      !> MPI_DOUBLE_PRECISION data
      real(WP), dimension(3) :: pos        !< Particle center coordinates
      real(WP), dimension(3) :: vel        !< Velocity of particle
      real(WP), dimension(3) :: acc        !< Acceleration of particle
      !> MPI_INTEGER data
      integer , dimension(3) :: ind        !< Index of cell containing particle center
      integer  :: flag                     !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=3
   integer           , dimension(part_nblock) :: part_lblock=[1,9,4]
   type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_INTEGER8,MPI_DOUBLE_PRECISION,MPI_INTEGER]
   !> MPI_PART derived datatype and size
   type(MPI_Datatype) :: MPI_PART
   integer :: MPI_PART_SIZE
   
   !> Fluid tracer solver object definition
   type :: tracer
      
      ! This is our underlying config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_TRACER'  !< Solver name (default=UNNAMED_TRACER)
      
      ! Particle data
      integer :: np                                       !< Global number of particles
      integer :: np_                                      !< Local number of particles
      integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
      type(part), dimension(:), allocatable :: p          !< Array of particles of type part
      
      ! CFL numbers
      real(WP) :: CFLp_x,CFLp_y,CFLp_z                    !< CFL numbers
      
      ! Injection parameters
      real(WP) :: inj_rate                                !< Rate of particle injection (#/time)
      real(WP), dimension(3) :: inj_pos                   !< Center location to inject particles
      
      ! Monitoring info
      real(WP) :: Umin,Umax,Umean,Uvar                    !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean,Vvar                    !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean,Wvar                    !< W velocity info
      integer  :: np_new,np_out                           !< Number of new and removed particles
      
   contains
      procedure :: update_partmesh                        !< Update a partmesh object using current particles
      procedure :: advance                                !< Step forward the particle ODEs
      procedure :: resize                                 !< Resize particle array to given size
      procedure :: recycle                                !< Recycle particle array by removing flagged particles
      procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
      procedure :: read                                   !< Parallel read particles from file
      procedure :: write                                  !< Parallel write particles to file
      procedure :: get_max                                !< Extract various monitoring data
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: inject                                 !< Inject particles at a prescribed boundary
   end type tracer
   
   
   !> Declare tracer solver constructor
   interface tracer
      procedure constructor
   end interface tracer
   
   
contains
   
   
   !> Default constructor for tracer solver
   function constructor(cfg,name) result(self)
      implicit none
      type(tracer) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Allocate variables
      allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
      self%np_=0; self%np=0
      call self%resize(0)
      self%np_new=0; self%np_out=0
      
      ! Initialize MPI derived datatype for a particle
      call prepare_mpi_part()
      
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
      use mpi_f08, only : MPI_SUM,MPI_INTEGER
      use mathtools, only: Pi
      implicit none
      class(tracer), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      type(part) :: pold
      integer :: i,ierr
      
      ! Zero out number of particles removed
      this%np_out=0
      
      ! Advance the equations
      do i=1,this%np_
         ! Avoid particles with id=0
         if (this%p(i)%id.eq.0) cycle
         ! Remember the particle
         pold=this%p(i)
         ! Interpolate fluid velocity to particle location
         this%p(i)%vel=this%cfg%get_velocity(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),U=U,V=V,W=W)
         ! Advance with Euler prediction
         this%p(i)%pos=pold%pos+0.5_WP*dt*this%p(i)%vel
         ! Correct with midpoint rule
         this%p(i)%vel=this%cfg%get_velocity(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),U=U,V=V,W=W)
         this%p(i)%pos=pold%pos+dt*this%p(i)%vel
         ! Update acceleration
         this%p(i)%acc=(this%p(i)%vel-pold%vel)/dt
         ! Relocalize the particle
         this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
         ! Handle particles in walls
         if (this%cfg%VF(this%p(i)%ind(1),this%p(i)%ind(2),this%p(i)%ind(3)).le.0.0_WP) this%p(i)%flag=1
         ! Correct the position to take into account periodicity
         if (this%cfg%xper) this%p(i)%pos(1)=this%cfg%x(this%cfg%imin)+modulo(this%p(i)%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
         if (this%cfg%yper) this%p(i)%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(this%p(i)%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
         if (this%cfg%zper) this%p(i)%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(this%p(i)%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
         ! Relocalize the particle
         this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
         ! Handle particles that have left the domain
         if (this%p(i)%pos(1).lt.this%cfg%x(this%cfg%imin).or.this%p(i)%pos(1).gt.this%cfg%x(this%cfg%imax+1)) this%p(i)%flag=1
         if (this%p(i)%pos(2).lt.this%cfg%y(this%cfg%jmin).or.this%p(i)%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) this%p(i)%flag=1
         if (this%p(i)%pos(3).lt.this%cfg%z(this%cfg%kmin).or.this%p(i)%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) this%p(i)%flag=1
         if (this%p(i)%flag.eq.1) then
            ! Count number of particles removed
            this%np_out=this%np_out+1
         else
            ! Interpolate fluid quantities to particle location
            this%p(i)%vel=this%cfg%get_velocity(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),U=U,V=V,W=W)
         end if
      end do
      
      ! Communicate particles
      call this%sync()
      
      ! Sum up particles removed
      call MPI_ALLREDUCE(this%np_out,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=i
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (this%cfg%amRoot) then
            write(message,'("Tracer solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end subroutine advance
   
   
   !> Inject particles from a prescribed location with given injection rate
   !> Requires injection parameters to be set beforehand
   subroutine inject(this,dt,U,V,W)
      use mpi_f08
      implicit none
      class(tracer), intent(inout) :: this
      real(WP), intent(inout) :: dt                  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: Ngoal,Nadded                       !< Injection rate parameters
      real(WP), save :: previous_error=0.0_WP        !< Store number of particles left over from previous timestep
      integer(kind=8) :: maxid_,maxid                !< Keep track of maximum particle id
      integer :: i,np0_,np_tmp,count,buf,ierr
      
      ! Initial number of particles
      np0_=this%np_
      this%np_new=0
      
      ! Get the particle number density that should be added to the system
      Ngoal =this%inj_rate*dt+previous_error
      Nadded=0.0_WP
      
      ! Determine id to assign to particle
      maxid_=0
      do i=1,this%np_
         maxid_=max(maxid_,this%p(i)%id)
      end do
      call MPI_ALLREDUCE(maxid_,maxid,1,MPI_INTEGER8,MPI_MAX,this%cfg%comm,ierr)
      
      ! Add new particles until desired number is achieved
      do while (Nadded.lt.Ngoal)
         if (this%cfg%amRoot) then
            ! Initialize parameters
            np_tmp = 0
            ! Loop while the added volume is not sufficient
            do while (real(np_tmp,WP).lt.Ngoal-Nadded)
               ! Increment counter
               np_tmp=np_tmp+1
               count = np0_+np_tmp
               ! Create space for new particle
               call this%resize(count)
               ! Set particle ID
               this%p(count)%id=maxid+int(np_tmp,8)
               ! Give a position at the injector to the particle
               this%p(count)%pos=get_position()
               ! Localize the particle
               this%p(count)%ind(1)=this%cfg%imin; this%p(count)%ind(2)=this%cfg%jmin; this%p(count)%ind(3)=this%cfg%kmin
               this%p(count)%ind=this%cfg%get_ijk_global(this%p(count)%pos,this%p(count)%ind)
               ! Interpolate fluid quantities to particle location
               this%p(count)%vel=this%cfg%get_velocity(pos=this%p(count)%pos,i0=this%p(count)%ind(1),j0=this%p(count)%ind(2),k0=this%p(count)%ind(3),U=U,V=V,W=W)
               ! Make it an "official" particle
               this%p(count)%flag=0
            end do
         end if
         ! Communicate particles
         call this%sync()
         ! Loop through newly created particles
         buf=0
         do i=np0_+1,this%np_
            ! Remove if out of bounds
            if (this%cfg%VF(this%p(i)%ind(1),this%p(i)%ind(2),this%p(i)%ind(3)).le.0.0_WP) this%p(i)%flag=1
            if (this%p(i)%flag.eq.0) then
               ! Update the added number for the timestep
               buf=buf+1
               ! Update the max particle id
               maxid=max(maxid,this%p(i)%id)
               ! Increment counter
               this%np_new=this%np_new+1
            end if
         end do
         ! Total particles added
         call MPI_ALLREDUCE(buf,np_tmp,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); Nadded=Nadded+real(np_tmp,WP)
         ! Clean up particles
         call this%recycle()
         ! Update initial npart
         np0_=this%np_
         ! Maximum particle id
         call MPI_ALLREDUCE(maxid,maxid_,1,MPI_INTEGER8,MPI_MAX,this%cfg%comm,ierr); maxid=maxid_
      end do
      
      ! Remember the error
      previous_error=Ngoal-Nadded
      
      ! Sum up injected particles
      call MPI_ALLREDUCE(this%np_new,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_new=i
      
   contains
   
      ! Position for bulk injection of particles
      function get_position() result(pos)
         use random, only: random_uniform
         implicit none
         real(WP), dimension(3) :: pos
         ! Set x position
         pos(1) = this%inj_pos(1)
         ! Random y & z position within a circular region
         pos(2)=random_uniform(lo=this%cfg%y(this%cfg%jmin),hi=this%cfg%y(this%cfg%jmax+1))
         if (this%cfg%nz.eq.1) then
            pos(3) = this%cfg%zm(this%cfg%kmin_)
         else
            pos(3)=random_uniform(lo=this%cfg%z(this%cfg%kmin),hi=this%cfg%z(this%cfg%kmax+1))
         end if
      end function get_position
      
   end subroutine inject
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tracer), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      integer :: i,ierr
      real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z
      
      ! Set the CFLs to zero
      my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP
      do i=1,this%np_
         my_CFLp_x=max(my_CFLp_x,abs(this%p(i)%vel(1))*this%cfg%dxi(this%p(i)%ind(1)))
         my_CFLp_y=max(my_CFLp_y,abs(this%p(i)%vel(2))*this%cfg%dyi(this%p(i)%ind(2)))
         my_CFLp_z=max(my_CFLp_z,abs(this%p(i)%vel(3))*this%cfg%dzi(this%p(i)%ind(3)))
      end do
      my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Return the maximum convective CFL
      cfl=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)
      
   end subroutine get_cfl
   
   
   !> Extract various monitoring data from particle field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tracer), intent(inout) :: this
      real(WP) :: buf,safe_np
      integer :: i,ierr
      
      ! Create safe np
      safe_np=real(max(this%np,1),WP)
      
      ! Diameter and velocity min/max/mean
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); this%Umean=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); this%Vmean=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); this%Wmean=0.0_WP
      do i=1,this%np_
         this%Umin=min(this%Umin,this%p(i)%vel(1)); this%Umax=max(this%Umax,this%p(i)%vel(1)); this%Umean=this%Umean+this%p(i)%vel(1)
         this%Vmin=min(this%Vmin,this%p(i)%vel(2)); this%Vmax=max(this%Vmax,this%p(i)%vel(2)); this%Vmean=this%Vmean+this%p(i)%vel(2)
         this%Wmin=min(this%Wmin,this%p(i)%vel(3)); this%Wmax=max(this%Wmax,this%p(i)%vel(3)); this%Wmean=this%Wmean+this%p(i)%vel(3)
      end do
      call MPI_ALLREDUCE(this%Umin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Umin =buf
      call MPI_ALLREDUCE(this%Umax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Umax =buf
      call MPI_ALLREDUCE(this%Umean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Umean=buf/safe_np
      call MPI_ALLREDUCE(this%Vmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vmin =buf
      call MPI_ALLREDUCE(this%Vmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vmax =buf
      call MPI_ALLREDUCE(this%Vmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vmean=buf/safe_np
      call MPI_ALLREDUCE(this%Wmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wmin =buf
      call MPI_ALLREDUCE(this%Wmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wmax =buf
      call MPI_ALLREDUCE(this%Wmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wmean=buf/safe_np
      
      ! Diameter and velocity variance
      this%Uvar=0.0_WP
      this%Vvar=0.0_WP
      this%Wvar=0.0_WP
      do i=1,this%np_
         this%Uvar=this%Uvar+(this%p(i)%vel(1)-this%Umean)**2.0_WP
         this%Vvar=this%Vvar+(this%p(i)%vel(2)-this%Vmean)**2.0_WP
         this%Wvar=this%Wvar+(this%p(i)%vel(3)-this%Wmean)**2.0_WP
      end do
      call MPI_ALLREDUCE(this%Uvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Uvar=buf/safe_np
      call MPI_ALLREDUCE(this%Vvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vvar=buf/safe_np
      call MPI_ALLREDUCE(this%Wvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wvar=buf/safe_np
      
   end subroutine get_max
   
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(tracer), intent(inout) :: this
      class(partmesh), intent(inout) :: pmesh
      integer :: i
      ! Reset particle mesh storage
      call pmesh%reset()
      ! Nothing else to do if no particle is present
         if (this%np_.eq.0) return
         ! Copy particle info
         call pmesh%set_size(this%np_)
         do i=1,this%np_
            pmesh%pos(:,i)=this%p(i)%pos
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
      if (ierr.ne.0) call die('[tracer prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Synchronize particle arrays across processors
   subroutine sync(this)
      use mpi_f08
      implicit none
      class(tracer), intent(inout) :: this
      integer, dimension(0:this%cfg%nproc-1) :: nsend_proc,nrecv_proc
      integer, dimension(0:this%cfg%nproc-1) :: nsend_disp,nrecv_disp
      integer :: n,prank,ierr
      type(part), dimension(:), allocatable :: buf_send
      ! Recycle first to minimize communication load
      call this%recycle()
      ! Prepare information about what to send
      nsend_proc=0
      do n=1,this%np_
         prank=this%cfg%get_rank(this%p(n)%ind)
         nsend_proc(prank)=nsend_proc(prank)+1
      end do
      nsend_proc(this%cfg%rank)=0
      ! Inform processors of what they will receive
      call MPI_ALLtoALL(nsend_proc,1,MPI_INTEGER,nrecv_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
      ! Prepare displacements for all-to-all
      nsend_disp(0)=0
      nrecv_disp(0)=this%np_   !< Directly add particles at the end of main array
      do n=1,this%cfg%nproc-1
         nsend_disp(n)=nsend_disp(n-1)+nsend_proc(n-1)
         nrecv_disp(n)=nrecv_disp(n-1)+nrecv_proc(n-1)
      end do
      ! Allocate buffer to send particles
      allocate(buf_send(sum(nsend_proc)))
      ! Pack the particles in the send buffer
      nsend_proc=0
      do n=1,this%np_
         ! Get the rank
         prank=this%cfg%get_rank(this%p(n)%ind)
         ! Skip particles still inside
         if (prank.eq.this%cfg%rank) cycle
         ! Pack up for sending
         nsend_proc(prank)=nsend_proc(prank)+1
         buf_send(nsend_disp(prank)+nsend_proc(prank))=this%p(n)
         ! Flag particle for removal
         this%p(n)%flag=1
      end do
      ! Allocate buffer for receiving particles
      call this%resize(this%np_+sum(nrecv_proc))
      ! Perform communication
      call MPI_ALLtoALLv(buf_send,nsend_proc,nsend_disp,MPI_PART,this%p,nrecv_proc,nrecv_disp,MPI_PART,this%cfg%comm,ierr)
      ! Deallocate buffer
      deallocate(buf_send)
      ! Recycle to remove duplicate particles
      call this%recycle()
   end subroutine sync
   
   
   !> Adaptation of particle array size
   subroutine resize(this,n)
      implicit none
      class(tracer), intent(inout) :: this
      integer, intent(in) :: n
      type(part), dimension(:), allocatable :: tmp
      integer :: size_now,size_new
      ! Resize particle array to size n
      if (.not.allocated(this%p)) then
         ! Allocate directly to size n
         allocate(this%p(n))
         this%p(1:n)%flag=1
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
      implicit none
      class(tracer), intent(inout) :: this
      integer :: new_size,i,ierr
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
      this%np_=new_size
      call MPI_ALLGATHER(this%np_,1,MPI_INTEGER,this%np_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
      this%np=sum(this%np_proc)
   end subroutine recycle
   
   
   !> Parallel write particles to file
   subroutine write(this,filename)
      use mpi_f08
      use messager, only: die
      use parallel, only: info_mpiio
      implicit none
      class(tracer), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: i,ierr,iunit

      ! Root serial-writes the file header
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[tracer write] Problem encountered while serial-opening data file: '//trim(filename))
         ! Number of particles and particle object size
         write(iunit) this%np,MPI_PART_SIZE
         ! Done with the header
         close(iunit)
      end if
      
      ! The rest is done in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[tracer write] Problem encountered while parallel-opening data file: '//trim(filename))
      
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
      class(tracer), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
      integer :: i,j,ierr,npadd,psize,nchunk,cnt
      integer, dimension(:,:), allocatable :: ppp
      
      ! First open the file in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[tracer read] Problem encountered while reading data file: '//trim(filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
      
      ! Remember current position
      call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)
      
      ! Check compatibility of particle type
      if (psize.ne.MPI_PART_SIZE) call die('[tracer read] Particle type unreadable')
      
      ! Naively share reading task among all processors
      nchunk=int(npadd/(this%cfg%nproc*part_chunk_size))+1
      allocate(ppp(this%cfg%nproc,nchunk))
      ppp=int(npadd/(this%cfg%nproc*nchunk))
      cnt=0
      out:do j=1,nchunk
         do i=1,this%cfg%nproc
            cnt=cnt+1
            if (cnt.gt.mod(npadd,this%cfg%nproc*nchunk)) exit out
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
   
   
end module tracer_class
