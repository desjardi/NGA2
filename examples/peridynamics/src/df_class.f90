!> Basic direct forcing IBM class:
!> Provides support for Lagrangian marker particles
module df_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: dfibm
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   
   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing
   
   !> Basic marker particle definition
   type :: part
      !> MPI_DOUBLE_PRECISION data
      real(WP) :: dV                       !< Element volume
      real(WP), dimension(3) :: pos        !< Particle center coordinates
      real(WP), dimension(3) :: vel        !< Velocity of particle
      !> MPI_INTEGER data
      integer :: id                        !< ID the object is associated with
      integer , dimension(3) :: ind        !< Index of cell containing particle center
      integer :: flag                      !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=2
   integer           , dimension(part_nblock) :: part_lblock=[7,5]
   type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_DOUBLE_PRECISION,MPI_INTEGER]
   !> MPI_PART derived datatype and size
   type(MPI_Datatype) :: MPI_PART
   integer :: MPI_PART_SIZE
   
   !> Basic ibm object definition
   type :: obj
      !> MPI_DOUBLE_PRECISION data
      real(WP) :: vol                      !< Object volume
      real(WP), dimension(3) :: pos        !< Center of mass
      real(WP), dimension(3) :: vel        !< Translational velocity of the object
      real(WP), dimension(3) :: angVel     !< Angular velocity of the object
      real(WP), dimension(3) :: F          !< Hydrodynamic force
   end type obj
   
   !> Direct forcing IBM solver object definition
   type :: dfibm
   
      ! This is our underlying config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_DFIBM'   !< Solver name (default=UNNAMED_DFIBM
      
      ! Marker particle data
      integer :: np                                       !< Global number of particles
      integer :: np_                                      !< Local number of particles
      integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
      type(part), dimension(:), allocatable :: p          !< Array of particles of type part
      
      ! Object data
      integer :: nobj                                     !< Global number of objects
      type(obj), dimension(:), allocatable :: o           !< Array of objects of type obj
      
      ! CFL numbers
      real(WP) :: CFLp_x,CFLp_y,CFLp_z                    !< CFL numbers
      
      ! Solver parameters
      real(WP) :: nstep=1                                 !< Number of substeps (default=1)
      logical :: can_move                                 !< Flag to allow moving IBM objects
      
      ! Monitoring info
      real(WP) :: VFmin,VFmax,VFmean                 !< Volume fraction info
      real(WP) :: Umin,Umax,Umean                    !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean                    !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean                    !< W velocity info
      real(WP) :: Fx,Fy,Fz                           !< Total force
      
      ! Volume fraction associated with IBM projection
      real(WP), dimension(:,:,:), allocatable :: VF       !< Volume fraction, cell-centered
      
      ! Momentum source
      real(WP), dimension(:,:,:), allocatable :: srcU     !< U momentum source on mesh, cell-centered
      real(WP), dimension(:,:,:), allocatable :: srcV     !< V momentum source on mesh, cell-centered
      real(WP), dimension(:,:,:), allocatable :: srcW     !< W momentum source on mesh, cell-centered
      
   contains
      procedure :: update_partmesh                        !< Update a partmesh object using current particles
      procedure :: setup_obj                              !< Setup IBM object
      procedure :: get_source                             !< Compute direct forcing source
      procedure :: resize                                 !< Resize particle array to given size
      procedure :: recycle                                !< Recycle particle array by removing flagged particles
      procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
      procedure :: read                                   !< Parallel read particles from file
      procedure :: write                                  !< Parallel write particles to file
      procedure :: get_max                                !< Extract various monitoring data
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: update_VF                              !< Compute volume fraction
      procedure :: get_delta                              !< Compute regularized delta function
      procedure :: interpolate                            !< Interpolation routine from mesh=>marker
      procedure :: extrapolate                            !< Extrapolation routine from marker=>mesh
   end type dfibm
   
   
   !> Declare df solver constructor
   interface dfibm
      procedure constructor
   end interface dfibm
   
contains


   !> Default constructor for direct forcing solver
   function constructor(cfg,name) result(self)
      implicit none
      type(dfibm) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Allocate variables
      allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
      self%np_=0; self%np=0
      call self%resize(0)
      
      ! Initialize MPI derived datatype for a particle
      call prepare_mpi_part()
      
      ! Allocate VF and src arrays on cfg mesh
      allocate(self%VF  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF  =0.0_WP
      allocate(self%srcU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcU=0.0_WP
      allocate(self%srcV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcV=0.0_WP
      allocate(self%srcW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcW=0.0_WP
      
      ! Initialize object
      self%can_move=.false.
      self%nobj=0
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%cfg%amRoot) then
            write(message,'("IBM solver [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function constructor
   
   
   !> Setup IBM objects, each processors own all objects
   subroutine setup_obj(this)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(dfibm), intent(inout) :: this
      integer :: i,j,n,ibuf,ierr
      real(WP) :: myVol,dV
      real(WP), dimension(3) :: pos0,dist
      ! Determine number of objects based on marker ID
      n=1
      do i=1,this%np_
         n=max(n,this%p(i)%id)
      end do
      call MPI_ALLREDUCE(n,this%nobj,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      ! Allocate and zero out the object array
      allocate(this%o(1:this%nobj))
      do i=1,this%nobj
         ! Zero-out object properties
         this%o(i)%pos=0.0_WP
         this%o(i)%vel=0.0_WP
         this%o(i)%angVel=0.0_WP
         this%o(i)%F=0.0_WP
         ! Compute center of mass and volume
         myVol=0.0_WP
         do j=1,this%np_
            ! Determine if particle associated with object
            if (this%p(j)%id.eq.i) then
               myVol=myVol+this%p(j)%dV
               this%o(i)%pos=this%o(i)%pos+this%p(j)%pos*this%p(j)%dV
            end if
         end do
         call MPI_ALLREDUCE(myVol,dV,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%o(i)%vol=dV
         call MPI_ALLREDUCE(this%o(i)%pos,pos0,3,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%o(i)%pos=pos0/dV
      end do
   end subroutine setup_obj
   
   
   !> Compute direct forcing source by a specified time step dt
   subroutine get_source(this,dt,U,V,W,rho)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(dfibm), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: dti,rho_
      real(WP), dimension(3) :: vel,src
      
      ! Zero out source term arrays
      this%srcU=0.0_WP
      this%srcV=0.0_WP
      this%srcW=0.0_WP
      
      ! Zero-out forces on objects
      do i=1,this%nobj
         this%o(i)%F=0.0_WP
      end do
      
      ! Advance the equations
      dti=1.0_WP/dt
      do i=1,this%np_
         ! Interpolate the velocity to the particle location
         vel=0.0_WP
         if (this%cfg%nx.gt.1) vel(1)=this%interpolate(A=U,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                                 ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='U')
         if (this%cfg%ny.gt.1) vel(2)=this%interpolate(A=V,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                                 ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='V')
         if (this%cfg%nz.gt.1) vel(3)=this%interpolate(A=W,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                                 ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='W')
         rho_=this%interpolate(A=rho,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                           ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='SC')
         ! Compute the source term
         src=(this%p(i)%vel-vel)*this%p(i)%dV
         ! Send source term back to the mesh
         if (this%cfg%nx.gt.1) call this%extrapolate(Ap=src(1),xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                                     ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%srcU,dir='U')
         if (this%cfg%ny.gt.1) call this%extrapolate(Ap=src(2),xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                                     ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%srcV,dir='V')
         if (this%cfg%nz.gt.1) call this%extrapolate(Ap=src(3),xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                                     ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%srcW,dir='W')
         ! Sum up force on object (Newton's 3rd law)
         j=max(this%p(i)%id,1)
         this%o(j)%F=this%o(j)%F-rho_*src*this%p(i)%dV*dti
      end do
      
      ! Sum at boundaries
      call this%cfg%syncsum(this%srcU)
      call this%cfg%syncsum(this%srcV)
      call this%cfg%syncsum(this%srcW)
      
      ! Sum over each object
      do j=1,this%nobj
         call MPI_ALLREDUCE(this%o(j)%F,src,3,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%o(j)%F=src
      end do
      
      ! Recompute volume fraction
      call this%update_VF()
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (this%cfg%amRoot) then
            write(message,'("IBM solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end subroutine get_source
   
   
   !> Update particle volume fraction using our current particles
   subroutine update_VF(this)
      use mathtools, only: Pi
      implicit none
      class(dfibm), intent(inout) :: this
      integer :: i
      ! Reset volume fraction
      this%VF=0.0_WP
      ! Transfer particle volume
      do i=1,this%np_
         ! Skip inactive particle
         if (this%p(i)%flag.eq.1) cycle
         ! Transfer volume to mesh
         call this%extrapolate(Ap=this%p(i)%dV,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),&
         &                                     ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%VF,dir='SC')
      end do
      ! Sum at boundaries
      call this%cfg%syncsum(this%VF)
   end subroutine update_VF
   
   
   !> Compute regularized delta function
   subroutine get_delta(this,delta,ic,jc,kc,xp,yp,zp,dir)
      implicit none
      class(dfibm), intent(inout) :: this
      real(WP), intent(out) :: delta   !< Return delta function
      integer, intent(in) :: ic,jc,kc  !< Cell index
      real(WP), intent(in) :: xp,yp,zp !< Position of marker 
      character(len=*) :: dir
      real(WP) :: deltax,deltay,deltaz,r
      
      ! Compute in X
      if (trim(adjustl(dir)).eq.'U') then
         r=(xp-this%cfg%x(ic))*this%cfg%dxmi(ic)
         deltax=roma_kernel(r)*this%cfg%dxmi(ic)
      else
         r=(xp-this%cfg%xm(ic))*this%cfg%dxi(ic)
         deltax=roma_kernel(r)*this%cfg%dxi(ic)
      end if
      
      ! Compute in Y
      if (trim(adjustl(dir)).eq.'V') then
         r=(yp-this%cfg%y(jc))*this%cfg%dymi(jc)
         deltay=roma_kernel(r)*this%cfg%dymi(jc)
      else
         r=(yp-this%cfg%ym(jc))*this%cfg%dyi(jc)
         deltay=roma_kernel(r)*this%cfg%dyi(jc)
      end if
      
      ! Compute in Z
      if (trim(adjustl(dir)).eq.'W') then
         r=(zp-this%cfg%z(kc))*this%cfg%dzmi(kc)
         deltaz=roma_kernel(r)*this%cfg%dzmi(kc)
      else
         r=(zp-this%cfg%zm(kc))*this%cfg%dzi(kc)
         deltaz=roma_kernel(r)*this%cfg%dzi(kc)
      end if
      
      ! Put it all together
      delta=deltax*deltay*deltaz
      
   contains
      ! Mollification kernel
      ! Roma A, Peskin C and Berger M 1999 J. Comput. Phys. 153 509–534
      function roma_kernel(r) result(phi)
         implicit none
         real(WP), intent(in) :: r
         real(WP)             :: phi
         if (abs(r).le.0.5_WP) then
            phi=1.0_WP/3.0_WP*(1.0_WP+sqrt(-3.0_WP*r**2+1.0_WP))
         else if (abs(r).gt.0.5_WP .and. abs(r).le.1.5_WP) then
            phi=1.0_WP/6.0_WP*(5.0_WP-3.0_WP*abs(r)-sqrt(-3.0_WP*(1.0_WP-abs(r))**2+1.0_WP))
         else
            phi=0.0_WP
         end if
      end function roma_kernel
      
   end subroutine get_delta
   
   
   !> Interpolation routine
   function interpolate(this,A,xp,yp,zp,ip,jp,kp,dir) result(Ap)
      implicit none
      class(dfibm), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: A         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(in) :: xp,yp,zp
      integer, intent(in) :: ip,jp,kp
      character(len=*) :: dir
      real(WP) :: Ap
      integer :: di,dj,dk
      integer :: i1,i2,j1,j2,k1,k2
      real(WP), dimension(-2:+2,-2:+2,-2:+2) :: delta
      ! Get the interpolation points
      i1=ip-2; i2=ip+2
      j1=jp-2; j2=jp+2
      k1=kp-2; k2=kp+2
      ! Loop over neighboring cells and compute regularized delta function
      do dk=-2,+2
         do dj=-2,+2
            do di=-2,+2
               call this%get_delta(delta=delta(di,dj,dk),ic=ip+di,jc=jp+dj,kc=kp+dk,xp=xp,yp=yp,zp=zp,dir=trim(dir))
            end do
         end do
      end do
      ! Perform the actual interpolation on Ap
      Ap = sum(delta*A(i1:i2,j1:j2,k1:k2))*this%cfg%vol(ip,jp,kp)
   end function interpolate
   
   
   !> Extrapolation routine
   subroutine extrapolate(this,Ap,xp,yp,zp,ip,jp,kp,A,dir)
      use messager, only: die
      implicit none
      class(dfibm), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:) :: A         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(in) :: xp,yp,zp
      integer, intent(in) :: ip,jp,kp
      real(WP), intent(in) :: Ap
      character(len=*) :: dir
      real(WP), dimension(-2:+2,-2:+2,-2:+2) :: delta
      integer  :: di,dj,dk
      ! If particle has left processor domain or reached last ghost cell, kill job
      if ( ip.lt.this%cfg%imin_-1.or.ip.gt.this%cfg%imax_+1.or.&
      &    jp.lt.this%cfg%jmin_-1.or.jp.gt.this%cfg%jmax_+1.or.&
      &    kp.lt.this%cfg%kmin_-1.or.kp.gt.this%cfg%kmax_+1) then
         write(*,*) ip,jp,kp,xp,yp,zp
         call die('[df extrapolate] Particle has left the domain')
      end if
      ! Loop over neighboring cells and compute regularized delta function
      do dk=-2,+2
         do dj=-2,+2
            do di=-2,+2
               call this%get_delta(delta=delta(di,dj,dk),ic=ip+di,jc=jp+dj,kc=kp+dk,xp=xp,yp=yp,zp=zp,dir=trim(dir))
            end do
         end do
      end do
      ! Perform the actual extrapolation on A
      A(ip-2:ip+2,jp-2:jp+2,kp-2:kp+2)=A(ip-2:ip+2,jp-2:jp+2,kp-2:kp+2)+delta*Ap
   end subroutine extrapolate
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(dfibm), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      integer :: i,ierr
      real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z
      
      ! Return if not used
      if (.not.this%can_move) then
         cfl=0.0_WP
         return
      end if
      
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
      
      ! Return the maximum CFL
      cfl=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)
      
   end subroutine get_cfl
   
   
   !> Extract various monitoring data from particle field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(dfibm), intent(inout) :: this
      real(WP) :: buf,safe_np
      integer :: i,j,k,ierr
      
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
      
      ! Get mean, max, and min volume fraction and total force
      this%VFmean=0.0_WP
      this%VFmax =-huge(1.0_WP)
      this%VFmin =+huge(1.0_WP)
      this%Fx=0.0_WP; this%Fy=0.0_WP; this%Fz=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%VFmean=this%VFmean+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*this%VF(i,j,k)
               this%VFmax =max(this%VFmax,this%VF(i,j,k))
               this%VFmin =min(this%VFmin,this%VF(i,j,k))
               this%Fx=this%Fx+sum(this%o(:)%F(1))*this%cfg%vol(i,j,k)
               this%Fy=this%Fy+sum(this%o(:)%F(2))*this%cfg%vol(i,j,k)
               this%Fz=this%Fz+sum(this%o(:)%F(3))*this%cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(this%VFmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFmean=buf/this%cfg%vol_total
      call MPI_ALLREDUCE(this%VFmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%VFmax =buf
      call MPI_ALLREDUCE(this%VFmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%VFmin =buf
      call MPI_ALLREDUCE(this%Fx    ,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Fx=buf/this%cfg%vol_total
      call MPI_ALLREDUCE(this%Fy    ,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Fy=buf/this%cfg%vol_total
      call MPI_ALLREDUCE(this%Fz    ,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Fz=buf/this%cfg%vol_total
      
   end subroutine get_max
   
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(dfibm), intent(inout) :: this
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
      if (ierr.ne.0) call die('[dfibm prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Synchronize particle arrays across processors
   subroutine sync(this)
      use mpi_f08
      implicit none
      class(dfibm), intent(inout) :: this
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
      class(dfibm), intent(inout) :: this
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
      class(dfibm), intent(inout) :: this
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
      class(dfibm), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: i,ierr,iunit
      
      ! Root serial-writes the file header
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[dfibm write] Problem encountered while serial-opening data file: '//trim(filename))
         ! Number of particles and particle object size
         write(iunit) this%np,MPI_PART_SIZE
         ! Done with the header
         close(iunit)
      end if
      
      ! The rest is done in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[dfibm write] Problem encountered while parallel-opening data file: '//trim(filename))
      
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
      class(dfibm), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
      integer :: i,j,ierr,npadd,psize,nchunk,cnt
      integer, dimension(:,:), allocatable :: ppp
      
      ! First open the file in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[dfibm read] Problem encountered while reading data file: '//trim(filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
      
      ! Remember current position
      call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)
      
      ! Check compatibility of particle type
      if (psize.ne.MPI_PART_SIZE) call die('[dfibm read] Particle type unreadable')
      
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
   
   
end module df_class
