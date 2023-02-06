!> Basic Lagrangian particle solver class:
!> Provides support for Lagrangian-transported objects
module lpt_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
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
      !> MPI_INTEGER8 data
      integer(kind=8) :: id                !< Particle ID
      !> MPI_DOUBLE_PRECISION data
      real(WP) :: d                        !< Particle diameter
      real(WP), dimension(3) :: pos        !< Particle center coordinates
      real(WP), dimension(3) :: vel        !< Velocity of particle
      real(WP), dimension(3) :: col        !< Collision force
      real(WP) :: dt                       !< Time step size for the particle
      !> MPI_INTEGER data
      integer , dimension(3) :: ind        !< Index of cell containing particle center
      integer  :: flag                     !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=3
   integer           , dimension(part_nblock) :: part_lblock=[1,11,4]
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
      
      ! Overlap particle (i.e., ghost) data
      integer :: ng_                                      !< Local number of ghosts
      type(part), dimension(:), allocatable :: g          !< Array of ghosts of type part
      
      ! Particle density
      real(WP) :: rho                                     !< Density of particle
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP            !< Acceleration of gravity
      
      ! Solver parameters
      real(WP) :: nstep=1                                 !< Number of substeps (default=1)
      
      ! Collisional parameters
      real(WP) :: Tcol                                    !< Characteristic collision time scale
      real(WP) :: e_n=0.9_WP                              !< Normal restitution coefficient
      real(WP) :: clip_col=0.2_WP                         !< Maximum allowable overlap
      real(WP), dimension(:,:,:),   allocatable :: Wdist  !< Signed wall distance - naive for now (could be redone with FMM)
      real(WP), dimension(:,:,:,:), allocatable :: Wnorm  !< Wall normal function - naive for now (could be redone with FMM)
      
      ! Monitoring info
      real(WP) :: VFmin,VFmax,VFmean,VFvar                !< Volume fraction info
      real(WP) :: dmin,dmax,dmean,dvar                    !< Diameter info
      real(WP) :: Umin,Umax,Umean,Uvar                    !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean,Vvar                    !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean,Wvar                    !< W velocity info
      integer  :: np_new,np_out                           !< Number of new and removed particles
      
      ! Particle volume fraction
      real(WP), dimension(:,:,:), allocatable :: VF       !< Particle volume fraction, cell-centered
      
      ! Momentum source
      real(WP), dimension(:,:,:), allocatable :: srcU     !< U momentum source on mesh, cell-centered
      real(WP), dimension(:,:,:), allocatable :: srcV     !< V momentum source on mesh, cell-centered
      real(WP), dimension(:,:,:), allocatable :: srcW     !< W momentum source on mesh, cell-centered
      
      ! Filtering operation
      real(WP) :: filter_width                            !< Characteristic filter width
      real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z    !< Divergence operator
      real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z    !< Gradient operator
      
   contains
      procedure :: update_partmesh                        !< Update a partmesh object using current particles
      procedure :: collide                                !< Evaluate interparticle collision force
      procedure :: advance                                !< Step forward the particle ODEs
      procedure :: get_rhs                                !< Compute rhs of particle odes
      procedure :: resize                                 !< Resize particle array to given size
      procedure :: resize_ghost                           !< Resize ghost array to given size
      procedure :: recycle                                !< Recycle particle array by removing flagged particles
      procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
      procedure :: share                                  !< Share particles across interprocessor boundaries
      procedure :: read                                   !< Parallel read particles from file
      procedure :: write                                  !< Parallel write particles to file
      procedure :: get_max                                !< Extract various monitoring data
      procedure :: update_VF                              !< Compute particle volume fraction
      procedure :: filter                                 !< Apply volume filtering to field
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
      
      ! Set default filter width to zero by default
      self%filter_width=0.0_WP
      
      ! Allocate finite volume divergence operators
      allocate(self%div_x(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      allocate(self%div_y(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      allocate(self%div_z(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               self%div_x(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< Divergence from [x ,ym,zm]
               self%div_y(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,y ,zm]
               self%div_z(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,ym,z ]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(self%grd_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
      allocate(self%grd_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
      allocate(self%grd_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               self%grd_x(:,i,j,k)=self%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< Gradient in x from [xm,ym,zm] to [x,ym,zm]
               self%grd_y(:,i,j,k)=self%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< Gradient in y from [xm,ym,zm] to [xm,y,zm]
               self%grd_z(:,i,j,k)=self%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< Gradient in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Loop over the domain and zero divergence in walls
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                  self%div_x(:,i,j,k)=0.0_WP
                  self%div_y(:,i,j,k)=0.0_WP
                  self%div_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Zero out gradient to wall faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i-1,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j-1,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j,k-1).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (self%cfg%nx.eq.1) then
         self%div_x=0.0_WP
         self%grd_x=0.0_WP
      end if
      if (self%cfg%ny.eq.1) then
         self%div_y=0.0_WP
         self%grd_y=0.0_WP
      end if
      if (self%cfg%nz.eq.1) then
         self%div_z=0.0_WP
         self%grd_z=0.0_WP
      end if
      
      ! Generate a wall distance/norm function
      allocate(self%Wdist(  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
      allocate(self%Wnorm(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
      ! First pass to set correct sign
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                  self%Wdist(i,j,k)=-sqrt(self%cfg%xL**2+self%cfg%yL**2+self%cfg%zL**2)
               else
                  self%Wdist(i,j,k)=+sqrt(self%cfg%xL**2+self%cfg%yL**2+self%cfg%zL**2)
               end if
               self%Wnorm(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      ! Second pass to compute local distance
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_+1,self%cfg%imaxo_
               if (self%Wdist(i,j,k)*self%Wdist(i-1,j,k).lt.0.0_WP) then
                  ! There is a wall at x(i)
                  if (abs(self%cfg%xm(i  )-self%cfg%x(i)).lt.abs(self%Wdist(i  ,j,k))) then
                     self%Wdist(i  ,j,k)=sign(self%cfg%xm(i  )-self%cfg%x(i),self%Wdist(i  ,j,k))
                     self%Wnorm(:,i  ,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-1,j,k),0.0_WP,0.0_WP]
                  end if
                  if (abs(self%cfg%xm(i-1)-self%cfg%x(i)).lt.abs(self%Wdist(i-1,j,k))) then
                     self%Wdist(i-1,j,k)=sign(self%cfg%xm(i-1)-self%cfg%x(i),self%Wdist(i-1,j,k))
                     self%Wnorm(:,i-1,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-1,j,k),0.0_WP,0.0_WP]
                  end if
               end if
            end do
         end do
      end do
      call self%cfg%sync(self%Wdist)
      call self%cfg%sync(self%Wnorm)
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_+1,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%Wdist(i,j,k)*self%Wdist(i,j-1,k).lt.0.0_WP) then
                  ! There is a wall at y(j)
                  if (abs(self%cfg%ym(j  )-self%cfg%y(j)).lt.abs(self%Wdist(i,j  ,k))) then
                     self%Wdist(i,j  ,k)=sign(self%cfg%ym(j  )-self%cfg%y(j),self%Wdist(i,j  ,k))
                     self%Wnorm(:,i,j  ,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-1,k),0.0_WP]
                  end if
                  if (abs(self%cfg%ym(j-1)-self%cfg%y(j)).lt.abs(self%Wdist(i,j-1,k))) then
                     self%Wdist(i,j-1,k)=sign(self%cfg%ym(j-1)-self%cfg%y(j),self%Wdist(i,j-1,k))
                     self%Wnorm(:,i,j-1,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-1,k),0.0_WP]
                  end if
               end if
            end do
         end do
      end do
      call self%cfg%sync(self%Wdist)
      call self%cfg%sync(self%Wnorm)
      do k=self%cfg%kmino_+1,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%Wdist(i,j,k)*self%Wdist(i,j,k-1).lt.0.0_WP) then
                  ! There is a wall at z(k)
                  if (abs(self%cfg%zm(k  )-self%cfg%z(k)).lt.abs(self%Wdist(i,j,k  ))) then
                     self%Wdist(i,j,k  )=sign(self%cfg%zm(k  )-self%cfg%z(k),self%Wdist(i,j,k  ))
                     self%Wnorm(:,i,j,k  )=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-1)]
                  end if
                  if (abs(self%cfg%zm(k-1)-self%cfg%z(k)).lt.abs(self%Wdist(i,j,k-1))) then
                     self%Wdist(i,j,k-1)=sign(self%cfg%zm(k-1)-self%cfg%z(k),self%Wdist(i,j,k-1))
                     self%Wnorm(:,i,j,k-1)=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-1)]
                  end if
               end if
            end do
         end do
      end do
      call self%cfg%sync(self%Wdist)
      call self%cfg%sync(self%Wnorm)
      
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
   
   
   !> Resolve collisional interaction between particles
   !> Requires Tcol and e_n to be set beforehand
   subroutine collide(this,dt)
      implicit none
      class(lpt), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
      
      ! Start by zeroing out the collision force
      zero_force: block
         integer :: i
         do i=1,this%np_
            this%p(i)%col=0.0_WP
         end do
      end block zero_force
      
      ! Then share particles across overlap
      call this%share()
      
      ! We can now assemble particle-in-cell information
      pic_prep: block
         use mpi_f08
         integer :: i,ip,jp,kp,ierr
         integer :: mymax_npic,max_npic
         
         ! Allocate number of particle in cell
         allocate(npic(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); npic=0
         
         ! Count particles and ghosts per cell
         do i=1,this%np_
            ip=this%p(i)%ind(1); jp=this%p(i)%ind(2); kp=this%p(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
         end do
         
         ! Get maximum number of particle in cell
         mymax_npic=maxval(npic); call MPI_ALLREDUCE(mymax_npic,max_npic,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         
         ! Allocate pic map
         allocate(ipic(1:max_npic,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); ipic=0
         
         ! Assemble pic map
         npic=0
         do i=1,this%np_
            ip=this%p(i)%ind(1); jp=this%p(i)%ind(2); kp=this%p(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=i
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=-i
         end do
         
      end block pic_prep
      
      ! Finally, calculate collision force
      collision_force: block
         use mathtools, only: Pi,normalize
         integer :: i1,i2,ii,jj,kk,nn
         real(WP) :: d1,m1,d2,m2,d12,m12
         real(WP), dimension(3) :: r1,v1,r2,v2,v12,n12,f_n
         real(WP) :: k_n,eta_n,lne,pilne2,rnv,r_influ,delta_n
         
         ! Precompute ln(e_n) and Pi^2+ln(e_n)^2
         lne=log(this%e_n); pilne2=Pi**2+lne**2
         
         ! Loop over all local particles
         do i1=1,this%np_
            
            ! Store particle data
            r1=this%p(i1)%pos
            v1=this%p(i1)%vel
            d1=this%p(i1)%d
            m1=this%rho*Pi/6.0_WP*d1**3
            
            ! First collide with walls
            d12=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=this%Wdist,bc='d')
            n12=this%Wnorm(:,this%p(i1)%ind(1),this%p(i1)%ind(2),this%p(i1)%ind(3))
            n12=-normalize(n12+[epsilon(1.0_WP),epsilon(1.0_WP),epsilon(1.0_WP)])
            rnv=dot_product(v1,n12)
            r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
            delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)
            
            ! Assess if there is collision
            if (delta_n.gt.0.0_WP) then
               ! Normal collision
               k_n=m1/(this%Tcol**2*pilne2)
               eta_n=-2.0_WP*lne*sqrt(m1*k_n)/pilne2
               f_n=-k_n*delta_n*n12-eta_n*rnv*n12
               ! Calculate collision force
               this%p(i1)%col=this%p(i1)%col+f_n/m1
            end if
            
            ! Loop over nearest cells
            do kk=this%p(i1)%ind(3)-1,this%p(i1)%ind(3)+1
               do jj=this%p(i1)%ind(2)-1,this%p(i1)%ind(2)+1
                  do ii=this%p(i1)%ind(1)-1,this%p(i1)%ind(1)+1
                     
                     ! Loop over particles in that cell
                     do nn=1,npic(ii,jj,kk)
                        
                        ! Get index of neighbor particle
                        i2=ipic(nn,ii,jj,kk)
                        
                        ! Get relevant data from correct storage
                        if (i2.gt.0) then
                           r2=this%p(i2)%pos
                           v2=this%p(i2)%vel
                           d2=this%p(i2)%d
                           m2=this%rho*Pi/6.0_WP*d2**3
                        else if (i2.lt.0) then
                           i2=-i2
                           r2=this%g(i2)%pos
                           v2=this%g(i2)%vel
                           d2=this%g(i2)%d
                           m2=this%rho*Pi/6.0_WP*d2**3
                        end if
                        
                        ! Compute relative information
                        d12=norm2(r1-r2)
                        if (d12.lt.10.0_WP*epsilon(d12)) cycle !< this should skip auto-collision
                        n12=(r2-r1)/d12
                        v12=v1-v2
                        rnv=dot_product(v12,n12)
                        r_influ=min(abs(rnv)*dt,0.1_WP*(d1+d2))
                        delta_n=min(0.5_WP*(d1+d2)+r_influ-d12,this%clip_col*0.5_WP*(d1+d2))
                        
                        ! Assess if there is collision
                        if (delta_n.gt.0.0_WP) then
                           ! Normal collision
                           m12=m1*m2/(m1+m2)
                           k_n=m12/(this%Tcol**2*pilne2)
                           eta_n=-2.0_WP*lne*sqrt(m12*k_n)/pilne2
                           f_n=-k_n*delta_n*n12-eta_n*rnv*n12
                           ! Calculate collision force
                           this%p(i1)%col=this%p(i1)%col+f_n/m1
                        end if
                        
                     end do
                     
                  end do
               end do
            end do
            
            ! Deal with dimensionality
            if (this%cfg%nx.eq.1) this%p(i1)%col(1)=0.0_WP
            if (this%cfg%ny.eq.1) this%p(i1)%col(2)=0.0_WP
            if (this%cfg%nz.eq.1) this%p(i1)%col(3)=0.0_WP
            
         end do
         
      end block collision_force
      
   end subroutine collide

   
   !> Advance the particle equations by a specified time step dt
   subroutine advance(this,dt,U,V,W,rho,visc)
      use mathtools, only: Pi
      implicit none
      class(lpt), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i
      real(WP) :: mydt,dt_done
      real(WP), dimension(3) :: acc,dmom
      type(part) :: pold
      
      ! Zero out source term arrays
      this%srcU=0.0_WP
      this%srcV=0.0_WP
      this%srcW=0.0_WP
      
      ! Advance the equations
      do i=1,this%np_
         ! Time-integrate until dt_done=dt
         dt_done=0.0_WP
         do while (dt_done.lt.dt)
            ! Decide the timestep size
            mydt=min(this%p(i)%dt,dt-dt_done)
            ! Remember the particle
            pold=this%p(i)
            ! Advance with Euler prediction
            call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=this%p(i),acc=acc,opt_dt=this%p(i)%dt)
            this%p(i)%pos=pold%pos+0.5_WP*mydt*this%p(i)%vel
            this%p(i)%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity+this%p(i)%col)
            ! Correct with midpoint rule
            call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=this%p(i),acc=acc,opt_dt=this%p(i)%dt)
            this%p(i)%pos=pold%pos+mydt*this%p(i)%vel
            this%p(i)%vel=pold%vel+mydt*(acc+this%gravity+this%p(i)%col)
            ! Relocalize
            this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
            ! Send source term back to the mesh
            dmom=mydt*acc*this%rho*Pi/6.0_WP*this%p(i)%d**3
            if (this%cfg%nx.gt.1) call this%cfg%set_scalar(Sp=-dmom(1),pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%srcU,bc='n')
            if (this%cfg%ny.gt.1) call this%cfg%set_scalar(Sp=-dmom(2),pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%srcV,bc='n')
            if (this%cfg%nz.gt.1) call this%cfg%set_scalar(Sp=-dmom(3),pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%srcW,bc='n')
            ! Increment
            dt_done=dt_done+mydt
         end do
         ! Correct the position to take into account periodicity
         if (this%cfg%xper) this%p(i)%pos(1)=this%cfg%x(this%cfg%imin)+modulo(this%p(i)%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
         if (this%cfg%yper) this%p(i)%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(this%p(i)%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
         if (this%cfg%zper) this%p(i)%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(this%p(i)%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
         ! Handle particles that have left the domain
         if (this%p(i)%pos(1).lt.this%cfg%x(this%cfg%imin).or.this%p(i)%pos(1).gt.this%cfg%x(this%cfg%imax+1)) this%p(i)%flag=1
         if (this%p(i)%pos(2).lt.this%cfg%y(this%cfg%jmin).or.this%p(i)%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) this%p(i)%flag=1
         if (this%p(i)%pos(3).lt.this%cfg%z(this%cfg%kmin).or.this%p(i)%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) this%p(i)%flag=1
         ! Relocalize the particle
         this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
      end do
      
      ! Communicate particles
      call this%sync()
      
      ! Divide source arrays by volume, sum at boundaries, and volume filter
      this%srcU=this%srcU/this%cfg%vol; call this%cfg%syncsum(this%srcU); call this%filter(this%srcU)
      this%srcV=this%srcV/this%cfg%vol; call this%cfg%syncsum(this%srcV); call this%filter(this%srcV)
      this%srcW=this%srcW/this%cfg%vol; call this%cfg%syncsum(this%srcW); call this%filter(this%srcW)
      
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
            write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end subroutine advance
   
   
   !> Calculate RHS of the particle ODEs
   subroutine get_rhs(this,U,V,W,rho,visc,p,acc,opt_dt)
      implicit none
      class(lpt), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      type(part), intent(in) :: p
      real(WP), dimension(3), intent(out) :: acc
      real(WP), intent(out) :: opt_dt
      real(WP) :: Re,corr,tau,b1,b2
      real(WP) :: fvisc,frho,pVF,fVF
      real(WP), dimension(3) :: fvel
      ! Interpolate the fluid phase velocity to the particle location
      fvel=this%cfg%get_velocity(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),U=U,V=V,W=W)
      ! Interpolate the fluid phase viscosity to the particle location
      fvisc=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=visc,bc='n')
      fvisc=fvisc+epsilon(1.0_WP)
      ! Interpolate the fluid phase density to the particle location
      frho=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=rho,bc='n')
      ! Interpolate the particle volume fraction to the particle location
      pVF=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=this%VF,bc='n')
      fVF=1.0_WP-pVF
      ! Particle Reynolds number
      Re=frho*norm2(p%vel-fvel)*p%d/fvisc+epsilon(1.0_WP)
      ! Stokes correlation
      !corr=1.0_WP
      ! Schiller Naumann correlation
      !corr=1.0_WP+0.15_WP*Re**(0.687_WP)
      ! Tenneti and Subramaniam (2011)
      b1=5.81_WP*pVF/fVF**3+0.48_WP*pVF**(1.0_WP/3.0_WP)/fVF**4
      b2=pVF**3*Re*(0.95_WP+0.61_WP*pVF**3/fVF**2)
      corr=fVF*(1.0_WP+0.15_WP*Re**(0.687_WP))/fVF**3+b1+b2
      ! Particle response time
      tau=this%rho*p%d**2/(18.0_WP*fvisc*corr)
      ! Return acceleration and optimal timestep size
      acc=fVF*(fvel-p%vel)/tau
      ! Add buoyancy
      acc=acc-frho/this%rho*this%gravity
      opt_dt=tau/real(this%nstep,WP)
   end subroutine get_rhs


   !> Update particle volume fraction using our current particles
   subroutine update_VF(this)
      use mathtools, only: Pi
      implicit none
      class(lpt), intent(inout) :: this
      integer :: i
      real(WP) :: Vp
      ! Reset volume fraction
      this%VF=0.0_WP
      ! Transfer particle volume
      do i=1,this%np_
         ! Skip inactive particle
         if (this%p(i)%flag.eq.1) cycle
         ! Transfer particle volume
         Vp=Pi/6.0_WP*this%p(i)%d**3
         call this%cfg%set_scalar(Sp=Vp,pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%VF,bc='n')
      end do
      this%VF=this%VF/this%cfg%vol
      ! Sum at boundaries
      call this%cfg%syncsum(this%VF)
      ! Apply volume filter
      call this%filter(this%VF)
   end subroutine update_VF


   !> Laplacian filtering operation
   subroutine filter(this,A)
      implicit none
      class(lpt), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: filter_coeff
      integer :: i,j,k,n,nstep
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      
      ! Return without filtering if filter width is zero
      if (this%filter_width.le.0.0_WP) return
      
      ! Recompute filter coeff and number of explicit steps needed
      filter_coeff=0.5_WP*(this%filter_width/(2.0_WP*sqrt(2.0_WP*log(2.0_WP))))**2
      nstep=ceiling(6.0_WP*filter_coeff/this%cfg%min_meshsize**2)
      filter_coeff=filter_coeff/real(nstep,WP)
      
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Apply filter
      do n=1,nstep
         ! Diffusive flux of A
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
                  FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
                  FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
               end do
            end do
         end do
         ! Divergence of fluxes
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  A(i,j,k)=A(i,j,k)+sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
               end do
            end do
         end do
         ! Sync A
         call this%cfg%sync(A)
      end do
      
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      
   end subroutine filter
   
   
   !> Extract various monitoring data from particle field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(lpt), intent(inout) :: this
      real(WP) :: buf,safe_np
      integer :: i,j,k,ierr
      
      ! Create safe np
      safe_np=real(max(this%np,1),WP)
      
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
      call MPI_ALLREDUCE(this%dmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dmean=buf/safe_np
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
      call MPI_ALLREDUCE(this%dvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dvar=buf/safe_np
      call MPI_ALLREDUCE(this%Uvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Uvar=buf/safe_np
      call MPI_ALLREDUCE(this%Vvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vvar=buf/safe_np
      call MPI_ALLREDUCE(this%Wvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wvar=buf/safe_np
      
      ! Get mean, max, and min volume fraction
      this%VFmean=0.0_WP
      this%VFmax =-huge(1.0_WP)
      this%VFmin =+huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%VFmean=this%VFmean+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*this%VF(i,j,k)
               this%VFmax =max(this%VFmax,this%VF(i,j,k))
               this%VFmin =min(this%VFmin,this%VF(i,j,k))
            end do
         end do
      end do
      call MPI_ALLREDUCE(this%VFmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFmean=buf/this%cfg%fluid_vol
      call MPI_ALLREDUCE(this%VFmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%VFmax =buf
      call MPI_ALLREDUCE(this%VFmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%VFmin =buf

      ! Get volume fraction variance
      this%VFvar=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%VFvar=this%VFvar+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*(this%VF(i,j,k)-this%VFmean)**2.0_WP
            end do
         end do
      end do
      call MPI_ALLREDUCE(this%VFvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFvar=buf/this%cfg%fluid_vol
      
   end subroutine get_max
   
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(lpt), intent(inout) :: this
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
      if (ierr.ne.0) call die('[lpt prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Synchronize particle arrays across processors
   subroutine sync(this)
      use mpi_f08
      implicit none
      class(lpt), intent(inout) :: this
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
   

   !> Share particles across processor boundaries
   subroutine share(this,nover)
      use mpi_f08
      use messager, only: warn,die
      implicit none
      class(lpt), intent(inout) :: this
      integer, optional :: nover
      type(part), dimension(:), allocatable :: tosend
      type(part), dimension(:), allocatable :: torecv
      integer :: no,nsend,nrecv
      type(MPI_Status) :: status
      integer :: icnt,isrc,idst,ierr
      integer :: i,n
      
      ! Check overlap size
      if (present(nover)) then
         no=nover
         if (no.gt.this%cfg%no) then
            call warn('[lpt_class share] Specified overlap is larger than that of cfg - reducing no')
            no=this%cfg%no
         else if (no.le.0) then
            call die('[lpt_class share] Specified overlap cannot be less or equal to zero')
         end if
      else
         no=1
      end if
      
      ! Clean up ghost array
      call this%resize_ghost(n=0); this%ng_=0
      
      ! Share ghost particles to the left in x
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).lt.this%cfg%imin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).lt.this%cfg%imin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%xper.and.tosend(nsend)%ind(1).lt.this%cfg%imin+no) then
               tosend(nsend)%pos(1)=tosend(nsend)%pos(1)+this%cfg%xL
               tosend(nsend)%ind(1)=tosend(nsend)%ind(1)+this%cfg%nx
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,0,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles to the right in x
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).gt.this%cfg%imax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).gt.this%cfg%imax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%xper.and.tosend(nsend)%ind(1).gt.this%cfg%imax-no) then
               tosend(nsend)%pos(1)=tosend(nsend)%pos(1)-this%cfg%xL
               tosend(nsend)%ind(1)=tosend(nsend)%ind(1)-this%cfg%nx
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,0,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles to the left in y
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).lt.this%cfg%jmin+no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)+this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)+this%cfg%ny
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,1,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles to the right in y
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).gt.this%cfg%jmax-no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)-this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)-this%cfg%ny
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,1,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles to the left in z
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%zper.and.tosend(nsend)%ind(3).lt.this%cfg%kmin+no) then
               tosend(nsend)%pos(3)=tosend(nsend)%pos(3)+this%cfg%zL
               tosend(nsend)%ind(3)=tosend(nsend)%ind(3)+this%cfg%nz
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,2,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles to the right in z
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%zper.and.tosend(nsend)%ind(3).gt.this%cfg%kmax-no) then
               tosend(nsend)%pos(3)=tosend(nsend)%pos(3)-this%cfg%zL
               tosend(nsend)%ind(3)=tosend(nsend)%ind(3)-this%cfg%nz
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,2,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)

   end subroutine share
   
   
   !> Adaptation of particle array size
   subroutine resize(this,n)
      implicit none
      class(lpt), intent(inout) :: this
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


   !> Adaptation of ghost array size
   subroutine resize_ghost(this,n)
      implicit none
      class(lpt), intent(inout) :: this
      integer, intent(in) :: n
      type(part), dimension(:), allocatable :: tmp
      integer :: size_now,size_new
      ! Resize ghost array to size n
      if (.not.allocated(this%g)) then
         ! Allocate directly to size n
         allocate(this%g(n))
         this%g(1:n)%flag=1
      else
         ! Update from a non-zero size to another non-zero size
         size_now=size(this%g,dim=1)
         if (n.gt.size_now) then
            size_new=max(n,int(real(size_now,WP)*coeff_up))
            allocate(tmp(size_new))
            tmp(1:size_now)=this%g
            tmp(size_now+1:)%flag=1
            call move_alloc(tmp,this%g)
         else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
            allocate(tmp(n))
            tmp(1:n)=this%g(1:n)
            call move_alloc(tmp,this%g)
         end if
      end if
   end subroutine resize_ghost
   
   
   !> Clean-up of particle array by removing flag=1 particles
   subroutine recycle(this)
      implicit none
      class(lpt), intent(inout) :: this
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
      integer :: i,j,ierr,npadd,psize,nchunk,cnt
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
   
   
end module lpt_class
