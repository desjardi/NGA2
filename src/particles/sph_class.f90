!> Smooth particle hydrodynamics solver
module sph_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER,MPI_DOUBLE_PRECISION
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: sph


   ! List of known available smoothing kernels
   integer, parameter, public :: cubic=1
   integer, parameter, public :: quintic=2
   integer, parameter, public :: gaussian=3
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   

   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing
   

   !> Particle definition
   type :: part
      !> MPI_DOUBLE_PRECISION data
      real(WP) :: dp                         !< Effective diameter
      real(WP) :: vol                        !< Particle volume
      real(WP) :: rho                        !< Particle density
      real(WP) :: P                          !< Particle pressure
      real(WP) :: dRhoDt                     !< Right-hand side for particle density
      real(WP), dimension(3) :: pos          !< Particle center coordinates
      real(WP), dimension(3) :: vel          !< Velocity of particle
      real(WP), dimension(3) :: dxdt         !< Right-hand side for particle position
      real(WP), dimension(3) :: dudt         !< Right-hand side for particle velocity
      !> MPI_INTEGER data
      integer :: id                          !< ID the object is associated with
      integer , dimension(3) :: ind          !< Index of cell containing particle center
      integer :: flag                        !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=2
   integer           , dimension(part_nblock) :: part_lblock=[17,5]
   type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_DOUBLE_PRECISION,MPI_INTEGER]
   !> MPI_PART derived datatype and size
   type(MPI_Datatype) :: MPI_PART
   integer :: MPI_PART_SIZE
   
   
   !> Smooth particle hydrodynamics object definition
   type :: sph
   
      ! This config is used for parallelization and for neighbor detection
      class(config), pointer :: cfg
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_SPH'
      
      ! Fluid properties
      real(WP) :: rho                                     !< Reference density
      real(WP) :: visc                                    !< Reference viscosity
      real(WP) :: c                                       !< Sound speed
      real(WP) :: cmin=0.0_WP                             !< Minimum sound speed
      real(WP) :: gamm=7.0_WP                             !< Specific heat ratio

      ! Kernel properties
      real(WP) :: h                                       !< Smoothing length
      real(WP) :: hi                                      !< Inverse of smoothing length
      real(WP) :: rmax                                    !< Cutoff radius
      
      ! Global and local particle data
      integer :: np                                       !< Global number of particles
      integer :: np_                                      !< Local number of particles
      integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
      type(part), dimension(:), allocatable :: p          !< Array of particles of type part
      
      ! Overlap particle (i.e., ghost) data
      integer :: ng_                                      !< Local number of ghosts
      type(part), dimension(:), allocatable :: g          !< Array of ghosts of type part
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP

      ! Wall detection
      real(WP), dimension(:,:,:),   allocatable :: Wdist  !< Signed wall distance - naive for now (could be redone with FMM)
      real(WP), dimension(:,:,:,:), allocatable :: Wnorm  !< Wall normal function - naive for now (could be redone with FMM)
      
      ! CFL numbers
      real(WP) :: CFLp_x,CFLp_y,CFLp_z,CFLp_c,CFLp_f
      
      ! Number of substeps for time integrator
      real(WP) :: nstep=1
      
      ! Monitoring info
      real(WP) :: RHOmin,RHOmax,RHOmean              !< RHO info
      real(WP) :: Umin,Umax,Umean                    !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean                    !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean                    !< W velocity info
      real(WP) :: dRHO                               !< Relative change in density
      integer  :: np_out                             !< Number of particles leaving the domain
      
    contains
      procedure :: sph_init                          !< Initialize the smoothing kernel and pressure
      procedure :: sph_update                        !< Update SPH dependent variables
      procedure :: get_rhs                           !< Compute rhs of particle odes
      procedure :: advance                           !< Step forward the particle ODEs
      procedure :: get_cfl                           !< Calculate maximum CFL
      procedure :: get_max                           !< Extract various monitoring data
      procedure :: update_partmesh                   !< Update a partmesh object using current particles
      procedure :: share                             !< Share particles across interprocessor boundaries
      procedure :: sync                              !< Synchronize particles across interprocessor boundaries
      procedure :: resize                            !< Resize particle array to given size
      procedure :: resize_ghost                      !< Resize ghost array to given size
      procedure :: recycle                           !< Recycle particle array by removing flagged particles
      procedure :: write                             !< Parallel write particles to file
      procedure :: read                              !< Parallel read particles from file
   end type sph
   
   
   !> Declare sph constructor
   interface sph
      procedure constructor
   end interface sph
   
contains
   
   
  ! Smoothing kernel (not normalized)
  function W(d,h,kernel,rmax)
    implicit none
    real(WP), intent(in) :: d,h
    integer , intent(in) :: kernel
    real(WP), intent(in), optional :: rmax
    real(WP) :: W
    real(WP) :: q
    q=d/h
    W = 0.0_WP
    select case (kernel)
    case (cubic)
       if (q.ge.0.0_WP .and. q.lt.1.0_WP) then
          W = 1.0_WP-1.5_WP*q**2+0.75_WP*q**3
       elseif (q.ge.1.0_WP .and. q.lt.2.0_WP) then
          W = 0.25_WP*(2.0_WP-q)**3
       end if
    case (quintic)
       if (q.ge.0.0_WP .and. q.lt.1.0_WP) then
          W = (3.0_WP-q)**5-6.0_WP*(2.0_WP-q)**5+15.0_WP*(1.0_WP-q)**5
       elseif (q.ge.1.0_WP .and. q.lt.2.0_WP) then
          W = (3.0_WP-q)**5-6.0_WP*(2.0_WP-q)**5
       elseif (q.ge.2.0_WP .and. q.lt.3.0_WP) then
          W = (3.0_WP-q)**5
       end if
    case (gaussian)
          W = max((exp(-q**2) - exp(-rmax)), 0.0_WP)
    end select
  end function W


  ! Gradient of the smoothing kernel (not normalized)
  function gradW(d,h,x1,x2,kernel)
    implicit none
    real(WP), dimension(3), intent(in) :: x1,x2
    real(WP), intent(in) :: d,h
    integer, intent(in) :: kernel
    real(WP), dimension(3) :: gradW
    real(WP) :: q,dWdr
    q = d/h
    select case (kernel)
    case(cubic)
       dWdr = 0.0_WP
       if (q.ge.0.0_WP .and. q.lt.1.0_WP) then
          dWdr = -0.75_WP/h*(4.0_WP*q-3.0_WP*q**2)
       elseif (q.ge.1.0_WP .and. q.lt.2.0_WP) then
          dWdr = -0.75_WP/h*(2.0_WP-q)**2
       end if
    case(quintic)
       dWdr = 0.0_WP
       if (q.ge.0.0_WP .and. q.lt.1.0_WP) then
          dWdr = -(5.0_WP*(3.0_WP-q)**4-30.0_WP*(2.0_WP-q)**4+75.0_WP*(1.0_WP-q)**4)/h
       elseif (q.ge.1.0_WP .and. q.lt.2.0_WP) then
          dWdr = -(5.0_WP*(3.0_WP-q)**4-30.0_WP*(2.0_WP-q)**4)/h
       elseif (q.ge.2.0_WP .and. q.lt.3.0_WP) then
          dWdr = -5.0_WP*(3.0_WP-q)**4/h
       end if
    case (gaussian)
       dWdr = -2.0_WP*d/h**2*exp(-q**2)
    end select
    gradW = dWdr*(x1-x2)/d
  end function gradW


   !> Default constructor for SPH solver
   function constructor(cfg,name) result(self)
      implicit none
      type(sph) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k,l
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Set default smoothing length based on underlying mesh
      self%h=self%cfg%min_meshsize/2.0_WP
      self%hi=1.0_WP/self%h
      self%rmax=2.0_WP*self%h
      
      ! Allocate variables
      allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
      self%np_=0; self%np=0
      call self%resize(0)
      
      ! Initialize MPI derived datatype for a particle
      call prepare_mpi_part()

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
    ! Second pass to compute local distance (get 2 closest cells)
    do l=1,2
       do k=self%cfg%kmino_,self%cfg%kmaxo_
          do j=self%cfg%jmino_,self%cfg%jmaxo_
             do i=self%cfg%imino_+l,self%cfg%imaxo_
                if (self%Wdist(i,j,k)*self%Wdist(i-l,j,k).lt.0.0_WP) then
                   ! There is a wall at x(i)
                   if (abs(self%cfg%xm(i  )-self%cfg%x(i-l+1)).lt.abs(self%Wdist(i  ,j,k))) then
                      self%Wdist(i  ,j,k)=sign(self%cfg%xm(i  )-self%cfg%x(i-l+1),self%Wdist(i  ,j,k))
                      self%Wnorm(:,i  ,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-l,j,k),0.0_WP,0.0_WP]
                   end if
                   if (abs(self%cfg%xm(i-l)-self%cfg%x(i)).lt.abs(self%Wdist(i-l,j,k))) then
                      self%Wdist(i-l,j,k)=sign(self%cfg%xm(i-l)-self%cfg%x(i),self%Wdist(i-l,j,k))
                      self%Wnorm(:,i-l,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-l,j,k),0.0_WP,0.0_WP]
                   end if
                end if
             end do
          end do
       end do
       do k=self%cfg%kmino_,self%cfg%kmaxo_
          do j=self%cfg%jmino_+l,self%cfg%jmaxo_
             do i=self%cfg%imino_,self%cfg%imaxo_
                if (self%Wdist(i,j,k)*self%Wdist(i,j-l,k).lt.0.0_WP) then
                   ! There is a wall at y(j)
                   if (abs(self%cfg%ym(j  )-self%cfg%y(j-l+1)).lt.abs(self%Wdist(i,j  ,k))) then
                      self%Wdist(i,j  ,k)=sign(self%cfg%ym(j  )-self%cfg%y(j-l+1),self%Wdist(i,j  ,k))
                      self%Wnorm(:,i,j  ,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-l,k),0.0_WP]
                   end if
                   if (abs(self%cfg%ym(j-l)-self%cfg%y(j)).lt.abs(self%Wdist(i,j-l,k))) then
                      self%Wdist(i,j-l,k)=sign(self%cfg%ym(j-l)-self%cfg%y(j),self%Wdist(i,j-l,k))
                      self%Wnorm(:,i,j-l,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-l,k),0.0_WP]
                   end if
                end if
             end do
          end do
       end do
       do k=self%cfg%kmino_+l,self%cfg%kmaxo_
          do j=self%cfg%jmino_,self%cfg%jmaxo_
             do i=self%cfg%imino_,self%cfg%imaxo_
                if (self%Wdist(i,j,k)*self%Wdist(i,j,k-l).lt.0.0_WP) then
                   ! There is a wall at z(k)
                   if (abs(self%cfg%zm(k  )-self%cfg%z(k-l+1)).lt.abs(self%Wdist(i,j,k  ))) then
                      self%Wdist(i,j,k  )=sign(self%cfg%zm(k  )-self%cfg%z(k-l+1),self%Wdist(i,j,k  ))
                      self%Wnorm(:,i,j,k  )=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-l)]
                   end if
                   if (abs(self%cfg%zm(k-l)-self%cfg%z(k)).lt.abs(self%Wdist(i,j,k-l))) then
                      self%Wdist(i,j,k-l)=sign(self%cfg%zm(k-l)-self%cfg%z(k),self%Wdist(i,j,k-l))
                      self%Wnorm(:,i,j,k-l)=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-l)]
                   end if
                end if
             end do
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
           write(message,'("SPH object [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
           if (verbose.gt.1) write(output_unit,'(a)') trim(message)
           if (verbose.gt.0) call log(message)
        end if
      end block logging
      
    end function constructor


    !> Initialize the smoothing kernel and pressure
    !> Density, velocity, and smoothing length must be defined already
    subroutine sph_init(this)
      use messager, only: die
      implicit none
      class(sph), intent(inout) :: this
      if (this%h.gt.this%cfg%min_meshsize) call die('[sph init_kernel] Smoothing length must be < dx')
      if (this%rmax.gt.this%cfg%min_meshsize) call die('[sph init_kernel] Cutoff must be < dx')
      this%hi=1.0_WP/this%h
      call this%sph_update()
    end subroutine sph_init


    !> Update dependent SPH variables
    !> Calculate pressure from equation of state
    !> Sound speed is set to 10x max velocity to ensure low density variation
    !> Effective diameter computed from particle mass and density
    subroutine sph_update(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(sph), intent(inout) :: this
      integer :: i,ndim,ierr
      real(WP) :: P0,Umax,buf
      Umax=0.0_WP
      do i=1,this%np_
         Umax=max(Umax,sqrt(sum(this%p(i)%vel**2)))
      end do
      call MPI_ALLREDUCE(Umax,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); Umax=buf
      this%c=max(this%cmin,10.0_WP*Umax)
      P0=this%rho*this%c**2/this%gamm
      ndim=count((/this%cfg%nx,this%cfg%ny,this%cfg%nz/).gt.1)
      do i=1,this%np_
         this%p(i)%P=P0*((this%p(i)%rho/this%rho)**this%gamm-1.0_WP)
         this%p(i)%dp=(this%p(i)%vol*this%rho/this%p(i)%rho)**(1.0_WP/real(ndim,WP))
      end do
    end subroutine sph_update


    !> Calculate RHS of the particle ODEs
    subroutine get_rhs(this)
      implicit none
      class(sph), intent(inout) :: this
      integer, dimension(:,:,:),   allocatable :: npic    !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
      
      ! Communicate particles in ghost cells
      call this%share()

      ! Create image and boundary particles and store in ghost array
      ! Image particles are mapped across the boundary with equal/opposite velocity
      ! Boundary particles are placed on the boundary with zero velocity
      boundary: block
        integer :: i,nimg
        real(WP) :: dist
        real(WP), dimension(3) :: n12
        type(part), dimension(:), allocatable :: ptmp
        ! Count number of image/boundary particles
        nimg=0
        do i=1,this%np_
           if (this%p(i)%id.eq.0) cycle
           dist=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%Wdist,bc='d')
           if (dist.lt.this%rmax) nimg=nimg+1
        end do
        do i=1,this%ng_
           if (this%g(i)%id.eq.0) cycle
           dist=this%cfg%get_scalar(pos=this%g(i)%pos,i0=this%g(i)%ind(1),j0=this%g(i)%ind(2),k0=this%g(i)%ind(3),S=this%Wdist,bc='d')
           if (dist.lt.this%rmax) nimg=nimg+1
        end do
        ! Image particles across the boundary
        if (nimg.gt.0) then
           allocate(ptmp(this%ng_)); ptmp=this%g
           deallocate(this%g); allocate(this%g(this%ng_+nimg))
           this%g(1:this%ng_)=ptmp
           nimg=this%ng_
           do i=1,this%np_
              if (this%p(i)%id.eq.0) cycle
              dist=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%Wdist,bc='d')
              if (dist.lt.this%rmax) then
                 n12=this%Wnorm(:,this%p(i)%ind(1),this%p(i)%ind(2),this%p(i)%ind(3))
                 n12=n12/(sum(n12)+epsilon(1.0_WP))
                 nimg=nimg+1
                 this%g(nimg)=this%p(i)
                 this%g(nimg)%id=0
                 this%g(nimg)%vel=this%g(nimg)%vel-2.0_WP*dot_product(n12,this%p(i)%vel)*n12
                 this%g(nimg)%pos=this%p(i)%pos-2.0_WP*n12*dist
                 this%g(nimg)%ind=this%cfg%get_ijk_global(this%g(nimg)%pos,this%p(i)%ind)
              end if
           end do
           ! Image ghost particles across the boundary
           do i=1,this%ng_
              if (ptmp(i)%id.eq.0) cycle
              dist=this%cfg%get_scalar(pos=ptmp(i)%pos,i0=ptmp(i)%ind(1),j0=ptmp(i)%ind(2),k0=ptmp(i)%ind(3),S=this%Wdist,bc='d')
              if (dist.lt.this%rmax) then
                 n12=this%Wnorm(:,ptmp(i)%ind(1),ptmp(i)%ind(2),ptmp(i)%ind(3))
                 n12=n12/(sum(n12)+epsilon(1.0_WP))
                 nimg=nimg+1
                 this%g(nimg)=ptmp(i)
                 this%g(nimg)%id=0
                 this%g(nimg)%vel=-2.0_WP*dot_product(n12,ptmp(i)%vel)*n12
                 this%g(nimg)%pos=ptmp(i)%pos-2.0_WP*n12*dist
                 this%g(nimg)%ind=this%cfg%get_ijk_global(this%g(nimg)%pos,ptmp(i)%ind)
              end if
           end do
           this%ng_=nimg
           deallocate(ptmp)
        end if
      end block boundary

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
            ipic(npic(ip,jp,kp),ip,jp,kp)=+i
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=-i
         end do
       end block pic_prep


       ! Update right-hand side
       update_rhs: block
         use mathtools, only: inverse_matrix
         integer :: i,j,k,n1,nn,n2,ndim
         type(part) :: p1,p2
         real(WP), dimension(3) :: rpos,dWdx,n12,buf
         real(WP), dimension(3,3) :: L
         real(WP), dimension(:,:), allocatable :: Ldim,Li
         real(WP) :: dist,sumW,stress,m2,dV,Dcol,Lp

         ! Get number of dimensions
         ndim=count((/this%cfg%nx,this%cfg%ny,this%cfg%nz/).gt.1)
         allocate(Ldim(ndim,ndim))
         allocate(Li(ndim,ndim))

         ! Collision coefficient
         Dcol=2.0_WP*(0.1_WP*this%c)**2

         ! Loop over particles
         do n1=1,this%np_
            ! Particles marked 0 do not update their forces
            if (this%p(n1)%id.eq.0) cycle
            ! Create copy of our particle
            p1=this%p(n1)
            ! Zero out kernel normalization variables
            sumW=0.0_WP; L=0.0_WP
            ! Loop over neighbor cells
            do k=p1%ind(3)-1,p1%ind(3)+1
               do j=p1%ind(2)-1,p1%ind(2)+1
                  do i=p1%ind(1)-1,p1%ind(1)+1
                     ! Loop over particles in that cell
                     do nn=1,npic(i,j,k)
                        ! Create copy of our neighbor
                        n2=ipic(nn,i,j,k)
                        if (n2.gt.0) then
                           p2=this%p(+n2)
                        else if (n2.lt.0) then
                           p2=this%g(-n2)
                        end if
                        ! Current distance
                        rpos=p1%pos-p2%pos
                        dist=sqrt(dot_product(rpos,rpos))
                        if (dist.lt.this%rmax .and. p1%id.ne.p2%id) then
                           m2=p2%vol*this%rho
                           dV=m2/p2%rho
                           sumW = sumW + W(dist,this%h,cubic)*dV
                           dWdx = gradW(dist,this%h,p1%pos,p2%pos,cubic)
                           L(1,1) = L(1,1) + dV*(p2%pos(1)-p1%pos(1))*dWdx(1)
                           L(1,2) = L(1,2) + dV*(p2%pos(1)-p1%pos(1))*dWdx(2)
                           L(1,3) = L(1,3) + dV*(p2%pos(1)-p1%pos(1))*dWdx(3)
                           L(2,1) = L(2,1) + dV*(p2%pos(2)-p1%pos(2))*dWdx(1)
                           L(2,2) = L(2,2) + dV*(p2%pos(2)-p1%pos(2))*dWdx(2)
                           L(2,3) = L(2,3) + dV*(p2%pos(2)-p1%pos(2))*dWdx(3)
                           L(3,1) = L(3,1) + dV*(p2%pos(3)-p1%pos(3))*dWdx(1)
                           L(3,2) = L(3,2) + dV*(p2%pos(3)-p1%pos(3))*dWdx(2)
                           L(3,3) = L(3,3) + dV*(p2%pos(3)-p1%pos(3))*dWdx(3)
                        end if
                     end do
                  end do
               end do
            end do
            ! Invert L matrix for gradient correction
            Ldim=L(1:ndim,1:ndim)
            Li=inverse_matrix(Ldim)
            L=0.0_WP
            L(1:ndim,1:ndim)=Li
            ! Zero-out RHS terms
            p1%dRhoDt=0.0_WP; p1%dudt=0.0_WP; p1%dxdt=0.0_WP
            ! Loop over neighbor cells again
            do k=p1%ind(3)-1,p1%ind(3)+1
               do j=p1%ind(2)-1,p1%ind(2)+1
                  do i=p1%ind(1)-1,p1%ind(1)+1
                     ! Loop over particles in that cell
                     do nn=1,npic(i,j,k)
                        ! Create copy of our neighbor
                        n2=ipic(nn,i,j,k)
                        if (n2.gt.0) then
                           p2=this%p(+n2)
                        else if (n2.lt.0) then
                           p2=this%g(-n2)
                        end if
                        ! Current distance
                        rpos=p1%pos-p2%pos
                        dist=sqrt(dot_product(rpos,rpos))
                        ! Sum up weight if particle is within cutoff
                        if (dist.lt.this%rmax .and. p1%id.ne.p2%id) then
                           ! Compute mass
                           m2=p2%vol*this%rho
                           ! Compute the gradient
                           buf=gradW(dist,this%h,p1%pos,p2%pos,cubic)
                           dWdx(1) = sum(L(1,:)*buf)
                           dWdx(2) = sum(L(2,:)*buf)
                           dWdx(3) = sum(L(3,:)*buf)
                           ! Pressure
                           stress=-(p1%p/p1%rho**2+p2%p/p2%rho**2)
                           ! Viscous stress
                           stress=stress+16.0_WP*this%visc/(p1%rho+p2%rho)*dot_product(p1%vel-p2%vel,rpos)/this%h/dist
                           ! Sum RHS terms
                           p1%dudt=p1%dudt+stress*m2*dWdx
                           p1%dRhoDt=p1%dRhoDt+m2*dot_product(p1%vel-p2%vel,dWdx)
                           p1%dxdt=p1%dxdt+m2*(p2%vel-p1%vel)/(p1%rho+p2%rho)*W(dist,this%h,cubic)/sumW
                        end if
                     end do
                  end do
               end do
            end do
            ! Collisions with walls
            dist=this%cfg%get_scalar(pos=p1%pos,i0=p1%ind(1),j0=p1%ind(2),k0=p1%ind(3),S=this%Wdist,bc='d')
            Lp=(p1%vol)**(1.0_WP/real(ndim,WP))
            if (dist.lt.Lp) then
               n12=this%Wnorm(:,p1%ind(1),p1%ind(2),p1%ind(3))
               n12=n12/(sum(n12)+epsilon(1.0_WP))
               dist=dist+1.0e-9_WP*Lp
               p1%dudt=p1%dudt+Dcol*((Lp/dist)**4-(Lp/dist)**2)*n12/dist
            end if
            ! Include gravity
            p1%dudt=p1%dudt+this%gravity
            ! Deal with dimensionality
            if (this%cfg%nx.eq.1) then
               p1%dxdt(1)=0.0_WP
               p1%dudt(1)=0.0_WP
            end if
            if (this%cfg%ny.eq.1) then
               p1%dxdt(2)=0.0_WP
               p1%dudt(2)=0.0_WP
            end if
            if (this%cfg%nz.eq.1) then
               p1%dxdt(3)=0.0_WP
               p1%dudt(3)=0.0_WP
            end if
            ! Copy back the particle
            this%p(n1)=p1
         end do
         deallocate(Ldim,Li)
       end block update_rhs

       ! Clean up
       if (allocated(npic)) deallocate(npic)
       if (allocated(ipic)) deallocate(ipic)
      
     end subroutine get_rhs

   
   !> Advance the particle equations by a specified time step dt
   !> p%id=0 => do not advance
   subroutine advance(this,dt)
      use mpi_f08, only : MPI_SUM,MPI_INTEGER
      implicit none
      class(sph), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      integer :: n,ierr
      
      ! Zero out number of particles removed
      this%np_out=0

      ! Compute right-hand side terms
      call this%get_rhs()

      do n=1,this%np_
         ! Skip over particles with id=0
         if (this%p(n)%id.eq.0) cycle
         this%p(n)%rho=this%p(n)%rho+dt*this%p(n)%dRhoDt
         this%p(n)%pos=this%p(n)%pos+dt*(this%p(n)%vel+this%p(n)%dxdt)
         this%p(n)%vel=this%p(n)%vel+dt*this%p(n)%dudt
         ! Correct the position to take into account periodicity
         if (this%cfg%xper) this%p(n)%pos(1)=this%cfg%x(this%cfg%imin)+modulo(this%p(n)%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
         if (this%cfg%yper) this%p(n)%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(this%p(n)%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
         if (this%cfg%zper) this%p(n)%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(this%p(n)%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
         ! Handle particles that have left the domain
         if (this%p(n)%pos(1).lt.this%cfg%x(this%cfg%imin).or.this%p(n)%pos(1).gt.this%cfg%x(this%cfg%imax+1)) this%p(n)%flag=1
         if (this%p(n)%pos(2).lt.this%cfg%y(this%cfg%jmin).or.this%p(n)%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) this%p(n)%flag=1
         if (this%p(n)%pos(3).lt.this%cfg%z(this%cfg%kmin).or.this%p(n)%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) this%p(n)%flag=1
         ! Relocalize the particle
         this%p(n)%ind=this%cfg%get_ijk_global(this%p(n)%pos,this%p(n)%ind)
         ! Count number of particles removed
         if (this%p(n)%flag.eq.1) this%np_out=this%np_out+1
      end do

      ! Update dependent variables
      call this%sph_update()

      ! Communicate particles
      call this%sync()

      ! Sum up particles removed
      call MPI_ALLREDUCE(this%np_out,n,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=n
      
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
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(sph), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      integer :: i,ierr
      real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z,my_CFLp_f
      
      ! Set the CFLs to zero
      my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP; my_CFLp_f=0.0_WP
      do i=1,this%np_
         my_CFLp_x=max(my_CFLp_x,abs(this%p(i)%vel(1))*this%hi)
         my_CFLp_y=max(my_CFLp_y,abs(this%p(i)%vel(2))*this%hi)
         my_CFLp_z=max(my_CFLp_z,abs(this%p(i)%vel(3))*this%hi)
         my_CFLp_f=max(my_CFLp_f,sqrt(sum(this%p(i)%dudt**2)))
      end do
      my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt
      my_CFLp_f=4.0_WP*sqrt(2.0_WP*my_CFLp_f*this%hi)*dt
      this%CFLp_c=this%c*dt*this%hi
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLp_f,this%CFLp_f,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Return the maximum CFL
      cfl=max(this%CFLp_c,this%CFLp_f)
      
   end subroutine get_cfl
   
   
   !> Extract various monitoring data from particle field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(sph), intent(inout) :: this
      real(WP) :: buf,safe_np
      integer :: i,ierr
      
      ! Create safe np
      safe_np=real(max(this%np,1),WP)
      
      ! Density and velocity min/max/mean
      this%RHOmin=huge(1.0_WP); this%RHOmax=-huge(1.0_WP); this%RHOmean=0.0_WP
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); this%Umean=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); this%Vmean=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); this%Wmean=0.0_WP
      this%dRHO=0.0_WP
      do i=1,this%np_
         this%RHOmin=min(this%RHOmin,this%p(i)%rho); this%RHOmax=max(this%RHOmax,this%p(i)%rho); this%RHOmean=this%RHOmean+this%p(i)%rho
         this%Umin=min(this%Umin,this%p(i)%vel(1)); this%Umax=max(this%Umax,this%p(i)%vel(1)); this%Umean=this%Umean+this%p(i)%vel(1)
         this%Vmin=min(this%Vmin,this%p(i)%vel(2)); this%Vmax=max(this%Vmax,this%p(i)%vel(2)); this%Vmean=this%Vmean+this%p(i)%vel(2)
         this%Wmin=min(this%Wmin,this%p(i)%vel(3)); this%Wmax=max(this%Wmax,this%p(i)%vel(3)); this%Wmean=this%Wmean+this%p(i)%vel(3)
         this%dRHO=this%dRHO+abs(this%p(i)%rho-this%rho)
      end do
      call MPI_ALLREDUCE(this%RHOmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%RHOmin =buf
      call MPI_ALLREDUCE(this%RHOmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%RHOmax =buf
      call MPI_ALLREDUCE(this%RHOmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%RHOmean=buf/safe_np
      call MPI_ALLREDUCE(this%Umin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Umin =buf
      call MPI_ALLREDUCE(this%Umax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Umax =buf
      call MPI_ALLREDUCE(this%Umean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Umean=buf/safe_np
      call MPI_ALLREDUCE(this%Vmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vmin =buf
      call MPI_ALLREDUCE(this%Vmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vmax =buf
      call MPI_ALLREDUCE(this%Vmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vmean=buf/safe_np
      call MPI_ALLREDUCE(this%Wmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wmin =buf
      call MPI_ALLREDUCE(this%Wmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wmax =buf
      call MPI_ALLREDUCE(this%Wmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wmean=buf/safe_np
      call MPI_ALLREDUCE(this%dRHO,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dRHO=sqrt(buf)/this%rho/safe_np
      
   end subroutine get_max
   
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(sph), intent(inout) :: this
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
      if (ierr.ne.0) call die('[sph prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Share particles across processor boundaries
   subroutine share(this,nover)
      use mpi_f08
      use messager, only: warn,die
      implicit none
      class(sph), intent(inout) :: this
      integer, optional :: nover
      type(part), dimension(:), allocatable :: tosend
      type(part), dimension(:), allocatable :: torecv
      integer :: no,nsend,nrecv
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr
      integer :: n
      
      ! Check overlap size
      if (present(nover)) then
         no=nover
         if (no.gt.this%cfg%no) then
            call warn('[sph share] Specified overlap is larger than that of cfg - reducing no')
            no=this%cfg%no
         else if (no.le.0) then
            call die('[sph share] Specified overlap cannot be less or equal to zero')
         end if
      else
         no=1
      end if
      
      ! Clean up ghost array
      call this%resize_ghost(n=0); this%ng_=0
      
      ! Share ghost particles in -x (no ghosts are sent here)
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
      
      ! Share ghost particles in +x (no ghosts are sent here)
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
      
      ! Share ghost particles in -y (ghosts need to be sent now)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
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
      do n=1,this%ng_
         if (this%g(n)%ind(2).lt.this%cfg%jmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
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
      
      ! Share ghost particles in +y (ghosts need to be sent now - but not newly received ghosts!)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
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
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(2).gt.this%cfg%jmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
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
      
      ! Share ghost particles in -z (ghosts need to be sent now)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
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
      do n=1,this%ng_
         if (this%g(n)%ind(3).lt.this%cfg%kmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
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
      
      ! Share ghost particles in +z (ghosts need to be sent now - but not newly received ghosts!)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
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
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(3).gt.this%cfg%kmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
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
   
   
   !> Synchronize particle arrays across processors
   subroutine sync(this)
      use mpi_f08
      implicit none
      class(sph), intent(inout) :: this
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
      class(sph), intent(inout) :: this
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
      class(sph), intent(inout) :: this
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
      class(sph), intent(inout) :: this
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
      class(sph), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: i,ierr,iunit
      
      ! Root serial-writes the file header
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[sph write] Problem encountered while serial-opening data file: '//trim(filename))
         ! Number of particles and particle object size
         write(iunit) this%np,MPI_PART_SIZE
         ! Done with the header
         close(iunit)
      end if
      
      ! The rest is done in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[sph write] Problem encountered while parallel-opening data file: '//trim(filename))
      
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
            write(message,'("[sph write] Wrote ",i0," particles to file [",a,"] on partitioned grid [",a,"]")') this%np,trim(filename),trim(this%cfg%name)
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
      class(sph), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
      integer :: i,j,ierr,npadd,psize,nchunk,cnt
      integer, dimension(:,:), allocatable :: ppp
      
      ! First open the file in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[sph read] Problem encountered while reading data file: '//trim(filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
      
      ! Remember current position
      call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)
      
      ! Check compatibility of particle type
      if (psize.ne.MPI_PART_SIZE) call die('[sph read] Particle type unreadable')
      
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
            write(message,'("[sph read] Read ",i0," particles from file [",a,"] on partitioned grid [",a,"]")') npadd,trim(filename),trim(this%cfg%name)
            if (verbose.gt.2) write(output_unit,'(a)') trim(message)
            if (verbose.gt.1) call log(message)
         end if
      end block logging
      
   end subroutine read
   
   
end module sph_class
