!> Lagrangian solid solver object
!> Implements peridynamics equations
module lss_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: lss
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   

   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing
   

   !> Maximum number of bonds per particle
   integer, parameter, public :: max_bond=200  !< Assumes something like a 7x7x7 stencil in 3D
   

   !> Bonded solid particle definition
   type :: part
      !> MPI_DOUBLE_PRECISION data
      real(WP) :: mw                         !< Weighted volume
      real(WP) :: dil                        !< Element dilatation
      real(WP), dimension(max_bond) :: dbond !< Length of initial bonds
      real(WP), dimension(3) :: pos          !< Particle center coordinates
      real(WP), dimension(3) :: vel          !< Velocity of particle
      real(WP), dimension(3) :: Abond        !< Bond acceleration for particle
      real(WP), dimension(3) :: Afluid       !< Fluid acceleration for particle
      !> MPI_INTEGER data
      integer :: id                          !< ID the object is associated with
      integer :: i                           !< Unique index of particle (assumed >0)
      integer :: nbond                           !< Number of initial bonds
      integer, dimension(max_bond) :: ibond  !< Indices of initially bonded particles (0 values ignored)
      integer , dimension(3) :: ind          !< Index of cell containing particle center
      integer :: flag                        !< Control parameter (0=normal, 1=done->will be removed)
   end type part
   !> Number of blocks, block length, and block types in a particle
   integer, parameter                         :: part_nblock=2
   integer           , dimension(part_nblock) :: part_lblock=[14+max_bond,7+max_bond]
   type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_DOUBLE_PRECISION,MPI_INTEGER]
   !> MPI_PART derived datatype and size
   type(MPI_Datatype) :: MPI_PART
   integer :: MPI_PART_SIZE
   
   
   !> Lagrangian solid solver object definition
   type :: lss
   
      ! This config is used for parallelization and for calculating bond/collision forces
      class(config), pointer :: cfg
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_LSS'
      
      ! Solid material properties
      real(WP) :: elastic_modulus                         !< Elastic modulus of the material
      real(WP) :: poisson_ratio                           !< Poisson's ratio of the material
      real(WP) :: rho                                     !< Density of the material
      real(WP) :: crit_energy                             !< Critical energy release
      real(WP) :: dV                                      !< Element volume
      
      ! Bonding parameters
      real(WP) :: delta                                   !< Bonding horizon (distance)
      integer :: nb                                       !< Cell-based horizon 
      
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
      
      ! Volume fraction associated with IBM projection
      real(WP), dimension(:,:,:), allocatable :: VF       !< Volume fraction, cell-centered
      
      ! Momentum source
      real(WP), dimension(:,:,:), allocatable :: srcU     !< U momentum source on mesh, cell-centered
      real(WP), dimension(:,:,:), allocatable :: srcV     !< V momentum source on mesh, cell-centered
      real(WP), dimension(:,:,:), allocatable :: srcW     !< W momentum source on mesh, cell-centered
      
      ! CFL numbers
      real(WP) :: CFLp_x,CFLp_y,CFLp_z
      
      ! Number of substeps for time integrator
      real(WP) :: nstep=1
      
      ! Monitoring info
      real(WP) :: Umin,Umax,Umean                    !< U velocity info
      real(WP) :: Vmin,Vmax,Vmean                    !< V velocity info
      real(WP) :: Wmin,Wmax,Wmean                    !< W velocity info
      integer  :: np_out                             !< Number of particles leaving the domain
      
   contains
      procedure :: bond_init                         !< Setup initial interparticle bonds
      procedure :: get_bond_force                    !< Compute interparticle bond force
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
      procedure :: get_source                             !< Compute direct forcing source
      procedure :: update_VF                              !< Compute volume fraction
      procedure :: get_delta                              !< Compute regularized delta function
      procedure :: interpolate                            !< Interpolation routine from mesh=>marker
      procedure :: extrapolate                            !< Extrapolation routine from marker=>mesh
   end type lss
   
   
   !> Declare lss constructor
   interface lss
      procedure constructor
   end interface lss
   
contains
   
   
   ! Quasi-Gaussian weighting function - h is the cut-off
   real(WP) function wgauss(d,h)
      implicit none
      real(WP), intent(in) :: d,h
      real(WP), parameter :: coeff=2.6_WP
      real(WP) :: hh
      hh=coeff*h
      if (d.ge.hh) then
         wgauss=0.0_WP
      else
         wgauss=(1.0_WP+4.0_WP*d/hh)*(1.0_WP-d/hh)**4
      end if
   end function wgauss
   
   
   !> Default constructor for Lagrangian solid solver
   function constructor(cfg,name) result(self)
      implicit none
      type(lss) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Set default bonding horizon based on underlying mesh
      self%delta=self%cfg%min_meshsize
      self%nb=1
      
      ! Allocate variables
      allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
      self%np_=0; self%np=0
      call self%resize(0)
      
      ! Initialize MPI derived datatype for a particle
      call prepare_mpi_part()
      
      ! Allocate levelset, VF and src arrays on cfg mesh
      allocate(self%VF  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF  =0.0_WP
      allocate(self%srcU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcU=0.0_WP
      allocate(self%srcV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcV=0.0_WP
      allocate(self%srcW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcW=0.0_WP
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%cfg%amRoot) then
            write(message,'("LSS object [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function constructor
   
   
   !> Initialize bond force between particles
   subroutine bond_init(this)
      use messager, only: die
      implicit none
      class(lss), intent(inout) :: this
      integer, dimension(:,:,:),   allocatable :: npic    !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
      
      ! Communicate particles in ghost cells
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
            ipic(npic(ip,jp,kp),ip,jp,kp)=+i
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=-i
         end do
         
      end block pic_prep
      
      ! Establish initial bonds
      create_bonds: block
         integer :: i,j,k,n1,nn,n2
         type(part) :: p1,p2
         real(WP), dimension(3) :: rpos
         real(WP) :: dist

         ! Loop over particles
         do n1=1,this%np_
            ! Create copy of our particle
            p1=this%p(n1)
            ! Zero out weighted volume
            p1%mw=0.0_WP
            ! Zero out bonds
            p1%ibond=0
            p1%nbond=0
            ! Loop over neighbor cells
            do k=p1%ind(3)-this%nb,p1%ind(3)+this%nb
               do j=p1%ind(2)-this%nb,p1%ind(2)+this%nb
                  do i=p1%ind(1)-this%nb,p1%ind(1)+this%nb
                     ! Loop over particles in that cell
                     do nn=1,npic(i,j,k)
                        ! Create copy of our neighbor
                        n2=ipic(nn,i,j,k)
                        if (n2.gt.0) then
                           p2=this%p(+n2)
                        else if (n2.lt.0) then
                           p2=this%g(-n2)
                        end if
                        ! Check interparticle distance
                        rpos=p2%pos-p1%pos
                        dist=sqrt(dot_product(rpos,rpos))
                        if (dist.lt.this%delta) then
                           ! Cannot self-bond
                           if (p1%i.eq.p2%i) cycle
                           ! Cannot bond with different id except <=0 (<=0 bonds with everyone)
                           if (p1%id.ne.p2%id.and.p1%id.ge.0.and.p2%id.ge.0) cycle
                           ! This particle is in horizon, create a bond
                           p1%nbond=p1%nbond+1
                           if (p1%nbond.gt.max_bond) call die('[lss_class bond_init] Number of detected bonds is larger than max allowed')
                           p1%ibond(p1%nbond)=p2%i
                           p1%dbond(p1%nbond)=dist
                           ! Increment weighted volume
                           p1%mw=p1%mw+wgauss(dist,this%delta)*dist**2*this%dV
                        end if
                     end do
                  end do
               end do
            end do
            ! Zero out initial dilatation
            p1%dil=0.0_WP
            ! Copy back the particle
            this%p(n1)=p1
         end do
      end block create_bonds
      
      ! Clean up
      if (allocated(npic)) deallocate(npic)
      if (allocated(ipic)) deallocate(ipic)

   end subroutine bond_init


   !> Calculate bond force between particles
   subroutine get_bond_force(this)
      implicit none
      class(lss), intent(inout) :: this
      integer, dimension(:,:,:),   allocatable :: npic    !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
      
      ! Communicate particles in ghost cells
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
            ipic(npic(ip,jp,kp),ip,jp,kp)=+i
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=-i
         end do
         
      end block pic_prep
      
      ! Update weighted volume and dilatation
      update_weighted_vol_and_dilatation: block
         integer :: i,j,k,n1,nn,n2
         type(part) :: p1,p2
         integer :: nb,nbond
         real(WP), dimension(3) :: rpos
         real(WP) :: dist
         
         ! Loop over particles
         do n1=1,this%np_
            ! Create copy of our particle
            p1=this%p(n1)
            ! Zero out weighted volume and dilatation
            p1%mw=0.0_WP
            p1%dil=0.0_WP
            ! Loop over neighbor cells
            do k=p1%ind(3)-this%nb,p1%ind(3)+this%nb
               do j=p1%ind(2)-this%nb,p1%ind(2)+this%nb
                  do i=p1%ind(1)-this%nb,p1%ind(1)+this%nb
                     ! Loop over particles in that cell
                     do nn=1,npic(i,j,k)
                        ! Create copy of our neighbor
                        n2=ipic(nn,i,j,k)
                        if (n2.gt.0) then
                           p2=this%p(+n2)
                        else if (n2.lt.0) then
                           p2=this%g(-n2)
                        end if
                        ! Check if a bond exists
                        do nb=1,max_bond
                           if (p1%ibond(nb).eq.p2%i) then
                              ! Increment weighted volume
                              p1%mw=p1%mw+wgauss(p1%dbond(nb),this%delta)*p1%dbond(nb)**2*this%dV
                              ! Get current distance
                              rpos=p2%pos-p1%pos
                              dist=sqrt(dot_product(rpos,rpos))
                              ! Increment dilatation
                              p1%dil=p1%dil+wgauss(p1%dbond(nb),this%delta)*p1%dbond(nb)*(dist-p1%dbond(nb))*this%dV
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            ! Rescale dilatation
            p1%dil=p1%dil*3.0_WP/p1%mw
            ! Copy back the particle
            this%p(n1)=p1
         end do
      end block update_weighted_vol_and_dilatation
      
      ! Re-communicate particles in ghost cells to update dil and mw
      call this%share()
      
      ! Update bond force, including collision force
      update_bond_force: block
         use mathtools, only: Pi
         integer :: i,j,k,n1,nn,n2
         type(part) :: p1,p2
         real(WP), dimension(3) :: rpos,t12,t21
         real(WP) :: dist,beta,alpha,ed
         real(WP) :: stretch,max_stretch,mu,kk
         real(WP) :: nc,rc,kc
         integer :: nb,nbond
         logical :: found_bond
         
         ! Recompute a few physical parameters
         mu=this%elastic_modulus/(2.0_WP+2.0_WP*this%poisson_ratio)
         kk=this%elastic_modulus/(3.0_WP-6.0_WP*this%poisson_ratio)
         max_stretch=sqrt(this%crit_energy/((3.0_WP*mu+(kk-5.0_WP*mu/3.0_WP)*0.75_WP**4)*this%delta))
         rc=this%dV**(1.0_WP/3.0_WP)
         nc=1.0_WP
         kc=15.0_WP*12.0_WP*this%elastic_modulus/(Pi*this%delta**4)
         
         ! Loop over particles
         do n1=1,this%np_
            ! Particles marked 0 do not update their forces
            if (this%p(n1)%id.eq.0) cycle
            ! Create copy of our particle
            p1=this%p(n1)
            ! Zero out bond force
            p1%Abond=0.0_WP
            ! Loop over neighbor cells
            do k=p1%ind(3)-this%nb,p1%ind(3)+this%nb
               do j=p1%ind(2)-this%nb,p1%ind(2)+this%nb
                  do i=p1%ind(1)-this%nb,p1%ind(1)+this%nb
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
                        rpos=p2%pos-p1%pos
                        dist=sqrt(dot_product(rpos,rpos))
                        ! Check if a bond exists
                        found_bond=.false.
                        do nb=1,max_bond
                           if (p1%ibond(nb).eq.p2%i) then
                              ! Check for breakage first
                              stretch=(dist-p1%dbond(nb))/p1%dbond(nb)
                              if (stretch.gt.max_stretch) then
                                 ! Remove the bond
                                 p1%ibond(nb)=0
                                 p1%dbond(nb)=0.0_WP
                                 cycle
                              end if
                              ! Beta1
                              beta=3.0_WP*kk*p1%dil
                              ! Alpha1
                              alpha=15.0_WP*mu/p1%mw
                              ! Extension1
                              ed=dist-p1%dbond(nb)*(1.0_WP+p1%dil/3.0_WP)
                              ! Force density 1->2
                              t12=+wgauss(p1%dbond(nb),this%delta)*(beta/p1%mw*p1%dbond(nb)+alpha*ed)*rpos/dist
                              ! Beta2
                              beta=3.0_WP*kk*p2%dil
                              ! Alpha2
                              alpha=15.0_WP*mu/p2%mw
                              ! Extension2
                              ed=dist-p1%dbond(nb)*(1.0_WP+p2%dil/3.0_WP)
                              ! Force density 2->1
                              t21=-wgauss(p1%dbond(nb),this%delta)*(beta/p2%mw*p1%dbond(nb)+alpha*ed)*rpos/dist
                              ! Increment bond force
                              p1%Abond=p1%Abond+(t12-t21)*this%dV/this%rho
                              ! If still here, we have an active bond
                              found_bond=.true.
                              cycle
                           end if
                        end do
                        ! Add collision force now
                        if (.not.found_bond.and.p1%i.ne.p2%i.and.dist.lt.rc) then
                           p1%Abond=p1%Abond-kc*((rc/dist)**nc-1.0_WP)*(rpos/dist)*this%dV/this%rho
                        end if
                     end do
                  end do
               end do
            end do
            ! Copy back the particle
            this%p(n1)=p1
         end do
      end block update_bond_force
      
      ! Clean up
      if (allocated(npic)) deallocate(npic)
      if (allocated(ipic)) deallocate(ipic)
      
   end subroutine get_bond_force

   
   !> Advance the particle equations by a specified time step dt
   !> p%id=-2 => do not solve for position nor velocity
   !> p%id=-1 => do not solve for velocity
   !> p%id= 0 => do not update force
   subroutine advance(this,dt)
      use mpi_f08,   only : MPI_SUM,MPI_INTEGER
      use mathtools, only: Pi
      implicit none
      class(lss), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      integer :: n,ierr
      type(part) :: myp,pold
      
      ! Zero out number of particles removed
      this%np_out=0
      
      ! Advance velocity based on old force and position based on mid-velocity
      do n=1,this%np_
         ! Advance with Verlet scheme
         if (this%p(n)%id.gt.-1) this%p(n)%vel=this%p(n)%vel+0.5_WP*dt*(this%gravity+this%p(n)%Abond+this%p(n)%Afluid)
         if (this%p(n)%id.gt.-2) this%p(n)%pos=this%p(n)%pos+dt*this%p(n)%vel
         ! Relocalize
         this%p(n)%ind=this%cfg%get_ijk_global(this%p(n)%pos,this%p(n)%ind)
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
      
      ! Communicate particles
      call this%sync()
      
      ! Sum up particles removed
      call MPI_ALLREDUCE(this%np_out,n,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=n
      
      ! Calculate bond force
      call this%get_bond_force()
      
      ! Advance velocity only based on new force
      do n=1,this%np_
         ! Advance with Verlet scheme
         if (this%p(n)%id.gt.-1) this%p(n)%vel=this%p(n)%vel+0.5_WP*dt*(this%gravity+this%p(n)%Abond+this%p(n)%Afluid)
      end do
      
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
   

   !> Compute direct forcing source by a specified time step dt
   subroutine get_source(this,dt,U,V,W,rho)
      implicit none
      class(lss), intent(inout) :: this
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
      
      ! Advance the equations
      dti=1.0_WP/dt
      do i=1,this%np_
         ! Zero out fluid force
         this%p(i)%Afluid=0.0_WP
         ! Interpolate the velocity to the particle location
         vel=0.0_WP
         if (this%cfg%nx.gt.1) vel(1)=this%interpolate(A=U,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='U')
         if (this%cfg%ny.gt.1) vel(2)=this%interpolate(A=V,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='V')
         if (this%cfg%nz.gt.1) vel(3)=this%interpolate(A=W,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='W')
         rho_=this%interpolate(A=rho,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),dir='SC')
         ! Compute the source term
         src=(this%p(i)%vel-vel)*this%dV
         ! Send source term back to the mesh
         if (this%cfg%nx.gt.1) call this%extrapolate(Ap=src(1),xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%srcU,dir='U')
         if (this%cfg%ny.gt.1) call this%extrapolate(Ap=src(2),xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%srcV,dir='V')
         if (this%cfg%nz.gt.1) call this%extrapolate(Ap=src(3),xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%srcW,dir='W')
         ! Form the force on particle
         this%p(i)%Afluid=-rho_*src*dti/(this%rho*this%dV)
      end do
      
      ! Sum at boundaries
      call this%cfg%syncsum(this%srcU)
      call this%cfg%syncsum(this%srcV)
      call this%cfg%syncsum(this%srcW)
      
      ! Recompute volume fraction
      call this%update_VF()
      
   end subroutine get_source

   
   !> Update particle volume fraction using our current particles
   subroutine update_VF(this)
      implicit none
      class(lss), intent(inout) :: this
      integer :: i
      ! Reset volume fraction
      this%VF=0.0_WP
      ! Transfer particle volume
      do i=1,this%np_
         ! Skip inactive particle
         if (this%p(i)%flag.eq.1) cycle
         ! Transfer volume to mesh
         call this%extrapolate(Ap=this%dV,xp=this%p(i)%pos(1),yp=this%p(i)%pos(2),zp=this%p(i)%pos(3),ip=this%p(i)%ind(1),jp=this%p(i)%ind(2),kp=this%p(i)%ind(3),A=this%VF,dir='SC')
      end do
      ! Sum at boundaries
      call this%cfg%syncsum(this%VF)
   end subroutine update_VF
   

   !> Compute regularized delta function
   subroutine get_delta(this,delta,ic,jc,kc,xp,yp,zp,dir)
      implicit none
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
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
      
      ! Return the maximum CFL
      cfl=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)
      
   end subroutine get_cfl
   
   
   !> Extract various monitoring data from particle field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(lss), intent(inout) :: this
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
      
   end subroutine get_max
   
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(lss), intent(inout) :: this
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
      if (ierr.ne.0) call die('[lss prepare_mpi_part] MPI Particle type creation failed')
      ! Get the size of this type
      call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   end subroutine prepare_mpi_part
   
   
   !> Share particles across processor boundaries
   subroutine share(this,nover)
      use mpi_f08
      use messager, only: warn,die
      implicit none
      class(lss), intent(inout) :: this
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
            call warn('[lss share] Specified overlap is larger than that of cfg - reducing no')
            no=this%cfg%no
         else if (no.le.0) then
            call die('[lss share] Specified overlap cannot be less or equal to zero')
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
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
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
      class(lss), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: i,ierr,iunit
      
      ! Root serial-writes the file header
      if (this%cfg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[lss write] Problem encountered while serial-opening data file: '//trim(filename))
         ! Number of particles and particle object size
         write(iunit) this%np,MPI_PART_SIZE
         ! Done with the header
         close(iunit)
      end if
      
      ! The rest is done in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[lss write] Problem encountered while parallel-opening data file: '//trim(filename))
      
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
            write(message,'("[lss write] Wrote ",i0," particles to file [",a,"] on partitioned grid [",a,"]")') this%np,trim(filename),trim(this%cfg%name)
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
      class(lss), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
      integer :: i,j,ierr,npadd,psize,nchunk,cnt
      integer, dimension(:,:), allocatable :: ppp
      
      ! First open the file in parallel
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[lss read] Problem encountered while reading data file: '//trim(filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
      
      ! Remember current position
      call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)
      
      ! Check compatibility of particle type
      if (psize.ne.MPI_PART_SIZE) call die('[lss read] Particle type unreadable')
      
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
            write(message,'("[lss read] Read ",i0," particles from file [",a,"] on partitioned grid [",a,"]")') npadd,trim(filename),trim(this%cfg%name)
            if (verbose.gt.2) write(output_unit,'(a)') trim(message)
            if (verbose.gt.1) call log(message)
         end if
      end block logging
      
   end subroutine read
   
   
end module lss_class
