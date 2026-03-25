!> AMR-aware Lagrangian particle tracking solver
!> Mirrors lpt_class capabilities on an AMReX AMR hierarchy.
!> Particle communication and sorting handled by NeighborParticleContainer<14,1>.
module amrlpt_class
   use precision, only: WP,I8
   use string, only: str_medium
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use iso_c_binding
   implicit none
   private

   ! Public exports
   public :: amrlpt,part
   public :: PART_MOVES,PART_COLLIDES,PART_EXCHANGES,PART_IS_DEAD

   ! Particle struct layout constants (must match #define in amrlpt_wrapper.cpp)
   integer, parameter, public :: AMRLPT_NREAL=14  !< extra reals per particle
   integer, parameter, public :: AMRLPT_NINT =1   !< extra ints  per particle

   ! Particle flags
   integer(c_int), parameter, public :: PART_MOVES     = 1  ! bit 0
   integer(c_int), parameter, public :: PART_COLLIDES  = 2  ! bit 1
   integer(c_int), parameter, public :: PART_EXCHANGES = 4  ! bit 2
   integer(c_int), parameter, public :: PART_IS_DEAD   = 0  ! no bits = remove

   ! Particle struct -- must match C++ Particle<14,1> memory layout exactly:
   ! - pos[3]    (pos, managed by AMReX)
   ! - rdata[14] (d, vel[3], angVel[3], Acol[3], Tcol[3], dt)
   ! - idcpu     (packed id+cpu, private)
   ! - idata[1]  (flag)
   type, bind(C), public :: part
      !> AMReX position (physical coordinates)
      real(c_double) :: pos(3)
      !> Extra reals (NStructReal=14)
      real(c_double) :: d               !< Particle diameter
      real(c_double) :: vel(3)          !< Particle velocity
      real(c_double) :: angVel(3)       !< Angular velocity
      real(c_double) :: Acol(3)         !< Collision acceleration
      real(c_double) :: Tcol(3)         !< Collision torque
      real(c_double) :: dt              !< Particle sub-timestep
      !> Packed id (39 bits, sign=valid) + cpu (24 bits)  [AMReX internal]
      integer(c_int64_t), private :: idcpu
      !> Extra int (NStructInt=1)
      integer(c_int) :: flag            !< 0=active, 1=mark for removal
   end type part

   !> C interface: bind(C) declarations to amrlpt_wrapper.cpp
   interface

      subroutine amrlpt_new_pc(pc, amrcore) bind(c)
         import :: c_ptr
         type(c_ptr) :: pc
         type(c_ptr), value :: amrcore
      end subroutine

      subroutine amrlpt_delete_pc(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_redistribute(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine

      subroutine amrlpt_fill_neighbors(pc,ngrow) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: ngrow
      end subroutine

      subroutine amrlpt_clear_neighbors(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_get_neighbor_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_build_neighbor_list(pc,rcrit) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: rcrit
      end subroutine

      subroutine amrlpt_get_neighbor_list_mfi(pc,lev,mfi,pairs,npairs) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: pairs
         integer(c_int64_t) :: npairs
      end subroutine

      subroutine amrlpt_get_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_add_particle_i(pc,lev,grid,tile,p) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,p
         integer(c_int), value :: lev,grid,tile
      end subroutine

      subroutine amrlpt_get_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine

      subroutine amrlpt_set_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t), value :: id
      end subroutine

      subroutine amrlpt_get_cpu(cpu) bind(c)
         import :: c_int
         integer(c_int) :: cpu
      end subroutine

      subroutine amrlpt_set_particle_id(id,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t), value :: id
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_set_particle_cpu(cpu,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int), value :: cpu
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_particle_is_valid(valid,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int) :: valid
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_total_np(pc,np) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_write(pc,path,is_chk) bind(c)
         import :: c_ptr,c_char,c_int
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
         integer(c_int), value :: is_chk
      end subroutine

      subroutine amrlpt_read(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
      end subroutine

   end interface

   !> AMR LPT solver type
   type :: amrlpt

      !> Associated AMR grid
      class(amrgrid), pointer :: amr=>null()

      !> AMReX NeighborParticleContainer<14,1> opaque handle
      type(c_ptr) :: pc=c_null_ptr

      !> Solver name
      character(len=str_medium) :: name='UNNAMED_AMRLPT'

      !> Global particle count
      integer(c_int64_t) :: np=0

      !> Physics parameters
      real(WP) :: rho                                   !< Particle material density
      real(WP), dimension(3) :: gravity=0.0_WP          !< Acceleration of gravity
      integer :: nstep=1                                !< Substeps per timestep
      character(len=str_medium) :: drag_model='Tenneti' !< Drag model

      !> Collision parameters
      real(WP) :: tau_col=1.0e-3_WP             !< Collision time scale
      real(WP) :: e_n=1.0_WP                    !< Normal restitution
      real(WP) :: e_w=1.0_WP                    !< Wall restitution
      real(WP) :: mu_f=0.0_WP                   !< Friction coefficient
      real(WP) :: clip_col=0.2_WP               !< Max overlap fraction

      !> CFL numbers
      real(WP) :: CFLp_x=0.0_WP,CFLp_y=0.0_WP,CFLp_z=0.0_WP
      real(WP) :: CFL_col=0.0_WP

      !> Monitoring data
      real(WP) :: dmin,dmax,dmean,dvar
      real(WP) :: Umin,Umax,Umean,Uvar
      real(WP) :: Vmin,Vmax,Vmean,Vvar
      real(WP) :: Wmin,Wmax,Wmean,Wvar
      real(WP) :: VFmin,VFmax,VFmean,VFvar
      integer  :: np_new=0,np_out=0
      real(WP) :: Vp_new=0.0_WP,Vp_out=0.0_WP,Vp_tot=0.0_WP
      integer  :: ncol=0

      !> Overlap size
      integer :: nover=2

      !> Filter width
      real(WP) :: filter_width=0.0_WP

      !> Particle volume fraction
      type(amrdata) :: VF

      !> Two-way coupling source terms
      type(amrdata) :: src

   contains
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Physics procedures
      procedure :: advance                !< Advance particle ODEs one timestep
      procedure :: update_VF              !< Compute particle volume fraction field
      procedure :: get_cfl                !< Compute particle CFL numbers
      ! Utilities
      procedure :: redistribute           !< Call AMReX redistribute
      procedure :: fill_ghosts            !< Fill ghost particle buffer
      procedure :: clear_ghosts           !< Release ghost particle buffer
      procedure :: build_neighbor_list    !< Build explicit pair list within rcrit
      procedure :: get_np                 !< Update global particle count
      procedure, private :: get_particles !< Get particle array for MFIter tile
      procedure :: interp                 !< Trilinear cell-centered interpolation
      procedure :: interp_face_velocities !< Trilinear face-centered interpolation
      procedure, private :: filter        !< Explicit diffusion filter
      ! Print solver info
      procedure :: get_info
      procedure :: print
      ! Checkpoint I/O
      procedure :: read
      procedure :: write
   end type amrlpt

contains

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize amrlpt solver
   subroutine initialize(this,amr,name)
      use amrdata_class, only: amrex_bc_foextrap
      implicit none
      class(amrlpt), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), optional :: name
      ! Assign solver name
      if (present(name)) this%name=trim(adjustl(name))
      ! Point to amr
      this%amr=>amr
      ! Create particle container
      call amrlpt_new_pc(this%pc,this%amr%amrcore)
      ! Initialize VF field
      call this%VF%initialize(amr=amr,name='VF',ncomp=1,ng=this%nover); call this%VF%register()
      where (.not.[this%amr%xper,this%amr%yper,this%amr%zper])
         this%VF%lo_bc(1:3,1)=amrex_bc_foextrap
         this%VF%hi_bc(1:3,1)=amrex_bc_foextrap
      end where
      ! Initialize source terms
      call this%src%initialize(amr=amr,name='src',ncomp=3,ng=this%nover); call this%src%register()
      where (.not.[this%amr%xper,this%amr%yper,this%amr%zper])
         this%src%lo_bc(1:3,:)=amrex_bc_foextrap
         this%src%hi_bc(1:3,:)=amrex_bc_foextrap
      end where
      ! Print out info
      call this%print()
   end subroutine initialize

   !> Finalize: destroy particle container and release grid pointer
   subroutine finalize(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call this%VF%finalize()
      call this%src%finalize()
      call amrlpt_delete_pc(this%pc)
      this%pc=c_null_ptr
      nullify(this%amr)
   end subroutine finalize

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Advance all particles by dt using Euler-midpoint substepping
   !> U,V,W: velocity amrdata (staggered MAC or collocated, with ghosts filled)
   !> rho,visc: cell-centered amrdata (with ghosts filled)
   !> Ucomp/Vcomp/Wcomp: component to use from each velocity field (optional, default 1)
   subroutine advance(this,dt,U,Ucomp,V,Vcomp,W,Wcomp,rho,visc)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mathtools, only: Pi
      use messager,  only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in) :: U,V,W,rho,visc
      integer, intent(in), optional :: Ucomp,Vcomp,Wcomp

      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRho,pVisc,pVolFrac,pSrc
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: lvl,i,uc,vc,wc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,mydt,dt_done,Ip
      real(WP), dimension(3) :: acc,dmom
      type(part) :: myp,pold
      logical :: is_stag

      ! Resolve optional component indices
      uc=1; if (present(Ucomp)) uc=Ucomp
      vc=1; if (present(Vcomp)) vc=Vcomp
      wc=1; if (present(Wcomp)) wc=Wcomp

      ! Check velocity nodal locations are consistent
      check_velocity: block
         logical, dimension(3) :: nU,nV,nW
         nU=U%nodal; nV=V%nodal; nW=W%nodal
         if (all(nU.eqv.[.true.,.false.,.false.]).and.&
         &   all(nV.eqv.[.false.,.true.,.false.]).and.&
         &   all(nW.eqv.[.false.,.false.,.true.])) then
            is_stag=.true.
         else if (.not.any(nU).and..not.any(nV).and..not.any(nW)) then
            is_stag=.false.
         else
            call die('[amrlpt advance] U/V/W must be staggered (face-centered) or collocated (cell-centered)')
         end if
      end block check_velocity

      ! Zero source mfabs
      call this%src%setval(0.0_WP)

      ! Track number of particles leaving domain
      this%np_out=0
      this%Vp_out=0.0_WP

      ! Loop over all AMR levels
      do lvl=0,this%amr%clvl()

         ! Get mesh size
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz

         ! MFIter over the level
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())

            ! Get pointers to data
            pU      =>U%mf(lvl)%dataptr(mfi)
            pV      =>V%mf(lvl)%dataptr(mfi)
            pW      =>W%mf(lvl)%dataptr(mfi)
            pRho    =>rho%mf(lvl)%dataptr(mfi)
            pVisc   =>visc%mf(lvl)%dataptr(mfi)
            pVolFrac=>this%VF%mf(lvl)%dataptr(mfi)
            pSrc    =>this%src%mf(lvl)%dataptr(mfi)

            ! Get particles on this tile
            call this%get_particles(lvl=lvl,mfi=mfi,p=p,np=np_)
            
            ! Loop over local particles
            do i=1,np_
               ! Skip particles that are not moving or exchanging
               if (IAND(p(i)%flag,PART_MOVES+PART_EXCHANGES).eq.0) cycle
               ! Create copy of particle
               myp=p(i)
               ! Time-integrate until dt_done=dt
               dt_done=0.0_WP
               do while (dt_done.lt.dt)
                  mydt=min(myp%dt,dt-dt_done)
                  if (mydt.le.0.0_WP) mydt=dt-dt_done
                  ! Remember the particle
                  pold=myp
                  ! Precompute moment of inertia for a sphere
                  Ip=0.1_WP*myp%d**2
                  ! Advance with Euler prediction
                  acc=get_rhs()
                  if (IAND(myp%flag,PART_MOVES).ne.0) then
                     myp%pos=pold%pos+0.5_WP*mydt*myp%vel
                     myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity+myp%Acol)
                     myp%angVel=pold%angVel+0.5_WP*mydt*myp%Tcol/Ip
                  end if
                  acc=get_rhs()
                  if (IAND(myp%flag,PART_MOVES).ne.0) then
                     myp%pos=pold%pos+mydt*myp%vel
                     myp%vel=pold%vel+mydt*(acc+this%gravity+myp%Acol)
                     myp%angVel=pold%angVel+mydt*myp%Tcol/Ip
                  end if
                  ! Transfer back to the mesh
                  if (IAND(myp%flag,PART_EXCHANGES).ne.0) then
                     dmom=mydt*acc*this%rho*Pi/6.0_WP*myp%d**3
                     call deposit(val=-dmom)
                  end if
                  ! Increment
                  dt_done=dt_done+mydt
               end do
               ! Track escape from non-periodic boundaries
               if ((.not.this%amr%xper.and.(myp%pos(1).lt.this%amr%xlo.or.myp%pos(1).gt.this%amr%xhi)).or.&
               &   (.not.this%amr%yper.and.(myp%pos(2).lt.this%amr%ylo.or.myp%pos(2).gt.this%amr%yhi)).or.&
               &   (.not.this%amr%zper.and.(myp%pos(3).lt.this%amr%zlo.or.myp%pos(3).gt.this%amr%zhi))) then
                  this%np_out=this%np_out+1
                  this%vp_out=this%vp_out+Pi/6.0_WP*myp%d**3
               end if
               ! Write back
               p(i)=myp
            end do
         end do
         call this%amr%mfiter_destroy(mfi)

      end do

      ! Redistribute particles
      call this%redistribute()

      ! Reduce info on particles leaving domain
      reduce_leaving_particles: block
         use mpi_f08,  only: MPI_INTEGER,MPI_IN_PLACE,MPI_SUM,MPI_ALLREDUCE
         use parallel, only: MPI_REAL_WP
         integer :: ierr
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_out,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_out,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block reduce_leaving_particles

      ! Accumulate ghost→valid, restrict-SUM fine into coarse, filter
      call this%src%syncsum(); call this%src%sum_down(); call this%src%average_down(); call this%src%fill(time=0.0_WP)
      call this%filter(this%src)

      ! Recompute particle volume fraction
      call this%update_VF()

   contains

      !> Calculate rhs of particle equations of motion
      function get_rhs() result(acc)
         real(WP), dimension(3) :: acc,fvel
         real(WP) :: pVF,fVF,frho,fvisc,Re,tau,corr,b1,b2
         real(WP) :: wx,wy,wz,fx,fy,fz
         integer  :: ic,jc,kc,ix,iy,iz
         ! Cell-centered stencil indices and weights
         ic=floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(myp%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
         jc=floor((myp%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(myp%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
         kc=floor((myp%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(myp%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
         ! Velocity interpolation
         if (is_stag) then
            ! Face-centered indices and weights
            ix=floor((myp%pos(1)-this%amr%xlo)*dxi); fx=(myp%pos(1)-this%amr%xlo)*dxi-real(ix,WP)
            iy=floor((myp%pos(2)-this%amr%ylo)*dyi); fy=(myp%pos(2)-this%amr%ylo)*dyi-real(iy,WP)
            iz=floor((myp%pos(3)-this%amr%zlo)*dzi); fz=(myp%pos(3)-this%amr%zlo)*dzi-real(iz,WP)
            fvel(1)=(1.0_WP-fx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix,jc,kc,uc)+fx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix+1,jc,kc,uc)+(1.0_WP-fx)*wy*(1.0_WP-wz)*pU(ix,jc+1,kc,uc)+fx*wy*(1.0_WP-wz)*pU(ix+1,jc+1,kc,uc)+(1.0_WP-fx)*(1.0_WP-wy)*wz*pU(ix,jc,kc+1,uc)+fx*(1.0_WP-wy)*wz*pU(ix+1,jc,kc+1,uc)+(1.0_WP-fx)*wy*wz*pU(ix,jc+1,kc+1,uc)+fx*wy*wz*pU(ix+1,jc+1,kc+1,uc)
            fvel(2)=(1.0_WP-wx)*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic,iy,kc,vc)+wx*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic+1,iy,kc,vc)+(1.0_WP-wx)*fy*(1.0_WP-wz)*pV(ic,iy+1,kc,vc)+wx*fy*(1.0_WP-wz)*pV(ic+1,iy+1,kc,vc)+(1.0_WP-wx)*(1.0_WP-fy)*wz*pV(ic,iy,kc+1,vc)+wx*(1.0_WP-fy)*wz*pV(ic+1,iy,kc+1,vc)+(1.0_WP-wx)*fy*wz*pV(ic,iy+1,kc+1,vc)+wx*fy*wz*pV(ic+1,iy+1,kc+1,vc)
            fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic,jc,iz,wc)+wx*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic+1,jc,iz,wc)+(1.0_WP-wx)*wy*(1.0_WP-fz)*pW(ic,jc+1,iz,wc)+wx*wy*(1.0_WP-fz)*pW(ic+1,jc+1,iz,wc)+(1.0_WP-wx)*(1.0_WP-wy)*fz*pW(ic,jc,iz+1,wc)+wx*(1.0_WP-wy)*fz*pW(ic+1,jc,iz+1,wc)+(1.0_WP-wx)*wy*fz*pW(ic,jc+1,iz+1,wc)+wx*wy*fz*pW(ic+1,jc+1,iz+1,wc)
         else
            fvel(1)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ic,jc,kc,uc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ic+1,jc,kc,uc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pU(ic,jc+1,kc,uc)+wx*wy*(1.0_WP-wz)*pU(ic+1,jc+1,kc,uc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pU(ic,jc,kc+1,uc)+wx*(1.0_WP-wy)*wz*pU(ic+1,jc,kc+1,uc)+(1.0_WP-wx)*wy*wz*pU(ic,jc+1,kc+1,uc)+wx*wy*wz*pU(ic+1,jc+1,kc+1,uc)
            fvel(2)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pV(ic,jc,kc,vc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pV(ic+1,jc,kc,vc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pV(ic,jc+1,kc,vc)+wx*wy*(1.0_WP-wz)*pV(ic+1,jc+1,kc,vc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pV(ic,jc,kc+1,vc)+wx*(1.0_WP-wy)*wz*pV(ic+1,jc,kc+1,vc)+(1.0_WP-wx)*wy*wz*pV(ic,jc+1,kc+1,vc)+wx*wy*wz*pV(ic+1,jc+1,kc+1,vc)
            fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pW(ic,jc,kc,wc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pW(ic+1,jc,kc,wc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pW(ic,jc+1,kc,wc)+wx*wy*(1.0_WP-wz)*pW(ic+1,jc+1,kc,wc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pW(ic,jc,kc+1,wc)+wx*(1.0_WP-wy)*wz*pW(ic+1,jc,kc+1,wc)+(1.0_WP-wx)*wy*wz*pW(ic,jc+1,kc+1,wc)+wx*wy*wz*pW(ic+1,jc+1,kc+1,wc)
         end if
         ! Density and viscosity
         frho =(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pRho (ic,jc,kc,1)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pRho (ic+1,jc,kc,1)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pRho (ic,jc+1,kc,1)+wx*wy*(1.0_WP-wz)*pRho (ic+1,jc+1,kc,1)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pRho (ic,jc,kc+1,1)+wx*(1.0_WP-wy)*wz*pRho (ic+1,jc,kc+1,1)+(1.0_WP-wx)*wy*wz*pRho (ic,jc+1,kc+1,1)+wx*wy*wz*pRho (ic+1,jc+1,kc+1,1)
         fvisc=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pVisc(ic,jc,kc,1)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pVisc(ic+1,jc,kc,1)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pVisc(ic,jc+1,kc,1)+wx*wy*(1.0_WP-wz)*pVisc(ic+1,jc+1,kc,1)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pVisc(ic,jc,kc+1,1)+wx*(1.0_WP-wy)*wz*pVisc(ic+1,jc,kc+1,1)+(1.0_WP-wx)*wy*wz*pVisc(ic,jc+1,kc+1,1)+wx*wy*wz*pVisc(ic+1,jc+1,kc+1,1)
         fvisc=fvisc+epsilon(1.0_WP)
         ! Volume fraction
         pVF=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pVolFrac(ic,jc,kc,1)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pVolFrac(ic+1,jc,kc,1)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pVolFrac(ic,jc+1,kc,1)+wx*wy*(1.0_WP-wz)*pVolFrac(ic+1,jc+1,kc,1)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pVolFrac(ic,jc,kc+1,1)+wx*(1.0_WP-wy)*wz*pVolFrac(ic+1,jc,kc+1,1)+(1.0_WP-wx)*wy*wz*pVolFrac(ic,jc+1,kc+1,1)+wx*wy*wz*pVolFrac(ic+1,jc+1,kc+1,1)
         fVF=1.0_WP-pVF
         ! Drag correction factor
         select case(trim(this%drag_model))
         case('None','none')
            corr=epsilon(1.0_WP)
         case('Stokes')
            corr=1.0_WP
         case('Schiller-Naumann','Schiller Naumann','SN')
            Re=frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            corr=1.0_WP+0.15_WP*Re**(0.687_WP)
         case('Tenneti') ! Tenneti and Subramaniam (2011)
            Re=fVF*frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            b1=5.81_WP*pVF/fVF**3+0.48_WP*pVF**(1.0_WP/3.0_WP)/fVF**4
            b2=pVF**3*Re*(0.95_WP+0.61_WP*pVF**3/fVF**2)
            corr=fVF*((1.0_WP+0.15_WP*Re**(0.687_WP))/fVF**3+b1+b2)
         case('Beetstra') ! Beetstra et al. (2007)
            Re=fVF*frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            b1=10.0_WP*pVF/fVF**2+fVF**2*(1.0_WP+1.5_WP*sqrt(pVF))
            b2=0.413_WP/24.0_WP*Re/fVF**2*(1.0_WP/fVF+3.0_WP*fVF*pVF+8.4_WP*Re**(-0.343_WP))/(1.0_WP+10.0_WP**(3.0_WP*pVF)*Re**(2.0_WP*fVF-2.5_WP))
            corr=b1+b2
         case default
            corr=1.0_WP
         end select
         ! Particle response time
         tau=this%rho*myp%d**2/(18.0_WP*fvisc*corr)
         ! Return acceleration and update particle timestep size
         acc=(fvel-myp%vel)/tau
         myp%dt=tau/real(this%nstep,WP)
      end function get_rhs

      !> Deposit a particle source term back to the mesh
      subroutine deposit(val)
         real(WP), dimension(3), intent(in) :: val
         real(WP) :: wx,wy,wz
         integer  :: ic,jc,kc
         ! Cell-centered stencil indices and weights
         ic=floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(myp%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
         jc=floor((myp%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(myp%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
         kc=floor((myp%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(myp%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
         pSrc(ic:ic+1,jc:jc+1,kc:kc+1,1)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,1)+reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])*val(1)
         pSrc(ic:ic+1,jc:jc+1,kc:kc+1,2)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,2)+reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])*val(2)
         pSrc(ic:ic+1,jc:jc+1,kc:kc+1,3)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,3)+reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])*val(3)
      end subroutine deposit

   end subroutine advance


   !> Update particle volume fraction field based on our current particles
   subroutine update_VF(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer(I8) :: np_
      integer :: lvl,i,ii,jj,kk
      real(WP) :: dxi,dyi,dzi,Vp,wx,wy,wz
      ! Zero VF on all levels
      call this%VF%setval(0.0_WP)
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         ! Get inverse of mesh size
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            ! Get pointer to data
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            ! Loop over particles
            call this%get_particles(lvl=lvl,mfi=mfi,p=p,np=np_)
            do i=1,np_
               ! Skip particles that don't exchange with the fluid
               if (IAND(p(i)%flag,PART_EXCHANGES).eq.0) cycle
               ! Get particle volume
               Vp=Pi/6.0_WP*p(i)%d**3
               ! Extrapolate
               ii=floor((p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ii,WP)
               jj=floor((p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jj,WP)
               kk=floor((p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kk,WP)
               pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)=pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)+Vp*reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])
            end do
         end do
         call amrex_mfiter_destroy(mfi)
         ! Divide by cell volume
         call this%VF%mf(lvl)%mult(1.0_WP/this%amr%cell_vol(lvl),0,1,this%nover)
      end do
      ! Sum overlap data across boxes and levels, sync
      call this%VF%syncsum(); call this%VF%sum_down(); call this%VF%average_down(); call this%VF%fill(time=0.0_WP)
      ! Filter
      call this%filter(this%VF)
   end subroutine update_VF

   !> CFL based on particle velocities relative to their local cell size
   subroutine get_cfl(this,dt,cflc,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional    :: cfl
      ! NOTE: without per-particle dx we can only report |vel|*dt.
      ! A proper implementation queries this%amr%dx(lev) at the particle cell.
      this%CFLp_x = 0.0_WP; this%CFLp_y = 0.0_WP; this%CFLp_z = 0.0_WP
      this%CFL_col = 0.0_WP
      cflc = 0.0_WP
      if (present(cfl)) cfl = 0.0_WP
      ! TODO: loop over tiles, query geom%dx at particle position
   end subroutine get_cfl

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Redistribute particles
   subroutine redistribute(this,minlvl,maxlvl,nover)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in), optional :: minlvl,maxlvl,nover
      integer :: lmin,lmax,no
      lmin= 0; if (present(minlvl)) lmin=minlvl
      lmax=-1; if (present(maxlvl)) lmax=maxlvl
      no  = 0; if (present(nover))  no  =nover
      call amrlpt_redistribute(this%pc,lmin,lmax,no)
      call amrlpt_total_np(this%pc,this%np)
   end subroutine redistribute

   !> Fill ghost particle buffer within no cells
   subroutine fill_ghosts(this,no)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in), optional :: no
      integer :: ng
      ng=this%nover; if (present(no)) ng=no
      call amrlpt_fill_neighbors(this%pc,ng)
   end subroutine fill_ghosts

   !> Release ghost particle buffer
   subroutine clear_ghosts(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call amrlpt_clear_neighbors(this%pc)
   end subroutine clear_ghosts

   !> Build explicit pair list within interaction radius rcrit
   subroutine build_neighbor_list(this,rcrit)
      use iso_c_binding, only: c_double
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: rcrit
      call amrlpt_build_neighbor_list(this%pc,real(rcrit,c_double))
   end subroutine build_neighbor_list

   !> Update global particle count
   subroutine get_np(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call amrlpt_total_np(this%pc,this%np)
   end subroutine get_np

   !> Get pointer to particle array for a given tile
   subroutine get_particles(this,lvl,mfi,p,np)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_c
      call amrlpt_get_particles_mfi(this%pc,lvl,mfi%p,dp,np_c)
      np=np_c
      if (np.gt.0) then
         call c_f_pointer(dp,p,[np])
      else
         nullify(p)
      end if
   end subroutine get_particles

   !> Interpolate from cell-centered data to a point
   function interp(this,lvl,pos,arr,comp) result(val)
      implicit none
      class(amrlpt), intent(in) :: this
      integer,  intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: arr
      integer,  intent(in) :: comp
      real(WP) :: val,wx,wy,wz
      integer  :: ii,jj,kk
      ! Get indices and weights
      ii=floor((pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP); wx=(pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP-real(ii,WP)
      jj=floor((pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP); wy=(pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP-real(jj,WP)
      kk=floor((pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP); wz=(pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP-real(kk,WP)
      ! Trilinear interpolation
      val=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*arr(ii  ,jj  ,kk  ,comp) &
      &  +        wx *(1.0_WP-wy)*(1.0_WP-wz)*arr(ii+1,jj  ,kk  ,comp) &
      &  +(1.0_WP-wx)*        wy *(1.0_WP-wz)*arr(ii  ,jj+1,kk  ,comp) &
      &  +        wx *        wy *(1.0_WP-wz)*arr(ii+1,jj+1,kk  ,comp) &
      &  +(1.0_WP-wx)*(1.0_WP-wy)*        wz *arr(ii  ,jj  ,kk+1,comp) &
      &  +        wx *(1.0_WP-wy)*        wz *arr(ii+1,jj  ,kk+1,comp) &
      &  +(1.0_WP-wx)*        wy *        wz *arr(ii  ,jj+1,kk+1,comp) &
      &  +        wx *        wy *        wz *arr(ii+1,jj+1,kk+1,comp)
   end function interp

   !> Interpolate face velocities to a point
   function interp_face_velocities(this,lvl,pos,pU,uc,pV,vc,pW,wc) result(fvel)
      implicit none
      class(amrlpt), intent(in) :: this
      integer,  intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: pU,pV,pW
      integer,  intent(in) :: uc,vc,wc
      real(WP), dimension(3) :: fvel
      real(WP) :: dxi,dyi,dzi,wx,wy,wz,fx,fy,fz
      integer  :: ic,jc,kc,ix,iy,iz
      ! Get cell-centered indices and weights
      dxi=1.0_WP/this%amr%dx(lvl); ic=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
      dyi=1.0_WP/this%amr%dy(lvl); jc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
      dzi=1.0_WP/this%amr%dz(lvl); kc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
      ! Get face-centered indices and weights
      ix=floor((pos(1)-this%amr%xlo)*dxi); fx=(pos(1)-this%amr%xlo)*dxi-real(ix,WP)
      iy=floor((pos(2)-this%amr%ylo)*dyi); fy=(pos(2)-this%amr%ylo)*dyi-real(iy,WP)
      iz=floor((pos(3)-this%amr%zlo)*dzi); fz=(pos(3)-this%amr%zlo)*dzi-real(iz,WP)
      ! U: staggered in x — cell indices jc/wy, kc/wz; face index ix/fx
      fvel(1)=(1.0_WP-fx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix  ,jc  ,kc  ,uc)+fx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix+1,jc  ,kc  ,uc) &
      &      +(1.0_WP-fx)*        wy *(1.0_WP-wz)*pU(ix  ,jc+1,kc  ,uc)+fx*        wy *(1.0_WP-wz)*pU(ix+1,jc+1,kc  ,uc) &
      &      +(1.0_WP-fx)*(1.0_WP-wy)*        wz *pU(ix  ,jc  ,kc+1,uc)+fx*(1.0_WP-wy)*        wz *pU(ix+1,jc  ,kc+1,uc) &
      &      +(1.0_WP-fx)*        wy *        wz *pU(ix  ,jc+1,kc+1,uc)+fx*        wy *        wz *pU(ix+1,jc+1,kc+1,uc)
      ! V: staggered in y — cell indices ic/wx, kc/wz; face index iy/fy
      fvel(2)=(1.0_WP-wx)*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic  ,iy  ,kc  ,vc)+wx*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic+1,iy  ,kc  ,vc) &
      &      +(1.0_WP-wx)*        fy *(1.0_WP-wz)*pV(ic  ,iy+1,kc  ,vc)+wx*        fy *(1.0_WP-wz)*pV(ic+1,iy+1,kc  ,vc) &
      &      +(1.0_WP-wx)*(1.0_WP-fy)*        wz *pV(ic  ,iy  ,kc+1,vc)+wx*(1.0_WP-fy)*        wz *pV(ic+1,iy  ,kc+1,vc) &
      &      +(1.0_WP-wx)*        fy *        wz *pV(ic  ,iy+1,kc+1,vc)+wx*        fy *        wz *pV(ic+1,iy+1,kc+1,vc)
      ! W: staggered in z — cell indices ic/wx, jc/wy; face index iz/fz
      fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic  ,jc  ,iz  ,wc)+wx*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic+1,jc  ,iz  ,wc) &
      &      +(1.0_WP-wx)*        wy *(1.0_WP-fz)*pW(ic  ,jc+1,iz  ,wc)+wx*        wy *(1.0_WP-fz)*pW(ic+1,jc+1,iz  ,wc) &
      &      +(1.0_WP-wx)*(1.0_WP-wy)*        fz *pW(ic  ,jc  ,iz+1,wc)+wx*(1.0_WP-wy)*        fz *pW(ic+1,jc  ,iz+1,wc) &
      &      +(1.0_WP-wx)*        wy *        fz *pW(ic  ,jc+1,iz+1,wc)+wx*        wy *        fz *pW(ic+1,jc+1,iz+1,wc)
   end function interp_face_velocities

   !> Explicit diffusion filter for a cell-centered amrdata field
   subroutine filter(this,A)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_multifab,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      real(WP) :: alpha,alpha_step
      integer  :: nstep,n,nc,lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pA,pFx,pFy,pFz
      real(WP) :: dxi,dyi,dzi

      ! Determine filter coefficient from finest level mesh size
      alpha=max(this%filter_width**2-this%amr%min_meshsize(this%amr%clvl())**2,0.0_WP)/(16.0_WP*log(2.0_WP))
      if (alpha.le.0.0_WP) return

      ! Number of explicit sub-steps for stability (dt < dx^2/6 in 3D)
      nstep=ceiling(6.0_WP*alpha/this%amr%min_meshsize(this%amr%clvl())**2)
      alpha_step=alpha/real(nstep,WP)

      ! Allocate face flux mfabs (no ghost needed for fluxes)
      allocate(Fx(0:this%amr%maxlvl),Fy(0:this%amr%maxlvl),Fz(0:this%amr%maxlvl))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=A%ncomp,nover=0,atface=[.true., .false.,.false.]); call Fx(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.true., .false.]); call Fy(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
      end do

      ! Explicit sub-steps
      do n=1,nstep
         ! Compute diffusive fluxes at each level
         do lvl=0,this%amr%clvl()
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pA =>A%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(1)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFx(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i-1,j,k,nc))*dxi
               end do; end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFy(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i,j-1,k,nc))*dyi
               end do; end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFz(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i,j,k-1,nc))*dzi
               end do; end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
         ! Average down face fluxes for C/F conservation (finest→coarsest)
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         end do
         ! Apply divergence to update A
         do lvl=0,this%amr%clvl()
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pA =>A%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%tilebox()
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pA(i,j,k,nc)=pA(i,j,k,nc)+dxi*(pFx(i+1,j,k,nc)-pFx(i,j,k,nc))+dyi*(pFy(i,j+1,k,nc)-pFy(i,j,k,nc))+dzi*(pFz(i,j,k+1,nc)-pFz(i,j,k,nc))
               end do; end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
         ! Restore coarse valid cells covered by fine
         call A%average_down()
         ! Fill ghosts at all levels
         call A%fill(time=0.0_WP)
      end do

      ! Destroy flux mfabs
      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(Fx(lvl))
         call amrex_multifab_destroy(Fy(lvl))
         call amrex_multifab_destroy(Fz(lvl))
      end do
      deallocate(Fx,Fy,Fz)

   end subroutine filter

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Compute particle statistics: np, d/vel min/max/mean/var.
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(c_ptr) :: dp
      type(part), pointer :: p(:)
      integer(c_int64_t) :: np_tile
      integer :: lev,n,ierr
      real(WP) :: d_sum,d_sq,vx_sum,vx_sq,vy_sum,vy_sq,vz_sum,vz_sq
      real(WP) :: inv_np
      ! Init per-rank accumulators
      this%np=0
      this%dmin=huge(1.0_WP); this%dmax=-huge(1.0_WP);  d_sum=0.0_WP;  d_sq=0.0_WP
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); vx_sum=0.0_WP; vx_sq=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); vy_sum=0.0_WP; vy_sq=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); vz_sum=0.0_WP; vz_sq=0.0_WP
      ! Loop over all AMR levels and tiles
      do lev=0,this%amr%clvl()
         call this%amr%mfiter_build(lev,mfi)
         do while (mfi%valid())
            call amrlpt_get_particles_mfi(this%pc,lev,mfi%p,dp,np_tile)
            if (np_tile.gt.0) then
               call c_f_pointer(dp,p,[np_tile])
               do n=1,int(np_tile)
                  if (p(n)%flag.eq.1) cycle
                  this%np=this%np+1
                  this%dmin=min(this%dmin,p(n)%d);      this%dmax=max(this%dmax,p(n)%d);       d_sum= d_sum+p(n)%d;       d_sq= d_sq+p(n)%d**2
                  this%Umin=min(this%Umin,p(n)%vel(1)); this%Umax=max(this%Umax,p(n)%vel(1)); vx_sum=vx_sum+p(n)%vel(1); vx_sq=vx_sq+p(n)%vel(1)**2
                  this%Vmin=min(this%Vmin,p(n)%vel(2)); this%Vmax=max(this%Vmax,p(n)%vel(2)); vy_sum=vy_sum+p(n)%vel(2); vy_sq=vy_sq+p(n)%vel(2)**2
                  this%Wmin=min(this%Wmin,p(n)%vel(3)); this%Wmax=max(this%Wmax,p(n)%vel(3)); vz_sum=vz_sum+p(n)%vel(3); vz_sq=vz_sq+p(n)%vel(3)**2
               end do
            end if
            call mfi%next()
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Global MPI reduce
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np,  1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sum,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sq,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sum,   1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sq,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sum,   1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sq,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sum,   1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sq,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      ! Derive mean and variance
      if (this%np.gt.0) then
         inv_np=1.0_WP/real(this%np,WP)
         this%dmean= d_sum*inv_np; this%dvar=max(0.0_WP, d_sq*inv_np-this%dmean**2)
         this%Umean=vx_sum*inv_np; this%Uvar=max(0.0_WP,vx_sq*inv_np-this%Umean**2)
         this%Vmean=vy_sum*inv_np; this%Vvar=max(0.0_WP,vy_sq*inv_np-this%Vmean**2)
         this%Wmean=vz_sum*inv_np; this%Wvar=max(0.0_WP,vz_sq*inv_np-this%Wmean**2)
      end if
   end subroutine get_info

   !> Print solver info
   subroutine print(this)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(amrlpt), intent(in) :: this
      character(len=str_long) :: message
      call log('AMR Lagrangian particle solver: '//trim(this%name))
      write(message,'("  Particle density  : ",ES12.5)') this%rho
      call log(trim(message))
      write(message,'("  Drag model        : ",a)') trim(this%drag_model)
      call log(trim(message))
      write(message,'("  ODE substeps      : ",i0)') this%nstep
      call log(trim(message))
      call log('  Grid: '//trim(this%amr%name))
   end subroutine print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Write AMReX checkpoint for particles to dirname
   subroutine write(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_write(this%pc,trim(dirname)//c_null_char,1_c_int)
   end subroutine write

   !> Read AMReX checkpoint for particles from dirname
   subroutine read(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_read(this%pc,trim(dirname)//c_null_char)
      call amrlpt_redistribute(this%pc,0,-1,0)
      call amrlpt_total_np(this%pc,this%np)
      call this%update_VF()
   end subroutine read

end module amrlpt_class
