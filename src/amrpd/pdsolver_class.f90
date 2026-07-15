!> Peridynamics solver: node-centered, CSR-based solid dynamics on flat
!> per-field arrays with persistent graph-halo communication. GRID-FREE:
!> no AMReX anywhere -- ownership follows the reference configuration
!> (Morton partition, motion-invariant), neighborhoods and communication
!> plans are built once and reused every substep.
!>
!> Physics: linear peridynamic solid (LPS, dimension-aware coefficients),
!> brittle stretch damage, per-side viscoelastic/viscoplastic flow with J2
!> (Mitchell OSB) yield, soft-sphere contact (walls + particle-particle via
!> a displacement-triggered spatial service), velocity-Verlet integration.
!> Checkpoint/restart is gid-space and rank-count portable, including all
!> bond damage and plastic history.
!>
!> Usage tiers (see amrpd_class for the grid-side face):
!>   1. pdsolver alone      -- standalone solid dynamics (this module only)
!>   2. pdsolver + amrpd    -- adds viz, mesh VF, AMR refinement, seeding
!>   3. ... + a flow solver -- two-way FSI via pdsolver%exchange
!>
!> Layout: owned nodes 1..nown; halo slots nown+1..ntot, keyed (gid, periodic
!> image offset) with shifts applied at exchange time. Each physical bond is
!> two CSR half-entries (one per endpoint row, Peridigm convention): kernels
!> compute each row's own force state -- ghost dilatation is never
!> communicated -- and a single halo reduce assembles cross-rank pairs.
!> Kernels are pure loops over owned nodes with no mutable module-level
!> state (OpenMP-ready by construction; threads deferred).
module pdsolver_class
   use precision,        only: WP,I8
   use string,           only: str_medium
   use pdhalo_class,     only: pddir,pdhalo,sort3_perm,PDHALO_KEY0
   use pdhash_class, only: gid_hash
   implicit none
   private

   public :: pdsolver,pd_partition
   public :: PDC_IS_DEAD,PDC_MOVES,PDC_INTEGRATES,PDC_BONDS

   ! Motion-control bit flags -- values MUST match amrpd's PART_* constants
   ! (handoff copies amrpd flags verbatim)
   integer, parameter :: PDC_IS_DEAD   =0
   integer, parameter :: PDC_MOVES     =1
   integer, parameter :: PDC_INTEGRATES=2
   integer, parameter :: PDC_BONDS     =4

   !> Graph-core PD solver
   type :: pdsolver
      character(len=str_medium) :: name='UNNAMED_PDSOLVER'

      ! Sizes
      integer :: nown=0                        !< owned nodes on this rank
      integer :: nhalo=0                       !< halo slots
      integer :: ntot=0                        !< nown+nhalo
      integer(I8) :: np=0                      !< global node count (get_info)
      integer(I8) :: nbond=0                   !< global bond count (half-entries/2, get_info)

      ! Material / discretization
      real(WP) :: rho            =0.0_WP       !< density
      real(WP) :: elastic_modulus=0.0_WP       !< Young's modulus
      real(WP) :: poisson_ratio  =0.0_WP       !< Poisson's ratio
      real(WP) :: delta          =0.0_WP       !< horizon
      real(WP) :: dV             =0.0_WP       !< nominal element volume (CFL length scale; kernels use per-node V)
      real(WP) :: s0             =huge(1.0_WP) !< critical bond stretch (huge = no damage)
      real(WP) :: dtcrit         =0.0_WP       !< Silling-Askari critical dt (diagnostic, stamped at connect)
      ! Viscoelastic / viscoplastic flow (PER-SIDE form: each half-entry evolves
      ! its own e_v with its own endpoint's dilatation and yield factor --
      ! exactly Peridigm's elastic_plastic.cxx, verified term-by-term against
      ! amrpd's J2 2026-07-14. This is the one INTENDED delta from amrpd, which
      ! averages the endpoints into a single per-bond e_v.)
      real(WP) :: tau            =huge(1.0_WP) !< Maxwell deviatoric relaxation time (huge = purely elastic)
      real(WP) :: visc_lambda    =1.0_WP       !< SLS relaxing fraction [0,1]
      real(WP) :: yield_stretch  =0.0_WP       !< legacy per-bond Perzyna yield strain (0 = pure Maxwell)
      real(WP) :: sigma_yield    =0.0_WP       !< J2 yield stress (Mitchell OSB family norm; overrides yield_stretch)
      real(WP), dimension(3) :: gravity=0.0_WP !< body acceleration
      logical,  dimension(3) :: collapsed=.false. !< collapsed (n==1) directions: velocity locked
      real(WP), dimension(3) :: Ldom=0.0_WP    !< domain lengths (image shifts)
      logical,  dimension(3) :: per=.false.    !< periodicity per direction
      real(WP), dimension(3) :: dom_lo=0.0_WP  !< domain lower bounds (wall contact)
      real(WP), dimension(3) :: dom_hi=0.0_WP  !< domain upper bounds (wall contact)

      ! Short-range soft-sphere contact (soft-sphere penalty + damping).
      ! Contact is a pure SPATIAL service, fully separate from the bond graph:
      ! candidates are (owned nodes + contact-halo slots) discovered by a
      ! displacement-triggered broad phase; the graph halo is never binned, so
      ! bonded remote partners arrive as contact slots when in range and
      ! double-counting is structurally impossible. The narrow phase is
      ! gather-only (each owned node accumulates from its candidates -- the
      ! partner gets its share from its own row), so no force reduction.
      logical  :: use_contact  =.false.
      real(WP) :: contact_dist =0.0_WP         !< d_c (p-p d_eff; wall d_eff = 0.5*d_c)
      real(WP) :: tau_col      =0.0_WP         !< collision duration (<=0 -> auto 5*dt)
      real(WP) :: e_n=0.7_WP,e_w=0.7_WP        !< restitution (p-p, wall)
      real(WP) :: clip_col     =0.2_WP         !< overlap clip fraction
      integer, dimension(3) :: lo_bc=0,hi_bc=0 !< per-face: 0=open, 1=wall (matches AMRPD_* values)
      real(WP) :: cskin        =0.0_WP         !< broad-phase skin (<=0 -> auto 0.5*contact_dist)
      type(pdhalo) :: chalo                    !< contact halo (rebuilt at trigger cadence; nown=ntot)
      integer :: nchalo=0                      !< contact slots (y/v extended to ntot+nchalo)
      integer, allocatable :: cptr(:),clst(:)  !< candidate CSR (owned rows; entries index owned+contact slots)
      real(WP), allocatable :: ylast(:,:)      !< (3,nown) positions at last broad-phase build

      ! Node state -- flat per-field arrays; owned first, halo slots appended.
      ! (3,:) fields are xyz-interleaved per node (Fortran-natural gather layout).
      integer(I8), allocatable :: gid(:)       !< (ntot) global id (halo slots carry partner gid)
      real(WP), allocatable :: x0(:,:)         !< (3,ntot) reference position; halo PRE-SHIFTED
      real(WP), allocatable :: y(:,:)          !< (3,ntot) current position; halo shifted at update
      real(WP), allocatable :: v(:,:)          !< (3,nown) velocity
      real(WP), allocatable :: f(:,:)          !< (3,ntot) bond force density (halo = scatter buffer)
      real(WP), allocatable :: ff(:,:)         !< (3,nown) external (fluid) force density
      real(WP), allocatable :: vol(:)          !< (ntot) per-node volume (reference; halo slots filled at connect)
      real(WP), allocatable :: mw(:)           !< (nown) weighted volume (reference, set at connect)
      real(WP), allocatable :: theta(:)        !< (nown) dilatation (recomputed each substep)
      real(WP), allocatable :: damage(:)       !< (nown) accumulated damage fraction (broken/reference bonds)
      real(WP), allocatable :: alive(:)        !< (ntot) 1=alive, 0=dead (exit through open face); halo-exchanged on death events only
      integer,  allocatable :: flag(:)         !< (nown) motion-control flags
      logical :: watch_exit=.false.            !< exit detection active (set at connect: domain set + any open non-periodic face)

      ! CSR families (built once at connect)
      integer, allocatable :: ptr(:)           !< (nown+1) row offsets
      integer, allocatable :: lst(:)           !< neighbor index (1..ntot) per half-entry
      integer(1), allocatable :: dmg(:)        !< per half-entry: 0 intact, 1 broken (irreversible)
      real(WP), allocatable :: e_v(:)          !< per half-entry: inelastic deviatoric stretch (per-side history)
      real(WP), allocatable :: td2(:),td2a(:)  !< (nown) J2 family deviatoric norm^2: previous substep / accumulator.
                                               !< Node-centered => pure own-row gather, NO communication (amrpd
                                               !< needed sum_ghosts_td2 + ghost refresh for the same quantity).

      ! Parallel machinery
      type(gid_hash) :: ohash                  !< gid -> owned index (built at set_nodes)
      type(pddir)    :: dir                    !< persistent gid directory (owner lookups; registered at connect/read_state)
      type(pdhalo)   :: halo                   !< persistent halo plan
      real(WP), allocatable :: rextra_tmp(:,:) !< read_state scratch (restart-field overlay across assemble)

      ! Monitoring
      real(WP) :: Umax=0.0_WP                  !< max |velocity component| (get_info)
      real(WP) :: CFLe=0.0_WP,CFLp=0.0_WP      !< elastic-wave / convective CFL (get_cfl)
      integer(I8) :: nbroken=0                 !< global broken half-entry count (internal)
      integer(I8) :: nb_broken=0               !< global broken BOND count (exact census, get_info)
      integer(I8) :: nb=0                      !< global bond count (exact census, stamped at assemble)
      integer(I8) :: nrebuild=0                !< broad-phase rebuild count (cumulative)
      integer(I8) :: nchalo_glob=0             !< global contact-slot count (get_info)
      integer(I8) :: ncand_glob=0              !< global contact-candidate count (get_info)

      ! Per-rank phase timers (accumulated in advance; reduced+reset in get_info)
      real(WP) :: wt_kick=0.0_WP,wt_halo=0.0_WP,wt_dil=0.0_WP,wt_force=0.0_WP,wt_reduce=0.0_WP
      real(WP) :: wt_contact=0.0_WP,wt_broad=0.0_WP
      real(WP) :: wtmax_kick=0.0_WP,wtmax_halo=0.0_WP,wtmax_dil=0.0_WP,wtmax_force=0.0_WP,wtmax_reduce=0.0_WP
      real(WP) :: wtmax_contact=0.0_WP,wtmax_broad=0.0_WP
      real(WP) :: wtmin_dil=0.0_WP,wtmin_force=0.0_WP

   contains
      procedure :: initialize
      procedure :: set_nodes
      procedure :: connect
      procedure :: detect_families
      procedure :: advance
      procedure :: exchange
      procedure :: query_owners
      procedure :: write_state
      procedure :: read_state
      procedure :: get_cfl
      procedure :: get_info
      procedure :: finalize
      procedure, private :: lps_coefs
      procedure, private :: compute_mw
      procedure, private :: contact_broadphase
      procedure, private :: contact_narrow
      procedure, private :: assemble
   end type pdsolver

contains


   !> Configure the solver (no allocation yet; set_nodes sizes the state)
   subroutine initialize(this,name,rho,elastic_modulus,poisson_ratio,delta,dV,gravity,collapsed,Ldom,per)
      implicit none
      class(pdsolver), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(in) :: rho,elastic_modulus,poisson_ratio,delta,dV
      real(WP), dimension(3), intent(in) :: gravity,Ldom
      logical,  dimension(3), intent(in) :: collapsed,per
      this%name=trim(adjustl(name))
      this%rho=rho
      this%elastic_modulus=elastic_modulus
      this%poisson_ratio=poisson_ratio
      this%delta=delta
      this%dV=dV
      this%gravity=gravity
      this%collapsed=collapsed
      this%Ldom=Ldom
      this%per=per
   end subroutine initialize

   !> Load this rank's owned nodes (any distribution; it becomes the static
   !> partition). Builds the gid->index hash used by connect and the halo plan.
   !> vol is the per-node volume (pass a constant-filled array for a uniform
   !> lattice; kernels use it per neighbor, Peridigm-style).
   subroutine set_nodes(this,n,gids,pos,vel,flags,vol)
      implicit none
      class(pdsolver), intent(inout) :: this
      integer, intent(in) :: n
      integer(I8), intent(in) :: gids(:)
      real(WP), intent(in) :: pos(:,:),vel(:,:)
      integer, intent(in) :: flags(:)
      real(WP), intent(in) :: vol(:)
      integer :: i
      this%nown=n
      this%nhalo=0
      this%ntot=n
      allocate(this%gid(max(n,1)),this%x0(3,max(n,1)),this%y(3,max(n,1)))
      allocate(this%v(3,max(n,1)),this%f(3,max(n,1)),this%ff(3,max(n,1)))
      allocate(this%mw(max(n,1)),this%theta(max(n,1)),this%flag(max(n,1)))
      allocate(this%vol(max(n,1)),this%damage(max(n,1)))
      do i=1,n
         this%gid(i) =gids(i)
         this%x0(:,i)=pos(:,i)
         this%y(:,i) =pos(:,i)
         this%v(:,i) =vel(:,i)
         this%flag(i)=flags(i)
         this%vol(i)   =vol(i)
      end do
      this%f=0.0_WP; this%ff=0.0_WP; this%mw=0.0_WP; this%theta=0.0_WP; this%damage=0.0_WP
      call this%ohash%build(n,gids(1:n))
   end subroutine set_nodes

   !> Build the static CSR families and the halo plan from a distributed bond
   !> list (this rank passes the bonds it holds; any distribution is fine --
   !> half-entries are routed to their node's owner through the gid directory).
   !> Collective. bkey packs the periodic image offset of the HI endpoint in
   !> amrpd's hist1 convention. Self-image bonds (gid_lo==gid_hi) yield ONE
   !> half-entry (the opposite-image bond exists separately in the input, exactly
   !> as amrpd stores them).
   subroutine connect(this,nb,bgid_lo,bgid_hi,bkey)
      use parallel, only: comm,nproc
      use messager, only: die
      use mpi_f08
      implicit none
      class(pdsolver), intent(inout) :: this
      integer, intent(in) :: nb
      integer(I8), intent(in) :: bgid_lo(:),bgid_hi(:)
      integer, intent(in) :: bkey(:)
      integer(I8), allocatable :: hnode(:),hnbr(:),rnode(:),rnbr(:)
      integer, allocatable :: hkey(:),howner(:),rkey(:)
      real(WP), allocatable :: rev(:)
      integer(1), allocatable :: rdmg(:)
      integer :: nhe,rn,i,ib,ierr

      ! Distributed gid directory over the node partition (persistent: also
      ! serves owner queries for face-tag restamping after restart)
      call this%dir%finalize()
      call this%dir%register(this%nown,this%gid(1:this%nown))

      ! Expand bonds into half-entries (one per endpoint row; one total for
      ! self-image bonds -- see header)
      nhe=0
      do ib=1,nb
         nhe=nhe+1
         if (bgid_lo(ib).ne.bgid_hi(ib)) nhe=nhe+1
      end do
      allocate(hnode(max(nhe,1)),hnbr(max(nhe,1)),hkey(max(nhe,1)),howner(max(nhe,1)))
      nhe=0
      do ib=1,nb
         nhe=nhe+1
         hnode(nhe)=bgid_lo(ib); hnbr(nhe)=bgid_hi(ib); hkey(nhe)=bkey(ib)
         if (bgid_lo(ib).ne.bgid_hi(ib)) then
            nhe=nhe+1
            hnode(nhe)=bgid_hi(ib); hnbr(nhe)=bgid_lo(ib); hkey(nhe)=negkey(bkey(ib))
         end if
      end do

      ! Route each half-entry to the rank owning its node
      call this%dir%query(nhe,hnode,howner)
      route_entries: block
         integer, dimension(0:nproc-1) :: sc,rc,sd,rd
         integer, allocatable :: pos(:)
         integer(I8), allocatable :: s8(:)
         integer, allocatable :: s4(:)
         integer :: r,h
         sc=0
         do i=1,nhe
            sc(howner(i))=sc(howner(i))+1
         end do
         call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
         sd(0)=0; rd(0)=0
         do r=1,nproc-1
            sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
         end do
         rn=sum(rc)
         allocate(rnode(max(rn,1)),rnbr(max(rn,1)),rkey(max(rn,1)))
         allocate(pos(0:nproc-1),s8(max(nhe,1)),s4(max(nhe,1)))
         ! node gids
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s8(pos(h))=hnode(i)
         end do
         call MPI_ALLTOALLV(s8,sc,sd,MPI_INTEGER8,rnode,rc,rd,MPI_INTEGER8,comm,ierr)
         ! neighbor gids
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s8(pos(h))=hnbr(i)
         end do
         call MPI_ALLTOALLV(s8,sc,sd,MPI_INTEGER8,rnbr,rc,rd,MPI_INTEGER8,comm,ierr)
         ! image keys
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s4(pos(h))=hkey(i)
         end do
         call MPI_ALLTOALLV(s4,sc,sd,MPI_INTEGER,rkey,rc,rd,MPI_INTEGER,comm,ierr)
         deallocate(pos,s8,s4)
      end block route_entries
      deallocate(hnode,hnbr,hkey,howner)

      deallocate(hnode,hnbr,hkey,howner)

      ! Fresh bonds carry zero inelastic state
      allocate(rev(max(rn,1)),rdmg(max(rn,1)))
      rev=0.0_WP; rdmg=0_1
      call this%assemble(rn,rnode,rnbr,rkey,rev,rdmg)
      deallocate(rnode,rnbr,rkey,rev,rdmg)
   end subroutine connect


   !> Build the bond families directly from the REFERENCE configuration -- no
   !> amrpd bond container, no bond expansion: distributed neighbor discovery
   !> at radius delta (bounds allgather + per-(rank,image-offset) offers of
   !> shifted x0, contact-broadphase pattern), then each owned row's
   !> half-entries are generated straight from the binned candidates and fed
   !> to assemble with zero inelastic state. Acceptance test r2 <= delta^2
   !> matches amrpd bond_init exactly. Collective; call after set_nodes.
   subroutine detect_families(this)
      use parallel, only: comm,rank,nproc,amRoot,MPI_REAL_WP
      use messager, only: log,die
      use string,   only: str_long
      use mpi_f08
      implicit none
      class(pdsolver), intent(inout) :: this
      real(WP), dimension(3) :: bl,bh,shift,pos_s,gl,hcell
      real(WP), allocatable :: allb(:,:),opos(:,:),rpos(:,:),cpos(:,:)
      integer(I8), allocatable :: ogid(:),rgid(:),cgid(:),hnode(:),hnbr(:)
      integer, allocatable :: okey(:),rkey(:),ckey(:),hkey(:),head(:),nxt(:)
      real(WP), allocatable :: rev(:)
      integer(1), allocatable :: rdmg(:)
      integer, dimension(0:nproc-1) :: sc,rc,sd,rd,sc3,rc3,sd3,rd3
      integer, dimension(3) :: nmax,nc
      integer :: d,r,n1,n2,n3,i,k,m,noff,nrecv,ncand,nhe,pass,ic,jc,kc,c1,c2,c3,ierr
      character(len=str_long) :: message

      ! Directory over the node partition (persistent)
      call this%dir%finalize()
      call this%dir%register(this%nown,this%gid(1:this%nown))

      ! Owned reference bounds, exchanged globally
      bl=huge(1.0_WP); bh=-huge(1.0_WP)
      do i=1,this%nown
         bl=min(bl,this%x0(:,i)); bh=max(bh,this%x0(:,i))
      end do
      allocate(allb(6,0:nproc-1))
      call MPI_ALLGATHER([bl,bh],6,MPI_REAL_WP,allb,6,MPI_REAL_WP,comm,ierr)
      do d=1,3
         nmax(d)=0
         if (this%per(d).and.this%Ldom(d).gt.0.0_WP) nmax(d)=min(4,int(this%delta/this%Ldom(d))+1)
      end do

      ! Offers of shifted reference positions (two passes: count, fill)
      do pass=1,2
         sc=0
         do r=0,nproc-1
            do n3=-nmax(3),nmax(3); do n2=-nmax(2),nmax(2); do n1=-nmax(1),nmax(1)
               if (r.eq.rank.and.n1.eq.0.and.n2.eq.0.and.n3.eq.0) cycle
               shift=[real(n1,WP)*this%Ldom(1),real(n2,WP)*this%Ldom(2),real(n3,WP)*this%Ldom(3)]
               if (any(bl+shift-this%delta.gt.allb(4:6,r)).or.any(bh+shift+this%delta.lt.allb(1:3,r))) cycle
               do i=1,this%nown
                  pos_s=this%x0(:,i)+shift
                  if (any(pos_s.lt.allb(1:3,r)-this%delta).or.any(pos_s.gt.allb(4:6,r)+this%delta)) cycle
                  sc(r)=sc(r)+1
                  if (pass.eq.2) then
                     ogid(sd(r)+sc(r))=this%gid(i)
                     okey(sd(r)+sc(r))=(n1+128)+(n2+128)*256+(n3+128)*65536
                     opos(:,sd(r)+sc(r))=pos_s
                  end if
               end do
            end do; end do; end do
         end do
         if (pass.eq.1) then
            sd(0)=0
            do r=1,nproc-1
               sd(r)=sd(r-1)+sc(r-1)
            end do
            noff=sum(sc)
            allocate(ogid(max(noff,1)),okey(max(noff,1)),opos(3,max(noff,1)))
         end if
      end do
      call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
      rd(0)=0
      do r=1,nproc-1
         rd(r)=rd(r-1)+rc(r-1)
      end do
      nrecv=sum(rc)
      allocate(rgid(max(nrecv,1)),rkey(max(nrecv,1)),rpos(3,max(nrecv,1)))
      call MPI_ALLTOALLV(ogid,sc,sd,MPI_INTEGER8,rgid,rc,rd,MPI_INTEGER8,comm,ierr)
      call MPI_ALLTOALLV(okey,sc,sd,MPI_INTEGER, rkey,rc,rd,MPI_INTEGER, comm,ierr)
      sc3=3*sc; sd3=3*sd; rc3=3*rc; rd3=3*rd
      call MPI_ALLTOALLV(opos,sc3,sd3,MPI_REAL_WP,rpos,rc3,rd3,MPI_REAL_WP,comm,ierr)
      deallocate(ogid,okey,opos)

      ! Candidate set = owned nodes (zero offset) + received offers
      ncand=this%nown+nrecv
      allocate(cgid(max(ncand,1)),ckey(max(ncand,1)),cpos(3,max(ncand,1)))
      do i=1,this%nown
         cgid(i)=this%gid(i); ckey(i)=PDHALO_KEY0; cpos(:,i)=this%x0(:,i)
      end do
      do i=1,nrecv
         cgid(this%nown+i)=rgid(i); ckey(this%nown+i)=rkey(i); cpos(:,this%nown+i)=rpos(:,i)
      end do
      deallocate(rgid,rkey,rpos)

      ! Bin candidates; generate each owned row directly (two passes)
      bl=huge(1.0_WP); bh=-huge(1.0_WP)
      do m=1,ncand
         bl=min(bl,cpos(:,m)); bh=max(bh,cpos(:,m))
      end do
      call setup_bins(bl,bh,this%delta,gl,hcell,nc)
      allocate(head(nc(1)*nc(2)*nc(3)),nxt(max(ncand,1)))
      head=0
      do m=1,ncand
         k=cell_of(cpos(:,m),gl,hcell,nc)
         nxt(m)=head(k); head(k)=m
      end do
      do pass=1,2
         nhe=0
         do i=1,this%nown
            ic=min(nc(1),max(1,int((this%x0(1,i)-gl(1))/hcell(1))+1))
            jc=min(nc(2),max(1,int((this%x0(2,i)-gl(2))/hcell(2))+1))
            kc=min(nc(3),max(1,int((this%x0(3,i)-gl(3))/hcell(3))+1))
            do c3=max(1,kc-1),min(nc(3),kc+1); do c2=max(1,jc-1),min(nc(2),jc+1); do c1=max(1,ic-1),min(nc(1),ic+1)
               m=head(c1+nc(1)*(c2-1)+nc(1)*nc(2)*(c3-1))
               do while (m.gt.0)
                  if (m.ne.i) then
                     if (sum((cpos(:,m)-this%x0(:,i))**2).le.this%delta**2) then
                        nhe=nhe+1
                        if (pass.eq.2) then
                           hnode(nhe)=this%gid(i)
                           hnbr(nhe) =cgid(m)
                           hkey(nhe) =ckey(m)
                        end if
                     end if
                  end if
                  m=nxt(m)
               end do
            end do; end do; end do
         end do
         if (pass.eq.1) allocate(hnode(max(nhe,1)),hnbr(max(nhe,1)),hkey(max(nhe,1)))
      end do
      deallocate(cgid,ckey,cpos,head,nxt,allb)

      ! Assemble with zero inelastic state (entries are already local rows)
      allocate(rev(max(nhe,1)),rdmg(max(nhe,1)))
      rev=0.0_WP; rdmg=0_1
      call this%assemble(nhe,hnode,hnbr,hkey,rev,rdmg)
      deallocate(hnode,hnbr,hkey,rev,rdmg)
      if (amRoot) then
         write(message,'("[",a,"] detect_families: ",i0," half-entries (~2x bonds)")') trim(this%name),this%nbond
         call log(message)
      end if
   end subroutine detect_families

   !> Assemble the CSR families, halo plan, and reference state from LOCAL
   !> half-entry arrays (already routed to this rank: every entry's node gid is
   !> owned here). Per-entry inelastic state (dmg, e_v) travels with the
   !> entries -- zeros for a fresh connect, loaded values on restart. Shared by
   !> connect and read_state; collective.
   subroutine assemble(this,rn,rnode,rnbr,rkey,rev,rdmg)
      use parallel, only: comm,nproc
      use messager, only: die
      use mpi_f08
      implicit none
      class(pdsolver), intent(inout) :: this
      integer, intent(in) :: rn
      integer(I8), intent(in) :: rnode(:),rnbr(:)
      integer, intent(in) :: rkey(:)
      real(WP), intent(in) :: rev(:)
      integer(1), intent(in) :: rdmg(:)
      integer, allocatable :: ridx(:),perm(:)
      integer :: i,s,ierr

      if (allocated(this%ptr)) deallocate(this%ptr)
      if (allocated(this%lst)) deallocate(this%lst)
      if (allocated(this%dmg)) deallocate(this%dmg)
      if (allocated(this%e_v)) deallocate(this%e_v)
      if (allocated(this%td2)) deallocate(this%td2)
      if (allocated(this%td2a)) deallocate(this%td2a)

      ! Resolve each received entry's node to an owned index
      allocate(ridx(max(rn,1)),perm(max(rn,1)))
      do i=1,rn
         ridx(i)=this%ohash%lookup(rnode(i))
         if (ridx(i).lt.1) call die('[pdsolver connect] half-entry routed to a rank that does not own its node')
         perm(i)=i
      end do

      ! Deterministic CSR order: sort by (node index, neighbor gid, image key)
      if (rn.gt.1) call sort3_perm(ridx,rnbr,rkey,perm,1,rn)

      ! Row pointers
      allocate(this%ptr(this%nown+1))
      row_pointers: block
         integer, allocatable :: cnt(:)
         allocate(cnt(this%nown)); cnt=0
         do i=1,rn
            cnt(ridx(i))=cnt(ridx(i))+1
         end do
         this%ptr(1)=1
         do i=1,this%nown
            this%ptr(i+1)=this%ptr(i)+cnt(i)
         end do
         deallocate(cnt)
      end block row_pointers

      ! Classify entries (owned direct vs halo reference), dedupe references,
      ! build the halo plan, and finalize the CSR neighbor indices
      build_refs_and_halo: block
         integer(I8), allocatable :: refgid(:),ugid(:)
         integer, allocatable :: refkey(:),refpos(:),rperm(:),zeros(:)
         integer, allocatable :: ukey(:),uowner(:),uslot(:)
         integer :: nref,nuniq,lid,u
         allocate(this%lst(max(rn,1)))
         allocate(this%dmg(max(rn,1))); this%dmg=0_1
         allocate(this%e_v(max(rn,1))); this%e_v=0.0_WP
         ! Per-entry inelastic state follows the deterministic CSR order
         do s=1,rn
            this%dmg(s)=rdmg(perm(s))
            this%e_v(s)=rev(perm(s))
         end do
         allocate(this%td2(max(this%nown,1)),this%td2a(max(this%nown,1)))
         this%td2=0.0_WP; this%td2a=0.0_WP
         allocate(refgid(max(rn,1)),refkey(max(rn,1)),refpos(max(rn,1)))
         nref=0
         do s=1,rn
            i=perm(s)
            if (rkey(i).eq.PDHALO_KEY0) then
               lid=this%ohash%lookup(rnbr(i))
               if (lid.ge.1) then
                  this%lst(s)=lid       ! owned, zero image offset: direct index
                  cycle
               end if
            end if
            nref=nref+1
            refgid(nref)=rnbr(i); refkey(nref)=rkey(i); refpos(nref)=s
         end do
         ! Unique (gid,key) references, deterministic order
         allocate(rperm(max(nref,1)),zeros(max(nref,1)))
         zeros=0
         do i=1,nref
            rperm(i)=i
         end do
         if (nref.gt.1) call sort3_perm(zeros,refgid,refkey,rperm,1,nref)
         allocate(ugid(max(nref,1)),ukey(max(nref,1)))
         nuniq=0
         do s=1,nref
            i=rperm(s)
            if (s.eq.1) then
               nuniq=1; ugid(1)=refgid(i); ukey(1)=refkey(i)
            else if (refgid(i).ne.refgid(rperm(s-1)).or.refkey(i).ne.refkey(rperm(s-1))) then
               nuniq=nuniq+1; ugid(nuniq)=refgid(i); ukey(nuniq)=refkey(i)
            end if
            this%lst(refpos(i))=-nuniq   ! provisional: -(unique ref id)
         end do
         ! Owners of the unique references, then the persistent halo plan
         allocate(uowner(max(nuniq,1)),uslot(max(nuniq,1)))
         call this%dir%query(nuniq,ugid,uowner)
         call this%halo%build(this%nown,this%ohash,nuniq,ugid,ukey,uowner,this%Ldom,this%per,uslot)
         this%nhalo=this%halo%nhalo
         this%ntot=this%nown+this%nhalo
         ! Finalize CSR: provisional negatives -> halo slot indices
         do s=1,rn
            if (this%lst(s).lt.0) this%lst(s)=this%nown+uslot(-this%lst(s))
         end do
         ! Extend node arrays to include halo slots; stamp halo gids
         extend_arrays: block
            integer(I8), allocatable :: g2(:)
            real(WP), allocatable :: a2(:,:)
            allocate(g2(max(this%ntot,1))); g2(1:this%nown)=this%gid(1:this%nown)
            do u=1,nuniq
               g2(this%nown+uslot(u))=ugid(u)
            end do
            call move_alloc(g2,this%gid)
            allocate(a2(3,max(this%ntot,1))); a2=0.0_WP; a2(:,1:this%nown)=this%x0(:,1:this%nown)
            call move_alloc(a2,this%x0)
            allocate(a2(3,max(this%ntot,1))); a2=0.0_WP; a2(:,1:this%nown)=this%y(:,1:this%nown)
            call move_alloc(a2,this%y)
            allocate(a2(3,max(this%ntot,1))); a2=0.0_WP; a2(:,1:this%nown)=this%f(:,1:this%nown)
            call move_alloc(a2,this%f)
            extend_volume: block
               real(WP), allocatable :: v2(:)
               allocate(v2(max(this%ntot,1))); v2=0.0_WP; v2(1:this%nown)=this%vol(1:this%nown)
               call move_alloc(v2,this%vol)
            end block extend_volume
         end block extend_arrays
         deallocate(refgid,refkey,refpos,rperm,zeros,ugid,ukey,uowner,uslot)
      end block build_refs_and_halo
      deallocate(ridx,perm)

      ! Fill halo reference positions ONCE, pre-shifted by the image offsets
      ! (x0 is static; this is the only x0 exchange of the entire run), and
      ! the halo per-node volumes (also static)
      call this%halo%update(this%x0,3,shifted=.true.)
      call this%halo%update1(this%vol)
      this%y(:,this%nown+1:this%ntot)=this%x0(:,this%nown+1:this%ntot)

      ! Life status (exit-through-open-face handling). Exchanged over the halo
      ! ONLY on substeps where a death occurs somewhere; steady state is free.
      if (allocated(this%alive)) deallocate(this%alive)
      allocate(this%alive(max(this%ntot,1))); this%alive=1.0_WP
      this%watch_exit=(this%dom_hi(1).gt.this%dom_lo(1)).and. &
      &  any((.not.this%per).and.(this%lo_bc.eq.0.or.this%hi_bc.eq.0))

      ! Stamp the reference weighted volume
      call this%compute_mw()

      ! Silling-Askari critical time step (Peridigm form, 3D bond-based
      ! micromodulus C = 18K/(pi*delta^4)):
      !   dt_crit_i = sqrt(2*rho / sum_family(V_j * C / zeta)), global min.
      ! DIAGNOSTIC only for now -- reported at init, does not bind dt. The
      ! micromodulus constant is 3D-based; in quasi-2D slabs treat it as
      ! indicative.
      critical_dt: block
         use mathtools, only: Pi
         use messager,  only: log
         use string,    only: str_long
         use parallel,  only: amRoot,MPI_REAL_WP
         real(WP) :: K_bulk,Cmicro,denom,zeta,dtc
         character(len=str_long) :: message
         integer :: i,e,j
         K_bulk=this%elastic_modulus/(3.0_WP*(1.0_WP-2.0_WP*this%poisson_ratio))
         Cmicro=18.0_WP*K_bulk/(Pi*this%delta**4)
         dtc=huge(1.0_WP)
         do i=1,this%nown
            denom=0.0_WP
            do e=this%ptr(i),this%ptr(i+1)-1
               j=this%lst(e)
               zeta=sqrt(sum((this%x0(:,j)-this%x0(:,i))**2))
               if (zeta.gt.0.0_WP) denom=denom+this%vol(j)*Cmicro/zeta
            end do
            if (denom.gt.0.0_WP) dtc=min(dtc,sqrt(2.0_WP*this%rho/denom))
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,dtc,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
         this%dtcrit=dtc
         if (amRoot) then
            write(message,'("[",a,"] Silling-Askari critical dt = ",es12.5," (diagnostic)")') trim(this%name),this%dtcrit
            call log(message)
         end if
      end block critical_dt

      ! Global half-entry count for logging (= 2*bonds - self-image bonds)
      count_bonds: block
         use parallel, only: comm
         integer(I8) :: nhe8
         nhe8=int(rn,I8)
         call MPI_ALLREDUCE(MPI_IN_PLACE,nhe8,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
         this%nbond=nhe8   ! total half-entries; = 2*bonds - self-image bonds
         ! Exact bond census (lower-gid rule; positive-offset self-images)
         count_nb: block
            integer :: i2,e2,j2
            this%nb=0_I8
            do i2=1,this%nown
               do e2=this%ptr(i2),this%ptr(i2+1)-1
                  j2=this%lst(e2)
                  if (this%gid(i2).lt.this%gid(j2)) then
                     this%nb=this%nb+1_I8
                  else if (this%gid(i2).eq.this%gid(j2).and.j2.gt.this%nown) then
                     if (shift_positive(this%halo%shift(:,j2-this%nown))) this%nb=this%nb+1_I8
                  end if
               end do
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%nb,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
         end block count_nb
      end block count_bonds
   end subroutine assemble

   !> Weighted volume: mw_i = sum_family w(zeta)*zeta^2*V_j (reference state;
   !> never updated by damage)
   subroutine compute_mw(this)
      implicit none
      class(pdsolver), intent(inout) :: this
      integer :: i,e,j
      real(WP) :: zeta
      do i=1,this%nown
         this%mw(i)=0.0_WP
         do e=this%ptr(i),this%ptr(i+1)-1
            j=this%lst(e)
            zeta=sqrt(sum((this%x0(:,j)-this%x0(:,i))**2))
            this%mw(i)=this%mw(i)+omega(zeta,this%delta)*zeta**2*this%vol(j)
         end do
      end do
   end subroutine compute_mw

   !> Dimension-aware LPS constitutive coefficients. psi_fac sets the J2 yield
   !> threshold on the family deviatoric force-state norm: yield when
   !> ||t_dev||^2 > psi_fac*sigma_yield^2/mw (Mitchell OSB).
   subroutine lps_coefs(this,fdim,coef_vol,coef_dev,psi_fac)
      implicit none
      class(pdsolver), intent(in) :: this
      real(WP), intent(out) :: fdim,coef_vol,coef_dev
      real(WP), intent(out), optional :: psi_fac
      real(WP) :: K_bulk,mu_shear
      integer :: ndim
      ndim=3-count(this%collapsed)
      K_bulk  =this%elastic_modulus/(3.0_WP*(1.0_WP-2.0_WP*this%poisson_ratio))
      mu_shear=this%elastic_modulus/(2.0_WP*(1.0_WP+this%poisson_ratio))
      select case (ndim)
      case (3)
         fdim=3.0_WP; coef_vol=3.0_WP*K_bulk;                   coef_dev=15.0_WP*mu_shear
         if (present(psi_fac)) psi_fac=5.0_WP
      case (2)
         fdim=2.0_WP; coef_vol=2.0_WP*(K_bulk+mu_shear/3.0_WP); coef_dev= 8.0_WP*mu_shear
         if (present(psi_fac)) psi_fac=8.0_WP/3.0_WP
      case default
         fdim=1.0_WP; coef_vol=this%elastic_modulus;            coef_dev= 0.0_WP
         if (present(psi_fac)) psi_fac=0.0_WP
      end select
   end subroutine lps_coefs

   !> Velocity-Verlet step: half-kick + drift, halo position update,
   !> dilatation gather, node-centered force sweep, halo force reduce,
   !> second half-kick. Mirrors amrpd%advance minus contact/VF (stage 1).
   subroutine advance(this,dt)
      use parallel, only: parallel_time
      implicit none
      class(pdsolver), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP) :: rho_inv,fdim,cvol,cdev,t0
      real(WP) :: zeta,dY,e_b,t,w
      real(WP) :: psi_fac,sY2,decay,e_d,td,beta,e_e,over
      logical :: plastic,do_j2
      real(WP), dimension(3) :: acc,dxv,fx
      integer :: i,e,j

      rho_inv=1.0_WP/this%rho
      call this%lps_coefs(fdim,cvol,cdev,psi_fac)
      ! Viscoplastic setup: decay is loop-invariant (exact exponential update,
      ! unconditionally stable -- no viscous CFL)
      sY2=this%sigma_yield**2
      plastic=(this%tau.gt.0.0_WP.and.this%tau.lt.huge(1.0_WP))
      do_j2=(this%sigma_yield.gt.0.0_WP)
      decay=0.0_WP
      if (plastic) decay=exp(-dt/this%tau)

      ! First half-kick and drift (owned nodes)
      t0=parallel_time()
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         acc=this%gravity+(this%f(:,i)+this%ff(:,i))*rho_inv
         if (iand(this%flag(i),PDC_INTEGRATES).ne.0) this%v(:,i)=this%v(:,i)+0.5_WP*dt*acc
         if (this%collapsed(1)) this%v(1,i)=0.0_WP
         if (this%collapsed(2)) this%v(2,i)=0.0_WP
         if (this%collapsed(3)) this%v(3,i)=0.0_WP
         if (iand(this%flag(i),PDC_MOVES).ne.0) this%y(:,i)=this%y(:,i)+dt*this%v(:,i)
      end do
      this%wt_kick=this%wt_kick+(parallel_time()-t0)

      ! Exit handling: nodes drifting out through an OPEN non-periodic face die
      ! (amrpd drops them at Redistribute; here they are flagged and muted).
      ! The death-count allreduce runs only when exits are possible at all, and
      ! the mute propagation only on substeps where a death actually occurred.
      if (this%watch_exit) then
         death_watch: block
            use parallel, only: comm
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_INTEGER
            integer :: nd,d,e,ierr
            logical :: out
            nd=0
            do i=1,this%nown
               if (this%flag(i).eq.PDC_IS_DEAD) cycle
               out=.false.
               do d=1,3
                  if (this%per(d)) cycle
                  if (this%lo_bc(d).eq.0.and.this%y(d,i).lt.this%dom_lo(d)) out=.true.
                  if (this%hi_bc(d).eq.0.and.this%y(d,i).gt.this%dom_hi(d)) out=.true.
               end do
               if (out) then
                  this%flag(i)=PDC_IS_DEAD
                  this%v(:,i)=0.0_WP
                  this%alive(i)=0.0_WP
                  nd=nd+1
               end if
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,nd,1,MPI_INTEGER,MPI_SUM,comm,ierr)
            if (nd.gt.0) then
               ! Propagate life status to halo slots, then permanently mute
               ! every entry touching a dead node (dmg=2: distinct from broken,
               ! so damage statistics stay honest -- amrpd does not count
               ! dropped-particle bonds as damage either)
               call this%halo%update1(this%alive)
               do i=1,this%nown
                  if (this%flag(i).eq.PDC_IS_DEAD) then
                     do e=this%ptr(i),this%ptr(i+1)-1
                        if (this%dmg(e).eq.0_1) this%dmg(e)=2_1
                     end do
                  else
                     do e=this%ptr(i),this%ptr(i+1)-1
                        if (this%dmg(e).eq.0_1.and.this%alive(this%lst(e)).lt.0.5_WP) this%dmg(e)=2_1
                     end do
                  end if
               end do
               ! Force a contact broad-phase rebuild so no candidate list
               ! references a corpse (dead nodes are excluded from offers/bins)
               if (this%use_contact.and.allocated(this%ylast)) deallocate(this%ylast)
            end if
         end block death_watch
      end if

      ! Refresh halo positions (owner y -> slots, with image shifts)
      t0=parallel_time()
      call this%halo%update(this%y,3,shifted=.true.)
      this%wt_halo=this%wt_halo+(parallel_time()-t0)

      ! Contact service: displacement-triggered broad phase (rebuilds contact
      ! halo + candidate list when cumulative drift exhausts the skin), then
      ! per-substep refresh of contact-slot positions AND velocities (the only
      ! place velocity crosses ranks; the graph halo never carries it)
      if (this%use_contact) then
         t0=parallel_time()
         call this%contact_broadphase()
         this%wt_broad=this%wt_broad+(parallel_time()-t0)
         t0=parallel_time()
         call this%chalo%update(this%y,3,shifted=.true.)
         call this%chalo%update(this%v,3,shifted=.false.)
         this%wt_halo=this%wt_halo+(parallel_time()-t0)
      end if

      ! Dilatation (pure gather; own family only; broken entries excluded --
      ! breaks happen in the force sweep AFTER this, matching amrpd's ordering)
      t0=parallel_time()
      do i=1,this%nown
         this%theta(i)=0.0_WP
         do e=this%ptr(i),this%ptr(i+1)-1
            if (this%dmg(e).ne.0_1) cycle
            j=this%lst(e)
            zeta=sqrt(sum((this%x0(:,j)-this%x0(:,i))**2))
            dY  =sqrt(sum((this%y(:,j) -this%y(:,i) )**2))
            e_b=dY-zeta
            this%theta(i)=this%theta(i)+omega(zeta,this%delta)*zeta*e_b*this%vol(j)
         end do
         if (this%mw(i).gt.0.0_WP) then
            this%theta(i)=fdim*this%theta(i)/this%mw(i)
         else
            this%theta(i)=0.0_WP
         end if
      end do
      this%wt_dil=this%wt_dil+(parallel_time()-t0)

      ! Node-centered force sweep: each row computes its OWN force state t
      ! (own theta, own mw) and scatters +t/-t; the neighbor's t arrives from
      ! the neighbor's own row (locally or via the halo reduce below).
      t0=parallel_time()
      this%f=0.0_WP
      do i=1,this%nown
         if (this%mw(i).le.0.0_WP) cycle
         do e=this%ptr(i),this%ptr(i+1)-1
            if (this%dmg(e).ne.0_1) cycle
            j=this%lst(e)
            zeta=sqrt(sum((this%x0(:,j)-this%x0(:,i))**2))
            dxv=this%y(:,j)-this%y(:,i)
            dY=sqrt(sum(dxv**2))
            if (dY.le.0.0_WP) cycle
            e_b=dY-zeta
            ! Brittle break on total stretch (e > s0*zeta), irreversible.
            ! Each row breaks its OWN half-entry and increments its OWN node's
            ! damage by 1/nb0 (nb0 = reference row length); the counterpart row
            ! breaks its half independently -- the criterion is symmetric in
            ! the endpoints, so both halves break in the same substep (for
            ! image bonds, up to shift-association roundoff: a 1-ulp-marginal
            ! bond may break one substep apart, a benign local transient --
            ! the intact half still applies its +/- pair, conserving momentum).
            if (e_b.gt.this%s0*zeta) then
               this%dmg(e)=1_1
               this%damage(i)=this%damage(i)+1.0_WP/real(this%ptr(i+1)-this%ptr(i),WP)
               cycle
            end if
            w=omega(zeta,this%delta)
            ! Deviatoric split: e_d carries this HALF-ENTRY's inelastic stretch
            ! e_v (per-side history: own theta, own mw -- Peridigm form; e_v=0
            ! recovers canonical elastic LPS bit-for-bit)
            e_d=e_b-this%theta(i)*zeta/fdim
            td=w/this%mw(i)*cdev*(e_d-this%visc_lambda*this%e_v(e))
            t =w/this%mw(i)*cvol*this%theta(i)*zeta+td
            ! J2 family norm: pure own-row gather (no communication)
            if (do_j2) this%td2a(i)=this%td2a(i)+td*td*this%vol(j)
            ! Pair contribution from THIS row's force state (Peridigm volumes:
            ! +t*V_j to self, -t*V_i to the neighbor)
            fx=t*dxv/dY
            this%f(:,i)=this%f(:,i)+fx*this%vol(j)
            this%f(:,j)=this%f(:,j)-fx*this%vol(i)
            ! Per-side viscoplastic flow of e_v (exact exponential). Two yield
            ! criteria, as in amrpd:
            !   sigma_yield>0: J2 radial return toward the yield surface when
            !     the PREVIOUS substep's family norm exceeds psi_fac*sY^2/mw,
            !     Perzyna-regularized by (1-decay); tau->0 recovers Peridigm's
            !     rate-independent return.
            !   else: per-bond overstress (yield_stretch=0 -> pure Maxwell).
            if (plastic) then
               if (do_j2) then
                  beta=1.0_WP
                  if (this%td2(i)*this%mw(i).gt.psi_fac*sY2) beta=sqrt(psi_fac*sY2/(this%td2(i)*this%mw(i)))
                  this%e_v(e)=this%e_v(e)+(1.0_WP-beta)*(e_d-this%e_v(e))*(1.0_WP-decay)
               else
                  e_e=e_d-this%e_v(e)
                  over=abs(e_e)-this%yield_stretch*zeta
                  if (over.gt.0.0_WP) this%e_v(e)=this%e_v(e)+sign(over*(1.0_WP-decay),e_e)
               end if
            end if
         end do
      end do
      ! Publish this substep's J2 norm (read by the NEXT substep's return)
      if (do_j2) then
         this%td2(1:this%nown)=this%td2a(1:this%nown)
         this%td2a(1:this%nown)=0.0_WP
      end if
      this%wt_force=this%wt_force+(parallel_time()-t0)

      ! Assemble cross-rank pair forces (halo slots -> owners, add)
      t0=parallel_time()
      call this%halo%reduce(this%f,3)
      this%wt_reduce=this%wt_reduce+(parallel_time()-t0)

      ! Short-range contact (walls + particle-particle), gather-only: adds
      ! into owned f, no reduction (amrpd ordering: after the bond force)
      if (this%use_contact) then
         t0=parallel_time()
         call this%contact_narrow(dt)
         this%wt_contact=this%wt_contact+(parallel_time()-t0)
      end if

      ! Second half-kick with the fresh force
      t0=parallel_time()
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         if (iand(this%flag(i),PDC_INTEGRATES).ne.0) then
            acc=this%gravity+(this%f(:,i)+this%ff(:,i))*rho_inv
            this%v(:,i)=this%v(:,i)+0.5_WP*dt*acc
         end if
         if (this%collapsed(1)) this%v(1,i)=0.0_WP
         if (this%collapsed(2)) this%v(2,i)=0.0_WP
         if (this%collapsed(3)) this%v(3,i)=0.0_WP
      end do
      this%wt_kick=this%wt_kick+(parallel_time()-t0)
   end subroutine advance

   !> Contact broad phase: displacement-triggered rebuild of the contact halo
   !> and the candidate CSR. The trigger is one scalar allreduce per substep so
   !> the (collective) rebuild decision is rank-consistent. rbuild =
   !> 1.2*contact_dist + 2*cskin: engagement reach is bounded by d_eff*(1+0.2)
   !> (the r_influ clip) and two nodes drifting cskin each can close 2*cskin
   !> between rebuilds, so the candidate set provably contains every pair that
   !> can produce force before the next rebuild.
   !>
   !> Discovery: allgather of per-rank owned-node bounds; for each (rank,
   !> periodic-image offset) whose shifted bounds approach mine within rbuild,
   !> OFFER my owned nodes in range as (gid, image key, shifted position). The
   !> receiver keeps offers with an owned node within rbuild (binned test) and
   !> builds the contact halo from the kept references via the standard pdhalo
   !> protocol (chalo%nown = ntot, so contact slots append after graph slots).
   !> Candidates are then binned over OWNED + CONTACT slots only -- the graph
   !> halo is never binned, so bonded remote partners arrive as contact slots
   !> when in range and double-counting is structurally impossible.
   subroutine contact_broadphase(this)
      use parallel, only: comm,rank,nproc,MPI_REAL_WP
      use mpi_f08
      use messager, only: die
      implicit none
      class(pdsolver), intent(inout) :: this
      real(WP) :: rbuild,drift
      integer :: i,ierr

      if (this%contact_dist.le.0.0_WP) call die('[pdsolver contact] use_contact requires contact_dist > 0')
      if (this%cskin.le.0.0_WP) this%cskin=0.5_WP*this%contact_dist
      rbuild=1.2_WP*this%contact_dist+2.0_WP*this%cskin

      ! Displacement trigger (collective decision)
      if (allocated(this%ylast)) then
         drift=0.0_WP
         do i=1,this%nown
            drift=max(drift,sum((this%y(:,i)-this%ylast(:,i))**2))
         end do
         drift=sqrt(drift)
      else
         drift=huge(1.0_WP)
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,drift,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      if (drift.le.this%cskin) return
      this%nrebuild=this%nrebuild+1_I8

      rebuild: block
         real(WP), dimension(3) :: bl,bh,shift,pos_s
         real(WP), allocatable :: allb(:,:),opos(:,:),rpos(:,:),kpos(:,:)
         integer(I8), allocatable :: ogid(:),rgid(:),kgid(:)
         integer, allocatable :: okey(:),rkey(:),kkey(:),kowner(:),slot(:)
         integer, dimension(0:nproc-1) :: sc,rc,sd,rd
         integer, dimension(0:nproc-1) :: sc3,rc3,sd3,rd3
         integer :: nmax(3),d,r,n1,n2,n3,noff,nrecv,nkeep,k,pass
         ! Binning workspace (owned nodes for offer filtering, then combined
         ! set for the candidate CSR)
         real(WP), dimension(3) :: gl,hcell
         integer, dimension(3) :: nc
         integer, allocatable :: head(:),nxt(:)

         ! Owned bounds and their global exchange
         bl=huge(1.0_WP); bh=-huge(1.0_WP)
         do i=1,this%nown
            bl=min(bl,this%y(:,i)); bh=max(bh,this%y(:,i))
         end do
         allocate(allb(6,0:nproc-1))
         call MPI_ALLGATHER([bl,bh],6,MPI_REAL_WP,allb,6,MPI_REAL_WP,comm,ierr)

         ! Admissible periodic-image offsets for contact range
         do d=1,3
            nmax(d)=0
            if (this%per(d).and.this%Ldom(d).gt.0.0_WP) nmax(d)=min(4,int(rbuild/this%Ldom(d))+1)
         end do

         ! Offers: two passes (count, then fill), grouped by destination rank
         do pass=1,2
            sc=0
            do r=0,nproc-1
               do n3=-nmax(3),nmax(3); do n2=-nmax(2),nmax(2); do n1=-nmax(1),nmax(1)
                  if (r.eq.rank.and.n1.eq.0.and.n2.eq.0.and.n3.eq.0) cycle
                  shift=[real(n1,WP)*this%Ldom(1),real(n2,WP)*this%Ldom(2),real(n3,WP)*this%Ldom(3)]
                  ! Shifted-bounds proximity prefilter
                  if (any(bl+shift-rbuild.gt.allb(4:6,r)).or.any(bh+shift+rbuild.lt.allb(1:3,r))) cycle
                  do i=1,this%nown
                     if (this%flag(i).eq.PDC_IS_DEAD) cycle
                     pos_s=this%y(:,i)+shift
                     if (any(pos_s.lt.allb(1:3,r)-rbuild).or.any(pos_s.gt.allb(4:6,r)+rbuild)) cycle
                     sc(r)=sc(r)+1
                     if (pass.eq.2) then
                        ogid(sd(r)+sc(r))=this%gid(i)
                        okey(sd(r)+sc(r))=(n1+128)+(n2+128)*256+(n3+128)*65536
                        opos(:,sd(r)+sc(r))=pos_s
                     end if
                  end do
               end do; end do; end do
            end do
            if (pass.eq.1) then
               sd(0)=0
               do r=1,nproc-1
                  sd(r)=sd(r-1)+sc(r-1)
               end do
               noff=sum(sc)
               allocate(ogid(max(noff,1)),okey(max(noff,1)),opos(3,max(noff,1)))
            end if
         end do

         ! Exchange offers (gid, key, shifted position)
         call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
         rd(0)=0
         do r=1,nproc-1
            rd(r)=rd(r-1)+rc(r-1)
         end do
         nrecv=sum(rc)
         allocate(rgid(max(nrecv,1)),rkey(max(nrecv,1)),rpos(3,max(nrecv,1)))
         call MPI_ALLTOALLV(ogid,sc,sd,MPI_INTEGER8,rgid,rc,rd,MPI_INTEGER8,comm,ierr)
         call MPI_ALLTOALLV(okey,sc,sd,MPI_INTEGER, rkey,rc,rd,MPI_INTEGER, comm,ierr)
         sc3=3*sc; sd3=3*sd; rc3=3*rc; rd3=3*rd
         call MPI_ALLTOALLV(opos,sc3,sd3,MPI_REAL_WP,rpos,rc3,rd3,MPI_REAL_WP,comm,ierr)
         deallocate(ogid,okey,opos)

         ! Filter offers: keep those with an owned node within rbuild.
         ! Bin owned nodes (cell size >= rbuild so a +/-1 cell sweep suffices;
         ! dims clamped so degenerate/huge extents stay bounded).
         call setup_bins(bl,bh,rbuild,gl,hcell,nc)
         allocate(head(nc(1)*nc(2)*nc(3)),nxt(max(this%nown,1)))
         head=0
         do i=1,this%nown
            if (this%flag(i).eq.PDC_IS_DEAD) cycle
            k=cell_of(this%y(:,i),gl,hcell,nc)
            nxt(i)=head(k); head(k)=i
         end do
         allocate(kgid(max(nrecv,1)),kkey(max(nrecv,1)),kowner(max(nrecv,1)),kpos(3,max(nrecv,1)))
         nkeep=0
         do r=0,nproc-1
            do i=rd(r)+1,rd(r)+rc(r)
               if (near_owned(rpos(:,i),rbuild,gl,hcell,nc,head,nxt)) then
                  nkeep=nkeep+1
                  kgid(nkeep)=rgid(i); kkey(nkeep)=rkey(i); kowner(nkeep)=r; kpos(:,nkeep)=rpos(:,i)
               end if
            end do
         end do
         deallocate(rgid,rkey,rpos,head,nxt)

         ! Rebuild the contact halo (slots append after graph slots: nown=ntot)
         call this%chalo%finalize()
         allocate(slot(max(nkeep,1)))
         call this%chalo%build(this%ntot,this%ohash,nkeep,kgid,kkey,kowner,this%Ldom,this%per,slot)
         this%nchalo=this%chalo%nhalo

         ! Extend y and v to cover contact slots; stamp slot positions from the
         ! kept offers (current values -- chalo%update refreshes each substep)
         resize_state: block
            real(WP), allocatable :: a2(:,:)
            integer :: ntc
            ntc=this%ntot+this%nchalo
            allocate(a2(3,max(ntc,1))); a2=0.0_WP
            a2(:,1:this%ntot)=this%y(:,1:this%ntot)
            call move_alloc(a2,this%y)
            allocate(a2(3,max(ntc,1))); a2=0.0_WP
            a2(:,1:this%nown)=this%v(:,1:this%nown)
            call move_alloc(a2,this%v)
            do k=1,nkeep
               this%y(:,this%ntot+slot(k))=kpos(:,k)
            end do
         end block resize_state
         deallocate(kgid,kkey,kowner,kpos,slot)

         ! Candidate CSR over the contact-visible set: owned nodes (indices
         ! 1..nown) + contact slots (ntot+1..ntot+nchalo). Two passes.
         candidates: block
            integer :: ns,m,jj,cnt,ic,jc,kc,c1,c2,c3
            integer, allocatable :: midx(:)
            real(WP), dimension(3) :: blc,bhc
            ns=this%nown+this%nchalo
            allocate(midx(max(ns,1)))
            do m=1,this%nown
               midx(m)=m
            end do
            do m=1,this%nchalo
               midx(this%nown+m)=this%ntot+m
            end do
            blc=bl; bhc=bh
            do m=this%nown+1,ns
               blc=min(blc,this%y(:,midx(m))); bhc=max(bhc,this%y(:,midx(m)))
            end do
            call setup_bins(blc,bhc,rbuild,gl,hcell,nc)
            allocate(head(nc(1)*nc(2)*nc(3)),nxt(max(ns,1)))
            head=0
            do m=1,ns
               if (m.le.this%nown) then
                  if (this%flag(m).eq.PDC_IS_DEAD) cycle
               end if
               k=cell_of(this%y(:,midx(m)),gl,hcell,nc)
               nxt(m)=head(k); head(k)=m
            end do
            if (allocated(this%cptr)) deallocate(this%cptr)
            if (allocated(this%clst)) deallocate(this%clst)
            allocate(this%cptr(this%nown+1))
            do pass=1,2
               do i=1,this%nown
                  cnt=0
                  if (this%flag(i).eq.PDC_IS_DEAD) then
                     if (pass.eq.1) this%cptr(i+1)=0
                     cycle
                  end if
                  ic=min(nc(1),max(1,int((this%y(1,i)-gl(1))/hcell(1))+1))
                  jc=min(nc(2),max(1,int((this%y(2,i)-gl(2))/hcell(2))+1))
                  kc=min(nc(3),max(1,int((this%y(3,i)-gl(3))/hcell(3))+1))
                  do c3=max(1,kc-1),min(nc(3),kc+1); do c2=max(1,jc-1),min(nc(2),jc+1); do c1=max(1,ic-1),min(nc(1),ic+1)
                     m=head(c1+nc(1)*(c2-1)+nc(1)*nc(2)*(c3-1))
                     do while (m.gt.0)
                        jj=midx(m)
                        if (jj.ne.i) then
                           if (sum((this%y(:,jj)-this%y(:,i))**2).le.rbuild**2) then
                              cnt=cnt+1
                              if (pass.eq.2) this%clst(this%cptr(i)+cnt-1)=jj
                           end if
                        end if
                        m=nxt(m)
                     end do
                  end do; end do; end do
                  if (pass.eq.1) this%cptr(i+1)=cnt   ! provisional count
               end do
               if (pass.eq.1) then
                  this%cptr(1)=1
                  do i=1,this%nown
                     this%cptr(i+1)=this%cptr(i)+this%cptr(i+1)
                  end do
                  allocate(this%clst(max(this%cptr(this%nown+1)-1,1)))
               end if
            end do
            deallocate(midx,head,nxt)
         end block candidates

         ! Snapshot positions for the drift trigger
         if (allocated(this%ylast)) deallocate(this%ylast)
         allocate(this%ylast(3,max(this%nown,1)))
         this%ylast(:,1:this%nown)=this%y(:,1:this%nown)
         deallocate(allb)
      end block rebuild

   contains

      !> Any owned node within r of position p? (binned +/-1 cell sweep)
      function near_owned(p,r,gl,h,nc,head,nxt) result(hit)
         real(WP), dimension(3), intent(in) :: p,gl,h
         real(WP), intent(in) :: r
         integer, dimension(3), intent(in) :: nc
         integer, intent(in) :: head(:),nxt(:)
         logical :: hit
         integer :: c(3),d,c1,c2,c3,m
         hit=.false.
         do d=1,3
            c(d)=min(nc(d),max(1,int((p(d)-gl(d))/h(d))+1))
         end do
         do c3=max(1,c(3)-1),min(nc(3),c(3)+1); do c2=max(1,c(2)-1),min(nc(2),c(2)+1); do c1=max(1,c(1)-1),min(nc(1),c(1)+1)
            m=head(c1+nc(1)*(c2-1)+nc(1)*nc(2)*(c3-1))
            do while (m.gt.0)
               if (sum((this%y(:,m)-p)**2).le.r**2) then
                  hit=.true.
                  return
               end if
               m=nxt(m)
            end do
         end do; end do; end do
      end function near_owned

   end subroutine contact_broadphase

   !> Contact narrow phase: soft-sphere walls + particle-particle over the
   !> candidate CSR, gather-only (soft-sphere penalty ported from amrlpt's collision model;
   !> IB contact arrives with the coupling layer). Adds force/volume into owned
   !> f. Walls use e_w with d_eff = 0.5*contact_dist and m_eff = m1; pairs use
   !> e_n with d_eff = contact_dist and m_eff = 0.5*m1 (m1 = rho*vol(i),
   !> matching amrpd's uniform rho*dV on a uniform lattice).
   subroutine contact_narrow(this,dt)
      use mathtools, only: Pi
      implicit none
      class(pdsolver), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP) :: tau,k_n,eta_n,k_w,eta_w,d_eff_w,m1
      real(WP), dimension(3) :: r1,v1,floc,r2
      real(WP), dimension(3), parameter :: vzero=[0.0_WP,0.0_WP,0.0_WP]
      integer :: i,k,j
      if (this%contact_dist.le.0.0_WP.or.dt.le.0.0_WP) return
      if (this%e_n.le.0.0_WP.or.this%e_w.le.0.0_WP) return
      if (this%tau_col.gt.0.0_WP) then
         tau=this%tau_col
      else
         tau=5.0_WP*dt
      end if
      k_n  =(Pi**2+log(this%e_n)**2)/tau**2
      eta_n=-2.0_WP*log(this%e_n)/tau
      k_w  =(Pi**2+log(this%e_w)**2)/tau**2
      eta_w=-2.0_WP*log(this%e_w)/tau
      d_eff_w=0.5_WP*this%contact_dist
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         m1=this%rho*this%vol(i)
         r1=this%y(:,i); v1=this%v(:,i)
         floc=0.0_WP
         ! Wall collisions on faces flagged as walls (virtual partner on the
         ! wall directly normal to the node)
         if (this%lo_bc(1).eq.1) then; r2=[this%dom_lo(1),r1(2),r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,vzero); end if
         if (this%hi_bc(1).eq.1) then; r2=[this%dom_hi(1),r1(2),r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,vzero); end if
         if (this%lo_bc(2).eq.1) then; r2=[r1(1),this%dom_lo(2),r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,vzero); end if
         if (this%hi_bc(2).eq.1) then; r2=[r1(1),this%dom_hi(2),r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,vzero); end if
         if (this%lo_bc(3).eq.1) then; r2=[r1(1),r1(2),this%dom_lo(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,vzero); end if
         if (this%hi_bc(3).eq.1) then; r2=[r1(1),r1(2),this%dom_hi(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,vzero); end if
         ! Particle-particle via the candidate CSR
         do k=this%cptr(i),this%cptr(i+1)-1
            j=this%clst(k)
            call apply_col(k_n,eta_n,this%contact_dist,0.5_WP*m1,this%y(:,j),this%v(:,j))
         end do
         ! Accumulate as force/volume (matches bond force units)
         this%f(:,i)=this%f(:,i)+floc/this%vol(i)
      end do

   contains

      !> Soft-sphere normal force from virtual partner (r2_in, v2_in) onto i.
      !> Host-associated r1, v1, dt, floc.
      subroutine apply_col(kk,ee,d_eff,m_eff,r2_in,v2_in)
         real(WP), intent(in) :: kk,ee,d_eff,m_eff
         real(WP), dimension(3), intent(in) :: r2_in,v2_in
         real(WP) :: d12,rnv,r_influ,delta_n
         real(WP), dimension(3) :: n12,v12,f_n
         d12=norm2(r2_in-r1)
         if (d12.lt.10.0_WP*epsilon(d12)) return   ! self-overlap guard
         n12=(r2_in-r1)/d12
         v12=v1-v2_in
         rnv=dot_product(v12,n12)
         r_influ=min(abs(rnv)*dt,0.2_WP*d_eff)
         delta_n=min(d_eff+r_influ-d12,this%clip_col*d_eff)
         if (delta_n.le.0.0_WP) return
         f_n=(-m_eff*kk*delta_n-m_eff*ee*rnv)*n12
         floc=floc+f_n
      end subroutine apply_col

   end subroutine contact_narrow

   !> Binding CFL: elastic wave + scaled convective (limits 0.5 / 0.1)
   subroutine get_cfl(this,dt,cfl)
      use parallel, only: comm,MPI_REAL_WP
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      implicit none
      class(pdsolver), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      real(WP), parameter :: CFL_scale_conv=5.0_WP
      real(WP) :: K_bulk,mu_shear,c_p,dp_inv,vmax
      integer :: i,ierr
      K_bulk  =this%elastic_modulus/(3.0_WP*(1.0_WP-2.0_WP*this%poisson_ratio))
      mu_shear=this%elastic_modulus/(2.0_WP*(1.0_WP+this%poisson_ratio))
      c_p     =sqrt((K_bulk+4.0_WP*mu_shear/3.0_WP)/this%rho)
      dp_inv  =1.0_WP/this%dV**(1.0_WP/3.0_WP)
      this%CFLe=c_p*dp_inv*dt
      this%CFLp=0.0_WP
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         vmax=max(abs(this%v(1,i)),abs(this%v(2,i)),abs(this%v(3,i)))
         this%CFLp=max(this%CFLp,vmax*dp_inv)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      this%CFLp=this%CFLp*dt
      cfl=max(CFL_scale_conv*this%CFLp,this%CFLe)
   end subroutine get_cfl

   !> Global counts, velocity max, and timer reduction (+reset). Collective.
   subroutine get_info(this)
      use parallel, only: comm,MPI_REAL_WP
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM,MPI_IN_PLACE,MPI_INTEGER8
      implicit none
      class(pdsolver), intent(inout) :: this
      integer :: i,ierr
      integer(I8) :: np_loc
      np_loc=0_I8
      this%Umax=0.0_WP
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         np_loc=np_loc+1_I8
         this%Umax=max(this%Umax,abs(this%v(1,i)),abs(this%v(2,i)),abs(this%v(3,i)))
      end do
      this%np=np_loc
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      ! Broken half-entry census (each broken bond counts twice, except
      ! self-image bonds which have a single half-entry)
      count_broken: block
         integer(I8) :: nb_loc
         integer :: e,i2,j2
         ! Half-entry count (internal) and EXACT broken-bond census: each bond
         ! is counted at exactly one of its two half-entries -- the one whose
         ! node gid is lower (ties = self-image bonds, counted at the
         ! positive-offset image so each appears once)
         nb_loc=0_I8; this%nb_broken=0_I8
         do i2=1,this%nown
            do e=this%ptr(i2),this%ptr(i2+1)-1
               if (this%dmg(e).eq.0_1) cycle
               nb_loc=nb_loc+1_I8
               j2=this%lst(e)
               if (this%gid(i2).lt.this%gid(j2)) then
                  this%nb_broken=this%nb_broken+1_I8
               else if (this%gid(i2).eq.this%gid(j2)) then
                  if (j2.gt.this%nown) then
                     if (shift_positive(this%halo%shift(:,j2-this%nown))) this%nb_broken=this%nb_broken+1_I8
                  end if
               end if
            end do
         end do
         this%nbroken=nb_loc
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nbroken,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nb_broken,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
      end block count_broken
      ! Timers: max (and min for the compute-heavy phases) across ranks, then reset
      call MPI_ALLREDUCE(this%wt_kick,  this%wtmax_kick,  1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(this%wt_halo,  this%wtmax_halo,  1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(this%wt_dil,   this%wtmax_dil,   1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(this%wt_dil,   this%wtmin_dil,   1,MPI_REAL_WP,MPI_MIN,comm,ierr)
      call MPI_ALLREDUCE(this%wt_force, this%wtmax_force, 1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(this%wt_force, this%wtmin_force, 1,MPI_REAL_WP,MPI_MIN,comm,ierr)
      call MPI_ALLREDUCE(this%wt_reduce,this%wtmax_reduce,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(this%wt_contact,this%wtmax_contact,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(this%wt_broad,  this%wtmax_broad,  1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      this%wt_kick=0.0_WP; this%wt_halo=0.0_WP; this%wt_dil=0.0_WP; this%wt_force=0.0_WP; this%wt_reduce=0.0_WP
      this%wt_contact=0.0_WP; this%wt_broad=0.0_WP
      ! Contact-service size census (visibility into the fragmentation-driven
      ! degradation mode of the static graph partition)
      contact_census: block
         integer(I8) :: tmp
         this%nchalo_glob=int(this%nchalo,I8)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nchalo_glob,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
         tmp=0_I8
         if (allocated(this%cptr)) tmp=int(this%cptr(this%nown+1)-1,I8)
         this%ncand_glob=tmp
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%ncand_glob,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
      end block contact_census
   end subroutine get_info

   !> Release all storage
   subroutine finalize(this)
      implicit none
      class(pdsolver), intent(inout) :: this
      if (allocated(this%gid))   deallocate(this%gid)
      if (allocated(this%x0))    deallocate(this%x0)
      if (allocated(this%y))     deallocate(this%y)
      if (allocated(this%v))     deallocate(this%v)
      if (allocated(this%f))     deallocate(this%f)
      if (allocated(this%ff))    deallocate(this%ff)
      if (allocated(this%vol))   deallocate(this%vol)
      if (allocated(this%mw))    deallocate(this%mw)
      if (allocated(this%theta)) deallocate(this%theta)
      if (allocated(this%damage))deallocate(this%damage)
      if (allocated(this%alive)) deallocate(this%alive)
      if (allocated(this%flag))  deallocate(this%flag)
      if (allocated(this%ptr))   deallocate(this%ptr)
      if (allocated(this%lst))   deallocate(this%lst)
      if (allocated(this%dmg))   deallocate(this%dmg)
      if (allocated(this%e_v))   deallocate(this%e_v)
      if (allocated(this%td2))   deallocate(this%td2)
      if (allocated(this%td2a))  deallocate(this%td2a)
      if (allocated(this%cptr))  deallocate(this%cptr)
      if (allocated(this%clst))  deallocate(this%clst)
      if (allocated(this%ylast)) deallocate(this%ylast)
      call this%ohash%finalize()
      call this%dir%finalize()
      call this%halo%finalize()
      call this%chalo%finalize()
      this%nown=0; this%nhalo=0; this%ntot=0; this%nchalo=0
   end subroutine finalize


   !> Mirror synchronization (the coupling bridge). Collective, once per FLUID
   !> step. The caller walks its face particles (AMReX container on the
   !> fluid decomposition) and passes per particle: gid, core owner rank (read
   !> from the face particle's repurposed flag tag), and the F_fluid it interpolated
   !> from the grid. This routine routes F_fluid to the owning nodes (held in
   !> ff across the subsequent PD subcycles) and replies with each node's
   !> current (pos, vel, damage, alive), returned aligned with the caller's
   !> input order for direct write-back into the face particles.
   subroutine exchange(this,nm,mgid,mowner,mff,mpos,mvel,mdmg,malive)
      use parallel, only: comm,nproc,MPI_REAL_WP
      use messager, only: die
      use mpi_f08
      implicit none
      class(pdsolver), intent(inout) :: this
      integer, intent(in) :: nm
      integer(I8), intent(in) :: mgid(:)
      integer, intent(in) :: mowner(:)
      real(WP), intent(in) :: mff(:,:)
      real(WP), intent(out) :: mpos(:,:),mvel(:,:)
      real(WP), intent(out) :: mdmg(:),malive(:)
      integer, dimension(0:nproc-1) :: sc,rc,sd,rd
      integer, dimension(0:nproc-1) :: scw,rcw,sdw,rdw
      integer, allocatable :: pos(:),qpos(:)
      integer(I8), allocatable :: sg(:),rg(:)
      real(WP), allocatable :: sff(:,:),rff(:,:),srep(:,:),rrep(:,:)
      integer :: i,r,nr,idx,ierr

      ! Count and pack by owner, remembering each entry's packed slot
      sc=0
      do i=1,nm
         sc(mowner(i))=sc(mowner(i))+1
      end do
      call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
      sd(0)=0; rd(0)=0
      do r=1,nproc-1
         sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
      end do
      allocate(pos(0:nproc-1),qpos(max(nm,1)))
      allocate(sg(max(nm,1)),sff(3,max(nm,1)))
      pos=sd
      do i=1,nm
         r=mowner(i); pos(r)=pos(r)+1
         sg(pos(r))=mgid(i); sff(:,pos(r))=mff(:,i); qpos(i)=pos(r)
      end do
      nr=sum(rc)
      allocate(rg(max(nr,1)),rff(3,max(nr,1)))
      call MPI_ALLTOALLV(sg,sc,sd,MPI_INTEGER8,rg,rc,rd,MPI_INTEGER8,comm,ierr)
      scw=3*sc; sdw=3*sd; rcw=3*rc; rdw=3*rd
      call MPI_ALLTOALLV(sff,scw,sdw,MPI_REAL_WP,rff,rcw,rdw,MPI_REAL_WP,comm,ierr)

      ! Owner side: ingest F_fluid, build the state reply in arrival order
      allocate(rrep(8,max(nr,1)))
      do i=1,nr
         idx=this%ohash%lookup(rg(i))
         if (idx.lt.1) call die('[pdsolver exchange] face gid not owned by tagged rank')
         this%ff(:,idx)=rff(:,i)
         rrep(1:3,i)=this%y(:,idx)
         rrep(4:6,i)=this%v(:,idx)
         rrep(7,i)  =this%damage(idx)
         rrep(8,i)  =this%alive(idx)
      end do

      ! Reply along the reverse route; unpack to the caller's original order
      allocate(srep(8,max(nm,1)))
      scw=8*rc; sdw=8*rd; rcw=8*sc; rdw=8*sd
      call MPI_ALLTOALLV(rrep,scw,sdw,MPI_REAL_WP,srep,rcw,rdw,MPI_REAL_WP,comm,ierr)
      do i=1,nm
         mpos(:,i) =srep(1:3,qpos(i))
         mvel(:,i) =srep(4:6,qpos(i))
         mdmg(i)   =srep(7,qpos(i))
         malive(i) =srep(8,qpos(i))
      end do
      deallocate(pos,qpos,sg,sff,rg,rff,rrep,srep)
   end subroutine exchange


   !> Owner-rank lookup for arbitrary node gids via the persistent directory.
   !> Collective. Drivers use it to re-stamp face routing tags after restart.
   subroutine query_owners(this,n,gids,owners)
      implicit none
      class(pdsolver), intent(inout) :: this
      integer, intent(in) :: n
      integer(I8), intent(in) :: gids(:)
      integer, intent(out) :: owners(:)
      call this%dir%query(n,gids,owners)
   end subroutine query_owners

   !> Checkpoint the core under <dirname>/: per-rank stream files + root
   !> header. Records are GID-SPACE (no local indices, no partition info) --
   !> nodes: (gid, flag, x0, y, v, f, vol, damage, td2); half-entries:
   !> (node_gid, nbr_gid, image_key, dmg, e_v), the image key reconstructed
   !> from the halo slot's shift. Rank-count portable on read.
   subroutine write_state(this,dirname)
      use parallel, only: rank,nproc,amRoot
      use messager, only: die
      use string,   only: str_medium
      implicit none
      class(pdsolver), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      character(len=str_medium) :: fname
      integer :: iunit,ios,i,e,j,nhe
      integer, allocatable :: hkey(:)
      integer(I8), allocatable :: hnode(:),hnbr(:)
      ! Half-entries in gid space
      nhe=this%ptr(this%nown+1)-1
      allocate(hnode(max(nhe,1)),hnbr(max(nhe,1)),hkey(max(nhe,1)))
      do i=1,this%nown
         do e=this%ptr(i),this%ptr(i+1)-1
            j=this%lst(e)
            hnode(e)=this%gid(i)
            hnbr(e) =this%gid(j)
            if (j.le.this%nown) then
               hkey(e)=PDHALO_KEY0
            else
               hkey(e)=key_of_shift(this%halo%shift(:,j-this%nown),this%Ldom)
            end if
         end do
      end do
      ! Per-rank stream file
      ! All solid state lives under <dirname>/pd/ (root creates it)
      make_dir: block
         use parallel, only: comm
         use mpi_f08,  only: MPI_BARRIER
         integer :: ierr2
         if (amRoot) call execute_command_line('mkdir -p '//trim(dirname)//'/pd')
         call MPI_BARRIER(comm,ierr2)
      end block make_dir
      write(fname,'(a,"/pd/pd_",i7.7,".bin")') trim(dirname),rank
      open(newunit=iunit,file=trim(fname),form='unformatted',access='stream',status='replace',iostat=ios)
      if (ios.ne.0) call die('[pdsolver write_state] cannot open '//trim(fname))
      write(iunit) this%nown,nhe
      write(iunit) this%gid(1:this%nown)
      write(iunit) this%flag(1:this%nown)
      write(iunit) this%x0(:,1:this%nown)
      write(iunit) this%y(:,1:this%nown)
      write(iunit) this%v(:,1:this%nown)
      write(iunit) this%f(:,1:this%nown)
      write(iunit) this%vol(1:this%nown)
      write(iunit) this%damage(1:this%nown)
      write(iunit) this%td2(1:this%nown)
      write(iunit) hnode(1:nhe)
      write(iunit) hnbr(1:nhe)
      write(iunit) hkey(1:nhe)
      write(iunit) this%dmg(1:nhe)
      write(iunit) this%e_v(1:nhe)
      close(iunit)
      deallocate(hnode,hnbr,hkey)
      ! Root header (file count for portable round-robin reads)
      if (amRoot) then
         open(newunit=iunit,file=trim(dirname)//'/pd/header',form='formatted',status='replace',iostat=ios)
         if (ios.ne.0) call die('[pdsolver write_state] cannot open header')
         write(iunit,'(a)') 'pdsolver checkpoint v1'
         write(iunit,'(i0)') nproc
         close(iunit)
      end if
   contains
      !> Reconstruct the packed image key from a slot shift vector
      pure function key_of_shift(s,L) result(k)
         implicit none
         real(WP), dimension(3), intent(in) :: s,L
         integer :: k,n1,n2,n3
         n1=0; n2=0; n3=0
         if (L(1).gt.0.0_WP) n1=nint(s(1)/L(1))
         if (L(2).gt.0.0_WP) n2=nint(s(2)/L(2))
         if (L(3).gt.0.0_WP) n3=nint(s(3)/L(3))
         k=(n1+128)+(n2+128)*256+(n3+128)*65536
      end function key_of_shift
   end subroutine write_state

   !> Restore the core from a checkpoint written by write_state. Collective;
   !> rank-count portable: files read round-robin, nodes re-partitioned by
   !> Morton order of the reference configuration, half-entries routed to
   !> their owners, CSR/halo rebuilt via assemble with the loaded per-entry
   !> state. The caller must configure the solver (initialize + material/
   !> contact/plastic component assignments) BEFORE calling this.
   subroutine read_state(this,dirname)
      use parallel, only: comm,rank,nproc,MPI_REAL_WP
      use messager, only: die
      use string,   only: str_medium
      use mpi_f08
      implicit none
      class(pdsolver), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      character(len=str_medium) :: fname,line
      integer :: nfiles,iunit,ios,f,i,r,ierr
      integer :: nn,nhe,nf,nhf
      integer(I8), allocatable :: gid(:),hnode(:),hnbr(:)
      integer, allocatable :: flag(:),hkey(:),owner(:)
      real(WP), allocatable :: x0(:,:),yy(:,:),vv(:,:),ffb(:,:),vol(:),dmgn(:),td2n(:)
      real(WP), allocatable :: hev(:)
      integer(1), allocatable :: hdmg(:)

      ! Header: number of files written
      nfiles=0
      if (rank.eq.0) then
         open(newunit=iunit,file=trim(dirname)//'/pd/header',form='formatted',status='old',iostat=ios)
         if (ios.ne.0) call die('[pdsolver read_state] no pd/header under '//trim(dirname))
         read(iunit,'(a)') line
         read(iunit,*) nfiles
         close(iunit)
      end if
      call MPI_BCAST(nfiles,1,MPI_INTEGER,0,comm,ierr)

      ! Read my round-robin share of the files, concatenating records
      nn=0; nhe=0
      do f=rank,nfiles-1,nproc
         write(fname,'(a,"/pd/pd_",i7.7,".bin")') trim(dirname),f
         open(newunit=iunit,file=trim(fname),form='unformatted',access='stream',status='old',iostat=ios)
         if (ios.ne.0) call die('[pdsolver read_state] cannot open '//trim(fname))
         read(iunit) nf,nhf
         call grow_i8(gid,nn,nf);   call grow_i4(flag,nn,nf)
         call grow_r2(x0,nn,nf);    call grow_r2(yy,nn,nf)
         call grow_r2(vv,nn,nf);    call grow_r2(ffb,nn,nf)
         call grow_r1(vol,nn,nf);   call grow_r1(dmgn,nn,nf); call grow_r1(td2n,nn,nf)
         read(iunit) gid(nn+1:nn+nf)
         read(iunit) flag(nn+1:nn+nf)
         read(iunit) x0(:,nn+1:nn+nf)
         read(iunit) yy(:,nn+1:nn+nf)
         read(iunit) vv(:,nn+1:nn+nf)
         read(iunit) ffb(:,nn+1:nn+nf)
         read(iunit) vol(nn+1:nn+nf)
         read(iunit) dmgn(nn+1:nn+nf)
         read(iunit) td2n(nn+1:nn+nf)
         call grow_i8(hnode,nhe,nhf); call grow_i8(hnbr,nhe,nhf)
         call grow_i4(hkey,nhe,nhf);  call grow_i1(hdmg,nhe,nhf); call grow_r1(hev,nhe,nhf)
         read(iunit) hnode(nhe+1:nhe+nhf)
         read(iunit) hnbr(nhe+1:nhe+nhf)
         read(iunit) hkey(nhe+1:nhe+nhf)
         read(iunit) hdmg(nhe+1:nhe+nhf)
         read(iunit) hev(nhe+1:nhe+nhf)
         close(iunit)
         nn=nn+nf; nhe=nhe+nhf
      end do
      if (.not.allocated(gid)) then   ! ranks with no files still join collectives
         allocate(gid(1),flag(1),x0(3,1),yy(3,1),vv(3,1),ffb(3,1),vol(1),dmgn(1),td2n(1))
         allocate(hnode(1),hnbr(1),hkey(1),hdmg(1),hev(1))
      end if

      ! Re-partition nodes by Morton order of the REFERENCE configuration and
      ! route the full records (pd_partition routes the set_nodes payload; the
      ! remaining fields ride a second, identically-ordered exchange)
      repartition: block
         integer(I8), allocatable :: rgid(:)
         real(WP), allocatable :: rx0(:,:),rvv(:,:),rvol(:),extra(:,:)
         integer, allocatable :: rflag(:)
         integer, dimension(0:nproc-1) :: sc,rc,sd,rd,scw,rcw,sdw,rdw
         integer, allocatable :: pos(:)
         integer :: nr
         allocate(owner(max(nn,1)))
         call pd_partition(nn,gid,x0,vv,flag,vol,owner,nr,rgid,rx0,rvv,rflag,rvol)
         ! Second exchange: (y, f, damage, td2) = 8 reals, packed in the same
         ! per-destination input order as pd_partition's own packing
         sc=0
         do i=1,nn
            sc(owner(i))=sc(owner(i))+1
         end do
         call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
         sd(0)=0; rd(0)=0
         do r=1,nproc-1
            sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
         end do
         allocate(pos(0:nproc-1),extra(8,max(nn,1)),this%rextra_tmp(8,max(nr,1)))
         pos=sd
         do i=1,nn
            r=owner(i); pos(r)=pos(r)+1
            extra(1:3,pos(r))=yy(:,i)
            extra(4:6,pos(r))=ffb(:,i)
            extra(7,pos(r))  =dmgn(i)
            extra(8,pos(r))  =td2n(i)
         end do
         scw=8*sc; sdw=8*sd; rcw=8*rc; rdw=8*rd
         call MPI_ALLTOALLV(extra,scw,sdw,MPI_REAL_WP,this%rextra_tmp,rcw,rdw,MPI_REAL_WP,comm,ierr)
         ! Load the routed nodes, then overlay the restart-only fields
         call this%set_nodes(nr,rgid,rx0,rvv,rflag,rvol)
         do i=1,nr
            this%y(:,i)   =this%rextra_tmp(1:3,i)
            this%f(:,i)   =this%rextra_tmp(4:6,i)
            this%damage(i)=this%rextra_tmp(7,i)
         end do
         deallocate(pos,extra,rgid,rx0,rvv,rflag,rvol)
      end block repartition

      ! Register the directory over the new partition, route half-entries to
      ! their owners (state travels along), and rebuild CSR/halo/reference
      route_and_assemble: block
         integer(I8), allocatable :: rnode(:),rnbr(:)
         integer, allocatable :: rkey(:),howner(:)
         real(WP), allocatable :: rev(:)
         integer(1), allocatable :: rdmg(:)
         integer, dimension(0:nproc-1) :: sc,rc,sd,rd
         integer, allocatable :: pos(:)
         integer(I8), allocatable :: s8(:)
         integer, allocatable :: s4(:)
         real(WP), allocatable :: sr(:)
         integer(1), allocatable :: s1(:)
         integer :: rn,h
         call this%dir%finalize()
         call this%dir%register(this%nown,this%gid(1:this%nown))
         allocate(howner(max(nhe,1)))
         call this%dir%query(nhe,hnode,howner)
         sc=0
         do i=1,nhe
            sc(howner(i))=sc(howner(i))+1
         end do
         call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
         sd(0)=0; rd(0)=0
         do r=1,nproc-1
            sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
         end do
         rn=sum(rc)
         allocate(rnode(max(rn,1)),rnbr(max(rn,1)),rkey(max(rn,1)),rev(max(rn,1)),rdmg(max(rn,1)))
         allocate(pos(0:nproc-1),s8(max(nhe,1)),s4(max(nhe,1)),sr(max(nhe,1)),s1(max(nhe,1)))
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s8(pos(h))=hnode(i)
         end do
         call MPI_ALLTOALLV(s8,sc,sd,MPI_INTEGER8,rnode,rc,rd,MPI_INTEGER8,comm,ierr)
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s8(pos(h))=hnbr(i)
         end do
         call MPI_ALLTOALLV(s8,sc,sd,MPI_INTEGER8,rnbr,rc,rd,MPI_INTEGER8,comm,ierr)
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s4(pos(h))=hkey(i)
         end do
         call MPI_ALLTOALLV(s4,sc,sd,MPI_INTEGER,rkey,rc,rd,MPI_INTEGER,comm,ierr)
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; sr(pos(h))=hev(i)
         end do
         call MPI_ALLTOALLV(sr,sc,sd,MPI_REAL_WP,rev,rc,rd,MPI_REAL_WP,comm,ierr)
         pos=sd
         do i=1,nhe
            h=howner(i); pos(h)=pos(h)+1; s1(pos(h))=hdmg(i)
         end do
         call MPI_ALLTOALLV(s1,sc,sd,MPI_INTEGER1,rdmg,rc,rd,MPI_INTEGER1,comm,ierr)
         call this%assemble(rn,rnode,rnbr,rkey,rev,rdmg)
         deallocate(rnode,rnbr,rkey,rev,rdmg,pos,s8,s4,sr,s1,howner)
      end block route_and_assemble

      ! Overlay td2 (assemble allocates it zeroed) and life status
      do i=1,this%nown
         this%td2(i)=this%rextra_tmp(8,i)
         if (this%flag(i).eq.PDC_IS_DEAD) this%alive(i)=0.0_WP
      end do
      deallocate(this%rextra_tmp)
      call this%halo%update1(this%alive)
      deallocate(gid,flag,x0,yy,vv,ffb,vol,dmgn,td2n,hnode,hnbr,hkey,hdmg,hev,owner)

   contains

      subroutine grow_i8(a,n,add)
         integer(I8), allocatable, intent(inout) :: a(:)
         integer, intent(in) :: n,add
         integer(I8), allocatable :: t(:)
         allocate(t(n+add)); if (n.gt.0) t(1:n)=a(1:n)
         call move_alloc(t,a)
      end subroutine grow_i8
      subroutine grow_i4(a,n,add)
         integer, allocatable, intent(inout) :: a(:)
         integer, intent(in) :: n,add
         integer, allocatable :: t(:)
         allocate(t(n+add)); if (n.gt.0) t(1:n)=a(1:n)
         call move_alloc(t,a)
      end subroutine grow_i4
      subroutine grow_i1(a,n,add)
         integer(1), allocatable, intent(inout) :: a(:)
         integer, intent(in) :: n,add
         integer(1), allocatable :: t(:)
         allocate(t(n+add)); if (n.gt.0) t(1:n)=a(1:n)
         call move_alloc(t,a)
      end subroutine grow_i1
      subroutine grow_r1(a,n,add)
         real(WP), allocatable, intent(inout) :: a(:)
         integer, intent(in) :: n,add
         real(WP), allocatable :: t(:)
         allocate(t(n+add)); if (n.gt.0) t(1:n)=a(1:n)
         call move_alloc(t,a)
      end subroutine grow_r1
      subroutine grow_r2(a,n,add)
         real(WP), allocatable, intent(inout) :: a(:,:)
         integer, intent(in) :: n,add
         real(WP), allocatable :: t(:,:)
         allocate(t(3,n+add)); if (n.gt.0) t(:,1:n)=a(:,1:n)
         call move_alloc(t,a)
      end subroutine grow_r2

   end subroutine read_state


   !> Static load-balancing partition of the reference configuration.
   !> Collective; called once at handoff, BEFORE set_nodes. Nodes are ordered
   !> by the Morton key of their reference position and split into equal-count
   !> contiguous ranges: on a uniform lattice family size is ~constant, so node
   !> count ~ bond work (a family-weighted split can substitute later), and
   !> bond work is motion-invariant -- this balance holds for the entire run
   !> regardless of deformation or flight, using ALL ranks even when the solid
   !> occupies a corner of the fluid domain.
   !> Inputs: this rank's extracted nodes (any distribution). Outputs: the
   !> nodes assigned to this rank, plus each INPUT node's assigned owner (for
   !> stamping the face particles' routing tags).
   subroutine pd_partition(n_in,gid_in,pos_in,vel_in,flag_in,vol_in,owner_out, &
   &                       n_out,gid_out,pos_out,vel_out,flag_out,vol_out)
      use parallel, only: comm,rank,nproc,amRoot,MPI_REAL_WP
      use pdhalo_class, only: sort3_perm
      use mpi_f08
      implicit none
      integer, intent(in) :: n_in
      integer(I8), intent(in) :: gid_in(:)
      real(WP), intent(in) :: pos_in(:,:),vel_in(:,:)
      integer, intent(in) :: flag_in(:)
      real(WP), intent(in) :: vol_in(:)
      integer, intent(out) :: owner_out(:)
      integer, intent(out) :: n_out
      integer(I8), allocatable, intent(out) :: gid_out(:)
      real(WP), allocatable, intent(out) :: pos_out(:,:),vel_out(:,:),vol_out(:)
      integer, allocatable, intent(out) :: flag_out(:)
      real(WP), dimension(3) :: blo,bhi,inv
      integer(I8), allocatable :: keys(:),splitters(:)
      integer :: i,r,d,ierr

      ! Global reference bounds
      blo=huge(1.0_WP); bhi=-huge(1.0_WP)
      do i=1,n_in
         blo=min(blo,pos_in(:,i)); bhi=max(bhi,pos_in(:,i))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,blo,3,MPI_REAL_WP,MPI_MIN,comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,bhi,3,MPI_REAL_WP,MPI_MAX,comm,ierr)
      do d=1,3
         inv(d)=0.0_WP
         if (bhi(d).gt.blo(d)) inv(d)=2097151.0_WP/(bhi(d)-blo(d))
      end do

      ! Morton keys of this rank's nodes
      allocate(keys(max(n_in,1)))
      do i=1,n_in
         keys(i)=morton(pos_in(:,i),blo,inv)
      end do

      ! Equal-weight splitters by SAMPLE SORT: each rank contributes a few
      ! evenly-spaced samples of its locally sorted keys, weighted by its node
      ! count -- root memory is O(nproc*S), not O(N_global), so this scales to
      ! very large rank counts (the old gather-all-keys approach walled at
      ! root memory and int32 N_global).
      allocate(splitters(max(nproc-1,1)))
      sample_splitters: block
         integer, parameter :: S=16
         integer(I8), allocatable :: lsamp(:),gsamp(:),gw(:)
         real(WP), allocatable :: w(:)
         integer, allocatable :: perm(:),zk(:),scnt(:),sdis(:)
         integer(I8) :: wtot,wcum,wtarg
         integer :: ns,j,r2,gtot
         ! Locally sort keys (permutation) and draw samples
         allocate(perm(max(n_in,1)),zk(max(n_in,1)))
         do i=1,n_in
            perm(i)=i
         end do
         zk=0
         if (n_in.gt.1) call sort3_perm(zk,keys(1:n_in),zk,perm,1,n_in)
         ns=min(S,n_in)
         allocate(lsamp(max(ns,1)))
         do j=1,ns
            lsamp(j)=keys(perm(min(n_in,int((real(j,WP)-0.5_WP)*real(n_in,WP)/real(ns,WP))+1)))
         end do
         deallocate(perm,zk)
         ! Gather samples (+ per-rank sample counts and node counts) on root
         allocate(scnt(nproc),sdis(nproc))
         call MPI_GATHER(ns,1,MPI_INTEGER,scnt,1,MPI_INTEGER,0,comm,ierr)
         gtot=0
         if (amRoot) then
            sdis(1)=0
            do r2=2,nproc
               sdis(r2)=sdis(r2-1)+scnt(r2-1)
            end do
            gtot=sum(scnt)
         end if
         allocate(gsamp(max(gtot,1)),gw(nproc))
         call MPI_GATHERV(lsamp,ns,MPI_INTEGER8,gsamp,scnt,sdis,MPI_INTEGER8,0,comm,ierr)
         call MPI_GATHER(int(n_in,I8),1,MPI_INTEGER8,gw,1,MPI_INTEGER8,0,comm,ierr)
         if (amRoot.and.gtot.gt.0) then
            root_split: block
               integer, allocatable :: p2(:),z2(:)
               real(WP), allocatable :: sw(:)
               ! Weight each sample by (its rank's node count)/(its rank's samples)
               allocate(sw(gtot),p2(gtot),z2(gtot))
               do r2=1,nproc
                  do j=sdis(r2)+1,sdis(r2)+scnt(r2)
                     sw(j)=real(gw(r2),WP)/real(max(scnt(r2),1),WP)
                  end do
               end do
               do j=1,gtot
                  p2(j)=j
               end do
               z2=0
               call sort3_perm(z2,gsamp(1:gtot),z2,p2,1,gtot)
               ! Single cumulative-weight pass placing all nproc-1 splitters
               wtot=sum(gw)
               wcum=0_I8; r2=1
               do j=1,gtot
                  if (r2.gt.nproc-1) exit
                  wcum=wcum+int(sw(p2(j)),I8)
                  do while (r2.le.nproc-1.and.wcum.ge.(int(r2,I8)*wtot)/int(nproc,I8))
                     splitters(r2)=gsamp(p2(j))
                     r2=r2+1
                  end do
               end do
               do while (r2.le.nproc-1)
                  splitters(r2)=huge(1_I8)   ! degenerate tail: empty upper buckets
                  r2=r2+1
               end do
               deallocate(sw,p2,z2)
            end block root_split
         end if
         deallocate(lsamp,gsamp,gw,scnt,sdis)
      end block sample_splitters
      if (nproc.gt.1) call MPI_BCAST(splitters,nproc-1,MPI_INTEGER8,0,comm,ierr)

      ! Assign owners: bucket = number of splitters <= key
      do i=1,n_in
         owner_out(i)=0
         do r=1,nproc-1
            if (keys(i).ge.splitters(r)) owner_out(i)=r
         end do
      end do
      deallocate(keys,splitters)

      ! Route node payloads to their owners
      route_nodes: block
         integer, dimension(0:nproc-1) :: sc,rc,sd,rd,scw,rcw,sdw,rdw
         integer, allocatable :: pos(:),sflag(:)
         integer(I8), allocatable :: sgid(:)
         real(WP), allocatable :: sdat(:,:),rdat(:,:)
         sc=0
         do i=1,n_in
            sc(owner_out(i))=sc(owner_out(i))+1
         end do
         call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
         sd(0)=0; rd(0)=0
         do r=1,nproc-1
            sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
         end do
         n_out=sum(rc)
         allocate(pos(0:nproc-1),sgid(max(n_in,1)),sflag(max(n_in,1)),sdat(7,max(n_in,1)))
         pos=sd
         do i=1,n_in
            r=owner_out(i); pos(r)=pos(r)+1
            sgid(pos(r))=gid_in(i)
            sflag(pos(r))=flag_in(i)
            sdat(1:3,pos(r))=pos_in(:,i)
            sdat(4:6,pos(r))=vel_in(:,i)
            sdat(7,pos(r))  =vol_in(i)
         end do
         allocate(gid_out(max(n_out,1)),flag_out(max(n_out,1)),rdat(7,max(n_out,1)))
         allocate(pos_out(3,max(n_out,1)),vel_out(3,max(n_out,1)),vol_out(max(n_out,1)))
         call MPI_ALLTOALLV(sgid,sc,sd,MPI_INTEGER8,gid_out,rc,rd,MPI_INTEGER8,comm,ierr)
         call MPI_ALLTOALLV(sflag,sc,sd,MPI_INTEGER,flag_out,rc,rd,MPI_INTEGER,comm,ierr)
         scw=7*sc; sdw=7*sd; rcw=7*rc; rdw=7*rd
         call MPI_ALLTOALLV(sdat,scw,sdw,MPI_REAL_WP,rdat,rcw,rdw,MPI_REAL_WP,comm,ierr)
         do i=1,n_out
            pos_out(:,i)=rdat(1:3,i)
            vel_out(:,i)=rdat(4:6,i)
            vol_out(i)  =rdat(7,i)
         end do
         deallocate(pos,sgid,sflag,sdat,rdat)
      end block route_nodes

   contains

      !> 63-bit Morton key: 21 bits per dimension, bit-interleaved
      pure function morton(p,lo,inv) result(key)
         implicit none
         real(WP), dimension(3), intent(in) :: p,lo,inv
         integer(I8) :: key
         integer(I8), dimension(3) :: ix
         integer :: b,d
         do d=1,3
            ix(d)=int(min(max((p(d)-lo(d))*inv(d),0.0_WP),2097151.0_WP),I8)
         end do
         key=0_I8
         do b=0,20
            do d=1,3
               if (btest(ix(d),b)) key=ibset(key,3*b+d-1)
            end do
         end do
      end function morton

   end subroutine pd_partition


   !> Bin geometry: cell size >= the search radius (so +/-1 cell sweeps are complete),
   !> dims clamped to keep total cell count bounded on huge/degenerate extents
   subroutine setup_bins(lo,hi,r,gl,h,nc)
      real(WP), dimension(3), intent(in) :: lo,hi
      real(WP), intent(in) :: r
      real(WP), dimension(3), intent(out) :: gl,h
      integer, dimension(3), intent(out) :: nc
      integer :: d
      do d=1,3
      gl(d)=lo(d)-0.5_WP*r
      nc(d)=max(1,min(256,int((hi(d)-lo(d)+r)/r)))
      h(d)=max((hi(d)+0.5_WP*r-gl(d))/real(nc(d),WP),r)
      end do
   end subroutine setup_bins

   !> Flattened cell index of a position (clamped into the grid)
   pure function cell_of(p,gl,h,nc) result(k)
      real(WP), dimension(3), intent(in) :: p,gl,h
      integer, dimension(3), intent(in) :: nc
      integer :: k,c(3),d
      do d=1,3
      c(d)=min(nc(d),max(1,int((p(d)-gl(d))/h(d))+1))
      end do
      k=c(1)+nc(1)*(c(2)-1)+nc(1)*nc(2)*(c(3)-1)
   end function cell_of


   !> Influence function (constant, Peridigm default).
   pure function omega(d,h) result(w)
      implicit none
      real(WP), intent(in) :: d,h
      real(WP) :: w
      w=1.0_WP
   end function omega

   !> Lexicographic sign of an image shift: .true. for the "positive" member
   !> of a self-image pair (first nonzero component positive), so each
   !> self-image bond is census-counted exactly once.
   pure function shift_positive(s) result(p)
      implicit none
      real(WP), dimension(3), intent(in) :: s
      logical :: p
      integer :: d
      p=.false.
      do d=1,3
         if (abs(s(d)).gt.0.0_WP) then
            p=(s(d).gt.0.0_WP)
            return
         end if
      end do
   end function shift_positive

   !> Negate a packed periodic image offset (amrpd hist1 convention)
   pure function negkey(key) result(nk)
      implicit none
      integer, intent(in) :: key
      integer :: nk,n1,n2,n3
      n1=mod(key,256)-128; n2=mod(key/256,256)-128; n3=key/65536-128
      nk=(-n1+128)+(-n2+128)*256+(-n3+128)*65536
   end function negkey

end module pdsolver_class
