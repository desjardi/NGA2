!> Persistent graph-halo communication for the peridynamics solver (pdsolver).
!>
!> Two objects:
!>   pddir  -- distributed GID directory. Owner-rank resolution for arbitrary
!>             global ids via a hashed home-rank protocol (Fibonacci-mixed:
!>             raw mod collapses on structured idcpu keys).
!>             Built once at init, used during plan construction, then discarded.
!>   pdhalo -- persistent halo exchange plan. A halo SLOT is a (gid, image-offset)
!>             pair: a node bonded to two periodic images of the same partner
!>             gets two slots with different shifts. Shifts are applied at
!>             unpack time on the receiver, so send buffers are pure copies and
!>             the same owned node can serve any number of slots/images.
!>             Two operations per substep:
!>               update(field) -- owner values -> halo slots (positions get +shift)
!>               reduce(field) -- halo-slot accumulations -> add back into owners
!>             Both are nonblocking isend/irecv with fixed, deterministic
!>             pack/unpack order (neighbor rank ascending, slot order within).
!>
!> Self-rank "neighbors" (periodic self-images or same-rank image bonds) are
!> handled uniformly through MPI self-messages -- no special-case code path.
module pdhalo_class
   use precision, only: WP,I8
   use mpi_f08
   implicit none
   private

   public :: pddir,pdhalo,sort3_perm,PDHALO_KEY0

   !> Packed zero image offset ((0+128) + (0+128)*256 + (0+128)*65536),
   !> matching amrpd's hist1 convention.
   integer, parameter :: PDHALO_KEY0=8421504

   !> Distributed GID directory (hashed home-rank protocol)
   type :: pddir
      integer :: n=0                          !< number of gids homed on this rank
      integer(I8), allocatable :: keys(:)     !< gids homed on this rank (sorted)
      integer,     allocatable :: owner(:)    !< owner rank per homed gid (aligned with keys)
   contains
      procedure :: register
      procedure :: query
      procedure :: finalize => dir_finalize
   end type pddir

   !> Persistent halo plan + exchange buffers
   type :: pdhalo
      integer :: nown=0                       !< owned nodes (halo slots are indexed nown+1..nown+nhalo)
      integer :: nhalo=0                      !< halo slot count
      ! Receive side: whom I receive halo data from (= owners of my slots)
      integer :: nrecv=0
      integer, allocatable :: nbr_recv(:)     !< source ranks, ascending
      integer, allocatable :: recv_ptr(:)     !< (nrecv+1) slot group offsets
      ! Send side: whom I send owned data to (= ranks holding slots of my nodes)
      integer :: nsend=0
      integer, allocatable :: nbr_send(:)     !< destination ranks, ascending
      integer, allocatable :: send_ptr(:)     !< (nsend+1) entry group offsets
      integer, allocatable :: send_idx(:)     !< owned node index per send entry (duplicates allowed: one per remote slot)
      ! Per-slot image shift (added to position components at unpack)
      real(WP), allocatable :: shift(:,:)     !< (3,nhalo)
      ! Persistent message buffers (grown on demand)
      real(WP), allocatable :: sbuf(:),rbuf(:)
   contains
      procedure :: build
      procedure :: update
      procedure :: update1
      procedure :: reduce
      procedure :: finalize => halo_finalize
   end type pdhalo

contains


   ! ===========================================================================
   ! Sorting utility: recursive quicksort of a permutation over a triple key
   ! (a int, g int64, k int), ordered lexicographically. Used for deterministic
   ! halo-slot and CSR ordering. a is typically an owner rank or a node index.
   ! ===========================================================================
   recursive subroutine sort3_perm(a,g,k,perm,lo,hi)
      implicit none
      integer,     intent(in)    :: a(:)
      integer(I8), intent(in)    :: g(:)
      integer,     intent(in)    :: k(:)
      integer,     intent(inout) :: perm(:)
      integer,     intent(in)    :: lo,hi
      integer :: i,j,tp,pv
      if (lo.ge.hi) return
      pv=perm((lo+hi)/2)
      i=lo; j=hi
      do
         do while (less3(perm(i),pv)); i=i+1; end do
         do while (less3(pv,perm(j))); j=j-1; end do
         if (i.le.j) then
            tp=perm(i); perm(i)=perm(j); perm(j)=tp
            i=i+1; j=j-1
         end if
         if (i.gt.j) exit
      end do
      call sort3_perm(a,g,k,perm,lo,j)
      call sort3_perm(a,g,k,perm,i,hi)
   contains
      logical function less3(p,q)
         integer, intent(in) :: p,q
         if (a(p).ne.a(q)) then
            less3=a(p).lt.a(q)
         else if (g(p).ne.g(q)) then
            less3=g(p).lt.g(q)
         else
            less3=k(p).lt.k(q)
         end if
      end function less3
   end subroutine sort3_perm


   ! ===========================================================================
   ! PDDIR -- distributed GID directory
   ! ===========================================================================

   !> Register this rank's owned gids with their home ranks. Collective.
   subroutine register(this,n,gids)
      use parallel, only: comm,rank,nproc
      use pdhash_class, only: gid_hash
      implicit none
      class(pddir), intent(inout) :: this
      integer, intent(in) :: n
      integer(I8), intent(in) :: gids(:)
      integer, dimension(0:nproc-1) :: sc,rc,sd,rd
      integer(I8), allocatable :: sg(:),rg(:)
      integer, allocatable :: pos(:)
      integer :: i,h,nr,r,ierr
      ! Count per home rank
      sc=0
      do i=1,n
         h=home(gids(i)); sc(h)=sc(h)+1
      end do
      call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
      sd(0)=0; rd(0)=0
      do r=1,nproc-1
         sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
      end do
      ! Pack and exchange gids
      allocate(sg(max(n,1)),pos(0:nproc-1))
      pos=sd
      do i=1,n
         h=home(gids(i)); pos(h)=pos(h)+1; sg(pos(h))=gids(i)
      end do
      nr=sum(rc)
      allocate(rg(max(nr,1)))
      call MPI_ALLTOALLV(sg,sc,sd,MPI_INTEGER8,rg,rc,rd,MPI_INTEGER8,comm,ierr)
      deallocate(sg,pos)
      ! Store: owner of each received gid = the rank it arrived from
      this%n=nr
      allocate(this%keys(max(nr,1)),this%owner(max(nr,1)))
      this%keys(1:nr)=rg(1:nr)
      do r=0,nproc-1
         do i=rd(r)+1,rd(r)+rc(r)
            this%owner(i)=r
         end do
      end do
      ! Sort keys with the owner array following (simple perm sort)
      sort_dir: block
         integer, allocatable :: perm(:),zk(:),ow(:)
         integer(I8), allocatable :: kk(:)
         integer :: m
         m=nr
         if (m.gt.0) then
            allocate(perm(m),zk(m),ow(m),kk(m))
            do i=1,m
               perm(i)=i
            end do
            zk=0
            call sort3_perm(zk,this%keys(1:m),zk,perm,1,m)
            kk=this%keys(1:m); ow=this%owner(1:m)
            do i=1,m
               this%keys(i) =kk(perm(i))
               this%owner(i)=ow(perm(i))
            end do
            deallocate(perm,zk,ow,kk)
         end if
      end block sort_dir
      deallocate(rg)
   end subroutine register

   !> Resolve owner ranks for m gids. Collective. Dies on unknown gid.
   subroutine query(this,m,gids,owners)
      use parallel, only: comm,nproc
      use messager, only: die
      implicit none
      class(pddir), intent(in) :: this
      integer, intent(in) :: m
      integer(I8), intent(in) :: gids(:)
      integer, intent(out) :: owners(:)
      integer, dimension(0:nproc-1) :: sc,rc,sd,rd
      integer(I8), allocatable :: sg(:),rg(:)
      integer, allocatable :: pos(:),qpos(:),rans(:),reply(:)
      integer :: i,h,nr,r,idx,ierr
      ! Count and pack queries by home rank; remember each query's packed slot
      sc=0
      do i=1,m
         h=home(gids(i)); sc(h)=sc(h)+1
      end do
      call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
      sd(0)=0; rd(0)=0
      do r=1,nproc-1
         sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
      end do
      allocate(sg(max(m,1)),pos(0:nproc-1),qpos(max(m,1)))
      pos=sd
      do i=1,m
         h=home(gids(i)); pos(h)=pos(h)+1; sg(pos(h))=gids(i); qpos(i)=pos(h)
      end do
      nr=sum(rc)
      allocate(rg(max(nr,1)))
      call MPI_ALLTOALLV(sg,sc,sd,MPI_INTEGER8,rg,rc,rd,MPI_INTEGER8,comm,ierr)
      ! Answer each received query by binary search of the sorted directory
      allocate(rans(max(nr,1)))
      do i=1,nr
         idx=dir_lookup(this,rg(i))
         if (idx.lt.1) call die('[pddir query] gid not found in directory')
         rans(i)=this%owner(idx)
      end do
      ! Send answers back along the reverse route (counts swapped)
      allocate(reply(max(m,1)))
      call MPI_ALLTOALLV(rans,rc,rd,MPI_INTEGER,reply,sc,sd,MPI_INTEGER,comm,ierr)
      do i=1,m
         owners(i)=reply(qpos(i))
      end do
      deallocate(sg,rg,pos,qpos,rans,reply)
   end subroutine query

   !> Binary search of the sorted directory keys. Returns index or -1.
   pure function dir_lookup(this,key) result(idx)
      implicit none
      class(pddir), intent(in) :: this
      integer(I8), intent(in) :: key
      integer :: idx,lo,hi,mid
      idx=-1
      if (.not.allocated(this%keys).or.this%n.eq.0) return
      lo=1; hi=this%n
      do while (lo.le.hi)
         mid=(lo+hi)/2
         if (this%keys(mid).lt.key) then
            lo=mid+1
         else if (this%keys(mid).gt.key) then
            hi=mid-1
         else
            idx=mid
            return
         end if
      end do
   end function dir_lookup

   !> Release directory storage
   subroutine dir_finalize(this)
      implicit none
      class(pddir), intent(inout) :: this
      if (allocated(this%keys))  deallocate(this%keys)
      if (allocated(this%owner)) deallocate(this%owner)
      this%n=0
   end subroutine dir_finalize

   !> Home rank of a gid. Keys are STRUCTURED (AMReX idcpu = id<<24|cpu: raw
   !> mod collapses onto few ranks -- all of them rank 0 for power-of-two
   !> nproc when cpu=0), so mix the bits first (Fibonacci hash; the multiply
   !> wraps by design, and the logical shift keeps the result nonnegative).
   pure function home(gid) result(h)
      use parallel, only: nproc
      implicit none
      integer(I8), intent(in) :: gid
      integer :: h
      integer(I8) :: k
      k=gid*(-7046029254386353131_I8)
      h=int(mod(ishft(k,-40),int(nproc,I8)))
   end function home


   ! ===========================================================================
   ! PDHALO -- persistent halo plan
   ! ===========================================================================

   !> Build the halo plan. Collective.
   !>   nown     : owned node count (slots index from nown+1)
   !>   ohash    : gid->owned-index hash over this rank's owned gids
   !>   nreq     : number of UNIQUE remote references (gid, image-key) pairs
   !>   rgid/rkey: the references (key packs the image offset, amrpd hist1 style)
   !>   rowner   : owner rank of each reference's gid (from pddir%query)
   !>   Ldom/per : domain lengths and periodicity (for shift vectors)
   !>   slot     : OUT -- final halo slot (1..nhalo) of each input reference
   subroutine build(this,nown,ohash,nreq,rgid,rkey,rowner,Ldom,per,slot)
      use parallel, only: comm,nproc
      use messager, only: die
      use pdhash_class, only: gid_hash
      implicit none
      class(pdhalo), intent(inout) :: this
      integer, intent(in) :: nown,nreq
      type(gid_hash), intent(in) :: ohash
      integer(I8), intent(in) :: rgid(:)
      integer, intent(in) :: rkey(:),rowner(:)
      real(WP), intent(in) :: Ldom(3)
      logical, intent(in) :: per(3)
      integer, intent(out) :: slot(:)
      integer, allocatable :: perm(:)
      integer, dimension(0:nproc-1) :: sc,rc,sd,rd
      integer :: i,s,r,n1,n2,n3,ierr,nr,lid
      integer(I8), allocatable :: sg(:),rg(:)

      this%nown=nown
      this%nhalo=nreq

      ! Deterministic slot order: sort references by (owner, gid, key)
      allocate(perm(max(nreq,1)))
      do i=1,nreq
         perm(i)=i
      end do
      if (nreq.gt.1) call sort3_perm(rowner,rgid,rkey,perm,1,nreq)
      do s=1,nreq
         slot(perm(s))=s
      end do

      ! Receive groups (one per distinct owner, ascending by construction)
      count_recv: block
         integer :: prev
         this%nrecv=0; prev=-1
         do s=1,nreq
            if (rowner(perm(s)).ne.prev) then
               this%nrecv=this%nrecv+1; prev=rowner(perm(s))
            end if
         end do
         allocate(this%nbr_recv(max(this%nrecv,1)),this%recv_ptr(this%nrecv+1))
         this%nrecv=0; prev=-1
         do s=1,nreq
            if (rowner(perm(s)).ne.prev) then
               this%nrecv=this%nrecv+1; prev=rowner(perm(s))
               this%nbr_recv(this%nrecv)=prev
               this%recv_ptr(this%nrecv)=s
            end if
         end do
         this%recv_ptr(this%nrecv+1)=nreq+1
      end block count_recv

      ! Per-slot shift vectors from the packed image key
      allocate(this%shift(3,max(nreq,1)))
      do s=1,nreq
         i=perm(s)
         n1=mod(rkey(i),256)-128; n2=mod(rkey(i)/256,256)-128; n3=rkey(i)/65536-128
         if ((n1.ne.0.and..not.per(1)).or.(n2.ne.0.and..not.per(2)).or.(n3.ne.0.and..not.per(3))) &
         &  call die('[pdhalo build] nonzero image offset along a non-periodic direction')
         this%shift(1,s)=real(n1,WP)*Ldom(1)
         this%shift(2,s)=real(n2,WP)*Ldom(2)
         this%shift(3,s)=real(n3,WP)*Ldom(3)
      end do

      ! Tell every owner which of its nodes we need (gids in slot order).
      ! Payload order within each destination = our slot order, and MPI
      ! preserves per-pair message order, so the owner's send list built in
      ! arrival order matches our slot order exactly.
      sc=0
      do s=1,nreq
         sc(rowner(perm(s)))=sc(rowner(perm(s)))+1
      end do
      call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
      sd(0)=0; rd(0)=0
      do r=1,nproc-1
         sd(r)=sd(r-1)+sc(r-1); rd(r)=rd(r-1)+rc(r-1)
      end do
      allocate(sg(max(nreq,1)))
      do s=1,nreq
         sg(s)=rgid(perm(s))     ! grouped by owner because slots are owner-sorted
      end do
      nr=sum(rc)
      allocate(rg(max(nr,1)))
      call MPI_ALLTOALLV(sg,sc,sd,MPI_INTEGER8,rg,rc,rd,MPI_INTEGER8,comm,ierr)

      ! Send groups: ranks that requested nodes from me
      count_send: block
         integer :: g
         this%nsend=count(rc.gt.0)
         allocate(this%nbr_send(max(this%nsend,1)),this%send_ptr(this%nsend+1))
         allocate(this%send_idx(max(nr,1)))
         g=0; this%send_ptr(1)=1
         do r=0,nproc-1
            if (rc(r).gt.0) then
               g=g+1
               this%nbr_send(g)=r
               this%send_ptr(g+1)=this%send_ptr(g)+rc(r)
               do i=rd(r)+1,rd(r)+rc(r)
                  lid=ohash%lookup(rg(i))
                  if (lid.lt.1) call die('[pdhalo build] halo request for a gid this rank does not own')
                  this%send_idx(this%send_ptr(g)+(i-rd(r)-1))=lid
               end do
            end if
         end do
      end block count_send

      deallocate(perm,sg,rg)
   end subroutine build

   !> Refresh halo slots with current owner values: field(:,1:nown) -> slots.
   !> field is (ncomp, nown+nhalo). If shifted, per-slot image shifts are added
   !> to components 1:3 (positions). Deterministic unpack order.
   subroutine update(this,field,ncomp,shifted)
      use parallel, only: comm,MPI_REAL_WP
      use messager, only: die
      implicit none
      class(pdhalo), intent(inout) :: this
      real(WP), intent(inout) :: field(:,:)
      integer, intent(in) :: ncomp
      logical, intent(in) :: shifted
      type(MPI_Request), allocatable :: reqs(:)
      integer :: i,g,s,off,cnt,nrq,ierr
      integer :: nsend_tot
      if (shifted.and.ncomp.lt.3) call die('[pdhalo update] shifted update requires ncomp>=3')
      nsend_tot=this%send_ptr(this%nsend+1)-1
      call ensure_buffers(this,ncomp*max(nsend_tot,1),ncomp*max(this%nhalo,1))
      allocate(reqs(this%nrecv+this%nsend))
      nrq=0
      ! Post receives (one message per source rank)
      do g=1,this%nrecv
         off=ncomp*(this%recv_ptr(g)-1)
         cnt=ncomp*(this%recv_ptr(g+1)-this%recv_ptr(g))
         nrq=nrq+1
         call MPI_IRECV(this%rbuf(off+1:off+cnt),cnt,MPI_REAL_WP,this%nbr_recv(g),101,comm,reqs(nrq),ierr)
      end do
      ! Pack and send (one message per destination rank)
      do g=1,this%nsend
         off=ncomp*(this%send_ptr(g)-1)
         do i=this%send_ptr(g),this%send_ptr(g+1)-1
            this%sbuf(off+ncomp*(i-this%send_ptr(g))+1:off+ncomp*(i-this%send_ptr(g))+ncomp)=field(1:ncomp,this%send_idx(i))
         end do
         cnt=ncomp*(this%send_ptr(g+1)-this%send_ptr(g))
         nrq=nrq+1
         call MPI_ISEND(this%sbuf(off+1:off+cnt),cnt,MPI_REAL_WP,this%nbr_send(g),101,comm,reqs(nrq),ierr)
      end do
      call MPI_WAITALL(nrq,reqs,MPI_STATUSES_IGNORE,ierr)
      ! Unpack into halo slots (slot s lives at field index nown+s)
      do s=1,this%nhalo
         field(1:ncomp,this%nown+s)=this%rbuf(ncomp*(s-1)+1:ncomp*(s-1)+ncomp)
      end do
      if (shifted) then
         do s=1,this%nhalo
            field(1:3,this%nown+s)=field(1:3,this%nown+s)+this%shift(1:3,s)
         end do
      end if
      deallocate(reqs)
   end subroutine update

   !> Scalar-field variant of update (no shift): owner values -> halo slots.
   !> Used for static per-node scalars (e.g., nodal volume) filled once at init.
   subroutine update1(this,field)
      use parallel, only: comm,MPI_REAL_WP
      implicit none
      class(pdhalo), intent(inout) :: this
      real(WP), intent(inout) :: field(:)
      type(MPI_Request), allocatable :: reqs(:)
      integer :: i,g,s,off,cnt,nrq,ierr
      integer :: nsend_tot
      nsend_tot=this%send_ptr(this%nsend+1)-1
      call ensure_buffers(this,max(nsend_tot,1),max(this%nhalo,1))
      allocate(reqs(this%nrecv+this%nsend))
      nrq=0
      do g=1,this%nrecv
         off=this%recv_ptr(g)-1
         cnt=this%recv_ptr(g+1)-this%recv_ptr(g)
         nrq=nrq+1
         call MPI_IRECV(this%rbuf(off+1:off+cnt),cnt,MPI_REAL_WP,this%nbr_recv(g),103,comm,reqs(nrq),ierr)
      end do
      do g=1,this%nsend
         off=this%send_ptr(g)-1
         do i=this%send_ptr(g),this%send_ptr(g+1)-1
            this%sbuf(i)=field(this%send_idx(i))
         end do
         cnt=this%send_ptr(g+1)-this%send_ptr(g)
         nrq=nrq+1
         call MPI_ISEND(this%sbuf(off+1:off+cnt),cnt,MPI_REAL_WP,this%nbr_send(g),103,comm,reqs(nrq),ierr)
      end do
      call MPI_WAITALL(nrq,reqs,MPI_STATUSES_IGNORE,ierr)
      do s=1,this%nhalo
         field(this%nown+s)=this%rbuf(s)
      end do
      deallocate(reqs)
   end subroutine update1

   !> Add halo-slot accumulations back into their owners: slots -> field(:,1:nown).
   !> Reverse of update: slot data flows to the owner, which adds it into the
   !> owned entries listed in send_idx. Deterministic add order (group order,
   !> then entry order within group).
   subroutine reduce(this,field,ncomp)
      use parallel, only: comm,MPI_REAL_WP
      implicit none
      class(pdhalo), intent(inout) :: this
      real(WP), intent(inout) :: field(:,:)
      integer, intent(in) :: ncomp
      type(MPI_Request), allocatable :: reqs(:)
      integer :: i,g,s,off,cnt,nrq,ierr
      integer :: nsend_tot
      nsend_tot=this%send_ptr(this%nsend+1)-1
      ! Buffers: sending nhalo slots, receiving nsend_tot contributions
      call ensure_buffers(this,ncomp*max(this%nhalo,1),ncomp*max(nsend_tot,1))
      allocate(reqs(this%nrecv+this%nsend))
      nrq=0
      ! Post receives along the send-plan links (contributions to my owned nodes)
      do g=1,this%nsend
         off=ncomp*(this%send_ptr(g)-1)
         cnt=ncomp*(this%send_ptr(g+1)-this%send_ptr(g))
         nrq=nrq+1
         call MPI_IRECV(this%rbuf(off+1:off+cnt),cnt,MPI_REAL_WP,this%nbr_send(g),102,comm,reqs(nrq),ierr)
      end do
      ! Pack halo slots and send to their owners along the recv-plan links
      do g=1,this%nrecv
         off=ncomp*(this%recv_ptr(g)-1)
         do s=this%recv_ptr(g),this%recv_ptr(g+1)-1
            this%sbuf(off+ncomp*(s-this%recv_ptr(g))+1:off+ncomp*(s-this%recv_ptr(g))+ncomp)=field(1:ncomp,this%nown+s)
         end do
         cnt=ncomp*(this%recv_ptr(g+1)-this%recv_ptr(g))
         nrq=nrq+1
         call MPI_ISEND(this%sbuf(off+1:off+cnt),cnt,MPI_REAL_WP,this%nbr_recv(g),102,comm,reqs(nrq),ierr)
      end do
      call MPI_WAITALL(nrq,reqs,MPI_STATUSES_IGNORE,ierr)
      ! Accumulate received contributions into owned nodes
      do i=1,nsend_tot
         field(1:ncomp,this%send_idx(i))=field(1:ncomp,this%send_idx(i))+this%rbuf(ncomp*(i-1)+1:ncomp*(i-1)+ncomp)
      end do
      deallocate(reqs)
   end subroutine reduce

   !> Grow persistent buffers on demand
   subroutine ensure_buffers(this,ns,nr)
      implicit none
      class(pdhalo), intent(inout) :: this
      integer, intent(in) :: ns,nr
      if (allocated(this%sbuf)) then
         if (size(this%sbuf).lt.ns) deallocate(this%sbuf)
      end if
      if (.not.allocated(this%sbuf)) allocate(this%sbuf(ns))
      if (allocated(this%rbuf)) then
         if (size(this%rbuf).lt.nr) deallocate(this%rbuf)
      end if
      if (.not.allocated(this%rbuf)) allocate(this%rbuf(nr))
   end subroutine ensure_buffers

   !> Release plan storage
   subroutine halo_finalize(this)
      implicit none
      class(pdhalo), intent(inout) :: this
      if (allocated(this%nbr_recv)) deallocate(this%nbr_recv)
      if (allocated(this%recv_ptr)) deallocate(this%recv_ptr)
      if (allocated(this%nbr_send)) deallocate(this%nbr_send)
      if (allocated(this%send_ptr)) deallocate(this%send_ptr)
      if (allocated(this%send_idx)) deallocate(this%send_idx)
      if (allocated(this%shift))    deallocate(this%shift)
      if (allocated(this%sbuf))     deallocate(this%sbuf)
      if (allocated(this%rbuf))     deallocate(this%rbuf)
      this%nown=0; this%nhalo=0; this%nrecv=0; this%nsend=0
   end subroutine halo_finalize


end module pdhalo_class
