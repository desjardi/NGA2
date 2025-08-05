!> Coupler concept is defined here: it takes in two pgrid
!> objects and builds the communication and interpolation
!> layer to exchange data between them.
module coupler_class
   use precision,      only: WP
   use string,         only: str_medium
   use pgrid_class,    only: pgrid
   use mpi_f08
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: coupler
   
   ! Interpolation data structure
   type :: interp_map
      integer , dimension(:,:), allocatable :: srcind                !< Src indices of dst points that this processor can interpolate
      integer , dimension(:,:), allocatable :: dstind                !< Dst indices of dst points that this processor will receive
      real(WP), dimension(:,:), allocatable :: w                     !< Interpolation weights for dst points that this processor can interpolate
      integer , dimension(:)  , allocatable :: nsend_proc,nsend_disp !< Number of points to send to each processor and corresponding displacement
      integer , dimension(:)  , allocatable :: nrecv_proc,nrecv_disp !< Number of points to receive from each processor and corresponding displacement
      !> Data to send
      real(WP), dimension(:), allocatable :: data_send
      integer :: nsend=0
      !> Data to receive
      real(WP), dimension(:), allocatable :: data_recv
      integer :: nrecv=0
   end type interp_map
   
   !> Coupler object definition
   type :: coupler
      
      ! This is the name of the coupler
      character(len=str_medium) :: name='UNNAMED_CPL'     !< Coupler name (default=UNNAMED_CPL)
      
      ! These are our two pgrids
      type(pgrid), pointer :: src=>NULL()                 !< Source grid
      type(pgrid), pointer :: dst=>NULL()                 !< Destination grid
      
      ! Logicals to help us know if we have received a src or dst grid
      logical :: got_src=.false.                          !< Were we given a src grid
      logical :: got_dst=.false.                          !< Were we given a dst grid
      
      ! This is our communication information
      type(MPI_Comm) :: comm                              !< Intracommunicator over the union of both groups
      type(MPI_Group) :: sgrp,dgrp,grp                    !< Source and destination groups and their union
      integer :: nproc                                    !< Number of processors
      integer :: rank                                     !< Processor grid rank
      logical :: amRoot                                   !< Am I root for the coupler?
      integer :: sroot                                    !< Rank of src grid root on union group
      integer :: droot                                    !< Rank of dst grid root on union group
      
      ! Destination grid partitioning
      integer :: dnproc,dnpx,dnpy,dnpz
      
      ! Interpolation map for each location type
      type(interp_map), pointer :: cmap,xmap,ymap,zmap    !< Interpolation and communication support for cell-center, x-face, y-face, and z-face data
      
      ! Track last push location
      character(len=1) :: last_push_loc=''

      ! Do we reuse buffers
      logical :: reuse_buffers=.true.
      
   contains
      procedure :: set_src                                !< Routine that sets the source grid
      procedure :: set_dst                                !< Routine that sets the destination grid
      procedure :: initialize                             !< Routine that prepares all interpolation metrics from src to dst
      procedure :: push                                   !< Src routine that pushes a field into our send data storage
      procedure :: pull                                   !< Dst routine that pulls a field from our received data storage
      procedure :: transfer                               !< Routine that performs the src->dst data transfer
      procedure :: finalize                               !< Finalize coupler object
   end type coupler
   
   
   !> Declare coupler constructor
   interface coupler
      procedure construct_from_two_groups
   end interface coupler
   
contains
   
   
   !> Coupler constructor from two groups
   function construct_from_two_groups(src_grp,dst_grp,name) result(self)
      use messager, only: die
      use parallel, only: comm
      implicit none
      type(coupler) :: self
      type(MPI_Group), intent(in) :: src_grp,dst_grp
      character(len=*), intent(in) :: name
      integer, dimension(1) :: rankin,rankout
      integer :: ierr
      ! Set name for the coupler
      self%name=trim(adjustl(name))
      ! Build group union
      self%sgrp=src_grp
      self%dgrp=dst_grp
      call MPI_GROUP_UNION(self%sgrp,self%dgrp,self%grp,ierr)
      ! Gather some info for communication
      call MPI_GROUP_SIZE(self%grp,self%nproc,ierr)
      if (self%nproc.eq.0) call die('[coupler constructor] Somehow the union of both groups is of size zero')
      call MPI_GROUP_RANK(self%grp,self%rank ,ierr)
      if (self%rank.eq.MPI_UNDEFINED) call die('[coupler constructor] All processors that call the constructor need to be in one of the two groups')
      ! Create intracommunicator for the new group
      call MPI_COMM_CREATE_GROUP(comm,self%grp,0,self%comm,ierr)
      ! Find roots for both grids on the shared communicator
      rankin=0; call MPI_GROUP_TRANSLATE_RANKS(self%sgrp,1,rankin,self%grp,rankout,ierr); self%sroot=rankout(1)
      rankin=0; call MPI_GROUP_TRANSLATE_RANKS(self%dgrp,1,rankin,self%grp,rankout,ierr); self%droot=rankout(1)
      ! Set coupler root to src root
      self%amRoot=(self%rank.eq.self%sroot)
   end function construct_from_two_groups
   
   
   !> Set the source grid - to be called by processors in src_group
   subroutine set_src(this,pg)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      ! Point to the grid
      this%src=>pg
      ! Set a flag
      this%got_src=.true.
   end subroutine set_src
   
   
   !> Set the destination grid - to be called by processors in dst_group
   subroutine set_dst(this,pg)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      ! Point to the grid
      this%dst=>pg
      ! Set a flag
      this%got_dst=.true.
   end subroutine set_dst
   
   
   !> Prepare interpolation metrics from src to dst
   subroutine initialize(this)
      implicit none
      class(coupler), intent(inout) :: this
      integer, dimension(:,:,:), allocatable :: rankmap
      
      ! Make destination grid available to all
      share_grid: block
         use sgrid_class, only: sgrid
         use parallel,    only: MPI_REAL_WP
         character(len=str_medium) :: simu_name
         real(WP), dimension(:), allocatable :: x
         real(WP), dimension(:), allocatable :: y
         real(WP), dimension(:), allocatable :: z
         logical :: xper,yper,zper
         integer :: no,nx,ny,nz,coord,ierr
         ! Destination root process extracts its own sgrid
         if (this%rank.eq.this%droot) then
            simu_name=this%dst%name
            coord=this%dst%coordsys
            xper=this%dst%xper
            yper=this%dst%yper
            zper=this%dst%zper
            nx=this%dst%nx
            ny=this%dst%ny
            nz=this%dst%nz
            no=this%dst%no
         end if
         ! Then it broadcasts it to our group
         call MPI_BCAST(simu_name,len(simu_name),MPI_CHARACTER,this%droot,this%comm,ierr)
         call MPI_BCAST(coord    ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(xper     ,1             ,MPI_LOGICAL  ,this%droot,this%comm,ierr)
         call MPI_BCAST(yper     ,1             ,MPI_LOGICAL  ,this%droot,this%comm,ierr)
         call MPI_BCAST(zper     ,1             ,MPI_LOGICAL  ,this%droot,this%comm,ierr)
         call MPI_BCAST(nx       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(ny       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(nz       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(no       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         ! Allocate x/y/z, fill it, and bcast
         allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1))
         if (this%rank.eq.this%droot) then
            x(1:nx+1)=this%dst%x(this%dst%imin:this%dst%imax+1)
            y(1:ny+1)=this%dst%y(this%dst%jmin:this%dst%jmax+1)
            z(1:nz+1)=this%dst%z(this%dst%kmin:this%dst%kmax+1)
         end if
         call MPI_BCAST(x,nx+1,MPI_REAL_WP,this%droot,this%comm,ierr)
         call MPI_BCAST(y,ny+1,MPI_REAL_WP,this%droot,this%comm,ierr)
         call MPI_BCAST(z,nz+1,MPI_REAL_WP,this%droot,this%comm,ierr)
         ! Finish creating the sgrid
         if (.not.this%got_dst) then
            allocate(this%dst)
            this%dst%sgrid=sgrid(coord,no,x,y,z,xper,yper,zper,trim(adjustl(simu_name)))
         end if
         ! Deallocate
         deallocate(x,y,z)
      end block share_grid
      
      ! Create rankmap, making destination partition map available to all
      share_partition: block
         integer :: ierr,n
         integer, dimension(:), allocatable :: diproc,djproc,dkproc
         ! Destination root process extracts partition
         if (this%rank.eq.this%droot) then
            this%dnproc=this%dst%nproc
            this%dnpx=this%dst%npx
            this%dnpy=this%dst%npy
            this%dnpz=this%dst%npz
         end if
         ! Broadcast it to our group
         call MPI_BCAST(this%dnpx,  1,MPI_INTEGER,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpy,  1,MPI_INTEGER,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpz,  1,MPI_INTEGER,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnproc,1,MPI_INTEGER,this%droot,this%comm,ierr)
         ! Prepare communication arrays
         allocate(diproc(0:this%nproc-1),djproc(0:this%nproc-1),dkproc(0:this%nproc-1))
         ! Provide a default iproc/jproc/kproc to processors without dst grid
         if (.not.this%got_dst) then
            this%dst%iproc=0
            this%dst%jproc=0
            this%dst%kproc=0
         end if
         ! Allgather the rank->(iproc,jproc,kproc) info
         call MPI_ALLGATHER(this%dst%iproc,1,MPI_INTEGER,diproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%dst%jproc,1,MPI_INTEGER,djproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%dst%kproc,1,MPI_INTEGER,dkproc,1,MPI_INTEGER,this%comm,ierr)
         ! Allocate the destination rankmap
         allocate(rankmap(this%dnpx,this%dnpy,this%dnpz))
         ! Finally, flip the rankmap data
         do n=0,this%nproc-1
            if (diproc(n).gt.0) then
               rankmap(diproc(n),djproc(n),dkproc(n))=n
            end if
         end do
         ! Deallocate communication arrays
         deallocate(diproc,djproc,dkproc)
      end block share_partition
      
      ! Setup interpolation maps
      allocate(this%cmap); call setup_interp_map(this%cmap,'c')
      allocate(this%xmap); call setup_interp_map(this%xmap,'x')
      allocate(this%ymap); call setup_interp_map(this%ymap,'y')
      allocate(this%zmap); call setup_interp_map(this%zmap,'z')
      
      ! We're done with rankmap, deallocate it
      deallocate(rankmap)
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (this%amRoot) then
            write(message,'("Coupler [",a,"] from pgrid [",a,"] to pgrid [",a,"]")') trim(this%name),trim(this%src%name),trim(this%dst%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   contains
      
      !> Prepare interpolation metrics from src to dst at loc
      subroutine setup_interp_map(map,loc)
         implicit none
         class(interp_map), intent(inout) :: map
         character(len=1) , intent(in)    :: loc
         integer, dimension(:)  , allocatable :: rk
         integer, dimension(:,:), allocatable :: dstind
         
         ! First, src processors identify all dst points that belong to them
         find_dst_points: block
            integer :: i,j,k,count,qx,rx,qy,ry,qz,rz
            real(WP), dimension(3) :: pt
            integer , dimension(3) :: coords
            
            ! Initialize counter
            map%nsend=0
            
            ! Only the src processors need to work here
            if (this%got_src) then
               
               ! Traverse the entire dst mesh and count points that can be interpolated
               do k=this%dst%kmino,this%dst%kmaxo
                  do j=this%dst%jmino,this%dst%jmaxo
                     do i=this%dst%imino,this%dst%imaxo
                        ! Choose appropriate location
                        select case (loc)
                        case ('c'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%zm(k)]
                        case ('x'); pt=[this%dst%x (i),this%dst%ym(j),this%dst%zm(k)]
                        case ('y'); pt=[this%dst%xm(i),this%dst%y (j),this%dst%zm(k)]
                        case ('z'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%z (k)]
                        end select
                        ! Skip grid points that lie outside our local domain (slight over-counting here due to .gt.)
                        if (pt(1).lt.this%src%x(this%src%imin_).or.pt(1).gt.this%src%x(this%src%imax_+1).or. &
                        &   pt(2).lt.this%src%y(this%src%jmin_).or.pt(2).gt.this%src%y(this%src%jmax_+1).or. &
                        &   pt(3).lt.this%src%z(this%src%kmin_).or.pt(3).gt.this%src%z(this%src%kmax_+1)) cycle
                        ! Increment our counter
                        map%nsend=map%nsend+1
                     end do
                  end do
               end do
               
               ! Continue only if points where found
               if (map%nsend.gt.0) then
                  
                  ! Allocate storage for both ind, rk, and w
                  allocate(map%srcind(3,map%nsend))
                  allocate(map%w(3,map%nsend))
                  allocate(rk(map%nsend))
                  allocate(dstind(3,map%nsend))
                  
                  ! Get ready to find the dst rank
                  qx=this%dst%nx/this%dnpx; rx=mod(this%dst%nx,this%dnpx)
                  qy=this%dst%ny/this%dnpy; ry=mod(this%dst%ny,this%dnpy)
                  qz=this%dst%nz/this%dnpz; rz=mod(this%dst%nz,this%dnpz)
                  
                  ! Traverse the entire dst mesh and identify points that can be interpolated
                  count=0
                  do k=this%dst%kmino,this%dst%kmaxo
                     do j=this%dst%jmino,this%dst%jmaxo
                        do i=this%dst%imino,this%dst%imaxo
                           ! Choose appropriate location
                           select case (loc)
                           case ('c'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%zm(k)]
                           case ('x'); pt=[this%dst%x (i),this%dst%ym(j),this%dst%zm(k)]
                           case ('y'); pt=[this%dst%xm(i),this%dst%y (j),this%dst%zm(k)]
                           case ('z'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%z (k)]
                           end select
                           ! Skip grid points that lie outside our local domain (slight over-counting here due to .gt.)
                           if (pt(1).lt.this%src%x(this%src%imin_).or.pt(1).gt.this%src%x(this%src%imax_+1).or. &
                           &   pt(2).lt.this%src%y(this%src%jmin_).or.pt(2).gt.this%src%y(this%src%jmax_+1).or. &
                           &   pt(3).lt.this%src%z(this%src%kmin_).or.pt(3).gt.this%src%z(this%src%kmax_+1)) cycle
                           ! The point is in our subdomain, so increment our counter
                           count=count+1
                           ! Locate point and store src index and interpolation weights with proper mesh location
                           select case (loc)
                           case ('c'); call get_weights_and_indices_c(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,map%w(:,count),map%srcind(:,count))
                           case ('x'); call get_weights_and_indices_x(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,map%w(:,count),map%srcind(:,count))
                           case ('y'); call get_weights_and_indices_y(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,map%w(:,count),map%srcind(:,count))
                           case ('z'); call get_weights_and_indices_z(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,map%w(:,count),map%srcind(:,count))
                           end select
                           ! Find coords of the dst processor
                           coords(1)=0; do while (i.ge.this%dst%imin+(coords(1)+1)*qx+min(coords(1)+1,rx).and.coords(1)+1.lt.this%dnpx); coords(1)=coords(1)+1; end do
                           coords(2)=0; do while (j.ge.this%dst%jmin+(coords(2)+1)*qy+min(coords(2)+1,ry).and.coords(2)+1.lt.this%dnpy); coords(2)=coords(2)+1; end do
                           coords(3)=0; do while (k.ge.this%dst%kmin+(coords(3)+1)*qz+min(coords(3)+1,rz).and.coords(3)+1.lt.this%dnpz); coords(3)=coords(3)+1; end do
                           ! Convert into a rank and store
                           rk(count)=rankmap(coords(1)+1,coords(2)+1,coords(3)+1)
                           ! Also store the dstind
                           dstind(:,count)=[i,j,k]
                        end do
                     end do
                  end do
                  
               end if
               
            end if
            
         end block find_dst_points
         
         ! Next step is to sort our data by recipient
         sort_communication: block
            integer :: n,ierr
            
            ! First brute-force quick-sort our data by dst recipient
            if (map%nsend.gt.0) call qs_commdata(rk,dstind,map%srcind,map%w)
            
            ! Allocate and zero out per processor counters
            allocate(map%nsend_proc(0:this%nproc-1)); map%nsend_proc=0
            allocate(map%nrecv_proc(0:this%nproc-1)); map%nrecv_proc=0
            
            ! Loop through identified points and count
            do n=1,map%nsend
               map%nsend_proc(rk(n))=map%nsend_proc(rk(n))+1
            end do
            
            ! We are done using rk
            if (map%nsend.gt.0) deallocate(rk)
            
            ! Prepare information about who receives what from whom
            do n=0,this%nproc-1
               call MPI_gather(map%nsend_proc(n),1,MPI_INTEGER,map%nrecv_proc,1,MPI_INTEGER,n,this%comm,ierr)
            end do
            
            ! Set size of receive buffer
            map%nrecv=sum(map%nrecv_proc)
            
            ! We need to generate displacements
            allocate(map%nsend_disp(0:this%nproc-1)); map%nsend_disp=0
            allocate(map%nrecv_disp(0:this%nproc-1)); map%nrecv_disp=0
            do n=1,this%nproc-1
               map%nsend_disp(n)=map%nsend_disp(n-1)+map%nsend_proc(n-1)
               map%nrecv_disp(n)=map%nrecv_disp(n-1)+map%nrecv_proc(n-1)
            end do
            
         end block sort_communication
         
         ! Communicate dstind to dst processors so they know what to expect
         share_dstind: block
            integer :: ierr
            integer, dimension(map%nsend) :: send_buffer
            integer, dimension(map%nrecv) :: recv_buffer
            ! Receivers allocate mapind
            if (map%nrecv.gt.0) allocate(map%dstind(3,map%nrecv))
            ! Communicate dstind(1)
            if (map%nsend.gt.0) send_buffer=dstind(1,:)
            call MPI_ALLtoALLv(send_buffer,map%nsend_proc,map%nsend_disp,MPI_INTEGER,recv_buffer,map%nrecv_proc,map%nrecv_disp,MPI_INTEGER,this%comm,ierr)
            if (map%nrecv.gt.0) map%dstind(1,:)=recv_buffer
            ! Communicate dstind(2)
            if (map%nsend.gt.0) send_buffer=dstind(2,:)
            call MPI_ALLtoALLv(send_buffer,map%nsend_proc,map%nsend_disp,MPI_INTEGER,recv_buffer,map%nrecv_proc,map%nrecv_disp,MPI_INTEGER,this%comm,ierr)
            if (map%nrecv.gt.0) map%dstind(2,:)=recv_buffer
            ! Communicate dstind(3)
            if (map%nsend.gt.0) send_buffer=dstind(3,:)
            call MPI_ALLtoALLv(send_buffer,map%nsend_proc,map%nsend_disp,MPI_INTEGER,recv_buffer,map%nrecv_proc,map%nrecv_disp,MPI_INTEGER,this%comm,ierr)
            if (map%nrecv.gt.0) map%dstind(3,:)=recv_buffer
            ! We are done using dstind
            if (map%nsend.gt.0) deallocate(dstind)
         end block share_dstind
         
      end subroutine setup_interp_map
      
      !> Private subroutine that finds weights w for the trilinear interpolation from cell centers
      !> to the provided position pos in the vicinity of cell i0,j0,k0 on pgrid pg
      subroutine get_weights_and_indices_c(pg,pos,i0,j0,k0,w,ind)
         implicit none
         class(pgrid), intent(in) :: pg
         real(WP), dimension(3), intent(in) :: pos
         integer, intent(in) :: i0,j0,k0
         integer :: i,j,k
         real(WP), dimension(3), intent(out) :: w
         integer , dimension(3), intent(out) :: ind
         ! Find right i index
         i=max(min(pg%imaxo_-1,i0),pg%imino_)
         do while (pos(1)-pg%xm(i  ).lt.0.0_WP.and.i  .gt.pg%imino_); i=i-1; end do
         do while (pos(1)-pg%xm(i+1).ge.0.0_WP.and.i+1.lt.pg%imaxo_); i=i+1; end do
         ! Find right j index
         j=max(min(pg%jmaxo_-1,j0),pg%jmino_)
         do while (pos(2)-pg%ym(j  ).lt.0.0_WP.and.j  .gt.pg%jmino_); j=j-1; end do
         do while (pos(2)-pg%ym(j+1).ge.0.0_WP.and.j+1.lt.pg%jmaxo_); j=j+1; end do
         ! Find right k index
         k=max(min(pg%kmaxo_-1,k0),pg%kmino_)
         do while (pos(3)-pg%zm(k  ).lt.0.0_WP.and.k  .gt.pg%kmino_); k=k-1; end do
         do while (pos(3)-pg%zm(k+1).ge.0.0_WP.and.k+1.lt.pg%kmaxo_); k=k+1; end do
         ! Return tri-linear interpolation coefficients
         w(1)=(pos(1)-pg%xm(i))/(pg%xm(i+1)-pg%xm(i))
         w(2)=(pos(2)-pg%ym(j))/(pg%ym(j+1)-pg%ym(j))
         w(3)=(pos(3)-pg%zm(k))/(pg%zm(k+1)-pg%zm(k))
         ! Return the indices too
         ind=[i,j,k]
      end subroutine get_weights_and_indices_c
      
      !> Private subroutine that finds weights w for the trilinear interpolation from x-faces
      !> to the provided position pos in the vicinity of cell i0,j0,k0 on pgrid pg
      subroutine get_weights_and_indices_x(pg,pos,i0,j0,k0,w,ind)
         implicit none
         class(pgrid), intent(in) :: pg
         real(WP), dimension(3), intent(in) :: pos
         integer, intent(in) :: i0,j0,k0
         integer :: i,j,k
         real(WP), dimension(3), intent(out) :: w
         integer , dimension(3), intent(out) :: ind
         ! Find right i index
         i=max(min(pg%imaxo_-1,i0),pg%imino_)
         do while (pos(1)-pg%x (i  ).lt.0.0_WP.and.i  .gt.pg%imino_); i=i-1; end do
         do while (pos(1)-pg%x (i+1).ge.0.0_WP.and.i+1.lt.pg%imaxo_); i=i+1; end do
         ! Find right j index
         j=max(min(pg%jmaxo_-1,j0),pg%jmino_)
         do while (pos(2)-pg%ym(j  ).lt.0.0_WP.and.j  .gt.pg%jmino_); j=j-1; end do
         do while (pos(2)-pg%ym(j+1).ge.0.0_WP.and.j+1.lt.pg%jmaxo_); j=j+1; end do
         ! Find right k index
         k=max(min(pg%kmaxo_-1,k0),pg%kmino_)
         do while (pos(3)-pg%zm(k  ).lt.0.0_WP.and.k  .gt.pg%kmino_); k=k-1; end do
         do while (pos(3)-pg%zm(k+1).ge.0.0_WP.and.k+1.lt.pg%kmaxo_); k=k+1; end do
         ! Return tri-linear interpolation coefficients
         w(1)=(pos(1)-pg%x (i))/(pg%x (i+1)-pg%x (i))
         w(2)=(pos(2)-pg%ym(j))/(pg%ym(j+1)-pg%ym(j))
         w(3)=(pos(3)-pg%zm(k))/(pg%zm(k+1)-pg%zm(k))
         ! Return the indices too
         ind=[i,j,k]
      end subroutine get_weights_and_indices_x
      
      !> Private subroutine that finds weights w for the trilinear interpolation from y-faces
      !> to the provided position pos in the vicinity of cell i0,j0,k0 on pgrid pg
      subroutine get_weights_and_indices_y(pg,pos,i0,j0,k0,w,ind)
         implicit none
         class(pgrid), intent(in) :: pg
         real(WP), dimension(3), intent(in) :: pos
         integer, intent(in) :: i0,j0,k0
         integer :: i,j,k
         real(WP), dimension(3), intent(out) :: w
         integer , dimension(3), intent(out) :: ind
         ! Find right i index
         i=max(min(pg%imaxo_-1,i0),pg%imino_)
         do while (pos(1)-pg%xm(i  ).lt.0.0_WP.and.i  .gt.pg%imino_); i=i-1; end do
         do while (pos(1)-pg%xm(i+1).ge.0.0_WP.and.i+1.lt.pg%imaxo_); i=i+1; end do
         ! Find right j index
         j=max(min(pg%jmaxo_-1,j0),pg%jmino_)
         do while (pos(2)-pg%y (j  ).lt.0.0_WP.and.j  .gt.pg%jmino_); j=j-1; end do
         do while (pos(2)-pg%y (j+1).ge.0.0_WP.and.j+1.lt.pg%jmaxo_); j=j+1; end do
         ! Find right k index
         k=max(min(pg%kmaxo_-1,k0),pg%kmino_)
         do while (pos(3)-pg%zm(k  ).lt.0.0_WP.and.k  .gt.pg%kmino_); k=k-1; end do
         do while (pos(3)-pg%zm(k+1).ge.0.0_WP.and.k+1.lt.pg%kmaxo_); k=k+1; end do
         ! Return tri-linear interpolation coefficients
         w(1)=(pos(1)-pg%xm(i))/(pg%xm(i+1)-pg%xm(i))
         w(2)=(pos(2)-pg%y (j))/(pg%y (j+1)-pg%y (j))
         w(3)=(pos(3)-pg%zm(k))/(pg%zm(k+1)-pg%zm(k))
         ! Return the indices too
         ind=[i,j,k]
      end subroutine get_weights_and_indices_y
      
      !> Private subroutine that finds weights w for the trilinear interpolation from z-faces
      !> to the provided position pos in the vicinity of cell i0,j0,k0 on pgrid pg
      subroutine get_weights_and_indices_z(pg,pos,i0,j0,k0,w,ind)
         implicit none
         class(pgrid), intent(in) :: pg
         real(WP), dimension(3), intent(in) :: pos
         integer, intent(in) :: i0,j0,k0
         integer :: i,j,k
         real(WP), dimension(3), intent(out) :: w
         integer , dimension(3), intent(out) :: ind
         ! Find right i index
         i=max(min(pg%imaxo_-1,i0),pg%imino_)
         do while (pos(1)-pg%xm(i  ).lt.0.0_WP.and.i  .gt.pg%imino_); i=i-1; end do
         do while (pos(1)-pg%xm(i+1).ge.0.0_WP.and.i+1.lt.pg%imaxo_); i=i+1; end do
         ! Find right j index
         j=max(min(pg%jmaxo_-1,j0),pg%jmino_)
         do while (pos(2)-pg%ym(j  ).lt.0.0_WP.and.j  .gt.pg%jmino_); j=j-1; end do
         do while (pos(2)-pg%ym(j+1).ge.0.0_WP.and.j+1.lt.pg%jmaxo_); j=j+1; end do
         ! Find right k index
         k=max(min(pg%kmaxo_-1,k0),pg%kmino_)
         do while (pos(3)-pg%z (k  ).lt.0.0_WP.and.k  .gt.pg%kmino_); k=k-1; end do
         do while (pos(3)-pg%z (k+1).ge.0.0_WP.and.k+1.lt.pg%kmaxo_); k=k+1; end do
         ! Return tri-linear interpolation coefficients
         w(1)=(pos(1)-pg%xm(i))/(pg%xm(i+1)-pg%xm(i))
         w(2)=(pos(2)-pg%ym(j))/(pg%ym(j+1)-pg%ym(j))
         w(3)=(pos(3)-pg%z (k))/(pg%z (k+1)-pg%z (k))
         ! Return the indices too
         ind=[i,j,k]
      end subroutine get_weights_and_indices_z
      
      !> Specialized quicksort driver for our communication data
      recursive subroutine qs_commdata(rk,dstind,srcind,w)
         implicit none
         integer , dimension(:)   :: rk
         integer , dimension(:,:) :: dstind
         integer , dimension(:,:) :: srcind
         real(WP), dimension(:,:) :: w
         integer :: imark
         if (size(rk).gt.1) then
            call qs_partition(rk,dstind,srcind,w,imark)
            call qs_commdata(rk(     :imark-1),dstind(:,     :imark-1),srcind(:,     :imark-1),w(:,     :imark-1))
            call qs_commdata(rk(imark:       ),dstind(:,imark:       ),srcind(:,imark:       ),w(:,imark:       ))
         end if
      end subroutine qs_commdata
      
      !> Specialized quicksort partitioning
      subroutine qs_partition(rk,dstind,srcind,w,marker)
         implicit none
         integer , dimension(:)   :: rk
         integer , dimension(:,:) :: dstind
         integer , dimension(:,:) :: srcind
         real(WP), dimension(:,:) :: w
         integer , intent(out) :: marker
         integer :: i,j,x,itmp
         integer , dimension(3) :: i3tmp
         real(WP), dimension(3) :: d3tmp
         x=rk(1)
         i=0
         j=size(rk)+1
         do
            j=j-1
            do
               if (rk(j).le.x) exit
               j=j-1
            end do
            i=i+1
            do
               if (rk(i).ge.x) exit
               i=i+1
            end do
            if (i.lt.j) then
               itmp =      rk(i);       rk(i)=      rk(j);       rk(j)= itmp  ! Swap rk(i) and rk(j)
               d3tmp=     w(:,i);      w(:,i)=     w(:,j);      w(:,j)=d3tmp  ! Swap w(:,i) and w(:,j)
               i3tmp=dstind(:,i); dstind(:,i)=dstind(:,j); dstind(:,j)=i3tmp  ! Swap dstind(:,i) and dstind(:,j)
               i3tmp=srcind(:,i); srcind(:,i)=srcind(:,j); srcind(:,j)=i3tmp  ! Swap srcind(:,i) and srcind(:,j)
            else if (i.eq.j) then
               marker=i+1
               return
            else
               marker=i
               return
            end if
         end do
      end subroutine qs_partition
      
   end subroutine initialize
   
   
   !> Routine that interpolates src data to the send storage - to be called by processors in src_group
   subroutine push(this,A,loc)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      real(WP), dimension(this%src%imino_:,this%src%jmino_:,this%src%kmino_:), intent(in) :: A !< Needs to be (src%imino_:src%imaxo_,src%jmino_:src%jmaxo_,src%kmino_:src%kmaxo_)
      character(len=1) :: loc
      type(interp_map), pointer :: map
      integer :: n
      ! Handle location
      select case (lowercase(loc))
      case ('c'); map=>this%cmap
      case ('x'); map=>this%xmap
      case ('y'); map=>this%ymap
      case ('z'); map=>this%zmap
      case default
         call die('[coupler push] Improper location value was provided')
      end select
      ! Remember push location
      this%last_push_loc=lowercase(loc)
      ! Allocate send buffer
      if (.not.allocated(map%data_send)) allocate(map%data_send(map%nsend))
      ! Fill buffer
      do n=1,map%nsend
         map%data_send(n)=(       map%w(3,n))*((       map%w(2,n))*((       map%w(1,n))*A(map%srcind(1,n)+1,map%srcind(2,n)+1,map%srcind(3,n)+1)  + &
         &                                                          (1.0_WP-map%w(1,n))*A(map%srcind(1,n)  ,map%srcind(2,n)+1,map%srcind(3,n)+1)) + &
         &                                     (1.0_WP-map%w(2,n))*((       map%w(1,n))*A(map%srcind(1,n)+1,map%srcind(2,n)  ,map%srcind(3,n)+1)  + &
         &                                                          (1.0_WP-map%w(1,n))*A(map%srcind(1,n)  ,map%srcind(2,n)  ,map%srcind(3,n)+1)))+ &
         &                (1.0_WP-map%w(3,n))*((       map%w(2,n))*((       map%w(1,n))*A(map%srcind(1,n)+1,map%srcind(2,n)+1,map%srcind(3,n)  )  + &
         &                                                          (1.0_WP-map%w(1,n))*A(map%srcind(1,n)  ,map%srcind(2,n)+1,map%srcind(3,n)  )) + &
         &                                     (1.0_WP-map%w(2,n))*((       map%w(1,n))*A(map%srcind(1,n)+1,map%srcind(2,n)  ,map%srcind(3,n)  )  + &
         &                                                          (1.0_WP-map%w(1,n))*A(map%srcind(1,n)  ,map%srcind(2,n)  ,map%srcind(3,n)  )))
      end do
   end subroutine push
   
   
   !> Routine that pulls dst data from the receive storage - to be called by processors in dst_group
   !> Free up recv buffer here
   subroutine pull(this,A,loc)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      real(WP), dimension(this%dst%imino_:,this%dst%jmino_:,this%dst%kmino_:), intent(out) :: A !< Needs to be (dst%imino_:dst%imaxo_,dst%jmino_:dst%jmaxo_,dst%kmino_:dst%kmaxo_)
      character(len=1) :: loc
      type(interp_map), pointer :: map
      integer :: n
      ! Check pull location is compatible with push location
      if (this%last_push_loc.ne.lowercase(loc)) call die('[coupler pull] Mismatched push/pull locations: push was '//this%last_push_loc//' but pull is '//lowercase(loc))
      ! Handle location
      select case (lowercase(loc))
      case ('c'); map=>this%cmap
      case ('x'); map=>this%xmap
      case ('y'); map=>this%ymap
      case ('z'); map=>this%zmap
      case default
         call die('[coupler pull] Improper location value was provided')
      end select
      ! Pull the data
      do n=1,map%nrecv
         A(map%dstind(1,n),map%dstind(2,n),map%dstind(3,n))=map%data_recv(n)
      end do
      ! Sync it before returning
      call this%dst%sync(A)
      ! Deallocate recv buffer
      if (.not.this%reuse_buffers.and.allocated(map%data_recv)) deallocate(map%data_recv)
   end subroutine pull
   
   
   !> Routine that transfers the data from src to dst - both src_group and dst_group processors need to call
   !> Allocate recv buffer and deallocate send buffer here
   subroutine transfer(this)
      use parallel, only: MPI_REAL_WP
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      integer :: ierr
      type(interp_map), pointer :: map
      ! Broadcast last_push_loc to all processes
      call MPI_BCAST(this%last_push_loc,1,MPI_CHARACTER,this%sroot,this%comm,ierr)
      ! Handle location
      select case (this%last_push_loc)
      case ('c'); map=>this%cmap
      case ('x'); map=>this%xmap
      case ('y'); map=>this%ymap
      case ('z'); map=>this%zmap
      case default
         call die('[coupler transfer] Improper location value was found')
      end select
      ! Allocate recv buffer (and maybe send buffer too)
      if (.not.allocated(map%data_recv)) allocate(map%data_recv(map%nrecv))
      if (.not.allocated(map%data_send)) allocate(map%data_send(map%nsend))
      ! Transfer data
      call MPI_ALLtoALLv(map%data_send,map%nsend_proc,map%nsend_disp,MPI_REAL_WP,map%data_recv,map%nrecv_proc,map%nrecv_disp,MPI_REAL_WP,this%comm,ierr)
      ! Deallocate send buffer
      if (.not.this%reuse_buffers.and.allocated(map%data_send)) deallocate(map%data_send)
   end subroutine transfer
   
   
   !> Finalize coupler object
   subroutine finalize(this)
      implicit none
      class(coupler), intent(inout) :: this
      integer :: ierr
      ! Finalize interp_map objects
      call finalize_interp_map(this%cmap); if (associated(this%cmap)) deallocate(this%cmap)
      call finalize_interp_map(this%xmap); if (associated(this%xmap)) deallocate(this%xmap)
      call finalize_interp_map(this%ymap); if (associated(this%ymap)) deallocate(this%ymap)
      call finalize_interp_map(this%zmap); if (associated(this%zmap)) deallocate(this%zmap)
      ! Nullify grid pointers
      nullify(this%src)
      nullify(this%dst)
      ! Free only the internally created communicator and group
      if (this%comm.ne.MPI_COMM_NULL ) call MPI_Comm_free(this%comm,ierr)
      if (this%grp .ne.MPI_GROUP_NULL) call MPI_Group_free(this%grp,ierr)
      ! Invalidate external groups
      this%sgrp=MPI_GROUP_NULL
      this%dgrp=MPI_GROUP_NULL
   contains
      !> Finalize interp_map object
      subroutine finalize_interp_map(map)
         type(interp_map), intent(inout) :: map
         if (allocated(map%srcind))     deallocate(map%srcind)
         if (allocated(map%dstind))     deallocate(map%dstind)
         if (allocated(map%w))          deallocate(map%w)
         if (allocated(map%nsend_proc)) deallocate(map%nsend_proc)
         if (allocated(map%nsend_disp)) deallocate(map%nsend_disp)
         if (allocated(map%nrecv_proc)) deallocate(map%nrecv_proc)
         if (allocated(map%nrecv_disp)) deallocate(map%nrecv_disp)
         map%nrecv=0; if (allocated(map%data_recv)) deallocate(map%data_recv)
         map%nsend=0; if (allocated(map%data_send)) deallocate(map%data_send)
      end subroutine finalize_interp_map   
   end subroutine finalize
   
   
end module coupler_class
