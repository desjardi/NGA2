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
   
   !> Coupler object definition
   type :: coupler
      
      ! This is the name of the coupler
      character(len=str_medium) :: name='UNNAMED_CPL'     !< Coupler name (default=UNNAMED_CPL)
      
      ! Overlap visualization
      !real(WP), dimension(:,:,:), allocatable :: overlap  !< Array that identifies overlap in the src and dst, on the dst grid (0=no overlap, 1=overlap)
      
      ! These are our two pgrids
      type(pgrid), pointer :: src=>NULL()                 !< Source grid
      type(pgrid), pointer :: dst=>NULL()                 !< Destination grid
      
      ! Cell-centered or face-centered
      character(len=1) :: sloc='c'                        !< Coupling location for src grid ('c'=cell center, 'x'=x-face, 'y'=y-face, 'z'=z-face)
      character(len=1) :: dloc='c'                        !< Coupling location for dst grid ('c'=cell center, 'x'=x-face, 'y'=y-face, 'z'=z-face)
      
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
      
      ! Rank map for dst grid
      integer :: dnproc,dnpx,dnpy,dnpz                    !< Destination grid partitioning
      integer, dimension(:,:,:), allocatable :: rankmap   !< Processor coordinate to union group rank map
      
      ! Interpolation support
      integer , dimension(:,:), allocatable :: srcind     !< Src indices of dst points that this processor can interpolate
      integer , dimension(:,:), allocatable :: dstind     !< Dst indices of dst points that this processor will receive
      real(WP), dimension(:,:), allocatable :: w          !< Interpolation weights for dst points that this processor can interpolate
      
      ! Communication support
      integer :: nsend                                    !< Total number of dst points that this processor can interpolate and will send out
      integer, dimension(:), allocatable :: nsend_proc    !< Number of points to send to each processor
      integer, dimension(:), allocatable :: nsend_disp    !< Data displacement when sending to each processor
      integer :: nrecv                                    !< Total number of dst points that this processor will receive
      integer, dimension(:), allocatable :: nrecv_proc    !< Number of points to receive from each processor
      integer, dimension(:), allocatable :: nrecv_disp    !< Data displacement when receiving from each processor
      
      ! Data storage
      real(WP), dimension(:), allocatable :: data_send    !< Data to send
      real(WP), dimension(:), allocatable :: data_recv    !< Received data
      
   contains
      procedure :: set_src                                !< Routine that sets the source grid
      procedure :: set_dst                                !< Routine that sets the destination grid
      procedure :: initialize                             !< Routine that prepares all interpolation metrics from src to dst
      procedure :: push                                   !< Src routine that pushes a field into our send data storage
      procedure :: pull                                   !< Dst routine that pulls a field from our received data storage
      procedure :: transfer                               !< Routine that performs the src->dst data transfer
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
   subroutine set_src(this,pg,loc)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      character(len=1), optional :: loc
      ! Point to the grid
      this%src=>pg
      ! Set a flag
      this%got_src=.true.
      ! Process location info
      if (present(loc)) this%sloc=lowercase(loc)
      if (this%sloc.ne.'c'.and.this%sloc.ne.'x'.and.this%sloc.ne.'y'.and.this%sloc.ne.'z') call die('[coupler set_src] Improper location value was provided')
   end subroutine set_src
   
   
   !> Set the destination grid - to be called by processors in dst_group
   subroutine set_dst(this,pg,loc)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      character(len=1), optional :: loc
      ! Point to the grid
      this%dst=>pg
      ! Set a flag
      this%got_dst=.true.
      ! Process location info
      if (present(loc)) this%dloc=lowercase(loc)
      if (this%dloc.ne.'c'.and.this%dloc.ne.'x'.and.this%dloc.ne.'y'.and.this%dloc.ne.'z') call die('[coupler set_dst] Improper location value was provided')
   end subroutine set_dst
   
   
   !> Prepare interpolation metrics from src to dst
   subroutine initialize(this)
      implicit none
      class(coupler), intent(inout) :: this
      integer , dimension(:)  , allocatable :: rk
      integer , dimension(:,:), allocatable :: dstind
      
      ! First step is to make destination grid available to all
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
      
      ! Second step is to make destination partition map available to all
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
         allocate(this%rankmap(this%dnpx,this%dnpy,this%dnpz))
         
         ! Finally, flip the rankmap data
         do n=0,this%nproc-1
            if (diproc(n).gt.0) then
               this%rankmap(diproc(n),djproc(n),dkproc(n))=n
            end if
         end do
         
         ! Deallocate communication arrays
         deallocate(diproc,djproc,dkproc)
         
      end block share_partition
      
      ! Now the src processors identify all dst points that belong to them
      find_dst_points: block
         integer :: i,j,k,count,qx,rx,qy,ry,qz,rz
         real(WP), dimension(3) :: pt
         integer , dimension(3) :: coords
         
         ! Initialize counter
         this%nsend=0
         
         ! Only the src processors need to work here
         if (this%got_src) then
            
            ! Traverse the entire dst mesh and count points that can be interpolated
            do k=this%dst%kmino,this%dst%kmaxo
               do j=this%dst%jmino,this%dst%jmaxo
                  do i=this%dst%imino,this%dst%imaxo
                     ! Choose appropriate location
                     select case (this%dloc)
                     case ('c'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%zm(k)]
                     case ('x'); pt=[this%dst%x (i),this%dst%ym(j),this%dst%zm(k)]
                     case ('y'); pt=[this%dst%xm(i),this%dst%y (j),this%dst%zm(k)]
                     case ('z'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%z (k)]
                     end select
                     ! Skip grid points that lie outside our local domain - allow for slight double-counting here with the .gt. ...
                     if (pt(1).lt.this%src%x(this%src%imin_).or.pt(1).gt.this%src%x(this%src%imax_+1).or. &
                     &   pt(2).lt.this%src%y(this%src%jmin_).or.pt(2).gt.this%src%y(this%src%jmax_+1).or. &
                     &   pt(3).lt.this%src%z(this%src%kmin_).or.pt(3).gt.this%src%z(this%src%kmax_+1)) cycle
                     ! Increment our counter
                     this%nsend=this%nsend+1
                  end do
               end do
            end do
            
            ! Continue only if points where found
            if (this%nsend.gt.0) then
               
               ! Allocate storage for both ind, rk, and w
               allocate(this%srcind(3,this%nsend))
               allocate(this%w(3,this%nsend))
               allocate(rk(this%nsend))
               allocate(dstind(3,this%nsend))
               
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
                        select case (this%dloc)
                        case ('c'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%zm(k)]
                        case ('x'); pt=[this%dst%x (i),this%dst%ym(j),this%dst%zm(k)]
                        case ('y'); pt=[this%dst%xm(i),this%dst%y (j),this%dst%zm(k)]
                        case ('z'); pt=[this%dst%xm(i),this%dst%ym(j),this%dst%z (k)]
                        end select
                        ! Skip grid points that lie outside our local domain - allow for slight double-counting here with the .gt. ...
                        if (pt(1).lt.this%src%x(this%src%imin_).or.pt(1).gt.this%src%x(this%src%imax_+1).or. &
                        &   pt(2).lt.this%src%y(this%src%jmin_).or.pt(2).gt.this%src%y(this%src%jmax_+1).or. &
                        &   pt(3).lt.this%src%z(this%src%kmin_).or.pt(3).gt.this%src%z(this%src%kmax_+1)) cycle
                        ! The point is in our subdomain, so increment our counter
                        count=count+1
                        ! Locate point and store src index and interpolation weights with proper mesh location
                        select case (this%sloc)
                        case ('c'); call get_weights_and_indices_c(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,this%w(:,count),this%srcind(:,count))
                        case ('x'); call get_weights_and_indices_x(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,this%w(:,count),this%srcind(:,count))
                        case ('y'); call get_weights_and_indices_y(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,this%w(:,count),this%srcind(:,count))
                        case ('z'); call get_weights_and_indices_z(this%src,pt,this%src%imin_,this%src%jmin_,this%src%kmin_,this%w(:,count),this%srcind(:,count))
                        end select
                        ! Find coords of the dst processor
                        coords(1)=0; do while (i.ge.this%dst%imin+(coords(1)+1)*qx+min(coords(1)+1,rx).and.coords(1)+1.lt.this%dnpx); coords(1)=coords(1)+1; end do
                        coords(2)=0; do while (j.ge.this%dst%jmin+(coords(2)+1)*qy+min(coords(2)+1,ry).and.coords(2)+1.lt.this%dnpy); coords(2)=coords(2)+1; end do
                        coords(3)=0; do while (k.ge.this%dst%kmin+(coords(3)+1)*qz+min(coords(3)+1,rz).and.coords(3)+1.lt.this%dnpz); coords(3)=coords(3)+1; end do
                        ! Convert into a rank and store
                        rk(count)=this%rankmap(coords(1)+1,coords(2)+1,coords(3)+1)
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
         if (this%nsend.gt.0) call qs_commdata(rk,dstind,this%srcind,this%w)
         
         ! Allocate and zero out per processor counters
         allocate(this%nsend_proc(0:this%nproc-1)); this%nsend_proc=0
         allocate(this%nrecv_proc(0:this%nproc-1)); this%nrecv_proc=0
         
         ! Loop through identified points and count
         do n=1,this%nsend
            this%nsend_proc(rk(n))=this%nsend_proc(rk(n))+1
         end do
         
         ! We are done using rk
         if (this%nsend.gt.0) deallocate(rk)
         
         ! Prepare information about who receives what from whom
         do n=0,this%nproc-1
            call MPI_gather(this%nsend_proc(n),1,MPI_INTEGER,this%nrecv_proc,1,MPI_INTEGER,n,this%comm,ierr)
         end do
         
         ! Set size of receive buffer
         this%nrecv=sum(this%nrecv_proc)
         
         ! We need to generate displacements
         allocate(this%nsend_disp(0:this%nproc-1)); this%nsend_disp=0
         allocate(this%nrecv_disp(0:this%nproc-1)); this%nrecv_disp=0
         do n=1,this%nproc-1
            this%nsend_disp(n)=this%nsend_disp(n-1)+this%nsend_proc(n-1)
            this%nrecv_disp(n)=this%nrecv_disp(n-1)+this%nrecv_proc(n-1)
         end do
         
      end block sort_communication
      
      ! Communicate dstind to dst processors so they know what to expect
      share_dstind: block
         integer :: ierr
         integer, dimension(this%nsend) :: send_buffer
         integer, dimension(this%nrecv) :: recv_buffer
         ! Receivers allocate mapind
         if (this%nrecv.gt.0) allocate(this%dstind(3,this%nrecv))
         ! Communicate dstind(1)
         if (this%nsend.gt.0) send_buffer=dstind(1,:)
         call MPI_ALLtoALLv(send_buffer,this%nsend_proc,this%nsend_disp,MPI_INTEGER,recv_buffer,this%nrecv_proc,this%nrecv_disp,MPI_INTEGER,this%comm,ierr)
         if (this%nrecv.gt.0) this%dstind(1,:)=recv_buffer
         ! Communicate dstind(2)
         if (this%nsend.gt.0) send_buffer=dstind(2,:)
         call MPI_ALLtoALLv(send_buffer,this%nsend_proc,this%nsend_disp,MPI_INTEGER,recv_buffer,this%nrecv_proc,this%nrecv_disp,MPI_INTEGER,this%comm,ierr)
         if (this%nrecv.gt.0) this%dstind(2,:)=recv_buffer
         ! Communicate dstind(3)
         if (this%nsend.gt.0) send_buffer=dstind(3,:)
         call MPI_ALLtoALLv(send_buffer,this%nsend_proc,this%nsend_disp,MPI_INTEGER,recv_buffer,this%nrecv_proc,this%nrecv_disp,MPI_INTEGER,this%comm,ierr)
         if (this%nrecv.gt.0) this%dstind(3,:)=recv_buffer
         ! We are done using dstind
         if (this%nsend.gt.0) deallocate(dstind)
      end block share_dstind
      
      ! For visualization, create coupling field=0 if not overlap was found, 1 if overlap was found
      !viz_overlap: block
      !   integer :: n
      !   if (this%got_dst) then
      !      ! Allocate the array
      !      allocate(this%overlap(this%dst%imino_:this%dst%imaxo_,this%dst%jmino_:this%dst%jmaxo_,this%dst%kmino_:this%dst%kmaxo_)); this%overlap=0.0_WP
      !      ! Fill it up with out dstind info
      !      do n=1,this%nrecv
      !         this%overlap(this%dstind(1,n),this%dstind(2,n),this%dstind(3,n))=1.0_WP
      !      end do
      !   end if
      !end block viz_overlap
      
      ! Final step is to allocate the data storage - this is done on the fly below
      !allocate(this%data_send(this%nsend))
      !allocate(this%data_recv(this%nrecv))
      
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
      
   end subroutine initialize
   
   
   !> Routine that interpolates src data to the send storage - to be called by processors in src_group
   !> Allocate send buffer here
   subroutine push(this,A)
      implicit none
      class(coupler), intent(inout) :: this
      real(WP), dimension(this%src%imino_:,this%src%jmino_:,this%src%kmino_:), intent(in) :: A !< Needs to be (src%imino_:src%imaxo_,src%jmino_:src%jmaxo_,src%kmino_:src%kmaxo_)
      integer :: n
      ! Allocate send buffer
      if (.not.allocated(this%data_send)) allocate(this%data_send(this%nsend))
      ! Fill buffer
      do n=1,this%nsend
         this%data_send(n)=(       this%w(3,n))*((       this%w(2,n))*((       this%w(1,n))*A(this%srcind(1,n)+1,this%srcind(2,n)+1,this%srcind(3,n)+1)  + &
         &                                                             (1.0_WP-this%w(1,n))*A(this%srcind(1,n)  ,this%srcind(2,n)+1,this%srcind(3,n)+1)) + &
         &                                       (1.0_WP-this%w(2,n))*((       this%w(1,n))*A(this%srcind(1,n)+1,this%srcind(2,n)  ,this%srcind(3,n)+1)  + &
         &                                                             (1.0_WP-this%w(1,n))*A(this%srcind(1,n)  ,this%srcind(2,n)  ,this%srcind(3,n)+1)))+ &
         &                 (1.0_WP-this%w(3,n))*((       this%w(2,n))*((       this%w(1,n))*A(this%srcind(1,n)+1,this%srcind(2,n)+1,this%srcind(3,n)  )  + &
         &                                                             (1.0_WP-this%w(1,n))*A(this%srcind(1,n)  ,this%srcind(2,n)+1,this%srcind(3,n)  )) + &
         &                                       (1.0_WP-this%w(2,n))*((       this%w(1,n))*A(this%srcind(1,n)+1,this%srcind(2,n)  ,this%srcind(3,n)  )  + &
         &                                                             (1.0_WP-this%w(1,n))*A(this%srcind(1,n)  ,this%srcind(2,n)  ,this%srcind(3,n)  )))
      end do
   end subroutine push
   
   
   !> Routine that pulls dst data from the receive storage - to be called by processors in dst_group
   !> Free up recv buffer here
   subroutine pull(this,A)
      implicit none
      class(coupler), intent(inout) :: this
      real(WP), dimension(this%dst%imino_:,this%dst%jmino_:,this%dst%kmino_:), intent(out) :: A !< Needs to be (dst%imino_:dst%imaxo_,dst%jmino_:dst%jmaxo_,dst%kmino_:dst%kmaxo_)
      integer :: n
      ! Pull the data
      do n=1,this%nrecv
         A(this%dstind(1,n),this%dstind(2,n),this%dstind(3,n))=this%data_recv(n)
      end do
      ! Sync it before returning
      call this%dst%sync(A)
      ! Deallocate recv buffer
      if (allocated(this%data_recv)) deallocate(this%data_recv)
   end subroutine pull
   
   
   !> Routine that transfers the data from src to dst - both src_group and dst_group processors need to call
   !> Allocate recv buffer and deallocate send buffer here
   subroutine transfer(this)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(coupler), intent(inout) :: this
      integer :: ierr
      ! Allocate recv buffer (and maybe send buffer too)
      if (.not.allocated(this%data_recv)) allocate(this%data_recv(this%nrecv))
      if (.not.allocated(this%data_send)) allocate(this%data_send(this%nsend))
      ! Transfer data
      call MPI_ALLtoALLv(this%data_send,this%nsend_proc,this%nsend_disp,MPI_REAL_WP,this%data_recv,this%nrecv_proc,this%nrecv_disp,MPI_REAL_WP,this%comm,ierr)
      ! Deallocate send buffer
      if (allocated(this%data_send)) deallocate(this%data_send)
   end subroutine transfer
   
   
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
   

end module coupler_class
