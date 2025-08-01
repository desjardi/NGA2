!> Interbox concept is defined here: it takes in two pgrid objects
!> and builds the facilty to exchange data between them.
!> The two meshes must be uniform and aligned, with integer refinement.
module interbox_class
   use precision,      only: WP
   use string,         only: str_medium
   use pgrid_class,    only: pgrid
   use mpi_f08
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: interbox,box
   
   ! Tolerance parameter
   real(WP), parameter :: tol=1.0e-12_WP
   
   ! Box type - simplified uniform Carteian mesh with domain decomposition
   type :: box
      integer , dimension(3) :: nx   !< Number of cells
      real(WP), dimension(3) :: dx   !< Mesh size
      real(WP), dimension(3) :: x0   !< Start location
      logical , dimension(3) :: per  !< Periodicity
      integer :: no                  !< Overlap size
      integer , dimension(3) :: np   !< Number of processors per direction
      integer , dimension(3) :: pi   !< Processor index per direction (pi in [1,np] if owner, 0 otherwise)
   end type box

   !> Interbox object definition
   type :: interbox
      
      ! These are our two pgrids
      type(pgrid), pointer :: crse=>NULL()                !< Coarse grid
      type(pgrid), pointer :: fine=>NULL()                !< Fine grid
      
      ! These are the two corresponding boxes
      type(box), pointer :: cb=>NULL()                    !< Coarse box
      type(box), pointer :: fb=>NULL()                    !< Fine box
      
      ! Logicals to help us know if we have received a coarse or a fine grid
      logical :: got_crse=.false.                         !< Were we given coarse grid
      logical :: got_fine=.false.                         !< Were we given fine grid
      
      ! This is our communication information
      type(MPI_Comm)  :: comm                             !< Intracommunicator over the union of both groups
      type(MPI_Group) :: cgrp,fgrp,grp                    !< Coarse and fine groups and their union
      integer :: nproc                                    !< Number of processors
      integer :: rank                                     !< Processor grid rank
      logical :: amRoot                                   !< Am I root for the interbox?
      integer :: croot                                    !< Rank of coarse grid root on union group
      integer :: froot                                    !< Rank of fine grid root on union group
      
      ! Rank maps
      integer, dimension(:,:,:), allocatable :: frankmap  !< Processor coordinate on fine mesh to union group rank map
      integer, dimension(:,:,:), allocatable :: crankmap  !< Processor coordinate on crse mesh to union group rank map
      
      ! ! Interpolation support
      ! integer , dimension(:,:), allocatable :: srcind     !< Src indices of dst points that this processor can interpolate
      ! integer , dimension(:,:), allocatable :: dstind     !< Dst indices of dst points that this processor will receive
      ! real(WP), dimension(:,:), allocatable :: w          !< Interpolation weights for dst points that this processor can interpolate
      
      ! ! Communication support
      ! integer :: nsend                                    !< Total number of dst points that this processor can interpolate and will send out
      ! integer, dimension(:), allocatable :: nsend_proc    !< Number of points to send to each processor
      ! integer, dimension(:), allocatable :: nsend_disp    !< Data displacement when sending to each processor
      ! integer :: nrecv                                    !< Total number of dst points that this processor will receive
      ! integer, dimension(:), allocatable :: nrecv_proc    !< Number of points to receive from each processor
      ! integer, dimension(:), allocatable :: nrecv_disp    !< Data displacement when receiving from each processor
      
      ! ! Data storage
      ! real(WP), dimension(:), allocatable :: data_send    !< Data to send
      ! real(WP), dimension(:), allocatable :: data_recv    !< Received data
      
   contains
      procedure :: create                                 !< Routine that creates an interbox object
      procedure :: set_fine_grid                          !< Routine that sets the fine grid
      procedure :: set_crse_grid                          !< Routine that sets the coarse grid
      procedure :: initialize                             !< Routine that initializes all interpolation metrics between 1 and 2
      procedure :: push                                   !< Src routine that pushes a field into our send data storage
      procedure :: pull                                   !< Dst routine that pulls a field from our received data storage
      procedure :: transfer                               !< Routine that performs the src->dst data transfer
      procedure :: finalize                               !< Finalize interbox object
   end type interbox
   
   
contains
   
   
   !> Interbox creation
   subroutine create(this,crse_grp,fine_grp,name)
      use messager, only: die
      use parallel, only: comm
      implicit none
      class(interbox), intent(inout) :: this
      type(MPI_Group), intent(in) :: crse_grp,fine_grp
      character(len=*), intent(in) :: name
      integer, dimension(1) :: rankin,rankout
      integer :: ierr
      ! Build group union
      this%cgrp=crse_grp; this%fgrp=fine_grp; call MPI_GROUP_UNION(this%cgrp,this%fgrp,this%grp,ierr)
      ! Gather some info for communication
      call MPI_GROUP_SIZE(this%grp,this%nproc,ierr); if (this%nproc.eq.0)            call die('[interbox initialize] Somehow the union of both groups is of size zero')
      call MPI_GROUP_RANK(this%grp,this%rank ,ierr); if (this%rank.eq.MPI_UNDEFINED) call die('[interbox initialize] All processors that call the constructor need to be in one of the two groups')
      ! Create intracommunicator for the new group
      call MPI_COMM_CREATE_GROUP(comm,this%grp,0,this%comm,ierr)
      ! Find roots for both grids on the shared communicator
      rankin=0; call MPI_GROUP_TRANSLATE_RANKS(this%cgrp,1,rankin,this%grp,rankout,ierr); this%croot=rankout(1)
      rankin=0; call MPI_GROUP_TRANSLATE_RANKS(this%fgrp,1,rankin,this%grp,rankout,ierr); this%froot=rankout(1)
      ! Set coupler root to fine root
      this%amRoot=(this%rank.eq.this%froot)
   end subroutine create
   
   
   !> Set the fine grid - to be called by processors in fine group
   subroutine set_fine_grid(this,pg)
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      ! Set pointer to fine grid
      this%fine=>pg; this%got_fine=.true.
      ! Check uniformity
      if (.not.this%fine%uniform_x) call die('[interbox] fine mesh is non-uniform in x')
      if (.not.this%fine%uniform_y) call die('[interbox] fine mesh is non-uniform in y')
      if (.not.this%fine%uniform_z) call die('[interbox] fine mesh is non-uniform in z')
      ! Create fb
      allocate(this%fb)
      this%fb%nx =[this%fine%nx                ,this%fine%ny                ,this%fine%nz                ]
      this%fb%dx =[this%fine%dx(this%fine%imin),this%fine%dy(this%fine%jmin),this%fine%dz(this%fine%kmin)]
      this%fb%x0 =[this%fine%x (this%fine%imin),this%fine%y (this%fine%jmin),this%fine%z (this%fine%kmin)]
      this%fb%per=[this%fine%xper              ,this%fine%yper              ,this%fine%zper              ]
      this%fb%no = this%fine%no
      this%fb%np =[this%fine%npx               ,this%fine%npy               ,this%fine%npz               ]
      this%fb%pi =[this%fine%iproc             ,this%fine%jproc             ,this%fine%kproc             ]
   end subroutine set_fine_grid
   
   
   !> Set the coarse grid - to be called by processors in coarse group
   subroutine set_crse_grid(this,pg)
      use messager, only: die
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      ! Set pointer to coarse grid
      this%crse=>pg; this%got_crse=.true.
      ! Check uniformity
      if (.not.this%crse%uniform_x) call die('[interbox] coarse mesh is non-uniform in x')
      if (.not.this%crse%uniform_y) call die('[interbox] coarse mesh is non-uniform in y')
      if (.not.this%crse%uniform_z) call die('[interbox] coarse mesh is non-uniform in z')
      ! Create cb
      allocate(this%cb)
      this%cb%nx =[this%crse%nx                ,this%crse%ny                ,this%crse%nz                ]
      this%cb%dx =[this%crse%dx(this%crse%imin),this%crse%dy(this%crse%jmin),this%crse%dz(this%crse%kmin)]
      this%cb%x0 =[this%crse%x (this%crse%imin),this%crse%y (this%crse%jmin),this%crse%z (this%crse%kmin)]
      this%cb%per=[this%crse%xper              ,this%crse%yper              ,this%crse%zper              ]
      this%cb%no = this%crse%no
      this%cb%np =[this%crse%npx               ,this%crse%npy               ,this%crse%npz               ]
      this%cb%pi =[this%crse%iproc             ,this%crse%jproc             ,this%crse%kproc             ]
   end subroutine set_crse_grid
   
   
   !> Initialize interbox 
   subroutine initialize(this)
      use messager, only: die
      implicit none
      class(interbox), intent(inout) :: this
      
      ! Make fb and cb grid available to all ranks
      share_grids: block
         use parallel, only: MPI_REAL_WP
         ! Coarse grid root broadcasts cb to the group
         if (.not.this%got_crse) allocate(this%cb)   ! Allocate it first
         call MPI_BCAST(this%cb%nx ,3,MPI_INTEGER,this%croot,this%comm,ierr)
         call MPI_BCAST(this%cb%dx ,3,MPI_REAL_WP,this%croot,this%comm,ierr)
         call MPI_BCAST(this%cb%x0 ,3,MPI_REAL_WP,this%croot,this%comm,ierr)
         call MPI_BCAST(this%cb%per,3,MPI_LOGICAL,this%croot,this%comm,ierr)
         call MPI_BCAST(this%cb%no ,1,MPI_INTEGER,this%croot,this%comm,ierr)
         call MPI_BCAST(this%cb%np ,3,MPI_INTEGER,this%croot,this%comm,ierr)
         if (.not.this%got_crse) this%cb%pi=[0,0,0]  ! Assign default proc index for non-owners
         ! Fine grid root broadcasts fb to the group
         if (.not.this%got_fine) allocate(this%fb)   ! Allocate it first
         call MPI_BCAST(this%fb%nx ,3,MPI_INTEGER,this%froot,this%comm,ierr)
         call MPI_BCAST(this%fb%dx ,3,MPI_REAL_WP,this%froot,this%comm,ierr)
         call MPI_BCAST(this%fb%x0 ,3,MPI_REAL_WP,this%froot,this%comm,ierr)
         call MPI_BCAST(this%fb%per,3,MPI_LOGICAL,this%froot,this%comm,ierr)
         call MPI_BCAST(this%fb%no ,1,MPI_INTEGER,this%froot,this%comm,ierr)
         call MPI_BCAST(this%fb%np ,3,MPI_INTEGER,this%froot,this%comm,ierr)
         if (.not.this%got_fine) this%fb%pi=[0,0,0]  ! Assign default proc index for non-owners
      end block share_grids
      
      ! Check refinement and alignement
      check_refinement_alignement: block
         real(WP) :: ratio,delta
         ! Ensure fine is finer than crse in all directions
         if (this%fb%dx(1).gt.this%cb%dx(1)-tol) call die('[interbox] fine grid is not finer than coarse grid in direction x')
         if (this%fb%dx(2).gt.this%cb%dx(2)-tol) call die('[interbox] fine grid is not finer than coarse grid in direction y')
         if (this%fb%dx(3).gt.this%cb%dx(3)-tol) call die('[interbox] fine grid is not finer than coarse grid in direction z')
         ! Ensure integer refinement ratio
         ratio=this%cb%dx(1)/this%fb%dx(1); this%ratio(1)=nint(ratio); if (abs(ratio-this%ratio(1)).gt.tol) call die('[interbox] refinement ratio in direction x is not an integer')
         ratio=this%cb%dx(2)/this%fb%dx(2); this%ratio(2)=nint(ratio); if (abs(ratio-this%ratio(2)).gt.tol) call die('[interbox] refinement ratio in direction y is not an integer')
         ratio=this%cb%dx(3)/this%fb%dx(3); this%ratio(3)=nint(ratio); if (abs(ratio-this%ratio(3)).gt.tol) call die('[interbox] refinement ratio in direction z is not an integer')
         ! Ensure grids are aligned
         delta=this%fb%x0(1)-this%cb%x0(1); if (abs(mod(delta,this%fb%dx(1))).gt.tol) call die('[interbox] fine and coarse grids are not aligned in direction x')
         delta=this%fb%x0(2)-this%cb%x0(2); if (abs(mod(delta,this%fb%dx(2))).gt.tol) call die('[interbox] fine and coarse grids are not aligned in direction y')
         delta=this%fb%x0(3)-this%cb%x0(3); if (abs(mod(delta,this%fb%dx(3))).gt.tol) call die('[interbox] fine and coarse grids are not aligned in direction z')
      end block check_refinement_alignement
      
      ! Make partition maps available to all
      share_partition: block
         integer :: ierr,n
         integer, dimension(:), allocatable :: fiproc,fjproc,fkproc
         integer, dimension(:), allocatable :: ciproc,cjproc,ckproc
         ! Prepare communication arrays
         allocate(fiproc(0:this%nproc-1),fjproc(0:this%nproc-1),fkproc(0:this%nproc-1))
         allocate(ciproc(0:this%nproc-1),cjproc(0:this%nproc-1),ckproc(0:this%nproc-1))
         ! Allgather the rank->(iproc,jproc,kproc) info
         call MPI_ALLGATHER(this%fine%pi(1),1,MPI_INTEGER,fiproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%fine%pi(2),1,MPI_INTEGER,fjproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%fine%pi(3),1,MPI_INTEGER,fkproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%crse%pi(1),1,MPI_INTEGER,ciproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%crse%pi(2),1,MPI_INTEGER,cjproc,1,MPI_INTEGER,this%comm,ierr)
         call MPI_ALLGATHER(this%crse%pi(3),1,MPI_INTEGER,ckproc,1,MPI_INTEGER,this%comm,ierr)
         ! Allocate rankmaps
         allocate(this%frankmap(this%fb%np(1),this%fb%np(2),this%fb%np(3)))
         allocate(this%crankmap(this%cb%np(1),this%cb%np(2),this%cb%np(3)))
         ! Finally, flip the rankmap data
         do n=0,this%nproc-1
            if (fiproc(n).gt.0) this%frankmap(fiproc(n),fjproc(n),fkproc(n))=n
            if (ciproc(n).gt.0) this%crankmap(ciproc(n),cjproc(n),ckproc(n))=n
         end do
         ! Deallocate communication arrays
         deallocate(fiproc,fjproc,fkproc)
         deallocate(ciproc,cjproc,ckproc)
      end block share_partition
      
      ! Fine ranks prepare range of coarse indices they are responsible for
      fine_to_coarse_prep: block
         if (this%got_fine) then
            ! Convert
            this%imin_crse=
            this%imax_crse=
            this%jmin_crse=
            this%jmax_crse=
            this%kmin_crse=
            this%kmax_crse=
            this%ncrse=
         end if
            i_dst_min = (this%src%imin_ - this%src%imin) * refx + this%dst%imin
            i_dst_max = (this%src%imax_ - this%src%imin) * refx + this%dst%imin + refx - 1
            j_dst_min = (this%src%jmin_ - this%src%jmin) * refy + this%dst%jmin
            j_dst_max = (this%src%jmax_ - this%src%jmin) * refy + this%dst%jmin + refy - 1
            k_dst_min = (this%src%kmin_ - this%src%kmin) * refz + this%dst%kmin
            k_dst_max = (this%src%kmax_ - this%src%kmin) * refz + this%dst%kmin + refz - 1
         end if

      end block fine_to_coarse_prep
      
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
      
   end subroutine initialize
   
   
   
   
   
   !> Prepare interpolation metrics from src to dst
   subroutine initialize(this)
      implicit none
      class(coupler), intent(inout) :: this
      integer , dimension(:)  , allocatable :: rk
      integer , dimension(:,:), allocatable :: dstind
      
      
      
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
   
   
   !> Finalize coupler object
   subroutine finalize(this)
      implicit none
      class(coupler), intent(inout) :: this
      integer :: ierr
      ! Deallocate allocatable arrays
      if (allocated(this%rankmap))     deallocate(this%rankmap)
      if (allocated(this%srcind))      deallocate(this%srcind)
      if (allocated(this%dstind))      deallocate(this%dstind)
      if (allocated(this%w))           deallocate(this%w)
      if (allocated(this%nsend_proc))  deallocate(this%nsend_proc)
      if (allocated(this%nsend_disp))  deallocate(this%nsend_disp)
      if (allocated(this%nrecv_proc))  deallocate(this%nrecv_proc)
      if (allocated(this%nrecv_disp))  deallocate(this%nrecv_disp)
      if (allocated(this%data_send))   deallocate(this%data_send)
      if (allocated(this%data_recv))   deallocate(this%data_recv)
      ! Nullify grid pointers
      nullify(this%src)
      nullify(this%dst)
      ! Free only the internally created communicator and group
      if (this%comm.ne.MPI_COMM_NULL ) call MPI_Comm_free(this%comm,ierr)
      if (this%grp .ne.MPI_GROUP_NULL) call MPI_Group_free(this%grp,ierr)
      ! Invalidate external groups
      this%sgrp=MPI_GROUP_NULL
      this%dgrp=MPI_GROUP_NULL
   end subroutine finalize
   
   
end module coupler_class
