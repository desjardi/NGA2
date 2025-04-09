!> Fast marching method class:
!> Provides support for creating a signed distance field to a distance Gmax from parallel (pgrid) array
module fmm_class
   use precision,   only: WP
   use string,      only: str_medium
   use config_class, only: config
   implicit none
   private

   ! Expose type/constructor/methods
   public :: fmm

   ! Communication frequency
   integer, parameter :: local_counter_max = 50

   type heap_type
      real(WP) :: G
      integer  :: i,j,k
   end type heap_type

   !> fmm object definition 
   type :: fmm 
      ! This is our config
      class(config), pointer :: cfg
      ! This is the name of the fmm
      character(len=str_medium) :: name='UNNAMED_FFM'
      ! i,j,k's for close nodes
      integer :: imin_close,imax_close
      integer :: jmin_close,jmax_close
      integer :: kmin_close,kmax_close
      ! i,j,k's on what needs to send to other procs
      integer :: i_passlo,i_passhi
      integer :: j_passlo,j_passhi
      integer :: k_passlo,k_passhi
      ! rank of proc to send to
      integer :: rank_x_lo,rank_x_hi
      integer :: rank_y_lo,rank_y_hi
      integer :: rank_z_lo,rank_z_hi
      ! MPI buffering
      integer(1), dimension(:), allocatable :: mpi_buffer
      integer :: mpi_buffer_size
   
   contains
      procedure :: initialize
      procedure :: build
   end type fmm

contains

   !> Initialize the fmm class 
   subroutine initialize(this,cfg,name)
      use mpi_f08
      implicit none
      class(fmm), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: isource,idest,ierr
      
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))

      ! Point to cfg object
      this%cfg=>cfg

      ! Determine the ranks of the procs that this proc should send to  
      call MPI_CART_SHIFT(this%cfg%comm,0,-1,isource,idest,ierr); this%rank_x_lo=idest
      call MPI_CART_SHIFT(this%cfg%comm,0,+1,isource,idest,ierr); this%rank_x_hi=idest
      call MPI_CART_SHIFT(this%cfg%comm,1,-1,isource,idest,ierr); this%rank_y_lo=idest
      call MPI_CART_SHIFT(this%cfg%comm,1,+1,isource,idest,ierr); this%rank_y_hi=idest
      call MPI_CART_SHIFT(this%cfg%comm,2,-1,isource,idest,ierr); this%rank_z_lo=idest
      call MPI_CART_SHIFT(this%cfg%comm,2,+1,isource,idest,ierr); this%rank_z_hi=idest
      
      ! Set bounds on what nodes to send to other procs
      if (this%cfg%nx.gt.1) then
         this%i_passlo=this%cfg%imin_+1
         this%i_passhi=this%cfg%imax_-1
      else
         this%i_passlo=this%cfg%imin_
         this%i_passhi=this%cfg%imax_
      end if
      if (this%cfg%ny.gt.1) then
         this%j_passlo=this%cfg%jmin_+1
         this%j_passhi=this%cfg%jmax_-1
      else
         this%j_passlo=this%cfg%jmin_
         this%j_passhi=this%cfg%jmax_
      end if
      if (this%cfg%nz.gt.1) then
         this%k_passlo=this%cfg%kmin_+1
         this%k_passhi=this%cfg%kmax_-1
      else
         this%k_passlo=this%cfg%kmin_
         this%k_passhi=this%cfg%kmax_
      end if

      ! Set bounds for the ghost nodes to consider given problem dimensions
      if (this%cfg%nx.gt.1) then
         this%imin_close=this%cfg%imin_-1
         this%imax_close=this%cfg%imax_+1
      else
         this%imin_close=this%cfg%imin_
         this%imax_close=this%cfg%imax_
      end if
      if (this%cfg%ny.gt.1) then
         this%jmin_close=this%cfg%jmin_-1
         this%jmax_close=this%cfg%jmax_+1
      else
         this%jmin_close=this%cfg%jmin_
         this%jmax_close=this%cfg%jmax_
      end if
      if (this%cfg%nz.gt.1) then
         this%kmin_close=this%cfg%kmin_-1
         this%kmax_close=this%cfg%kmax_+1
      else
         this%kmin_close=this%cfg%kmin_
         this%kmax_close=this%cfg%kmax_
      end if

      ! Attach an mpi buffer for parallelization
      this%mpi_buffer_size = (4*8+MPI_BSEND_OVERHEAD)*6*max(local_counter_max,20)
      allocate(this%mpi_buffer(this%mpi_buffer_size))
      call MPI_BUFFER_ATTACH(this%mpi_buffer,this%mpi_buffer_size,ierr)
      
   end subroutine initialize

   !> Update the distance field using estimate provided in G function out to Gmax
   subroutine build(this,G,Gmax)
      use messager, only: die
      implicit none
      class(fmm), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(inout) :: G
      real(WP), intent(in) :: Gmax
      integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: phi_flag
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: phi_fmm
      integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,3) :: stc_plus,stc_minus
      !integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: order_fmm
      integer :: n_plus,n_minus
      integer :: iter
      integer :: close_count,close_minus_count,close_plus_count
      integer :: fmm_accepted,fmm_close,fmm_far
      integer, dimension(this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,3), target :: close_minus_ijk,close_plus_ijk
      integer, dimension(:,:), pointer :: close_ijk
      ! Counter and mapping for accepted nodes
      integer, dimension(:,:), allocatable :: accepted_ijk
      ! Combined counter and mapping for all accepted nodes
      integer :: n_all_accepted
      integer, dimension(:,:), allocatable :: all_accepted_ijk
      ! Tags for cell status 
      integer, parameter :: fmm_far_plus       = 1
      integer, parameter :: fmm_far_minus      = 2
      integer, parameter :: fmm_close_plus     = 3
      integer, parameter :: fmm_close_minus    = 4
      integer, parameter :: fmm_accepted_plus  = 5
      integer, parameter :: fmm_accepted_minus = 6
      integer, parameter :: fmm_tmp            = 7
      ! Heap data
      type(heap_type), dimension(:),     allocatable :: heap  
      integer,         dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: heap_map
      integer :: nheap
      ! Communication
      integer, dimension(3) :: my_ibuf,ibuf
  
      ! Initialize the counters
      n_plus = 0
      n_minus = 0
      phi_flag = 0
      
      ! First tag all nodes and make plus and minus counts
      tag_nodes: block 
         integer :: i,j,k
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Cycle if too far or BC
                  !if (band(i,j,k).eq.0) cycle  !!!! without band then the n_plus & n_minus is large!!!
                  if (this%cfg%VF(i,j,k).eq.0.0_WP) cycle
                  ! Check with side
                  if (G(i,j,k).ge.0.0_WP) then
                     n_plus = n_plus + 1
                     phi_flag(i,j,k) = fmm_far_plus
                  else
                     n_minus = n_minus + 1
                     phi_flag(i,j,k) = fmm_far_minus
                  end if
               end do
            end do
         end do
      end block tag_nodes

      ! Initialize the close counters
      close_plus_count = 0
      close_minus_count = 0

      ! Set up list of initially close nodes
      setuplists: block 
         integer :: n,i,j,k
         integer :: ii,jj,kk
         do k=this%kmin_close,this%kmax_close
            do j=this%jmin_close,this%jmax_close
               do i=this%imin_close,this%imax_close
                  ! Cycle if too far or BC
                  !if (band(i,j,k).eq.0) cycle
                  if (this%cfg%VF(i,j,k).eq.0.0_WP) cycle
                  ! Check 6 direct neighbors
                  do n=1,6
                     select case(n)
                     case(1)
                        ii=i-1;jj=j;kk=k
                     case(2)
                        ii=i+1;jj=j;kk=k
                     case(3)
                        ii=i;jj=j-1;kk=k
                     case(4)
                        ii=i;jj=j+1;kk=k
                     case(5)
                        ii=i;jj=j;kk=k-1
                     case(6)
                        ii=i;jj=j;kk=k+1
                     end select
                     ! Don't add nodes that are outside the bands or BC
                     !if (band(ii,jj,kk).eq.0) cycle
                     if (this%cfg%VF(ii,jj,kk).eq.0.0_WP) cycle
                     ! Find interface crossing
                     if ((G(i,j,k)*G(ii,jj,kk)).le.0.0_WP) then
                        if (G(i,j,k).lt.0.0_WP) then
                           if (phi_flag(i,j,k).ne.fmm_close_minus) then
                              phi_flag(i,j,k) = fmm_close_minus
                              close_minus_count = close_minus_count + 1
                              close_minus_ijk(close_minus_count,1) = i
                              close_minus_ijk(close_minus_count,2) = j
                              close_minus_ijk(close_minus_count,3) = k
                           end if
                        else
                           if (phi_flag(i,j,k).ne.fmm_close_plus) then
                              phi_flag(i,j,k) = fmm_close_plus
                              close_plus_count = close_plus_count + 1
                              close_plus_ijk(close_plus_count,1) = i
                              close_plus_ijk(close_plus_count,2) = j
                              close_plus_ijk(close_plus_count,3) = k
                           end if
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end block setuplists

      !! Allocate a heap for sorting purposes
      nheap = 0
      if (allocated(heap)) deallocate(heap)
      allocate(heap(max(n_minus,n_plus)))
      
      ! Allocate a list of accepted nodes
      if (allocated(accepted_ijk)) deallocate(accepted_ijk)
      allocate(accepted_ijk(max(n_minus,n_plus),3))
      
      ! Allocate a list of all accepted nodes
      if (allocated(all_accepted_ijk)) deallocate(all_accepted_ijk)
      allocate(all_accepted_ijk(n_minus+n_plus,3))

      ! Zero the stencils
      stc_plus = 0
      stc_minus= 0

      ! Set up the temp level set field variable
      setuptemp: block
         integer :: i,j,k
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  phi_fmm(i,j,k) = Gmax
               end do
            end do
         end do
      end block setuptemp

      ! First negative, then positive side of the front
      do iter=2,1,-1
         
         ! Use the ibuf to store:
         ! 1: count of received messages
         ! 2: count of sent messages
         ! 3: our current nheap
         
         my_ibuf(1:2) = 0   ! clear the message counters
         
         ! Switch between sides
         if (iter.eq.1) then
            fmm_accepted = fmm_accepted_plus
            fmm_close    = fmm_close_plus
            fmm_far      = fmm_far_plus
            close_count  = close_plus_count
            close_ijk    => close_plus_ijk
         else
            fmm_accepted = fmm_accepted_minus
            fmm_close    = fmm_close_minus
            fmm_far      = fmm_far_minus
            close_count  = close_minus_count
            close_ijk    => close_minus_ijk
         end if
         
         ! Set up the initial close nodes list
         close_nodes: block
            integer :: n,i,j,k
            integer :: nn,ii,jj,kk 
            integer :: local_index, n_nbrs
            real(WP) :: local_phi
            real(WP), dimension(6) :: G_nbrs
            real(WP), dimension(3,6) :: dx_nbrs
            integer,  dimension(6) :: index_nbrs
            do nn = 1,close_count
               ! Get the index
               i = close_ijk(nn,1) 
               j = close_ijk(nn,2) 
               k = close_ijk(nn,3)
               
               ! Zero the number of neighbors
               n_nbrs = 0
               ! Loop over neighbors       
               do n = 1,6
                  select case(n)
                  case(1)
                     ii=i-1;jj=j;kk=k
                     local_index = -1
                  case(2) 
                     ii=i+1;jj=j;kk=k
                     local_index = +1
                  case(3)
                     ii=i;jj=j-1;kk=k
                     local_index = -2
                  case(4)
                     ii=i;jj=j+1;kk=k
                     local_index = +2
                  case(5)
                     ii=i;jj=j;kk=k-1
                     local_index = -3
                  case(6)
                     ii=i;jj=j;kk=k+1
                     local_index = +3
                  end select
                  ! Don't add nodes that are outside the BC
                  if (this%cfg%VF(ii,jj,kk).eq.0.0_WP) cycle
                  ! Form local metrics
                  if (G(i,j,k)*G(ii,jj,kk).le.0.0_WP) then
                     if (G(i,j,k).eq.G(ii,jj,kk)) cycle
                     n_nbrs = n_nbrs + 1
                     G_nbrs(   n_nbrs) = 0.0_WP
                     index_nbrs(n_nbrs) = local_index 
                     dx_nbrs(1,n_nbrs) = abs(-G(i,j,k)*(this%cfg%xm(ii)-this%cfg%xm(i))/(G(ii,jj,kk)-G(i,j,k)))
                     dx_nbrs(2,n_nbrs) = abs(-G(i,j,k)*(this%cfg%ym(jj)-this%cfg%ym(j))/(G(ii,jj,kk)-G(i,j,k)))
                     dx_nbrs(3,n_nbrs) = abs(-G(i,j,k)*(this%cfg%zm(kk)-this%cfg%zm(k))/(G(ii,jj,kk)-G(i,j,k)))
                  end if
               end do

               ! Reset the front itself
               if (n_nbrs.gt.0) then
                  local_phi = phi_calc(n_nbrs,G_nbrs,index_nbrs,dx_nbrs)
                  phi_fmm(i,j,k) = min(phi_fmm(i,j,k),local_phi)
               end if
               ! Add node to the heap
               if (abs(phi_fmm(i,j,k)).lt.Gmax) then
                  nheap           = nheap + 1
                  heap(nheap)%G   = phi_fmm(i,j,k)
                  heap(nheap)%i   = i
                  heap(nheap)%j   = j
                  heap(nheap)%k   = k
                  heap_map(i,j,k) = nheap
               else
                  phi_fmm(i,j,k) = sign(Gmax,phi_fmm(i,j,k))
                  phi_flag(i,j,k) = fmm_far
               end if
            end do
         end block close_nodes
         
         ! Heapify close nodes
         call heapify(nheap)

         ! Start FMM
         mainFMM: block
            integer :: n_accepted
            logical :: global_done,local_done
            integer :: local_counter
            n_accepted = 0
            global_done = .FALSE.
            global_loop: do while (.not.global_done)
               local_counter = 0
               local_done = .FALSE.     
               local_loop: do while (.not.local_done)
                  
                  !> Parallel message processing
                  message_processing: block
                     integer :: i,j,k 
                     integer :: ii,jj,kk
                     real(WP) :: this_phi
                     integer :: iheap
                     ! Check for new message from other processors
                     do while (multiphase_fmm_recv(i,j,k,this_phi))

                        if (this_phi.lt.phi_fmm(i,j,k)) then
                           
                           ! Check if this value is less than the currently
                           ! accepted value, and if it is, then roll back.
                           do while (n_accepted.gt.0)
                              ii = accepted_ijk(n_accepted,1)
                              jj = accepted_ijk(n_accepted,2)
                              kk = accepted_ijk(n_accepted,3)
                              
                              ! Check if we have rolled back enough...
                              if (this_phi.gt.phi_fmm(ii,jj,kk)) then
                                 if (ii.ge.this%cfg%imin_.and.ii.le.this%cfg%imax_) then
                                    if (jj.ge.this%cfg%jmin_.and.jj.le.this%cfg%jmax_) then
                                       if (kk.ge.this%cfg%kmin_.and.kk.le.this%cfg%kmax_) then
                                          exit
                                       end if
                                    end if
                                 end if
                              end if
                              
                              phi_flag(ii,jj,kk) = fmm_close
                              nheap = nheap + 1
                              heap(nheap)%G = phi_fmm(ii,jj,kk)
                              heap(nheap)%i = ii
                              heap(nheap)%j = jj
                              heap(nheap)%k = kk
                              heap_map(ii,jj,kk) = nheap
                              call upsift(nheap)
                              ! ... and reduce the accepted list.
                              n_accepted = n_accepted - 1
                              
                           end do
                           
                           ! Now take the sent value and take action - we must be either
                           ! far or close. If we were accepted, we should have been rolled back
                           if (phi_flag(i,j,k).eq.fmm_close) then
                              ! Reset node value and resort
                              phi_fmm(i,j,k) = this_phi
                              iheap = heap_map(i,j,k)
                              heap(iheap)%G = this_phi
                              call upsift(iheap)
                           elseif (phi_flag(i,j,k).eq.fmm_far) then
                              ! Add as close node and sort
                              if (abs(this_phi).lt.Gmax) then 
                                 phi_fmm(i,j,k) = this_phi
                                 phi_flag(i,j,k) = fmm_close
                                 nheap = nheap + 1
                                 heap(nheap)%G = this_phi
                                 heap(nheap)%i = i
                                 heap(nheap)%j = j
                                 heap(nheap)%k = k
                                 heap_map(i,j,k) = nheap
                                 call upsift(nheap)
                              else
                                 phi_fmm(i,j,k) = sign(Gmax,this_phi)
                                 phi_flag(i,j,k) = fmm_far
                              end if
                           end if
                           
                        end if
                        ! Add one to recv message counter...
                        my_ibuf(1) = my_ibuf(1) + 1
                     end do
                  end block message_processing
      
                  ! If heap not empty
                  nonemptyheap: if (nheap.gt.0) then

                     ! Accept the closest and add neighbors
                     accept_closest: block
                        integer :: i,j,k 
                        integer :: n_already_close,n_new_close
                        integer, dimension(6) :: already_close
                        integer, dimension(6) :: new_close
                        ! Get index
                        i = heap(1)%i 
                        j = heap(1)%j 
                        k = heap(1)%k
                        ! Set the closest as accepted              
                        phi_flag(i,j,k) = fmm_accepted
                        local_counter = local_counter + 1
                        n_accepted = n_accepted + 1
                        accepted_ijk(n_accepted,1) = i
                        accepted_ijk(n_accepted,2) = j
                        accepted_ijk(n_accepted,3) = k
                        ! Now delete the heap root
                        call delroot(nheap)

                        ! Add neighbors to the heap
                        add_neighbors: block 
                           integer :: ii,jj,kk
                           integer :: nn,n_nbrs
                           integer :: local_index
                           ! Loop over neighbors
                           n_already_close = 0
                           n_new_close = 0
                           n_nbrs = 0
                           do nn = 1,6
                              select case(nn)
                              case(1)
                                 ii=max(this%imin_close,i-1);jj=j;kk=k
                                 if (i.eq.ii) cycle
                                 local_index = -1
                              case(2)
                                 ii=min(this%imax_close,i+1);jj=j;kk=k
                                 if (i.eq.ii) cycle
                                 local_index = +1
                              case(3)
                                 ii=i;jj=max(this%jmin_close,j-1);kk=k
                                 if (j.eq.jj) cycle
                                 local_index = -2
                              case(4)
                                 ii=i;jj=min(this%jmax_close,j+1);kk=k
                                 if (j.eq.jj) cycle
                                 local_index = +2
                              case(5)
                                 ii=i;jj=j;kk=max(this%kmin_close,k-1)
                                 if (k.eq.kk) cycle
                                 local_index = -3
                              case(6)
                                 ii=i;jj=j;kk=min(this%kmax_close,k+1)
                                 if (k.eq.kk) cycle
                                 local_index = +3
                              end select
                              ! Don't add nodes that are outside the BC
                              if (this%cfg%VF(ii,jj,kk).eq.0.0_WP) cycle
                              ! Count the nodes to be used in extending the distance function
                              if (phi_flag(ii,jj,kk).eq.fmm_close) then
                                 n_already_close = n_already_close + 1
                                 already_close(n_already_close) = nn
                                 phi_flag(ii,jj,kk) = fmm_tmp
                              else if (phi_flag(ii,jj,kk).eq.fmm_far) then
                                 n_new_close = n_new_close + 1
                                 new_close(n_new_close) = nn
                                 phi_flag(ii,jj,kk) = fmm_tmp
                              end if
                           end do
                           end block add_neighbors
                     
                        ! Work on the physical domain only and let other processors 
                        ! know a boundary node has been accepted
                        if ( i.ge.this%cfg%imin_ .and. i.le.this%cfg%imax_ .and. &
                             j.ge.this%cfg%jmin_ .and. j.le.this%cfg%jmax_ .and. &
                             k.ge.this%cfg%kmin_ .and. k.le.this%cfg%kmax_ ) then
                           call multiphase_fmm_send(i,j,k,phi_fmm(i,j,k),my_ibuf(2))
                        end if
                  
                        ! Process already close nodes
                        process_close_nodes: block
                           integer :: ii,jj,kk
                           integer :: nnn,iii,jjj,kkk
                           integer :: m,nn,n_nbrs
                           integer :: local_index
                           real(WP) :: local_phi
                           real(WP), dimension(6) :: phi_nbrs
                           real(WP), dimension(3,6) :: dx_nbrs
                           integer,  dimension(6) :: index_nbrs
                           do nn = 1,n_already_close
                              m = already_close(nn)  
                              select case(m)
                              case(1)
                                 ii=max(this%imin_close,i-1);jj=j;kk=k
                                 if (i.eq.ii) cycle
                                 local_index = -1
                              case(2)
                                 ii=min(this%imax_close,i+1);jj=j;kk=k
                                 if (i.eq.ii) cycle
                                 local_index = +1
                              case(3)
                                 ii=i;jj=max(this%jmin_close,j-1);kk=k
                                 if (j.eq.jj) cycle
                                 local_index = -2
                              case(4)
                                 ii=i;jj=min(this%jmax_close,j+1);kk=k
                                 if (j.eq.jj) cycle
                                 local_index = +2
                              case(5)
                                 ii=i;jj=j;kk=max(this%kmin_close,k-1)
                                 if (k.eq.kk) cycle
                                 local_index = -3
                              case(6)
                                 ii=i;jj=j;kk=min(this%kmax_close,k+1)
                                 if (k.eq.kk) cycle
                                 local_index = +3
                              end select
                              ! Change the flag
                              phi_flag(ii,jj,kk) = fmm_close
                              ! Loop over neighbors
                              n_nbrs = 0
                              do nnn = 1,6
                                 select case(nnn)
                                 case(1)
                                    iii=ii-1;jjj=jj;kkk=kk
                                    local_index = -1
                                 case(2)
                                    iii=ii+1;jjj=jj;kkk=kk
                                    local_index = +1
                                 case(3)
                                    iii=ii;jjj=jj-1;kkk=kk
                                    local_index = -2
                                 case(4)
                                    iii=ii;jjj=jj+1;kkk=kk
                                    local_index = +2
                                 case(5)
                                    iii=ii;jjj=jj;kkk=kk-1
                                    local_index = -3
                                 case(6)
                                    iii=ii;jjj=jj;kkk=kk+1
                                    local_index = +3
                                 end select
                                 ! Don't add nodes that are outside the BC
                                 if (this%cfg%VF(iii,jjj,kkk).eq.0.0_WP) cycle
                                 ! Check for nbrs and look for...
                                 if ((G(ii,jj,kk)*G(iii,jjj,kkk)).le.0.0_WP) then
                                    if (G(ii,jj,kk).eq.G(iii,jjj,kkk)) cycle
                                    ! ... an intersection
                                    n_nbrs = n_nbrs + 1
                                    phi_nbrs(n_nbrs) = 0.0_WP
                                    index_nbrs(n_nbrs) = local_index
                                    dx_nbrs(1,n_nbrs) = abs(-G(ii,jj,kk)*(this%cfg%xm(iii)-this%cfg%xm(ii))/(G(iii,jjj,kkk)-G(ii,jj,kk)))
                                    dx_nbrs(2,n_nbrs) = abs(-G(ii,jj,kk)*(this%cfg%ym(jjj)-this%cfg%ym(jj))/(G(iii,jjj,kkk)-G(ii,jj,kk)))
                                    dx_nbrs(3,n_nbrs) = abs(-G(ii,jj,kk)*(this%cfg%zm(kkk)-this%cfg%zm(kk))/(G(iii,jjj,kkk)-G(ii,jj,kk)))
                                 else if (phi_flag(iii,jjj,kkk).eq.fmm_accepted) then
                                    ! ... an accepted nbr
                                    n_nbrs = n_nbrs + 1
                                    phi_nbrs(n_nbrs) = phi_fmm(iii,jjj,kkk)
                                    index_nbrs(n_nbrs) = local_index
                                    dx_nbrs(1,n_nbrs) = abs(this%cfg%xm(iii)-this%cfg%xm(ii))
                                    dx_nbrs(2,n_nbrs) = abs(this%cfg%ym(jjj)-this%cfg%ym(jj))
                                    dx_nbrs(3,n_nbrs) = abs(this%cfg%zm(kkk)-this%cfg%zm(kk))
                                 end if
                              end do
                              ! Recompute nodal values
                              if (n_nbrs.gt.0) then
                                 local_phi = phi_calc(n_nbrs,phi_nbrs,index_nbrs,dx_nbrs)
                                 phi_fmm(ii,jj,kk) = min(phi_fmm(ii,jj,kk),local_phi)
                              end if
                              ! Modify the heap
                              m = heap_map(ii,jj,kk)
                              if (phi_fmm(ii,jj,kk).gt.heap(m)%G) then
                                 heap(m)%G = phi_fmm(ii,jj,kk)
                                 call downsift(nheap,m)
                              else
                                 heap(m)%G = phi_fmm(ii,jj,kk)
                                 call upsift(m)
                              end if
                           end do
                        end block process_close_nodes
                  
                        ! Process new close nodes
                        process_new_nodes: block
                           integer :: ii,jj,kk
                           integer :: iii,jjj,kkk
                           integer :: m,nn,nnn,n_nbrs
                           integer :: local_index
                           real(WP) :: local_phi
                           real(WP), dimension(3,6) :: dx_nbrs
                           real(WP), dimension(6) :: phi_nbrs
                           integer,  dimension(6) :: index_nbrs
                           do nn = 1,n_new_close
                              m = new_close(nn)  
                              select case(m)
                              case(1)
                                 ii=max(this%imin_close,i-1);jj=j;kk=k
                                 if (i.eq.ii) cycle
                                 local_index = -1
                              case(2)
                                 ii=min(this%imax_close,i+1);jj=j;kk=k
                                 if (i.eq.ii) cycle
                                 local_index = +1
                              case(3)
                                 ii=i;jj=max(this%jmin_close,j-1);kk=k
                                 if (j.eq.jj) cycle
                                 local_index = -2
                              case(4)
                                 ii=i;jj=min(this%jmax_close,j+1);kk=k
                                 if (j.eq.jj) cycle
                                 local_index = +2
                              case(5)
                                 ii=i;jj=j;kk=max(this%kmin_close,k-1)
                                 if (k.eq.kk) cycle
                                 local_index = -3
                              case(6)
                                 ii=i;jj=j;kk=min(this%kmax_close,k+1)
                                 if (k.eq.kk) cycle
                                 local_index = +3
                              end select
                              ! Change the flag
                              phi_flag(ii,jj,kk) = fmm_close
                              ! Loop over neighbors
                              n_nbrs = 0
                              do nnn = 1,6
                                 select case(nnn)
                                 case(1)
                                    iii=ii-1;jjj=jj;kkk=kk
                                    local_index = -1
                                 case(2)
                                    iii=ii+1;jjj=jj;kkk=kk
                                    local_index = +1
                                 case(3)
                                    iii=ii;jjj=jj-1;kkk=kk
                                    local_index = -2
                                 case(4)
                                    iii=ii;jjj=jj+1;kkk=kk
                                    local_index = +2
                                 case(5)
                                    iii=ii;jjj=jj;kkk=kk-1
                                    local_index = -3
                                 case(6)
                                    iii=ii;jjj=jj;kkk=kk+1
                                    local_index = +3
                                 end select
                                 ! Don't add nodes that are outside the BC
                                 if (this%cfg%VF(iii,jjj,kkk).eq.0.0_WP) cycle
                                 ! Check for nbrs and look for...
                                 if ((G(ii,jj,kk)*G(iii,jjj,kkk)).le.0.0_WP) then
                                    if (G(ii,jj,kk).eq.G(iii,jjj,kkk)) cycle
                                    ! ... an intersection
                                    n_nbrs = n_nbrs + 1
                                    index_nbrs(n_nbrs) = local_index
                                    phi_nbrs(n_nbrs) = 0.0_WP
                                    dx_nbrs(1,n_nbrs) = abs(-G(ii,jj,kk)*(this%cfg%xm(iii)-this%cfg%xm(ii))/(G(iii,jjj,kkk)-G(ii,jj,kk)))
                                    dx_nbrs(2,n_nbrs) = abs(-G(ii,jj,kk)*(this%cfg%ym(jjj)-this%cfg%ym(jj))/(G(iii,jjj,kkk)-G(ii,jj,kk)))
                                    dx_nbrs(3,n_nbrs) = abs(-G(ii,jj,kk)*(this%cfg%zm(kkk)-this%cfg%zm(kk))/(G(iii,jjj,kkk)-G(ii,jj,kk)))
                                 else if (phi_flag(iii,jjj,kkk).eq.fmm_accepted) then
                                    ! ... an accepted nbr
                                    n_nbrs = n_nbrs + 1
                                    index_nbrs(n_nbrs) = local_index
                                    phi_nbrs(n_nbrs) = phi_fmm(iii,jjj,kkk)
                                    dx_nbrs(1,n_nbrs) = abs(this%cfg%xm(iii)-this%cfg%xm(ii))
                                    dx_nbrs(2,n_nbrs) = abs(this%cfg%ym(jjj)-this%cfg%ym(jj))
                                    dx_nbrs(3,n_nbrs) = abs(this%cfg%zm(kkk)-this%cfg%zm(kk))
                                 end if
                              end do
                              ! Recompute nodal value
                              if (n_nbrs.gt.0) then
                                 local_phi = phi_calc(n_nbrs,phi_nbrs,index_nbrs,dx_nbrs)
                                 phi_fmm(ii,jj,kk) = min(phi_fmm(ii,jj,kk),local_phi)
                              end if
                              ! Add to the heap
                              if (abs(phi_fmm(ii,jj,kk)).lt.Gmax) then
                                 nheap = nheap + 1
                                 heap(nheap)%G = phi_fmm(ii,jj,kk)
                                 heap(nheap)%i = ii
                                 heap(nheap)%j = jj
                                 heap(nheap)%k = kk
                                 heap_map(ii,jj,kk) = nheap
                                 call upsift(nheap)
                              else
                                 phi_fmm(ii,jj,kk) = sign(Gmax,phi_fmm(ii,jj,kk))
                                 phi_flag(ii,jj,kk) = fmm_far
                              end if
                           end do
                        end block process_new_nodes
                     end block accept_closest
                  end if nonemptyheap
                  
                  local_done = ((nheap.eq.0).or.(local_counter.eq.local_counter_max))
                  
               end do local_loop

               communcate_messages:block
                  use mpi_f08, only: MPI_ALLREDUCE, MPI_SUM, MPI_INTEGER
                  integer :: ierr
                  my_ibuf(3) = nheap
                  call MPI_ALLREDUCE(my_ibuf,ibuf,3,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
                  global_done = ((ibuf(1).eq.ibuf(2)).and.(ibuf(3).eq.0))
               end block communcate_messages
               
            end do global_loop
            
            ! Recover the accepted nodes list 1) Positive 2) Negative
            recover_accepted: block
               select case(iter)
               case(1)
                  all_accepted_ijk(n_all_accepted+1:n_all_accepted+n_accepted,:)=accepted_ijk(1:n_accepted,:)
                  n_all_accepted=n_all_accepted+n_accepted
               case(2)
                  all_accepted_ijk(1:n_accepted,:)=accepted_ijk(n_accepted:1:-1,:)
                  n_all_accepted=n_accepted
               end select
            end block recover_accepted

         end block mainFMM
      end do ! sides

      ! Now update level set
      update_levelset: block
         integer:: i,j,k
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (phi_flag(i,j,k).eq.fmm_accepted_plus) then
                     G(i,j,k) = +phi_fmm(i,j,k)
                  else if (phi_flag(i,j,k).eq.fmm_accepted_minus) then
                     G(i,j,k) = -phi_fmm(i,j,k)
                  end if
               end do
            end do
         end do
         ! Communicate level set
         call this%cfg%sync(G)
      end block update_levelset

         ! ! Store ordered list
         ! do n=1,n_all_accepted
         !    ! Get index
         !    i=all_accepted_ijk(n,1)
         !    j=all_accepted_ijk(n,2)
         !    k=all_accepted_ijk(n,3)
         !    order_fmm(i,j,k)=n
         ! end do

         ! ! Generate the stencils for variable extension
         ! do n=1,n_all_accepted
         !    ! Get index
         !    i=all_accepted_ijk(n,1)
         !    j=all_accepted_ijk(n,2)
         !    k=all_accepted_ijk(n,3)
         !    ! Find useable neighbors in plus direction
         !    n_nbrs = 0
         !    do nn = 1,6
         !       select case(nn)
         !       case(1)
         !          ii=max(imin_close,i-1);jj=j;kk=k
         !          if (i.eq.ii) cycle
         !          local_index = -1
         !       case(2)
         !          ii=min(imax_close,i+1);jj=j;kk=k
         !          if (i.eq.ii) cycle
         !          local_index = +1
         !       case(3)
         !          ii=i;jj=max(jmin_close,j-1);kk=k
         !          if (j.eq.jj) cycle
         !          local_index = -2
         !       case(4)
         !          ii=i;jj=min(jmax_close,j+1);kk=k
         !          if (j.eq.jj) cycle
         !          local_index = +2
         !       case(5)
         !          ii=i;jj=j;kk=max(kmin_close,k-1)
         !          if (k.eq.kk) cycle
         !          local_index = -3
         !       case(6)
         !          ii=i;jj=j;kk=min(kmax_close,k+1)
         !          if (k.eq.kk) cycle
         !          local_index = +3
         !       end select
         !       ! Don't consider BC nodes or outside bands
         !    !   if (band(ii,jj,kk).eq.0) cycle
         !       if (vol(ii,jj,kk).eq.0.0_WP) cycle
         !       ! Check if neighbor is closer
         !       if (order_fmm(ii,jj,kk).lt.order_fmm(i,j,k).and.G(ii,jj,kk).lt.G(i,j,k)) then
         !          n_nbrs = n_nbrs + 1
         !          phi_nbrs(n_nbrs) = G(ii,jj,kk)
         !          index_nbrs(n_nbrs) = local_index
         !          dx_nbrs(1,n_nbrs) = abs(xm(ii)-xm(i))
         !          dx_nbrs(2,n_nbrs) = abs(ym(jj)-ym(j))
         !          dx_nbrs(3,n_nbrs) = abs(zm(kk)-zm(k))
         !       end if
         !    end do
         !    ! Build the stencil
         !    if (n_nbrs.gt.0) then
         !       ! Check smallest distance
         !       if (icyl.eq.1) then 
         !          radius = abs(ym(j))
         !       else
         !          radius = 1.0_WP
         !       end if
         !       call phi_calc(local_phi,radius)
         !       ! Keep only the stencil
         !       stc_plus(i,j,k,:) = stc(:)
         !    end if
         !    ! Find useable neighbors in minus direction
         !    n_nbrs = 0
         !    do nn = 1,6
         !       select case(nn)
         !       case(1)
         !          ii=max(imin_close,i-1);jj=j;kk=k
         !          if (i.eq.ii) cycle
         !          local_index = -1
         !       case(2)
         !          ii=min(imax_close,i+1);jj=j;kk=k
         !          if (i.eq.ii) cycle
         !          local_index = +1
         !       case(3)
         !          ii=i;jj=max(jmin_close,j-1);kk=k
         !          if (j.eq.jj) cycle
         !          local_index = -2
         !       case(4)
         !          ii=i;jj=min(jmax_close,j+1);kk=k
         !          if (j.eq.jj) cycle
         !          local_index = +2
         !       case(5)
         !          ii=i;jj=j;kk=max(kmin_close,k-1)
         !          if (k.eq.kk) cycle
         !          local_index = -3
         !       case(6)
         !          ii=i;jj=j;kk=min(kmax_close,k+1)
         !          if (k.eq.kk) cycle
         !          local_index = +3
         !       end select
         !       ! Don't consider BC nodes or outside bands
         !    !   if (band(ii,jj,kk).eq.0) cycle
         !       if (vol(ii,jj,kk).eq.0.0_WP) cycle
         !       ! Check if neighbor is farther
         !       if (order_fmm(ii,jj,kk).gt.order_fmm(i,j,k).and.G(ii,jj,kk).gt.G(i,j,k)) then
         !          n_nbrs = n_nbrs + 1
         !          phi_nbrs(n_nbrs) = G(ii,jj,kk)
         !          index_nbrs(n_nbrs) = local_index
         !          dx_nbrs(1,n_nbrs) = abs(xm(ii)-xm(i))
         !          dx_nbrs(2,n_nbrs) = abs(ym(jj)-ym(j))
         !          dx_nbrs(3,n_nbrs) = abs(zm(kk)-zm(k))
         !       end if
            ! end do
         !    ! Build the stencil
         !    if (n_nbrs.gt.0) then
         !       ! Check smallest distance
         !       if (icyl.eq.1) then 
         !          radius = abs(ym(j))
         !       else
         !          radius = 1.0_WP
         !       end if
         !       call phi_calc(local_phi,radius)
         !       ! Keep only the stencil
         !       stc_minus(i,j,k,:) = stc(:)
         !    end if

   ! ! ====================================== !
   ! ! Constant scalar extension based on FMM !
   ! ! ====================================== !
   ! subroutine extend(A,Aflag,dir)
   !    use multiphase_fmm
   !    use parallel
   !    implicit none
      
   !    ! Input data
   !    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
   !    integer,  dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: Aflag
   !    character(len=*), intent(in) :: dir
      
   !    ! Work variables
   !    integer  :: n,nn,ierr
   !    integer  :: i,j,k,i0,j0,k0
   !    integer  :: ii,jj,kk
   !    real(WP) :: radius,local_sc,this_sc
   !    integer, dimension(2) :: my_ibuf
   !    integer, dimension(2) :: ibuf
      
   !    ! Initialize
   !    my_ibuf=0
   !    ibuf=0
      
   !    ! Barrier
   !    call MPI_BARRIER(comm,ierr)
      
   !    ! Select work direction
   !    select case (trim(dir))
   !    case ('+')
         
   !       do n=1,n_all_accepted
            
   !          ! Get index
   !          i=all_accepted_ijk(n,1)
   !          j=all_accepted_ijk(n,2)
   !          k=all_accepted_ijk(n,3)
            
   !          ! If already set, skip cell
   !          if (Aflag(i,j,k).eq.1) cycle
            
   !          ! If outside, skip cell
   !          if (i.lt.imin_.or.i.gt.imax_.or.j.lt.jmin_.or.j.gt.jmax_.or.k.lt.kmin_.or.k.gt.kmax_) cycle
            
   !          ! If wall, skip cell
   !          if (vol(i,j,k).eq.0.0_WP) cycle
            
   !          ! Receive messages if present
   !          do while (multiphase_fmm_recv(i0,j0,k0,this_sc))
   !             my_ibuf(1)=my_ibuf(1)+1
   !             A(i0,j0,k0)=this_sc
   !             Aflag(i0,j0,k0)=1
   !          end do
            
   !          ! Use the stencil to find neighbors
   !          n_nbrs = 0
   !          neighbors_plus: do nn=1,3
   !             select case(stc_plus(i,j,k,nn))
   !             case(-1)
   !                ii=i-1;jj=j;kk=k
   !             case(+1)
   !                ii=i+1;jj=j;kk=k
   !             case(-2)
   !                ii=i;jj=j-1;kk=k
   !             case(+2)
   !                ii=i;jj=j+1;kk=k
   !             case(-3)
   !                ii=i;jj=j;kk=k-1
   !             case(+3)
   !                ii=i;jj=j;kk=k+1
   !             case (0)
   !                cycle neighbors_plus
   !             end select
   !             ! Check if the data is available
   !             do while (Aflag(ii,jj,kk).eq.0)
   !                call multiphase_fmm_brecv(i0,j0,k0,this_sc)
   !                my_ibuf(1)=my_ibuf(1)+1
   !                A(i0,j0,k0)=this_sc
   !                Aflag(i0,j0,k0)=1
   !             end do
   !             n_nbrs = n_nbrs + 1
   !             sc_nbrs(n_nbrs) = A(ii,jj,kk)
   !             phi_nbrs(n_nbrs) = G(ii,jj,kk)
   !             index_nbrs(n_nbrs) = stc_plus(i,j,k,nn)
   !             dx_nbrs(1,n_nbrs) = abs(xm(ii)-xm(i))
   !             dx_nbrs(2,n_nbrs) = abs(ym(jj)-ym(j))
   !             dx_nbrs(3,n_nbrs) = abs(zm(kk)-zm(k))
   !          end do neighbors_plus
            
   !          ! Extend the quantity
   !          phi_me = G(i,j,k)
   !          if (icyl.eq.1) then 
   !             radius = abs(ym(j))
   !          else
   !             radius = 1.0_WP
   !          end if
   !          call sc_calc(local_sc,radius)
   !          A(i,j,k) = local_sc
   !          Aflag(i,j,k) = 1
            
   !          ! Communicate the value if necessary
   !          call multiphase_fmm_send(i,j,k,A(i,j,k),my_ibuf(2))
            
   !       end do
         
   !       ! Cleaning up remaining messages
   !       call MPI_BARRIER(comm,ierr)
   !       call MPI_ALLREDUCE(my_ibuf,ibuf,2,MPI_INTEGER,MPI_SUM,comm,ierr)
   !       do while (ibuf(1).ne.ibuf(2))
   !          do while (multiphase_fmm_recv(i0,j0,k0,this_sc))
   !             my_ibuf(1)=my_ibuf(1)+1
   !             A(i0,j0,k0)=this_sc
   !             Aflag(i0,j0,k0)=1
   !          end do
   !          call MPI_BARRIER(comm,ierr)
   !          call MPI_ALLREDUCE(my_ibuf,ibuf,2,MPI_INTEGER,MPI_SUM,comm,ierr)
   !       end do
         
   !    case ('-')
         
   !       do n=n_all_accepted,1,-1
            
   !          ! Get index
   !          i=all_accepted_ijk(n,1)
   !          j=all_accepted_ijk(n,2)
   !          k=all_accepted_ijk(n,3)
            
   !          ! If already set, skip cell
   !          if (Aflag(i,j,k).eq.1) cycle
            
   !          ! If outside, skip cell
   !          if (i.lt.imin_.or.i.gt.imax_.or.j.lt.jmin_.or.j.gt.jmax_.or.k.lt.kmin_.or.k.gt.kmax_) cycle
            
   !          ! If wall, skip cell
   !          if (vol(i,j,k).eq.0.0_WP) cycle
            
   !          ! Receive messages if present
   !          do while (multiphase_fmm_recv(i0,j0,k0,this_sc))
   !             my_ibuf(1)=my_ibuf(1)+1
   !             A(i0,j0,k0)=this_sc
   !             Aflag(i0,j0,k0)=1
   !          end do
            
   !          ! Use the stencil to find neighbors
   !          n_nbrs = 0
   !          neighbors_minus: do nn=1,3
   !             select case(stc_minus(i,j,k,nn))
   !             case(-1)
   !                ii=i-1;jj=j;kk=k
   !             case(+1)
   !                ii=i+1;jj=j;kk=k
   !             case(-2)
   !                ii=i;jj=j-1;kk=k
   !             case(+2)
   !                ii=i;jj=j+1;kk=k
   !             case(-3)
   !                ii=i;jj=j;kk=k-1
   !             case(+3)
   !                ii=i;jj=j;kk=k+1
   !             case (0)
   !                cycle neighbors_minus
   !             end select
   !             ! Check if the data is available
   !             do while (Aflag(ii,jj,kk).eq.0)
   !                call multiphase_fmm_brecv(i0,j0,k0,this_sc)
   !                my_ibuf(1)=my_ibuf(1)+1
   !                A(i0,j0,k0)=this_sc
   !                Aflag(i0,j0,k0)=1
   !             end do
   !             n_nbrs = n_nbrs + 1
   !             sc_nbrs(n_nbrs) = A(ii,jj,kk)
   !             phi_nbrs(n_nbrs) = G(ii,jj,kk)
   !             index_nbrs(n_nbrs) = stc_minus(i,j,k,nn)
   !             dx_nbrs(1,n_nbrs) = abs(xm(ii)-xm(i))
   !             dx_nbrs(2,n_nbrs) = abs(ym(jj)-ym(j))
   !             dx_nbrs(3,n_nbrs) = abs(zm(kk)-zm(k))
   !          end do neighbors_minus
            
   !          ! Extend the quantity
   !          phi_me = G(i,j,k)
   !          if (icyl.eq.1) then 
   !             radius = abs(ym(j))
   !          else
   !             radius = 1.0_WP
   !          end if
   !          call sc_calc(local_sc,radius)
   !          A(i,j,k) = local_sc
   !          Aflag(i,j,k) = 1
            
   !          ! Communicate the value if necessary
   !          call multiphase_fmm_send(i,j,k,A(i,j,k),my_ibuf(2))
            
   !       end do
         
   !       ! Cleaning up remaining messages
   !       call MPI_BARRIER(comm,ierr)
   !       call MPI_ALLREDUCE(my_ibuf,ibuf,2,MPI_INTEGER,MPI_SUM,comm,ierr)
   !       do while (ibuf(1).ne.ibuf(2))
   !          do while (multiphase_fmm_recv(i0,j0,k0,this_sc))
   !             my_ibuf(1)=my_ibuf(1)+1
   !             A(i0,j0,k0)=this_sc
   !             Aflag(i0,j0,k0)=1
   !          end do
   !          call MPI_BARRIER(comm,ierr)
   !          call MPI_ALLREDUCE(my_ibuf,ibuf,2,MPI_INTEGER,MPI_SUM,comm,ierr)
   !       end do
         
   !    end select
      
   !    ! Barrier
   !    call MPI_BARRIER(comm,ierr)
      
   !    return
   ! end subroutine extend

   contains

      ! ===================================== !
      ! Make the subtree of heap starting in  !
      ! parent node fulfil the heap property  !
      ! -nh  : heap size                      !
      ! -node: index of subtree of interest   !
      !                                       !
      ! Usage: if parent value has increased, !
      !        it needs to be sifted down     !
      ! ===================================== !
      subroutine downsift(nh,node)
         implicit none
         
         integer, intent(in) :: nh,node
         integer :: parent,child
         type(heap_type) :: myheap
         
         ! Copy parent to bring down heap
         parent=node
         myheap=heap(parent)
         
         ! First child of parent
         child=2*parent
            
         ! Continue until bottom of heap
         down_sift: do while (child.le.nh)
            
            ! Is there a second child?
            if (child+1.le.nh) then
               ! Pick closest of two children
               if (heap(child+1)%G.lt.heap(child)%G) child=child+1
            end if
            
            ! Check heap property
            if (myheap%G.le.heap(child)%G) exit down_sift
            
            ! Otherwise need to swap child and parent
            heap(parent)=heap(child); heap_map(heap(parent)%i,heap(parent)%j,heap(parent)%k)=parent
            
            ! Child is our new parent
            parent=child
            
            ! First child of new parent
            child=2*parent
            
         end do down_sift
         
         ! Finish up by putting the start node in place of last child (stored in parent)
         heap(parent)=myheap; heap_map(heap(parent)%i,heap(parent)%j,heap(parent)%k)=parent
         
         return
      end subroutine downsift
      
      ! ======================================= !
      ! Make the child node in heap is at the   !
      ! right location by moving it up the tree !
      ! -node: index of node of interest        !
      !                                         !
      ! Usage: if child value has decreased, it !
      !        needs to be sifted up            !
      ! ======================================= !
      subroutine upsift(node)
         implicit none
         
         integer, intent(in) :: node
         integer :: parent,child
         type(heap_type) :: myheap
         
         ! Copy child to bring up heap
         child=node
         myheap=heap(child)
         
         ! Continue until top of heap
         up_sift: do while (child.gt.1) 
            
            ! Get parent
            parent=child/2
            
            ! If parent is smaller than start value minval is already low, so exit
            if (heap(parent)%G.le.myheap%G) exit up_sift
            
            ! Otherwise need to swap child with parent
            heap(child)=heap(parent); heap_map(heap(child)%i,heap(child)%j,heap(child)%k)=child
            
            ! Parent is our new child
            child=parent
            
         end do up_sift
         
         ! Put the stored start value into whatever parent node swapped with it
         heap(child)=myheap
         heap_map(heap(child)%i,heap(child)%j,heap(child)%k)=child
         
         return
      end subroutine upsift
      
      ! ===================== !
      ! Heap creation routine !
      ! ===================== !
      subroutine heapify(nh)
         implicit none
         
         integer, intent(in) :: nh
         integer :: i
         
         ! Traverse the array and build heap
         do i=nh/2,1,-1
            call downsift(nh,i)
         end do
         
         return
      end subroutine heapify
      
      ! ===================== !
      ! Root deletion routine !
      ! ===================== !
      subroutine delroot(nh)
         implicit none
         
         integer, intent(inout) :: nh
         
         ! Check that heap is not empty
         if (nh.gt.1) then
            
            ! Move bottom of heap to top of heap
            heap(1)=heap(nh)
            
            ! Sift top value down the heap
            call downsift(nh-1,1)
            
         end if
         
         ! Reduce heap size
         nh=nh-1
         
         return
      end subroutine delroot
      
      !> Compute distance of a close cell to neighboring accepted cells 
      function phi_calc(n_nbrs,G_nbrs,index_nbrs,dx_nbrs) result(G_loc)
         implicit none
         integer, intent(in) :: n_nbrs 
         real(WP), dimension(  n_nbrs), intent(in) :: G_nbrs 
         integer,  dimension(  n_nbrs), intent(in) :: index_nbrs 
         real(WP), dimension(3,n_nbrs), intent(in) :: dx_nbrs
         integer :: n
         real(WP) :: G_loc, G_tmp
         integer :: local_1,local_2,local_3
         integer :: loctmp1,loctmp2,loctmp3,loctmp4,loctmp5
         integer, dimension(3)    :: stc

         ! Compute distance using the number of accepted neighbors
         select case (n_nbrs)
         case(0)
            ! Do nothing: G_loc already set
            
         case(1)
            G_loc = G_nbrs(1) + sqrt(sum(dx_nbrs(1:3,1)**2))
            stc(1) = index_nbrs(1)
            
         case(2) ! 2 possible cases - 1) the two nbrs either lie across node from each other or 2) they don't
            ! Nbrs on both sides of the node
            if (abs(index_nbrs(1)).eq.abs(index_nbrs(2))) then
               G_loc = G_nbrs(1) + sqrt(sum(dx_nbrs(1:3,1)**2))
               stc(1) = index_nbrs(1)
               G_tmp = G_nbrs(2) + sqrt(sum(dx_nbrs(1:3,2)**2))
               if (G_tmp.lt.G_loc) then
                  G_loc = G_tmp
                  stc(1) = index_nbrs(2)
               end if
            else
               ! Nbrs only on one side of the node
               call phi_calc_2D(G_nbrs(1),G_nbrs(2),dx_nbrs(1:3,1),dx_nbrs(1:3,2),G_loc)
               stc(1) = index_nbrs(1); stc(2) = index_nbrs(2)
            end if
            
         case(3) ! 2 possible cases - 1) any two of the nbrs lie across node from each other or 2) no two nbrs do
            ! First check for the crossing cases
            if (abs(index_nbrs(1)).eq.abs(index_nbrs(2))) then
               ! Two cases
               ! #1 --------------------------------------------
               local_1 = 1
               local_2 = 3
               call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_loc)
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
               ! #2 --------------------------------------------
               local_1 = 2
               local_2 = 3
               call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_tmp)
               if (G_tmp.lt.G_loc) then
                  G_loc = G_tmp
                  stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
               end if
            else if (abs(index_nbrs(2)).eq.abs(index_nbrs(3))) then
               ! Two cases
               ! #1 --------------------------------------------
               local_1 = 1
               local_2 = 2
               call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_loc)
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
               ! #2 --------------------------------------------
               local_1 = 1
               local_2 = 3
               call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_tmp)
               if (G_tmp.lt.G_loc) then
                  G_loc = G_tmp
                  stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
               end if
            else
               ! And now the non-crossed case: three nbrs all along different directions
               call phi_calc_3D(G_nbrs(1),G_nbrs(2),G_nbrs(3),dx_nbrs(1:3,1),dx_nbrs(1:3,2),dx_nbrs(1:3,3),G_loc)
               stc(1) = index_nbrs(1); stc(2) = index_nbrs(2); stc(3) = index_nbrs(3)
            end if
            
         case(4) ! 2 possible cases: either 1) both sets of nbrs lie across node from each other, or 2) only one set of nbrs does
            
            if (abs(index_nbrs(1)).eq.abs(index_nbrs(2))) then
               if (abs(index_nbrs(3)).eq.abs(index_nbrs(4))) then
                  ! Try 4 different 2D cases
                  ! #1 --------------------------------------------
                  local_1 = 1
                  local_2 = 3
                  call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_loc)
                  stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
                  ! #2 --------------------------------------------
                  local_1 = 1
                  local_2 = 4
                  call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_tmp)
                  if (G_tmp.lt.G_loc) then
                     G_loc = G_tmp
                     stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
                  end if
                  ! #3 --------------------------------------------
                  local_1 = 2
                  local_2 = 3
                  call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_tmp)
                  if (G_tmp.lt.G_loc) then
                     G_loc = G_tmp
                     stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
                  end if
                  ! #4 --------------------------------------------
                  local_1 = 2
                  local_2 = 4
                  call phi_calc_2D(G_nbrs(local_1),G_nbrs(local_2),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),G_tmp)
                  if (G_tmp.lt.G_loc) then
                     G_loc = G_tmp
                     stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2)
                  end if
               else
                  ! Try 2 different 3D cases
                  ! #1 --------------------------------------------
                  local_1 = 1
                  local_2 = 3
                  local_3 = 4
                  call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_loc)
                  stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
                  ! #2 --------------------------------------------
                  local_1 = 2
                  local_2 = 3
                  local_3 = 4
                  call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
                  if (G_tmp.lt.G_loc) then
                     G_loc = G_tmp
                     stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
                  end if
               end if
            else if (abs(index_nbrs(2)).eq.abs(index_nbrs(3))) then
               ! Try 2 different 3D cases
               ! #1 --------------------------------------------
               local_1 = 1
               local_2 = 2
               local_3 = 4
               call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_loc)
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
               ! #2 --------------------------------------------
               local_1 = 1
               local_2 = 3
               local_3 = 4
               call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
               if (G_tmp.lt.G_loc) then
                  G_loc = G_tmp
                  stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
               end if
            else 
               ! Try 2 different 3D cases
               ! #1 --------------------------------------------
               local_1 = 1
               local_2 = 2
               local_3 = 3
               call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_loc)
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
               ! #2 --------------------------------------------
               local_1 = 1
               local_2 = 2
               local_3 = 4
               call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
               if (G_tmp.lt.G_loc) then
                  G_loc = G_tmp
                  stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
               end if
            end if
            
         case(5) ! Only 1 case - test using the 4 possible planes and take smallest solution
            
            ! Renumber so that the 5th node has no one across from him
            if (abs(index_nbrs(1)).ne.abs(index_nbrs(2))) then
               loctmp1 = 2
               loctmp2 = 3
               loctmp3 = 4
               loctmp4 = 5
               loctmp5 = 1
            else if (abs(index_nbrs(3)).ne.abs(index_nbrs(4))) then
               loctmp1 = 4
               loctmp2 = 5
               loctmp3 = 1
               loctmp4 = 2
               loctmp5 = 3
            else 
               loctmp1 = 1
               loctmp2 = 2
               loctmp3 = 3
               loctmp4 = 4
               loctmp5 = 5
            end if
            ! #1 --------------------------------------------
            local_1 = loctmp1
            local_2 = loctmp3
            local_3 = loctmp5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_loc)
            stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            ! #2 --------------------------------------------
            local_1 = loctmp1
            local_2 = loctmp4
            local_3 = loctmp5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #3 --------------------------------------------
            local_1 = loctmp2
            local_2 = loctmp3
            local_3 = loctmp5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #4 --------------------------------------------
            local_1 = loctmp2
            local_2 = loctmp4
            local_3 = loctmp5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            
         case(6) ! Only 1 case - test using the 8 possible planes and take smallest solution
            
            ! #1 --------------------------------------------
            local_1 = 1
            local_2 = 3
            local_3 = 5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_loc)
            stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            ! #2 --------------------------------------------
            local_1 = 1
            local_2 = 3
            local_3 = 6
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #3 --------------------------------------------
            local_1 = 1
            local_2 = 4
            local_3 = 5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #4 --------------------------------------------
            local_1 = 1
            local_2 = 4
            local_3 = 6
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #5 --------------------------------------------
            local_1 = 2
            local_2 = 3
            local_3 = 5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #6 --------------------------------------------
            local_1 = 2
            local_2 = 4
            local_3 = 5
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #7 --------------------------------------------
            local_1 = 2
            local_2 = 3
            local_3 = 6
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            ! #8 --------------------------------------------
            local_1 = 2
            local_2 = 4
            local_3 = 6
            call phi_calc_3D(G_nbrs(local_1),G_nbrs(local_2),G_nbrs(local_3),dx_nbrs(1:3,local_1),dx_nbrs(1:3,local_2),dx_nbrs(1:3,local_3),G_tmp)
            if (G_tmp.lt.G_loc) then
               G_loc = G_tmp
               stc(1) = index_nbrs(local_1); stc(2) = index_nbrs(local_2); stc(3) = index_nbrs(local_3)
            end if
            
         end select

         ! Machine accuracy fix to guaranty ordering
         do n=1,n_nbrs
            if (G_loc.lt.G_nbrs(n)) G_loc=G_nbrs(n)
         end do
      end function phi_calc
    
      subroutine phi_calc_2D(G1,G2,dx1,dx2,G_out)
         implicit none
         real(WP), intent(in)               :: G1,G2
         real(WP), dimension(3), intent(in) :: dx1,dx2
         real(WP), intent(out)              :: G_out
         real(WP)                           :: a,b,c,d1,d2,root1,root2
         real(WP)                           :: d1i,d2i
         
         ! Square of the distance between nodes
         d1 = sum(dx1(1:3)**2)
         d2 = sum(dx2(1:3)**2)
         
         ! Take the inverse
         d1i = 1.0_WP/d1
         d2i = 1.0_WP/d2
         
         ! Value of G given by root of quadratic equation
         a = d1i + d2i
         b = -2.0_WP*G1*d1i -2.0_WP*G2*d2i
         c = (G1**2)*d1i + (G2**2)*d2i - 1.0_WP
         
         root1 = (-b + sqrt(max(b**2 - 4.0_WP*a*c,0.0_WP)))/(2.0_WP*a)
         root2 = (-b - sqrt(max(b**2 - 4.0_WP*a*c,0.0_WP)))/(2.0_WP*a)
         
         ! Have to be careful about selecting these roots...
         ! Don't want negative values
         if (root1.le.0.0_WP) root1 = huge(1.0_WP)
         if (root2.le.0.0_WP) root2 = huge(1.0_WP)
         ! Also don't want values smaller than our original two G values
         if (root1.le.G1)     root1 = huge(1.0_WP)
         if (root1.le.G2)     root1 = huge(1.0_WP)
         if (root2.le.G1)     root2 = huge(1.0_WP)
         if (root2.le.G2)     root2 = huge(1.0_WP)
         
         G_out =  min(root1,root2)
         
         return
      end subroutine phi_calc_2D
      
      subroutine phi_calc_3D(G1,G2,G3,dx1,dx2,dx3,G_out)
         implicit none
         real(WP), intent(in)               :: G1,G2,G3
         real(WP), dimension(3), intent(in) :: dx1,dx2,dx3
         real(WP), intent(out)              :: G_out
         real(WP)                           :: a,b,c,d1,d2,d3,root1,root2
         
         ! Square of the distance between nodes
         d1 = sum(dx1(1:3)**2)
         d2 = sum(dx2(1:3)**2)
         d3 = sum(dx3(1:3)**2)

         ! Take the inverse
         d1 = 1.0_WP/d1
         d2 = 1.0_WP/d2
         d3 = 1.0_WP/d3  
         
         ! Value of G given by root of quadratic eqn
         a = d1 + d2 + d3
         b = -2.0_WP*G1*d1 -2.0_WP*G2*d2 -2.0_WP*G3*d3
         c = (G1**2)*d1 + (G2**2)*d2 + (G3**2)*d3 - 1.0_WP
         
         root1 = (-b + sqrt(max(b**2 - 4.0_WP*a*c,0.0_WP)))/(2.0_WP*a)
         root2 = (-b - sqrt(max(b**2 - 4.0_WP*a*c,0.0_WP)))/(2.0_WP*a)
         
         ! Have to be careful about selecting these roots...
         ! Don't want negative values
         if (root1.le.0.0_WP) root1 = huge(1.0_WP)
         if (root2.le.0.0_WP) root2 = huge(1.0_WP)
         ! Also don't want values smaller than our original two G values
         if (root1.le.G1)     root1 = huge(1.0_WP)
         if (root1.le.G2)     root1 = huge(1.0_WP)
         if (root1.le.G3)     root1 = huge(1.0_WP)
         if (root2.le.G1)     root2 = huge(1.0_WP)
         if (root2.le.G2)     root2 = huge(1.0_WP)
         if (root2.le.G3)     root2 = huge(1.0_WP)
         
         G_out =  min(root1,root2)
         
         return
      end subroutine phi_calc_3D

      ! =================================== !
      ! Parallel receive of single quantity !
      ! =================================== !
      function multiphase_fmm_recv(i,j,k,phi_value) result(imessage)
         use mpi_f08, only : MPI_Iprobe, MPI_Recv, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_STATUS, MPI_SOURCE
         use parallel, only: MPI_REAL_WP
         implicit none
         
         integer,   intent(out) :: i,j,k
         real(WP),  intent(out) :: phi_value
         real(WP), dimension(4) :: val_recv
         type(MPI_Status)   :: status
         integer :: ierr,isource
         logical :: imessage
         
         ! Probe for message
         call MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,this%cfg%comm,imessage,status,ierr)
         
         ! If message is present, receive it
         if (imessage) then
            isource=status%MPI_SOURCE
            call MPI_Recv(val_recv,4,MPI_REAL_WP,isource,0,this%cfg%comm,status,ierr)
            i=nint(val_recv(1))
            j=nint(val_recv(2))
            k=nint(val_recv(3))
            phi_value=val_recv(4)
         end if
         
         return
      end function multiphase_fmm_recv

      ! ====================================== !
      ! Blocking receive for extension routine !
      ! ====================================== !
      subroutine multiphase_fmm_brecv(i,j,k,phi_value)
         use mpi_f08, only: MPI_Recv, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_Status
         use parallel, only: MPI_REAL_WP
         implicit none
         
         integer,  intent(out) :: i,j,k
         real(WP), intent(out) :: phi_value
         real(WP), dimension(4) :: recv_buf
         type(MPI_Status)   :: status
         integer :: ierr
         
         ! Blocking receive for extension
         call MPI_Recv(recv_buf,4,MPI_REAL_WP,MPI_ANY_SOURCE,MPI_ANY_TAG,this%cfg%comm,status,ierr)
         i=nint(recv_buf(1))
         j=nint(recv_buf(2))
         k=nint(recv_buf(3))
         phi_value=recv_buf(4)
         
         return
      end subroutine multiphase_fmm_brecv
      
      
      ! ================================ !
      ! Parallel send of single quantity !
      ! ================================ !
      subroutine multiphase_fmm_send(i,j,k,value,counter)
         implicit none
         
         integer,  intent(in) :: i,j,k
         real(WP), intent(in) :: value
         integer, intent(inout) :: counter
         integer :: i0,j0,k0
         
         ! Communicate if necessary
         if (this%rank_x_lo.ge.0 .and. i.lt.this%i_passlo) then
            i0=i;j0=j;k0=k
            if (this%cfg%xper .and.this%cfg%iproc.eq.1  ) i0=i+this%cfg%nx
            call multiphase_fmm_parallel_send(this%rank_x_lo,i0,j0,k0,value)
            counter=counter+1
         end if
         if (this%rank_x_hi.ge.0 .and. i.gt.this%i_passhi) then
            i0=i;j0=j;k0=k
            if (this%cfg%xper .and. this%cfg%iproc.eq.this%cfg%npx) i0=i-this%cfg%nx
            call multiphase_fmm_parallel_send(this%rank_x_hi,i0,j0,k0,value)
            counter=counter+1
         end if
         if (this%rank_y_lo.ge.0 .and. j.lt.this%j_passlo) then
            i0=i;j0=j;k0=k
            if (this%cfg%yper .and. this%cfg%jproc.eq.1  ) j0=j+this%cfg%ny
            call multiphase_fmm_parallel_send(this%rank_y_lo,i0,j0,k0,value)
            counter=counter+1
         end if
         if (this%rank_y_hi.ge.0 .and. j.gt.this%j_passhi) then
            i0=i;j0=j;k0=k
            if (this%cfg%yper .and. this%cfg%jproc.eq.this%cfg%npy) j0=j-this%cfg%ny
            call multiphase_fmm_parallel_send(this%rank_y_hi,i0,j0,k0,value)
            counter=counter+1
         end if
         if (this%rank_z_lo.ge.0 .and. k.lt.this%k_passlo) then
            i0=i;j0=j;k0=k
            if (this%cfg%zper .and. this%cfg%kproc.eq.1  ) k0=k+this%cfg%nz
            call multiphase_fmm_parallel_send(this%rank_z_lo,i0,j0,k0,value)
            counter=counter+1
         end if
         if (this%rank_z_hi.ge.0 .and. k.gt.this%k_passhi) then
            i0=i;j0=j;k0=k
            if (this%cfg%zper .and. this%cfg%kproc.eq.this%cfg%npz) k0=k-this%cfg%nz
            call multiphase_fmm_parallel_send(this%rank_z_hi,i0,j0,k0,value)
            counter=counter+1
         end if
      end subroutine multiphase_fmm_send

      
      subroutine multiphase_fmm_parallel_send(idest,i,j,k,phi_value)
         use mpi_f08, only: MPI_Bsend
         use parallel, only: MPI_REAL_WP
         implicit none
         
         integer,  intent(in) :: idest
         integer,  intent(in) :: i,j,k
         real(WP), intent(in) :: phi_value
         real(WP), dimension(4) :: buffer
         integer :: ierr
         
         buffer(1) = real(i,WP)
         buffer(2) = real(j,WP)
         buffer(3) = real(k,WP)
         buffer(4) = phi_value
         
         ! Use a buffered send
         call MPI_Bsend(buffer,4,MPI_REAL_WP,idest,0,this%cfg%comm,ierr)
         
         return
      end subroutine multiphase_fmm_parallel_send
         
   end subroutine build
   
end module fmm_class