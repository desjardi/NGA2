!> Iterator concept is defined here.
!> The intention is to have the ability to work in a sub-region of a pgrid.
!> This is directly useful for boundary conditions, but probably for more.
module itr_class
   use precision,   only: WP
   use string,      only: str_medium
   use pgrid_class, only: pgrid
   use mpi_f08,     only: MPI_Comm
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: itr
   
   !> itr object definition
   type :: itr
      ! Parallel info for the iterator
      type(MPI_Comm) :: comm                              !< Communicator for the iterator
      integer :: nproc                                    !< Number of processes in that iterator
      integer :: rank                                     !< Rank for that iterator
      logical :: amIn                                     !< Am I working in this iterator?
      logical :: amRoot                                   !< Am I the root for this iterator
      ! This is our partitioned grid
      class(pgrid), pointer :: pg                         !< This is the pgrid the iterator is build for
      ! This is the name of the iterator
      character(len=str_medium) :: name='UNNAMED_ITR'     !< Iterator name (default=UNNAMED_ITR)
      ! This is the unstructured mapping
      integer :: no_,n_                                   !< Total number of local cells in iterator, with and without overlap
      integer, dimension(:,:), allocatable :: map         !< This is our map for the iterator
      
   contains
      procedure :: print=>itr_print                       !< Output iterator to the screen
   end type itr
   
   
   !> Declare single-grid iterator constructor
   interface itr
      procedure construct_from_function
   end interface itr
   
contains

   !> Iterator constructor from a tester function
   function construct_from_function(pg,name,test_if_inside) result(self)
      use mpi_f08, only: MPI_COMM_SPLIT,MPI_COMM_NULL,MPI_UNDEFINED
      implicit none
      type(itr) :: self
      class(pgrid), target, intent(in) :: pg
      character(len=*), intent(in) :: name
      interface
         logical function test_if_inside(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function test_if_inside
      end interface
      integer :: i,j,k,cnt
      integer :: color,key,ierr
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%pg=>pg
      
      ! Loop over local domain with overlap to count iterator cells
      self%n_=0; self%no_=0
      do k=self%pg%kmino_,self%pg%kmaxo_
         do j=self%pg%jmino_,self%pg%jmaxo_
            do i=self%pg%imino_,self%pg%imaxo_
               if (test_if_inside(self%pg,i,j,k)) then
                  self%no_=self%no_+1
                  if (i.ge.self%pg%imin_.and.i.le.self%pg%imax_.and.&
                  &   j.ge.self%pg%jmin_.and.j.le.self%pg%jmax_.and.&
                  &   k.ge.self%pg%kmin_.and.k.le.self%pg%kmax_) self%n_=self%n_+1
               end if
            end do
         end do
      end do
      
      ! Create sub-communicator
      if (self%no_.eq.0) then
         color=MPI_UNDEFINED
         key=0
         self%amIn=.false.
      else
         color=0
         key=0
         self%amIn=.true.
      end if
      call MPI_COMM_SPLIT(self%pg%comm,color,key,self%comm,ierr)
      
      ! If we are not in the iterator, we are done
      if (.not.self%amIn) then
         self%nproc=0
         self%rank=0
         self%amRoot=.false.
         return
      end if
      
      ! Get nproc, rank, and root
      call MPI_COMM_SIZE(self%comm,self%nproc,ierr)
      call MPI_COMM_RANK(self%comm,self%rank,ierr)
      self%amRoot=(self%rank.eq.0)
      
      ! Create unstructured mapping to itr cells - first inside cells then overlap
      allocate(self%map(1:3,1:self%no_))
      cnt=0
      do k=self%pg%kmin_,self%pg%kmax_
         do j=self%pg%jmin_,self%pg%jmax_
            do i=self%pg%imin_,self%pg%imax_
               if (test_if_inside(self%pg,i,j,k)) then
                  cnt=cnt+1
                  self%map(1:3,cnt)=[i,j,k]
               end if
            end do
         end do
      end do
      do k=self%pg%kmino_,self%pg%kmaxo_
         do j=self%pg%jmino_,self%pg%jmaxo_
            do i=self%pg%imino_,self%pg%imaxo_
               ! Skip inside cells
               if (i.ge.self%pg%imin_.and.i.le.self%pg%imax_.and.&
               &   j.ge.self%pg%jmin_.and.j.le.self%pg%jmax_.and.&
               &   k.ge.self%pg%kmin_.and.k.le.self%pg%kmax_) cycle
               ! Only consider overlap cells
               if (test_if_inside(self%pg,i,j,k)) then
                  cnt=cnt+1
                  self%map(1:3,cnt)=[i,j,k]
               end if
            end do
         end do
      end do
            
   end function construct_from_function
   
   !> Basic printing of iterator
   subroutine itr_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER,MPI_MIN,MPI_MAX
      implicit none
      class(itr), intent(in) :: this
      integer :: i,ierr,ntot,maxi,mini,maxj,minj,maxk,mink
      
      ! Return if not in
      if (.not.this%amIn) return
      
      ! Do some processing on the iterator
      call MPI_ALLREDUCE(       this%n_       ,ntot,1,MPI_INTEGER,MPI_SUM,this%comm,ierr)
      call MPI_ALLREDUCE(minval(this%map(1,:)),mini,1,MPI_INTEGER,MPI_MIN,this%comm,ierr)
      call MPI_ALLREDUCE(maxval(this%map(1,:)),maxi,1,MPI_INTEGER,MPI_MAX,this%comm,ierr)
      call MPI_ALLREDUCE(minval(this%map(2,:)),minj,1,MPI_INTEGER,MPI_MIN,this%comm,ierr)
      call MPI_ALLREDUCE(maxval(this%map(2,:)),maxj,1,MPI_INTEGER,MPI_MAX,this%comm,ierr)
      call MPI_ALLREDUCE(minval(this%map(3,:)),mink,1,MPI_INTEGER,MPI_MIN,this%comm,ierr)
      call MPI_ALLREDUCE(maxval(this%map(3,:)),maxk,1,MPI_INTEGER,MPI_MAX,this%comm,ierr)
      
      ! Output
      if (this%amRoot) then
         write(output_unit,'("Iterator [",a,"] for pgrid [",a,"]")') trim(this%name),trim(this%pg%name)
         write(output_unit,'(" >>> Interior cells = ",i0)') ntot
         write(output_unit,'(" >>> Index extent = [",i0,",",i0,"]x[",i0,",",i0,"]x[",i0,",",i0,"]")') mini,maxi,minj,maxj,mink,maxk
      end if
      do i=0,this%nproc-1
         ! Block for clean output
         call MPI_BARRIER(this%comm,ierr)
         ! Output info
         if (this%rank.eq.i) then
            write(output_unit,'(" >>> Rank ",i0," local cells = ",i0," local cells with overlap = ",i0)') this%pg%rank,this%n_,this%no_
         end if
      end do
   end subroutine itr_print
   
end module itr_class
