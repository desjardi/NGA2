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
      ! This is our communication information
      type(MPI_Comm) :: comm                              !< Intracommunicator over the union of both groups
      type(MPI_Group) :: sgrp,dgrp,grp                    !< Source and destination groups and their union
      integer :: nproc                                    !< Number of processors
      integer :: rank                                     !< Processor grid rank
      logical :: amRoot                                   !< Am I root for the coupler?
   contains
      
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
      self%amRoot=(self%rank.eq.0)
      
      ! Create intracommunicator for the new group
      call MPI_COMM_CREATE_GROUP(comm,self%grp,0,self%comm,ierr)
      
      
      
      
      
      ! Log/screen output
      !logging: block
      !   use, intrinsic :: iso_fortran_env, only: output_unit
      !   use param,    only: verbose
      !   use messager, only: log
      !   use string,   only: str_long
      !   character(len=str_long) :: message
      !   if (self%src%amRoot) then
      !      write(message,'("Coupler [",a,"] from pgrid [",a,"] to pgrid [",a,"]")') trim(self%name),trim(self%src%name),trim(self%dst%name)
      !      if (verbose.gt.1) write(output_unit,'(a)') trim(message)
      !      if (verbose.gt.0) call log(message)
      !   end if
      !end block logging
      
   end function construct_from_two_groups
   
   
end module coupler_class
