!> Various definitions and tools for running an NGA2 simulation
module simulation
   use periodicpipe_class, only: periodicpipe
   use roundjet_class,     only: roundjet
   use coupler_class,      only: coupler
   implicit none
   private
   
   !> Periodic pipe simulation
   logical :: in_pipe_group
   type(periodicpipe) :: pipe
   
   !> Round jet simulation
   type(roundjet) :: jet
   
   !> Couplers from pipe to round jet
   type(coupler) :: xcpl,ycpl,zcpl
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: pipe_group
      
      ! Create an MPI group for the pipe simulation
      create_pipe_group: block
         use parallel, only: group
         ! Assume that all cores work on the pipe
         in_pipe_group=.true.
         pipe_group=group
      end block create_pipe_group
      
      ! Initialize periodic pipe simulation
      if (in_pipe_group) call pipe%init(group=pipe_group)
      
      ! Initialize round jet simulation
      call jet%init()
      
      ! If restarting, the domains could be out of sync, so resync
      ! time by forcing pipe to be at same time as jet
      pipe%time%t=jet%time%t
      
      ! Initialize couplers from injector to atomization
      create_coupler: block
         use parallel, only: group
         xcpl=coupler(src_grp=pipe_group,dst_grp=group,name='pipe2jet')
         ycpl=coupler(src_grp=pipe_group,dst_grp=group,name='pipe2jet')
         zcpl=coupler(src_grp=pipe_group,dst_grp=group,name='pipe2jet')
         if (in_pipe_group) call xcpl%set_src(pipe%cfg,'x')
         if (in_pipe_group) call ycpl%set_src(pipe%cfg,'y')
         if (in_pipe_group) call zcpl%set_src(pipe%cfg,'z')
         call xcpl%set_dst(jet%cfg,'x'); call xcpl%initialize()
         call ycpl%set_dst(jet%cfg,'y'); call ycpl%initialize()
         call zcpl%set_dst(jet%cfg,'z'); call zcpl%initialize()
      end block create_coupler
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Jet drives overall time integration
      do while (.not.jet%time%done())
         
         ! Advance pipe simulation until it's caught up
         if (in_pipe_group) then
            do while (pipe%time%t.le.jet%time%t)
               call pipe%step()
            end do
         end if
         
         ! Handle coupling between pipe and jet
         coupling: block
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cpl12x/y/z couplers
            if (in_pipe_group) call xcpl%push(pipe%fs%U)
            if (in_pipe_group) call ycpl%push(pipe%fs%V)
            if (in_pipe_group) call zcpl%push(pipe%fs%W)
            call xcpl%transfer(); call xcpl%pull(jet%resU)
            call ycpl%transfer(); call ycpl%pull(jet%resV)
            call zcpl%transfer(); call zcpl%pull(jet%resW)
            ! Apply time-varying Dirichlet conditions
            call jet%fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               jet%fs%U(i  ,j,k)=jet%resU(i  ,j,k)
               jet%fs%V(i-1,j,k)=jet%resV(i-1,j,k)
               jet%fs%W(i-1,j,k)=jet%resW(i-1,j,k)
            end do
         end block coupling
         
         ! Advance jet simulation
         call jet%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize pipe simulation
      if (in_pipe_group) call pipe%final()
      
      ! Finalize jet simulation
      call jet%final()
      
   end subroutine simulation_final
   
   
end module simulation
