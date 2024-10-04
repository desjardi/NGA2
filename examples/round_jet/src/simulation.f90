!> Various definitions and tools for running an NGA2 simulation
module simulation
   use periodicpipe_class, only: periodicpipe
   use roundjet_class,     only: roundjet
   use coupler_class,      only: coupler
   implicit none
   private
   
   !> Periodic pipe simulation
   type(periodicpipe) :: pipe
   
   !> Round jet simulation
   type(roundjet) :: jet
   
   !> Couplers from pipe to round jet
   type(coupler) :: xcpl,ycpl,zcpl
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize periodic pipe simulation
      call pipe%init()
      
      ! Initialize round jet simulation
      call jet%init()
      
      ! If restarting, the domains could be out of sync, so resync
      ! time by forcing pipe to be at same time as jet
      pipe%time%t=jet%time%t
      
      ! Initialize couplers from injector to atomization
      create_coupler: block
         use parallel, only: group
         xcpl=coupler(src_grp=group,dst_grp=group,name='pipe_to_jet'); call xcpl%set_src(pipe%cfg,'x'); call xcpl%set_dst(jet%cfg,'x'); call xcpl%initialize()
         ycpl=coupler(src_grp=group,dst_grp=group,name='pipe_to_jet'); call ycpl%set_src(pipe%cfg,'y'); call ycpl%set_dst(jet%cfg,'y'); call ycpl%initialize()
         zcpl=coupler(src_grp=group,dst_grp=group,name='pipe_to_jet'); call zcpl%set_src(pipe%cfg,'z'); call zcpl%set_dst(jet%cfg,'z'); call zcpl%initialize()
      end block create_coupler
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Jet drives overall time integration
      do while (.not.jet%time%done())
         
         ! Advance pipe simulation until it's caught up
         do while (pipe%time%t.le.jet%time%t)
            call pipe%step()
         end do
         
         ! Handle coupling between pipe and jet
         coupling: block
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cpl12x/y/z couplers
            call xcpl%push(pipe%fs%U); call xcpl%transfer(); call xcpl%pull(jet%resU)
            call ycpl%push(pipe%fs%V); call ycpl%transfer(); call ycpl%pull(jet%resV)
            call zcpl%push(pipe%fs%W); call zcpl%transfer(); call zcpl%pull(jet%resW)
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
      call pipe%final()
      
      ! Finalize jet simulation
      call jet%final()
      
   end subroutine simulation_final
   
   
end module simulation
