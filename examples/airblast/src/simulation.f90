!> Various definitions and tools for running an NGA2 simulation
module simulation
   use postproc_class, only: postproc
   use nozzle_class,   only: nozzle
   use atom_class,     only: atom
   use coupler_class,  only: coupler
   implicit none
   private
   
   !> Injector simulation
   type(nozzle) :: injector
   
   !> Atomization simulation
   type(atom) :: atomization
   
   !> Couplers from injector to atomization
   type(coupler) :: xcpl_i2a,ycpl_i2a,zcpl_i2a
   
   !> Postprocessing tool
   !type(postproc) :: pproc
   !logical :: only_pproc=.false.
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Postproc handling
      !call pproc%analyze(only_pproc)
      !if (only_pproc) return
      
      ! Initialize injector simulation
      call injector%init()
      
      ! Initialize atomization simulation
      call atomization%init()
      
      ! If restarting, the domains could be out of sync, so resync
      ! time by forcing injector to be at same time as atomization
      injector%time%t=atomization%time%t
      
      ! Initialize couplers from injector to atomization
      create_coupler_i2a: block
         use parallel, only: group
         xcpl_i2a=coupler(src_grp=group,dst_grp=group,name='nozzle_to_atom'); call xcpl_i2a%set_src(injector%cfg,'x'); call xcpl_i2a%set_dst(atomization%cfg,'x'); call xcpl_i2a%initialize()
         ycpl_i2a=coupler(src_grp=group,dst_grp=group,name='nozzle_to_atom'); call ycpl_i2a%set_src(injector%cfg,'y'); call ycpl_i2a%set_dst(atomization%cfg,'y'); call ycpl_i2a%initialize()
         zcpl_i2a=coupler(src_grp=group,dst_grp=group,name='nozzle_to_atom'); call zcpl_i2a%set_src(injector%cfg,'z'); call zcpl_i2a%set_dst(atomization%cfg,'z'); call zcpl_i2a%initialize()
      end block create_coupler_i2a
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Postproc handling
      !if (only_pproc) return
      
      ! Atomization drives overall time integration
      do while (.not.atomization%time%done())
         
         ! Advance injector simulation until it's caught up
         do while (injector%time%t.le.atomization%time%t)
            call injector%step()
         end do
         
         ! Handle coupling between injector and atomization
         coupling_i2a: block
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cpl12x/y/z couplers
            call xcpl_i2a%push(injector%fs%U); call xcpl_i2a%transfer(); call xcpl_i2a%pull(atomization%resU)
            call ycpl_i2a%push(injector%fs%V); call ycpl_i2a%transfer(); call ycpl_i2a%pull(atomization%resV)
            call zcpl_i2a%push(injector%fs%W); call zcpl_i2a%transfer(); call zcpl_i2a%pull(atomization%resW)
            ! Apply time-varying Dirichlet conditions
            call atomization%fs%get_bcond('gas_inlet',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               atomization%fs%U(i  ,j,k)=atomization%resU(i  ,j,k)*sum(atomization%fs%itpr_x(:,i  ,j,k)*atomization%cfg%VF(i-1:i,    j,    k))
               atomization%fs%V(i-1,j,k)=atomization%resV(i-1,j,k)*sum(atomization%fs%itpr_y(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j-1:j,    k))
               atomization%fs%W(i-1,j,k)=atomization%resW(i-1,j,k)*sum(atomization%fs%itpr_z(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j    ,k-1:k))
            end do
         end block coupling_i2a
         
         ! Advance atomization simulation
         call atomization%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Postproc handling
      !if (only_pproc) return

      ! Finalize injector simulation
      call injector%final()
      
      ! Finalize atomization simulation
      call atomization%final()
      
   end subroutine simulation_final
   

end module simulation
