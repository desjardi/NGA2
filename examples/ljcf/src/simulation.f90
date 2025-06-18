!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,      only: WP
   use hit_class,      only: hit
   use ljcf_class,     only: ljcf
   use coupler_class,  only: coupler
   implicit none
   private
   
   !> HIT simulation
   type(hit) :: turb
   logical :: isInHITGrp
   
   !> LJCF atomization simulation
   type(ljcf) :: atom
   
   !> Coupler from turb to atom
   type(coupler) :: xcpl,ycpl,zcpl
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: hit_group
      
      ! Initialize atomization simulation
      call atom%init()
      
      ! Create an MPI group using leftmost processors only
      create_hit_group: block
         use parallel, only: group,comm
         use mpi_f08,  only: MPI_Group_incl
         integer, dimension(:), allocatable :: ranks
         integer, dimension(3) :: coord
         integer :: n,ngrp,ierr,ny,nz
         ngrp=atom%cfg%npy*atom%cfg%npz
         allocate(ranks(ngrp))
         ngrp=0
         do nz=1,atom%cfg%npz
            do ny=1,atom%cfg%npy
               ngrp=ngrp+1
               coord=[0,ny-1,nz-1]
               call MPI_CART_RANK(atom%cfg%comm,coord,ranks(ngrp),ierr)
            end do
         end do
         call MPI_Group_incl(group,ngrp,ranks,hit_group,ierr)
         if (atom%cfg%iproc.eq.1) then
            isInHITGrp=.true.
         else
            isInHITGrp=.false.
         end if
      end block create_hit_group
      
      ! Initialize HIT simulation
      if (isInHITGrp) call turb%init(group=hit_group,xend=atom%cfg%x(atom%cfg%imin))
      
      ! If restarting, the domains could be out of sync, so resync
      ! time by forcing HIT to be at same time as jet
      if (isInHITGrp) then
         turb%time%t=atom%time%t
         turb%time%told=turb%time%t-turb%time%dt
      end if
      
      ! Initialize couplers from turb to atom
      create_coupler: block
         use parallel, only: group
         xcpl=coupler(src_grp=hit_group,dst_grp=group,name='turb2atom')
         ycpl=coupler(src_grp=hit_group,dst_grp=group,name='turb2atom')
         zcpl=coupler(src_grp=hit_group,dst_grp=group,name='turb2atom')
         if (isInHITGrp) call xcpl%set_src(turb%cfg,'x')
         if (isInHITGrp) call ycpl%set_src(turb%cfg,'y')
         if (isInHITGrp) call zcpl%set_src(turb%cfg,'z')
         call xcpl%set_dst(atom%cfg,'x'); call xcpl%initialize()
         call ycpl%set_dst(atom%cfg,'y'); call ycpl%initialize()
         call zcpl%set_dst(atom%cfg,'z'); call zcpl%initialize()
      end block create_coupler
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      do while (.not.atom%time%done())
         
         ! Advance HIT simulation and transfer velocity info
         if (isInHITGrp) then
            ! Advance HIT with maximum stable dt until caught up
            advance_hit: block
               real(WP) :: dt
               dt=0.15_WP*turb%cfg%min_meshsize/turb%Urms_tgt
               do while (turb%time%t.le.atom%time%t)
                  call turb%step(dt)
               end do
            end block advance_hit
         end if
         
         ! Handle coupling between HIT and atomization simulation
         coupling: block
            ! Push data from HIT simulation
            if (isInHITGrp) then
               push_velocity: block
                  real(WP) :: rescaling,tinterp
                  rescaling=turb%ti/turb%Urms_tgt
                  tinterp=(turb%time%t-atom%time%t)/(turb%time%t-turb%time%told)
                  turb%resU=rescaling*((1.0_WP-tinterp)*turb%fs%U+tinterp*turb%fs%Uold); call xcpl%push(turb%resU)
                  turb%resV=rescaling*((1.0_WP-tinterp)*turb%fs%V+tinterp*turb%fs%Vold); call ycpl%push(turb%resV)
                  turb%resW=rescaling*((1.0_WP-tinterp)*turb%fs%W+tinterp*turb%fs%Wold); call zcpl%push(turb%resW)
                  end block push_velocity
            end if
            ! Transfer and pull
            call xcpl%transfer(); call xcpl%pull(atom%resU)
            call ycpl%transfer(); call ycpl%pull(atom%resV)
            call zcpl%transfer(); call zcpl%pull(atom%resW)
            ! Apply time-dependent Dirichlet condition
            apply_boundary_condition: block
               use tpns_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call atom%fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  atom%fs%U(i  ,j,k)=1.0_WP !atom%resU(i  ,j,k)+1.0_WP
                  atom%fs%V(i-1,j,k)=0.0_WP !atom%resV(i-1,j,k)
                  atom%fs%W(i-1,j,k)=0.0_WP !atom%resW(i-1,j,k)
               end do
            end block apply_boundary_condition
         end block coupling
         
         ! Advance atomization simulation
         call atom%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize atomization simulation
      call atom%final()
      
      ! Finalize HIT simulation
      if (isInHITGrp) call turb%final()
      
   end subroutine simulation_final
   
   
end module simulation