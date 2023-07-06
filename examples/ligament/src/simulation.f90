!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,      only: WP
   use hit_class,      only: hit
   use ligament_class, only: ligament
   implicit none
   private
   
   !> HIT simulation
   type(hit) :: turb
   logical :: isInHITGrp
   
   !> Ligament atomization simulation
   type(ligament) :: atom
   
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
      
      ! Prepare HIT simulation
      if (isInHITGrp) then
         prepare_hit: block
            real(WP) :: dt
            ! Initialize HIT
            call turb%init(group=hit_group)
            ! Run HIT until t/tau_eddy=20
            dt=0.15_WP*turb%cfg%min_meshsize/turb%Urms_tgt !< Estimate maximum stable dt
            !do while (turb%time%t.lt.20.0_WP*turb%tau_tgt); call turb%step(dt); end do
         end block prepare_hit

         print*,'rank in turb:',turb%cfg%rank,'loc in turb:',turb%cfg%iproc,turb%cfg%jproc,turb%cfg%kproc,'rank in atom:',atom%cfg%rank,'loc in atom:',atom%cfg%iproc,atom%cfg%jproc,atom%cfg%kproc
      end if
      
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      do while (.not.atom%time%done())
         ! Advance atomization simulation
         call atom%step()
         ! Advance HIT simulation and transfer velocity info
         if (isInHITGrp) then
            ! Advance HIT with maximum stable dt until caught up
            advance_hit: block
               real(WP) :: dt
               dt=0.15_WP*turb%cfg%min_meshsize/turb%Urms_tgt
               !do while (turb%time%t.lt.atom%time%t+20.0_WP*turb%tau_tgt); call turb%step(dt); end do
            end block advance_hit
            ! Transfer turbulent velocity from hit to atomization region
            apply_boundary_condition: block
               use tpns_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,ihit
               real(WP) :: rescaling,tinterp
               rescaling=turb%ti/turb%Urms_tgt
               tinterp=(turb%time%t-(atom%time%t+20.0_WP*turb%tau_tgt))/(turb%time%t-turb%time%told)
               call atom%fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%n_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  ihit=i-atom%fs%cfg%imin+turb%fs%cfg%imax+1
                  if (j.lt.turb%cfg%jmin_.or.j.gt.turb%cfg%jmax_.or.k.lt.turb%cfg%kmin_.or.k.gt.turb%cfg%kmax_) then
                     print*,'ranks:',turb%cfg%rank,atom%cfg%rank,'i,j,k',i,j,k,ihit,'range:',turb%cfg%jmin_,turb%cfg%jmax_,turb%cfg%kmin_,turb%cfg%kmax_
                  end if
                  atom%fs%U(i  ,j,k)=1.0_WP+rescaling*((1.0_WP-tinterp)*turb%fs%U(ihit  ,j,k)+tinterp*turb%fs%Uold(ihit  ,j,k))
                  atom%fs%V(i-1,j,k)=       rescaling*((1.0_WP-tinterp)*turb%fs%V(ihit-1,j,k)+tinterp*turb%fs%Vold(ihit-1,j,k))
                  atom%fs%W(i-1,j,k)=       rescaling*((1.0_WP-tinterp)*turb%fs%W(ihit-1,j,k)+tinterp*turb%fs%Wold(ihit-1,j,k))
               end do
            end block apply_boundary_condition
         end if
         ! Communicate atom velocity field
         call atom%cfg%sync(atom%fs%U)
         call atom%cfg%sync(atom%fs%V)
         call atom%cfg%sync(atom%fs%W)
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