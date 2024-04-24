!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use geometry,      only: cfg1,grp1,isInGrp1
   use geometry,      only: cfg2,grp2,isInGrp2
   use ensight_class, only: ensight
   use coupler_class, only: coupler
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Ensight postprocessing
   type(ensight) :: ens1,ens2
   
   !> Give ourselves two work arrays
   real(WP), dimension(:,:,:), allocatable :: U1,U2
   
   !> Coupler
   type(coupler) :: cpl
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      
      ! Both groups prepare the coupler
      if (isInGrp1.or.isInGrp2) then
         ! Create the coupler
         cpl=coupler(src_grp=grp1,dst_grp=grp2,name='test')
         ! Set the grids
         if (isInGrp1) call cpl%set_src(cfg1)
         if (isInGrp2) call cpl%set_dst(cfg2)
         ! Initialize the metrics
         call cpl%initialize()
      end if
      
      
      ! Group1 does its initialization work here
      if (isInGrp1) then
         
         ! Group1 allocates and initializes the field to transfer
         create_field1: block
            use mathtools, only: twoPi
            integer :: i,j,k
            ! Allocate array
            allocate(U1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
            ! Initialize to solid body rotation
            do k=cfg1%kmino_,cfg1%kmaxo_
               do j=cfg1%jmino_,cfg1%jmaxo_
                  do i=cfg1%imino_,cfg1%imaxo_
                     U1(i,j,k)=-twoPi*cfg1%ym(j)
                  end do
               end do
            end do
         end block create_field1
         
         ! Group1 also outputs to Ensight
         ens1=ensight(cfg=cfg1,name='grid1')
         call ens1%add_scalar('U',U1)
         call ens1%write_data(0.0_WP)
         
      end if
      
      
      ! Group2 allocates the field to receive
      if (isInGrp2) then
         allocate(U2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         U2=-1.0_WP
      end if
      
      
      ! Both groups work on the coupling
      coupling_step: block
         if (isInGrp1) call cpl%push(U1)
         call cpl%transfer()
         if (isInGrp2) call cpl%pull(U2)
      end block coupling_step
      
      
      ! Group2 outputs its received data
      if (isInGrp2) then
         ens2=ensight(cfg=cfg2,name='grid2')
         call ens2%add_scalar('U',U2)
         !call ens2%add_scalar('overlap',cpl%overlap)
         call ens2%write_data(0.0_WP)
      end if
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      ! Nothing here
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Nothing here
   end subroutine simulation_final
   
end module simulation
