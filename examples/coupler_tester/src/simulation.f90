!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use geometry,      only: cfg1,cfg2
   use geometry,      only: isInGroup1,isInGroup2
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
      
      ! Group1's initial work
      if (isInGroup1) then
         
         ! Group1 prepares the field to transfer
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
         create_ensight1: block
            ens1=ensight(cfg=cfg1,name='grid1')
            call ens1%add_scalar('U',U1)
            call ens1%write_data(0.0_WP)
         end block create_ensight1
         
      end if
      
      ! Group2's initial work
      if (isInGroup2) then
         
         ! Group1 prepares the field to receive
         create_field2: block
            allocate(U2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); U2=-1.0_WP
         end block create_field2
         
         ! Group1 also outputs to Ensight
         create_ensight2: block
            ens2=ensight(cfg=cfg2,name='grid2')
            call ens2%add_scalar('U',U2)
            call ens2%write_data(0.0_WP)
         end block create_ensight2
         
      end if
      
      
      ! Now we create the coupler
      
      
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
