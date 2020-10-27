!> Various definitions and tools for running an NGA2 simulation
module simulation
   use incomp_class,  only: incomp
   use geometry,      only: cfg
   use ensight_class, only: ensight
   implicit none
   private
   
   !> Single flow solver
   type(incomp), public :: fs
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   
   public :: simulation_init,simulation_run
   
contains
   
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param,     only: param_read
      use precision, only: WP
      use ils_class, only: rbgs
      use ensight_class, only: ensight
      implicit none
      
      ! Create an incompressible flow solver
      create_solver: block
         ! Create solver
         fs=incomp(cfg,'Bob')
         ! Assign constant fluid properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',fs%visc)
         ! Configure pressure solver
         fs%psolv%maxit=100
         fs%psolv%acvg=1.0e-4_WP
         fs%psolv%rcvg=1.0e-4_WP
         fs%psolv%method=rbgs
         ! Check solver objects
         call fs%print()
         call fs%psolv%print()
      end block create_solver
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg,'test')
         ! Add variables to output
         call ens_out%add_scalar('P',fs%P)
         ! Try outputting ensight data
         call ens_out%write_data(0.0_WP)
      end block create_ensight
      
      
      ! Try to use the pressure solver
      test_pressure_solver: block
         real(WP), dimension(:,:,:), allocatable :: rhs
         ! Set Neumann top and bottom
         
         ! Create a RHS
         allocate(rhs(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         
         !call fs%psolv%set_rhs(rhs)
         ! Call the solver
         !call fs%psolv%solve(fs%P)
         
         
      end block test_pressure_solver
      
      
   end subroutine simulation_init
   
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      
      
      
   end subroutine simulation_run
   
   
end module simulation
