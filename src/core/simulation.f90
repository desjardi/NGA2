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
      use ils_class, only: rbgs,amg
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
         ! Initialize solver
         call fs%psolv%init_solver(amg)
         ! Check solver objects
         call fs%print()
      end block create_solver
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg,'test')
         ! Add variables to output
         call ens_out%add_scalar('Pressure',fs%P)
         call ens_out%add_vector('Velocity',fs%U,fs%V,fs%W)
         ! Try outputting ensight data
         !call ens_out%write_data(0.0_WP)
      end block create_ensight
      
      
      ! Try to use the pressure solver
      test_pressure_solver: block
         ! Create a scaled RHS and output it
         fs%psolv%rhs=0.0_WP
         if (fs%cfg%jproc.eq.         1) fs%psolv%rhs(:,fs%cfg%jmin_,:)=+1.0_WP
         if (fs%cfg%jproc.eq.fs%cfg%npy) fs%psolv%rhs(:,fs%cfg%jmax_,:)=-1.0_WP
         fs%psolv%rhs=-fs%cfg%vol*fs%psolv%rhs
         call ens_out%add_scalar('RHS',fs%psolv%rhs)
         ! Set initial guess to zero
         fs%psolv%sol=0.0_WP
         ! Call the solver
         call fs%psolv%solve()
         ! Copy back to pressure
         fs%P=fs%psolv%sol
         call fs%psolv%print()
         ! Output to ensight
         call ens_out%write_data(0.0_WP)
      end block test_pressure_solver
      
      
   end subroutine simulation_init
   
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      
      
      
   end subroutine simulation_run
   
   
end module simulation
