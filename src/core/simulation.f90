!> Various definitions and tools for running an NGA2 simulation
module simulation
   use incomp_class,  only: incomp,dirichlet
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
   
   !> Function that localizes the top of the domain
   function yplus_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yplus_locator
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param,     only: param_read
      use precision, only: WP
      use ils_class, only: rbgs,amg
      !use ensight_class, only: ensight
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
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         fs%U=0.0_WP
         fs%V=0.0_WP
         fs%W=0.0_WP
      end block initialize_velocity
      
      
      ! Initialize boundary conditions
      initialize_bc: block
         call fs%add_bcond('stokes',dirichlet,yplus_locator)
      end block initialize_bc
      
      
   end subroutine simulation_init
   
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      !logical :: done
      
      ! Implement an Euler explicit integration of NS
      !do while (.not.done)
         
         
         
         
      !end do
      
   end subroutine simulation_run
   
   
end module simulation
