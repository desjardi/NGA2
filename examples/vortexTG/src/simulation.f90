!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use matm_class,        only: matm
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver, volume fraction solver, and material model set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,kefile
   
   public :: simulation_init,simulation_run,simulation_final
   
   ! Ad-hoc monitor variable and function
   real(WP) :: KE
   public   :: current_KE
   
contains
  
  function current_KE(this) result(buf)
    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
    use parallel,  only: MPI_REAL_WP
    implicit none
    type(mast), intent(in) :: this
    real(WP) :: buf, KE
    integer  :: ierr,i,j,k
    
    KE = 0.0_WP
    
    do k=this%cfg%kmin_,this%cfg%kmax_
      do j=this%cfg%jmin_,this%cfg%jmax_
        do i=this%cfg%imin_,this%cfg%imax_
          KE = KE + 0.5_WP*this%RHO(i,j,k)* &
          (this%Ui(i,j,k)**2+this%Vi(i,j,k)**2+this%Wi(i,j,k)**2)
        end do
      end do
    end do
    call MPI_ALLREDUCE(KE,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
    
  end function
  
  
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         !call param_read('Max steps',time%nmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Single-phase simulation, but MAST requires vfs solver
      create_and_initialize_vof: block
         use vfs_class, only: lvira
         integer :: i,j,k
         real(WP) :: vol,area
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Single phase
         vf%VF = 0.0_WP
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a compressible two-phase flow solver
      create_and_initialize_flow_solver: block
         use mast_class, only: clipped_neumann,dirichlet,bc_scope,bcond
         use matm_class, only: none
         use ils_class,  only: pcg_bbox,pcg_amg
         use mathtools,  only: Pi
         character(len=1) :: orientation
         integer :: i,j,k,n
         real(WP), dimension(3) :: xyz
         real(WP) :: gamm_l,Pref_l,gamm_g,rho_g0,visc_g,psol,p0
         type(bcond), pointer :: mybc
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get EOS parameters from input
         call param_read('Liquid gamma',gamm_l)
         call param_read('Liquid Pref', Pref_l)
         call param_read('Gas gamma',gamm_g)
         ! Register equations of state
         call matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
         call matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         
         ! Assign constant viscosity to each phase and turn off heat diffusion
         call param_read('Gas dynamic viscosity',visc_g)
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g,hdffmodel_gas=none)
         ! Surface tension should be set to 0, it is a single-phase case
         fs%sigma = 0.0_WP
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_bbox,implicit_ils=pcg_amg)

         ! Initial conditions
         call param_read('Gas density',rho_g0)
         fs%Grho = rho_g0
         ! Liquid quantities are irrelevant, set to unity
         fs%Lrho = 1.0_WP; fs%LrhoE = 1.0_WP
         call param_read('Background pressure',p0)
         ! Orientation of flow field
         if (fs%cfg%nz.eq.1) then
           orientation = 'z'
         elseif (fs%cfg%ny.eq.1) then
           orientation = 'y'
         elseif (fs%cfg%nx.eq.1) then
           orientation = 'x'
         else
           call param_read('Orientation',orientation)
         end if
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
           do j=fs%cfg%jmino_,fs%cfg%jmaxo_
             do i=fs%cfg%imino_,fs%cfg%imaxo_
               ! Set up flow field according to designated orientation
               select case(orientation)
               case('x')
                 fs%Ui(i,j,k) =  0.0_WP
                 fs%Vi(i,j,k) =  cos(fs%cfg%ym(j))*sin(fs%cfg%zm(k))
                 fs%Wi(i,j,k) = -sin(fs%cfg%ym(j))*cos(fs%cfg%zm(k))
                 psol = p0 - 0.25_WP*(cos(2.0_WP*fs%cfg%ym(j))+cos(2.0_WP*fs%cfg%zm(k)))
               case('y')
                 fs%Ui(i,j,k) = -sin(fs%cfg%zm(k))*cos(fs%cfg%xm(i))
                 fs%Vi(i,j,k) =  0.0_WP
                 fs%Wi(i,j,k) =  cos(fs%cfg%zm(k))*sin(fs%cfg%xm(i))
                 psol = p0 - 0.25_WP*(cos(2.0_WP*fs%cfg%zm(k))+cos(2.0_WP*fs%cfg%xm(i)))
               case('z') ! This is the default
                 fs%Ui(i,j,k) = -cos(fs%cfg%xm(i))*sin(fs%cfg%ym(j))
                 fs%Vi(i,j,k) =  sin(fs%cfg%xm(i))*cos(fs%cfg%ym(j))
                 fs%Wi(i,j,k) =  0.0_WP
                 psol = p0 - 0.25_WP*(cos(2.0_WP*fs%cfg%xm(i))+cos(2.0_WP*fs%cfg%ym(j)))
               end select
               ! Calculate total energy
               fs%GrhoE(i,j,k) = matmod%EOS_energy(psol,rho_g0,fs%Ui(i,j,k),fs%Vi(i,j,k),fs%Wi(i,j,k),'gas')
             end do 
           end do
         end do

         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         ! Apply face BC - this is for the sync operations within
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt,bc_scope)
         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         ! Perform initial pressure relax
         call fs%pressure_relax(vf,matmod)
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA
         
         ! Update KE
         KE = current_KE(fs)

      end block create_and_initialize_flow_solver
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='TaylorGreenVortex')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('Density',fs%RHO)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl,matmod=matmod)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create kinetic energy monitor
         kefile=monitor(fs%cfg%amRoot,'kinetic_energy')
         call kefile%add_column(time%n,'Timestep number')
         call kefile%add_column(time%t,'Time')
         call kefile%add_column(KE,'Kinetic energy')
         call kefile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl,matmod=matmod)
         call time%adjust_dt()
         call time%increment()
         
         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()
         
         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Insert viscous step here, or possibly incorporate into predictor above
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)
            
            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            fs%psolv%sol=fs%PA-fs%P
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Single-phase, no need for relaxation
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Update KE
         KE = current_KE(fs)
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call kefile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays - none
      
   end subroutine simulation_final
   
end module simulation
