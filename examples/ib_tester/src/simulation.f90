!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use df_class,          only: dfibm
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a direct forcing solver, an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp),      public :: fs
   type(dfibm),       public :: df
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,ibmfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP) :: inlet_velocity

contains
   
   
   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg,gmres_pfmg
         use incomp_class, only: dirichlet,clipped_neumann
         real(WP) :: visc
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='inflow', type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=gmres_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our direct forcing solver
      initialize_df: block
         use mathtools, only: twoPi,arctan
         integer :: i,j,k,np
         real(WP) :: Dcyl,Xcyl,amp,freq,theta,r
         ! Create solver
         df=dfibm(cfg=cfg,name='IBM')
         ! Read cylinder properties
         call param_read('Number of markers',np)
         call param_read('Diameter',Dcyl)
         call param_read('Position',Xcyl)
         call param_read('Perturbation amp',amp,default=0.05_WP*Dcyl)
         call param_read('Perturbation freq',freq,default=0.0_WP*Dcyl)
         ! Root process initializes marker particles
         if (df%cfg%amRoot) then
            call df%resize(np)
            ! Distribute marker particles
            do i=1,np
               ! Set various parameters for the marker
               df%p(i)%id  =1
               df%p(i)%vel =0.0_WP
               ! Set position
               theta=real(i-1,WP)*twoPi/real(np,WP)
               r=0.5_WP*Dcyl+amp*sin(freq*theta)
               df%p(i)%pos(1)=r*cos(theta)
               df%p(i)%pos(2)=r*sin(theta)
               df%p(i)%pos(3)=0.0_WP
               ! Assign element area
               df%p(i)%dA=twoPi*r/real(np,WP)*df%cfg%zL
               ! Assign outward normal vector
               df%p(i)%norm = df%p(i)%pos/r
               ! Shift cylinder
               df%p(i)%pos(1)=df%p(i)%pos(1)+Xcyl
               ! Locate the particle on the mesh
               df%p(i)%ind=df%cfg%get_ijk_global(df%p(i)%pos,[df%cfg%imin,df%cfg%jmin,df%cfg%kmin])
               ! Activate the particle
               df%p(i)%flag=0
            end do
         end if
         call df%sync()
         
         ! Get initial volume fraction
         call df%update_VF()
         
         ! All processes initialize IBM objects
         call df%setup_obj()
         
         ! Define levelset (only used for visualization)
         do k=df%cfg%kmin_,df%cfg%kmax_
            do j=df%cfg%jmin_,df%cfg%jmax_
               do i=df%cfg%imin_,df%cfg%imax_
                  theta=arctan(df%cfg%xm(i),df%cfg%ym(j))
                  r=0.5_WP*Dcyl+amp*sin(freq*theta)
                  df%G(i,j,k)=sqrt((df%cfg%xm(i)-Xcyl)**2+df%cfg%ym(j)**2)-r
               end do
            end do
         end do
         call df%cfg%sync(df%G)
         
         if (df%cfg%amRoot) then
            print*,"===== Direct Forcing Setup Description ====="
            print*,'Number of marker particles', df%np
            print*,'Number of IBM objects', df%nobj
            print*,'Particle spacing / dx', twoPi*r/real(np,WP)/df%cfg%xL*real(df%cfg%nx,WP)
         end if
      end block initialize_df
      
      
      ! Create partmesh object for marker particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=1,name='ibm')
         pmesh%varname(1)='area'
         pmesh%vecname(1)='velocity'
         call df%update_partmesh(pmesh)
         do i=1,df%np_
            pmesh%var(1,i)=df%p(i)%dA
            pmesh%vec(:,1,i)=df%p(i)%vel
         end do
      end block create_pmesh
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Initial fields
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! Set inflow velocity/momentum
         call param_read('Inlet velocity',inlet_velocity)
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=inlet_velocity
         end do
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='ibm')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('markers',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VF',df%VF)
         call ens_out%add_scalar('levelset',df%G)
         call ens_out%add_scalar('pressure',fs%P)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call df%get_max()
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
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
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
         ! Create IBM monitor
         ibmfile=monitor(amroot=df%cfg%amRoot,name='ibm')
         call ibmfile%add_column(time%n,'Timestep number')
         call ibmfile%add_column(time%t,'Time')
         call ibmfile%add_column(df%VFmin,'VF min')
         call ibmfile%add_column(df%VFmax,'VF max')
         call ibmfile%add_column(df%Fx,'Fx')
         call ibmfile%add_column(df%Fy,'Fy')
         call ibmfile%add_column(df%Fz,'Fz')
         call ibmfile%add_column(df%np,'Marker number')
         call ibmfile%add_column(df%Umin,'Marker Umin')
         call ibmfile%add_column(df%Umax,'Marker Umax')
         call ibmfile%add_column(df%Vmin,'Marker Vmin')
         call ibmfile%add_column(df%Vmax,'Marker Vmax')
         call ibmfile%add_column(df%Wmin,'Marker Wmin')
         call ibmfile%add_column(df%Wmax,'Marker Wmax')
         call ibmfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Add momentum source term from direct forcing
            ibm_correction: block
               integer :: i,j,k
               resU=fs%rho
               call df%get_source(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU)
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%U(i,j,k)=fs%U(i,j,k)+sum(fs%itpr_x(:,i,j,k)*df%srcU(i-1:i,j,k))
                        fs%V(i,j,k)=fs%V(i,j,k)+sum(fs%itpr_y(:,i,j,k)*df%srcV(i,j-1:j,k))
                        fs%W(i,j,k)=fs%W(i,j,k)+sum(fs%itpr_z(:,i,j,k)*df%srcW(i,j,k-1:k))
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
            end block ibm_correction
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Reset Dirichlet BCs
            dirichlet_velocity: block
               use incomp_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)=inlet_velocity
               end do
            end block dirichlet_velocity
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call df%update_partmesh(pmesh)
               do i=1,df%np_
                  pmesh%var(1,i)=df%p(i)%dA
                  pmesh%vec(:,1,i)=df%p(i)%vel
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call df%get_max()
         call mfile%write()
         call cflfile%write()
         call ibmfile%write()
         
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
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
end module simulation
