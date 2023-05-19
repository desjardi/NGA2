!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Dpipe
   use lpt_class,         only: lpt
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(lpt),         public :: lp
   type(incomp),      public :: fs
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,partfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP) :: visc,Uin,dpart
   
   
contains
   
   
   !> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1.and.sqrt(pg%ym(j)**2+pg%zm(k)**2).le.0.5_WP*Dpipe) isIn=.true.
   end function right_boundary


   !> Function that localizes the left domain boundary
   function left_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_boundary
   

   !> Initialization of problem driver
   subroutine simulation_init
      implicit none
      

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         use param, only: param_read
         use incomp_class, only: clipped_neumann,dirichlet,bcond
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Inflow on the left
         call fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=left_boundary)
         ! Outflow on the right
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary)
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero out velocity
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Setup inflow velocity
         call param_read('Inflow velocity',Uin)
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=cfg%VF(i,j,k)*Uin
         end do
         ! Correct MFR
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block create_flow_solver
      
      
      ! Initialize our LPT
      initialize_lpt: block
         use param, only: param_read
         use random, only: random_uniform
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Initialize with zero particles
         call lp%resize(0)
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Injection parameters
         lp%inj_d=Dpipe
         lp%inj_pos=[0.0_WP,0.0_WP,0.0_WP]
         lp%inj_vel=[Uin,0.0_WP,0.0_WP]
         call param_read('Particle diameter',dpart)
         call param_read('Particle mass flow rate',lp%mfr)
      end block initialize_lpt
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=1,name='lpt')
         pmesh%varname(1)='diameter'
         pmesh%vecname(1)='velocity'
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(1,i)=lp%p(i)%d
            pmesh%vec(:,1,i)=lp%p(i)%vel
         end do
      end block create_pmesh
      
      
      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pipe')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('levelset',cfg%Gib)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call lp%get_max()
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
         ! LPT monitor file
         partfile=monitor(amroot=lp%cfg%amRoot,name='particles')
         call partfile%add_column(time%n,'Timestep number')
         call partfile%add_column(time%t,'Time')
         call partfile%add_column(time%dt,'Timestep size')
         call partfile%add_column(lp%np,'Particle number')
         call partfile%add_column(lp%np_new,'Npart new')
         call partfile%add_column(lp%np_out,'Npart removed')
         call partfile%add_column(lp%VFmax,'Max VF')
         call partfile%add_column(lp%Umin,'Particle Umin')
         call partfile%add_column(lp%Umax,'Particle Umax')
         call partfile%add_column(lp%Vmin,'Particle Vmin')
         call partfile%add_column(lp%Vmax,'Particle Vmax')
         call partfile%add_column(lp%Wmin,'Particle Wmin')
         call partfile%add_column(lp%Wmax,'Particle Wmax')
         call partfile%add_column(lp%dmin,'Particle dmin')
         call partfile%add_column(lp%dmax,'Particle dmax')
         call partfile%write()
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
         
         ! Inject particles
         call inject(dt=time%dt)
         
         ! Advance particles by dt
         resU=fs%rho
         call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=fs%visc)
         
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
            
            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                        fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                        fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
            end block ibforcing
           
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
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
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(1,i)=lp%p(i)%d
                  pmesh%vec(:,1,i)=lp%p(i)%vel
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call lp%get_max()
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call partfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Inject particles from a prescribed location with given mass flowrate
   !> Requires injection parameters to be set beforehand
   subroutine inject(dt)
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      use random,    only: random_uniform
      implicit none
      real(WP), intent(inout) :: dt                  !< Timestep size over which to advance
      real(WP) :: Mgoal,Madded,Mtmp,buf              !< Mass flow rate parameters
      real(WP), save :: previous_error=0.0_WP        !< Store mass left over from previous timestep
      integer :: i,np0_,np_tmp,count,ierr
      
      ! Initial number of particles
      np0_=lp%np_
      lp%np_new=0
      
      ! Get the particle mass that should be added to the system
      Mgoal =lp%mfr*dt+previous_error
      Madded=0.0_WP
      
      ! Add new particles until desired mass is achieved
      do while (Madded.lt.Mgoal)
         
         if (lp%cfg%amRoot) then
            ! Initialize parameters
            Mtmp=0.0_WP
            np_tmp=0
            ! Loop while the added volume is not sufficient
            do while (Mtmp.lt.Mgoal-Madded)
               ! Increment counter
               np_tmp=np_tmp+1
               count=np0_+np_tmp
               ! Create space for new particle
               call lp%resize(count)
               ! Assign diameter
               lp%p(count)%d=dpart
               ! Set various parameters for the particle
               lp%p(count)%id    =int(1,8)
               lp%p(count)%dt    =0.0_WP
               lp%p(count)%Acol  =0.0_WP
               lp%p(count)%Tcol  =0.0_WP
               lp%p(count)%angVel=0.0_WP
               ! Give a position at the injector to the particle
               buf=huge(1.0_WP)
               do while (buf.gt.0.5_WP*lp%inj_d)
                  lp%p(count)%pos=[lp%inj_pos(1),&
                  &                random_uniform(lp%inj_pos(2)-0.5_WP*lp%inj_d,lp%inj_pos(2)+0.5_WP*lp%inj_d),&
                  &                random_uniform(lp%inj_pos(3)-0.5_WP*lp%inj_d,lp%inj_pos(3)+0.5_WP*lp%inj_d)]
                  buf=sqrt(lp%p(count)%pos(2)**2+lp%p(count)%pos(3)**2)
               end do
               if (lp%cfg%ny.eq.1)  lp%p(count)%pos(2)=0.0_WP
               if (lp%cfg%nz.eq.1)  lp%p(count)%pos(3)=0.0_WP
               ! Localize the particle
               lp%p(count)%ind=[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin]
               lp%p(count)%ind=lp%cfg%get_ijk_global(lp%p(count)%pos,lp%p(count)%ind)
               ! Give it a velocity
               lp%p(count)%vel=lp%inj_vel
               ! Make it an "official" particle
               lp%p(count)%flag=0
               ! Update the added mass for the timestep
               Mtmp=Mtmp+lp%rho*Pi/6.0_WP*lp%p(count)%d**3
            end do
         end if
         ! Communicate particles
         call lp%sync()
         ! Loop through newly created particles
         buf=0.0_WP
         do i=np0_+1,lp%np_
            ! Remove if out of bounds
            if (lp%cfg%VF(lp%p(i)%ind(1),lp%p(i)%ind(2),lp%p(i)%ind(3)).le.0.0_WP) lp%p(i)%flag=1
            if (lp%p(i)%flag.eq.0) then
               ! Update the added mass for the timestep
               buf=buf+lp%rho*Pi/6.0_WP*lp%p(i)%d**3
               ! Increment counter
               lp%np_new=lp%np_new+1
            end if
         end do
         ! Total mass added
         call MPI_ALLREDUCE(buf,Mtmp,1,MPI_REAL_WP,MPI_SUM,lp%cfg%comm,ierr); Madded=Madded+Mtmp
         ! Clean up particles
         call lp%recycle()
         ! Update initial npart
         np0_=lp%np_
      end do
      
      ! Remember the error
      previous_error=Mgoal-Madded

      ! Sum up injected particles
      call MPI_ALLREDUCE(lp%np_new,i,1,MPI_INTEGER,MPI_SUM,lp%cfg%comm,ierr); lp%np_new=i
      
   end subroutine inject
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
end module simulation
