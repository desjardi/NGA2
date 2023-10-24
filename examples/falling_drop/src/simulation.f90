!> Various definitions and tools for running an NGA2 simulation
module simulation
	use precision,         only: WP
	use geometry,          only: cfg
	use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
	use tpns_class,        only: tpns
	use vfs_class,         only: vfs
   use tpscalar_class,    only: tpscalar
	use timetracker_class, only: timetracker
	use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
	use event_class,       only: event
	use monitor_class,     only: monitor
	implicit none
	private
	
	!> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
	type(hypre_str),   public :: ps
	type(ddadi),       public :: vs
	type(tpns),        public :: fs
	type(vfs),         public :: vf
   type(tpscalar),    public :: sc
	type(timetracker), public :: time
	
	!> Ensight postprocessing
   type(surfmesh) :: smesh
	type(ensight)  :: ens_out
	type(event)    :: ens_evt
	
	!> Simulation monitor file
	type(monitor) :: mfile,cflfile,scfile
	
	public :: simulation_init,simulation_run,simulation_final
	
	!> Private work arrays
	real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
	real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
	
	!> Problem definition
	real(WP), dimension(3) :: center
	real(WP) :: radius,depth
	
contains


   !> Function that defines a level set function for a falling drop problem
	function levelset_falling_drop(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G
		! Create the droplet
	   G=radius-sqrt(sum((xyz-center)**2))
	   ! Add the pool
	   G=max(G,depth-xyz(2))
	end function levelset_falling_drop
	
	
	!> Initialization of problem solver
	subroutine simulation_init
		use param, only: param_read
		implicit none
		
		
		! Allocate work arrays
	   allocate_work_arrays: block
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,2))
		   allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
		end block allocate_work_arrays
		
		
		! Initialize time tracker with 2 subiterations
	   initialize_timetracker: block
		   time=timetracker(amRoot=cfg%amRoot)
			call param_read('Max timestep size',time%dtmax)
			call param_read('Max cfl number',time%cflmax)
			call param_read('Max time',time%tmax)
			time%dt=time%dtmax
			time%itmax=2
		end block initialize_timetracker
		
      
		! Initialize our VOF solver and field
	   create_and_initialize_vof: block
		   use mms_geom, only: cube_refine_vol
		   use vfs_class, only: lvira,VFhi,VFlo
			integer :: i,j,k,n,si,sj,sk
			real(WP), dimension(3,8) :: cube_vertex
			real(WP), dimension(3) :: v_cent,a_cent
			real(WP) :: vol,area
			integer, parameter :: amr_ref_lvl=4
			! Create a VOF solver
		   !vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF',store_detailed_flux=.true.)
		   ! Initialize to a droplet and a pool
		   center=[0.0_WP,0.01_WP,0.0_WP]
		   radius=0.002_WP
		   depth =0.005_WP
			do k=vf%cfg%kmino_,vf%cfg%kmaxo_
				do j=vf%cfg%jmino_,vf%cfg%jmaxo_
					do i=vf%cfg%imino_,vf%cfg%imaxo_
						! Set cube vertices
				      n=0
						do sk=0,1
							do sj=0,1
								do si=0,1
									n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
								end do
							end do
						end do
						! Call adaptive refinement code to get volume and barycenters recursively
				      vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
						call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_falling_drop,0.0_WP,amr_ref_lvl)
						vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
						if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
							vf%Lbary(:,i,j,k)=v_cent
							vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
						else
							vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
							vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
						end if
					end do
				end do
			end do
			! Update the band
		   call vf%update_band()
			! Perform interface reconstruction from VOF field
		   call vf%build_interface()
			! Set interface planes at the boundaries
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
		
		
		! Create a two-phase flow solver without bconds
	   create_and_initialize_flow_solver: block
		   use hypre_str_class, only: pcg_pfmg
         use mathtools,       only: Pi
			! Create flow solver
		   fs=tpns(cfg=cfg,name='Two-phase NS')
		   ! Assign constant viscosity to each phase
		   call param_read('Liquid dynamic viscosity',fs%visc_l)
			call param_read('Gas dynamic viscosity',fs%visc_g)
			! Assign constant density to each phase
		   call param_read('Liquid density',fs%rho_l)
			call param_read('Gas density',fs%rho_g)
			! Read in surface tension coefficient
		   call param_read('Surface tension coefficient',fs%sigma)
         call param_read('Static contact angle',fs%contact_angle)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
		   ! Assign acceleration of gravity
		   call param_read('Gravity',fs%gravity)
			! Configure pressure solver
			ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
		   call fs%setup(pressure_solver=ps,implicit_solver=vs)
		   ! Zero initial field
		   fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
		   ! Calculate cell-centered velocities and divergence
		   call fs%interp_vel(Ui,Vi,Wi)
		   call fs%get_div()
	   end block create_and_initialize_flow_solver
	   
      
      ! Create a liquid scalar solver
      create_scalar: block
         integer :: i,j,k
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=2,name='tpscalar_test')
         ! Make it liquid and give it a name
         sc%SCname=['Zl','Zg']
         sc%phase =[  0 ,  1 ]
         ! Assign zero diffusivity
         sc%diff=0.0_WP
         ! Setup without an implicit solver
         call sc%setup()
         ! Initialize scalar fields
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Liquid scalar
                  if (vf%VF(i,j,k).gt.0.0_WP) then
                     ! We are in the liquid
                     if (cfg%ym(j).gt.depth) then
                        ! We are above the pool
                        sc%SC(i,j,k,1)=1.0_WP
                     else
                        ! We are in the pool
                        sc%SC(i,j,k,1)=2.0_WP
                     end if
                  end if
                  ! Gas scalar
                  if (vf%VF(i,j,k).lt.1.0_WP) then
                     ! We are in the gas
                     sc%SC(i,j,k,2)=(cfg%ym(j)-depth)/(cfg%yL-depth)
                  end if
               end do
            end do
         end do
      end block create_scalar
      
      
	   ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


	   ! Add Ensight output
	   create_ensight: block
         integer :: nsc
		   ! Create Ensight output from cfg
		   ens_out=ensight(cfg=cfg,name='FallingDrop')
			! Create event for Ensight output
		   ens_evt=event(time=time,name='Ensight output')
		   call param_read('Ensight output period',ens_evt%tper)
		   ! Add variables to output
		   call ens_out%add_vector('velocity',Ui,Vi,Wi)
		   call ens_out%add_scalar('VOF',vf%VF)
		   call ens_out%add_scalar('pressure',fs%P)
		   call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         do nsc=1,sc%nscalar
            call ens_out%add_scalar(trim(sc%SCname(nsc)),sc%SC(:,:,:,nsc))
         end do
		   ! Output to ensight
		   if (ens_evt%occurs()) call ens_out%write_data(time%t)
	   end block create_ensight
	   
	   
	   ! Create a monitor file
	   create_monitor: block
         integer :: nsc
		   ! Prepare some info about fields
		   call fs%get_cfl(time%dt,time%cfl)
		   call fs%get_max()
		   call vf%get_max()
         call sc%get_max(VF=vf%VF)
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
		   call mfile%add_column(vf%VFmax,'VOF maximum')
		   call mfile%add_column(vf%VFmin,'VOF minimum')
		   call mfile%add_column(vf%VFint,'VOF integral')
		   call mfile%add_column(fs%divmax,'Maximum divergence')
		   call mfile%add_column(fs%psolv%it,'Pressure iteration')
		   call mfile%add_column(fs%psolv%rerr,'Pressure error')
		   call mfile%write()
		   ! Create CFL monitor
		   cflfile=monitor(fs%cfg%amRoot,'cfl')
		   call cflfile%add_column(time%n,'Timestep number')
		   call cflfile%add_column(time%t,'Time')
		   call cflfile%add_column(fs%CFLst,'STension CFL')
		   call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
		   call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
		   call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
		   call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
		   call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
		   call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
		   call cflfile%write()
         ! Create scalar monitor
         scfile=monitor(sc%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         do nsc=1,sc%nscalar
            call scfile%add_column(sc%SCmin(nsc),trim(sc%SCname(nsc))//'_min')
            call scfile%add_column(sc%SCmax(nsc),trim(sc%SCname(nsc))//'_max')
            call scfile%add_column(sc%SCint(nsc),trim(sc%SCname(nsc))//'_int')
         end do
         call scfile%write()
	   end block create_monitor
	   
	   
	end subroutine simulation_init
	
	
	!> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
	subroutine simulation_run
      use tpns_class, only: static_contact,arithmetic_visc
		implicit none
		
		! Perform time integration
	   do while (.not.time%done())
			
			! Increment time
		   call fs%get_cfl(time%dt,time%cfl)
		   call time%adjust_dt()
			call time%increment()
			
			! Remember old VOF
		   vf%VFold=vf%VF
         
         ! Remember old SC
         sc%SCold=sc%SC
         
			! Remember old velocity
		   fs%Uold=fs%U
			fs%Vold=fs%V
			fs%Wold=fs%W
			
			! Apply time-varying Dirichlet conditions
		   ! This is where time-dpt Dirichlet would be enforced
		   
		   ! Prepare old staggered density (at n)
		   call fs%get_olddensity(vf=vf)
			
			! VOF solver step
		   call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
			
			! Prepare new staggered viscosity (at n+1)
		   call fs%get_viscosity(vf=vf,strat=arithmetic_visc)

         ! Now transport our phase-specific scalars
         advance_scalar: block
            use vfs_class, only: VFlo,VFhi
            integer :: nsc
            ! Explicit calculation of dSC/dt from scalar equation
            call sc%get_dSCdt(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)
            ! Advance scalar fields
            do nsc=1,sc%nscalar
               if (sc%phase(nsc).eq.0) then
                  ! Liquid scalar
                  where (vf%VF.gt.VFlo)
                     sc%SC(:,:,:,nsc)=(vf%VFold*sc%SCold(:,:,:,nsc)+time%dt*resSC(:,:,:,nsc))/vf%VF
                  else where
                     sc%SC(:,:,:,nsc)=0.0_WP
                  end where
               else if (sc%phase(nsc).eq.1) then
                  ! Gas scalar
                  where (vf%VF.lt.VFhi)
                     sc%SC(:,:,:,nsc)=((1.0_WP-vf%VFold)*sc%SCold(:,:,:,nsc)+time%dt*resSC(:,:,:,nsc))/(1.0_WP-vf%VF)
                  else where
                     sc%SC(:,:,:,nsc)=0.0_WP
                  end where
               end if
            end do
            ! Apply boundary conditions
            call sc%apply_bcond(time%t,time%dt)
         end block advance_scalar
         
		   ! Perform sub-iterations
		   do while (time%it.le.time%itmax)
				
				! Build mid-time velocity
			   fs%U=0.5_WP*(fs%U+fs%Uold)
				fs%V=0.5_WP*(fs%V+fs%Vold)
				fs%W=0.5_WP*(fs%W+fs%Wold)
				
				! Preliminary mass and momentum transport step at the interface
			   call fs%prepare_advection_upwind(dt=time%dt)
				
				! Explicit calculation of drho*u/dt from NS
			   call fs%get_dmomdt(resU,resV,resW)
				
				! Add momentum source terms
			   call fs%addsrc_gravity(resU,resV,resW)
				
				! Assemble explicit residual
			   resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
				resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
				resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
				
				! Form implicit residuals
			   call fs%solve_implicit(time%dt,resU,resV,resW)
				
				! Apply these residuals
			   fs%U=2.0_WP*fs%U-fs%Uold+resU
				fs%V=2.0_WP*fs%V-fs%Vold+resV
				fs%W=2.0_WP*fs%W-fs%Wold+resW
				
				! Apply other boundary conditions
			   call fs%apply_bcond(time%t,time%dt)
				
				! Solve Poisson equation
			   call fs%update_laplacian()
				call fs%correct_mfr()
				call fs%get_div()
				call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
				fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
				fs%psolv%sol=0.0_WP
				call fs%psolv%solve()
				call fs%shift_p(fs%psolv%sol)
				
				! Correct velocity
			   call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
				fs%P=fs%P+fs%psolv%sol
				fs%U=fs%U-time%dt*resU/fs%rho_U
				fs%V=fs%V-time%dt*resV/fs%rho_V
				fs%W=fs%W-time%dt*resW/fs%rho_W
				
				! Increment sub-iteration counter
			   time%it=time%it+1
				
			end do
			
			! Recompute interpolated velocity and divergence
		   call fs%interp_vel(Ui,Vi,Wi)
			call fs%get_div()
			
			! Output to ensight
		   if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
			
			! Perform and output monitoring
		   call fs%get_max()
			call vf%get_max()
         call sc%get_max(VF=vf%VF)
			call mfile%write()
			call cflfile%write()
         call scfile%write()
			
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
	   deallocate(resSC,resU,resV,resW,Ui,Vi,Wi)
	   
	end subroutine simulation_final
	
	
end module simulation