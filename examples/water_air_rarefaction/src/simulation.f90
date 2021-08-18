!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Problem definition
   real(WP), dimension(3) :: center1,center2,vel1,vel2
   real(WP) :: radius1,radius2
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
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
         use vfs_class, only: r2p,lvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize two droplets
         call param_read('Droplet 1 diameter',radius1); radius1=0.5_WP*radius1
         call param_read('Droplet 1 position',center1)
         call param_read('Droplet 2 diameter',radius2); radius2=0.5_WP*radius2
         call param_read('Droplet 2 position',center2)
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
                  !call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_colliding_drops,0.0_WP,amr_ref_lvl)
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
         use ils_class, only: gmres_amg
         use mathtools, only: Pi
         integer :: i,j,k
         real(WP), dimension(3) :: xyz
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l0)
         call param_read('Gas dynamic viscosity',fs%visc_g0)
         ! Assign constant density to each phase
         !call param_read('Liquid density',fs%rho_l)
         !call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=gmres_amg,implicit_ils=gmres_amg)
         ! Initial droplet velocity
         call param_read('Droplet 1 velocity',vel1)
         call param_read('Droplet 2 velocity',vel2)
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  ! U velocity
                  xyz=[fs%cfg%x(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  if (radius1-sqrt(sum((xyz-center1)**2)).ge.0.0_WP) fs%U(i,j,k)=vel1(1)
                  if (radius2-sqrt(sum((xyz-center2)**2)).ge.0.0_WP) fs%U(i,j,k)=vel2(1)
                  ! V velocity
                  xyz=[fs%cfg%xm(i),fs%cfg%y(j),fs%cfg%zm(k)]
                  if (radius1-sqrt(sum((xyz-center1)**2)).ge.0.0_WP) fs%V(i,j,k)=vel1(2)
                  if (radius2-sqrt(sum((xyz-center2)**2)).ge.0.0_WP) fs%V(i,j,k)=vel2(2)
                  ! V velocity
                  xyz=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%z(k)]
                  if (radius1-sqrt(sum((xyz-center1)**2)).ge.0.0_WP) fs%W(i,j,k)=vel1(3)
                  if (radius2-sqrt(sum((xyz-center2)**2)).ge.0.0_WP) fs%W(i,j,k)=vel2(3)
               end do
            end do
         end do
         ! Initialize flow properties and therm variables
         !call fs%therm_init()
         ! Calculate face velocities
         !call fs%interp_vel()
         ! Perform initial pressure relax
         !call fs%pressure_relax()

      end block create_and_initialize_flow_solver
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='CollidingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',vf%surfgrid)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
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
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         !call mfile%add_column(fs%divmax,'Maximum divergence')
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
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old flow variables (mixture)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoUold=fs%GrhoU; fs%LrhoUold=fs%LrhoU
         fs%GrhoVold=fs%GrhoV; fs%LrhoVold=fs%LrhoV
         fs%GrhoWold=fs%GrhoW; fs%LrhoWold=fs%LrhoW
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()
         
         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Apply boundary conditions

         ! Other routines to add later: sgs, lpt, prescribe

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(vf)

         ! Initialize loop quantities
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf)

            ! Insert viscous step here, or possibly incorporate into predictor above

            ! Insert sponge step here
            
            ! Helmholtz solve
            ! call fs%calcH_LHS(fs%psolv%op)
            ! call fs%interp_vel(fs%U,fs%V,fs%W)
            ! call fs%calcH_RHS(fs%psolv%rhs)
            ! call fs%psolv%setup()
            ! fs%psolv%sol=0.0_WP
            ! call fs%psolv%solve()
            
            ! ! Corrector step
            ! fs%P=fs%P+fs%psolv%sol
            ! call fs%correct_face(fs%psolv%sol)
            ! call fs%correct_cellcenter()
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do

         ! Pressure relaxation
         !call fs%pressure_relax()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
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
      
      ! Deallocate work arrays
      
   end subroutine simulation_final
   
end module simulation
