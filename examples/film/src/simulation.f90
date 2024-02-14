!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use tpns_class,        only: tpns
   use hypre_str_class,   only: hypre_str
   use vfs_class,         only: vfs
   !use ccl_class,         only: ccl
   use timetracker_class, only: timetracker
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Two-phase flow solver, linear solver for pressure, volume fraction solver, and a time tracker
   type(tpns)        :: fs
   type(hypre_str)   :: ps
   type(vfs)         :: vf
   type(timetracker) :: time

   !> CCL framework
   !type(ccl),         public :: cc
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Problem definition
   real(WP) :: rho_ratio,visc_ratio,Oh,dh
   integer :: nh
   real(WP), dimension(:,:), allocatable :: ph
   
   !> Post-processing info
   !integer :: ndrop
   !real(WP), dimension(:), allocatable :: drop_diam
   !real(WP) :: mean_diam,min_diam,max_diam
   !type(event) :: drop_evt
   !type(monitor) :: dropfile

contains
   
   
   !> Function that defines a level set function for a perforated film
   function levelset_perf_film(xyz,t) result(G)
      use mathtools, only: twoPi
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,dist
      real(WP), dimension(3) :: c,pos,myh
      integer :: n
      ! Create the film of unity thickness
      G=0.5_WP-abs(xyz(2))
      ! Add the holes
      do n=1,nh
         ! Store 2D positions
         pos=[xyz(1) ,0.0_WP,xyz(3) ]
         myh=[ph(1,n),0.0_WP,ph(2,n)]
         ! Account for periodicity
         if (myh(1)-pos(1).gt.+0.5_WP*cfg%xL) myh(1)=myh(1)-cfg%xL
         if (myh(1)-pos(1).lt.-0.5_WP*cfg%xL) myh(1)=myh(1)+cfg%xL
         if (myh(3)-pos(3).gt.+0.5_WP*cfg%zL) myh(3)=myh(3)-cfg%zL
         if (myh(3)-pos(3).lt.-0.5_WP*cfg%zL) myh(3)=myh(3)+cfg%zL
         ! Check hole is close enough
         dist=norm2(myh-pos)
         if (dist.lt.0.5_WP*dh) then
            ! Get closest point on torus
            c=myh+0.5_WP*dh*(pos-myh)/(dist+tiny(dist))
            G=0.5_WP-norm2(xyz-c)
         end if
      end do
   end function levelset_perf_film
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
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
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      

      ! Prepare random holes for film perforation
      initialize_holes: block
         use random,   only: random_uniform
         use parallel, only: MPI_REAL_WP
         use mpi_f08,  only: MPI_BCAST
         integer  :: n,nn,ierr
         real(WP) :: xh,zh,xxh,zzh
         logical  :: is_overlap
         real(WP), parameter :: safety_margin=1.5_WP
         ! Read in the hole size
         call param_read('Size of holes',dh)
         ! Read in the number of holes
         call param_read('Number of holes',nh)
         ! Allocate hole position array
         allocate(ph(2,nh))
         ! Root assigns positions to holes
         if (cfg%amRoot) then
            n=0
            do while (n.lt.nh)
               ! Draw a random position on the film
               xh=random_uniform(lo=cfg%x(cfg%imin),hi=cfg%x(cfg%imax+1))
               zh=random_uniform(lo=cfg%z(cfg%kmin),hi=cfg%z(cfg%kmax+1))
               ! Compare to all previous holes
               is_overlap=.false.
               do nn=1,n
                  ! Get the position of the other hole
                  xxh=ph(1,nn); zzh=ph(2,nn)
                  ! Account for periodicity
                  if (xxh-xh.gt.+0.5_WP*cfg%xL) xxh=xxh-cfg%xL
                  if (xxh-xh.lt.-0.5_WP*cfg%xL) xxh=xxh+cfg%xL
                  if (zzh-zh.gt.+0.5_WP*cfg%zL) zzh=zzh-cfg%zL
                  if (zzh-zh.lt.-0.5_WP*cfg%zL) zzh=zzh+cfg%zL
                  ! Check for overlap
                  if (norm2([xxh-xh,zzh-zh]).lt.safety_margin*dh) is_overlap=.true.
               end do
               ! If no overlap was found, add the hole to the list
               if (.not.is_overlap) then
                  n=n+1
                  ph(1,n)=xh
                  ph(2,n)=zh
               end if
            end do
         end if
         ! Broadcoast the hole positions
         call MPI_BCAST(ph,2*nh,MPI_REAL_WP,0,cfg%comm,ierr)
      end block initialize_holes

      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: swartz,lvira,elvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
         ! Initialize to a film
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_perf_film,0.0_WP,amr_ref_lvl)
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
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant density and viscosity to each phase
         fs%rho_l =1.0_WP; call param_read('Density ratio'  ,  rho_ratio); fs%rho_g =fs%rho_l / rho_ratio
         fs%visc_l=1.0_WP; call param_read('Viscosity ratio', visc_ratio); fs%visc_g=fs%visc_l/visc_ratio
         ! Read in Ohnesorge number and assign surface tension coefficient
         call param_read('Oh',Oh); fs%sigma=Oh**(-2)
         ! Prepare and configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=16
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='Film')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Prepare drop post-processing
      !postprocess_drop: block
      !   use vfs_class, only: VFlo
      !   ! Creat CCL object
      !   cc=ccl(cfg=cfg,name='CCL')
      !   cc%max_interface_planes=1
      !   cc%VFlo=VFlo
      !   ! Create event for drop size output
      !   drop_evt=event(time=time,name='Drop output')
      !   call param_read('Drop output period',drop_evt%tper)
      !   ! Prepare directory for drop size output
      !   if (vf%cfg%amRoot) call execute_command_line('mkdir -p drops')
      !   ! Perform first analysis
      !   call analyze_drops()
      !end block postprocess_drop
      
      
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
         ! Create drop monitor
         !dropfile=monitor(vf%cfg%amRoot,'drop')
         !call dropfile%add_column(time%n,'Timestep number')
         !call dropfile%add_column(time%t,'Time')
         !call dropfile%add_column(ndrop,'Ndrop')
         !call dropfile%add_column( min_diam, 'Min diameter')
         !call dropfile%add_column(mean_diam,'Mean diameter')
         !call dropfile%add_column( max_diam, 'Max diameter')
         !call dropfile%write()
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
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
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
         call fs%get_viscosity(vf=vf)
         
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
            !call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho_U
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho_V
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho_W
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
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
         call vf%update_surfmesh(smesh)
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()

         ! Count drops and get mean diameter
         !call analyze_drops()
         !call dropfile%write()

      end do
      
      
   end subroutine simulation_run
   
   
   !> Analyze VOF field for drops
   !subroutine analyze_drops()
   !   use mathtools, only: Pi
   !   use string,    only: str_medium
   !   implicit none
   !   integer :: n,iunit,ierr
   !   character(len=str_medium) :: filename,timestamp
   !   
   !   ! Perform CCL on VOF field
   !   call cc%build_lists(VF=vf%VF,U=fs%U,V=fs%V,W=fs%W)
   !   
   !   ! Store number of droplets and allocate diameter
   !   ndrop=cc%n_meta_struct
   !   if (allocated(drop_diam)) deallocate(drop_diam)
   !   allocate(drop_diam(ndrop))
   !   
   !   ! Loops over identified structures and get equivalent diameter
   !   mean_diam=0.0_WP
   !   do n=1,ndrop
   !      drop_diam(n)=(6.0_WP*cc%meta_structures_list(n)%vol/Pi)**(1.0_WP/3.0_WP)
   !      mean_diam=mean_diam+drop_diam(n)
   !   end do
   !   if (ndrop.gt.0) mean_diam=mean_diam/real(ndrop,WP)
   !   min_diam=minval(drop_diam)
   !   max_diam=maxval(drop_diam)
   !   
   !   ! Clean up CCL
   !   call cc%deallocate_lists()
   !   
   !   ! If root and if event triggers, print out drop sizes
   !   if (vf%cfg%amRoot.and.drop_evt%occurs()) then
   !      filename='diameter_'
   !      write(timestamp,'(es12.5)') time%t
   !      open(newunit=iunit,file='drops/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
   !      do n=1,ndrop
   !         write(iunit,'(es12.5)') drop_diam(n)
   !      end do
   !      close(iunit)
   !   end if
   !
   !end subroutine analyze_drops
   
   
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
