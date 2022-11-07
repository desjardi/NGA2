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
   type(monitor) :: mfile,dfile,cflfile,consfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Choice of relaxation model
   integer  :: relax_model
   !> Problem definition
   real(WP) :: ddrop
   logical  :: yes_Temp
   
   !> Monitor output variables
   real(WP) :: xd,yd,zd, vmag2, vmaginf
   real(WP) :: LV0,mass0,mom0_x,mom0_y,mom0_z,totegy0
   real(WP) :: LVdiff,massdiff,momdiff_x,momdiff_y,momdiff_z,totegydiff

contains
  
   !> Updates quantities specific to the monitor of this case
   subroutine mast_drop_update()
     use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX
     use parallel,  only: MPI_REAL_WP
     implicit none
     real(WP) :: md, vol_tot, buf
     integer  :: ierr,i,j,k
     
     ! Initialize
     vol_tot = 0.0_WP
     md = 0.0_WP
     xd = 0.0_WP; yd = 0.0_WP; zd = 0.0_WP
     vmag2 = 0.0_WP; vmaginf = 0.0_WP
     do k=fs%cfg%kmin_,fs%cfg%kmax_
       do j=fs%cfg%jmin_,fs%cfg%jmax_
         do i=fs%cfg%imin_,fs%cfg%imax_
           ! Volume (for normalizing)
           vol_tot = vol_tot + fs%cfg%vol(i,j,k)
           ! Mass (for normalizing)
           md = md + fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*fs%Lrho(i,j,k)
           ! Position of liquid (multiplied by mass)
           xd = xd + vf%Lbary(1,i,j,k)*fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*fs%Lrho(i,j,k)
           yd = yd + vf%Lbary(2,i,j,k)*fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*fs%Lrho(i,j,k)
           zd = zd + vf%Lbary(3,i,j,k)*fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*fs%Lrho(i,j,k)
           ! Velocity norms
           vmag2 = vmag2 + fs%cfg%vol(i,j,k)*(fs%Ui(i,j,k)**2+fs%Vi(i,j,k)**2+fs%Wi(i,j,k)**2)
           vmaginf = max(vmaginf,sqrt(fs%Ui(i,j,k)**2+fs%Vi(i,j,k)**2+fs%Wi(i,j,k)**2))
         end do
       end do
     end do
     ! Normalizing quantities
     call MPI_ALLREDUCE(vol_tot,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); vol_tot = buf;
     call MPI_ALLREDUCE(md,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); md = buf;
     ! Droplet position
     call MPI_ALLREDUCE(xd,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); xd = buf/md;
     call MPI_ALLREDUCE(yd,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); yd = buf/md;
     call MPI_ALLREDUCE(zd,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); zd = buf/md;
     ! Velocity norms
     call MPI_ALLREDUCE(vmag2,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); vmag2 = sqrt(buf/vol_tot);
     call MPI_ALLREDUCE(vmaginf,buf,1,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr); vmaginf = buf;
     
   end subroutine mast_drop_update
   
   !> Updates quantities specific to the monitor of this case
   subroutine conservation_update()
     use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX
     use parallel,  only: MPI_REAL_WP
     implicit none
     real(WP) :: vol_tot, buf
     integer  :: ierr,i,j,k
     
     ! Initialize
     vol_tot = 0.0_WP
     LVdiff = 0.0_WP; massdiff = 0.0_WP; totegydiff = 0.0_WP
     momdiff_x = 0.0_WP; momdiff_y = 0.0_WP; momdiff_z = 0.0_WP
     do k=fs%cfg%kmin_,fs%cfg%kmax_
       do j=fs%cfg%jmin_,fs%cfg%jmax_
         do i=fs%cfg%imin_,fs%cfg%imax_
           ! Volume (for normalizing)
           vol_tot = vol_tot + fs%cfg%vol(i,j,k)
           ! Liquid volume
           LVdiff = LVdiff + fs%cfg%vol(i,j,k)*vf%VF(i,j,k)
           ! Mass
           massdiff = massdiff + fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)
           ! Momentum
           momdiff_x = momdiff_x + fs%cfg%vol(i,j,k)*fs%rhoUi(i,j,k)
           momdiff_y = momdiff_y + fs%cfg%vol(i,j,k)*fs%rhoVi(i,j,k)
           momdiff_z = momdiff_z + fs%cfg%vol(i,j,k)*fs%rhoWi(i,j,k)
           ! Total energy
           totegydiff = totegydiff + fs%cfg%vol(i,j,k)*(vf%VF(i,j,k)*fs%LrhoE(i,j,k)+(1.0-vf%VF(i,j,k))*fs%GrhoE(i,j,k))
         end do
       end do
     end do
     
     ! Parallel sums
     call MPI_ALLREDUCE(vol_tot,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); vol_tot = buf;
     ! Normalizing quantities is not necessary, but makes things more understandable
     call MPI_ALLREDUCE(LVdiff,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); LVdiff = buf/vol_tot;
     call MPI_ALLREDUCE(massdiff,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); massdiff = buf/vol_tot;
     call MPI_ALLREDUCE(momdiff_x,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); momdiff_x = buf/vol_tot;
     call MPI_ALLREDUCE(momdiff_y,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); momdiff_y = buf/vol_tot;
     call MPI_ALLREDUCE(momdiff_z,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); momdiff_z = buf/vol_tot;
     call MPI_ALLREDUCE(totegydiff,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); totegydiff = buf/vol_tot;
     ! Subtract for differences
     LVdiff = LVdiff - LV0
     massdiff = massdiff - mass0
     momdiff_x = momdiff_x - mom0_x; momdiff_y = momdiff_y - mom0_y; momdiff_z = momdiff_z - mom0_z
     totegydiff = totegydiff - totegy0
     
   end subroutine conservation_update
   
   !> Function that defines a level set function for a spherical droplet at center
   function levelset_drop_center(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)/(ddrop/2.0)
   end function levelset_drop_center
   
   function levelset_cyl_center(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP-sqrt((xyz(1)/ddrop*2.0_WP)**2+(xyz(2)/ddrop*2.0_WP)**2)
   end function levelset_cyl_center

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
         call param_read('Max steps',time%nmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker


      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,elvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
         ! Initialize liquid at left
         call param_read('Droplet diameter',ddrop)
         ! Droplet always initialized at center
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
                  if (vf%cfg%nz.eq.1) then
                    call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cyl_center,0.0_WP,amr_ref_lvl)
                  else
                    call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_drop_center,0.0_WP,amr_ref_lvl)
                  end if
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                    vf%Lbary(:,i,j,k)=v_cent
                    vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                    if (vf%cfg%nz.eq.1) vf%Gbary(3,i,j,k)=v_cent(3);
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
         use mast_class, only: clipped_neumann,dirichlet,bc_scope,bcond,mech_egy_mech_hhz,thermmech_egy_mech_hhz
         use ils_class,  only: pcg_bbox,pcg_amg
         use mathtools,  only: Pi
         use parallel,   only: amRoot
         use string,     only: str_medium
         integer :: i,j,k,impl_option
         real(WP) :: gamm_l,Pref_l,gamm_g,visc_l,visc_g,hdff_l,hdff_g,cv_l,cv_g,b_l,q_l
         real(WP) :: GSS, GP0, LP0, Grho0, Lrho0, GTemp0, LTemp0
         real(WP), dimension(3) :: u_g,u_l,u_mix
         character(len=str_medium) :: impl_str
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get EOS parameters from input
         call param_read('Liquid gamma',gamm_l)
         call param_read('Liquid Pref', Pref_l)
         call param_read('Gas gamma',gamm_g)
         ! Check if NASG EOS is being used
         b_l = 0.0_WP
         if (param_exists('Liquid b')) then
           call param_read('Liquid b', b_l)
           call param_read('Liquid q', q_l)
         end if
         ! Register equations of state
         if (b_l .eq. 0.0_WP) then
           call matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
         else
           call matmod%register_NobleAbelstiffenedgas('liquid',gamm_l,Pref_l,q_l,b_l)
         end if
         call matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',visc_l)
         call param_read('Gas dynamic viscosity',visc_g)
         ! Set heat diffusion to 0 unless specified
         hdff_l = 0.0_WP; hdff_g = 0.0_WP
         if (param_exists('Liquid therm cond')) then
           call param_read('Liquid therm cond', hdff_l)
           call param_read('Gas therm cond', hdff_g)
         end if
         ! Set specific heat to default unless specified
         cv_l = matmod%cv_l0; cv_g = matmod%cv_g0
         if (param_exists('Liquid spec heat')) then
           call param_read('Liquid spec heat', cv_l)
           call param_read('Gas spec heat', cv_g)
         end if
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g,&
            viscconst_liquid=visc_l, hdffconst_gas=hdff_g, hdffconst_liquid=hdff_l, &
            sphtconst_gas=cv_g, sphtconst_liquid=cv_l)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver, default is amg
         impl_option = pcg_amg
         if (param_exists('Implicit solver')) call param_read('Implicit solver',impl_str)
         if (trim(impl_str).eq.'blackbox') impl_option = pcg_bbox
         call fs%setup(pressure_ils=pcg_bbox,implicit_ils=impl_option)

         ! Set up problem: velocity, density, pressure
         ! Velocity is quiescent unless specified
         u_g = 0.0_WP; u_l = 0.0_WP
         yes_Temp = .false.
         if (param_exists('Liquid velocity')) then
           call param_read('Liquid velocity',u_l)
         end if
         if (param_exists('Gas velocity')) then
           call param_read('Gas velocity',u_g)
         end if
         if (param_exists('Background pressure')) then
           ! Set background pressure if specified
           call param_read('Background pressure',GP0)
           ! Calculate liquid pressure
           if (fs%cfg%nz.eq.1) then
             ! Cylinder configuration, curv = 1/r
             LP0 = GP0 + 2.0/ddrop*fs%sigma
           else
             ! Sphere configuration, curv = 1/r + 1/r
             LP0 = GP0 + 4.0/ddrop*fs%sigma
           end if
           if (param_exists('Liquid density')) then
             ! Set phase density if specified
             call param_read('Liquid density',Lrho0)
             call param_read('Gas density',Grho0)
           else
             ! Calculate phase density from specified temperature
             yes_Temp = .true.
             call param_read('Liquid temperature',LTemp0)
             call param_read('Gas temperature',GTemp0)
             Lrho0 = matmod%EOS_density(LP0,LTemp0,cv_l,'liquid')
             Grho0 = matmod%EOS_density(GP0,GTemp0,cv_g,'gas')
           end if
         else
           ! Get specified sound speed in gas
           call param_read('Gas sound speed',GSS)
           ! Get phase density
           call param_read('Liquid density',Lrho0)
           call param_read('Gas density',Grho0)
           ! Calculate background pressure
           GP0 = Grho0*GSS**2/gamm_g
           ! Calculate liquid pressure
           if (fs%cfg%nz.eq.1) then
             ! Cylinder configuration, curv = 1/r
             LP0 = GP0 + 2.0/ddrop*fs%sigma
           else
             ! Sphere configuration, curv = 1/r + 1/r
             LP0 = GP0 + 4.0/ddrop*fs%sigma
           end if
         end if
         ! Get gravity if specified
         if (param_exists('Gravity')) then
            call param_read('Gravity',fs%gravity)
         end if
         if (amRoot) then
           print*,"===== Problem Setup Non-dimensional Numbers ====="
           if (visc_g.ne.0.0_WP) then
             print*,'Reynolds (gas)    = ', Grho0*sqrt(u_g(1)**2+u_g(2)**2+u_g(3)**2)*ddrop/visc_g
             print*,'Reynolds (liquid) = ', Lrho0*sqrt(u_l(1)**2+u_l(2)**2+u_l(3)**2)*ddrop/visc_l
             print*,'Acoustic Reynolds = ', Grho0*sqrt(gamm_g*GP0/Grho0)*ddrop/visc_g
             print*,'Laplace           = ', fs%sigma*Lrho0*ddrop/visc_l**2
           else
             print*,'Reynolds (gas)    = Infinity'
             print*,'Reynolds (liquid) = Infinity'
             print*,'Acoustic Reynolds = Infinity'
             if (fs%sigma.ne.0.0_WP) then
               print*,'Laplace           = Infinity'
             else
               print*,'Laplace           = Undefined'
             end if
           end if
           if (hdff_g.ne.0.0_WP) then
             print*,'Prandtl (gas)     = ', visc_g*gamm_g*cv_g/hdff_g
             print*,'Prandtl (liquid)  = ', visc_l*gamm_l*cv_l/hdff_l
           else
             if (visc_g.ne.0.0_WP) then
               print*,'Prandtl (gas)     = Infinity'
               print*,'Prandtl (liquid)  = Infinity'
             else
               print*,'Prandtl (gas)     = Undefined'
               print*,'Prandtl (liquid)  = Undefined'
             end if
           end if
           if (fs%sigma.ne.0.0_WP) then
             print*,'Weber             = ', Grho0*((u_g(1)-u_l(1))**2+(u_g(2)-u_l(2))**2+(u_g(3)-u_l(3))**2)*ddrop/fs%sigma
             print*,'Acoustic Weber    = ', gamm_g*GP0*ddrop/fs%sigma
           else
             print*,'Weber             = Infinity'
             print*,'Acoustic Weber    = Infinity'
           end if
           if (sum(abs(fs%gravity)).ne.0.0_WP) then
             print*,'Froude            = ', sqrt(u_l(1)**2+u_l(2)**2+u_l(3)**2)/(fs%gravity(1)**2+fs%gravity(2)**2+fs%gravity(3)**2)/sqrt(ddrop)
           else if (sqrt(u_l(1)**2+u_l(2)**2+u_l(3)**2).ne.0.0_WP) then
             print*,'Froude            = Infinity'
           else
             print*,'Froude            = Undefined'
           end if
         end if

         ! Initialize flow variables
         fs%Lrho = Lrho0
         fs%Grho = Grho0
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
           do j=vf%cfg%jmino_,vf%cfg%jmaxo_
             do i=vf%cfg%imino_,vf%cfg%imaxo_
               ! Get local velocity based on mixture momentum, density
               u_mix = ((vf%VF(i,j,k))*Lrho0*u_l+(1.0-vf%VF(i,j,k))*Grho0*u_g) / &
                       ((vf%VF(i,j,k))*Lrho0+(1.0-vf%VF(i,j,k))*Grho0)
               ! Set velocity and energy
               fs%Ui(i,j,k) = u_mix(1)
               fs%Vi(i,j,k) = u_mix(2)
               fs%Wi(i,j,k) = u_mix(3)
               fs%GrhoE(i,j,k) = matmod%EOS_energy(GP0,Grho0,u_mix(1),u_mix(2),u_mix(3),'gas')
               fs%LrhoE(i,j,k) = matmod%EOS_energy(LP0,Lrho0,u_mix(1),u_mix(2),u_mix(3),'liquid')
             end do
           end do
         end do

         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         ! Need to fill ghost cells
         bc_scope='velocity'
         call fs%apply_bcond(time%dt,bc_scope)
         
         ! Choose relaxation model
         if (abs(cv_l)+abs(cv_g).gt.0.0_WP) then
           ! Use thermo-mechanical if heat conduction is considered
           relax_model = thermmech_egy_mech_hhz
         else
           ! Use mechanical otherwise
           relax_model = mech_egy_mech_hhz
         end if

         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         ! Perform initial pressure relax
         call fs%pressure_relax(vf,matmod,relax_model)
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA
         ! Calculate initial temperature
         call matmod%update_temperature(vf,fs%Tmptr)
         ! Initialize first guess for pressure (0 works best)
         fs%psolv%sol=0.0_WP

      end block create_and_initialize_flow_solver


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='MASTdrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('Peos',fs%PA)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_scalar('Temperature',fs%Tmptr)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
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
         call mfile%add_column(fs%Tmax,'Tmax')
         call mfile%write()
         ! Create drop monitor
         dfile=monitor(fs%cfg%amRoot,'drop')
         call dfile%add_column(time%n,'Timestep number')
         call dfile%add_column(time%t,'Time')
         call dfile%add_column(xd,'xdrop')
         call dfile%add_column(yd,'ydrop')
         call dfile%add_column(zd,'zdrop')
         call dfile%add_column(vmag2,'Vmag2Norm')
         call dfile%add_column(vmaginf,'VmagInfNorm')
         call mast_drop_update()
         call dfile%write()
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
         ! Create convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'cvg')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(fs%impl_it_x,'Impl_x iteration')
         call cvgfile%add_column(fs%impl_rerr_x,'Impl_x error')
         call cvgfile%add_column(fs%impl_it_y,'Impl_y iteration')
         call cvgfile%add_column(fs%impl_rerr_y,'Impl_y error')
         call cvgfile%add_column(fs%implicit%it,'Impl_z iteration')
         call cvgfile%add_column(fs%implicit%rerr,'Impl_z error')
         call cvgfile%add_column(fs%psolv%it,'Pressure iteration')
         call cvgfile%add_column(fs%psolv%rerr,'Pressure error')
         
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(LVdiff,'D_LiqVol')
         call consfile%add_column(massdiff,'D_Mass')
         call consfile%add_column(momdiff_x,'D_XMomentum')
         call consfile%add_column(momdiff_y,'D_YMomentum')
         call consfile%add_column(momdiff_z,'D_ZMomentum')
         call consfile%add_column(totegydiff,'D_TotEnergy')
         ! Set initial quantities to zero
         LV0 = 0.0_WP; mass0 = 0.0_WP; totegy0 = 0.0_WP
         mom0_x = 0.0_WP; mom0_y = 0.0_WP; mom0_z = 0.0_WP
         call conservation_update()
         ! Copy current quantities to initial
         LV0 = LVdiff; mass0 = massdiff; totegy0 = totegydiff
         mom0_x = momdiff_x; mom0_y = momdiff_y; mom0_z = momdiff_z
         ! Reset differences to zero, by definition
         LVdiff = 0.0_WP; massdiff = 0.0_WP; totegydiff = 0.0_WP
         momdiff_x = 0.0_WP; momdiff_y = 0.0_WP; momdiff_z = 0.0_WP
         call consfile%write()
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

            ! Diffusion and source term step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)
            
            ! Record convergence monitor
            call cvgfile%write()

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         fs%PA = matmod%EOS_all(vf);
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call mast_drop_update()
         call dfile%write()
         call cflfile%write()
         call conservation_update()
         call consfile%write()

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      use string,    only: str_long
      implicit none
      integer,          parameter :: col_len=14
      character(len=*), parameter :: aformat='(a12)'
      character(len=*), parameter :: rformat='(es12.5)'
      integer  :: i,j,k,n,ierr,iunit
      real(WP) :: dxr, r, buf0, buf1, l_left, l_right
      real(WP), dimension(:), allocatable :: Tprof, RHOprof, normprof, xmprof
      character(len=str_long) :: line, header
      
      ! Output temperature and density profiles 
      if (yes_Temp) then 
        
        ! Resolution of final profile
        dxr = (2.0_WP/real(fs%cfg%nx,WP))
        ! Allocate and zero profile arrays
        allocate(Tprof(fs%cfg%nx/2),RHOprof(fs%cfg%nx/2),normprof(fs%cfg%nx/2),xmprof(fs%cfg%nx/2))
        Tprof = 0.0_WP; RHOprof = 0.0_WP; normprof = 0.0_WP
        ! Populate the x profile
        do i=1,fs%cfg%nx/2
          xmprof(i) = 0.5_WP*dxr + real(i-1,WP)*dxr
        end do
        
        do k=fs%cfg%kmin_,fs%cfg%kmax_
          do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
              ! Calculate current radius
              r = sqrt(fs%cfg%xm(i)**2+fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
              
              !! == For plotting profile == !!
              ! Determine which increment of the domain this belongs to
              n = 1+floor((r-0.5_WP*dxr)/dxr)
              ! Contribute to the left
              if (n.le.fs%cfg%nx/2) then
                l_left = real(n,WP)*dxr-(r-0.5_WP*dxr)
                Tprof(n) = Tprof(n) + l_left*fs%Tmptr(i,j,k)
                RHOprof(n) = RHOprof(n) + l_left*fs%RHO(i,j,k)
                normprof(n) = normprof(n) + l_left
              end if
              ! Contribute to the right
              if (n+1.le.fs%cfg%nx/2) then
                l_right = (r+0.5_WP*dxr)-real(n,WP)*dxr
                Tprof(n+1) = Tprof(n+1) + l_right*fs%Tmptr(i,j,k)
                RHOprof(n+1) = RHOprof(n+1) + l_right*fs%RHO(i,j,k)
                normprof(n+1) = normprof(n+1) + l_right
              end if
              
            end do
          end do
        end do
        ! Create final profile array using info from every proc
        do n=1,fs%cfg%nx/2
          ! Sum normalization at current location
          call MPI_ALLREDUCE(normprof(n),buf0,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr);
          ! Sum temperature
          call MPI_ALLREDUCE(Tprof(n),buf1,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr);
          ! Update temperature
          Tprof(n) = buf1/buf0
          ! Sum density
          call MPI_ALLREDUCE(RHOprof(n),buf1,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr);
          ! Update density
          RHOprof(n) = buf1/buf0
        end do
        ! Write profile to file
        if (fs%cfg%amRoot) then
          ! Open file
          if (fs%sigma.eq.0.0_WP) then
            open(newunit=iunit,file='profile_noST.txt',form="formatted",iostat=ierr,status="REPLACE")
          else
            open(newunit=iunit,file='profile_ST.txt',form="formatted",iostat=ierr,status="REPLACE")
          end if
          ! Write header
          write (header(1+0*col_len:),aformat) 'xm'
          write (header(1+1*col_len:),aformat) 'T'
          write (header(1+2*col_len:),aformat) 'RHO'
          write(iunit,'(a)') trim(header)
          ! Write data
          do n=1,fs%cfg%nx/2
            ! Create the line to dump
            write (line(1+0*col_len:),rformat) xmprof(n)
            write (line(1+1*col_len:),rformat) Tprof(n)
            write (line(1+2*col_len:),rformat) RHOprof(n)
            ! Dump the line
            write(iunit,'(a)') trim(line)
          end do
          ! Close file
          close(iunit)
        end if
      end if

      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays - none

   end subroutine simulation_final

end module simulation
