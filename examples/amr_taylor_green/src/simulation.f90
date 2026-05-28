!> AMR Taylor-Green Vortex test case for compressible solver
module simulation
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrcomp_class,     only: amrcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible solver
   type(timetracker) :: time
   type(amrcomp), target :: fs
   type(amrdata) :: dQdt,k1,k2,k3,k4
   
   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Rec_tag=huge(1.0_WP)
   real(WP) :: Res_tag=huge(1.0_WP)
   
   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile
   
   !> Flow parameters
   real(WP) :: Mach,Reynolds,Prandtl

   !> Stiffened gas EOS parameters
   real(WP) :: Gamma,Pinf,Cv
   
contains


   !> P=EOS(RHO,I) - Stiffened gas
   pure real(WP) function get_P(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_P=RHO*I*(Gamma-1.0_WP)-Gamma*Pinf
   end function get_P
   
   !> T=f(RHO,P)
   pure real(WP) function get_T(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_T=(P+Pinf)/(Cv*RHO*(Gamma-1.0_WP))
   end function get_T
   
   !> C=f(RHO,P) - Speed of sound
   pure real(WP) function get_C(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_C=sqrt(Gamma*(P+Pinf)/RHO)
   end function get_C

   !> Compute viscosity using Sutherland's law, zero bulk viscosity, and set diffusivity based on Prandtl number
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pT,pVisc,pBeta,pDiff
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pT=>fs%T%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Sutherland's law
               pVisc(i,j,k,1)=Reynolds**(-1.0_WP)*(1.4042_WP*pT(i,j,k,1)**1.5_WP)/(pT(i,j,k,1)+0.4042_WP)
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Heat diffusivity: k = Cp*mu/Pr = Cv*Gamma*mu/Pr
               pDiff(i,j,k,1)=Cv*Gamma*pVisc(i,j,k,1)/Prandtl
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities
   
   !> User init callback - set initial conditions for Q
   subroutine user_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,x,y,z,rho,P,U,V,W,IE
      integer :: i,j,k   
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         ! Get pointer to Q
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Get tilebox with overlap
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get physical coordinates
            x=(real(i,WP)+0.5_WP)*dx
            y=(real(j,WP)+0.5_WP)*dy
            z=(real(k,WP)+0.5_WP)*dz
            ! TGV velocity
            U=+sin(x)*cos(y)*cos(z)
            V=-cos(x)*sin(y)*cos(z)
            W=0.0_WP
            ! TGV pressure
            P=1.0_WP/(Gamma*Mach**2)+(cos(2.0_WP*x)+cos(2.0_WP*y))*(cos(2.0_WP*z)+2.0_WP)/16.0_WP
            ! Internal energy with T=1 normalization
            IE=Cv
            ! Density from EOS inversion
            rho=(P+Gamma*Pinf)/(IE*(Gamma-1.0_WP))
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=rho*U
            pQ(i,j,k,3)=rho*V
            pQ(i,j,k,4)=rho*W
            pQ(i,j,k,5)=rho*IE
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine user_init

   !> Tagger for this case based on velocity gradient (turbulence) and divergence (shocks)
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVisc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real(WP) :: vort_mag,div_neg,rho,mu,Rec,Res
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pVisc=>solver%visc%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Local density and viscosity
            rho=pQ(i,j,k,1)
            mu=pVisc(i,j,k,1)
            if (mu.le.0.0_WP) mu=1.0_WP/Reynolds ! visc may not have been set yet
            ! Centered velocity gradients
            dudx=0.5_WP*dxi*(pQ(i+1,j,k,2)/max(pQ(i+1,j,k,1),solver%rho_floor)-pQ(i-1,j,k,2)/max(pQ(i-1,j,k,1),solver%rho_floor))
            dudy=0.5_WP*dyi*(pQ(i,j+1,k,2)/max(pQ(i,j+1,k,1),solver%rho_floor)-pQ(i,j-1,k,2)/max(pQ(i,j-1,k,1),solver%rho_floor))
            dudz=0.5_WP*dzi*(pQ(i,j,k+1,2)/max(pQ(i,j,k+1,1),solver%rho_floor)-pQ(i,j,k-1,2)/max(pQ(i,j,k-1,1),solver%rho_floor))
            dvdx=0.5_WP*dxi*(pQ(i+1,j,k,3)/max(pQ(i+1,j,k,1),solver%rho_floor)-pQ(i-1,j,k,3)/max(pQ(i-1,j,k,1),solver%rho_floor))
            dvdy=0.5_WP*dyi*(pQ(i,j+1,k,3)/max(pQ(i,j+1,k,1),solver%rho_floor)-pQ(i,j-1,k,3)/max(pQ(i,j-1,k,1),solver%rho_floor))
            dvdz=0.5_WP*dzi*(pQ(i,j,k+1,3)/max(pQ(i,j,k+1,1),solver%rho_floor)-pQ(i,j,k-1,3)/max(pQ(i,j,k-1,1),solver%rho_floor))
            dwdx=0.5_WP*dxi*(pQ(i+1,j,k,4)/max(pQ(i+1,j,k,1),solver%rho_floor)-pQ(i-1,j,k,4)/max(pQ(i-1,j,k,1),solver%rho_floor))
            dwdy=0.5_WP*dyi*(pQ(i,j+1,k,4)/max(pQ(i,j+1,k,1),solver%rho_floor)-pQ(i,j-1,k,4)/max(pQ(i,j-1,k,1),solver%rho_floor))
            dwdz=0.5_WP*dzi*(pQ(i,j,k+1,4)/max(pQ(i,j,k+1,1),solver%rho_floor)-pQ(i,j,k-1,4)/max(pQ(i,j,k-1,1),solver%rho_floor))
            ! Vorticity magnitude
            vort_mag=sqrt((dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2)
            ! Divergence (only compression, i.e., negative divergence)
            div_neg=min(dudx+dvdy+dwdz,0.0_WP)
            ! Local cell Reynolds for turbulence
            Rec=rho*vort_mag*solver%amr%min_meshsize(lvl)**2/mu
            ! Local cell Reynolds for shocks
            Res=rho*abs(div_neg)*solver%amr%min_meshsize(lvl)**2/mu
            ! Tag based on either criterion
            if (Rec.gt.Rec_tag.or.Res.gt.Res_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read EoS and flow parameters
      init_eos_and_flow: block
         call param_read('Gamma',Gamma)
         call param_read('Pinf',Pinf)
         call param_read('Mach',Mach)
         call param_read('Reynolds',Reynolds)
         call param_read('Prandtl',Prandtl)
         Cv=1.0_WP/(Gamma*(Gamma-1.0_WP)*Mach**2)
      end block init_eos_and_flow
      
      ! Initialize AMR grid
      create_amrgrid: block
         use mathtools, only: twoPi
         amr%name='amrtgv'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=0.0_WP; amr%xhi=twoPi
         amr%ylo=0.0_WP; amr%yhi=twoPi
         amr%zlo=0.0_WP; amr%zhi=twoPi
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max levels',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible solver
      create_solver: block
         ! Create flow solver
         call fs%initialize(amr=amr)
         ! Provide thermodynamic model
         fs%getP=>get_P
         fs%getC=>get_C
         fs%getT=>get_T
         ! Set initial conditions
         fs%user_init=>user_init
      end block create_solver
      
      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=amrex_interp_none); call dQdt%register()
         call k1%initialize(amr,name='k1',ncomp=5,ng=0,interp=amrex_interp_none); call k1%register()
         call k2%initialize(amr,name='k2',ncomp=5,ng=0,interp=amrex_interp_none); call k2%register()
         call k3%initialize(amr,name='k3',ncomp=5,ng=0,interp=amrex_interp_none); call k3%register()
         call k4%initialize(amr,name='k4',ncomp=5,ng=0,interp=amrex_interp_none); call k4%register()
      end block create_workspace
      
      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Rec',Rec_tag)
         call param_read('Tagging Res',Res_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Compute viscosities
         call get_viscosities()
         ! Add artificial bulk viscosity
         call fs%add_viscartif(dt=time%dt)
      end block init_regridding
      
      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'TGV')
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%U,1,'U')
         call viz%add_scalar(fs%V,1,'V')
         call viz%add_scalar(fs%W,1,'W')
         call viz%add_scalar(fs%beta,1,'beta')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz
      
      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Qmin(1),'RHOmin')
         call mfile%add_column(fs%Qmax(1),'RHOmax')
         call mfile%add_column(fs%Imin,'Imin')
         call mfile%add_column(fs%Imax,'Imax')
         call mfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Mass')
         call consfile%add_column(fs%Qint(2),'U Momentum')
         call consfile%add_column(fs%Qint(3),'V Momentum')
         call consfile%add_column(fs%Qint(4),'W Momentum')
         call consfile%add_column(fs%Qint(5),'Internal energy')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitors

   end subroutine simulation_init
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old conserved variables
         call fs%Qold%copy(src=fs%Q)
         
         ! ===== RK4 Stage 1: k1 = f(t, Q) =====
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t); call k1%copy(src=dQdt)
         
         ! ===== RK4 Stage 2: k2 = f(t+dt/2, Q+k1/2) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=k1)  ! Q=Qold+dt/2*k1
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t+0.5_WP*time%dt); call k2%copy(src=dQdt)
         
         ! ===== RK4 Stage 3: k3 = f(t+dt/2, Q+k2/2) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=k2)  ! Q=Qold+dt/2*k2
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t+0.5_WP*time%dt); call k3%copy(src=dQdt)
         
         ! ===== RK4 Stage 4: k4 = f(t+dt, Q+k3) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=time%dt,src=k3)  ! Q=Qold+dt*k3
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+time%dt)
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t+time%dt); call k4%copy(src=dQdt)
         
         ! ===== RK4 Combination: Q = Qold + (k1 + 2*k2 + 2*k3 + k4)/6 =====
         call fs%Q%copy(src=fs%Qold)
         call fs%Q%saxpy(a=time%dt/6.0_WP,src=k1)
         call fs%Q%saxpy(a=time%dt/3.0_WP,src=k2)
         call fs%Q%saxpy(a=time%dt/3.0_WP,src=k3)
         call fs%Q%saxpy(a=time%dt/6.0_WP,src=k4)

         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         
         ! Recompute primitive variables
         call fs%get_primitive(fs%Q)
         
         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add artificial bulk viscosity
         call fs%add_viscartif(dt=time%dt)
         
         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call k1%finalize()
      call k2%finalize()
      call k3%finalize()
      call k4%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final
   
end module simulation
