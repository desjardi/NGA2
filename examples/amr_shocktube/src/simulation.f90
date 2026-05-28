!> 1D multiphase shock tube simulation with amr code
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use amrgrid_class,     only: amrgrid
   use amrmpcomp_class,   only: amrmpcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> Stiffened gas EOS parameters (liquid and gas)
   real(WP) :: GammaL,PinfL,CvL
   real(WP) :: GammaG,PinfG,CvG

   !> Flow parameters
   real(WP) :: rhoG,pG_init            !< Gas state (right side)
   real(WP) :: rhoL,pL_init            !< Liquid state (left side)
   real(WP) :: x_int                   !< Interface position

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for planar interface at x = x_int
   !> Returns positive for liquid (x < x_int), negative for gas (x >= x_int)
   function planar_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=x_int-xyz(1)
   end function planar_levelset

   !> Liquid EOS: P=f(RHO,I) - Stiffened gas
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> Liquid EOS: T=f(RHO,P)
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> Liquid EOS: C=f(RHO,P)
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(max(0.0_WP,GammaL*(P+PinfL)/RHO))
   end function get_CL
   !> Liquid EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IL=(P+GammaL*PinfL)/(RHO*(GammaL-1.0_WP))
   end function get_IL

   !> Gas EOS: P=f(RHO,I) - Ideal gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> Gas EOS: T=f(RHO,P)
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> Gas EOS: C=f(RHO,P)
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(max(0.0_WP,GammaG*(P+PinfG)/RHO))
   end function get_CG
   !> Gas EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IG=(P+GammaG*PinfG)/(RHO*(GammaG-1.0_WP))
   end function get_IG

   !> Generalized mechanical relaxation for stiffened gas EOS pair
   !> Solves quadratic for equilibrium pressure Peq where PL+Pjump=PG=Peq,
   !> then adjusts VF and internal energies via p*dV work exchange.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7)
   !> Enforces pressure jump provided in Pjump
   subroutine P_relax_generalized(VF,Q,Pjump)
      use amrmpcomp_class, only: VFlo,VFhi
      implicit none
      real(WP), intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP), intent(in) :: Pjump
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      real(WP), parameter :: RHOGmin=1.0e-2_WP
      real(WP), parameter :: phist=1.0_WP,phi0=0.0_WP   !< Temporal weighting, phist=1 should yield best results
      ! Skip relaxation for the first few timesteps to allow pressure to stabilize (since IC are not in mechanical equilibrium)
      if (time%t.lt.5.0e-6_WP) return
      ! Skip if any conserved quantity is non-positive (EOS undefined)
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2))
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)**2
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG)**2
      cJ=ZL/(ZG+ZL)
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      n1=VF*phist
      n0=VF*(phi0*Pint-phist*cJ*pjump)+Q(3)
      d1=phist+1.0_WP/(GammaL-1.0_WP)
      d0=phi0*Pint-phist*cJ*pjump+GammaL/(GammaL-1.0_WP)*PinfL
      a=d1*(1.0_WP/(GammaG-1.0_WP)+phist*VF)+n1*(-1.0_WP/(GammaG-1.0_WP)-phist)
      b=d1*((GammaG*PinfG-pjump)/(GammaG-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump))+n1*(-(GammaG*PinfG-pjump)/(GammaG-1.0_WP)-phi0*Pint+phist*cJ*pjump)+d0*(1.0_WP/(GammaG-1.0_WP)+phist*VF)+n0*(-1.0_WP/(GammaG-1.0_WP)-phist)
      d=d0*((GammaG*PinfG-pjump)/(GammaG-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump))+n0*(-(GammaG*PinfG-pjump)/(GammaG-1.0_WP)-phi0*Pint+phist*cJ*pjump)
      ! Get equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.-PinfL.or.Peq-Pjump.le.-PinfG) return
      ! Get equilibrium volume fraction
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.VFlo.or.VFeq.gt.VFhi) return
      ! Adjust conserved quantities
      Q(3)=Q(3)-(phi0*Pint+phist*Peq)*(VFeq-VF)
      Q(4)=Q(4)+(phi0*Pint+phist*Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_generalized

   !> Set viscosities to zero
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=0.0_WP)
      call fs%beta%setval(val=0.0_WP)
      call fs%diff%setval(val=0.0_WP)
   end subroutine get_viscosities

   !> User init callback - set Q and VF for 1D shocktube
   !> Liquid on left (x < x_int), gas on right (x >= x_int)
   subroutine shocktube_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,IEG
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energies
      IEL=get_IL(rhoL,pL_init)
      IEG=get_IG(rhoG,pG_init)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF from planar levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=planar_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=myVF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Set conserved variables: Q=(VF*rhoL,(1-VF)*rhoG,VF*rhoL*IL,(1-VF)*rhoG*IG,0,0,0)
            pQ(i,j,k,1)=(       myVF)*rhoL
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*IEG
            pQ(i,j,k,5)=0.0_WP
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shocktube_init

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         ! Gas EoS parameters (ideal gas: Pinf=0)
         call param_read('GammaG',GammaG)
         PinfG=0.0_WP
         ! Liquid EoS parameters (stiffened gas)
         call param_read('GammaL',GammaL)
         call param_read('PinfL',PinfL)
         ! Specific heats
         call param_read('CvG',CvG)
         call param_read('CvL',CvL)
         ! Gas state (right side)
         call param_read('Gas density',rhoG)
         call param_read('Gas pressure',pG_init)
         ! Liquid state (left side)
         call param_read('Liquid density',rhoL)
         call param_read('Liquid pressure',pL_init)
         ! Interface location
         call param_read('Interface location',x_int)
         ! Log
         write(message,'("[Gas]    GammaG=",es12.5," PinfG=",es12.5," CvG=",es12.5)') GammaG,PinfG,CvG; call log(message)
         write(message,'("[Liquid] GammaL=",es12.5," PinfL=",es12.5," CvL=",es12.5)') GammaL,PinfL,CvL; call log(message)
         write(message,'("[Gas]    rhoG=",es12.5," pG=",es12.5)') rhoG,pG_init; call log(message)
         write(message,'("[Liquid] rhoL=",es12.5," pL=",es12.5)') rhoL,pL_init; call log(message)
         write(message,'("[Interface] x_int=",es12.5)') x_int; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid (1D domain: [-1.5, 1.5])
      create_amrgrid: block
         amr%name='shocktube'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-1.5_WP; amr%xhi=+1.5_WP
         amr%ylo=-1.0_WP; amr%yhi=+1.0_WP
         amr%zlo=-1.0_WP; amr%zhi=+1.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-1D: set y extent to match one cell at finest level
         if (amr%ny.eq.1) then
            amr%ylo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%yhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
         ! Enable quasi-2D/1D: set z extent to match one cell at finest level
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
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

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS,BC_LIQ
         use amrdata_class,    only: interp_face_lin
         ! Create flow solver
         call fs%initialize(amr=amr,name='shocktube')
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         if (amr%ny.eq.1) fs%interp_vel=interp_face_lin
         ! Provide thermodynamic model (6 EOS pointers)
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getTG=>get_TG
         ! Provide pressure relaxation model
         fs%relax=>P_relax_generalized
         ! Set initial conditions
         fs%user_init=>shocktube_init
         ! Set BCs: outflow (foextrap) on both x boundaries
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_LIQ
            fs%hi_bc(1)=BC_GAS
            fs%Q%lo_bc(1,:)=amrex_bc_foextrap
            fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%U%lo_bc(1,:)=amrex_bc_foextrap; fs%U%hi_bc(1,:)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)=amrex_bc_foextrap; fs%V%hi_bc(1,:)=amrex_bc_foextrap
            fs%W%lo_bc(1,:)=amrex_bc_foextrap; fs%W%hi_bc(1,:)=amrex_bc_foextrap
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=7,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize grid (no regridding for uniform mesh)
      init_grid: block
         ! Fresh start
         call amr%init_from_scratch(time=time%t)
         ! Build PLIC
         call fs%build_plic(time%t)
         call fs%build_subVF()
         ! Initialize primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! Initialize face velocities
         call fs%get_face_velocity()
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Set viscosities to zero
         call get_viscosities()
         ! Add artificial viscosity only
         call fs%add_viscartif(dt=time%dt,Cartif=5.0_WP)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_grid

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'shocktube',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%C,1,'C')
         call viz%add_surfmesh(fs%smesh,'plic')
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
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%write()
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
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
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
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_sl,'sl_max')
         call tfile%add_column(fs%wtmin_sl,'sl_min')
         call tfile%add_column(fs%wtmax_fv,'fv_max')
         call tfile%add_column(fs%wtmin_fv,'fv_min')
         call tfile%add_column(fs%wtmax_div,'div_max')
         call tfile%add_column(fs%wtmin_div,'div_min')
         call tfile%add_column(fs%wtmax_plicnet,'plicnet_max')
         call tfile%add_column(fs%wtmin_plicnet,'plicnet_min')
         call tfile%add_column(fs%wtmax_polygon,'polygon_max')
         call tfile%add_column(fs%wtmin_polygon,'polygon_min')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%add_column(fs%nmixed_min,'mixed_min')
         call tfile%write()
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

         ! Remember old state
         call fs%store_old()

         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Get most up-to-date pressure
         call fs%apply_relax(time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities
         call fs%get_face_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         ! Add surface tension term
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Get most up-to-date pressure
         call fs%apply_relax(time=time%t)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities
         call fs%get_face_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=time%dt)
         ! Add surface tension term
         call fs%add_surface_tension(scale=time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================================================================================

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cartif=5.0_WP)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

      ! Extract mixture data at final time for x in [0, 1.5]
      mixture_data: block
         use mpi_f08
         use parallel,         only: MPI_REAL_WP
         use string,           only: str_medium
         use filesys,          only: makedir,isdir
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         integer :: i,j,k,lvl,iunit,ierr,nlocal,ntotal
         integer, dimension(:), allocatable :: recvcount,displs
         real(WP), dimension(:), allocatable :: local_x,local_VF,local_RHO,local_P,local_UVW
         real(WP), dimension(:), allocatable :: global_x,global_VF,global_RHO,global_P,global_UVW
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pPL,pPG,pUVW
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: x_cc,vf_val,dx
         character(len=str_medium) :: filename,timestamp
         lvl=0 ! Single level, no AMR
         dx=amr%dx(lvl)
         ! First pass: count local cells in range [0, 1.5]
         nlocal=0
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            do i=bx%lo(1),bx%hi(1)
               x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x_cc.ge.0.0_WP.and.x_cc.le.1.5_WP) nlocal=nlocal+1
            end do
         end do
         call amr%mfiter_destroy(mfi)
         ! Allocate local arrays
         allocate(local_x(nlocal),local_VF(nlocal),local_RHO(nlocal),local_P(nlocal),local_UVW(nlocal))
         ! Second pass: fill local arrays
         nlocal=0
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
            pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
            pPL  =>fs%PL%mf(lvl)%dataptr(mfi)
            pPG  =>fs%PG%mf(lvl)%dataptr(mfi)
            pUVW   =>fs%UVW%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            j=bx%lo(2); k=bx%lo(3)
            do i=bx%lo(1),bx%hi(1)
               x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x_cc.ge.0.0_WP.and.x_cc.le.1.5_WP) then
                  nlocal=nlocal+1
                  vf_val=pVF(i,j,k,1)
                  local_x(nlocal)  =x_cc
                  local_VF(nlocal) =vf_val
                  local_RHO(nlocal)=vf_val*pRHOL(i,j,k,1)+(1.0_WP-vf_val)*pRHOG(i,j,k,1)
                  local_P(nlocal)  =vf_val*pPL(i,j,k,1)  +(1.0_WP-vf_val)*pPG(i,j,k,1)
                  local_UVW(nlocal)  =pUVW(i,j,k,1)
               end if
            end do
         end do
         call amr%mfiter_destroy(mfi)
         ! Gather counts from all ranks
         allocate(recvcount(amr%nproc),displs(amr%nproc))
         call MPI_GATHER(nlocal,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,0,amr%comm,ierr)
         ! Compute displacements and total count on root
         ntotal=0
         if (amr%amRoot) then
            displs(1)=0
            do i=2,amr%nproc
               displs(i)=displs(i-1)+recvcount(i-1)
            end do
            ntotal=displs(amr%nproc)+recvcount(amr%nproc)
         end if
         ! Allocate global arrays on root
         allocate(global_x(ntotal),global_VF(ntotal),global_RHO(ntotal),global_P(ntotal),global_UVW(ntotal))
         ! Gather data to root
         call MPI_GATHERV(local_x,  nlocal,MPI_REAL_WP,global_x,  recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_VF, nlocal,MPI_REAL_WP,global_VF, recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_RHO,nlocal,MPI_REAL_WP,global_RHO,recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_P,  nlocal,MPI_REAL_WP,global_P,  recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_UVW,  nlocal,MPI_REAL_WP,global_UVW,  recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         ! Sort by position on root (insertion sort — small array)
         if (amr%amRoot) then
            do i=2,ntotal
               if (global_x(i).lt.global_x(i-1)) then
                  ! Find insertion point
                  j=i-1
                  do while (j.ge.1.and.global_x(j).gt.global_x(i))
                     j=j-1
                  end do
                  j=j+1
                  ! Rotate element i into position j
                  vf_val=global_x(i);   global_x(j+1:i)  =global_x(j:i-1);   global_x(j)  =vf_val
                  vf_val=global_VF(i);  global_VF(j+1:i) =global_VF(j:i-1);  global_VF(j) =vf_val
                  vf_val=global_RHO(i); global_RHO(j+1:i)=global_RHO(j:i-1); global_RHO(j)=vf_val
                  vf_val=global_P(i);   global_P(j+1:i)  =global_P(j:i-1);   global_P(j)  =vf_val
                  vf_val=global_UVW(i);   global_UVW(j+1:i)  =global_UVW(j:i-1);   global_UVW(j)  =vf_val
               end if
            end do
         end if
         ! Only root outputs to file
         if (amr%amRoot) then
            if (.not.isdir('data')) call makedir('data')
            filename='data_'; write(timestamp,'(es12.5)') time%t
            open(newunit=iunit,file='data/'//trim(adjustl(filename))//trim(adjustl(timestamp))//'.csv', &
            &    form='formatted',status='replace',iostat=ierr)
            write(iunit,'(a)') 'Position,VF,RHO,P,Velocity'
            do i=1,ntotal
               write(iunit,'(g0.17,",",g0.17,",",g0.17,",",g0.17,",",g0.17)') &
               &    global_x(i),global_VF(i),global_RHO(i),global_P(i),global_UVW(i)
            end do
            close(iunit)
         end if
         ! Clean up
         deallocate(local_x,local_VF,local_RHO,local_P,local_UVW)
         deallocate(global_x,global_VF,global_RHO,global_P,global_UVW)
         deallocate(recvcount,displs)
      end block mixture_data

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
