!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mpcomp_class,      only: mpcomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Multiphase compressible flow solver and corresponding time tracker
   type(mpcomp),      public :: fs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ma
   
   !> Equation of state and shock properties
   real(WP) :: Gamma
   real(WP) :: M2,p2,rho2,a2,u2
   real(WP) :: M1,p1,rho1,a1,u1
      
contains
   
   
   !> Function that returns a smooth Heaviside representation of exact shock
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock
   
   
   !> P=EOS(RHO,I)
   pure real(WP) function get_P(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_P=RHO*I*(Gamma-1.0_WP)
   end function get_P
   !> C=f(RHO,P)
   pure real(WP) function get_C(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_C=sqrt(Gamma*P/RHO)
   end function get_C
   
   
   !> Mechanical relaxation model
   subroutine P_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: E
      ! Adjust conserved quantities to match single-phase model
      VF=Q(1)/(Q(1)+Q(2))
      E=Q(3)+Q(4)
      Q(3)=E*VF
      Q(4)=E*(1.0_WP-VF)
   end subroutine P_relax
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Prepare EoS and flow conditions
      initialize_eos: block
         use string,   only: str_long
         use messager, only: log
         character(str_long) :: message
         ! Read in gamma
         call param_read('Gamma',Gamma)
         ! Read in minimal shock information
         call param_read('Mach',M1)
         ! Generate preshock conditions
         rho1=1.0_WP
         p1=0.25_WP*rho1/gamma*((gamma+1.0_WP)*M1/(M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
         a1=sqrt(gamma*p1/rho1)
         ! Generate pre-shock velocity
         u1=M1*a1
         ! Also generate post-shock conditions
         p2=p1*(2.0_WP*gamma/(gamma+1.0_WP)*(M1**2-1.0_WP)+1.0_WP)
         rho2=rho1*(gamma+1.0_WP)*M1**2/((gamma-1.0_WP)*M1**2+2.0_WP)
         u2=u1*rho1/rho2
         a2=sqrt(gamma*p2/rho2)
         M2=u2/a2
         ! Output shock info
         if (cfg%amRoot) then
            write(message,'("[Pre -shock conditions] =>  rho1=",es12.5)')  rho1; call log(message)
            write(message,'("[Pre -shock conditions] =>    p1=",es12.5)')    p1; call log(message)
            write(message,'("[Pre -shock conditions] =>    u1=",es12.5)')    u1; call log(message)
            write(message,'("[Pre -shock conditions] =>    a1=",es12.5)')    a1; call log(message)
            write(message,'("[Pre -shock conditions] =>    M1=",es12.5)')    M1; call log(message)
            write(message,'("[Post-shock conditions] =>  rho2=",es12.5)')  rho2; call log(message)
            write(message,'("[Post-shock conditions] =>    p2=",es12.5)')    p2; call log(message)
            write(message,'("[Post-shock conditions] =>    u2=",es12.5)')    u2; call log(message)
            write(message,'("[Post-shock conditions] =>    a2=",es12.5)')    a2; call log(message)
            write(message,'("[Post-shock conditions] =>    M2=",es12.5)')    M2; call log(message)
         end if
      end block initialize_eos
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker
      
      ! Create multipgase compressible flow solver
      create_velocity_solver: block
         ! Initialize solver with required thermodynamic functions
         call fs%initialize(cfg=cfg,getPL=get_P,getCL=get_C,getPG=get_P,getCG=get_C,name='Compressible NS')
         ! Provide pressure relaxation model
         fs%Prelax=>P_relax
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dQdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
         allocate(Ma  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize our initial conditions
      initial_conditions: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use mms_geom,              only: initialize_volume_moments
         use mpcomp_class,          only: VFlo
         integer :: i,j,k
         ! Initialize gas fields
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Set VF to zero
                  fs%VF(i,j,k)=0.0_WP
                  fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  ! Set corresponding PLIC
                  call setNumberOfPlanes(fs%PLIC(i,j,k),1); call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,fs%VF(i,j,k)-0.5_WP))
                  ! Volume moments for a slab
                  call initialize_volume_moments(lo=[cfg%x(i),cfg%y(j),cfg%z(k)],hi=[cfg%x(i+1),cfg%y(j+1),cfg%z(k+1)],&
                  levelset=levelset_slab,time=0.0_WP,level=4,VFlo=VFlo,VF=fs%VF(i,j,k),BL=fs%BL(:,i,j,k),BG=fs%BG(:,i,j,k))
                  ! Setup normal shock
                  fs%RHOG(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i),delta=0.5_WP)
                  fs%RHOL(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i),delta=0.5_WP)
                  fs%U   (i,j,k)=u1  +(u2  -u1  )*Hshock(cfg%xm(i),delta=0.5_WP)
                  fs%PG  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i),delta=0.5_WP)
                  fs%PL  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i),delta=0.5_WP)
                  fs%IG  (i,j,k)=fs%PG(i,j,k)/(fs%RHOG(i,j,k)*(Gamma-1.0_WP))
                  fs%IL  (i,j,k)=fs%PL(i,j,k)/(fs%RHOL(i,j,k)*(Gamma-1.0_WP))
               end do
            end do
         end do
         ! Uncomment this for moving shock
         fs%U=fs%U-u1
         ! Build total energy
         fs%EL=fs%IL+0.5_WP*fs%U**2
         fs%EG=fs%IG+0.5_WP*fs%U**2
         ! Initialize conserved variables
         fs%Q(:,:,:,1)=(       fs%VF)*fs%RHOL
         fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
         fs%Q(:,:,:,3)=fs%Q(:,:,:,1)*fs%EL
         fs%Q(:,:,:,4)=fs%Q(:,:,:,2)*fs%EG
         fs%RHO=fs%Q(:,:,:,1)+fs%Q(:,:,:,2)
         fs%Q(:,:,:,5)=fs%RHO*fs%U
         fs%Q(:,:,:,6)=fs%RHO*fs%V
         fs%Q(:,:,:,7)=fs%RHO*fs%W
         ! Rebuild primitive variables
         call fs%get_primitive()
         ! Compute local Mach number
         Ma=sqrt(fs%U**2+fs%V**2+fs%W**2)/fs%C
      end block initial_conditions
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='staticshock')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('VOF',fs%VF)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_vector('U',fs%U,fs%V,fs%W)
         call ens_out%add_scalar('RHO' ,fs%RHO)
         call ens_out%add_scalar('E',fs%E)
         call ens_out%add_scalar('Mach',Ma)

         call ens_out%add_scalar('PL',fs%PL)
         call ens_out%add_scalar('PG',fs%PG)
         call ens_out%add_scalar('RHOL' ,fs%RHOL)
         call ens_out%add_scalar('RHOG' ,fs%RHOG)
         call ens_out%add_scalar('EL',fs%EL)
         call ens_out%add_scalar('EG',fs%EG)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call fs%get_info()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOGmax,'max(RHO)')
         call mfile%add_column(fs%RHOGmin,'min(RHO)')
         call mfile%add_column(fs%EGmax  ,'max(E)'  )
         call mfile%add_column(fs%EGmin  ,'min(E)'  )
         call mfile%add_column(fs%PGmax  ,'max(P)'  )
         call mfile%add_column(fs%PGmin  ,'min(P)'  )
         call mfile%add_column(fs%TGmax  ,'max(T)'  )
         call mfile%add_column(fs%TGmin  ,'min(T)'  )
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(2),'Mass')
         call consfile%add_column(fs%Qint(4),'Total energy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%RHOKGint,'Kinetic energy')
         call consfile%write()
      end block create_monitor
      
   contains
      
      !> Function that defines a level set function for a slab
      function levelset_slab(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=20.0_WP-abs(xyz(1)+40.0_WP)
      end function levelset_slab
      
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
         
         ! Remember conserved variables
         fs%Qold=fs%Q
         
         ! Remember phasic quantities
         fs%RHOLold=fs%RHOL; fs%ILold=fs%IL
         fs%RHOGold=fs%RHOG; fs%IGold=fs%IG
         
         ! Remember volume moments and interface
         fs%VFold=fs%VF
         fs%BLold=fs%BL
         fs%BGold=fs%BG
         copy_plic_to_old: block
            use irl_fortran_interface, only: copy
            integer :: i,j,k
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
               call copy(fs%PLICold(i,j,k),fs%PLIC(i,j,k))
            end do; end do; end do
         end block copy_plic_to_old
         
         ! Tag cells for semi-Lagrangian transport
         call fs%SLtag()
         
         ! Perform first semi-Lagrangian transport step =====================================================
         call fs%SLstep(dt=0.5_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         call fs%build_interface()
         
         ! First RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,1))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Enforce mechanical equilibrium
         call fs%relax_pressure()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Second RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,2))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Enforce mechanical equilibrium
         call fs%relax_pressure()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Perform second semi-Lagrangian transport step ====================================================
         call fs%SLstep(dt=1.0_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         call fs%build_interface()
         
         ! Third RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,3))
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Enforce mechanical equilibrium
         call fs%relax_pressure()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Fourth RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,4))
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Enforce mechanical equilibrium
         call fs%relax_pressure()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Compute local Mach number
         Ma=sqrt(fs%U**2+fs%V**2+fs%W**2)/fs%C
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(dQdt,Ma)
   end subroutine simulation_final
   
   
end module simulation
