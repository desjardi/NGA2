!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Only get a LPT solver and corresponding time tracker
   type(lpt),         public :: lp
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Fluid info
   real(WP), dimension(:,:,:), allocatable :: U,V,W
   real(WP), dimension(:,:,:), allocatable :: rho,visc

   !> Injection parameters
   real(WP) :: mfr,inj_dmean,inj_dsd,inj_dmin,inj_dmax,inj_dshift,inj_duration
   real(WP), dimension(3) :: inj_vel,inj_pos1,inj_pos2

   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker
      
      
      ! Initialize our LPT solver
      initialize_lpt: block
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Turn off drag
         lp%drag_model='Stokes'
         ! Initialize with zero particles
         call lp%resize(0)
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n,default=0.7_WP)
         call param_read('Wall restitution',lp%e_w,default=lp%e_n)
         call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
      end block initialize_lpt
      
      
      ! Set injection parameters
      injection_parameters: block
         real(WP) :: pmin,pmax
         call param_read('Particle inject duration'   ,inj_duration,default=huge(1.0_WP))
         call param_read('Particle inject point 1'    ,inj_pos1)
         call param_read('Particle inject point 2'    ,inj_pos2)
         pmin=min(inj_pos1(1),inj_pos2(1)); pmax=max(inj_pos1(1),inj_pos2(1)); inj_pos1(1)=pmin; inj_pos2(1)=pmax
         pmin=min(inj_pos1(2),inj_pos2(2)); pmax=max(inj_pos1(2),inj_pos2(2)); inj_pos1(2)=pmin; inj_pos2(2)=pmax
         pmin=min(inj_pos1(3),inj_pos2(3)); pmax=max(inj_pos1(3),inj_pos2(3)); inj_pos1(3)=pmin; inj_pos2(3)=pmax
         call param_read('Particle mass flow rate'    ,mfr)
         call param_read('Particle velocity'          ,inj_vel)
         call param_read('Particle mean diameter'     ,inj_dmean)
         call param_read('Particle standard deviation',inj_dsd   ,default=0.0_WP)
         call param_read('Particle min diameter'      ,inj_dmin  ,default=tiny(1.0_WP))
         call param_read('Particle max diameter'      ,inj_dmax  ,default=huge(1.0_WP))
         call param_read('Particle diameter shift'    ,inj_dshift,default=0.0_WP)
         if (inj_dsd.eq.0.0_WP) then
            inj_dmin=inj_dmean-epsilon(1.0_WP)
            inj_dmax=inj_dmean+epsilon(1.0_WP)
         end if
      end block injection_parameters
      
      
      ! Restart bed here
      read_bed: block
         integer :: i
         ! Reset injection duration to zero
         inj_duration=0.0_WP
         ! Restart from bed file
         call lp%read(filename='bed.file')
         ! Loop through particles
         do i=1,lp%np_
            ! Remove top of bed
            !if (lp%p(i)%pos(2).gt.0.4_WP) lp%p(i)%flag=1
            ! Zero out velocity
            lp%p(i)%vel=0.0_WP
         end do
         call lp%sync()
         ! Recalculate VF
         call lp%update_VF()
      end block read_bed
      
      
      ! Prepare our fluid phase info based on a Taylor vortex
      initialize_fluid: block
         real(WP) :: rhof,viscf
         ! Allocate arrays
         allocate(rho (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); rho =1.0_WP
         allocate(visc(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); visc=0.0_WP
         allocate(U   (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); U   =0.0_WP
         allocate(V   (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); V   =0.0_WP
         allocate(W   (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); W   =0.0_WP
         ! Set constant density and viscosity
         call param_read('Density'  ,rhof ); rho =rhof
         call param_read('Viscosity',viscf); visc=viscf
      end block initialize_fluid
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=2,nvec=2,name='lpt')
         pmesh%varname(1)='radius'
         pmesh%varname(2)='id'
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='Fcol'
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(  1,i)=0.5_WP*lp%p(i)%d
            pmesh%var(  2,i)=real(lp%p(i)%id,WP)
            pmesh%vec(:,1,i)=lp%p(i)%vel
            pmesh%vec(:,2,i)=lp%p(i)%Acol
         end do
      end block create_pmesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=lp%cfg,name='particles')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call lp%get_max()
         ! Create simulation monitor
         mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(lp%np,'Particle number')
         call mfile%add_column(lp%np_new,'Npart new')
         call mfile%add_column(lp%np_out,'Npart removed')
         call mfile%add_column(lp%ncol,'Particle collisions')
         call mfile%add_column(lp%VFmax,'Max VF')
         call mfile%add_column(lp%Umin,'Particle Umin')
         call mfile%add_column(lp%Umax,'Particle Umax')
         call mfile%add_column(lp%Vmin,'Particle Vmin')
         call mfile%add_column(lp%Vmax,'Particle Vmax')
         call mfile%add_column(lp%Wmin,'Particle Wmin')
         call mfile%add_column(lp%Wmax,'Particle Wmax')
         call mfile%add_column(lp%dmin,'Particle dmin')
         call mfile%add_column(lp%dmax,'Particle dmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(lp%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(lp%CFLp_x,'Particle xCFL')
         call cflfile%add_column(lp%CFLp_y,'Particle yCFL')
         call cflfile%add_column(lp%CFLp_z,'Particle zCFL')
         call cflfile%add_column(lp%CFL_col,'Collision CFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Inject particles
         if (time%t.lt.inj_duration) call inject_particles(dt=time%dt,avoid_overlap=.true.)
         
         ! Collide particles
         call lp%collide(dt=time%dt)
         
         ! Advance particles by dt
         call lp%advance(dt=time%dt,U=U,V=V,W=W,rho=rho,visc=visc)
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(  1,i)=0.5_WP*lp%p(i)%d
                  pmesh%var(  2,i)=real(lp%p(i)%id,WP)
                  pmesh%vec(:,1,i)=lp%p(i)%vel
                  pmesh%vec(:,2,i)=lp%p(i)%Acol
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         
      end do
      
      ! Write out the particles
      call lp%write(filename='bed.file')
      
   end subroutine simulation_run
   
   
   !> Inject particles from a prescribed location with a given mass flowrate
   !> Requires injection parameters to be set
   subroutine inject_particles(dt,avoid_overlap)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX,MPI_SUM,MPI_IN_PLACE,MPI_INTEGER,MPI_INTEGER8,MPI_Status
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      use lpt_class, only: part,MPI_PART
      use precision
      implicit none
      real(WP), intent(inout) :: dt                  !< Timestep size over which to advance
      logical, intent(in), optional :: avoid_overlap !< Option to avoid overlap during injection
      real(WP), dimension(3) :: inj_min,inj_max      !< Min/max extents of injection
      real(WP) :: Mgoal,Madded,Mtmp                  !< Mass flow rate parameters
      real(WP), save :: previous_error=0.0_WP        !< Store mass left over from previous timestep
      integer(kind=I8) :: maxid                      !< Keep track of maximum particle id
      integer :: i,j,np0_,np2,np_tmp,count,ierr
      integer, dimension(:), allocatable :: nrecv
      type(part), dimension(:), allocatable :: p2
      type(MPI_Status) :: status
      logical :: avoid_overlap_,overlap
      
      ! Initial number of particles
      np0_=lp%np_
      lp%np_new=0
      
      ! Get the particle mass that should be added to the system
      Mgoal =mfr*dt+previous_error
      Madded=0.0_WP
      
      ! Determine id to assign to particle
      maxid=0
      do i=1,lp%np_
         maxid=max(maxid,lp%p(i)%id)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,maxid,1,MPI_INTEGER8,MPI_MAX,lp%cfg%comm,ierr)
      
      ! Communicate nearby particles to check for overlap
      avoid_overlap_=.false.
      if (present(avoid_overlap)) avoid_overlap_=avoid_overlap
      if (avoid_overlap_) then
         allocate(nrecv(lp%cfg%nproc))
         count=0
         inj_min=inj_pos1-inj_dmax
         inj_max=inj_pos2+inj_dmax
         do i=1,lp%np_
            if (lp%p(i)%pos(1).gt.inj_min(1).and.lp%p(i)%pos(1).lt.inj_max(1).and.&
            &   lp%p(i)%pos(2).gt.inj_min(2).and.lp%p(i)%pos(2).lt.inj_max(2).and.&
            &   lp%p(i)%pos(3).gt.inj_min(3).and.lp%p(i)%pos(3).lt.inj_max(3)) count=count+1
         end do
         call MPI_GATHER(count,1,MPI_INTEGER,nrecv,1,MPI_INTEGER,0,lp%cfg%comm,ierr)
         if (lp%cfg%amRoot) then
            np2=sum(nrecv)
            allocate(p2(np2))
         else
            allocate(p2(count))
         end if
         count=0
         do i=1,lp%np_
            if (lp%p(i)%pos(1).gt.inj_min(1).and.lp%p(i)%pos(1).lt.inj_max(1).and.&
            &   lp%p(i)%pos(2).gt.inj_min(2).and.lp%p(i)%pos(2).lt.inj_max(2).and.&
            &   lp%p(i)%pos(3).gt.inj_min(3).and.lp%p(i)%pos(3).lt.inj_max(3)) then
               count=count+1
               p2(count)=lp%p(i)
            end if
         end do
         if (lp%cfg%amRoot) then
            do i=2,lp%cfg%nproc
               if (nrecv(i).gt.0) call MPI_recv(p2(sum(nrecv(1:i-1))+1:sum(nrecv(1:i))),nrecv(i),MPI_PART,i-1,0,lp%cfg%comm,status,ierr)
            end do
         else
            if (count.gt.0) call MPI_send(p2,count,MPI_PART,0,0,lp%cfg%comm,ierr)
         end if
         deallocate(nrecv)
      end if
      
      ! Add new particles until desired mass is achieved
      do while (Madded.lt.Mgoal)
         ! Root creates generates new particles
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
               ! Generate a diameter
               lp%p(count)%d=get_diameter()
               ! Set various parameters for the particle
               lp%p(count)%id    =maxid+int(np_tmp,8)
               lp%p(count)%dt    =0.0_WP
               lp%p(count)%Acol  =0.0_WP
               lp%p(count)%Tcol  =0.0_WP
               lp%p(count)%angVel=0.0_WP
               ! Give a position at the injector to the particle
               lp%p(count)%pos=get_position(lp%p(count)%d)
               overlap=.false.
               ! Check overlap with particles recently injected
               if (avoid_overlap_) then
                  do j=1,np_tmp-1
                     if (norm2(lp%p(count)%pos-lp%p(j)%pos).lt.0.5_WP*(lp%p(count)%d+lp%p(j)%d)) overlap=.true.
                  end do
                  ! Check overlap with all other particles
                  if (.not.overlap) then
                     do j=1,np2
                        if (norm2(lp%p(count)%pos-p2(j)%pos).lt.0.5_WP*(lp%p(count)%d+p2(j)%d)) overlap=.true.
                     end do
                  end if
               end if
               if (overlap) then
                  ! Try again
                  np_tmp=np_tmp-1
               else
                  ! Localize the particle
                  lp%p(count)%ind(1)=lp%cfg%imin; lp%p(count)%ind(2)=lp%cfg%jmin; lp%p(count)%ind(3)=lp%cfg%kmin
                  lp%p(count)%ind=lp%cfg%get_ijk_global(lp%p(count)%pos,lp%p(count)%ind)
                  ! Give it a velocity
                  lp%p(count)%vel=inj_vel
                  ! Make it an "official" particle
                  lp%p(count)%flag=0
                  ! Update the added mass for the timestep
                  Mtmp=Mtmp+lp%rho*Pi/6.0_WP*lp%p(count)%d**3
               end if
            end do
         end if
         ! Communicate particles
         call lp%sync()
         ! Loop through newly created particles
         Mtmp=0.0_WP
         do i=np0_+1,lp%np_
            ! Remove if out of bounds
            if (lp%cfg%VF(lp%p(i)%ind(1),lp%p(i)%ind(2),lp%p(i)%ind(3)).le.0.0_WP) lp%p(i)%flag=1
            if (lp%p(i)%flag.eq.0) then
               ! Update the added mass for the timestep
               Mtmp=Mtmp+lp%rho*Pi/6.0_WP*lp%p(i)%d**3
               ! Update the max particle id
               maxid=max(maxid,lp%p(i)%id)
               ! Increment counter
               lp%np_new=lp%np_new+1
            end if
         end do
         ! Total mass added
         call MPI_ALLREDUCE(MPI_IN_PLACE,Mtmp,1,MPI_REAL_WP,MPI_SUM,lp%cfg%comm,ierr); Madded=Madded+Mtmp
         ! Clean up particles
         call lp%recycle()
         ! Update initial npart
         np0_=lp%np_
         ! Maximum particle id
         call MPI_ALLREDUCE(MPI_IN_PLACE,maxid,1,MPI_INTEGER8,MPI_MAX,lp%cfg%comm,ierr)
      end do
      
      ! Remember the error
      previous_error=Mgoal-Madded
      
      ! Sum up injected particles
      call MPI_ALLREDUCE(MPI_IN_PLACE,lp%np_new,1,MPI_INTEGER,MPI_SUM,lp%cfg%comm,ierr)
      
      ! Deallocate
      if (allocated(p2)) deallocate(p2)
      
   contains
      
      ! Compute particle diameter from a lognormal distribution
      function get_diameter() result(dp)
         use random, only: random_lognormal
         implicit none
         real(WP) :: dp
         dp=-1.0_WP
         do while (dp.lt.inj_dmin.or.dp.gt.inj_dmax)
            dp=random_lognormal(m=inj_dmean-inj_dshift,sd=inj_dsd)+inj_dshift
         end do
      end function get_diameter
      
      ! Position for bulk injection of particles
      function get_position(d) result(pos)
         use random,    only: random_uniform
         use mathtools, only: twoPi
         implicit none
         real(WP), intent(in)   :: d
         real(WP), dimension(3) :: pos
         pos(1)=random_uniform(lo=inj_pos1(1)+0.5_WP*d,hi=inj_pos2(1)-0.5_WP*d)
         pos(2)=inj_pos1(2)
         pos(3)=random_uniform(lo=inj_pos1(3)+0.5_WP*d,hi=inj_pos2(3)-0.5_WP*d)
      end function get_position
      
   end subroutine inject_particles
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(rho,visc,U,V,W)
      
   end subroutine simulation_final
   
   
end module simulation
