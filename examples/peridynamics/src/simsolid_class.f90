!> Definition of a solid simulation superclass that drives an lss_class object
module simsolid_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use lss_class,         only: lss
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   
   public :: simsolid
   

   !> Solid simulation object
   type :: simsolid
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config
      type(config) :: cfg
      
      !> Get Lagrangian solid solver and timetracker
      type(lss),         public :: ls
      type(timetracker), public :: time
      
      !> Ensight postprocessing
      type(partmesh) :: pmesh
      type(ensight)  :: ens_out
      type(event)    :: ens_evt
      
      !> Monitor file
      type(monitor) :: mfile
      
   contains
      procedure :: init                    !< Initialize solid simulation
      procedure :: step                    !< Advance solid simulation by one time step
      procedure :: final                   !< Finalize solid simulation
   end type simsolid
   
   
contains
   
   
   !> Initialization of solid simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(simsolid), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_solid')
      
      ! Create a config
      initialize_cfg: block
         !this%cfg=config
      end block initialize_cfg
      
      ! Initialize time tracker
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
      end block initialize_timetracker

      ! Initialize Lagrangian solid solver
      initialize_lss: block
         use mathtools, only: twoPi,Pi
         use random,    only: random_uniform
         integer :: i,j,k,np
         real(WP) :: radius
         
         ! Create solver
         this%ls=lss(cfg=this%cfg,name='solid')
         
         ! Set material properties
         this%ls%elastic_modulus=1.0e5_WP
         this%ls%poisson_ratio=0.333_WP
         this%ls%rho=1000.0_WP
         
         ! Read sphere properties
         call this%input%read('Number of markers',np)
         ! Root process initializes marker particles
         if (this%ls%cfg%amRoot) then
            call this%ls%resize(np)
            ! Distribute marker particles
            do i=1,np
               ! Set object id
               this%ls%p(i)%id=1
               ! Give zero dt
               this%ls%p(i)%dt=1.0e-3_WP
               ! Set position
               radius=huge(1.0_WP)
               do while (radius.gt.0.5_WP)
                  this%ls%p(i)%pos=[random_uniform(-0.5_WP,+0.5_WP),random_uniform(-0.5_WP,+0.5_WP),random_uniform(-0.5_WP,+0.5_WP)]
                  radius=sqrt(dot_product(this%ls%p(i)%pos,this%ls%p(i)%pos))
               end do
               ! Assign velocity
               this%ls%p(i)%vel=[0.1_WP*this%ls%p(i)%pos(1),0.0_WP,0.0_WP]
               ! Assign element volume
               this%ls%p(i)%dV=Pi/(6.0_WP*real(np,WP))
               ! Locate the particle on the mesh
               this%ls%p(i)%ind=this%ls%cfg%get_ijk_global(this%ls%p(i)%pos,[this%ls%cfg%imin,this%ls%cfg%jmin,this%ls%cfg%kmin])
               ! Assign a unique integer to particle
               this%ls%p(i)%i=i
               ! Activate the particle
               this%ls%p(i)%flag=0
            end do
         end if
         call this%ls%sync()
         
         ! Initalize bonds
         call this%ls%bond_init()
         
      end block initialize_lss


      ! Create partmesh object for visualizing Lagrangian particles
      create_pmesh: block
         use lss_class, only: max_bond
         integer :: i,n
         this%pmesh=partmesh(nvar=2,nvec=2,name='solid')
         this%pmesh%varname(1)='nbond'
         this%pmesh%varname(2)='dilatation'
         this%pmesh%vecname(1)='velocity'
         this%pmesh%vecname(2)='bond_force'
         call this%ls%update_partmesh(this%pmesh)
         do i=1,this%ls%np_
            this%pmesh%var(1,i)=0.0_WP
            do n=1,max_bond
               if (this%ls%p(i)%ibond(n).gt.0) this%pmesh%var(1,i)=this%pmesh%var(1,i)+1.0_WP
            end do
            this%pmesh%var(2,i)  =this%ls%p(i)%dil
            this%pmesh%vec(:,1,i)=this%ls%p(i)%vel
            this%pmesh%vec(:,2,i)=this%ls%p(i)%Abond
         end do
      end block create_pmesh


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='solid')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         call this%ens_out%add_particle('particles',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%ls%get_cfl(this%time%dt,this%time%cfl)
         call this%ls%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%cfg%amRoot,'simulation_solid')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%ls%np,'Particle number')
         call this%mfile%add_column(this%ls%Umin,'Particle Umin')
         call this%mfile%add_column(this%ls%Umax,'Particle Umax')
         call this%mfile%add_column(this%ls%Vmin,'Particle Vmin')
         call this%mfile%add_column(this%ls%Vmax,'Particle Vmax')
         call this%mfile%add_column(this%ls%Wmin,'Particle Wmin')
         call this%mfile%add_column(this%ls%Wmax,'Particle Wmax')
         call this%mfile%write()
      end block create_monitor

      
   end subroutine init
   
   
   !> Take one time step of the solid simulation
	subroutine step(this)
		implicit none
		class(simsolid), intent(inout) :: this
      
      ! Advance solid solver
      call this%ls%advance(this%time%dt)
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         update_pmesh: block
            use lss_class, only: max_bond
            integer :: i,n
            call this%ls%update_partmesh(this%pmesh)
            do i=1,this%ls%np_
               this%pmesh%var(1,i)=0.0_WP
               do n=1,max_bond
                  if (this%ls%p(i)%ibond(n).gt.0) this%pmesh%var(1,i)=this%pmesh%var(1,i)+1.0_WP
               end do
               this%pmesh%var(2,i)  =this%ls%p(i)%dil
               this%pmesh%vec(:,1,i)=this%ls%p(i)%vel
               this%pmesh%vec(:,2,i)=this%ls%p(i)%Abond
            end do
         end block update_pmesh
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%ls%get_max()
      call this%mfile%write()
      
   end subroutine step
   
   
   !> Finalize solid simulation
   subroutine final(this)
		implicit none
		class(simsolid), intent(inout) :: this
		
		! Deallocate work arrays
		
	end subroutine final
   
   
end module simsolid_class