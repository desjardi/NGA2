!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrmpinc_class,    only: amrmpinc
   use amrdata_class,     only: amrdata
   use amrio_class,       only: amrio
   use amrcclabel_class,  only: amrcclabel
   use amrvof_class,      only: amrvof
   use monitor_class,     only: monitor

   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! VOF solver
   type(amrvof), target :: vof

   ! Ellipsoid parameters
   integer :: nEllipsoid
   real(WP), dimension(:,:), allocatable :: ellipsoid_center
   real(WP), dimension(:,:), allocatable :: ellipsoid_radius

   ! Visualization 
   type(amrviz) :: viz 

   ! I/O 
   type(amrio) :: io

   ! CCLabel
   type(amrcclabel) :: cclabel

   ! Monitoring
   type(monitor) :: gridfile

contains

   !> Function that identifies cells within a structure
   logical function make_label(pVF,lo,i,j,k)
      implicit none
      real(WP), dimension(:,:,:,:), intent(in) :: pVF
      integer, dimension(3), intent(in) :: lo
      integer, intent(in) :: i,j,k
      integer :: il,jl,kl
      il = i - lo(1) + 1
      jl = j - lo(2) + 1
      kl = k - lo(3) + 1
      if (pVF(il,jl,kl,1).gt.0.0_WP) then
         make_label=.true.
      else
         make_label=.false.
      end if
   end function make_label

   !> Function that identifies if neighbors are within the same structure
   logical function same_label(pVF,lo,i,j,k,ii,jj,kk)
      implicit none
      real(WP), dimension(:,:,:,:), intent(in) :: pVF
      integer, dimension(3), intent(in) :: lo
      integer, intent(in) :: i,j,k,ii,jj,kk
      integer :: il,jl,kl,iil,jjl,kkl
      il  = i  - lo(1) + 1
      jl  = j  - lo(2) + 1
      kl  = k  - lo(3) + 1
      iil = ii - lo(1) + 1
      jjl = jj - lo(2) + 1
      kkl = kk - lo(3) + 1
      if (pVF(il,jl,kl,1).gt.0.0_WP .and. pVF(iil,jjl,kkl,1).gt.0.0_WP) then
         same_label=.true.
      else
         same_label=.false.
      end if
   end function same_label

   !> Function that identifies cells within a structure on coarse level
   logical function coarse_make_label(pVF,lo,i,j,k)
      use amrmpinc_class,   only: VFhi
      implicit none
      real(WP), dimension(:,:,:,:), intent(in) :: pVF
      integer, dimension(3), intent(in) :: lo
      integer, intent(in) :: i,j,k
      integer :: il,jl,kl
      il = i - lo(1) + 1
      jl = j - lo(2) + 1
      kl = k - lo(3) + 1
      if (pVF(il,jl,kl,1).gt.VFhi) then
         coarse_make_label=.true.
      else
         coarse_make_label=.false.
      end if
   end function coarse_make_label

   !> Function that identifies if neighbors are within the same structure on coarse level
   logical function coarse_same_label(pVF,lo,i,j,k,ii,jj,kk)
      use amrmpinc_class,   only: VFhi
      implicit none
      real(WP), dimension(:,:,:,:), intent(in) :: pVF
      integer, dimension(3), intent(in) :: lo
      integer, intent(in) :: i,j,k,ii,jj,kk
      integer :: il,jl,kl,iil,jjl,kkl
      il  = i  - lo(1) + 1
      jl  = j  - lo(2) + 1
      kl  = k  - lo(3) + 1
      iil = ii - lo(1) + 1
      jjl = jj - lo(2) + 1
      kkl = kk - lo(3) + 1
      if (pVF(il,jl,kl,1).gt.VFhi .and. pVF(iil,jjl,kkl,1).gt.VFhi) then
         coarse_same_label=.true.
      else
         coarse_same_label=.false.
      end if
   end function coarse_same_label

   !> Ellipsoids levelset function with periodicity
   function Ellipsoids_levelset(xyz,t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,phi
      real(WP), dimension(3) :: d,L
      integer :: n
      ! Distance to nearest Ellipsoid 
      G = -huge(1.0_WP)
      do n=1,nEllipsoid
         d=xyz-ellipsoid_center(:,n)
         L=[amr%xhi-amr%xlo,amr%yhi-amr%ylo,amr%zhi-amr%zlo]
         d=d-L*nint(d/L)  ! Nearest image
         phi = 1.0_WP - sqrt( &
              (d(1)/ellipsoid_radius(1,n))**2 + &
              (d(2)/ellipsoid_radius(2,n))**2 + &
              (d(3)/ellipsoid_radius(3,n))**2 )
         G=max(G,phi)
      end do
   end function Ellipsoids_levelset

   !> Initialize VF field with Ellipsoids using levelset-based moments
   subroutine Ellipsoids_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom,         only: initialize_volume_moments
      use amrmpinc_class,   only: VFlo,VFhi
      implicit none
      class(amrvof), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG
      real(WP) :: dx,dy,dz
      real(WP), dimension(3) :: BL,BG
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl)
      dy=solver%amr%dy(lvl)
      dz=solver%amr%dz(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Get tile box with ghost cells
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF and barycenters from levelset with 3 levels of refinement
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=Ellipsoids_levelset,time=time,level=nref,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine Ellipsoids_init
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         amr%name='cclabel_tester'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=0.0_WP; amr%xhi=1.0_WP
         amr%ylo=0.0_WP; amr%yhi=1.0_WP
         amr%zlo=0.0_WP; amr%zhi=1.0_WP
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid

      ! Setup Ellipsoids parameters
      setup_Ellipsoids: block
         use random, only: random_uniform
         use string,   only: str_medium
         use messager,       only: die
         integer :: nD,nseed
         real(WP), dimension(3) :: center,radius
         real(WP) :: radius_scale
         integer :: myseed
         integer, dimension(:), allocatable :: seed
         character(len=str_medium) :: case
         call param_read('Droplet case',case,default='Random')
         if (case == 'Random') then 
            call param_read('Number of ellipsoids',nEllipsoid,default=4)
            call param_read('Random seed',myseed,default=3)
            call param_read('Radius scale',radius_scale,default=0.5_WP)
            ! Provide seed for random number generator
            call random_seed(size=nseed)
            allocate(seed(nseed))
            seed(:)=myseed
            call random_seed(put=seed)
         else if (case == 'Cylinder') then
            nEllipsoid=1
         end if
         
         ! Allocate arrays
         allocate(ellipsoid_center(3,nEllipsoid))
         allocate(ellipsoid_radius(3,nEllipsoid))

         ! Define centers and radii of ellipsoids
         do nD=1,nEllipsoid
            if (case == 'Random') then
               ! Random center and radius
               center=[random_uniform(amr%xlo, amr%xhi), &
                       random_uniform(amr%ylo, amr%yhi), &
                       random_uniform(amr%zlo, amr%zhi)  ]
               radius=[radius_scale*random_uniform(amr%xlo, amr%xhi), &
                       radius_scale*random_uniform(amr%ylo, amr%yhi), &
                       radius_scale*random_uniform(amr%zlo, amr%zhi)  ]
            else if (case == 'Cylinder') then
               ! Large cylinder for testing multiple levels representing one structure
               center=[0.5_WP, 0.5_WP,0.5_WP]
               radius=[0.4_WP,10.0_WP,0.4_WP]
            else
               call die('Unknown droplet case')
            end if

            ellipsoid_center(:,nD)=center
            ellipsoid_radius(:,nD)=radius
         end do
      end block setup_Ellipsoids

      ! Initialize our VOF field
      create_and_initialize_vof: block
         call vof%initialize(amr,name='Ellipsoids_vof')
         vof%user_vof_init=>Ellipsoids_init
      end block create_and_initialize_vof

      ! Initialize CCLabel
      create_and_initialize_cclabel: block
         call cclabel%initialize(amr,name='Ellipsoids_cclabel')
      end block create_and_initialize_cclabel

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Fresh start
         call amr%init_from_scratch(time=0.0_WP)
         ! Build PLIC
         call vof%build_plic(0.0_WP)
      end block init_regridding

      ! Create visualization
      create_visualization: block
         ! Create amrviz output
         call viz%initialize(amr,'amrcclabel',use_hdf5=.false.)
         call viz%add_scalar(vof%VF,1,'VF')
         call viz%add_scalar(cclabel%id,1,'ID')
         call viz%add_surfmesh(vof%smesh,'plic')
      end block create_visualization

      ! Create monitor
      create_monitor: block
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitor

      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
     
      ! Compute CCLabel
      call cclabel%build(make_label,same_label,coarse_make_label,coarse_same_label,vof%VF)

      ! Write visualization with IDs
      call viz%write(time=0.0_WP)
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays
      call io%finalize()
      call amr%finalize()
      call vof%finalize()
      
   end subroutine simulation_final
   
end module simulation
