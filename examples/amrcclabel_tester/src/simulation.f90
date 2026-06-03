!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrmpinc_class,    only: amrmpinc
   use amrdata_class,     only: amrdata
   use amrio_class,       only: amrio
   use amrcclabel_class,  only: amrcclabel

   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! VOF solver
   type(amrvof), target :: vof

   ! Sphere parameters
   integer :: nSphere
   real(WP), dimension(:,:), allocatable :: sphere_center
   real(WP), dimension(:), allocatable :: sphere_radius

   ! Visualization 
   type(amrviz) :: viz 

   ! I/O 
   type(amrio) :: io

   ! CCLabel
   type(amrcclabel) :: cclabel

contains

   
   !> Function that identifies cells within a structure
   logical function make_label(pVF,i,j,k)
      implicit none
      real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pVF
      integer, intent(in) :: i,j,k
      if (pVF(i,j,k,1).gt.0.0_WP) then
         make_label=.true.
      else
         make_label=.false.
      end if
   end function make_label

    !> Function that identifies if neighbors are within the same structure
   logical function same_label(pVF,i,j,k,ii,jj,kk)
       implicit none
       real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pVF
       integer, intent(in) :: i,j,k,ii,jj,kk
       if (pVF(i,j,k,1).gt.0.0_WP .and. pVF(ii,jj,kk,1).gt.0.0_WP) then
          same_label=.true.
       else
          same_label=.false.
       end if
   end function same_label

     !> Spheres levelset function with periodicity
   function spheres_levelset(xyz,t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP), dimension(3) :: d,L
      ! Distance to nearest sphere 
      do n=1,nSphere
         d=xyz-sphere_center(:,n)
         L=[amr%xhi-amr%xlo,amr%yhi-amr%ylo,amr%zhi-amr%zlo]
         d=d-L*nint(d/L)  ! Nearest image
         G=min(G,sphere_radius(n)-sqrt(sum(d**2)))
      end do
   end function spheres_levelset

   !> Initialize VF field with sphere using levelset-based moments
   subroutine spheres_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom,         only: initialize_volume_moments
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
            &                              levelset=spheres_levelset,time=time,level=nref,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine spheres_init
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         amr%name='vof_advect'
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

      ! Setup spheres parameters
      setup_spheres: block
         use random, only: random_uniform
         call param_read('Sphere diameter',radius); radius=radius/2.0_WP
         call param_read('Number of spheres',nSphere);
         ! Allocate arrays
         allocate(sphere_center(3,nSphere))
         allocate(sphere_radius(nSphere))
         ! Provide seed for random number generator
         call random_seed(size=nseed)
         allocate(seed(nseed))
         seed(:)=1
         call random_seed(put=seed)
         do nD=1,nSphere
            center=[random_uniform(amr%xlo, amr%xhi), &
                    random_uniform(amr%ylo, amr%yhi), &
                    random_uniform(amr%zlo, amr%zhi)  ]
            sphere_center(:,nD)=center
            sphere_radius(nD)=radius
         end do
      end block setup_spheres

      ! Initialize our VOF field
      create_and_initialize_vof: block
         call vof%initialize(amr,name='spheres_vof')
         vof%user_vof_init=>spheres_init
      end block create_and_initialize_vof


      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Fresh start
         call amr%init_from_scratch(time=0.0_WP)
         ! Build PLIC
         call vof%build_plic(time%t)
      end block init_regridding

      ! Create visualization
      create_visualization: block
         ! Create amrviz output
         call viz%initialize(amr,'amrcclabel',use_hdf5=.false.)
         call viz%add_scalar(vof%VF,1,'VF')
         call viz%add_scalar(cclabel%id,1,'ID')
         call viz%add_surfmesh(vof%smesh,'plic')
      end block create_visualization

      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
     
      

      ! Write visualization with IDs
      call viz%write(time=0.0_WP)
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays
      call io%finalize()
      call amr%finalize()
      call vf%finalize()
      
   end subroutine simulation_final
   
end module simulation
