!> Incompressible scalar solver class built onto an amrconfig
module amrscalar_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrconfig_class,  only: amrconfig
   use amrex_amr_module, only: amrex_multifab,amrex_fluxregister
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrscalar
   
   !> Constant density scalar solver object definition
   type :: amrscalar
      
      ! This is our amrconfig
      class(amrconfig), pointer :: amr                               !< This is the amrconfig the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_AMRSCALAR'          !< Solver name (default=UNNAMED_AMRSCALAR)
      
      ! Scalar variable definition
      integer :: nscalar                                             !< Number of scalars
      character(len=str_medium), dimension(:), allocatable :: SCname !< Names of scalars
      
      ! Scalar and old scalar data
      type(amrex_multifab), dimension(:), allocatable :: SC          !< SC multifab array
      type(amrex_multifab), dimension(:), allocatable :: SCold       !< Old SC multifab array
      
      ! Flux register for interlevel conservation
      type(amrex_fluxregister), dimension(:), allocatable :: Freg    !< Flux register

      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), allocatable :: lo_bc                  !< Boundary condition descriptor in minus direction
      integer, dimension(:,:), allocatable :: hi_bc                  !< Boundary condition descriptor in plus direction
      
      ! Interpolation method
      integer :: interp                                              !< Interpolation method

      ! Overlap size
      integer :: nover                                               !< Size of the overlap/ghost

      ! Monitoring quantities
      real(WP), dimension(:), allocatable :: SCmax,SCmin,SCint       !< Maximum and minimum, integral scalar
      
   contains
      ! Basic procedures -------------------------------------------------------------------------------------
      procedure :: initialize       !< Initialize scalar solver
      procedure :: finalize         !< Finalize scalar solver
      procedure :: get_info         !< Calculate various information on our amrscalar object
      !procedure :: print            !< Output solver info to the screen
      ! Regrid procedures ------------------------------------------------------------------------------------
      procedure :: delete           !< Delete data at level (lvl)
      procedure :: create           !< Create data at level (lvl) and leave it uninitialized
      procedure :: refine           !< Refine data at level (lvl) using cfill procedure
      procedure :: remake           !< Remake data at level (lvl) using  fill procedure
      ! Fill procedures --------------------------------------------------------------------------------------
      procedure ::  cfill           !< Fill provided mfab at level (lvl) from this%SC at level (lvl-1)           - involves boundary conditions - done at a single time using this%SC
      procedure ::   fill           !< Fill provided mfab at level (lvl) from this%SC at level (lvl-1) and (lvl) - involves boundary conditions - done at two times using this%SC and this%SCold
      ! Advance procedures -----------------------------------------------------------------------------------
      procedure :: get_dSCdt        !< Calculate dSC/dt at level (lvl)
      procedure :: copy2old         !< Copy SC in SCold
      procedure :: reflux_avg_lvl   !< Perform refluxing and averaging at a given level
      procedure :: reflux_avg       !< Perform successive refluxing and averaging at all levels
   end type amrscalar
   
   
contains
   
   
   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,nscalar,name)
      use messager,         only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_bc_reflect_even,amrex_interp_cell_cons
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrconfig), target, intent(in) :: amr
      integer, intent(in) :: nscalar
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to amrconfig object
      this%amr=>amr
      
      ! Allocate variables
      allocate(this%SC   (0:this%amr%nlvl))
      allocate(this%SCold(0:this%amr%nlvl))
      allocate(this%Freg (0:this%amr%nlvl))
      
      ! Set the number of scalars
      this%nscalar=nscalar
      if (this%nscalar.le.0) call die('[amrscalar initialize] At least 1 scalar is required')
      
      ! Initialize scalar names
      allocate(this%SCname(1:this%nscalar))
      this%SCname='' ! User will set names
      
      ! Initialize info storage
      allocate(this%SCmin(1:this%nscalar)); this%SCmin=+huge(1.0_WP)
      allocate(this%SCmax(1:this%nscalar)); this%SCmax=-huge(1.0_WP)
      allocate(this%SCint(1:this%nscalar)); this%SCint= 0.0_WP
      
      ! Allocate boundary condition descriptor
      allocate(this%lo_bc(1:3,1:this%nscalar))
      allocate(this%hi_bc(1:3,1:this%nscalar))
      
      ! Initialize bc based on periodicity info - user can change later
      this%lo_bc=amrex_bc_reflect_even
      if (this%amr%xper) this%lo_bc(1,:)=amrex_bc_int_dir
      if (this%amr%yper) this%lo_bc(2,:)=amrex_bc_int_dir
      if (this%amr%zper) this%lo_bc(3,:)=amrex_bc_int_dir
      this%hi_bc=this%lo_bc
      
      ! Set overlap size for QUICK scheme
      this%nover=2
      
      ! Assume conservative interpolation - user can change later
      this%interp=amrex_interp_cell_cons
      
   end subroutine initialize
   
   
   !> Finalization for amrscalar solver
   impure elemental subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy,amrex_fluxregister_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl
      do lvl=0,this%amr%nlvl
         call amrex_multifab_destroy(this%SC   (lvl))
         call amrex_multifab_destroy(this%SCold(lvl))
         call amrex_fluxregister_destroy(this%Freg(lvl))
      end do
      deallocate(this%SC,this%SCold,this%Freg)
      deallocate(this%SCmax,this%SCmin,this%SCint,this%lo_bc,this%hi_bc,this%SCname)
      nullify(this%amr)
   end subroutine finalize
   
   
   !> Calculate various information on our amrscalar object
   subroutine get_info(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl,nsc
      
      ! Reset info
      this%SCmin=+huge(1.0_WP)
      this%SCmax=-huge(1.0_WP)
      this%SCint= 0.0_WP
      
      ! Loop over scalars
      do nsc=1,this%nscalar
         ! Loop over all levels
         do lvl=0,this%amr%clvl()
            ! Get min and max at that level
            this%SCmin(nsc)=min(this%SCmin(nsc),this%SC(lvl)%min(comp=nsc))
            this%SCmax(nsc)=max(this%SCmax(nsc),this%SC(lvl)%max(comp=nsc))
         end do
         ! Get int at level 0
         this%SCint(nsc)=this%SC(0)%sum(comp=nsc)*(this%amr%dx(0)*this%amr%dy(0)*this%amr%dz(0))/this%amr%vol
      end do
      
   end subroutine get_info
   
   
   !> Delete solver data at level lvl
   subroutine delete(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy,amrex_fluxregister_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      call amrex_multifab_destroy(this%SC   (lvl))
      call amrex_multifab_destroy(this%SCold(lvl))
      call amrex_fluxregister_destroy(this%Freg(lvl))
   end subroutine delete
   
   
   !> Create solver data at level lvl
   subroutine create(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap,amrex_fluxregister_build
      implicit none
      class(amrscalar), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Delete level
      call this%delete(lvl)
      ! Rebuild level
      call amrex_multifab_build(this%SC   (lvl),ba,dm,this%nscalar,0)
      call amrex_multifab_build(this%SCold(lvl),ba,dm,this%nscalar,0)
      ! Create flux register
      if (lvl.gt.0) call amrex_fluxregister_build(this%Freg(lvl),ba,dm,this%amr%rref(lvl-1),lvl,this%nscalar)
   end subroutine create
   
   
   !> Refine solver data at level lvl
   subroutine refine(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrscalar), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Recreate level
      call this%create(lvl,time,ba,dm)
      ! Populate SC from coarse level
      call this%cfill(lvl,time,this%SC(lvl))
   end subroutine refine
   
   
   !> Remake solver data at level lvl
   subroutine remake(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
      implicit none
      class(amrscalar), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab) :: SCnew
      ! Create SCnew and populate it from current data
      call amrex_multifab_build(SCnew,ba,dm,this%nscalar,0)
      call this%fill(lvl,time,SCnew)
      ! Recreate level
      call this%create(lvl,time,ba,dm)
      ! Copy SCnew to SC
      call this%SC(lvl)%copy(SCnew,1,1,this%nscalar,0)
      ! Destroy SCnew
      call amrex_multifab_destroy(SCnew)
   end subroutine remake
   
   
   !> Fill provided mfab at level (lvl) from this%SC at level (lvl-1)
   subroutine cfill(this,lvl,time,SC)
      use amrex_amr_module, only: amrex_multifab,amrex_fillcoarsepatch
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: SC
      call amrex_fillcoarsepatch(          SC,&  !< fine data being filled...
      &                   time,this%SC(lvl-1),&  !< using coarse data at old time...
      &                   time,this%SC(lvl-1),&  !<   and coarse data at new time...
      &           this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
      &           this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
      &                 time,1,1,this%nscalar,&  !< time when we want the data, scomp, dcomp, ncomp...
      &                  this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
      &                           this%interp,&  !< interpolation strategy...
      &                 this%lo_bc,this%hi_bc)   !< domain bconds
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc,this%hi_bc)
            end if
         end do
         ! This will need hooks for user-provided BCs
      end subroutine fillbc
   end subroutine cfill
   
   
   !> Fill provided mfab at level (lvl) from this%SC at level (lvl-1) and (lvl)
   subroutine fill(this,lvl,time,SC,tnew,told)
      use messager,         only: die
      use amrex_amr_module, only: amrex_multifab,amrex_fillpatch
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: SC
      real(WP), intent(in), optional :: tnew,told
      real(WP) :: new_time,old_time
      ! Handle time info
      if (present(told).and.present(tnew)) then
         ! Use both SC and SCold to populate mfab
         new_time=tnew
         old_time=told
      else if (.not.present(told).and..not.present(tnew)) then
         ! Only use SC to populate mfab
         new_time=time
         old_time=time-1.0e200_WP
      else
         ! One is provided but not the other...
         call die('[amrscalar fill] tnew and told must be either both provided or both ignored')
      end if
      ! Fill mfab
      if (lvl.eq.0) then
         ! Fill without interpolation, just direct copy and bconds
         call amrex_fillpatch(             SC,&  !< base data being filled...
         &           old_time,this%SCold(lvl),&  !< using base data at old time...
         &           new_time,this%SC   (lvl),&  !<   and base data at new time...
         &          this%amr%geom(lvl),fillbc,&  !< base geometry with function to apply bconds...
         &              time,1,1,this%nscalar)   !< time when we want the data, scomp, dcomp, ncomp
         ! Unclear why lo_bc and hi_bc aren't involved here...
      else
         ! Fill with a mix of interpolation, direct copy and bconds
         call amrex_fillpatch(             SC,&  !< fine data being filled...
         &         old_time,this%SCold(lvl-1),&  !< using coarse data at old time...
         &         new_time,this%SC   (lvl-1),&  !<   and coarse data at new time...
         &        this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
         &         old_time,this%SCold(lvl  ),&  !<     and fine data at old time...
         &         new_time,this%SC   (lvl  ),&  !<     and fine data at new time...
         &        this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
         &              time,1,1,this%nscalar,&  !< time when we want the data, scomp, dcomp, ncomp...
         &               this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
         &                        this%interp,&  !< interpolation strategy...
         &              this%lo_bc,this%hi_bc)   !< domain bconds
      end if
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc,this%hi_bc)
            end if
         end do
         ! This will need hooks for user-provided BCs
      end subroutine fillbc
   end subroutine fill
   
   
   !> Calculate dSC/dt at level (lvl)
   !! dSCdt needs to be ncomp=nscalar,nover=0    ,atface=[.false.,.false.,.false.]
   !! SC    needs to be ncomp=nscalar,nover=nover,atface=[.false.,.false.,.false.]
   !! U     needs to be ncomp=1      ,nover=0    ,atface=[.true. ,.false.,.false.]
   !! V     needs to be ncomp=1      ,nover=0    ,atface=[.false.,.true. ,.false.]
   !! W     needs to be ncomp=1      ,nover=0    ,atface=[.false.,.false.,.true. ]
   subroutine get_dSCdt(this,lvl,dSCdt,SC,U,V,W)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_fab,amrex_fab_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: dSCdt
      type(amrex_multifab), intent(in) :: SC,U,V,W
      type(amrex_mfiter)    :: mfi
      type(amrex_box)       :: bx,tbx
      type(amrex_multifab), dimension(3) :: flux
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs,pSC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: FX,FY,FZ
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      integer :: i,j,k,nsc
      
      ! Prepare flux multifabs
      call this%amr%mfab_build(lvl=lvl,mfab=flux(1),ncomp=this%nscalar,nover=0,atface=[.true. ,.false.,.false.])
      call this%amr%mfab_build(lvl=lvl,mfab=flux(2),ncomp=this%nscalar,nover=0,atface=[.false.,.true. ,.false.])
      call this%amr%mfab_build(lvl=lvl,mfab=flux(3),ncomp=this%nscalar,nover=0,atface=[.false.,.false.,.true. ])
      
      ! Loop over boxes - no tiling for now
      call this%amr%mfiter_build(lvl,mfi)
      do while(mfi%next())
         bx=mfi%tilebox()
         
         ! Get SC and dSCdt data
         pSC=>    SC%dataptr(mfi)
         rhs=> dSCdt%dataptr(mfi)
         
         ! Get velocity data
         pU=>U%dataptr(mfi)
         pV=>V%dataptr(mfi)
         pW=>W%dataptr(mfi)
         
         ! Get face flux
         FX=>flux(1)%dataptr(mfi)
         FY=>flux(2)%dataptr(mfi)
         FZ=>flux(3)%dataptr(mfi)
         
         ! Calculate fluxes for all components
         do nsc=1,this%nscalar
            ! Convective fluxes in x
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)+1
                     FX(i,j,k,nsc)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i-2,j,k,nsc)+5.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+2.0_WP/6.0_WP*pSC(i  ,j,k,nsc)) &
                     &             -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+5.0_WP/6.0_WP*pSC(i  ,j,k,nsc)-1.0_WP/6.0_WP*pSC(i+1,j,k,nsc))
                  end do
               end do
            end do
            ! Convective fluxes in y
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)+1
                  do i=bx%lo(1),bx%hi(1)
                     FY(i,j,k,nsc)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j-2,k,nsc)+5.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+2.0_WP/6.0_WP*pSC(i,j  ,k,nsc)) &
                     &             -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+5.0_WP/6.0_WP*pSC(i,j  ,k,nsc)-1.0_WP/6.0_WP*pSC(i,j+1,k,nsc))
                  end do
               end do
            end do
            ! Convective fluxes in z
            do k=bx%lo(3),bx%hi(3)+1
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     FZ(i,j,k,nsc)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j,k-2,nsc)+5.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+2.0_WP/6.0_WP*pSC(i,j,k  ,nsc)) &
                     &             -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+5.0_WP/6.0_WP*pSC(i,j,k  ,nsc)-1.0_WP/6.0_WP*pSC(i,j,k+1,nsc))
                  end do
               end do
            end do
            ! Calculate rhs for all components
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     rhs(i,j,k,nsc)=(FX(i+1,j,k,nsc)-FX(i,j,k,nsc))/this%amr%dx(lvl)+&
                     &              (FY(i,j+1,k,nsc)-FY(i,j,k,nsc))/this%amr%dy(lvl)+&
                     &              (FZ(i,j,k+1,nsc)-FZ(i,j,k,nsc))/this%amr%dz(lvl)
                  end do
               end do
            end do
         end do
         ! Rescale fluxes by face area for later refluxing
         FX=FX*this%amr%dy(lvl)*this%amr%dz(lvl)
         FY=FY*this%amr%dz(lvl)*this%amr%dx(lvl)
         FZ=FZ*this%amr%dx(lvl)*this%amr%dy(lvl)
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Handle refluxing
      if (lvl.gt.0)               call this%Freg(lvl  )%fineadd (flux,-1.0_WP)
      if (lvl.lt.this%amr%clvl()) call this%Freg(lvl+1)%crseinit(flux,+1.0_WP)
      
      ! Deallocate flux storage
      call this%amr%mfab_destroy(flux(1))
      call this%amr%mfab_destroy(flux(2))
      call this%amr%mfab_destroy(flux(3))
      
   end subroutine get_dSCdt
   
   
   !> Copy SC in SCold
   subroutine copy2old(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl
      ! Loop over all levels
      do lvl=0,this%amr%clvl()
         call this%SCold(lvl)%copy(srcmf=this%SC(lvl),srccomp=1,dstcomp=1,nc=this%nscalar,ng=0)
      end do
   end subroutine copy2old
   
   
   !> Perform successive refluxing and averaging at a level
   subroutine reflux_avg_lvl(this,lvl,dt)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl
      interface
         subroutine amrex_fi_fluxregister_reflux (fr,mf,scale,geom) bind(c)
            import :: c_ptr,WP
            implicit none
            type(c_ptr), value :: fr,mf,geom
            real(WP), value :: scale
         end subroutine amrex_fi_fluxregister_reflux
      end interface
      call amrex_fi_fluxregister_reflux(this%Freg(lvl+1)%p,this%SC(lvl)%p,dt,this%amr%geom(lvl)%p)
      call this%amr%average_downto(this%SC,lvl)
   end subroutine reflux_avg_lvl
   
   
   !> Perform successive refluxing and averaging at all levels
   subroutine reflux_avg(this,dt)
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl
      ! Loop over all levels except last one
      do lvl=this%amr%clvl()-1,0,-1
         call this%reflux_avg_lvl(lvl,dt)
      end do
   end subroutine reflux_avg
   
   
end module amrscalar_class
