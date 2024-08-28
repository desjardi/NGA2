!> Incompressible scalar solver class built onto an amrcore
module amrscalar_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrconfig_class,  only: amrconfig,amrdata
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
      type(amrdata) :: SC,SCold          !< SC amrdata object
      
      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), pointer :: lo_bc                      !< Boundary condition descriptor in minus direction
      integer, dimension(:,:), pointer :: hi_bc                      !< Boundary condition descriptor in plus direction
      
      ! Interpolation method
      integer :: interp                                              !< Interpolation method

      ! Overlap size
      integer :: nover                                               !< Size of the overlap/ghost
      
   contains
      procedure :: initialize       !< Initialize scalar solver
      procedure :: finalize         !< Finalize scalar solver
      procedure :: get_dSCdt_lvl    !< Calculate dSC/dt at level (lvl)
      procedure :: get_info         !< Calculate various information on our amrscalar object
   end type amrscalar
   
   
contains
   
   
   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,nscalar,name)
      use messager,         only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_bc_reflect_even,amrex_interp_cell_cons
      use amrconfig_class,  only: resize,interp
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrconfig), target, intent(in) :: amr
      integer, intent(in) :: nscalar
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to amrconfig object
      this%amr=>amr
      
      ! Set the number of scalars
      this%nscalar=nscalar
      if (this%nscalar.le.0) call die('[amrscalar initialize] At least 1 scalar is required')
      
      ! Initialize scalar names
      allocate(this%SCname(1:this%nscalar))
      this%SCname='' ! User will set names
      
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

      ! Initialize amrdata variables
      call this%SC   %initialize(amr=this%amr,ncomp=this%nscalar,nover=this%nover,reg=interp,name='SC')
      call this%SCold%initialize(amr=this%amr,ncomp=this%nscalar,nover=this%nover,reg=resize,name='SCold')
      
      ! Transfer solver boundary conditions via pointer so the user can change later
      deallocate(this%SC   %lo_bc); this%SC   %lo_bc=>this%lo_bc; this%SC   %bc_ptr=.true.
      deallocate(this%SC   %hi_bc); this%SC   %hi_bc=>this%hi_bc; this%SC   %bc_ptr=.true.
      deallocate(this%SCold%lo_bc); this%SCold%lo_bc=>this%lo_bc; this%SCold%bc_ptr=.true.
      deallocate(this%SCold%hi_bc); this%SCold%hi_bc=>this%hi_bc; this%SCold%bc_ptr=.true.

   end subroutine initialize
   
   
   !> Finalization for amrscalar solver
   impure elemental subroutine finalize(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      call this%SC%finalize()
      call this%SCold%finalize()
      this%amr=>NULL()
      deallocate(this%SCname)
      if (associated(this%lo_bc)) deallocate(this%lo_bc); nullify(this%lo_bc)
      if (associated(this%hi_bc)) deallocate(this%hi_bc); nullify(this%hi_bc)
   end subroutine finalize
   
   
   !> Calculate dSC/dt at level (lvl)
   subroutine get_dSCdt_lvl(this,lvl,dSCdt,U,V,W)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box,amrex_fab,amrex_fab_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrdata), intent(inout) :: dSCdt
      type(amrdata), intent(in)    :: U,V,W
      type(amrex_mfiter)    :: mfi
      type(amrex_box)       :: bx,tbx
      type(amrex_fab)       :: xflux,yflux,zflux
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs,SC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: FX,FY,FZ
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      integer :: i,j,k,nsc
      
      ! Loop over boxes - no tiling for now
      call this%amr%mfiter_build(lvl,mfi)
      do while(mfi%next())
         bx=mfi%tilebox()
         
         ! Get SC and dSCdt data
         SC =>this%SC%data(lvl)%dataptr(mfi)
         rhs=>  dSCdt%data(lvl)%dataptr(mfi)
         
         ! Get velocity data
         pU=>U%data(lvl)%dataptr(mfi)
         pV=>V%data(lvl)%dataptr(mfi)
         pW=>W%data(lvl)%dataptr(mfi)
         
         ! Prepare face flux storage
         tbx=bx; call tbx%nodalize(1); call xflux%resize(tbx,1); FX=>xflux%dataptr()
         tbx=bx; call tbx%nodalize(2); call yflux%resize(tbx,1); FY=>yflux%dataptr()
         tbx=bx; call tbx%nodalize(3); call zflux%resize(tbx,1); FZ=>zflux%dataptr()
         
         ! Calculate fluxes for all components
         do nsc=1,this%nscalar
            ! Convective fluxes in x
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)+1
                     FX(i,j,k,1)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*(-1.0_WP/6.0_WP*SC(i-2,j,k,nsc)+5.0_WP/6.0_WP*SC(i-1,j,k,nsc)+2.0_WP/6.0_WP*SC(i  ,j,k,nsc)) &
                     &           -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*(+2.0_WP/6.0_WP*SC(i-1,j,k,nsc)+5.0_WP/6.0_WP*SC(i  ,j,k,nsc)-1.0_WP/6.0_WP*SC(i+1,j,k,nsc))
                  end do
               end do
            end do
            ! Convective fluxes in y
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)+1
                  do i=bx%lo(1),bx%hi(1)
                     FY(i,j,k,1)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*(-1.0_WP/6.0_WP*SC(i,j-2,k,nsc)+5.0_WP/6.0_WP*SC(i,j-1,k,nsc)+2.0_WP/6.0_WP*SC(i,j  ,k,nsc)) &
                     &           -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*(+2.0_WP/6.0_WP*SC(i,j-1,k,nsc)+5.0_WP/6.0_WP*SC(i,j  ,k,nsc)-1.0_WP/6.0_WP*SC(i,j+1,k,nsc))
                  end do
               end do
            end do
            ! Convective fluxes in z
            do k=bx%lo(3),bx%hi(3)+1
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     FZ(i,j,k,1)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*(-1.0_WP/6.0_WP*SC(i,j,k-2,nsc)+5.0_WP/6.0_WP*SC(i,j,k-1,nsc)+2.0_WP/6.0_WP*SC(i,j,k  ,nsc)) &
                     &           -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*(+2.0_WP/6.0_WP*SC(i,j,k-1,nsc)+5.0_WP/6.0_WP*SC(i,j,k  ,nsc)-1.0_WP/6.0_WP*SC(i,j,k+1,nsc))
                  end do
               end do
            end do
            ! Calculate rhs for all components
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     rhs(i,j,k,nsc)=(FX(i+1,j,k,1)-FX(i,j,k,1))/this%amr%dx(lvl)+&
                     &              (FY(i,j+1,k,1)-FY(i,j,k,1))/this%amr%dy(lvl)+&
                     &              (FZ(i,j,k+1,1)-FZ(i,j,k,1))/this%amr%dz(lvl)
                  end do
               end do
            end do
         end do
         
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Deallocate flux storage
      call amrex_fab_destroy(xflux)
      call amrex_fab_destroy(yflux)
      call amrex_fab_destroy(zflux)
      
   end subroutine get_dSCdt_lvl
   
   
   !> Calculate various information on our amrscalar object
   subroutine get_info(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      ! Generate info on SC object
      call this%SC%get_info()
   end subroutine get_info
   
   
end module amrscalar_class
