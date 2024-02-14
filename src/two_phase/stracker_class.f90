!> Two-phase structure tracking class
!> Augment cclabel object with scalar solver persistent indexing of structures
module stracker_class
   use precision,      only: WP
   use vfs_class,      only: vfs
   use cclabel_class,  only: cclabel
   use tpscalar_class, only: tpscalar
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: stracker
   
   
   !> Structure tracker object definition
   type, extends(cclabel) :: stracker
      ! Pointer to our VF solver
      class(vfs), pointer :: vf
      ! Include a two-phase scalar transport solver
      type(tpscalar) :: sc
      ! Pointer to make_label function
      !procedure(make_label_ftype), pointer :: make_label=>null()
   contains
      procedure :: init                                    !< Initialization of stracker class (different name is used because of extension...)
      procedure :: advance                                 !< Advance labels based on VOF transport
      !procedure, private :: same_label                     !< Checks if two points belong to same structure
   end type stracker
   

   !> Type of the make_label function used to generate a structure - user-provided
   !abstract interface
   !   logical function make_label_ftype(ind1,ind2,ind3)
   !      integer, intent(in) :: ind1,ind2,ind3
   !   end function make_label_ftype
   !end interface
   

contains
   
   
   !> Function that identifies if cell pairs have same label
   !> This implements merging and break-up of structures
   !logical function same_label(this,i1,j1,k1,i2,j2,k2)
   !   implicit none
   !   class(stracker), intent(in) :: this
   !   integer, intent(in) :: i1,j1,k1,i2,j2,k2
   !   if (this%sc%SC(i1,j1,k1,1).eq.this%sc%SC(i2,j2,k2,1)) then
   !      same_label=.true.
   !   else
   !      same_label=.false.
   !   end if
   !end function same_label
   
   
   !> Structure tracker initialization
   subroutine init(this,vf,phase,name)!,make_label)
      use messager, only: die
      implicit none
      class(stracker), intent(inout) :: this
      class(vfs), target, intent(in) :: vf
      integer, intent(in) :: phase
      !procedure(make_label_ftype) :: make_label
      character(len=*), optional :: name
      ! Point to our vfs object
      this%vf=>vf
      ! Check that vf stores detailed face flux
      if (.not.this%vf%store_detailed_flux) call die('[stracker init] vfs object needs to store detailed face fluxes')
      ! Create a cclabel object to get SIN
      call this%cclabel%initialize(pg=this%vf%cfg%pgrid)
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Create a scalar solver to transport LIN
      call this%sc%initialize(cfg=this%vf%cfg,nscalar=1)
      this%sc%phase=phase; this%sc%SCname(1)='LIN'
      ! Store pointer to make_label
      !this%make_label=>make_label
      ! Perform an initial CCL
      !call this%build(this%make_label,this%same_label)
   end subroutine init
   
   
   !> Advance structure tracker data by one time step
   subroutine advance(this,dt,U,V,W)
      implicit none
      class(stracker), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%vf%cfg%imino_:,this%vf%cfg%jmino_:,this%vf%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%vf%cfg%imino_:,this%vf%cfg%jmino_:,this%vf%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%vf%cfg%imino_:,this%vf%cfg%jmino_:,this%vf%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: p,q
      real(WP), dimension(:,:,:,:), allocatable :: resSC
      
      ! Advance LIN
      this%sc%SCold=this%sc%SC
      allocate(resSC(this%sc%cfg%imino_:this%sc%cfg%imaxo_,this%sc%cfg%jmino_:this%sc%cfg%jmaxo_,this%sc%cfg%kmino_:this%sc%cfg%kmaxo_,1:1))
      call this%sc%get_dSCdt(dSCdt=resSC,U=U,V=V,W=W,VFold=this%vf%VFold,VF=this%vf%VF,detailed_face_flux=this%vf%detailed_face_flux,dt=dt)
      p=real(this%sc%phase(1),WP); q=1.0_WP-2.0_WP*p
      where (this%sc%mask.eq.0.and.this%vf%VF.ne.p) this%sc%SC(:,:,:,1)=((p+q*this%vf%VFold)*this%sc%SCold(:,:,:,1)+dt*resSC(:,:,:,1))/(p+q*this%vf%VF)
      where (this%vf%VF.eq.p) this%sc%SC(:,:,:,1)=0.0_WP
      call this%sc%cfg%sync(this%sc%SC(:,:,:,1))
      deallocate(resSC)
      
      ! Perform CCL
      !call this%build(this%make_label,this%same_label)
      
   end subroutine advance


   
   
end module stracker_class