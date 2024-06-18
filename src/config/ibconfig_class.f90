!> Single-grid IB-capable config concept is defined here: this is a partitioned grid
!> as well as geometry based on VF augmented with a level set and its normal and a few methods
module ibconfig_class
   use precision,    only: WP
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ibconfig
   
   ! List of known available methods for calculating VF from G
   integer, parameter, public :: bigot=1
   integer, parameter, public :: sharp=2
   
   !> Config object definition as an extension of config
   type, extends(config) :: ibconfig
      
      real(WP), dimension(:,:,:),   allocatable :: Gib             !< Level set function, negative in fluid and positive in solid
      real(WP), dimension(:,:,:,:), allocatable :: Nib             !< IB normal vector, oriented into the solid
      
   contains
      procedure :: calculate_normal          !< Calculate robustly the normal vector from G field
      procedure :: calculate_vf              !< Calculate fluid volume fraction
   end type ibconfig
   
   
   !> Declare single-grid ibconfig constructor
   interface ibconfig
      procedure construct_from_sgrid
   end interface ibconfig
   
   
contains
   
   
   !> Single-grid config constructor from a serial grid
   function construct_from_sgrid(grp,decomp,grid) result(self)
      use sgrid_class, only: sgrid
      use string,      only: str_medium
      use mpi_f08,     only: MPI_Group
      implicit none
      type(ibconfig) :: self
      type(sgrid), intent(in) :: grid
      type(MPI_Group), intent(in) :: grp
      integer, dimension(3), intent(in) :: decomp
      ! Create a config with the provided group and decomposition
      self%config=config(grp,decomp,grid)
      ! Allocate IB info
      allocate(self%Gib(    self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_)); self%Gib=-1.0_WP
      allocate(self%Nib(1:3,self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_)); self%Nib=+0.0_WP
   end function construct_from_sgrid
   
   
   !> Specialized normal vector calculation from G
   !> Well suited to handle G discontinuities
   subroutine calculate_normal(this)
      implicit none
      class(ibconfig), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,stx,sty,stz
      real(WP), dimension(3) :: Gx,Gy,Gz
      real(WP) :: norm,norm_test
      ! Get the normal vector
      do k=this%kmino_+1,this%kmaxo_-1
         do j=this%jmino_+1,this%jmaxo_-1
            do i=this%imino_+1,this%imaxo_-1
               ! Consider 3 stencils in x
               Gx(1)=(this%Gib(i+1,j,k)-this%Gib(i-1,j,k))/(this%xm(i+1)-this%xm(i-1))
               Gx(2)=(this%Gib(i  ,j,k)-this%Gib(i-1,j,k))/(this%xm(i  )-this%xm(i-1))
               Gx(3)=(this%Gib(i+1,j,k)-this%Gib(i  ,j,k))/(this%xm(i+1)-this%xm(i  ))
               ! Consider 3 stencils in y
               Gy(1)=(this%Gib(i,j+1,k)-this%Gib(i,j-1,k))/(this%ym(j+1)-this%ym(j-1))
               Gy(2)=(this%Gib(i,j  ,k)-this%Gib(i,j-1,k))/(this%ym(j  )-this%ym(j-1))
               Gy(3)=(this%Gib(i,j+1,k)-this%Gib(i,j  ,k))/(this%ym(j+1)-this%ym(j  ))
               ! Consider 3 stencils in z
               Gz(1)=(this%Gib(i,j,k+1)-this%Gib(i,j,k-1))/(this%zm(k+1)-this%zm(k-1))
               Gz(2)=(this%Gib(i,j,k  )-this%Gib(i,j,k-1))/(this%zm(k  )-this%zm(k-1))
               Gz(3)=(this%Gib(i,j,k+1)-this%Gib(i,j,k  ))/(this%zm(k+1)-this%zm(k  ))
               ! Find stencil combination closest to unity norm
               norm=huge(1.0_WP)
               do kk=1,3
                  do jj=1,3
                     do ii=1,3
                        norm_test=sqrt(Gx(ii)**2+Gy(jj)**2+Gz(kk)**2)
                        if (abs(norm_test-1.0_WP).lt.abs(norm-1.0_WP)) then
                           norm=norm_test; stx=ii; sty=jj; stz=kk
                        end if
                     end do
                  end do
               end do
               norm=max(norm,epsilon(1.0_WP))
               this%Nib(1:3,i,j,k)=[Gx(stx),Gy(sty),Gz(stz)]/norm
            end do
         end do
      end do
      ! Synchronize it
      call this%sync(this%Nib)
      ! Extend to non-periodic edges
      if (.not.this%xper) then
         if (this%iproc.eq.1) then
            this%Nib(:,this%imino,:,:)=this%Nib(:,this%imino+1,:,:)
         else if (this%iproc.eq.this%npx) then
            this%Nib(:,this%imaxo,:,:)=this%Nib(:,this%imaxo-1,:,:)
         end if
      end if
      if (.not.this%yper) then
         if (this%jproc.eq.1) then
            this%Nib(:,:,this%jmino,:)=this%Nib(:,:,this%jmino+1,:)
         else if (this%jproc.eq.this%npy) then
            this%Nib(:,:,this%jmaxo,:)=this%Nib(:,:,this%jmaxo-1,:)
         end if
      end if
      if (.not.this%zper) then
         if (this%kproc.eq.1) then
            this%Nib(:,:,:,this%kmino)=this%Nib(:,:,:,this%kmino+1)
         else if (this%kproc.eq.this%npz) then
            this%Nib(:,:,:,this%kmaxo)=this%Nib(:,:,:,this%kmaxo-1)
         end if
      end if
   end subroutine calculate_normal
   
   
   !> Calculation of VF from G
   subroutine calculate_vf(this,method,allow_zero_vf)
      use messager, only: die
      implicit none
      class(ibconfig), intent(inout) :: this
      integer, intent(in) :: method
      logical, intent(in), optional :: allow_zero_vf
      
      ! Select the method
      select case (method)
      case(bigot)
         ! Get fluid volume fraction from G using Bigot's expression
         bigot: block
            integer :: i,j,k
            real(WP) :: lambda,eta
            do k=this%kmino_,this%kmaxo_
               do j=this%jmino_,this%jmaxo_
                  do i=this%imino_,this%imaxo_
                     lambda=sum(abs(this%Nib(:,i,j,k)))
                     eta=0.065_WP*(1.0_WP-lambda**2)+0.39_WP
                     this%VF(i,j,k)=0.5_WP*(1.0_WP-tanh(this%Gib(i,j,k)/(lambda*eta*this%meshsize(i,j,k))))
                  end do
               end do
            end do
         end block bigot
      case (sharp)
         ! Get fluid volume fraction from G using marching tets
         sharp: block
            use mms_geom, only: marching_tets
            integer :: i,j,k,n,si,sj,sk
            real(WP), dimension(3,8) :: hex_vertex
            real(WP), dimension(8) :: G
            real(WP), dimension(3) :: v_cent,a_cent
            real(WP) :: vol,area
            do k=this%kmino_,this%kmaxo_
               do j=this%jmino_,this%jmaxo_
                  do i=this%imino_,this%imaxo_
                     ! Form an hexahedron with corresponding distance info
                     n=0
                     do sk=0,1
                        do sj=0,1
                           do si=0,1
                              n=n+1
                              hex_vertex(:,n)=[this%x(i+si),this%y(j+sj),this%z(k+sk)]
                              G(n)=-this%get_scalar(hex_vertex(:,n),i,j,k,this%Gib,'d')
                           end do
                        end do
                     end do
                     ! Rescale by min_meshsize and shift
                     !do n=1,8
                     !   hex_vertex(:,n)=(hex_vertex(:,n)-[this%x(i),this%y(j),this%z(k)])/this%min_meshsize
                     !   G(n)=G(n)/this%min_meshsize
                     !end do
                     ! Use marching tets to get volume fraction
                     call marching_tets(hex_vertex,G,vol,area,v_cent,a_cent)
                     this%VF(i,j,k)=vol/this%vol(i,j,k)!*(this%min_meshsize)**3
                  end do
               end do
            end do
            call this%sync(this%VF) !< Sync needed because of the Gib interpolation above
         end block sharp
      case default
         call die('[ibconfig calculate_vf] Unknown method to calculate VF')
      end select
      
      ! Check if VF=0 is allowed
      if (present(allow_zero_vf)) then
         if (.not.allow_zero_vf) this%VF=max(this%VF,epsilon(1.0_WP))
      end if
      
      ! Update total fluid volume
      call this%calc_fluid_vol()

   end subroutine calculate_vf


end module ibconfig_class
