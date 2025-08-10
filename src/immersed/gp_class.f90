!> Ghost point immersed boundary class
!> Provides support for Dirichle, Neumann, and Robin boundary conditions
module gp_class
  use precision,      only: WP
  use ibconfig_class, only: ibconfig
  implicit none
  private

  ! Expose type/constructor/methods
  public :: gpibm

  ! List of known available bcond types for this solver
  integer, parameter, public :: dirichlet=1         !< Dirichlet condition
  integer, parameter, public :: neumann=2           !< Normal gradient
  integer, parameter, public :: robin=3             !< Robin condition

  !> Basic image point definition
  type :: image
     integer, dimension(3,2) :: st                       !< Interpolation stencil extent
     integer , dimension(3) :: ind                       !< Index of cell containing image point
     real(WP), dimension(3) :: pos                       !< Coordinates of image point
     real(WP), dimension(2,2,2) :: itp                   !< Interpolation weights
  end type image

  !> Basic ghost point definition
  type :: ghost
     integer, dimension(3) :: ind                        !< Index of cell containing ghost point
     real(WP) :: dist                                    !< Distance to surface
     type(image) :: im                                   !< Associated image point
  end type ghost

  !> SGS model object definition
  type :: gpibm

     ! This is our config
     class(ibconfig), pointer :: cfg                     !< This is the config the model is build for

     ! These are the ghost points
     integer :: no                                       !< Number of overlapping ghost points
     integer :: ngp,ngpx,ngpy,ngpz                       !< Number of ghost points for our solver
     type(ghost), dimension(:), allocatable :: gp        !< Array of ghost points at cell center
     type(ghost), dimension(:), allocatable :: gpx       !< Array of ghost points at X-face
     type(ghost), dimension(:), allocatable :: gpy       !< Array of ghost points at Y-face
     type(ghost), dimension(:), allocatable :: gpz       !< Array of ghost points at Z-face
     integer, dimension(:,:,:), allocatable :: label     !< Integer array used for labeling ghost/image points

   contains

     procedure :: update                                 !< Update ghost point information
     procedure :: apply_bcond                            !< Apply specified boundary conditions
     procedure :: apply_robin                            !< Apply Robin boundary condition
     !generic :: apply_bcond=>apply_bcond_scalar,apply_bcond_array   !< Apply specified boundary conditions
     !procedure, private :: apply_bcond_scalar,apply_bcond_array

  end type gpibm


  !> Declare model constructor
  interface gpibm
     procedure constructor
  end interface gpibm

contains


  !> Default constructor for model
  function constructor(cfg,no) result(self)
    implicit none
    type(gpibm) :: self
    class(ibconfig), target, intent(in) :: cfg
    integer, intent(in) :: no

    ! Point to config object
    self%cfg=>cfg

    ! Reset ghost point list
    self%no=no
    self%ngp=0; self%ngpx=0; self%ngpy=0; self%ngpz=0
    if (allocated(self%gp)) deallocate(self%gp)
    if (allocated(self%gpx)) deallocate(self%gpx)
    if (allocated(self%gpy)) deallocate(self%gpy)
    if (allocated(self%gpz)) deallocate(self%gpz)

    ! Allocate label array (0=fluid cell, +1=ghost point, -1=image point)
    allocate(self%label(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%label=0

  end function constructor


  !> Updates list of ghost points and associated image point data
  subroutine update(this)
    use messager, only: die
    use param, only: verbose
    implicit none
    class(gpibm), intent(inout) :: this
    integer :: i,j,k,ii,jj,kk,i1,i2,j1,j2,k1,k2,n
    real(WP) :: buf
    real(WP), dimension(3) :: pos,pos_im
    real(WP), dimension(2,2,2) :: alpha3D,delta3D,dist3D,eta3D
    real(WP), parameter :: eps=1.0e-10_WP
    real(WP), dimension(:,:,:), allocatable :: Gx,Gy,Gz
    real(WP), dimension(:,:,:,:), allocatable :: Nx,Ny,Nz
    real(WP), dimension(:,:,:,:), allocatable :: itpr_x,itpr_y,itpr_z
    logical :: success

    ! Reset ghost points
    this%ngp=0; this%ngpx=0; this%ngpy=0; this%ngpz=0
    if (allocated(this%gp)) deallocate(this%gp)
    if (allocated(this%gpx)) deallocate(this%gpx)
    if (allocated(this%gpy)) deallocate(this%gpy)
    if (allocated(this%gpz)) deallocate(this%gpz)

    ! Create temporary metrics
    allocate(itpr_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
    allocate(itpr_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
    allocate(itpr_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
    ! Create density (or other things) interpolation coefficients to cell face
    do k=this%cfg%kmin_,this%cfg%kmax_+1
       do j=this%cfg%jmin_,this%cfg%jmax_+1
          do i=this%cfg%imin_,this%cfg%imax_+1
             itpr_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
             itpr_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
             itpr_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
          end do
       end do
    end do

    ! Store Gib and Nib at cell faces
    allocate(Gx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(Gy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(Gz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(Nx(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(Ny(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(Nz(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             Gx(i,j,k)=sum(itpr_x(:,i,j,k)*this%cfg%Gib(i-1:i,j,k))
             Gy(i,j,k)=sum(itpr_y(:,i,j,k)*this%cfg%Gib(i,j-1:j,k))
             Gz(i,j,k)=sum(itpr_z(:,i,j,k)*this%cfg%Gib(i,j,k-1:k))
             Nx(1,i,j,k)=sum(itpr_x(:,i,j,k)*this%cfg%Nib(1,i-1:i,j,k))
             Nx(2,i,j,k)=sum(itpr_x(:,i,j,k)*this%cfg%Nib(2,i-1:i,j,k))
             Nx(3,i,j,k)=sum(itpr_x(:,i,j,k)*this%cfg%Nib(3,i-1:i,j,k))
             Ny(1,i,j,k)=sum(itpr_y(:,i,j,k)*this%cfg%Nib(1,i,j-1:j,k))
             Ny(2,i,j,k)=sum(itpr_y(:,i,j,k)*this%cfg%Nib(2,i,j-1:j,k))
             Ny(3,i,j,k)=sum(itpr_y(:,i,j,k)*this%cfg%Nib(3,i,j-1:j,k))
             Nz(1,i,j,k)=sum(itpr_z(:,i,j,k)*this%cfg%Nib(1,i,j,k-1:k))
             Nz(2,i,j,k)=sum(itpr_z(:,i,j,k)*this%cfg%Nib(2,i,j,k-1:k))
             Nz(3,i,j,k)=sum(itpr_z(:,i,j,k)*this%cfg%Nib(3,i,j,k-1:k))
          end do
       end do
    end do
    call this%cfg%sync(Gx)
    call this%cfg%sync(Gy)
    call this%cfg%sync(Gz)
    call this%cfg%sync(Nx)
    call this%cfg%sync(Ny)
    call this%cfg%sync(Nz)

    ! Cell center
    !========================================================================================
    ! Identify ghost points
    this%ngp=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (this%cfg%Gib(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (this%cfg%Gib(ii,jj,kk).gt.0.0_WP.and.this%cfg%VF(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) this%ngp=this%ngp+1
          end do
       end do
    end do
    if (this%ngp.gt.0) allocate(this%gp(this%ngp))
    ! Store ghost points and associated image points
    n=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (this%cfg%Gib(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (this%cfg%Gib(ii,jj,kk).gt.0.0_WP.and.this%cfg%VF(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) then
                n=n+1
                ! Store index of ghost points
                this%gp(n)%ind(1)=i; this%gp(n)%ind(2)=j; this%gp(n)%ind(3)=k
                pos=(/this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)/)
                ! Get position and index of image point
                this%gp(n)%dist=abs(this%cfg%Gib(i,j,k))
                pos_im=pos+2.0_WP*this%gp(n)%dist*this%cfg%Nib(:,i,j,k)
                ! Find right i index
                i1=i
                do while (pos_im(1)-this%cfg%xm(i1  ).lt.0.0_WP.and.i1  .gt.this%cfg%imino_); i1=i1-1; end do
                do while (pos_im(1)-this%cfg%xm(i1+1).ge.0.0_WP.and.i1+1.lt.this%cfg%imaxo_); i1=i1+1; end do
                i2=i1+1     
                ! Find right j index
                j1=j
                do while (pos_im(2)-this%cfg%ym(j1  ).lt.0.0_WP.and.j1  .gt.this%cfg%jmino_); j1=j1-1; end do
                do while (pos_im(2)-this%cfg%ym(j1+1).ge.0.0_WP.and.j1+1.lt.this%cfg%jmaxo_); j1=j1+1; end do
                j2=j1+1      
                ! Find right k index
                k1=k
                do while (pos_im(3)-this%cfg%zm(k1  ).lt.0.0_WP.and.k1  .gt.this%cfg%kmino_); k1=k1-1; end do
                do while (pos_im(3)-this%cfg%zm(k1+1).ge.0.0_WP.and.k1+1.lt.this%cfg%kmaxo_); k1=k1+1; end do
                k2=k1+1      
                ! Check if the image point is inside the levelset
                alpha3D = 1.0_WP
                if (this%cfg%Gib(i1,j1,k1).le.0.0_WP) alpha3D(1,1,1) = 0.0_WP
                if (this%cfg%Gib(i2,j1,k1).le.0.0_WP) alpha3D(2,1,1) = 0.0_WP
                if (this%cfg%Gib(i1,j2,k1).le.0.0_WP) alpha3D(1,2,1) = 0.0_WP
                if (this%cfg%Gib(i2,j2,k1).le.0.0_WP) alpha3D(2,2,1) = 0.0_WP
                if (this%cfg%Gib(i1,j1,k2).le.0.0_WP) alpha3D(1,1,2) = 0.0_WP
                if (this%cfg%Gib(i2,j1,k2).le.0.0_WP) alpha3D(2,1,2) = 0.0_WP
                if (this%cfg%Gib(i1,j2,k2).le.0.0_WP) alpha3D(1,2,2) = 0.0_WP
                if (this%cfg%Gib(i2,j2,k2).le.0.0_WP) alpha3D(2,2,2) = 0.0_WP
                buf = sum(alpha3D)
                if (buf.gt.0.0_WP) then
                   ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
                   dist3D(1,1,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,1,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,2,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,2,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,1,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,1,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(1,2,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,2,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   delta3D=0.0_WP
                   if (dist3D(1,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,1)=1.0_WP
                   else if (dist3D(2,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,1)=1.0_WP
                   else if (dist3D(1,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,1)=1.0_WP
                   else if (dist3D(2,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,1)=1.0_WP
                   else if (dist3D(1,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,2)=1.0_WP
                   else if (dist3D(2,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,2)=1.0_WP
                   else if (dist3D(1,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,2)=1.0_WP
                   else if (dist3D(2,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,2)=1.0_WP
                   else
                      eta3D=1.0_WP/dist3D**2
                      buf=sum(eta3D*alpha3D)
                      delta3D=alpha3D*eta3D/buf
                   end if
                   ! Store image point data
                   this%gp(n)%im%pos=pos_im
                   this%gp(n)%im%ind(1)=i1; this%gp(n)%im%ind(2)=j1; this%gp(n)%im%ind(3)=k1
                   this%gp(n)%im%ind=this%cfg%get_ijk_global(this%gp(n)%im%pos,this%gp(n)%im%ind)
                   this%gp(n)%im%st(1,1)=i1; this%gp(n)%im%st(1,2)=i2
                   this%gp(n)%im%st(2,1)=j1; this%gp(n)%im%st(2,2)=j2
                   this%gp(n)%im%st(3,1)=k1; this%gp(n)%im%st(3,2)=k2
                   this%gp(n)%im%itp=delta3D
                else
                   call die('[gpibm update] Problem computing interpolation weight for image point (centered)')
                end if
             end if
          end do
       end do
    end do

    ! Update label array
    this%label=0
    do n=1,this%ngp
       i=this%gp(n)%ind(1); j=this%gp(n)%ind(2); k=this%gp(n)%ind(3)
       this%label(i,j,k)=+1 !< Ghost point
       i=this%gp(n)%im%ind(1); j=this%gp(n)%im%ind(2); k=this%gp(n)%im%ind(3)
       this%label(i,j,k)=-1 !< Image points
    end do

    ! X-face
    !========================================================================================
    ! Identify ghost points
    this%ngpx=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (Gx(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (Gx(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) this%ngpx=this%ngpx+1
          end do
       end do
    end do
    if (this%ngpx.gt.0) allocate(this%gpx(this%ngpx))
    ! Store ghost points and associated image points
    n=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (Gx(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (Gx(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) then
                n=n+1
                ! Store index of ghost points
                this%gpx(n)%ind(1)=i; this%gpx(n)%ind(2)=j; this%gpx(n)%ind(3)=k
                pos=(/this%cfg%x(i),this%cfg%ym(j),this%cfg%zm(k)/)
                ! Get position and index of image point
                this%gpx(n)%dist=abs(Gx(i,j,k))
                pos_im=pos+2.0_WP*this%gpx(n)%dist*Nx(:,i,j,k)
                ! Find right i index
                i1=i
                do while (pos_im(1)-this%cfg%x(i1  ).lt.0.0_WP.and.i1  .gt.this%cfg%imino_); i1=i1-1; end do
                do while (pos_im(1)-this%cfg%x(i1+1).ge.0.0_WP.and.i1+1.lt.this%cfg%imaxo_); i1=i1+1; end do
                i2=i1+1     
                ! Find right j index
                j1=j
                do while (pos_im(2)-this%cfg%ym(j1  ).lt.0.0_WP.and.j1  .gt.this%cfg%jmino_); j1=j1-1; end do
                do while (pos_im(2)-this%cfg%ym(j1+1).ge.0.0_WP.and.j1+1.lt.this%cfg%jmaxo_); j1=j1+1; end do
                j2=j1+1      
                ! Find right k index
                k1=k
                do while (pos_im(3)-this%cfg%zm(k1  ).lt.0.0_WP.and.k1  .gt.this%cfg%kmino_); k1=k1-1; end do
                do while (pos_im(3)-this%cfg%zm(k1+1).ge.0.0_WP.and.k1+1.lt.this%cfg%kmaxo_); k1=k1+1; end do
                k2=k1+1      
                ! Check if the image point is inside the levelset
                alpha3D = 1.0_WP
                if (Gx(i1,j1,k1).le.0.0_WP) alpha3D(1,1,1) = 0.0_WP
                if (Gx(i2,j1,k1).le.0.0_WP) alpha3D(2,1,1) = 0.0_WP
                if (Gx(i1,j2,k1).le.0.0_WP) alpha3D(1,2,1) = 0.0_WP
                if (Gx(i2,j2,k1).le.0.0_WP) alpha3D(2,2,1) = 0.0_WP
                if (Gx(i1,j1,k2).le.0.0_WP) alpha3D(1,1,2) = 0.0_WP
                if (Gx(i2,j1,k2).le.0.0_WP) alpha3D(2,1,2) = 0.0_WP
                if (Gx(i1,j2,k2).le.0.0_WP) alpha3D(1,2,2) = 0.0_WP
                if (Gx(i2,j2,k2).le.0.0_WP) alpha3D(2,2,2) = 0.0_WP
                buf = sum(alpha3D)
                if (buf.gt.0.0_WP) then
                   ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
                   dist3D(1,1,1)=sqrt((this%cfg%x(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,1,1)=sqrt((this%cfg%x(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,2,1)=sqrt((this%cfg%x(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,2,1)=sqrt((this%cfg%x(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,1,2)=sqrt((this%cfg%x(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,1,2)=sqrt((this%cfg%x(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(1,2,2)=sqrt((this%cfg%x(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,2,2)=sqrt((this%cfg%x(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   delta3D=0.0_WP
                   if (dist3D(1,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,1)=1.0_WP
                   else if (dist3D(2,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,1)=1.0_WP
                   else if (dist3D(1,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,1)=1.0_WP
                   else if (dist3D(2,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,1)=1.0_WP
                   else if (dist3D(1,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,2)=1.0_WP
                   else if (dist3D(2,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,2)=1.0_WP
                   else if (dist3D(1,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,2)=1.0_WP
                   else if (dist3D(2,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,2)=1.0_WP
                   else
                      eta3D=1.0_WP/dist3D**2
                      buf=sum(eta3D*alpha3D)
                      delta3D=alpha3D*eta3D/buf
                   end if
                   ! Store image point data
                   this%gpx(n)%im%pos=pos_im
                   this%gpx(n)%im%ind(1)=i1; this%gpx(n)%im%ind(2)=j1; this%gpx(n)%im%ind(3)=k1
                   this%gpx(n)%im%ind=this%cfg%get_ijk_global(this%gpx(n)%im%pos,this%gpx(n)%im%ind)
                   this%gpx(n)%im%st(1,1)=i1; this%gpx(n)%im%st(1,2)=i2
                   this%gpx(n)%im%st(2,1)=j1; this%gpx(n)%im%st(2,2)=j2
                   this%gpx(n)%im%st(3,1)=k1; this%gpx(n)%im%st(3,2)=k2
                   this%gpx(n)%im%itp=delta3D
                else
                   call die('[gpibm update] Problem computing interpolation weight for image point (X-face)')
                end if
             end if
          end do
       end do
    end do

    ! Y-face
    !========================================================================================
    ! Identify ghost points
    this%ngpy=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (Gy(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (Gy(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) this%ngpy=this%ngpy+1
          end do
       end do
    end do
    if (this%ngpy.gt.0) allocate(this%gpy(this%ngpy))
    ! Store ghost points and associated image points
    n=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (Gy(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (Gy(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) then
                n=n+1
                ! Store index of ghost points
                this%gpy(n)%ind(1)=i; this%gpy(n)%ind(2)=j; this%gpy(n)%ind(3)=k
                pos=(/this%cfg%xm(i),this%cfg%y(j),this%cfg%zm(k)/)
                ! Get position and index of image point
                this%gpy(n)%dist=abs(Gy(i,j,k))
                pos_im=pos+2.0_WP*this%gpy(n)%dist*Ny(:,i,j,k)
                ! Find right i index
                i1=i
                do while (pos_im(1)-this%cfg%xm(i1  ).lt.0.0_WP.and.i1  .gt.this%cfg%imino_); i1=i1-1; end do
                do while (pos_im(1)-this%cfg%xm(i1+1).ge.0.0_WP.and.i1+1.lt.this%cfg%imaxo_); i1=i1+1; end do
                i2=i1+1
                ! Find right j index
                j1=j
                do while (pos_im(2)-this%cfg%y(j1  ).lt.0.0_WP.and.j1  .gt.this%cfg%jmino_); j1=j1-1; end do
                do while (pos_im(2)-this%cfg%y(j1+1).ge.0.0_WP.and.j1+1.lt.this%cfg%jmaxo_); j1=j1+1; end do
                j2=j1+1
                ! Find right k index
                k1=k
                do while (pos_im(3)-this%cfg%zm(k1  ).lt.0.0_WP.and.k1  .gt.this%cfg%kmino_); k1=k1-1; end do
                do while (pos_im(3)-this%cfg%zm(k1+1).ge.0.0_WP.and.k1+1.lt.this%cfg%kmaxo_); k1=k1+1; end do
                k2=k1+1      
                ! Check if the image point is inside the levelset
                alpha3D = 1.0_WP
                if (Gy(i1,j1,k1).le.0.0_WP) alpha3D(1,1,1) = 0.0_WP
                if (Gy(i2,j1,k1).le.0.0_WP) alpha3D(2,1,1) = 0.0_WP
                if (Gy(i1,j2,k1).le.0.0_WP) alpha3D(1,2,1) = 0.0_WP
                if (Gy(i2,j2,k1).le.0.0_WP) alpha3D(2,2,1) = 0.0_WP
                if (Gy(i1,j1,k2).le.0.0_WP) alpha3D(1,1,2) = 0.0_WP
                if (Gy(i2,j1,k2).le.0.0_WP) alpha3D(2,1,2) = 0.0_WP
                if (Gy(i1,j2,k2).le.0.0_WP) alpha3D(1,2,2) = 0.0_WP
                if (Gy(i2,j2,k2).le.0.0_WP) alpha3D(2,2,2) = 0.0_WP
                buf = sum(alpha3D)
                if (buf.gt.0.0_WP) then
                   ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
                   dist3D(1,1,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%y(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,1,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%y(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,2,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%y(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,2,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%y(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,1,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%y(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,1,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%y(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(1,2,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%y(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,2,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%y(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   delta3D=0.0_WP
                   if (dist3D(1,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,1)=1.0_WP
                   else if (dist3D(2,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,1)=1.0_WP
                   else if (dist3D(1,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,1)=1.0_WP
                   else if (dist3D(2,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,1)=1.0_WP
                   else if (dist3D(1,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,2)=1.0_WP
                   else if (dist3D(2,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,2)=1.0_WP
                   else if (dist3D(1,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,2)=1.0_WP
                   else if (dist3D(2,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,2)=1.0_WP
                   else
                      eta3D=1.0_WP/dist3D**2
                      buf=sum(eta3D*alpha3D)
                      delta3D=alpha3D*eta3D/buf
                   end if
                   ! Store image point data
                   this%gpy(n)%im%pos=pos_im
                   this%gpy(n)%im%ind(1)=i1; this%gpy(n)%im%ind(2)=j1; this%gpy(n)%im%ind(3)=k1
                   this%gpy(n)%im%ind=this%cfg%get_ijk_global(this%gpy(n)%im%pos,this%gpy(n)%im%ind)
                   this%gpy(n)%im%st(1,1)=i1; this%gpy(n)%im%st(1,2)=i2
                   this%gpy(n)%im%st(2,1)=j1; this%gpy(n)%im%st(2,2)=j2
                   this%gpy(n)%im%st(3,1)=k1; this%gpy(n)%im%st(3,2)=k2
                   this%gpy(n)%im%itp=delta3D
                else
                   print *, this%cfg%xm(i),this%cfg%y(j),this%cfg%zm(k)
                   call die('[gpibm update] Problem computing interpolation weight for image point (Y-face)')
                end if
             end if
          end do
       end do
    end do

    ! Z-face
    !========================================================================================
    ! Identify ghost points
    this%ngpz=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (Gz(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (Gz(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) this%ngpz=this%ngpz+1
          end do
       end do
    end do
    if (this%ngpz.gt.0) allocate(this%gpz(this%ngpz))
    ! Store ghost points and associated image points
    n=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (Gz(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (Gz(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) then
                n=n+1
                ! Store index of ghost points
                this%gpz(n)%ind(1)=i; this%gpz(n)%ind(2)=j; this%gpz(n)%ind(3)=k
                pos=(/this%cfg%xm(i),this%cfg%ym(j),this%cfg%z(k)/)
                ! Get position and index of image point
                this%gpz(n)%dist=abs(Gz(i,j,k))
                pos_im=pos+2.0_WP*this%gpz(n)%dist*Nz(:,i,j,k)
                ! Find right i index
                i1=i
                do while (pos_im(1)-this%cfg%xm(i1  ).lt.0.0_WP.and.i1  .gt.this%cfg%imino_); i1=i1-1; end do
                do while (pos_im(1)-this%cfg%xm(i1+1).ge.0.0_WP.and.i1+1.lt.this%cfg%imaxo_); i1=i1+1; end do
                i2=i1+1     
                ! Find right j index
                j1=j
                do while (pos_im(2)-this%cfg%ym(j1  ).lt.0.0_WP.and.j1  .gt.this%cfg%jmino_); j1=j1-1; end do
                do while (pos_im(2)-this%cfg%ym(j1+1).ge.0.0_WP.and.j1+1.lt.this%cfg%jmaxo_); j1=j1+1; end do
                j2=j1+1      
                ! Find right k index
                k1=k
                do while (pos_im(3)-this%cfg%z(k1  ).lt.0.0_WP.and.k1  .gt.this%cfg%kmino_); k1=k1-1; end do
                do while (pos_im(3)-this%cfg%z(k1+1).ge.0.0_WP.and.k1+1.lt.this%cfg%kmaxo_); k1=k1+1; end do
                k2=k1+1      
                ! Check if the image point is inside the levelset
                alpha3D = 1.0_WP
                if (Gz(i1,j1,k1).le.0.0_WP) alpha3D(1,1,1) = 0.0_WP
                if (Gz(i2,j1,k1).le.0.0_WP) alpha3D(2,1,1) = 0.0_WP
                if (Gz(i1,j2,k1).le.0.0_WP) alpha3D(1,2,1) = 0.0_WP
                if (Gz(i2,j2,k1).le.0.0_WP) alpha3D(2,2,1) = 0.0_WP
                if (Gz(i1,j1,k2).le.0.0_WP) alpha3D(1,1,2) = 0.0_WP
                if (Gz(i2,j1,k2).le.0.0_WP) alpha3D(2,1,2) = 0.0_WP
                if (Gz(i1,j2,k2).le.0.0_WP) alpha3D(1,2,2) = 0.0_WP
                if (Gz(i2,j2,k2).le.0.0_WP) alpha3D(2,2,2) = 0.0_WP
                buf = sum(alpha3D)
                if (buf.gt.0.0_WP) then
                   ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
                   dist3D(1,1,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%z(k1)-pos_im(3))**2)
                   dist3D(2,1,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%z(k1)-pos_im(3))**2)
                   dist3D(1,2,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%z(k1)-pos_im(3))**2)
                   dist3D(2,2,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%z(k1)-pos_im(3))**2)
                   dist3D(1,1,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%z(k2)-pos_im(3))**2)
                   dist3D(2,1,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%z(k2)-pos_im(3))**2)
                   dist3D(1,2,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%z(k2)-pos_im(3))**2)
                   dist3D(2,2,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%z(k2)-pos_im(3))**2)
                   delta3D=0.0_WP
                   if (dist3D(1,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,1)=1.0_WP
                   else if (dist3D(2,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,1)=1.0_WP
                   else if (dist3D(1,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,1)=1.0_WP
                   else if (dist3D(2,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,1)=1.0_WP
                   else if (dist3D(1,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,2)=1.0_WP
                   else if (dist3D(2,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,2)=1.0_WP
                   else if (dist3D(1,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,2)=1.0_WP
                   else if (dist3D(2,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,2)=1.0_WP
                   else
                      eta3D=1.0_WP/dist3D**2
                      buf=sum(eta3D*alpha3D)
                      delta3D=alpha3D*eta3D/buf
                   end if
                   ! Store image point data
                   this%gpz(n)%im%pos=pos_im
                   this%gpz(n)%im%ind(1)=i1; this%gpz(n)%im%ind(2)=j1; this%gpz(n)%im%ind(3)=k1
                   this%gpz(n)%im%ind=this%cfg%get_ijk_global(this%gpz(n)%im%pos,this%gpz(n)%im%ind)
                   this%gpz(n)%im%st(1,1)=i1; this%gpz(n)%im%st(1,2)=i2
                   this%gpz(n)%im%st(2,1)=j1; this%gpz(n)%im%st(2,2)=j2
                   this%gpz(n)%im%st(3,1)=k1; this%gpz(n)%im%st(3,2)=k2
                   this%gpz(n)%im%itp=delta3D
                else
                   call die('[gpibm update] Problem computing interpolation weight for image point (Z-face)')
                end if
             end if
          end do
       end do
    end do

    ! Clean up
    deallocate(Gx,Gy,Gz,Nx,Ny,Nz,itpr_x,itpr_y,itpr_z)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Ghost point solver on partitioned grid [",a,"]: ",i0," ghost points found")') trim(this%cfg%name),this%ngp
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine update


  !> Enforce boundary condition at ghost points of A using scalar boundary point BP
  subroutine apply_bcond(this,type,BP,A,dir)
    use messager,       only: die
    implicit none
    class(gpibm), intent(inout) :: this
    integer, intent(in) :: type
    real(WP), intent(in) :: BP
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    character(len=*) :: dir
    type(ghost), dimension(:), allocatable :: gp
    integer :: i,j,k,n,ngp
    real(WP) :: IP,dist
    select case(trim(adjustl(dir)))
    case ('U','u','X','x')
       ngp=this%ngpx
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gpx
       end if
    case ('V','v','Y','y')
       ngp=this%ngpy
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gpy
       end if
    case ('W','w','Z','z')
       ngp=this%ngpz
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gpz
       end if
    case default
       ngp=this%ngp
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gp
       end if
    end select
    select case (type)
    case (dirichlet)
       do n=1,ngp
          i=gp(n)%ind(1); j=gp(n)%ind(2); k=gp(n)%ind(3)
          !IP=sum(gp(n)%im%itp*A(gp(n)%im%st(1,1:2),gp(n)%im%st(2,1:2),gp(n)%im%st(3,1:2)))
          A(i,j,k)=BP!2.0_WP*BP-IP
       end do
    case (neumann)
       do n=1,ngp
          i=gp(n)%ind(1); j=gp(n)%ind(2); k=gp(n)%ind(3)
          IP=sum(gp(n)%im%itp*A(gp(n)%im%st(1,1:2),gp(n)%im%st(2,1:2),gp(n)%im%st(3,1:2)))
          A(i,j,k)=IP-2.0_WP*gp(n)%dist*BP
       end do
    case default
       call die('[gpibm apply_bcond] Unknown bcond type')
    end select
    call this%cfg%sync(A)
    if (allocated(gp)) deallocate(gp)
  end subroutine apply_bcond


  !> Enforce Robin boundary condition of the form alpha*A + beta*dA/dn = gamma
  !> Dirichlet: alpha=1, beta=0
  !> Neumann: alpha=0, beta=1
  subroutine apply_robin(this,alpha,beta,gamma,A,dir)
    use messager,       only: die
    implicit none
    class(gpibm), intent(inout) :: this
    real(WP), intent(inout), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:) :: alpha,beta,gamma
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    character(len=*) :: dir
    type(ghost), dimension(:), allocatable :: gp
    integer :: i,j,k,n,ngp
    real(WP) :: IP,dist,alpha_,beta_,gamma_
    integer, dimension(3) :: ind
    real(WP), dimension(3) :: pos
    select case(trim(adjustl(dir)))
    case ('U')
       ngp=this%ngpx
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gpx
       end if
    case ('V')
       ngp=this%ngpy
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gpy
       end if
    case ('W')
       ngp=this%ngpz
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gpz
       end if
    case default
       ngp=this%ngp
       if (ngp.gt.0) then
          allocate(gp(ngp)); gp=this%gp
       end if
    end select
    do n=1,ngp
       i=gp(n)%ind(1); j=gp(n)%ind(2); k=gp(n)%ind(3)
       ! Interpolate coefficients to the boundary
       if (trim(adjustl(dir)).eq.'U') then
          pos(1)=0.5_WP*(this%cfg%x(i)+gp(n)%im%pos(1))
       else
          pos(1)=0.5_WP*(this%cfg%xm(i)+gp(n)%im%pos(1))
       end if
       if (trim(adjustl(dir)).eq.'V') then
          pos(2)=0.5_WP*(this%cfg%y(j)+gp(n)%im%pos(2))
       else
          pos(2)=0.5_WP*(this%cfg%ym(j)+gp(n)%im%pos(2))
       end if
       if (trim(adjustl(dir)).eq.'W') then
          pos(3)=0.5_WP*(this%cfg%z(k)+gp(n)%im%pos(3))
       else
          pos(3)=0.5_WP*(this%cfg%zm(k)+gp(n)%im%pos(3))
       end if
       ind=this%cfg%get_ijk_global(pos,gp(n)%ind)
       alpha_=this%cfg%get_scalar(pos=pos,i0=ind(1),j0=ind(2),k0=ind(3),S=alpha,bc='n')
       beta_ =this%cfg%get_scalar(pos=pos,i0=ind(1),j0=ind(2),k0=ind(3),S=beta ,bc='n')
       gamma_=this%cfg%get_scalar(pos=pos,i0=ind(1),j0=ind(2),k0=ind(3),S=gamma,bc='n')
       ! Get image point
       IP=sum(gp(n)%im%itp*A(gp(n)%im%st(1,1:2),gp(n)%im%st(2,1:2),gp(n)%im%st(3,1:2)))
       A(i,j,k)=(2.0_WP*gamma_-IP*(alpha_+beta_/gp(n)%dist))/(alpha_-beta_/gp(n)%dist)
    end do
    call this%cfg%sync(A)
    if (allocated(gp)) deallocate(gp)
  end subroutine apply_robin


end module gp_class
