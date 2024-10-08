!> Turbulent subgrid scale model class
!> Provides support for evaluating turbulent viscosity, diffusivity, and variance
module sgsmodel_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: sgsmodel
   
   ! List of SGS LES models available
   integer, parameter, public :: dynamic_smag =1    !< Dynamic Smagorinsky
   integer, parameter, public :: constant_smag=2    !< Constant Smagorinsky 
   integer, parameter, public :: vreman       =3    !< Vreman 2004 
   
   !> SGS model object definition
   type :: sgsmodel
      
      ! This is our config
      class(config), pointer :: cfg                             !< This is the config the model is build for
      
      ! Safe index bounds
      integer :: imin_in,jmin_in,kmin_in                        !< Safe min in each direction
      integer :: imax_in,jmax_in,kmax_in                        !< Safe max in each direction

      ! Some clipping parameters
      real(WP) :: Cs_ref=0.17_WP
      
      ! LM and MM tensor norms and eddy viscosity
      real(WP), dimension(:,:,:), allocatable :: LM,MM          !< LM and MM tensor norms
      real(WP), dimension(:,:,:), allocatable :: visc           !< Turbulent eddy viscosity
      
      ! Some information of the fields
      real(WP) :: max_visc                                      !< Maximum eddy viscosity
      real(WP) :: min_visc                                      !< Minimum eddy viscosity
      
      ! Filtering scheme
      real(WP), dimension(:,:,:,:,:,:), allocatable :: filterd  !< Filtering operator with Dirichlet at walls
      real(WP), dimension(:,:,:,:,:,:), allocatable :: filtern  !< Filtering operator with Neumann at walls
      real(WP), dimension(:,:,:), allocatable :: ratio          !< Filter ratio
      real(WP), dimension(:,:,:), allocatable :: delta          !< Filter size
      
   contains
      
      procedure :: log=>sgs_log                                 !< Log SGS info
      procedure :: print=>sgs_print                             !< Output SGS info to the screen
      procedure :: get_visc                                     !< Calls appropriate eddy viscosity subroutine
      procedure :: visc_dynamic                                 !< Calculate the SGS viscosity (Dynamic Smag)
      procedure :: visc_cst                                     !< Calculate the SGS viscosity (Constant Smag)
      procedure :: visc_vreman                                  !< Calculate the SGS viscosity (Vreman 2004)      
      
   end type sgsmodel
   
   
   !> Declare model constructor
   interface sgsmodel
      procedure constructor
   end interface sgsmodel
   
contains
   
   
   !> Default constructor for model
   function constructor(cfg,umask,vmask,wmask) result(self)
      implicit none
      type(sgsmodel) :: self
      class(config), target, intent(in) :: cfg
      integer, dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: umask  !< U-face masks from a velocity solver (needed to differentiate a wall from a Dirichlet)
      integer, dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: vmask  !< V-face masks from a velocity solver
      integer, dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: wmask  !< W-face masks from a velocity solver
      integer :: i,j,k,stx,sty,stz
      real(WP), dimension(:,:,:,:), allocatable :: xfilter,yfilter,zfilter
      
      ! Point to config object
      self%cfg=>cfg
      
      ! Allocate storage for 1D filters
      allocate(xfilter(-1:+1,self%cfg%imin_  :self%cfg%imax_  ,self%cfg%jmin_-1:self%cfg%jmax_+1,self%cfg%kmin_-1:self%cfg%kmax_+1))
      allocate(yfilter(-1:+1,self%cfg%imin_-1:self%cfg%imax_+1,self%cfg%jmin_  :self%cfg%jmax_  ,self%cfg%kmin_-1:self%cfg%kmax_+1))
      allocate(zfilter(-1:+1,self%cfg%imin_-1:self%cfg%imax_+1,self%cfg%jmin_-1:self%cfg%jmax_+1,self%cfg%kmin_  :self%cfg%kmax_  ))
      
      ! Build Simpson's rule 1D filter in 1D with zero gradient assumption at the walls
      do k=self%cfg%kmin_-1,self%cfg%kmax_+1
         do j=self%cfg%jmin_-1,self%cfg%jmax_+1
            do i=self%cfg%imin_  ,self%cfg%imax_
               xfilter(-1,i,j,k)=(-0.5_WP*(self%cfg%dxm(i+1)-self%cfg%dxm(i))-(self%cfg%dxm(i+1)-self%cfg%dxm(i))**2/(6.0_WP*self%cfg%dxm(i  ))+self%cfg%dxm(i+1)/3.0_WP)/(self%cfg%dxm(i+1)+self%cfg%dxm(i))
               xfilter( 0,i,j,k)=2.0_WP/3.0_WP+(self%cfg%dxm(i+1)-self%cfg%dxm(i))**2/(6.0_WP*self%cfg%dxm(i+1)*self%cfg%dxm(i))
               xfilter(+1,i,j,k)=(+0.5_WP*(self%cfg%dxm(i+1)-self%cfg%dxm(i))-(self%cfg%dxm(i+1)-self%cfg%dxm(i))**2/(6.0_WP*self%cfg%dxm(i+1))+self%cfg%dxm(i  )/3.0_WP)/(self%cfg%dxm(i+1)+self%cfg%dxm(i))
               if (self%cfg%VF(i-1,j,k).eq.0.0_WP) then
                  xfilter( 0,i,j,k)=xfilter( 0,i,j,k)+xfilter(-1,i,j,k)
                  xfilter(-1,i,j,k)=0.0_WP
               end if
               if (self%cfg%VF(i+1,j,k).eq.0.0_WP) then
                  xfilter( 0,i,j,k)=xfilter( 0,i,j,k)+xfilter(+1,i,j,k)
                  xfilter(+1,i,j,k)=0.0_WP
               end if
               if (self%cfg%VF(i,j,k).eq.0.0_WP) xfilter(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      do k=self%cfg%kmin_-1,self%cfg%kmax_+1
         do j=self%cfg%jmin_  ,self%cfg%jmax_
            do i=self%cfg%imin_-1,self%cfg%imax_+1
               yfilter(-1,i,j,k)=(-0.5_WP*(self%cfg%dym(j+1)-self%cfg%dym(j))-(self%cfg%dym(j+1)-self%cfg%dym(j))**2/(6.0_WP*self%cfg%dym(j  ))+self%cfg%dym(j+1)/3.0_WP)/(self%cfg%dym(j+1)+self%cfg%dym(j))
               yfilter( 0,i,j,k)=2.0_WP/3.0_WP+(self%cfg%dym(j+1)-self%cfg%dym(j))**2/(6.0_WP*self%cfg%dym(j+1)*self%cfg%dym(j))
               yfilter(+1,i,j,k)=(+0.5_WP*(self%cfg%dym(j+1)-self%cfg%dym(j))-(self%cfg%dym(j+1)-self%cfg%dym(j))**2/(6.0_WP*self%cfg%dym(j+1))+self%cfg%dym(j  )/3.0_WP)/(self%cfg%dym(j+1)+self%cfg%dym(j))
               if (self%cfg%VF(i,j-1,k).eq.0.0_WP) then
                  yfilter( 0,i,j,k)=yfilter( 0,i,j,k)+yfilter(-1,i,j,k)
                  yfilter(-1,i,j,k)=0.0_WP
               end if
               if (self%cfg%VF(i,j+1,k).eq.0.0_WP) then
                  yfilter( 0,i,j,k)=yfilter( 0,i,j,k)+yfilter(+1,i,j,k)
                  yfilter(+1,i,j,k)=0.0_WP
               end if
               if (self%cfg%VF(i,j,k).eq.0.0_WP) yfilter(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      do k=self%cfg%kmin_  ,self%cfg%kmax_
         do j=self%cfg%jmin_-1,self%cfg%jmax_+1
            do i=self%cfg%imin_-1,self%cfg%imax_+1
               zfilter(-1,i,j,k)=(-0.5_WP*(self%cfg%dzm(k+1)-self%cfg%dzm(k))-(self%cfg%dzm(k+1)-self%cfg%dzm(k))**2/(6.0_WP*self%cfg%dzm(k  ))+self%cfg%dzm(k+1)/3.0_WP)/(self%cfg%dxm(k+1)+self%cfg%dzm(k))
               zfilter( 0,i,j,k)=2.0_WP/3.0_WP+(self%cfg%dzm(k+1)-self%cfg%dzm(k))**2/(6.0_WP*self%cfg%dzm(k+1)*self%cfg%dzm(k))
               zfilter(+1,i,j,k)=(+0.5_WP*(self%cfg%dzm(k+1)-self%cfg%dzm(k))-(self%cfg%dzm(k+1)-self%cfg%dzm(k))**2/(6.0_WP*self%cfg%dzm(k+1))+self%cfg%dzm(k  )/3.0_WP)/(self%cfg%dzm(k+1)+self%cfg%dzm(k))
               if (self%cfg%VF(i,j,k-1).eq.0.0_WP) then
                  zfilter( 0,i,j,k)=zfilter( 0,i,j,k)+zfilter(-1,i,j,k)
                  zfilter(-1,i,j,k)=0.0_WP
               end if
               if (self%cfg%VF(i,j,k+1).eq.0.0_WP) then
                  zfilter( 0,i,j,k)=zfilter( 0,i,j,k)+zfilter(+1,i,j,k)
                  zfilter(+1,i,j,k)=0.0_WP
               end if
               if (self%cfg%VF(i,j,k).eq.0.0_WP) zfilter(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Define the 3D filtering operator as the product of 1D filter operators
      allocate(self%filtern(-1:+1,-1:+1,-1:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               do stz=-1,+1
                  do sty=-1,+1
                     do stx=-1,+1
                        self%filtern(stx,sty,stz,i,j,k)=xfilter(stx,i,j+sty,k+stz)*yfilter(sty,i+stx,j,k+stz)*zfilter(stz,i+stx,j+sty,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      
      ! Build Simpson's rule 1D filter in 1D with zero value assumption at the walls
      do k=self%cfg%kmin_-1,self%cfg%kmax_+1
         do j=self%cfg%jmin_-1,self%cfg%jmax_+1
            do i=self%cfg%imin_  ,self%cfg%imax_
               xfilter(-1,i,j,k)=(-0.5_WP*(self%cfg%dxm(i+1)-self%cfg%dxm(i))-(self%cfg%dxm(i+1)-self%cfg%dxm(i))**2/(6.0_WP*self%cfg%dxm(i  ))+self%cfg%dxm(i+1)/3.0_WP)/(self%cfg%dxm(i+1)+self%cfg%dxm(i))
               xfilter( 0,i,j,k)=2.0_WP/3.0_WP+(self%cfg%dxm(i+1)-self%cfg%dxm(i))**2/(6.0_WP*self%cfg%dxm(i+1)*self%cfg%dxm(i))
               xfilter(+1,i,j,k)=(+0.5_WP*(self%cfg%dxm(i+1)-self%cfg%dxm(i))-(self%cfg%dxm(i+1)-self%cfg%dxm(i))**2/(6.0_WP*self%cfg%dxm(i+1))+self%cfg%dxm(i  )/3.0_WP)/(self%cfg%dxm(i+1)+self%cfg%dxm(i))
               if (umask(i,j,k).eq.1) then
                  xfilter( 0,i,j,k)=xfilter( 0,i,j,k)-xfilter(-1,i,j,k)
                  xfilter(-1,i,j,k)=0.0_WP
               end if
               if (umask(i+1,j,k).eq.1) then
                  xfilter( 0,i,j,k)=xfilter( 0,i,j,k)-xfilter(+1,i,j,k)
                  xfilter(+1,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      do k=self%cfg%kmin_-1,self%cfg%kmax_+1
         do j=self%cfg%jmin_  ,self%cfg%jmax_
            do i=self%cfg%imin_-1,self%cfg%imax_+1
               yfilter(-1,i,j,k)=(-0.5_WP*(self%cfg%dym(j+1)-self%cfg%dym(j))-(self%cfg%dym(j+1)-self%cfg%dym(j))**2/(6.0_WP*self%cfg%dym(j  ))+self%cfg%dym(j+1)/3.0_WP)/(self%cfg%dym(j+1)+self%cfg%dym(j))
               yfilter( 0,i,j,k)=2.0_WP/3.0_WP+(self%cfg%dym(j+1)-self%cfg%dym(j))**2/(6.0_WP*self%cfg%dym(j+1)*self%cfg%dym(j))
               yfilter(+1,i,j,k)=(+0.5_WP*(self%cfg%dym(j+1)-self%cfg%dym(j))-(self%cfg%dym(j+1)-self%cfg%dym(j))**2/(6.0_WP*self%cfg%dym(j+1))+self%cfg%dym(j  )/3.0_WP)/(self%cfg%dym(j+1)+self%cfg%dym(j))
               if (vmask(i,j,k).eq.1) then
                  yfilter( 0,i,j,k)=yfilter( 0,i,j,k)-yfilter(-1,i,j,k)
                  yfilter(-1,i,j,k)=0.0_WP
               end if
               if (vmask(i,j+1,k).eq.1) then
                  yfilter( 0,i,j,k)=yfilter( 0,i,j,k)-yfilter(+1,i,j,k)
                  yfilter(+1,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      do k=self%cfg%kmin_  ,self%cfg%kmax_
         do j=self%cfg%jmin_-1,self%cfg%jmax_+1
            do i=self%cfg%imin_-1,self%cfg%imax_+1
               zfilter(-1,i,j,k)=(-0.5_WP*(self%cfg%dzm(k+1)-self%cfg%dzm(k))-(self%cfg%dzm(k+1)-self%cfg%dzm(k))**2/(6.0_WP*self%cfg%dzm(k  ))+self%cfg%dzm(k+1)/3.0_WP)/(self%cfg%dxm(k+1)+self%cfg%dzm(k))
               zfilter( 0,i,j,k)=2.0_WP/3.0_WP+(self%cfg%dzm(k+1)-self%cfg%dzm(k))**2/(6.0_WP*self%cfg%dzm(k+1)*self%cfg%dzm(k))
               zfilter(+1,i,j,k)=(+0.5_WP*(self%cfg%dzm(k+1)-self%cfg%dzm(k))-(self%cfg%dzm(k+1)-self%cfg%dzm(k))**2/(6.0_WP*self%cfg%dzm(k+1))+self%cfg%dzm(k  )/3.0_WP)/(self%cfg%dzm(k+1)+self%cfg%dzm(k))
               if (wmask(i,j,k).eq.1) then
                  zfilter( 0,i,j,k)=zfilter( 0,i,j,k)-zfilter(-1,i,j,k)
                  zfilter(-1,i,j,k)=0.0_WP
               end if
               if (wmask(i,j,k+1).eq.1) then
                  zfilter( 0,i,j,k)=zfilter( 0,i,j,k)-zfilter(+1,i,j,k)
                  zfilter(+1,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Define the 3D filtering operator as the product of 1D filter operators
      allocate(self%filterd(-1:+1,-1:+1,-1:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               do stz=-1,+1
                  do sty=-1,+1
                     do stx=-1,+1
                        self%filterd(stx,sty,stz,i,j,k)=xfilter(stx,i,j+sty,k+stz)*yfilter(sty,i+stx,j,k+stz)*zfilter(stz,i+stx,j+sty,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      
      ! Define filter characteristics
      allocate(self%ratio(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
      allocate(self%delta(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               self%delta(i,j,k)=(self%cfg%dx(i)*self%cfg%dy(j)*self%cfg%dz(k))**(1.0_WP/3.0_WP)
               self%ratio(i,j,k)=((self%cfg%xm(i+1)-self%cfg%xm(i-1))*(self%cfg%ym(j+1)-self%cfg%ym(j-1))*(self%cfg%zm(k+1)-self%cfg%zm(k-1))/(self%cfg%dx(i)*self%cfg%dy(j)*self%cfg%dz(k)))**(1.0_WP/3.0_WP)
            end do
         end do
      end do
      
      ! Deallocate 1D filters
      deallocate(xfilter,yfilter,zfilter)
      
      ! Allocate LM and MM
      allocate(self%LM(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LM=0.0_WP
      allocate(self%MM(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%MM=0.0_WP
      
      ! Allocate visc
      allocate(self%visc(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc=0.0_WP
      
      ! Safe loop extents
      self%imin_in=self%cfg%imino_; if (self%cfg%iproc.eq.1           .and..not.self%cfg%xper) self%imin_in=self%cfg%imin
      self%imax_in=self%cfg%imaxo_; if (self%cfg%iproc.eq.self%cfg%npx.and..not.self%cfg%xper) self%imax_in=self%cfg%imax
      self%jmin_in=self%cfg%jmino_; if (self%cfg%jproc.eq.1           .and..not.self%cfg%yper) self%jmin_in=self%cfg%jmin
      self%jmax_in=self%cfg%jmaxo_; if (self%cfg%jproc.eq.self%cfg%npy.and..not.self%cfg%yper) self%jmax_in=self%cfg%jmax
      self%kmin_in=self%cfg%kmino_; if (self%cfg%kproc.eq.1           .and..not.self%cfg%zper) self%kmin_in=self%cfg%kmin
      self%kmax_in=self%cfg%kmaxo_; if (self%cfg%kproc.eq.self%cfg%npz.and..not.self%cfg%zper) self%kmax_in=self%cfg%kmax
      
   end function constructor
   
   
   !> Calls appropriate eddy viscosity subroutine according to SGSmodel%type
   subroutine get_visc(this,type,dt,rho,Ui,Vi,Wi,SR,gradu)
     use messager, only: die
     use param, only: verbose
     implicit none
     class(sgsmodel), intent(inout) :: this
     integer, intent(in) :: type !< Model type
     real(WP), intent(in) :: dt !< dt since the last call to the model
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho !< Density including all ghosts
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: Ui  !< Interpolated velocities including all ghosts
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: Vi  !< Interpolated velocities including all ghosts
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: Wi  !< Interpolated velocities including all ghosts
     real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: SR  !< Strain rate tensor
     real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: gradu !< Velocity gradient tensor

     select case(type)
     case(dynamic_smag)
        if (.not.present(Ui).or..not.(present(Vi)).or..not.present(Wi).or..not.present(SR)) &
             call die('[sgs get_visc] Dynamic Smagorinsky model requires Ui, Vi, Wi, and SR')
        call this%visc_dynamic(dt,rho,Ui,Vi,Wi,SR,gradu)
     case(constant_smag)
        if (.not.present(SR)) call die('[sgs get_visc] Constant Smagorinsky model requires SR')
        call this%visc_cst(rho,SR)
     case(vreman)
        if (.not.present(gradu)) call die('[sgs get_visc] Vreman model requires gradu')
        call this%visc_vreman(rho,gradu)
     end select

     ! Calculate some info on the model
     calc_info: block
       use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
       use parallel, only: MPI_REAL_WP
       integer :: ierr
       call MPI_ALLREDUCE(maxval(this%visc),this%max_visc,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
       call MPI_ALLREDUCE(minval(this%visc),this%min_visc,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
     end block calc_info
      
     ! Output info about model
     if (verbose.gt.0) call this%log()
     if (verbose.gt.1) call this%print()

   end subroutine get_visc
   
   
   !> Get subgrid scale dynamic viscosity - Dynamic
   subroutine visc_dynamic(this,dt,rho,Ui,Vi,Wi,SR,gradu)
      implicit none
      class(sgsmodel), intent(inout) :: this
      real(WP), intent(in) :: dt !< dt since the last call to the model
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho !< Density including all ghosts
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Ui  !< Interpolated velocities including all ghosts
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Vi  !< Interpolated velocities including all ghosts
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Wi  !< Interpolated velocities including all ghosts
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR  !< Strain rate tensor
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradu !< Velocity gradient tensor
      integer :: i,j,k
      real(WP) :: Frho,FU,FV,FW,FS_
      real(WP) :: Cs,tau,alpha,interp
      real(WP), dimension(3) :: pos
      real(WP), dimension(6) :: FSR,FrhoS_SR,FrhoUU,Mij,Lij
      real(WP), dimension(:,:,:), allocatable :: S_,LMold,MMold
      
      ! Remember the previous LM and MM
      allocate(LMold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); LMold=this%LM
      allocate(MMold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); MMold=this%MM
      
      ! Prepare magnitude of SR tensor
      allocate(S_(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      S_=sqrt(SR(1,:,:,:)**2+SR(2,:,:,:)**2+SR(3,:,:,:)**2+2.0_WP*(SR(4,:,:,:)**2+SR(5,:,:,:)**2+SR(6,:,:,:)**2))
      
      ! Compute the favre-filtered LM and MM tensors
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall cells
               if (this%cfg%VF(i,j,k).lt.epsilon(1.0_WP)) then
                  this%LM(i,j,k)=0.0_WP
                  this%MM(i,j,k)=0.0_WP
                  cycle
               end if
               ! Filtered density (Neumann)
               Frho=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1))
               ! Favre-filtered velocity (Dirichlet)
               FU=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Ui(i-1:i+1,j-1:j+1,k-1:k+1))/Frho
               FV=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Vi(i-1:i+1,j-1:j+1,k-1:k+1))/Frho
               FW=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Wi(i-1:i+1,j-1:j+1,k-1:k+1))/Frho
               ! Filtered strain rate tensor (Neumann)
               FSR(1)=sum(this%filtern(:,:,:,i,j,k)*SR(1,i-1:i+1,j-1:j+1,k-1:k+1))
               FSR(2)=sum(this%filtern(:,:,:,i,j,k)*SR(2,i-1:i+1,j-1:j+1,k-1:k+1))
               FSR(3)=sum(this%filtern(:,:,:,i,j,k)*SR(3,i-1:i+1,j-1:j+1,k-1:k+1))
               FSR(4)=sum(this%filtern(:,:,:,i,j,k)*SR(4,i-1:i+1,j-1:j+1,k-1:k+1))
               FSR(5)=sum(this%filtern(:,:,:,i,j,k)*SR(5,i-1:i+1,j-1:j+1,k-1:k+1))
               FSR(6)=sum(this%filtern(:,:,:,i,j,k)*SR(6,i-1:i+1,j-1:j+1,k-1:k+1))
               FS_=sqrt(sum(FSR(1:3)**2+2.0_WP*FSR(4:6)**2))
               ! Filtered rho*S_*SR (Neumann)
               FrhoS_SR(1)=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*S_(i-1:i+1,j-1:j+1,k-1:k+1)*SR(1,i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoS_SR(2)=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*S_(i-1:i+1,j-1:j+1,k-1:k+1)*SR(2,i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoS_SR(3)=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*S_(i-1:i+1,j-1:j+1,k-1:k+1)*SR(3,i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoS_SR(4)=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*S_(i-1:i+1,j-1:j+1,k-1:k+1)*SR(4,i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoS_SR(5)=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*S_(i-1:i+1,j-1:j+1,k-1:k+1)*SR(5,i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoS_SR(6)=sum(this%filtern(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*S_(i-1:i+1,j-1:j+1,k-1:k+1)*SR(6,i-1:i+1,j-1:j+1,k-1:k+1))
               ! Filtered rho*U*U (Dirichlet)
               FrhoUU(1)=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Ui(i-1:i+1,j-1:j+1,k-1:k+1)*Ui(i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoUU(2)=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Vi(i-1:i+1,j-1:j+1,k-1:k+1)*Vi(i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoUU(3)=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Wi(i-1:i+1,j-1:j+1,k-1:k+1)*Wi(i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoUU(4)=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Ui(i-1:i+1,j-1:j+1,k-1:k+1)*Vi(i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoUU(5)=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Vi(i-1:i+1,j-1:j+1,k-1:k+1)*Wi(i-1:i+1,j-1:j+1,k-1:k+1))
               FrhoUU(6)=sum(this%filterd(:,:,:,i,j,k)*rho(i-1:i+1,j-1:j+1,k-1:k+1)*Wi(i-1:i+1,j-1:j+1,k-1:k+1)*Ui(i-1:i+1,j-1:j+1,k-1:k+1))
               ! Compute Mij=<Frho*S_*SR>-ratio^2<Frho><FS_><FSR>
               Mij=2.0_WP*this%delta(i,j,k)**2*(FrhoS_SR-this%ratio(i,j,k)**2*Frho*FS_*FSR)
               ! Compute Lij=<rho*U*U>-<rho><U><U>
               Lij(1)=FrhoUU(1)-Frho*FU*FU
               Lij(2)=FrhoUU(2)-Frho*FV*FV
               Lij(3)=FrhoUU(3)-Frho*FW*FW
               Lij(4)=FrhoUU(4)-Frho*FU*FV
               Lij(5)=FrhoUU(5)-Frho*FV*FW
               Lij(6)=FrhoUU(6)-Frho*FW*FU
               ! Compute LM=sum(Lij Mij) and MM=sum(Mij Mij)
               this%LM(i,j,k)=sum(Lij(1:3)*Mij(1:3)+2.0_WP*Lij(4:6)*Mij(4:6))
               this%MM(i,j,k)=sum(Mij(1:3)*Mij(1:3)+2.0_WP*Mij(4:6)*Mij(4:6))
            end do
         end do
      end do
      
      ! Apply Meneveau's Lagrangian averaging strategy
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Simple Lagrangian backtracking
               pos(1)=this%cfg%xm(i)-Ui(i,j,k)*dt
               pos(2)=this%cfg%ym(j)-Vi(i,j,k)*dt
               pos(3)=this%cfg%zm(k)-Wi(i,j,k)*dt
               ! Advance LM and MM
               tau=dt*(LMold(i,j,k)*MMold(i,j,k))**0.125_WP/(1.5_WP*this%delta(i,j,k))
               alpha=tau/(1.0_WP+tau)
               interp=this%cfg%get_scalar(pos=pos,i0=i,j0=j,k0=k,S=LMold,bc='n')
               this%LM(i,j,k)=alpha*this%LM(i,j,k)+(1.0_WP-alpha)*interp
               interp=this%cfg%get_scalar(pos=pos,i0=i,j0=j,k0=k,S=MMold,bc='n')
               this%MM(i,j,k)=alpha*this%MM(i,j,k)+(1.0_WP-alpha)*interp
               ! Safe limits
               this%LM(i,j,k)=max(this%LM(i,j,k),100.0_WP*epsilon(1.0_WP))
               this%MM(i,j,k)=max(this%MM(i,j,k),100.0_WP*epsilon(1.0_WP)/this%Cs_ref**2)
            end do
         end do
      end do
      
      ! Synchronize LM and MM
      call this%cfg%sync(this%LM)
      call this%cfg%sync(this%MM)
      
      ! Compute the eddy viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (abs(this%MM(i,j,k)).gt.100.0_WP*epsilon(1.0_WP)/this%Cs_ref**2) then
                  Cs=max(this%LM(i,j,k)/this%MM(i,j,k),0.0_WP)
               else
                  Cs=0.0_WP
               end if
               this%visc(i,j,k)=rho(i,j,k)*S_(i,j,k)*Cs*this%delta(i,j,k)**2
            end do
         end do
      end do
      ! Synchronize visc
      call this%cfg%sync(this%visc)
      
      ! Deallocate work arrays
      deallocate(LMold,MMold,S_)
      
   end subroutine visc_dynamic
   
   
   !> Get subgrid scale viscosity - Constant coefficient
   subroutine visc_cst(this,rho,SR)
      implicit none
      class(sgsmodel), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho   !< Density including all ghosts
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR !< Strain rate tensor
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: S_
      
      ! Prepare magnitude of SR tensor
      allocate(S_(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      S_=sqrt(SR(1,:,:,:)**2+SR(2,:,:,:)**2+SR(3,:,:,:)**2+2.0_WP*(SR(4,:,:,:)**2+SR(5,:,:,:)**2+SR(6,:,:,:)**2))
      
      ! Compute the eddy viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%visc(i,j,k)=rho(i,j,k)*S_(i,j,k)*(this%Cs_ref*this%delta(i,j,k))**2
            end do
         end do
      end do
      
      ! Synchronize visc
      call this%cfg%sync(this%visc)
      
      ! Deallocate work arrays
      deallocate(S_)
      
   end subroutine visc_cst
   
   
   !> Get subgrid scale dynamic viscosity - Vreman
   subroutine visc_vreman(this,rho,gradu)
      implicit none
      class(sgsmodel), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho         !< Density including all ghosts
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradu !< Velocity gradient
      real(WP) :: A,B,C
      real(WP), dimension(1:3,1:3) :: beta
      integer :: i,j,k,ii,jj
      
      ! Model constant is c=2.5*Cs_ref**2
      ! Vreman uses c=0.07 which corresponds to Cs_ref=0.17
      C=2.5_WP*this%Cs_ref**2

      ! Compute the eddy viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Compute A=gradu_ij*gradu_ij invariant
               A=gradu(1,1,i,j,k)**2+gradu(1,2,i,j,k)**2+gradu(1,3,i,j,k)**2+&
               & gradu(2,1,i,j,k)**2+gradu(2,2,i,j,k)**2+gradu(2,3,i,j,k)**2+&
               & gradu(3,1,i,j,k)**2+gradu(3,2,i,j,k)**2+gradu(3,3,i,j,k)**2
               ! Compute beta_ij=dx_m*dx_m*gradu_mi*gradu_mj
               do jj=1,3
                  do ii=1,3
                     beta(ii,jj)=this%cfg%dx(i)**2*gradu(1,ii,i,j,k)*gradu(1,jj,i,j,k)&
                     &          +this%cfg%dy(j)**2*gradu(2,ii,i,j,k)*gradu(2,jj,i,j,k)&
                     &          +this%cfg%dz(k)**2*gradu(3,ii,i,j,k)*gradu(3,jj,i,j,k)
                  end do
               end do
               ! Compute B invariant
               B=beta(1,1)*beta(2,2)-beta(1,2)**2&
               &+beta(1,1)*beta(3,3)-beta(1,3)**2&
               &+beta(2,2)*beta(3,3)-beta(2,3)**2
               ! Assemble algebraic eddy viscosity model
               if (B.lt.1.0e-8_WP) then
                  this%visc(i,j,k)=0.0_WP
               else
                  this%visc(i,j,k)=rho(i,j,k)*C*sqrt(B/A)
               end if
            end do
         end do
      end do
      
      ! Synchronize visc
      call this%cfg%sync(this%visc)
      
   end subroutine visc_vreman
   
   
   !> Log info for model
   subroutine sgs_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(sgsmodel), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("SGS turbulence modeling for config [",a,"] - visc range = [",es12.5,",",es12.5,"]")') trim(this%cfg%name),this%min_visc,this%max_visc; call log(message)
      end if
   end subroutine sgs_log
   
   
   !> Print out info for model
   subroutine sgs_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(sgsmodel), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("SGS turbulence modeling for config [",a,"] - visc range = [",es12.5,",",es12.5,"]")') trim(this%cfg%name),this%min_visc,this%max_visc
      end if
      
   end subroutine sgs_print
   
   
end module sgsmodel_class
