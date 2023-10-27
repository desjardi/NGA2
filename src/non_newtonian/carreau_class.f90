!> Calculates shear rate-dependent viscosity using Carreau model
module carreau_class
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: carreau
   
   !> Carreau model object definition
   type :: carreau
      ! This is our config
      class(config), pointer :: cfg                        !< This is the config the model is build for
      ! Model parameters
      real(WP) :: ncoeff                                   !< Powerlaw coefficient
      real(WP) :: tref                                     !< Reference fluid timescale
      real(WP) :: visc_zero                                !< Zero shear rate  viscosity
      real(WP) :: visc_inf=0.0_WP                          !< Infinite shear rate viscosity (0 by default)
      ! Viscosity
      real(WP), dimension(:,:,:), allocatable :: visc      !< Carreau viscosity
      ! Monitoring info
      real(WP) :: visc_min,visc_max                        !< Min and max viscosity
   contains
      procedure :: initialize                              !< Initialization of Carreau fluid class
      procedure :: update_visc                             !< Update viscosity given strain rate tensor
   end type carreau
   
contains
   
   
   !> Carreau model initialization
   subroutine initialize(this,cfg)
      implicit none
      class(carreau), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      ! Store pointer to cfg
      this%cfg=>cfg
      ! Allocate storage for rate dependent viscosity
      allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc=0.0_WP
   end subroutine initialize
   

   !> Compute visc from SR using Carreau model
   subroutine update_visc(this,SR)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(carreau), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR
      real(WP) :: SRmag,my_visc_max,my_visc_min
      integer :: i,j,k,ierr
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Compute second invariant of strain rate tensor = sqrt(2*SR**2)
               SRmag=sqrt(2.0_WP*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2)))
               ! Compute viscosity
               this%visc(i,j,k)=this%visc_inf+(this%visc_zero-this%visc_inf)*(1.0_WP+(this%tref*SRmag)**2)**(0.5_WP*this%ncoeff-0.5_WP)
            end do
         end do
      end do
      my_visc_max=maxval(this%visc); call MPI_ALLREDUCE(my_visc_max,this%visc_max,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_visc_min=minval(this%visc); call MPI_ALLREDUCE(my_visc_min,this%visc_min,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine update_visc
   
   
end module carreau_class