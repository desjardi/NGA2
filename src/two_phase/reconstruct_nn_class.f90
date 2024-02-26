!> Neural Network Class for interface reconstructions
module reconstruct_nn_class
    use precision,         only: WP
    use config_class,      only: config
    use multimatrix_class, only: multimatrix
    implicit none
    private

    public :: reconstruct_nn
    public :: ReLU

    type, extends(multimatrix) :: reconstruct_nn
        ! Indices of bias vectors
        integer :: ivec_lay0_bias, ivec_lay1_bias, ivec_lay2_bias, ivec_outp_bias
        ! Indices of weight matrices
        integer :: imat_lay0_weight, imat_lay1_weight, imat_lay2_weight, imat_outp_weight
        ! Transposed weight matrices
        real(WP), dimension(:,:), allocatable :: lay0_weight_T,lay1_weight_T,lay2_weight_T
        real(WP), dimension(:,:), allocatable :: outp_weight_T

        contains
        procedure :: get_normal                                     

    end type reconstruct_nn

    interface reconstruct_nn
      procedure constructor
    end interface reconstruct_nn

contains
    function constructor(cfg,fdata,name) result(self)
        use messager, only: die
        implicit none
        type(reconstruct_nn)                   :: self
        class(config), target, intent(in)      :: cfg
        character(len=*), intent(in)           :: fdata
        character(len=*), intent(in), optional :: name
        integer :: ivector,imatrix

        ! load in weights and biases
        self%multimatrix=multimatrix(cfg=cfg,fdata=fdata,name=name)

        ! Get bias vectors
        do ivector=1,self%nvector
            select case (trim(self%vectors(ivector)%name))
                case ('layers_0_bias')
                    self%ivec_lay0_bias=ivector
                case ('layers_1_bias')
                    self%ivec_lay1_bias=ivector
                case ('layers_2_bias')
                    self%ivec_lay2_bias=ivector
                case ('output_bias')
                    self%ivec_outp_bias=ivector
            end select
        end do
        if (max(self%ivec_lay0_bias,self%ivec_lay1_bias,self%ivec_lay2_bias,self%ivec_outp_bias).ne.self%nvector) then
            call die('[Neural Network] Inconsistent number of vectors')
        end if

        ! Get weight matrices
        do imatrix=1,self%nmatrix
            select case (trim(self%matrices(imatrix)%name))
            case ('layers_0_weight')
                self%imat_lay0_weight=imatrix
            case ('layers_1_weight')
                self%imat_lay1_weight=imatrix
            case ('layers_2_weight')
                self%imat_lay2_weight=imatrix
            case ('output_weight')
                self%imat_outp_weight=imatrix
            end select
        end do
        if(max(self%imat_lay0_weight,self%imat_lay1_weight,self%imat_lay2_weight,self%imat_outp_weight).ne.self%nmatrix) then
            call die('[Neural Network] Inconsistent number of matrices')
        end if

        ! Allocate transposed weight matrices
        allocate(self%lay0_weight_T(size(self%matrices(self%imat_lay0_weight)%matrix,dim=2),size(self%matrices(self%imat_lay0_weight)%matrix,dim=1)))
        allocate(self%lay1_weight_T(size(self%matrices(self%imat_lay1_weight)%matrix,dim=2),size(self%matrices(self%imat_lay1_weight)%matrix,dim=1)))
        allocate(self%lay2_weight_T(size(self%matrices(self%imat_lay2_weight)%matrix,dim=2),size(self%matrices(self%imat_lay2_weight)%matrix,dim=1)))
        allocate(self%outp_weight_T(size(self%matrices(self%imat_outp_weight)%matrix,dim=2),size(self%matrices(self%imat_outp_weight)%matrix,dim=1)))

        ! Transpose weight matrices
        self%lay0_weight_T=transpose(self%matrices(self%imat_lay0_weight)%matrix)
        self%lay1_weight_T=transpose(self%matrices(self%imat_lay1_weight)%matrix)
        self%lay2_weight_T=transpose(self%matrices(self%imat_lay2_weight)%matrix)
        self%outp_weight_T=transpose(self%matrices(self%imat_outp_weight)%matrix)
    end function constructor

    !> Forward pass of Neural Network to return PLIC normal vector
    subroutine get_normal(this,fractions,normal)
        implicit none
        class(reconstruct_nn), intent(inout)   :: this
        real(WP), dimension(:), intent(in)  :: fractions
        real(WP), dimension(:), intent(out) :: normal
        real(WP), dimension(:), allocatable :: tmparr

        tmparr=ReLU(matmul(fractions,this%lay0_weight_T)+this%vectors(this%ivec_lay0_bias)%vector)
        tmparr=ReLU(matmul(tmparr,this%lay1_weight_T)+this%vectors(this%ivec_lay1_bias)%vector)
        tmparr=ReLU(matmul(tmparr,this%lay2_weight_T)+this%vectors(this%ivec_lay2_bias)%vector)
        normal=     matmul(tmparr,this%outp_weight_T)+this%vectors(this%ivec_outp_bias)%vector
        
        deallocate(tmparr)
     end subroutine get_normal

    !> ReLU activation function
    elemental pure function ReLU(x)
        real(WP), intent(in) :: x
        real(WP) :: ReLU
        ReLU=max(0.0_WP,x)
    end function ReLU
    
end module reconstruct_nn_class
  

