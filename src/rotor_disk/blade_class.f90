module blade_class
    use precision, only: WP
    implicit none
    private

    public :: blade

    type :: blade

        integer :: Nr                               ! Number of radial points
        integer :: Na                               ! Number of azimuthal points

        real(WP), dimension(:), allocatable :: radius 
        real(WP), dimension(:), allocatable :: chord
        real(WP), dimension(:), allocatable :: twist            ! Twist angle in degrees

        real(WP), dimension(:), allocatable :: aoa              ! Angle of attack in degrees
        real(WP), dimension(:), allocatable :: cl
        real(WP), dimension(:), allocatable :: cd


    contains

        procedure :: read

        procedure :: interpolate_r     ! Interpolate blade data in radial direction, given radius, return chord and twist
        procedure :: interpolate_a      ! Interpolate blade data in azimuthal direction, given aoa, return cl and cd



    end type blade

    interface blade
        procedure blade_constructor
    end interface blade
    
contains


    function blade_constructor(Nr, Na) result(self)
        implicit none
        type(blade) :: self
        integer, intent(in) :: Nr, Na

        self%Nr = Nr
        self%Na = Na

        ! allocate arrays
        allocate(self%radius(Nr))
        allocate(self%chord(Nr))
        allocate(self%twist(Nr))

        allocate(self%aoa(Na))
        allocate(self%cl(Na))
        allocate(self%cd(Na))

    end function blade_constructor


    subroutine read(self, radius, chord, twist, aoa, cl, cd)
        implicit none
        class(blade), intent(inout) :: self
        
        real(WP), dimension(:), intent(in) :: radius, chord, twist, aoa, cl, cd

        self%radius = radius
        self%chord = chord
        self%twist = twist

        self%aoa = aoa
        self%cl = cl
        self%cd = cd
    end subroutine read



    subroutine interpolate_r(self, r, c, t)
        implicit none
        class(blade), intent(inout) :: self
        real(WP), intent(in) :: r
        real(WP), intent(out) :: c
        real(WP), intent(out) :: t

        integer :: i
        real(WP) :: r_corrected

        r_corrected = r
        if (r < minval(self%radius)) r_corrected = minval(self%radius)
        if (r > maxval(self%radius)) r_corrected = maxval(self%radius)


        do i = 1, self%Nr-1
            if ((self%radius(i) <= r_corrected) .and. (r_corrected <= self%radius(i+1))) then
                c = self%chord(i) + (self%chord(i+1) - self%chord(i)) * (r_corrected - self%radius(i)) / (self%radius(i+1) - self%radius(i))
                t = self%twist(i) + (self%twist(i+1) - self%twist(i)) * (r_corrected - self%radius(i)) / (self%radius(i+1) - self%radius(i))
                return
            end if
        end do

    end subroutine interpolate_r
    

    subroutine interpolate_a(self, a, cl, cd)
        implicit none
        class(blade), intent(inout) :: self
        real(WP), intent(in) :: a               ! Angle of attack in degrees
        real(WP), intent(out) :: cl
        real(WP), intent(out) :: cd

        integer :: i
        real(WP) :: a_corrected

        a_corrected = a
        if (a < minval(self%aoa)) a_corrected = minval(self%aoa)
        if (a > maxval(self%aoa)) a_corrected = maxval(self%aoa)



        do i=1, self%Na-1
            if ((self%aoa(i) <= a_corrected) .and. (a_corrected) <= (self%aoa(i+1))) then
                cl = self%cl(i) + (self%cl(i+1) - self%cl(i)) * (a_corrected - self%aoa(i)) / (self%aoa(i+1) - self%aoa(i))
                cd = self%cd(i) + (self%cd(i+1) - self%cd(i)) * (a_corrected - self%aoa(i)) / (self%aoa(i+1) - self%aoa(i))
                return
            end if
        end do

        
    end subroutine interpolate_a



end module blade_class