!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use pfft3d_class,      only: pfft3d
  use ensight_class,     only: ensight
  implicit none
  private

  type(pfft3d), public :: ps

  !> Ensight postprocessing
  type(ensight)  :: ens_out

  public :: simulation_init, simulation_run, simulation_final

  !> Problem constants
  integer :: npts
  real(WP), dimension(:,:), allocatable :: ptsrcs

contains

  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none

    ! Read params
    read_input: block
      use messager, only: die
      integer :: n

      call param_read('Number of point sources', npts)
      allocate(ptsrcs(npts,4))
      call param_read('Point source x', ptsrcs(:,1))
      call param_read('Point source y', ptsrcs(:,2))
      call param_read('Point source z', ptsrcs(:,3))
      call param_read('Point source strengths', ptsrcs(:,4))

      do n = 1, npts
        if (ptsrcs(n,1) .lt. 0.0_WP .or. ptsrcs(n,1) .gt. cfg%xL .or.         &
            ptsrcs(n,2) .lt. 0.0_WP .or. ptsrcs(n,2) .gt. cfg%yL .or.         &
            ptsrcs(n,3) .lt. 0.0_WP .or. ptsrcs(n,3) .gt. cfg%zL) then
          call die("point not in domain")
        end if
      end do

    end block read_input

    ! Create a single-phase flow solver without bconds
    create_and_initialize_flow_solver: block
      integer :: i, j, k

      ! create solver
      ps = pfft3d(cfg=cfg, nst=7, name='P3dfft Test')
      call ps%init()

      ! set stencil
      ps%stc(1,:) = (/  0,  0,  0 /)
      ps%stc(2,:) = (/ +1,  0,  0 /)
      ps%stc(3,:) = (/ -1,  0,  0 /)
      ps%stc(4,:) = (/  0, +1,  0 /)
      ps%stc(5,:) = (/  0, -1,  0 /)
      ps%stc(6,:) = (/  0,  0, +1 /)
      ps%stc(7,:) = (/  0,  0, -1 /)
      do k = cfg%kmin_, cfg%kmax_
        do j = cfg%jmin_, cfg%jmax_
          do i = cfg%imin_, cfg%imax_
            ps%opr(2,i,j,k) = cfg%dxmi(i)**2
            ps%opr(3,i,j,k) = cfg%dxmi(i)**2
            ps%opr(4,i,j,k) = cfg%dymi(j)**2
            ps%opr(5,i,j,k) = cfg%dymi(j)**2
            ps%opr(6,i,j,k) = cfg%dzmi(k)**2
            ps%opr(7,i,j,k) = cfg%dzmi(k)**2
            ps%opr(1,i,j,k) = - sum(ps%opr(2:7,i,j,k))
          end do
        end do
      end do

      ! setup
      call ps%setup()

    end block create_and_initialize_flow_solver

    ! Add Ensight output
    create_ensight: block

      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg, name='PoissionTest')

      ! Add variables to output
      call ens_out%add_scalar('rhs',ps%rhs)
      call ens_out%add_scalar('sol',ps%sol)

    end block create_ensight

  end subroutine simulation_init

  !> Solve problem
  subroutine simulation_run
    implicit none

    ! set up rhs
    ! this sucks, but it doesn't matter
    rhs_setup: block
      integer :: n, i, j, k, ic, jc, kc
      real(WP) :: d2, d2c

      ps%rhs(:,:,:) = 0.0_WP

      do n = 1, npts
        d2c = huge(d2c)
        do k = cfg%kmin, cfg%kmax
          do j = cfg%jmin, cfg%jmax
            do i = cfg%imin, cfg%imax
              d2 =      (cfg%xm(i) - ptsrcs(n,1))**2
              d2 = d2 + (cfg%ym(j) - ptsrcs(n,2))**2
              d2 = d2 + (cfg%zm(k) - ptsrcs(n,3))**2
              if (d2c .gt. d2) then
                d2c = d2; ic = i; jc = j; kc = k;
              end if
            end do
          end do
        end do
        if (cfg%imin_ .le. ic .and. ic .le. cfg%imax_ .and.                   &
            cfg%jmin_ .le. jc .and. jc .le. cfg%jmax_ .and.                   &
            cfg%kmin_ .le. kc .and. kc .le. cfg%kmax_) then
          ps%rhs(ic,jc,kc) = ps%rhs(ic,jc,kc) + ptsrcs(n,4)
        end if
      end do

    end block rhs_setup

    ! Solve Poisson equation
    call ps%solve()

    ! Output to ensight
    call ens_out%write_data(0.0_WP)

  end subroutine simulation_run

  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Get rid of all objects - need destructors
    ! monitor
    ! ensight
    ! bcond
    ! timetracker

    ! Deallocate work arrays
    deallocate(ptsrcs)

  end subroutine simulation_final

end module simulation

