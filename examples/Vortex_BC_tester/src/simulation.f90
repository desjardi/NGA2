!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use block1_class,  only: block1
   !use block2_class,  only: block2
   !use coupler_class, only: coupler
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> Block 1 and 2 objects
   type(block1) :: b1
   !type(block2) :: b2

   !> Couplers between blocks
   !type(coupler) :: cpl12x,cpl12y,cpl12z        !> block 1 to block 2 couplers
   !type(coupler) :: cpl21x,cpl21y,cpl21z        !> block 2 to block 1 couplers

   !> Storage for coupled fields
   !real(WP), dimension(:,:,:), allocatable :: U1on2,V1on2,W1on2
   !real(WP), dimension(:,:,:), allocatable :: U2on1,V2on1,W2on1

contains


   !> Initialization of problem solver
   subroutine simulation_init
      use geometry, only: cfg1
      implicit none

      ! Initialize both blocks
      b1%cfg=>cfg1; call b1%init()
      !b2%cfg=>cfg2; call b2%init()

      ! Initialize the couplers
      !coupler_prep: block
         !use parallel, only: group
         ! Create block1-to-block2 couplers
         !cpl12x=coupler(src_grp=group,dst_grp=group,name='in_to_out_x'); call cpl12x%set_src(cfg1,'x'); call cpl12x%set_dst(cfg2,'x'); call cpl12x%initialize()
         !cpl12y=coupler(src_grp=group,dst_grp=group,name='in_to_out_y'); call cpl12y%set_src(cfg1,'y'); call cpl12y%set_dst(cfg2,'y'); call cpl12y%initialize()
         !cpl12z=coupler(src_grp=group,dst_grp=group,name='in_to_out_z'); call cpl12z%set_src(cfg1,'z'); call cpl12z%set_dst(cfg2,'z'); call cpl12z%initialize()
         !allocate(U1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         !allocate(V1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         !allocate(W1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         ! Create block2-to-block1 couplers
         !cpl21x=coupler(src_grp=group,dst_grp=group,name='out_to_in_x'); call cpl21x%set_src(cfg2,'x'); call cpl21x%set_dst(cfg1,'x'); call cpl21x%initialize()
         !cpl21y=coupler(src_grp=group,dst_grp=group,name='out_to_in_y'); call cpl21y%set_src(cfg2,'y'); call cpl21y%set_dst(cfg1,'y'); call cpl21y%initialize()
         !cpl21z=coupler(src_grp=group,dst_grp=group,name='out_to_in_z'); call cpl21z%set_src(cfg2,'z'); call cpl21z%set_dst(cfg1,'z'); call cpl21z%initialize()
         !allocate(U2on1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         !allocate(V2on1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         !allocate(W2on1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
      !end block coupler_prep

      ! Setup nudging region in block 2
      !b2%nudge_trans=10.0_WP*b2%cfg%min_meshsize
      !b2%nudge_xmin =-1.0_WP !b1%cfg%x(b1%cfg%imin)
      !b2%nudge_xmax =b1%cfg%x(b1%cfg%imax+1)
      !b2%nudge_ymin =b1%cfg%y(b1%cfg%jmin)
      !b2%nudge_ymax =b1%cfg%y(b1%cfg%jmax+1)
      !b2%nudge_zmin =-1.0_WP !b1%cfg%z(b1%cfg%kmin)
      !b2%nudge_zmax =+1.0_WP !b1%cfg%z(b1%cfg%kmax+1)

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      integer :: i,j,k

      ! Perform time integration - block 1 is the main driver here
      do while (.not.b1%time%done())

         ! Advance block 1
         call b1%step()


         ! Advance block 2 until we've caught up
         !do while (b2%time%t.lt.b1%time%t)

            ! Exchange data using cpl12x/y/z couplers and the most recent velocity
            !U1on2=0.0_WP; call cpl12x%push(b1%fs%U); call cpl12x%transfer(); call cpl12x%pull(U1on2)
            !V1on2=0.0_WP; call cpl12y%push(b1%fs%V); call cpl12y%transfer(); call cpl12y%pull(V1on2)
            !W1on2=0.0_WP; call cpl12z%push(b1%fs%W); call cpl12z%transfer(); call cpl12z%pull(W1on2)

            ! Advance block 2
            !call b2%step(U1on2,V1on2,W1on2)

         !end do

      end do

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Deallocate work arrays for both blocks
      call b1%final()
      !call b2%final()

   end subroutine simulation_final


end module simulation
