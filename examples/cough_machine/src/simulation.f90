!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use block1_class,  only: block1,filmthickness_over_dx,min_filmthickness,diam_over_filmthickness,max_eccentricity,d_threshold
   use block2_class,  only: block2
   use coupler_class, only: coupler
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> Block 1 and 2 objects
   type(block1) :: b1
   type(block2) :: b2

   !> Couplers between blocks
   type(coupler) :: cpl12x,cpl12y,cpl12z
   type(coupler) :: cpl21x,cpl21y,cpl21z
   
   !> Storage for coupled fields
   real(WP), dimension(:,:,:), allocatable :: U1on2,V1on2,W1on2
   real(WP), dimension(:,:,:), allocatable :: U2on1,V2on1,W2on1

contains


   !> Initialization of problem solver
   subroutine simulation_init
      use geometry, only: cfg1,cfg2
      implicit none

      ! Initialize both blocks
      b1%cfg=>cfg1; call b1%init()
      b2%cfg=>cfg2; call b2%init()

      ! Initialize the couplers
      coupler_prep: block
         use parallel, only: group
         ! Block 1 to block 2
         cpl12x=coupler(src_grp=group,dst_grp=group,name='in_to_out_x'); call cpl12x%set_src(cfg1,'x'); call cpl12x%set_dst(cfg2,'x'); call cpl12x%initialize()
         cpl12y=coupler(src_grp=group,dst_grp=group,name='in_to_out_y'); call cpl12y%set_src(cfg1,'y'); call cpl12y%set_dst(cfg2,'y'); call cpl12y%initialize()
         cpl12z=coupler(src_grp=group,dst_grp=group,name='in_to_out_z'); call cpl12z%set_src(cfg1,'z'); call cpl12z%set_dst(cfg2,'z'); call cpl12z%initialize()
         allocate(U1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); U1on2=0.0_WP
         allocate(V1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); V1on2=0.0_WP
         allocate(W1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); W1on2=0.0_WP
         ! Block 2 to block 1
         cpl21x=coupler(src_grp=group,dst_grp=group,name='out_to_in_x'); call cpl21x%set_src(cfg2,'x'); call cpl21x%set_dst(cfg1,'x'); call cpl21x%initialize()
         cpl21y=coupler(src_grp=group,dst_grp=group,name='out_to_in_y'); call cpl21y%set_src(cfg2,'y'); call cpl21y%set_dst(cfg1,'y'); call cpl21y%initialize()
         cpl21z=coupler(src_grp=group,dst_grp=group,name='out_to_in_z'); call cpl21z%set_src(cfg2,'z'); call cpl21z%set_dst(cfg1,'z'); call cpl21z%initialize()
         allocate(U2on1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); U2on1=0.0_WP
         allocate(V2on1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); V2on1=0.0_WP
         allocate(W2on1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); W2on1=0.0_WP
      end block coupler_prep

      ! Setup nudging region in block 2
      b2%nudge_trans=20.0_WP*b2%cfg%min_meshsize
      b2%nudge_xmin =-1.0_WP !b1%cfg%x(b1%cfg%imin)
      b2%nudge_xmax =b1%cfg%x(b1%cfg%imax+1)
      b2%nudge_ymin =b1%cfg%y(b1%cfg%jmin)
      b2%nudge_ymax =b1%cfg%y(b1%cfg%jmax+1)
      b2%nudge_zmin =-1.0_WP !b1%cfg%z(b1%cfg%kmin)
      b2%nudge_zmax =+1.0_WP !b1%cfg%z(b1%cfg%kmax+1)

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      integer :: i,j,k

      ! Perform time integration - block 1 is the main driver here
      do while (.not.b1%time%done())

         ! Advance block 1
         call b1%step(U2on1,V2on1,W2on1)

         ! ###############################################
         ! ####### TRANSFER DROPLETS FROM 1->2 HERE ######
         ! ###############################################
         ! Perform CCL and transfer
         call transfer_vf_to_drops()
         ! After we're done clip all VOF at the exit area and along the sides - hopefully nothing's left
         !do k=b1%fs%cfg%kmino_,b1%fs%cfg%kmaxo_
         !   do j=b1%fs%cfg%jmino_,b1%fs%cfg%jmaxo_
         !      do i=b1%fs%cfg%imino_,b1%fs%cfg%imaxo_
         !         if (i.ge.b1%vf%cfg%imax-5) b1%vf%VF(i,j,k)=0.0_WP
         !         if (j.ge.b1%vf%cfg%jmax-5) b1%vf%VF(i,j,k)=0.0_WP
         !         if (j.le.b1%vf%cfg%jmin+5) b1%vf%VF(i,j,k)=0.0_WP
         !         if (k.ge.b1%vf%cfg%kmax-5) b1%vf%VF(i,j,k)=0.0_WP
         !         if (k.le.b1%vf%cfg%kmin+5) b1%vf%VF(i,j,k)=0.0_WP
         !      end do
         !   end do
         !end do

         ! Advance block 2 until we've caught up
         do while (b2%time%t.lt.b1%time%t)

            ! Exchange data using cpl12x/y/z couplers and the most recent velocity
            U1on2=0.0_WP; call cpl12x%push(b1%fs%U); call cpl12x%transfer(); call cpl12x%pull(U1on2)
            V1on2=0.0_WP; call cpl12y%push(b1%fs%V); call cpl12y%transfer(); call cpl12y%pull(V1on2)
            W1on2=0.0_WP; call cpl12z%push(b1%fs%W); call cpl12z%transfer(); call cpl12z%pull(W1on2)

            ! Advance block 2
            call b2%step(U1on2,V1on2,W1on2)
            
            ! Exchange data using cpl21x/y/z couplers and the most recent velocity
            U2on1=0.0_WP; call cpl21x%push(b2%fs%U); call cpl21x%transfer(); call cpl21x%pull(U2on1)
            V2on1=0.0_WP; call cpl21y%push(b2%fs%V); call cpl21y%transfer(); call cpl21y%pull(V2on1)
            W2on1=0.0_WP; call cpl21z%push(b2%fs%W); call cpl21z%transfer(); call cpl21z%pull(W2on1)
            
         end do

      end do

   end subroutine simulation_run

   !> Transfer vf to drops
      subroutine transfer_vf_to_drops()
         implicit none

         ! Perform a first pass with simplest CCL
         call b1%cc1%build_lists(VF=b1%vf%VF,U=b1%fs%U,V=b1%fs%V,W=b1%fs%W)

         ! Loop through identified detached structs and remove those that are spherical enough
         remove_struct: block
            use mathtools, only: pi
            integer :: m,n,l,i,j,k,np
            real(WP) :: lmin,lmax,eccentricity,diam

            ! Loops over film segments contained locally
            do m=1,b1%cc1%n_meta_struct

               ! Test if sphericity is compatible with transfer
               lmin=b1%cc1%meta_structures_list(m)%lengths(3)
               if (lmin.eq.0.0_WP) lmin=b1%cc1%meta_structures_list(m)%lengths(2) ! Handle 2D case
               lmax=b1%cc1%meta_structures_list(m)%lengths(1)
               eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
               if (eccentricity.gt.max_eccentricity) cycle

               ! Test if diameter is compatible with transfer
               diam=(6.0_WP*b1%cc1%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
               if (diam.eq.0.0_WP.or.diam.gt.d_threshold) cycle

               ! Create drop from available liquid volume - only one root does that
               if (b1%cc1%cfg%amRoot) then
                  ! Make room for new drop
                  np=b2%lp%np_+1; call b2%lp%resize(np)
                  ! Add the drop
                  b2%lp%p(np)%id  =int(0,8)                                                                                               !< Give id (maybe based on break-up model?)
                  b2%lp%p(np)%dt  =0.0_WP                                                                                                 !< Let the drop find it own integration time
                  b2%lp%p(np)%d   =diam                                                                                                   !< Assign diameter to account for full volume
                  b2%lp%p(np)%pos =[b1%cc1%meta_structures_list(m)%x,b1%cc1%meta_structures_list(m)%y,b1%cc1%meta_structures_list(m)%z]   !< Place the drop at the liquid barycenter
                  b2%lp%p(np)%vel =[b1%cc1%meta_structures_list(m)%u,b1%cc1%meta_structures_list(m)%v,b1%cc1%meta_structures_list(m)%w]   !< Assign mean structure velocity as drop velocity
                  b2%lp%p(np)%ind =b2%lp%cfg%get_ijk_global(b2%lp%p(np)%pos,[b2%lp%cfg%imin,b2%lp%cfg%jmin,b2%lp%cfg%kmin])               !< Place the drop in the proper cell for the lp%cfg
                  b2%lp%p(np)%flag=0                                                                                                      !< Activate it
                  ! Increment particle counter
                  b2%lp%np_=np
               end if

               ! Find local structs with matching id
               do n=b1%cc1%sync_offset+1,b1%cc1%sync_offset+b1%cc1%n_struct
                  if (b1%cc1%struct_list(b1%cc1%struct_map_(n))%parent.ne.b1%cc1%meta_structures_list(m)%id) cycle
                  ! Remove liquid in meta-structure cells
                  do l=1,b1%cc1%struct_list(b1%cc1%struct_map_(n))%nnode ! Loops over cells within local
                     i=b1%cc1%struct_list(b1%cc1%struct_map_(n))%node(1,l)
                     j=b1%cc1%struct_list(b1%cc1%struct_map_(n))%node(2,l)
                     k=b1%cc1%struct_list(b1%cc1%struct_map_(n))%node(3,l)
                     ! Remove liquid in that cell
                     b1%vf%VF(i,j,k)=0.0_WP
                  end do
               end do

            end do

         end block remove_struct

         ! Sync VF and clean up IRL and band
         call b1%vf%cfg%sync(b1%vf%VF)
         call b1%vf%clean_irl_and_band()

         ! Clean up CCL
         call b1%cc1%deallocate_lists()

         ! Perform more detailed CCL in a second pass
         b1%cc1%max_interface_planes=2
         call b1%cc1%build_lists(VF=b1%vf%VF,poly=b1%vf%interface_polygon,U=b1%fs%U,V=b1%fs%V,W=b1%fs%W)
         call b1%cc1%get_min_thickness()
         call b1%cc1%sort_by_thickness()

         ! Loop through identified films and remove those that are thin enough
         remove_film: block
            use mathtools, only: pi
            integer :: m,n,i,j,k,np,ip,np_old
            real(WP) :: Vt,Vl,Hl,Vd

            ! Loops over film segments contained locally
            do m=b1%cc1%film_sync_offset+1,b1%cc1%film_sync_offset+b1%cc1%n_film

               ! Skip non-liquid films
               if (b1%cc1%film_list(b1%cc1%film_map_(m))%phase.ne.1) cycle

               ! Skip films that are still thick enough
               if (b1%cc1%film_list(b1%cc1%film_map_(m))%min_thickness.gt.min_filmthickness) cycle

               ! We are still here: transfer the film to drops
               Vt=0.0_WP       ! Transferred volume
               Vl=0.0_WP       ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
               np_old=b2%lp%np_  ! Remember old number of particles
               do n=1,b1%cc1%film_list(b1%cc1%film_map_(m))%nnode ! Loops over cells within local film segment
                  i=b1%cc1%film_list(b1%cc1%film_map_(m))%node(1,n)
                  j=b1%cc1%film_list(b1%cc1%film_map_(m))%node(2,n)
                  k=b1%cc1%film_list(b1%cc1%film_map_(m))%node(3,n)
                  ! Increment liquid volume to remove
                  Vl=Vl+b1%vf%VF(i,j,k)*b1%vf%cfg%vol(i,j,k)
                  ! Estimate drop size based on local film thickness in current cell
                  Hl=max(b1%cc1%film_thickness(i,j,k),min_filmthickness)
                  Vd=pi/6.0_WP*(diam_over_filmthickness*Hl)**3
                  ! Create drops from available liquid volume
                  do while (Vl-Vd.gt.0.0_WP)
                     ! Make room for new drop
                     np=b2%lp%np_+1; call b2%lp%resize(np)
                     ! Add the drop
                     b2%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                     b2%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                     b2%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                     b2%lp%p(np)%pos =b1%vf%Lbary(:,i,j,k)                         !< Place the drop at the liquid barycenter
                     b2%lp%p(np)%vel =b1%fs%cfg%get_velocity(pos=b2%lp%p(np)%pos,i0=i,j0=j,k0=k,U=b1%fs%U,V=b1%fs%V,W=b1%fs%W) !< Interpolate local cell velocity as drop velocity
                     b2%lp%p(np)%ind =b2%lp%cfg%get_ijk_global(b2%lp%p(np)%pos,[b2%lp%cfg%imin,b2%lp%cfg%jmin,b2%lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                     b2%lp%p(np)%flag=0                                          !< Activate it
                     ! Increment particle counter
                     b2%lp%np_=np
                     ! Update tracked volumes
                     Vl=Vl-Vd
                     Vt=Vt+Vd
                  end do
                  ! Remove liquid in that cell
                  b1%vf%VF(i,j,k)=0.0_WP
               end do

               ! Based on how many particles were created, decide what to do with left-over volume
               if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                  ! Add one last drop for remaining liquid volume
                  np=b2%lp%np_+1; call b2%lp%resize(np)
                  ! Add the drop
                  b2%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                  b2%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  b2%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
                  b2%lp%p(np)%pos =b1%vf%Lbary(:,i,j,k)                         !< Place the drop at the liquid barycenter
                  b2%lp%p(np)%vel =b1%fs%cfg%get_velocity(pos=b2%lp%p(np)%pos,i0=i,j0=j,k0=k,U=b1%fs%U,V=b1%fs%V,W=b1%fs%W) !< Interpolate local cell velocity as drop velocity
                  b2%lp%p(np)%ind =b2%lp%cfg%get_ijk_global(b2%lp%p(np)%pos,[b2%lp%cfg%imin,b2%lp%cfg%jmin,b2%lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                  b2%lp%p(np)%flag=0                                          !< Activate it
                  ! Increment particle counter
                  b2%lp%np_=np
               else ! Some particles were created, make them all larger
                  do ip=np_old+1,b2%lp%np_
                     b2%lp%p(ip)%d=b2%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
                  end do
               end if
            end do

         end block remove_film

         ! Sync VF and clean up IRL and band
         call b1%vf%cfg%sync(b1%vf%VF)
         call b1%vf%clean_irl_and_band()

         ! Clean up CCL
         call b1%cc1%deallocate_lists()

         ! Resync the spray
         call b2%lp%sync()

      end subroutine transfer_vf_to_drops


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Deallocate work arrays for both blocks
      call b1%final()
      call b2%final()

   end subroutine simulation_final


end module simulation
