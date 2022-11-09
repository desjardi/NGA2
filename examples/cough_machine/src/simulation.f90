!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,             only: str_medium
   use precision,          only: WP
   use block1_class,       only: block1
   use block2_class,       only: block2,filmthickness_over_dx,min_filmthickness,diam_over_filmthickness,max_eccentricity,d_threshold
   use block3_class,       only: block3
   use coupler_class,      only: coupler
   use event_class,        only: event
   use timetracker_class,  only: timetracker
   use datafile_class,     only: datafile
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> Block 1, 2 and 3 objects
   type(block1) :: b1
   type(block2) :: b2
   type(block3) :: b3

   !> Couplers between blocks
   type(coupler) :: cpl12x,cpl12y,cpl12z
   type(coupler) :: cpl23x,cpl23y,cpl23z
   type(coupler) :: cpl32x,cpl32y,cpl32z

   !> Event tracker for saving restarts and logical constant for checking restarts
   type(event)   :: save_evt
   logical, public :: restarted
   
   !> Storage for coupled fields
   real(WP), dimension(:,:,:), allocatable :: U1on2,V1on2,W1on2
   real(WP), dimension(:,:,:), allocatable :: U2on3,V2on3,W2on3
   real(WP), dimension(:,:,:), allocatable :: U3on2,V3on2,W3on2

   !> How much volume is clipped
   real(WP) :: Vclipped, V_before, V_after

   !>Parameter to compare diameters of junk
   real(WP) :: diam_thresh=1e-6

contains


   !> Initialization of problem solver
   subroutine simulation_init
      use geometry, only: cfg1,cfg2,cfg3
      use param, only: param_read
      implicit none

      ! Handle restart/saves here
      restart_and_save: block
         character(len=str_medium) :: dir_restart
         if (cfg2%amRoot) call execute_command_line('mkdir -p restart')
         ! CAREFUL - WE NEED TO CREATE THE TIMETRACKER BEFORE THE EVENT !
         b2%time=timetracker(cfg2%amRoot,name='cough_machine_in')
         ! Create event for saving restart files
         save_evt=event(b2%time,'Restart output')
         call param_read('Restart output period',save_evt%tper)
         ! Check if we are restarting
         call param_read(tag='Restart from',val=dir_restart,short='r',default='')
         restarted=.false.; if (len_trim(dir_restart).gt.0) restarted=.true.
         if (restarted) then
            ! If we are, read the name of the directory
            call param_read('Restart from',dir_restart,'r')
            ! Read in filenames
            b1%df=datafile(pg=cfg1,fdata='./restart/data1_'//trim(adjustl(dir_restart)))
            b2%df=datafile(pg=cfg2,fdata='./restart/data2_'//trim(adjustl(dir_restart)))
            b3%df=datafile(pg=cfg3,fdata='./restart/data3_'//trim(adjustl(dir_restart)))
            b3%lpt_file='./restart/datalpt_'//trim(adjustl(dir_restart))
         else
            ! If we are not restarting, we will still need datafiles for saving restart files
            b1%df=datafile(pg=cfg1,filename=trim(cfg1%name),nval=2,nvar=6)
            b1%df%valname(1)='t'
            b1%df%valname(2)='dt'
            b1%df%varname(1)='U'
            b1%df%varname(2)='V'
            b1%df%varname(3)='W'
            b1%df%varname(4)='P'
            b1%df%varname(5)='LM'
            b1%df%varname(6)='MM'
            b2%df=datafile(pg=cfg2,filename=trim(cfg2%name),nval=2,nvar=7)
            b2%df%valname(1)='t'
            b2%df%valname(2)='dt'
            b2%df%varname(1)='U'
            b2%df%varname(2)='V'
            b2%df%varname(3)='W'
            b2%df%varname(4)='P'
            b2%df%varname(5)='LM'
            b2%df%varname(6)='MM'
            b2%df%varname(7)='VF'
            b3%df=datafile(pg=cfg3,filename=trim(cfg3%name),nval=2,nvar=6)
            b3%df%valname(1)='t'
            b3%df%valname(2)='dt'
            b3%df%varname(1)='U'
            b3%df%varname(2)='V'
            b3%df%varname(3)='W'
            b3%df%varname(4)='P'
            b3%df%varname(5)='LM'
            b3%df%varname(6)='MM'
         end if
      end block restart_and_save

      ! Initialize both blocks
      b1%cfg=>cfg1; call b1%init(restarted)
      b2%cfg=>cfg2; call b2%init(restarted)
      b3%cfg=>cfg3; call b3%init(restarted)

      ! Initialize the couplers
      coupler_prep: block
         use parallel, only: group
         ! Block 1 to block 2
         cpl12x=coupler(src_grp=group,dst_grp=group,name='duct_to_in_x'); call cpl12x%set_src(cfg1,'x'); call cpl12x%set_dst(cfg2,'x'); call cpl12x%initialize()
         cpl12y=coupler(src_grp=group,dst_grp=group,name='duct_to_in_y'); call cpl12y%set_src(cfg1,'y'); call cpl12y%set_dst(cfg2,'y'); call cpl12y%initialize()
         cpl12z=coupler(src_grp=group,dst_grp=group,name='duct_to_in_z'); call cpl12z%set_src(cfg1,'z'); call cpl12z%set_dst(cfg2,'z'); call cpl12z%initialize()
         allocate(U1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); U1on2=0.0_WP
         allocate(V1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); V1on2=0.0_WP
         allocate(W1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); W1on2=0.0_WP
         ! Block 2 to block 3
         cpl23x=coupler(src_grp=group,dst_grp=group,name='in_to_out_x'); call cpl23x%set_src(cfg2,'x'); call cpl23x%set_dst(cfg3,'x'); call cpl23x%initialize()
         cpl23y=coupler(src_grp=group,dst_grp=group,name='in_to_out_y'); call cpl23y%set_src(cfg2,'y'); call cpl23y%set_dst(cfg3,'y'); call cpl23y%initialize()
         cpl23z=coupler(src_grp=group,dst_grp=group,name='in_to_out_z'); call cpl23z%set_src(cfg2,'z'); call cpl23z%set_dst(cfg3,'z'); call cpl23z%initialize()
         allocate(U2on3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_)); U2on3=0.0_WP
         allocate(V2on3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_)); V2on3=0.0_WP
         allocate(W2on3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_)); W2on3=0.0_WP
         ! Block 3 to block 2
         cpl32x=coupler(src_grp=group,dst_grp=group,name='out_to_in_x'); call cpl32x%set_src(cfg3,'x'); call cpl32x%set_dst(cfg2,'x'); call cpl32x%initialize()
         cpl32y=coupler(src_grp=group,dst_grp=group,name='out_to_in_y'); call cpl32y%set_src(cfg3,'y'); call cpl32y%set_dst(cfg2,'y'); call cpl32y%initialize()
         cpl32z=coupler(src_grp=group,dst_grp=group,name='out_to_in_z'); call cpl32z%set_src(cfg3,'z'); call cpl32z%set_dst(cfg2,'z'); call cpl32z%initialize()
         allocate(U3on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); U3on2=0.0_WP
         allocate(V3on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); V3on2=0.0_WP
         allocate(W3on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); W3on2=0.0_WP
      end block coupler_prep


      ! Setup nudging region in block 2
      b2%nudge_trans=20.0_WP*b2%cfg%min_meshsize
      b2%nudge_xmin =b1%cfg%x(3*b1%cfg%nx/4)
      b2%nudge_xmax =b1%cfg%x(3*b1%cfg%nx/4)
      b2%nudge_ymin =b1%cfg%y(b2%cfg%jmin)
      b2%nudge_ymax =b1%cfg%y(b2%cfg%jmax+1)
      b2%nudge_zmin =b1%cfg%z(b2%cfg%kmin)
      b2%nudge_zmax =b1%cfg%z(b2%cfg%kmax+1)

      ! Setup nudging region in block 3
      b3%nudge_trans=20.0_WP*b3%cfg%min_meshsize
      b3%nudge_xmin =-1.0_WP !b2%cfg%x(b2%cfg%imin)
      b3%nudge_xmax =b2%cfg%x(b2%cfg%imax+1)
      b3%nudge_ymin =b2%cfg%y(b2%cfg%jmin)
      b3%nudge_ymax =b2%cfg%y(b2%cfg%jmax+1)
      b3%nudge_zmin =b2%cfg%z(b2%cfg%kmin)
      b3%nudge_zmax =b2%cfg%z(b2%cfg%kmax+1)
      
   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      integer :: i,j,k
      character(len=str_medium) :: dirname,timestamp

      ! Perform time integration - block 2 is the main driver here
      do while (.not.b2%time%done())

         ! Advance block 2
         call b2%step(U1on2,V1on2,W1on2)

         ! ###############################################
         ! ####### TRANSFER DROPLETS FROM 2->3 HERE ######
         ! ###############################################
         ! Perform CCL and transfer
         ! call transfer_vf_to_drops()

         ! Preclip
         call b2%vf%cfg%integrate(b2%vf%VF,V_before)

         ! After we're done clip all VOF at the exit area and along the sides - hopefully nothing's left
         do k=b2%fs%cfg%kmino_,b2%fs%cfg%kmaxo_
           do j=b2%fs%cfg%jmino_,b2%fs%cfg%jmaxo_
              do i=b2%fs%cfg%imino_,b2%fs%cfg%imaxo_
                 if (i.ge.b2%vf%cfg%imax-5) b2%vf%VF(i,j,k)=0.0_WP
                 if (j.ge.b2%vf%cfg%jmax-5) b2%vf%VF(i,j,k)=0.0_WP
                 if (j.le.b2%vf%cfg%jmin+5) b2%vf%VF(i,j,k)=0.0_WP
                 if (k.ge.b2%vf%cfg%kmax-5) b2%vf%VF(i,j,k)=0.0_WP
                 if (k.le.b2%vf%cfg%kmin+5) b2%vf%VF(i,j,k)=0.0_WP
              end do
           end do
         end do

         ! Preclip
         call b2%vf%cfg%integrate(b2%vf%VF,V_after)
         Vclipped = V_after - V_before

         ! Advance block 1 until we've caught up
         do while (b1%time%t.lt.b2%time%t)

            ! Advance block 1
            call b1%step()

            ! Exchange data using cpl12x/y/z couplers and the most recent velocity
            U1on2=0.0_WP; call cpl12x%push(b1%fs%U); call cpl12x%transfer(); call cpl12x%pull(U1on2)
            V1on2=0.0_WP; call cpl12y%push(b1%fs%V); call cpl12y%transfer(); call cpl12y%pull(V1on2)
            W1on2=0.0_WP; call cpl12z%push(b1%fs%W); call cpl12z%transfer(); call cpl12z%pull(W1on2)

         end do

         ! ! Advance block 3 until we've caught up
         ! do while (b3%time%t.lt.b2%time%t)

         !    ! Exchange data using cpl23x/y/z couplers and the most recent velocity
         !    U2on3=0.0_WP; call cpl23x%push(b2%fs%U); call cpl23x%transfer(); call cpl23x%pull(U2on3)
         !    V2on3=0.0_WP; call cpl23y%push(b2%fs%V); call cpl23y%transfer(); call cpl23y%pull(V2on3)
         !    W2on3=0.0_WP; call cpl23z%push(b2%fs%W); call cpl23z%transfer(); call cpl23z%pull(W2on3)

         !    ! Advance block 3
         !    call b3%step(U2on3,V2on3,W2on3)
            
         !    ! Exchange data using cpl32x/y/z couplers and the most recent velocity
         !    U3on2=0.0_WP; call cpl32x%push(b3%fs%U); call cpl32x%transfer(); call cpl32x%pull(U3on2)
         !    V3on2=0.0_WP; call cpl32y%push(b3%fs%V); call cpl32y%transfer(); call cpl32y%pull(V3on2)
         !    W3on2=0.0_WP; call cpl32z%push(b3%fs%W); call cpl32z%transfer(); call cpl32z%pull(W3on2)
            
         ! end do

         ! ###############################################
         ! ############## PERFORM I/O HERE ###############
         ! ###############################################
         
         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            ! Prefix for files
            write(timestamp,'(es12.5)') b2%time%t
            ! Populate df1 and write it
            call b1%df%pushval(name=   't',val=b1%time%t )
            call b1%df%pushval(name=  'dt',val=b1%time%dt)
            call b1%df%pushvar(name=   'U',var=b1%fs%U   )
            call b1%df%pushvar(name=   'V',var=b1%fs%V   )
            call b1%df%pushvar(name=   'W',var=b1%fs%W   )
            call b1%df%pushvar(name=   'P',var=b1%fs%P   )
            call b1%df%pushvar(name=  'LM',var=b1%sgs%LM )
            call b1%df%pushvar(name=  'MM',var=b1%sgs%MM )
            call b1%df%write(fdata='./restart/data1_'//trim(adjustl(timestamp)))
            ! Populate df2 and write it
            call b2%df%pushval(name=  't',val=b2%time%t )
            call b2%df%pushval(name= 'dt',val=b2%time%dt)
            call b2%df%pushvar(name=  'U',var=b2%fs%U   )
            call b2%df%pushvar(name=  'V',var=b2%fs%V   )
            call b2%df%pushvar(name=  'W',var=b2%fs%W   )
            call b2%df%pushvar(name=  'P',var=b2%fs%P   )
            call b2%df%pushvar(name= 'LM',var=b2%sgs%LM )
            call b2%df%pushvar(name= 'MM',var=b2%sgs%MM )
            call b2%df%pushvar(name= 'VF',var=b2%vf%VF  )
            call b2%df%write(fdata='./restart/data2_'//trim(adjustl(timestamp)))
            ! Populate df3 and write it
            call b3%df%pushval(name=   't',val=b3%time%t )
            call b3%df%pushval(name=  'dt',val=b3%time%dt)
            call b3%df%pushvar(name=   'U',var=b3%fs%U   )
            call b3%df%pushvar(name=   'V',var=b3%fs%V   )
            call b3%df%pushvar(name=   'W',var=b3%fs%W   )
            call b3%df%pushvar(name=   'P',var=b3%fs%P   )
            call b3%df%pushvar(name=  'LM',var=b3%sgs%LM )
            call b3%df%pushvar(name=  'MM',var=b3%sgs%MM )
            call b3%df%write(fdata='./restart/data3_'//trim(adjustl(timestamp)))
            ! Also output particles
            call b3%lp%write(filename='./restart/datalpt_'//trim(adjustl(timestamp)))
         end if
         
      end do

   end subroutine simulation_run

   !> Transfer vf to drops
      subroutine transfer_vf_to_drops()
         implicit none

         ! Perform a first pass with simplest CCL
         call b2%cc2%build_lists(VF=b2%vf%VF,U=b2%fs%U,V=b2%fs%V,W=b2%fs%W)

         ! Loop through identified detached structs and remove those that are spherical enough
         remove_struct: block
            use mathtools, only: pi
            integer :: m,n,l,i,j,k,np
            real(WP) :: lmin,lmax,eccentricity,diam         

            ! Loops over film segments contained locally
            do m=1,b2%cc2%n_meta_struct

               ! Skip 0 volume structs
               if(b2%cc2%meta_structures_list(m)%vol.lt.(((pi*diam_thresh)**3)/6.0_WP)) then
                     ! Find local structs with matching id
                     do n=b2%cc2%sync_offset+1,b2%cc2%sync_offset+b2%cc2%n_struct
                        if (b2%cc2%struct_list(b2%cc2%struct_map_(n))%parent.ne.b2%cc2%meta_structures_list(m)%id) cycle
                           ! Record 0 vol struc timestep, time, ID and volume
                           b2%cc2%zero_struct_id =b2%cc2%meta_structures_list(m)%id
                           b2%cc2%zero_struct_vol=b2%cc2%meta_structures_list(m)%vol
                           call b2%volfile%write()
                        ! Remove liquid in meta-structure cells
                        do l=1,b2%cc2%struct_list(b2%cc2%struct_map_(n))%nnode ! Loops over cells within local
                           i=b2%cc2%struct_list(b2%cc2%struct_map_(n))%node(1,l)
                           j=b2%cc2%struct_list(b2%cc2%struct_map_(n))%node(2,l)
                           k=b2%cc2%struct_list(b2%cc2%struct_map_(n))%node(3,l)
                           ! Remove liquid in that cell
                           b2%vf%VF(i,j,k)=0.0_WP
                        end do
                     end do
                     cycle
               end if
               
               ! Test if sphericity is compatible with transfer
               lmin=b2%cc2%meta_structures_list(m)%lengths(3)
               if (lmin.eq.0.0_WP) lmin=b2%cc2%meta_structures_list(m)%lengths(2) ! Handle 2D case
               lmax=b2%cc2%meta_structures_list(m)%lengths(1)+epsilon(0.0_WP)
               eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
               if (eccentricity.gt.max_eccentricity) cycle

               ! Test if diameter is compatible with transfer
               diam=(6.0_WP*b2%cc2%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
               if (diam.eq.0.0_WP.or.diam.gt.d_threshold) cycle

               ! Create drop from available liquid volume - only one root does that
               if (b2%cc2%cfg%amRoot) then
                  ! Make room for new drop
                  np=b3%lp%np_+1; call b3%lp%resize(np)
                  ! Add the drop
                  b3%lp%p(np)%id  =int(0,8)                                                                                               !< Give id (maybe based on break-up model?)
                  b3%lp%p(np)%dt  =0.0_WP                                                                                                 !< Let the drop find it own integration time
                  b3%lp%p(np)%d   =diam                                                                                                   !< Assign diameter to account for full volume
                  b3%lp%p(np)%pos =[b2%cc2%meta_structures_list(m)%x,b2%cc2%meta_structures_list(m)%y,b2%cc2%meta_structures_list(m)%z]   !< Place the drop at the liquid barycenter
                  b3%lp%p(np)%vel =[b2%cc2%meta_structures_list(m)%u,b2%cc2%meta_structures_list(m)%v,b2%cc2%meta_structures_list(m)%w]   !< Assign mean structure velocity as drop velocity
                  b3%lp%p(np)%ind =b3%lp%cfg%get_ijk_global(b3%lp%p(np)%pos,[b3%lp%cfg%imin,b3%lp%cfg%jmin,b3%lp%cfg%kmin])               !< Place the drop in the proper cell for the lp%cfg
                  b3%lp%p(np)%flag=0                                                                                                      !< Activate it
                  ! Increment particle counter
                  b3%lp%np_=np
               end if

               ! Find local structs with matching id
               do n=b2%cc2%sync_offset+1,b2%cc2%sync_offset+b2%cc2%n_struct
                  if (b2%cc2%struct_list(b2%cc2%struct_map_(n))%parent.ne.b2%cc2%meta_structures_list(m)%id) cycle
                  ! Remove liquid in meta-structure cells
                  do l=1,b2%cc2%struct_list(b2%cc2%struct_map_(n))%nnode ! Loops over cells within local
                     i=b2%cc2%struct_list(b2%cc2%struct_map_(n))%node(1,l)
                     j=b2%cc2%struct_list(b2%cc2%struct_map_(n))%node(2,l)
                     k=b2%cc2%struct_list(b2%cc2%struct_map_(n))%node(3,l)
                     ! Remove liquid in that cell
                     b2%vf%VF(i,j,k)=0.0_WP
                  end do
               end do

            end do


         end block remove_struct

         ! Sync VF and clean up IRL and band
         call b2%vf%cfg%sync(b2%vf%VF)
         call b2%vf%clean_irl_and_band()

         ! Clean up CCL
         call b2%cc2%deallocate_lists()

         ! Perform more detailed CCL in a second pass
         b2%cc2%max_interface_planes=2
         call b2%cc2%build_lists(VF=b2%vf%VF,poly=b2%vf%interface_polygon,U=b2%fs%U,V=b2%fs%V,W=b2%fs%W)
         call b2%cc2%get_min_thickness()
         call b2%cc2%sort_by_thickness()

         ! Loop through identified films and remove those that are thin enough
         remove_film: block
            use mathtools, only: pi
            integer :: m,n,i,j,k,np,ip,np_old
            real(WP) :: Vt,Vl,Hl,Vd,diameter

            ! Loops over film segments contained locally
            do m=b2%cc2%film_sync_offset+1,b2%cc2%film_sync_offset+b2%cc2%n_film

               ! Skip non-liquid films
               if (b2%cc2%film_list(b2%cc2%film_map_(m))%phase.ne.1) cycle

               ! Skip films that are still thick enough
               if (b2%cc2%film_list(b2%cc2%film_map_(m))%min_thickness.gt.min_filmthickness) cycle

               ! We are still here: transfer the film to drops
               Vt=0.0_WP       ! Transferred volume
               Vl=0.0_WP       ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
               np_old=b3%lp%np_  ! Remember old number of particles
               do n=1,b2%cc2%film_list(b2%cc2%film_map_(m))%nnode ! Loops over cells within local film segment
                  i=b2%cc2%film_list(b2%cc2%film_map_(m))%node(1,n)
                  j=b2%cc2%film_list(b2%cc2%film_map_(m))%node(2,n)
                  k=b2%cc2%film_list(b2%cc2%film_map_(m))%node(3,n)
                  ! Increment liquid volume to remove
                  Vl=Vl+b2%vf%VF(i,j,k)*b2%vf%cfg%vol(i,j,k)
                  ! Estimate drop size based on local film thickness in current cell
                  Hl=max(b2%cc2%film_thickness(i,j,k),min_filmthickness)
                  diameter=max(diam_over_filmthickness*Hl,diam_thresh)
                  Vd=pi/6.0_WP*diameter**3
                  !Vd=pi/6.0_WP*(diam_over_filmthickness*Hl)**3
                  ! Create drops from available liquid volume
                  do while (Vl-Vd.gt.0.0_WP)
                     ! Make room for new drop
                     np=b3%lp%np_+1; call b3%lp%resize(np)
                     ! Add the drop
                     b3%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                     b3%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                     b3%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                     b3%lp%p(np)%pos =b2%vf%Lbary(:,i,j,k)                       !< Place the drop at the liquid barycenter
                     b3%lp%p(np)%vel =b2%fs%cfg%get_velocity(pos=b3%lp%p(np)%pos,i0=i,j0=j,k0=k,U=b2%fs%U,V=b2%fs%V,W=b2%fs%W) !< Interpolate local cell velocity as drop velocity
                     b3%lp%p(np)%ind =b3%lp%cfg%get_ijk_global(b3%lp%p(np)%pos,[b3%lp%cfg%imin,b3%lp%cfg%jmin,b3%lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                     b3%lp%p(np)%flag=0                                          !< Activate it
                     ! Increment particle counter
                     b3%lp%np_=np
                     ! Update tracked volumes
                     Vl=Vl-Vd
                     Vt=Vt+Vd
                  end do
                  ! Remove liquid in that cell
                  b2%vf%VF(i,j,k)=0.0_WP
               end do

               ! Based on how many particles were created, decide what to do with left-over volume
               if (Vt.eq.0.0_WP.and.Vl.gt.((pi*diam_thresh)**3)/6.0_WP) then ! No particle was created, we need one...
                  ! Add one last drop for remaining liquid volume
                  np=b3%lp%np_+1; call b3%lp%resize(np)
                  ! Add the drop
                  b3%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                  b3%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  b3%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
                  b3%lp%p(np)%pos =b2%vf%Lbary(:,i,j,k)                       !< Place the drop at the liquid barycenter
                  b3%lp%p(np)%vel =b2%fs%cfg%get_velocity(pos=b3%lp%p(np)%pos,i0=i,j0=j,k0=k,U=b2%fs%U,V=b2%fs%V,W=b2%fs%W) !< Interpolate local cell velocity as drop velocity
                  b3%lp%p(np)%ind =b3%lp%cfg%get_ijk_global(b3%lp%p(np)%pos,[b3%lp%cfg%imin,b3%lp%cfg%jmin,b3%lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                  b3%lp%p(np)%flag=0                                          !< Activate it
                  ! Increment particle counter
                  b3%lp%np_=np
               else ! Some particles were created, make them all larger
                  do ip=np_old+1,b3%lp%np_
                     b3%lp%p(ip)%d=b3%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
                  end do
               end if
            end do


         end block remove_film

         ! Sync VF and clean up IRL and band
         call b2%vf%cfg%sync(b2%vf%VF)
         call b2%vf%clean_irl_and_band()

         ! Clean up CCL
         call b2%cc2%deallocate_lists()

         ! Resync the spray
         call b3%lp%sync()

      end subroutine transfer_vf_to_drops


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      character(len=str_medium) :: timestamp

      write(timestamp,'(es12.5)') b2%time%t
      ! Populate df2 and write it
      call b2%df%pushval(name=  't',val=b2%time%t )
      call b2%df%pushval(name= 'dt',val=b2%time%dt)
      call b2%df%pushvar(name=  'U',var=b2%fs%U   )
      call b2%df%pushvar(name=  'V',var=b2%fs%V   )
      call b2%df%pushvar(name=  'W',var=b2%fs%W   )
      call b2%df%pushvar(name=  'P',var=b2%fs%P   )
      call b2%df%pushvar(name= 'LM',var=b2%sgs%LM )
      call b2%df%pushvar(name= 'MM',var=b2%sgs%MM )
      call b2%df%pushvar(name= 'VF',var=b2%vf%VF  )
      call b2%df%write(fdata='data2_'//trim(adjustl(timestamp)))
      ! Populate df3 and write it
      call b3%df%pushval(name=   't',val=b3%time%t )
      call b3%df%pushval(name=  'dt',val=b3%time%dt)
      call b3%df%pushvar(name=   'U',var=b3%fs%U   )
      call b3%df%pushvar(name=   'V',var=b3%fs%V   )
      call b3%df%pushvar(name=   'W',var=b3%fs%W   )
      call b3%df%pushvar(name=   'P',var=b3%fs%P   )
      call b3%df%pushvar(name=  'LM',var=b3%sgs%LM )
      call b3%df%pushvar(name=  'MM',var=b3%sgs%MM )
      call b3%df%write(fdata='data3_'//trim(adjustl(timestamp)))
      ! Also output particles
      call b3%lp%write(filename='datalpt_'//trim(adjustl(timestamp)))

      ! Deallocate work arrays for both blocks
      call b2%final()
      call b3%final()

   end subroutine simulation_final


end module simulation
