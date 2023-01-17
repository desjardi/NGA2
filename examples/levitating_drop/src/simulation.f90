!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs, VFhi, VFlo
   use matm_class,        only: matm
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single two-phase flow solver, volume fraction solver, and material model set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(timetracker), public :: time

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,dfile,cflfile,consfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Choice of relaxation model
   integer  :: relax_model
   !> Problem definition
   real(WP), dimension(3) :: dctr
   real(WP) :: ddrop, Pref, lam_wave, Grho0, gamm_g, Pmax, Ls, f_wave, c0
   
   !> Monitor input variables
   real(WP), dimension(3) :: prob_loc = (/0.0_WP,0.25_WP,0.0_WP/)
   real(WP) :: prob_diam = 1.0e-3_WP
   !> Monitor output variables
   real(WP) :: r_drop, Prms_probe, y_drop, Fx_avg, Fy_avg, Fz_avg
   !> Monitor working variables
   integer :: n_sample



contains

   !> Function that localizes the top (y+) of the domain
   function top_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain

   !> Function that localizes the bottom (y-) of the domain
   function btm_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.eq.pg%jmin-1) isIn=.true.
   end function btm_of_domain

   !> Function that localizes top sponge
   function top_sponge(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (pg%y(pg%jmax+1)-pg%ym(j).le.Ls) isIn=.true.
   end function top_sponge

   !> Function that localizes bottom sponge
   function btm_sponge(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (pg%ym(j)-pg%y(pg%jmin).le.Ls) isIn=.true.
   end function btm_sponge

   subroutine apply_sponges(t)
     use mathtools,  only: Pi
     implicit none
     real(WP), intent(in) :: t
     integer :: i,j,k
     logical :: in_sponge
     real(WP) :: swt, ps, rhos, us, vs, ws

     do k=fs%cfg%kmin_,fs%cfg%kmax_
       do j=fs%cfg%jmin_,fs%cfg%jmax_
         do i=fs%cfg%imin_,fs%cfg%imax_
           in_sponge = .false.
           ! Check top sponge
           if (top_sponge(fs%cfg,i,j,k)) then
             in_sponge = .true.
             ! Get sponge coefficient
             swt = min(1.0_WP,max(0.0_WP,1.0_WP/Ls*(fs%cfg%ym(j)-Ls)+0.5_WP))**2
           end if
           ! Check bottom sponge
           if (btm_sponge(fs%cfg,i,j,k)) then
             in_sponge = .true.
             ! Get sponge coefficient
             swt = min(1.0_WP,max(0.0_WP,1.0_WP/Ls*(-fs%cfg%ym(j)-Ls)+0.5_WP))**2
           end if
           if (in_sponge) then
             ! Get sponge solution if within a sponge
             ps = Pref + Pmax*sin(2.0_WP*pi*fs%cfg%ym(j)/(c0/f_wave))*cos(2.0_WP*pi*f_wave*t)
             rhos = Grho0 + (ps-Pref)/(gamm_g*ps/Grho0)
             us = 0.0_WP; ws = 0.0_WP
             vs = 0.0_WP - Pmax/sqrt(rhos*gamm_g*ps)*cos(2.0_WP*pi*fs%cfg%ym(j)/(c0/f_wave))*sin(2.0_WP*pi*f_wave*t)
             ! Apply changes to variables
             fs%Grho (i,j,k) = fs%Grho (i,j,k) - swt*(fs%Grho(i,j,k)-rhos)
             fs%Ui   (i,j,k) = fs%Ui   (i,j,k) - swt*(fs%Ui(i,j,k)-us)
             fs%Vi   (i,j,k) = fs%Vi   (i,j,k) - swt*(fs%Vi(i,j,k)-vs)
             fs%Wi   (i,j,k) = fs%Wi   (i,j,k) - swt*(fs%Wi(i,j,k)-ws)
             fs%GP   (i,j,k) = fs%GP   (i,j,k) - swt*(fs%GP(i,j,k)-ps)
             fs%GrhoE(i,j,k) = fs%GrhoE(i,j,k) - swt*(fs%GrhoE(i,j,k)-matmod%EOS_energy(ps,rhos,us,vs,ws,'gas'))
             ! Update related quantities
             fs%RHO  (i,j,k) = fs%Grho(i,j,k)
             fs%rhoUi(i,j,k) = fs%RHO(i,j,k)*fs%Ui(i,j,k)
             fs%rhoVi(i,j,k) = fs%RHO(i,j,k)*fs%Vi(i,j,k)
             fs%rhoWi(i,j,k) = fs%RHO(i,j,k)*fs%Wi(i,j,k)
             ! Remove liquid
             vf%VF   (i,j,k) = 0.0_WP
             fs%Lrho (i,j,k) = 0.0_WP
             fs%LP   (i,j,k) = 0.0_WP
             fs%LrhoE(i,j,k) = 0.0_WP
           end if

         end do
       end do
     end do

     ! Reset interface-based quantities due to VF changes
     call vf%remove_flotsams()
     call vf%advect_interface(0.0_WP,fs%U,fs%V,fs%W)
     call vf%build_interface()
     call vf%remove_sheets()
     call vf%polygonalize_interface()
     call vf%subcell_vol()
     call vf%get_curvature()
     call vf%reset_moments()

   end subroutine
  
   !> Updates quantities specific to the monitor of this case
   subroutine lev_drop_update()
     use irl_fortran_interface, only: calculateNormal, calculateVolume, getNumberOfPlanes, getNumberOfVertices, calculateCentroid
     use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX
     use parallel,  only: MPI_REAL_WP
     implicit none
     real(WP) :: buf,Lvol_buf,p_buf,A_buf,A_tmp,G,r_eq,pgas,area,x_drop,z_drop
     real(WP), dimension(3) :: pctr, norm, f_vec
     integer, parameter :: nF = 20
     integer  :: ierr,i,j,k,ii,jj,kk,nplane,n
     
     ! Initialize
     Lvol_buf = 0.0_WP; r_drop = 0.0_WP
     p_buf = 0.0_WP; A_buf = 0.0_WP
     x_drop = 0.0_WP; y_drop = 0.0_WP; z_drop = 0.0_WP
     f_vec = 0.0_WP

     ! Increment number of samples
     n_sample = n_sample+1
     do k=fs%cfg%kmin_,fs%cfg%kmax_
       do j=fs%cfg%jmin_,fs%cfg%jmax_
         do i=fs%cfg%imin_,fs%cfg%imax_

           ! Calculate the mean pressure measurement at the probe for the current timestep
           if (fs%cfg%y(j+1).ge.prob_loc(2).and.fs%cfg%y(j).le.prob_loc(2)) then ! check cells in y
              ! get interpolation quantities for y
              if (fs%cfg%ym(j).gt.prob_loc(2)) then
                 jj = j-1
              else
                 jj = j
              end if
              if (fs%cfg%nz.eq.1) then
                 ! Check if x is within specs
                 if (fs%cfg%x(i).le.prob_loc(1)+prob_diam/2.0_WP.and.fs%cfg%x(i+1).ge.prob_loc(1)-prob_diam/2.0_WP) then
                    A_tmp = fs%cfg%dz(k)*(min(fs%cfg%x(i+1),prob_loc(1)+prob_diam/2.0_WP) &
                          - max(fs%cfg%x(i),prob_loc(1)-prob_diam/2.0_WP))
                    p_buf = p_buf+A_tmp*(fs%P(i,jj,k)-Pref+(fs%P(i,jj+1,k)-fs%P(i,jj,k)) / &
                            (fs%cfg%ym(jj+1)-fs%cfg%ym(jj))*(prob_loc(2)-fs%cfg%ym(jj)))
                    A_buf = A_buf+A_tmp
                 end if
              else
                 ! Check if x,z are near enough to or within area
                 G = sqrt((fs%cfg%xm(i)-prob_loc(1))**2+(fs%cfg%zm(k)-prob_loc(3))**2)-prob_diam/2.0_WP
                 if (G.lt.2.0_WP*max(fs%cfg%dx(i),fs%cfg%dz(k))) then
                    ! Create inner fine grid to get area within probe
                    A_tmp = 0.0_WP
                    do ii=1,nF
                       do kk=1,nF
                          G = sqrt((fs%cfg%x(i)+fs%cfg%dx(i)/real(nF,WP)*(real(ii,WP)-0.5_WP)-prob_loc(1))**2 &
                               +(fs%cfg%z(k)+fs%cfg%dz(k)/real(nF,WP)*(real(kk,WP)-0.5_WP)-prob_loc(3))**2) &
                               -prob_diam/2.0_WP
                          if (G.le.0.0_WP) then
                             A_tmp = A_tmp+fs%cfg%dx(i)*fs%cfg%dz(k)/real(nF,WP)**2
                          end if
                       end do
                    end do
                    ! Accumulate pressure and area
                    p_buf = p_buf+A_tmp*(fs%P(i,jj,k)-Pref+(fs%P(i,jj+1,k)-fs%P(i,jj,k)) / &
                            (fs%cfg%ym(jj+1)-fs%cfg%ym(jj))*(prob_loc(2)-fs%cfg%ym(jj)))
                    A_buf = A_buf+A_tmp
                 end if
              end if
           end if

           ! Calculate the location of drop center of gravity
           Lvol_buf = Lvol_buf + fs%cfg%vol(i,j,k)*vf%VF(i,j,k)
           if (vf%VF(i,j,k).gt.VFhi) then
              ! add entire volume of cell to calculation of center of mass
              x_drop = x_drop + fs%cfg%vol(i,j,k)*fs%cfg%xm(i)
              y_drop = y_drop + fs%cfg%vol(i,j,k)*fs%cfg%ym(j)
              z_drop = z_drop + fs%cfg%vol(i,j,k)*fs%cfg%zm(k)
           elseif (vf%VF(i,j,k).ge.VFlo) then
              ! add liquid volume multiplied by liquid barycenter location
              x_drop = x_drop + fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*vf%Lbary(1,i,j,k)
              y_drop = y_drop + fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*vf%Lbary(2,i,j,k)
              z_drop = z_drop + fs%cfg%vol(i,j,k)*vf%VF(i,j,k)*vf%Lbary(3,i,j,k)

              ! Calculate the radiation force by integrating over the droplet surface
              ! get local gas pressure
              pgas = fs%P(i,j,k) - vf%VF(i,j,k)*fs%hpjump(i,j,k)
              ! multiply by local plane area vector, add to force
              nplane = getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
              do n=1,nplane
                 if (getNumberOfVertices(vf%interface_polygon(n,i,j,k)).gt.0) then
                    norm = calculateNormal(vf%interface_polygon(n,i,j,k))
                    area = abs(calculateVolume(vf%interface_polygon(n,i,j,k)))
                    f_vec = f_vec - pgas*norm*area
                    ! negative because vector direction should be opposite of outward normal
                 end if
              end do
           end if
         end do
       end do
     end do

     ! Get x and z centers for radius calculation
     call MPI_ALLREDUCE(Lvol_buf,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Lvol_buf = buf
     call MPI_ALLREDUCE(x_drop,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); x_drop = buf/(Lvol_buf+tiny(1.0_WP))
     call MPI_ALLREDUCE(z_drop,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); z_drop = buf/(Lvol_buf+tiny(1.0_WP))

     do k=fs%cfg%kmin_,fs%cfg%kmax_
       do j=fs%cfg%jmin_,fs%cfg%jmax_
         do i=fs%cfg%imin_,fs%cfg%imax_
           if ((vf%VF(i,j,k).ge.VFlo).and.(vf%VF(i,j,k).le.VFhi)) then
              ! Get plane information
              nplane = getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
              do n=1,nplane
                 if (getNumberOfVertices(vf%interface_polygon(n,i,j,k)).gt.0) then
                   ! Get intersection of plane with line from center
                   pctr = calculateCentroid(vf%interface_polygon(n,i,j,k))
                   ! Calculate radius in xz plane
                   r_eq = sqrt((pctr(1)-x_drop)**2+(pctr(3)-z_drop)**2)
                   r_drop = max(r_drop,r_eq)
                 end if
              end do
           end if
         end do
       end do
     end do
     !> Outputs: <!
     ! Normalized radius, RMS pressure at probe, droplet y position, time-averaged force

     ! Get average pressure on microphone (probe) surface for current time
     call MPI_ALLREDUCE(A_buf,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); A_buf = buf
     call MPI_ALLREDUCE(p_buf,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); p_buf = buf
     p_buf = p_buf/A_buf
     ! Contribute to RMS value over time
     Prms_probe = sqrt((real(n_sample-1,WP)*Prms_probe**2+p_buf**2)/real(n_sample,WP))

     ! Get current radius of droplet
     call MPI_ALLREDUCE(r_drop,buf,1,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr); r_drop = buf;
     ! Normalize by the initial radius
     r_drop = 2.0_WP*r_drop/ddrop
     
     ! Get mean y position of droplet
     call MPI_ALLREDUCE(y_drop,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); y_drop = buf/(Lvol_buf+tiny(1.0_WP))

     ! Sum forces across processors
     call MPI_ALLREDUCE(f_vec(1),buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); f_vec(1) = buf
     call MPI_ALLREDUCE(f_vec(2),buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); f_vec(2) = buf
     call MPI_ALLREDUCE(f_vec(3),buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); f_vec(3) = buf
     ! Contribute to avg value over time (assumes constant dt)
     Fx_avg = (real(n_sample-1,WP)*Fx_avg+f_vec(1))/real(n_sample,WP)
     Fy_avg = (real(n_sample-1,WP)*Fy_avg+f_vec(2))/real(n_sample,WP)
     Fz_avg = (real(n_sample-1,WP)*Fz_avg+f_vec(3))/real(n_sample,WP)
     
   end subroutine lev_drop_update
   
   !> Function that defines a level set function for a spherical droplet at center
   function levelset_drop_center(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP-sqrt((xyz(1)-dctr(1))**2+(xyz(2)-dctr(2))**2+(xyz(3)-dctr(3))**2)/(ddrop/2.0)
   end function levelset_drop_center
   
   function levelset_cyl_center(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP-sqrt((xyz(1)-dctr(1))**2+(xyz(2)-dctr(2))**2)/(ddrop/2.0)
   end function levelset_cyl_center

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Max steps',time%nmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker


      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,elvira,VFhi,VFlo
         use parallel,    only: amRoot
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
         ! Droplet information
         call param_read('Droplet diameter',ddrop)
         call param_read('Droplet relative location',dctr)
         ! Set the energy based on background pressure
         call param_read('Background pressure', Pref)
         ! Calculate the corresponding speed of sound
         call param_read('Gas gamma',gamm_g)
         call param_read('Gas density',Grho0)
         c0 = sqrt(gamm_g*Pref/Grho0)
         ! Get wave frequency
         call param_read('Wave frequency',f_wave)
         ! Calculate wavelength
         lam_wave = c0/f_wave
         if (amRoot) print*, 'standing wavelength', lam_wave
         ! Reference center scaled by wavelength
         dctr = dctr*lam_wave
         ! Droplet always initialized at center
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  if (vf%cfg%nz.eq.1) then
                    call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cyl_center,0.0_WP,amr_ref_lvl)
                  else
                    call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_drop_center,0.0_WP,amr_ref_lvl)
                  end if
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                    vf%Lbary(:,i,j,k)=v_cent
                    vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                    if (vf%cfg%nz.eq.1) vf%Gbary(3,i,j,k)=v_cent(3);
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a compressible two-phase flow solver
      create_and_initialize_flow_solver: block
         use mast_class, only: clipped_neumann,bc_scope,bcond,thermmech_egy_mech_hhz,mech_egy_mech_hhz
         use ils_class,  only: pcg_bbox,gmres_smg
         use mathtools,  only: Pi
         use parallel,   only: amRoot
         use string,     only: str_medium
         integer :: i,j,k
         real(WP) :: gamm_l,Pref_l,visc_l,visc_g,hdff_l,hdff_g,cv_l,cv_g,b_l,q_l
         real(WP) :: LP0, Lrho0, ps, rhos
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get EOS parameters from input
         call param_read('Liquid gamma',gamm_l)
         call param_read('Liquid Pref', Pref_l)
         call param_read('Liquid b', b_l)
         call param_read('Liquid q', q_l)
         ! Register equations of state
         call matmod%register_NobleAbelstiffenedgas('liquid',gamm_l,Pref_l,q_l,b_l)
         call matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',visc_l)
         call param_read('Gas dynamic viscosity',visc_g)
         ! Set heat diffusion to 0 unless specified
         hdff_l = 0.0_WP; hdff_g = 0.0_WP
         if (param_exists('Liquid therm cond')) then
           call param_read('Liquid therm cond', hdff_l)
           call param_read('Gas therm cond', hdff_g)
         end if
         ! Set specific heat to default
         cv_l = matmod%cv_l0; cv_g = matmod%cv_g0
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g,&
            viscconst_liquid=visc_l, hdffconst_gas=hdff_g, hdffconst_liquid=hdff_l, &
            sphtconst_gas=cv_g, sphtconst_liquid=cv_l)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_bbox,implicit_ils=gmres_smg)

         ! Calculate liquid pressure
         if (fs%cfg%nz.eq.1) then
           ! Cylinder configuration, curv = 1/r
           LP0 = Pref + 2.0/ddrop*fs%sigma
         else
           ! Sphere configuration, curv = 1/r + 1/r
           LP0 = Pref + 4.0/ddrop*fs%sigma
         end if
         ! Set liquid density
         call param_read('Liquid density',Lrho0)
         fs%Lrho = Lrho0
         ! Set liquid energy
         fs%LrhoE = matmod%EOS_energy(LP0,Lrho0,0.0_WP,0.0_WP,0.0_WP,'liquid')

         ! Get sponge parameters
         call param_read('Sponge max pressure',Pmax)
         call param_read('Sponge length',Ls)

         ! Initialize velocity and gas phase
         fs%Ui = 0.0_WP; fs%Vi = 0.0_WP; fs%Wi = 0.0_WP
         fs%U = 0.0_WP; fs%V = 0.0_WP; fs%W = 0.0_WP
         do j=vf%cfg%jmino_,vf%cfg%jmaxo_
           ! standing wave, node at center
           ps = Pref + Pmax*sin(2.0_WP*pi*fs%cfg%ym(j)/lam_wave)
           rhos = Grho0 + (ps-Pref)/(gamm_g*ps/Grho0)
           ! store values
           fs%Grho (:,j,:) = rhos
           fs%GrhoE(:,j,:) = matmod%EOS_energy(ps,rhos,0.0_WP,0.0_WP,0.0_WP,'gas')
         end do

         ! Define BCs at top and bottom - though sponges will determine bdy behavior
         call fs%add_bcond(name='top_y'   ,type=clipped_neumann,locator=top_of_domain  ,celldir='yp')
         call fs%add_bcond(name='btm_y'   ,type=clipped_neumann,locator=btm_of_domain  ,celldir='ym')
         
         ! Choose relaxation model
         if (abs(cv_l)+abs(cv_g).gt.0.0_WP) then
           ! Use thermo-mechanical if heat conduction is considered
           relax_model = thermmech_egy_mech_hhz
         else
           ! Use mechanical otherwise
           relax_model = mech_egy_mech_hhz
         end if

         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         ! Perform initial pressure relax
         call fs%pressure_relax(vf,matmod,relax_model)
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA
         ! Calculate initial temperature
         call matmod%update_temperature(vf,fs%Tmptr)
         ! Initialize first guess for pressure (0 works best)
         fs%psolv%sol=0.0_WP

      end block create_and_initialize_flow_solver


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='lev_drop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('Peos',fs%PA)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_scalar('Temperature',fs%Tmptr)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Tmax,'Tmax')
         call mfile%write()
         ! Check input file for monitor specs
         if (param_exists('Probe relative location')) then
            call param_read('Probe relative location',prob_loc)
         end if
         if (param_exists('Probe diameter')) then
            call param_read('Probe diameter',prob_diam)
         end if
         prob_loc = lam_wave * prob_loc ! scale by wavelength
         ! Initialize monoitor variables
         Prms_probe = 0.0_WP
         Fx_avg = 0.0_WP; Fy_avg = 0.0_WP; Fz_avg = 0.0_WP
         n_sample = 0
         ! Create levitating drop monitor
         ! Normalized radius, RMS pressure at probe, droplet y position, time-averaged force
         dfile=monitor(fs%cfg%amRoot,'levdrop')
         call dfile%add_column(time%n,'Timestep number')
         call dfile%add_column(time%t,'Time')
         call dfile%add_column(r_drop,'R^*')
         call dfile%add_column(Prms_probe,'Prms_p')
         call dfile%add_column(y_drop,'y_drop')
         call dfile%add_column(Fx_avg,'Fx_avg')
         call dfile%add_column(Fy_avg,'Fy_avg')
         call dfile%add_column(Fz_avg,'Fz_avg')
         call lev_drop_update()
         call dfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'cvg')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(fs%impl_it_x,'Impl_x iteration')
         call cvgfile%add_column(fs%impl_rerr_x,'Impl_x error')
         call cvgfile%add_column(fs%impl_it_y,'Impl_y iteration')
         call cvgfile%add_column(fs%impl_rerr_y,'Impl_y error')
         call cvgfile%add_column(fs%implicit%it,'Impl_z iteration')
         call cvgfile%add_column(fs%implicit%rerr,'Impl_z error')
         call cvgfile%add_column(fs%psolv%it,'Pressure iteration')
         call cvgfile%add_column(fs%psolv%rerr,'Pressure error')
         
      end block create_monitor


   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()

         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Diffusion and built-in source term (gravity) step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

            ! Perform sponge forcing
            call apply_sponges(time%t)

            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)
            
            ! Record convergence monitor
            call cvgfile%write()

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         fs%PA = matmod%EOS_all(vf);
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call lev_drop_update()
         call dfile%write()
         call cflfile%write()

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
     ! Nothing to do
   end subroutine simulation_final

end module simulation
