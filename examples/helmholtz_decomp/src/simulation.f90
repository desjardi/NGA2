!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use incomp_class,      only: incomp
   use ensight_class,     only: ensight
   implicit none
   private

   !> Solvers
   type(fft2d),       public :: ps !< Poisson solver
   type(incomp),      public :: fs !< Flow solver

   !> Ensight postprocessing
   type(ensight) :: ens_out

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,Udi,Vdi,Wdi,Usi,Vsi,Wsi
   real(WP), dimension(:,:,:,:), allocatable :: vort,buf3D

contains


   !> Initialization of problem solver
   subroutine simulation_init
     use param, only: param_read
     use mathtools, only : Pi
      implicit none


      ! Allocate work arrays
      allocate_work_arrays: block
        allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Udi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vdi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wdi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Usi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vsi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wsi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(vort(3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(buf3d(3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays

      ! Initialize solver
      create_and_initialize_flow_solver: block
        integer :: i,j,k
        ! Create flow solver with default parameters
        fs=incomp(cfg=cfg,name='NS solver')
        fs%visc=0.0_WP; fs%rho=1.0_WP        
        ! Configure pressure solver
        ps=fft2d(cfg=cfg,name='Pressure',nst=7)
        call fs%setup(pressure_solver=ps)
        fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
        do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 fs%U(i,j,k)=sin(fs%cfg%ym(j)/fs%cfg%yL*2.0_WP*Pi)*sin(fs%cfg%x(i)/fs%cfg%xL*2.0_WP*Pi)
              end do
           end do
        end do
        ! Store cell-centered velocity
       call fs%interp_vel(Ui,Vi,Wi)
      end block create_and_initialize_flow_solver


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='helmholtz')
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('solenoidal',Usi,Vsi,Wsi)
         call ens_out%add_vector('dilatational',Udi,Vdi,Wdi)
       end block create_ensight

   end subroutine simulation_init


   !> The main solver
   subroutine simulation_run
     use mpi_f08,  only: MPI_MAX,MPI_ALLREDUCE
     use parallel, only: MPI_REAL_WP
     implicit none
     integer :: i,j,k,ierror
     real(WP) :: my_max,max_val

     ! Compute divergence and curl
     call fs%get_div()
     call fs%get_vorticity(vort)
     my_max=maxval(abs(fs%div(cfg%imin_+1:cfg%imax_-1,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)))
     call MPI_ALLREDUCE(my_max,max_val,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierror)
     if (cfg%amRoot) print *, 'Max div(U)=',max_val
     
     ! Solve Poisson equation for velocity potential
     fs%psolv%rhs=-fs%cfg%vol*fs%div
     fs%psolv%sol=0.0_WP
     call fs%psolv%solve()

     ! Get dilataiotnal velocity: Ud=grad(phi)
     call fs%get_pgrad(fs%psolv%sol,fs%U,fs%V,fs%W)
     call fs%interp_vel(Udi,Vdi,Wdi)

     ! Verify Ud is irrotational
     call fs%get_vorticity(buf3d)
     my_max=maxval(abs(buf3d(:,cfg%imin_+1:cfg%imax_-1,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)))
     call MPI_ALLREDUCE(my_max,max_val,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierror)
     if (cfg%amRoot) print *, 'Max Curl(Ud)=',max_val

     ! Solve Poisson equation for 'A' in each direction
     do i=1,3
        ps%rhs=fs%cfg%vol*vort(i,:,:,:)
        ps%sol=0.0_WP
        call ps%solve()
        buf3d(i,:,:,:)=ps%sol
     end do

     ! Interpolate to cell faces
     do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
           do i=fs%cfg%imin_,fs%cfg%imax_
              fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*buf3d(1,i-1:i,j,k))
              fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*buf3d(2,i,j-1:j,k))
              fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*buf3d(3,i,j,k-1:k))
           end do
        end do
     end do
     call fs%cfg%sync(fs%U);call fs%cfg%sync(fs%V); call fs%cfg%sync(fs%W)

     ! Get solenoidal velocity: Us=curl(A)
     call fs%get_vorticity(buf3d)
     Usi=buf3d(1,:,:,:); Vsi=buf3d(2,:,:,:); Wsi=buf3d(3,:,:,:)

     ! Interpolate to cell faces
     do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
           do i=fs%cfg%imin_,fs%cfg%imax_
              fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*Usi(i-1:i,j,k))
              fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*Vsi(i,j-1:j,k))
              fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*Wsi(i,j,k-1:k))
           end do
        end do
     end do
     call fs%cfg%sync(fs%U);call fs%cfg%sync(fs%V); call fs%cfg%sync(fs%W)

     ! Verify Us is solenoidal
     call fs%get_div()
     my_max=maxval(abs(fs%div(cfg%imin_+1:cfg%imax_-1,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)))
     call MPI_ALLREDUCE(my_max,max_val,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierror)
     if (cfg%amRoot) print *, 'Max div(Us)=',max_val

     ! Output to ensight
     call ens_out%write_data(0.0_WP)

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Deallocate work arrays
      deallocate(Ui,Vi,Wi,Udi,Vdi,Wdi,Usi,Vsi,Wsi,vort,buf3d)

    end subroutine simulation_final

end module simulation
