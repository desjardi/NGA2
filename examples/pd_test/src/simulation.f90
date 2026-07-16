!> Tier-1 pdsolver regression: grid-free solid dynamics, no AMReX objects.
!>
!> A tilted cube of nodes with an initial contraction velocity field
!> oscillates elastically. Exercises: pd_partition, detect_families, the
!> halo machinery, the LPS kernels, and checkpoint/restart (write_state at
!> 'Checkpoint time', restartable via 'Restart from' for a continuity check
!> against an uninterrupted run's monitor trace -- also at a different
!> rank count, which validates the gid-space re-partitioning).
module simulation
   use precision,         only: WP,I8
   use pdsolver_class,    only: pdsolver,pd_partition,PDC_MOVES,PDC_INTEGRATES,PDC_BONDS
   use timetracker_class, only: timetracker
   use monitor_class,     only: monitor
   use messager,          only: log
   use string,            only: str_medium
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final

   type(pdsolver) :: pd
   type(timetracker) :: time
   type(monitor) :: mfile
   real(WP) :: ckpt_time=-1.0_WP
   logical  :: ckpt_done=.false.
   character(len=str_medium) :: restart_dir=''

contains

   subroutine simulation_init()
      use param,    only: param_read
      use parallel, only: amRoot
      implicit none
      real(WP) :: rho,E,nu,elem,delta,cube_L,contract,tilt(3)
      integer :: cube_n

      time=timetracker(amRoot=amRoot)
      call param_read('Max time',time%tmax)
      call param_read('Max dt',  time%dtmax)
      time%dt=time%dtmax
      call param_read('Checkpoint time',ckpt_time,default=-1.0_WP)
      call param_read('Restart from',restart_dir,default='')

      call param_read('Material density',rho)
      call param_read('Elastic modulus', E)
      call param_read('Poisson ratio',   nu)
      call param_read('Element size',    elem)
      call param_read('Horizon',         delta,default=3.0125_WP*elem)
      ! Configure by field assignment (grid-free: no domain, no periodicity)
      pd%name='pd'
      pd%rho=rho; pd%elastic_modulus=E; pd%poisson_ratio=nu
      pd%delta=delta; pd%dV=elem**3

      if (len_trim(restart_dir).gt.0) then
         call pd%read_state(trim(restart_dir))
         if (amRoot) call log('[pd_test] restarted from '//trim(restart_dir))
      else
         seed: block
            use mathtools, only: Pi
            integer(I8), allocatable :: gids(:),rgid(:)
            real(WP), allocatable :: pos(:,:),vel(:,:),voll(:),rpos(:,:),rvel(:,:),rvol(:)
            integer, allocatable :: flags(:),owner(:),rflag(:)
            integer :: i,j,k,n,nn,nr
            real(WP) :: dx,x0,y0,z0,x1,y1,z1,x2,y2,z2,cx,sx,cy,sy,cz,sz
            call param_read('Cube size',cube_L,default=0.5_WP)
            call param_read('Cube tilt',tilt,default=[20.0_WP,30.0_WP,0.0_WP])
            call param_read('Contraction rate',contract,default=10.0_WP)
            cube_n=nint(cube_L/elem); cube_L=real(cube_n,WP)*elem
            tilt=tilt*Pi/180.0_WP
            cx=cos(tilt(1)); sx=sin(tilt(1)); cy=cos(tilt(2)); sy=sin(tilt(2))
            cz=cos(tilt(3)); sz=sin(tilt(3))
            dx=cube_L/real(cube_n,WP)
            ! Root builds the whole lattice; pd_partition routes it (gids are
            ! simply 1..n -- any unique positive keys work)
            nn=0
            if (amRoot) nn=cube_n**3
            allocate(gids(max(nn,1)),pos(3,max(nn,1)),vel(3,max(nn,1)),flags(max(nn,1)),voll(max(nn,1)),owner(max(nn,1)))
            n=0
            do k=0,cube_n-1; do j=0,cube_n-1; do i=0,cube_n-1
               if (.not.amRoot) exit
               n=n+1
               x0=-0.5_WP*cube_L+(real(i,WP)+0.5_WP)*dx
               y0=-0.5_WP*cube_L+(real(j,WP)+0.5_WP)*dx
               z0=-0.5_WP*cube_L+(real(k,WP)+0.5_WP)*dx
               x1=x0;          y1=cx*y0-sx*z0; z1=sx*y0+cx*z0
               x2=cy*x1+sy*z1; y2=y1;          z2=-sy*x1+cy*z1
               pos(:,n)=[cz*x2-sz*y2, sz*x2+cz*y2, z2]
               vel(:,n)=-contract*pos(:,n)
               gids(n)=int(n,I8)
               flags(n)=PDC_MOVES+PDC_INTEGRATES+PDC_BONDS
               voll(n)=elem**3
            end do; end do; end do
            call pd_partition(nn,gids,pos,vel,flags,voll,owner,nr,rgid,rpos,rvel,rflag,rvol)
            call pd%set_nodes(nr,rgid,rpos,rvel,rflag,rvol)
            call pd%detect_families()
         end block seed
      end if

      call pd%get_info()
      mfile=monitor(amRoot=amRoot,name='simulation')
      call mfile%add_column(time%n,'Timestep')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(pd%np,'Nodes')
      call mfile%add_column(pd%nb,'Bonds')
      call mfile%add_column(pd%nb_broken,'Bonds broken')
      call mfile%add_column(pd%Umax,'Umax')
      call mfile%add_column(pd%wtmax_force,'force_max')
      call mfile%write()
   end subroutine simulation_init

   subroutine simulation_run()
      use parallel, only: amRoot
      implicit none
      do while (.not.time%done())
         call time%increment()
         call pd%advance(time%dt)
         if (mod(time%n,10).eq.0) then
            call pd%get_info()
            call mfile%write()
         end if
         if (ckpt_time.gt.0.0_WP.and.time%t.ge.ckpt_time.and..not.ckpt_done) then
            save_state: block
               use parallel, only: comm,rank
               use mpi_f08,  only: MPI_BARRIER
               integer :: ierr
               if (rank.eq.0) call execute_command_line('mkdir -p restart/pd_test')
               call MPI_BARRIER(comm,ierr)
               call pd%write_state('restart/pd_test')
               ckpt_done=.true.
               if (amRoot) call log('[pd_test] checkpoint written at restart/pd_test')
            end block save_state
         end if
      end do
      if (amRoot) call log('[pd_test] run complete')
   end subroutine simulation_run

   subroutine simulation_final()
      implicit none
      call mfile%finalize()
      call pd%finalize()
      call time%finalize()
   end subroutine simulation_final

end module simulation
