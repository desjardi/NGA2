!> Various definitions and tools for running an NGA2 simulation
module simulation
   implicit none
   public :: simulation_init,simulation_run,simulation_final
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use precision, only: WP
      use param,     only: param_read
      use string,    only: str_medium,str_long
      use lpt_class, only: prepare_mpi_part,part,MPI_PART,MPI_PART_SIZE
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      character(len=str_long) :: message
      type(part), dimension(:), allocatable :: p_in,p_out
      integer :: npart
      
      
      ! Prepare MPI_PART datatype
      call prepare_mpi_part()
      
      
      ! Read in the particle file
      check_serial: block
         use parallel, only: nproc
         use messager, only: die
         if (nproc.ne.1) call die('[part_sorter] This code assumes serial execution')
      end block check_serial
      
      
      ! Read in the particle file
      read_particles: block
         use mpi_f08
         use messager, only: die
         use parallel, only: comm,info_mpiio
         character(len=str_medium) :: file_in
         type(MPI_File) :: ifile
         type(MPI_Status):: status
         integer :: ierr,psize
         integer(kind=MPI_OFFSET_KIND) :: offset
         ! Read in the name of the particle file to read
         call param_read('Part file to read',file_in,short='r')
         ! Read the header of the particle file
         call MPI_FILE_OPEN(comm,trim(file_in),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[part_sorter] Problem encountered while reading particle file: '//trim(file_in))
         call MPI_FILE_READ_ALL(ifile,npart,1,MPI_INTEGER,status,ierr)
         call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
         ! Allocate storage and check compatibility of particle type
         allocate(p_in(npart))
         if (psize.ne.MPI_PART_SIZE) call die('[part_sorter] Particle type unreadable')
         ! Read the rest of the particle file
         call MPI_FILE_GET_POSITION(ifile,offset,ierr)
         call MPI_FILE_READ_AT(ifile,offset,p_in,npart,MPI_PART,status,ierr)
         call MPI_FILE_CLOSE(ifile,ierr)
         ! Some output
         write(message,'("Particle sorter: ",i0," particles were read")') npart
         write(output_unit,'(a)') trim(message)
      end block read_particles
      
      
      ! Sort particles
      sort_particles: block
         use quicksort, only: quick_sort
         integer :: i
         real(WP), dimension(:), allocatable :: pos
         integer , dimension(:), allocatable :: id
         ! Prepare sorting data
         allocate(pos(npart),id(npart))
         do i=1,npart
            id(i)=i
            pos(i)=p_in(i)%pos(2)
         end do
         ! Apply quicksort
         call quick_sort(A=pos,B=id)
         ! Loop through particles and sort them
         allocate(p_out(npart))
         do i=1,npart
            p_out(i)=p_in(id(i))
            p_out(i)%id=i
         end do
         ! Some output
         write(message,'("Particle sorter: ",i0," particles were sorted")') npart
         write(output_unit,'(a)') trim(message)
      end block sort_particles
      
      
      ! Write out the particles
      write_particles: block
         use mpi_f08
         use parallel, only: comm,info_mpiio
         use messager, only: die
         use string,   only: str_medium
         character(len=str_medium) :: file_out
         type(MPI_File) :: ifile
         type(MPI_Status):: status
         integer :: ierr,iunit
         integer(kind=MPI_OFFSET_KIND) :: offset
         ! Read in the name of the particle file to read
         call param_read('Part file to write',file_out,short='w')
         ! Write the header of the particle file
         open(newunit=iunit,file=trim(file_out),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[part_sorter] Problem encountered while writing particle file: '//trim(file_out))
         write(iunit) npart,MPI_PART_SIZE
         close(iunit)
         ! Write the rest of the particle file
         call MPI_FILE_OPEN(comm,trim(file_out),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[part_sorter] Problem encountered while writing particle file: '//trim(file_out))
         call MPI_FILE_GET_POSITION(ifile,offset,ierr)
         call MPI_FILE_WRITE_AT(ifile,offset,p_out,npart,MPI_PART,status,ierr)
         call MPI_FILE_CLOSE(ifile,ierr)
         ! Some output
         write(message,'("Particle sorter: ",i0," particles were written")') npart
         write(output_unit,'(a)') trim(message)
      end block write_particles
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
   end subroutine simulation_final
   
   
end module simulation
