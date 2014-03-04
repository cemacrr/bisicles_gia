
#include "cdriverconstants.h"

program fwrapper
  ! Want to know how to run BISICLES from FORTRAN 90? Avoid learning
  ! C++ ? Follow the example below, but hang your head in shame
#ifdef CH_MPI
  use mpi
#endif	
  implicit none
  character(len=25) :: file

  integer instance_id,  it, nt, max_step
  integer, dimension(1:2) :: dims, boxlo, boxhi
  real(kind=8) :: dx, max_time 
  integer, parameter :: nx = 64, ny = 96
  real(kind=8), dimension(:,:), allocatable :: smb, bmbf, bmbg, usrf, seb
  integer ixlo,ixhi,iylo,iyhi
  
  integer :: rank,nproc,ierr

#ifdef CH_MPI
  call MPI_Init ( ierr )
  call MPI_Comm_rank ( mpi_comm_world, rank, ierr )
  call MPI_Comm_size ( mpi_comm_world, nproc, ierr )
#else
  rank = 0
  nproc = 1
#endif

  !example domain decompistion scheme : in serial, one block of data and in parallel two blocks
  !owned by ranks 0,1, with all the other processors contributing nothing 

  if (nproc.eq.1) then
     ixlo = 0
     ixhi = nx - 1
     iylo = 0
     iyhi = ny - 1
  else 
     if (rank.eq.0) then
         ixlo = 0
         ixhi = nx/2 - 1
         iylo = 0
         iyhi = ny - 1
      else if (rank.eq.1) then
         ixlo = nx/2
         ixhi = nx - 1
         iylo = 0
         iyhi = ny - 1
     else
        !no data : specify a block outside the ice sheet domain, which bisicles then ignores
        ixlo = nx + 1
        iylo = ny + 1
        ixhi = nx + 1
        iyhi = ny + 1
     end if

  end if	

  !if (rank.le.1) then
     allocate (smb(ixlo:ixhi,iylo:iyhi))
     allocate (seb(ixlo:ixhi,iylo:iyhi))
     allocate (usrf(ixlo:ixhi,iylo:iyhi))
     allocate (bmbf(ixlo:ixhi,iylo:iyhi))
     allocate (bmbg(ixlo:ixhi,iylo:iyhi))
  !end if

  !need to leave at least one character beyond the filename to null-terminate on the C-side
  file = "inputs.pigv5.1km.l1l2.l1" 
  
  !create an instance
  call f_bisicles_new_instance(instance_id, file, 25)
 
  !if (rank.le.1) then
     !now set up some uniform mesh data to read from / write to the interface
     !mesh spacing
     dx = 4.0e+3
     !data dimensions
     dims(1) = nx
     dims(2) = ny
     !lower left corner on the BISICLES level with resolution dx
     boxlo(1) = ixlo
     boxlo(2) = iylo
     !top right corner
     boxhi(1) = ixhi
     boxhi(2) = iyhi
     
     
     !tell BISICLES to read a surface flux from smb, and a basal fluxes from bmbf and bmbg
     smb = 1.0d0 +  dble(mod(rank,2))
     call f_bisicles_set_2d_data(instance_id, smb, BISICLES_FIELD_SURFACE_FLUX, dx, dims, boxlo, boxhi)
     bmbf = -0.0d0
     call f_bisicles_set_2d_data(instance_id, bmbf, BISICLES_FIELD_FLOATING_ICE_BASAL_FLUX, dx, dims, boxlo, boxhi)
     seb = 1.0d+7
     call f_bisicles_set_2d_data(instance_id, seb, BISICLES_FIELD_SURFACE_HEAT_FLUX, dx, dims, boxlo, boxhi)
  !end if

  !At this point, we have given BISICLES as much data as it needs to run. So, initialize it here, at which point the AMR velocity problem gets solved (or a checkpoint loaded)
  !After this step we can still change the data in (say) smb but we cannot change its location
  call f_bisicles_init_instance(instance_id)

  !now do some time stepping
  nt = 3
  max_step = 0
  max_time = 0.0
  do it = 1, nt
     !if (rank.le.1) then
        !read the surface elevation
        call f_bisicles_get_2d_data(instance_id, usrf, BISICLES_FIELD_SURFACE_ELEVATION, dx, dims, boxlo, boxhi)
        
        !re-compute the smb - just change the data. this is obviously a weird field that is supposed 
        !to show the domain decomposed reads/writes are working (look at surfaceThicknessBalance in the output)
        smb = usrf / 1000.0 +  dble(mod(rank,2))
     !end if
     !advance in time
     max_time = max_time + 1.0d0; !advance by one year
     max_step = max_step + 10; !unless it takes to long, in which case give up
     call f_bisicles_advance(instance_id, max_time, max_step)
     
  end do

  !free any memory allocated on the C++ side
  call f_bisicles_free_instance(instance_id)
#ifdef CH_MPI
  call MPI_Finalize ( ierr )
#endif

end program fwrapper
