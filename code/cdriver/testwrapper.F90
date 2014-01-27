
#include "cdriverconstants.h"

program fwrapper
  ! Want to know how to run BISICLES from FORTRAN 90? Avoid learning
  ! C++ ? Follow the example below, but hang your head in shame
  
  implicit none
  character(len=25) :: file

  integer instance_id,  it, nt, max_step
  integer, dimension(1:2) :: dims, boxlo, boxhi
  real(kind=8) :: dx, max_time 
  integer, parameter :: nx = 64, ny = 96
  real(kind=8), dimension(1:nx,1:ny) :: smb, bmbf, bmbg, usrf

  !need to leave at least one character beyond the filename to null-terminate on the C-side
  file = "inputs.pigv5.1km.l1l2.l1" 
  
  !create an instance
  call bisicles_new_instance(instance_id, file, 25)
  
  !now set up some uniform mesh data to read from / write to the interface
  !mesh spacing
  dx = 4.0e+3
  !data dimensions
  dims(1) = nx
  dims(2) = ny
  !lower left corner on the BISICLES level with resolution dx
  boxlo(1) = 0
  boxlo(2) = 0
  !top right corner
  boxhi(1) = dims(1) - 1
  boxhi(2) = dims(2) - 1


  !tell BISICLES to read a surface flux from smb, and a basal fluxes from bmbf and bmbg
  smb = 0.0d0
  call bisicles_set_2d_data(instance_id, smb, BISICLES_FIELD_SURFACE_FLUX, dx, dims, boxlo, boxhi)
  bmbf = -100.0d0
  call bisicles_set_2d_data(instance_id, bmbf, BISICLES_FIELD_FLOATING_ICE_BASAL_FLUX, dx, dims, boxlo, boxhi)


  !At this point, we have given BISICLES as much data as it needs to run. So, initialize it here, at which point the AMR velocity problem gets solved (or a checkpoint loaded)
  !After this step we can still change the data in (say) smb but we cannot change its location
  call bisicles_init_instance(instance_id)

  !now do some time stepping
  nt = 3
  max_step = 0
  max_time = 0.0
  do it = 1, nt
     !read the surface elevation
     call bisicles_get_2d_data(instance_id, usrf, BISICLES_FIELD_SURFACE_ELEVATION, dx, dims, boxlo, boxhi)

     !re-compute the smb - just change the data
     smb = smb + 1.0d0
     
     !advance in time
     max_time = max_time + 1.0d0; !advance by one year
     max_step = max_step + 10; !unless it takes to long, in which case give up
     call bisicles_advance(instance_id, max_time, max_step)
     
  end do

  !free any memory allocated on the C++ side
  call bisicles_free_instance(instance_id)

end program fwrapper
