#include <iostream>
#include "cdriverconstants.h"
#include "cwrapper.H"

int main(int argc, char* argv[])
{
  char fname[] = "inputs.pigv5.1km.l1l2.l1";
  int instance_id;

  bisicles_new_instance(&instance_id, fname);

  std::cout << instance_id << std::endl;

  
  //just a test...
    
  //base grid is 64*96 4km cells
#define NCELL  6144
    
  double smb[NCELL];
  for (int i =0; i < NCELL; i++)
    {
      smb[i] =  double(i)/double(NCELL);
    }
       
  double usrf[NCELL];

  double dx = 4.0e+3;
  int dims[2] = {64,96};
  int boxlo[2] = {0,0};
  int boxhi[2] = {63,95};

  int surface_flux_id = BISICLES_FIELD_SURFACE_FLUX;
  int surface_elevation_id = BISICLES_FIELD_SURFACE_ELEVATION;


  bisicles_set_2d_data(&instance_id, smb, &surface_flux_id, &dx , dims, boxlo, boxhi);
    

  bisicles_init_instance(&instance_id);


  int max_step = 0;
  double max_time = 0.0;
  int n_ice_years = 3;
  for (int n = 0; n < n_ice_years; n++)
    {
      max_step += 10;
      max_time += 1.0;

      bisicles_advance(&instance_id, &max_time, &max_step);

      for (int i =0; i < NCELL; i++)
	{
	  smb[i] +=  1.0;
	}
      
      bisicles_get_2d_data(&instance_id, usrf, 
			   &surface_elevation_id, &dx , 
			   dims, boxlo, boxhi);

      double usrf_mean = 0.0;
      for (int i = 0; i < NCELL; i++)
	{
	  usrf_mean += usrf[i];
	}
      usrf_mean /= double(NCELL);
      std::cout << "n = " << n << ", usrf_mean = " << usrf_mean << std::endl;

    }
  bisicles_free_instance(&instance_id);
}
