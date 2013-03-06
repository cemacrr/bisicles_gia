//===========================================================================
// bike_driver.H
//===========================================================================

#ifndef BIKE_DRIVER
#define BIKE_DRIVER

#include "DyCoreModelRegistry.H"

class DycoreToGlimmer;
class BisiclesToGlimmer;

void bike_driver_init(int argc, int exec_mode,  BisiclesToGlimmer * btg_ptr, const char * input_fname);

void bike_driver_run(BisiclesToGlimmer * btg_ptr, float& cur_time_yr, float time_inc);

void bike_driver_finalize(int amr_obj_index);

class AmrIce;
class ParmParse;

// A struct to accumulate all the objects that must perist
// for an entire cycle.
struct Bike
{
  AmrIce *amrIce;
  ParmParse *parmParse;
  
  Bike()
    :amrIce(NULL),parmParse(NULL)
  {
    ;
  }

  ~Bike();
};

int
bike_store(int model_index, Bike ** amr_object, int mode);

#endif
