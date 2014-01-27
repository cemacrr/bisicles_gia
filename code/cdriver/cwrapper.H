#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef CWRAPPER_H
#define CWRAPPER_H
#include "NamespaceHeader.H"
extern "C"
{
  ///create a new BISICLES instance, return an ID 
  void bisicles_new_instance(int *instance_id, const char *input_fname);

  void bisicles_new_instance_(int *instance_id,  char *input_fname, const int *len_fname);
  

  ///free the BISICLES instance specified by instance ID 
  void bisicles_free_instance(int *instance_id);

  void bisicles_free_instance_(int *instance_id);
  /// instruct BISICLES to read 2d data from an array at data_ptr.
  /**
     
   */
  

  void bisicles_set_2d_data(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi);

  void bisicles_set_2d_data_(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);
 

  void bisicles_get_2d_data(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			    const int *boxlo, const int *boxhi);

  void bisicles_get_2d_data_(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);


  /// initialize bisicles 
  void bisicles_init_instance(int *instance_id);
  void bisicles_init_instance_(int *instance_id);


  ///step BISICLES forward in time by up to max_time or max_steps (whichever is sooner)
  void bisicles_advance(int *instance_id, double *max_time, int *max_step);
  void bisicles_advance_(int *instance_id, double *max_time, int *max_step);

}
#include "NamespaceFooter.H"
#endif
