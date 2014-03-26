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
#ifdef CH_MPI
  void bisicles_new_instance(int *instance_id, const char *input_fname, MPI_Comm comm);
#else
  void bisicles_new_instance(int *instance_id, const char *input_fname, int comm);
#endif

  ///free the BISICLES instance specified by instance ID 
  void bisicles_free_instance(int *instance_id);
  /// instruct BISICLES to read 2d data from an array at data_ptr.
  /**
     
   */
  void bisicles_set_2d_data(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi);
  /// instruct BISICLES to read 2d geometry (ice thickness and bedrock eleveation) data from an arrays at thck_data_ptr, topg_data_ptr.
  /**
     This special case exists because thck_data_ptr and *topg_data_ptr must be defined on the same BoxLayout
  */
  void bisicles_set_2d_geometry(int *instance_id,  double *thck_data_ptr, double *topg_data_ptr,
				const double *dx, const int *dims, 
				const int *boxlo, const int *boxhi);
  
  void bisicles_get_2d_data(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			    const int *boxlo, const int *boxhi);
  /// initialize bisicles 
  void bisicles_init_instance(int *instance_id);
  ///step BISICLES forward in time by up to max_time or max_steps (whichever is sooner)
  void bisicles_advance(int *instance_id, double *max_time, int *max_step);
  
  
  ///write a checkpoint file 
  void bisicles_write_checkpoint(int *instance_id, const char *checkpoint_fname);
  ///read a checkpoint file 
  void bisicles_read_checkpoint(int *instance_id, const char *checkpoint_fname);
  

  // FORTRAN compatibility stuff : function foo gets a f_foo_ and f_foo
  void f_bisicles_new_instance_(int *instance_id,  char *input_fname, const int *len_fname, const int *mpi_comm);

  void f_bisicles_free_instance_(int *instance_id);

  void f_bisicles_set_2d_data_(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);

  void f_bisicles_set_2d_geometry_(int *instance_id,  double *thck_data_ptr, double *topg_data_ptr,
				 const double *dx, const int *dims, 
				 const int *boxlo, const int *boxhi);

  void f_bisicles_get_2d_data_(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);

  void f_bisicles_init_instance_(int *instance_id);

  void f_bisicles_advance_(int *instance_id, double *max_time, int *max_step);

  void f_bisicles_write_checkpoint_(int *instance_id, char *checkpoint_fname, const int *len_fname);

  void f_bisicles_read_checkpoint_(int *instance_id, char *checkpoint_fname, const int *len_fname);

  void f_bisicles_new_instance(int *instance_id,  char *input_fname, const int *len_fname, const int *mpi_comm);

  void f_bisicles_free_instance(int *instance_id);

  void f_bisicles_set_2d_data(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);

  void f_bisicles_set_2d_geometry(int *instance_id,  double *thck_data_ptr, double *topg_data_ptr,
				 const double *dx, const int *dims, 
				 const int *boxlo, const int *boxhi);

  void f_bisicles_get_2d_data(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);

  void f_bisicles_init_instance(int *instance_id);

  void f_bisicles_advance(int *instance_id, double *max_time, int *max_step);

  void f_bisicles_write_checkpoint_(int *instance_id, char *checkpoint_fname, const int *len_fname);

  void f_bisicles_read_checkpoint_(int *instance_id, char *checkpoint_fname, const int *len_fname);

}
#include "NamespaceFooter.H"
#endif
