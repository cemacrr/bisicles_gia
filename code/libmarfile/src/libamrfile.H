
 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

/*===========================================================================
 libamrfile.cpp
 C-interface functions to read Chombo's amr files and access their data
==========================================================================*/

#include "NamespaceHeader.H"


extern "C"
  {
    void amr_read_file(int *status, int *amr_id, const char *file);
    void amr_read_file_R(int *status, int *amr_id, char **file);
    
    void amr_free(int *status, int *amr_id);
    void amr_query_n_level(int *status, int *n_level, const int *amr_id);
    void amr_query_n_fab(int *status, int *n_fab, const int *amr_id, const int *level_id);
#if CH_SPACEDIM == 2
    void amr_query_fab_dimensions_2d(int *status,
				     int *nx, 
				     int *ny, 
				     int *ncomp, 
				     const int *amr_id, 
				     const int *level_id, 
				     const int* fab_id);
    void amr_read_fab_data_2d(int *status, 
			      double *fab_data, 
			      double *x_data, 
			      double *y_data, 
			      const int *amr_id, 
			      const int *level_id, 
			      const int* fab_id,
			      const int* comp_id,
			      const int* nghost);
#endif 
  }

#include "NamespaceFooter.H"
