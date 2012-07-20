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
 testlibamrfile.c
 test the libamrfile functions
==========================================================================*/
#include <iostream>
#include "libamrfile.H"
#include "NamespaceHeader.H"

int main(int argc, char* argv[]) {

#ifdef CH_MPI
#error MPI not supported
#endif


  int amr_id;
  int status;
  char file[32] = "testfile.2d.hdf5";

  amr_read_file(&status, &amr_id, file);

  if (status != 0)
    {
      exit(status);
    }
  
   int n_level;
  amr_query_n_level(&status, &n_level, &amr_id);
  if (status != 0)
    {
      exit(status);
    }

  std::cout << " n_level = " <<  n_level << std::endl;

  for (int lev = 0; lev < n_level; lev++)
    {

      int n_fab;
      amr_query_n_fab(&status, &n_fab, &amr_id, &lev);
      if (status != 0)
	{
	  exit(status);
	}
  
      std::cout << " level = " << lev << " n_fab = " <<  n_fab << std::endl;

      for (int i = 0; i < n_fab; i++)
	{
	  int nx, ny, ncomp;
	  amr_query_fab_dimensions_2d(&status, &nx, &ny, &ncomp, &amr_id, &lev, &i);

	  if (status != 0)
	    {
	      exit(status);
	    } 

	  

	  double *fab_data = new double[nx*ny];
	  double *x_data = new double[nx];
	  double *y_data = new double[ny];
	  int comp = 0;
	  int nghost = 1;
	  amr_read_fab_data_2d(&status, fab_data , x_data, y_data, &amr_id, &lev, &i, &comp, &nghost);


	  std::cout << " level = " << lev << " fab = " <<  i 
		    << " nx = " << nx << " ny = " << ny << " ncomp = " << ncomp
	            << " " << x_data[0] << " <= x <= " << x_data[nx-1] 
		    << " "  << y_data[0] << " <= y <= " << y_data[ny-1] 
		    << std::endl;

	  if (status != 0)
	    {
	      exit(status);
	    } 

	  delete fab_data; 
	  delete x_data;
	  delete y_data;

	}
    }


  amr_free(&status, &amr_id);
  if (status != 0)
    {
      exit(status);
    } 


  return 0;
}

#include "NamespaceFooter.H"
