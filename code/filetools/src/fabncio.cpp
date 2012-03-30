 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
//===========================================================================
// ncio.cpp
// read/write FABs from netcdf files
//===========================================================================

#include "fabncio.H"
#include "netcdf.h"
#include "NamespaceHeader.H"

void writeNetCDF(const std::string& a_file,
		 const Vector<std::string>& a_names,
		 const FArrayBox& a_fab, 
		 const Real& a_dx)
{
  int rc; int ncID;
  if ( (rc = nc_create(a_file.c_str(), NC_CLOBBER, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  std::string xname[SpaceDim] = {"x","y"};
  int dimID[SpaceDim];

  int bufSize = 1;
  size_t start[SpaceDim] = {0,0};
  size_t count[SpaceDim];
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      int n = 1 + a_fab.box().bigEnd()[dir] - a_fab.box().smallEnd()[dir];
      count[dir] = n;
      bufSize*=n;
      if ( (rc = nc_def_dim(ncID, xname[dir].c_str(), n, &dimID[dir])) != NC_NOERR)
	{
	  MayDay::Error("failed to define dimension");
	}
    }

  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    
    int varID;
    int one = 1;
    if ( (rc = nc_def_var(ncID, xname[dir].c_str(), NC_DOUBLE,
			  one, &dimID[dir], &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define x");
      }
  }
  int FORTRAN_dimID[SpaceDim];
  size_t FORTRAN_start[SpaceDim];
  size_t FORTRAN_count[SpaceDim];
  for (int i = 0 ; i < SpaceDim ; i++)
    {
      FORTRAN_dimID[i] = dimID[SpaceDim-i-1];
      FORTRAN_count[i] = count[SpaceDim-i-1];
      FORTRAN_start[i] = start[SpaceDim-i-1];
    }
  for (int i =0; i < a_names.size(); i++)
    {
      int varID;
      if ( (rc = nc_def_var(ncID, a_names[i].c_str(), NC_DOUBLE,
			    SpaceDim, FORTRAN_dimID, &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to define variable");
	}
    }

  if ( (rc = nc_enddef(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to define netcdf file");
    }
  
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      int varID;
      if ( (rc = nc_inq_varid(ncID, xname[dir].c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
      double *xptr = new double[count[dir]];

      xptr[0] =  (Real(a_fab.box().smallEnd()[dir])+0.5)*a_dx;
      for (int i = 1; i < count[dir]; ++i)
	{
	  xptr[i] = xptr[i-1] + a_dx; 
	}
      if ( (rc = nc_put_vara_double(ncID, varID, &start[dir], &count[dir], xptr)) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}
      
      delete xptr;
    }

  double* dptr = new double[bufSize];
  for (int i =0; i < a_names.size(); i++)
    {
      int varID;
      if ( (rc = nc_inq_varid(ncID, a_names[i].c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
       
      Interval ivl(i,i);
      a_fab.linearOut((void*)dptr,a_fab.box(),ivl);
      if ( (rc = nc_put_vara_double(ncID, varID, FORTRAN_start, FORTRAN_count, dptr)) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}

    }

  if ( (rc = nc_close(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to close netcdf file");
    }
  if (dptr != NULL)
    delete[] dptr;
 
}

//construct a fab with n=a_var.size() components,
//filled with data loaded from a netcdf file
//with a one-dimensional variable x and 
//SpaceDim-dimensional variables var[0]-var[n-1]
void readNetCDF(const std::string& a_file,
		const Vector<std::string>& a_var,
		FArrayBox& a_fab,
		Real& a_dx)
{

  int rc; int ncID;
  if ( (rc = nc_open(a_file.c_str(), NC_NOWRITE, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  {
    //determine a_dx
    int xID;
    if ( (rc = nc_inq_varid(ncID, "x", &xID) ) != NC_NOERR)
      {
	MayDay::Warning("nc_inq_varid failed to find id for x : setting dx = 5.0km");
	a_dx = 5.0e+3;
      }
    else
      {
	int nxdims;
	if ( (rc = nc_inq_varndims(ncID, xID, &nxdims) )  != NC_NOERR)
	  {
	    MayDay::Error("nc_inq_varndims failed to find ndims for x");
	  }
	if (nxdims != 1)
	  {
	    MayDay::Error("wrong dimensions in x");
	  }
	double x[2];
	size_t start[2] = {0,0};
	size_t count[2] = {2,2};
	if ( ( rc =  nc_get_vara_double(ncID, xID, start, count, x) ) != NC_NOERR)
	  {
	    MayDay::Error("nc_get_vara_double failed to read first two values of x");
	  }
	a_dx = x[1] - x[0];
      }
    if (a_dx <= 0.0)
      MayDay::Error("a_dx <= 0.0");
  }

  double* dptr = NULL;
  Vector<int> varID(a_var.size());
  size_t dimLength[SpaceDim];
  for (int i =0; i < a_var.size(); ++i)
    {
      if ( (rc = nc_inq_varid(ncID, a_var[i].c_str(), &varID[i]) ) != NC_NOERR)
	{
	  MayDay::Error("nc_inq_varid failed");
	}

      int ndims;
      if ( (rc = nc_inq_varndims(ncID, varID[i], &ndims ) ) != NC_NOERR)
	{
	  MayDay::Error("nc_inq_varndims failed");
	}
      if (ndims != SpaceDim)
	{
	  MayDay::Error("wrong dimensions in variable");
	}

      if (i == 0)
	{
	  int dimID[SpaceDim];

	  if ( (rc = nc_inq_vardimid(ncID, varID[i], dimID ) ) != NC_NOERR)
	    {
	      MayDay::Error("nc_inq_vardimid  failed");
	    }
	  
	  
	  
	  IntVect hi; int bufSize = 1;
	  for (int dir = 0; dir < SpaceDim; ++dir)
	    {
	      
	      if ( (rc =  nc_inq_dimlen  (ncID, dimID[dir], &dimLength[dir]) ) != NC_NOERR)
		{
		  MayDay::Error("nc_inq_dimlen failed");
		}
	     
	      bufSize *= dimLength[dir];
	    }
	  if (SpaceDim == 2)
	    {
	      hi[0] = dimLength[1]-1;
	      hi[1] = dimLength[0]-1;
	    }
	  else
	    {
	      MayDay::Error("2d only for now");
	    }
	  a_fab.define(Box(IntVect::Zero,hi), a_var.size());
	  
	  dptr = new double[bufSize];

	 

	} //end if (i == 0)
      size_t start[SpaceDim] = {0,0};
      if ( ( rc =  nc_get_vara_double(ncID, varID[i], start, dimLength, dptr) ) != NC_NOERR)
	{
	  MayDay::Error("nc_get_vara_double failed");
	}
      a_fab.linearIn(dptr,a_fab.box(),Interval(i,i));

    }

  if (dptr != NULL)
    delete[] dptr;

 if ( (rc = nc_close(ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to close netcdf file");
    }
}



#include "NamespaceFooter.H"
