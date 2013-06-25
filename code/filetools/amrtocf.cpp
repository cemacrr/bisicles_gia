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
// amrtocf.cpp
// Read data from hdf5 files containing Chombo block-structured AMR hierachies
// Write unstructurd data *plus* the grid data needed to reconstruct
// the block structured data to a CF compliant netcdf file
// Or the other way round
//===========================================================================

#include <iostream>
#include "AMRIO.H"
#include "netcdf.h"

void AMRtoCF(const std::string& ifile, const std::string& ofile);
void CFtoAMR(const std::string& ifile, const std::string& ofile);

class ValidData
{

  int m_nComp;
  Vector<RealVect> m_dx;
  int m_nLevel;

  Vector<int> m_level;
  Vector<IntVect> m_iv;
  Vector<Vector< Real > > m_field;
  Vector<Vector<Real > > m_x;

  ValidData()
  {
  }

public:
  ValidData(int a_nComp, Vector<RealVect>& a_dx)
  {
    m_nComp = a_nComp;
    m_dx = a_dx;
    m_nLevel = m_dx.size();
    m_field.resize(m_nComp);
    m_x.resize(SpaceDim);
  }

  void append(const int a_lev, const IntVect& a_iv, const Vector<Real>& a_data)
  {
    CH_assert(a_lev < m_nLevel);
    CH_assert(a_data.size() == m_field.size());

    m_level.push_back(a_lev);
    m_iv.push_back(a_iv);
    for (int ic = 0; ic < m_nComp; ic++)
      {
	m_field[ic].push_back(a_data[ic]);
      }
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	m_x[dir].push_back(m_dx[a_lev][dir]*(Real(a_iv[dir]) + 0.5));  
      }
  }
  const int nComp() const
  {
    return m_nComp;
  }


  const size_t size() const
  {
    return m_level.size();
  }

  const double* field(int a_comp) const
  {
    return &m_field[a_comp][0];
  }
  const double* x(int a_dir) const
  {
    return &m_x[a_dir][0];
  }
  

};


//typedef std::pair<int,IntVect> ValidDataKey;
//typedef Vector<Real>  ValidDataValue;
//typedef std::pair<ValidDataKey, ValidDataValue > ValidDataElement;
//typedef Vector<ValidDataElement> ValidData  ;

//given block structured data produce a list of valid data
//on all levels
void BStoValid ( ValidData & , 
		 const Vector<LevelData<FArrayBox>*>&, 
		 const Vector<int>& a_ratio );

void writeCF ( const std::string&,  const ValidData& , const Vector<std::string>& );


int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif
  

    if(argc < 3) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> " 
		  << std::endl; 
	exit(0); 
      }
    
    char* in_file = argv[1];
    char* out_file = argv[2];
    
  
    std::string ifile(in_file);
    std::string ofile(out_file);
 
    if ( (ifile.size() >= 5) && (ifile.find(".hdf5") == ifile.size()-5))
      {
	AMRtoCF(ifile,ofile);
      }
    else if ((ifile.size() >= 3) && (ifile.find(".nc") == ifile.size()-3))
      {
        CFtoAMR(ifile,ofile);
      }
    else
      {
	pout() << "unknown input file type [" <<  ifile << "]" << std::endl;
	MayDay::Error("unknown input file type");
      }


   

   
    
  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
void AMRtoCF(const std::string& ifile, const std::string& ofile)
{
  

  //read the AMR data
  Vector<std::string> names;
  Vector<LevelData<FArrayBox>* > data;
  Vector<DisjointBoxLayout> grids;
  Vector<int> ratio;
  Real crseDx = 0.0, dt = 0.0, time = 0.0;
  Box domBox;
  int numLevels;
  int status = ReadAMRHierarchyHDF5
    (ifile,grids,data,names,domBox,crseDx,dt,time,
     ratio,numLevels);
  if (status != 0)
    {
      MayDay::Error("failed to read AMR hierarchy");
    }

  Vector<RealVect> dx(numLevels);
  dx[0] = crseDx * IntVect::Unit;
  for (int lev = 1; lev < numLevels; lev++)
    {
      dx[lev] = dx[lev-1]/Real(ratio[lev-1]);
    }

  RealVect x0 = RealVect::Zero;

  //extract valid data
  ValidData validData(names.size(), dx);
  BStoValid(validData, data, ratio);

  //write to CF file
  writeCF ( ofile,  validData , names );



  for (int lev = 0; lev < numLevels; lev++)
    {
      if (data[lev] != NULL)
	{
	  delete data[lev];data[lev]=NULL;
	}
    }
}
void CFtoAMR(const std::string& ifile, const std::string& ofile)
{
 MayDay::Error("CFtoAMR not implemented");
}

//given block structured data and a valid region mask, produce unstructured data 
void BStoValid ( ValidData& a_validData , 
		 const Vector<LevelData<FArrayBox>*>& a_bsData, 
		 const Vector<int>& a_ratio)
{
  
  //build valid data mask
  int numLevels = a_bsData.size();
  Vector<LevelData<BaseFab<int> >* > mask(numLevels,NULL);
  for (int lev = numLevels - 1; lev >= 0; lev--)
	  {
	    const DisjointBoxLayout& grids = a_bsData[lev]->disjointBoxLayout();
	    mask[lev] = new LevelData<BaseFab<int> >(grids,1,IntVect::Zero);
	    for (DataIterator dit(grids); dit.ok(); ++dit)
	      {
		(*mask[lev])[dit].setVal(1);
		
		if (lev < numLevels - 1)
		  {
		    const DisjointBoxLayout& fgrids = a_bsData[lev+1]->disjointBoxLayout();
		    for (DataIterator fit(fgrids) ; fit.ok(); ++fit)
		      {
			Box covered = fgrids[fit];
			covered.coarsen(a_ratio[lev]);
			covered &= grids[dit];
			if (!covered.isEmpty())
			  {
			    (*mask[lev])[dit].setVal(0,covered,0,1);
			  }
		      }
		  }
	      }
	    
	  }


  for (int lev = 0; lev < a_bsData.size(); lev++)
    {
      LevelData<FArrayBox>& levelData = *a_bsData[lev];
      LevelData<BaseFab<int> >& levelMask = *mask[lev];
      for (DataIterator dit = levelData.dataIterator(); dit.ok(); ++dit)
	{
	  const Box& box = levelData.disjointBoxLayout()[dit];
	  for (BoxIterator bit(box);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	     
	      if (levelMask[dit](iv) == 1)
		{
		  Vector<Real> v(levelData.nComp());
		  for (int ic = 0; ic < v.size(); ic++)
		    v[ic] = levelData[dit](iv,ic);
		  a_validData.append(lev,iv,v);
		  
		}
	    }
	}

      

    }
  
  //clean up mask
  for (int lev = 0; lev < mask.size(); lev++)
    {
      if (mask[lev] != NULL)
	{
	  delete mask[lev];mask[lev]=NULL;
	}
    }


}

/// write valid data to a NetCDF-CF compliant file, 
/// along with the mesh data needed to reconstruct 
/// a Chombo AMR hierarchy  from it
void writeCF ( const std::string& a_file,  const ValidData& a_validData , 
	       const Vector<std::string>& a_names )
{
  int rc; int ncID; int varID;
  //create new file
  if ( (rc = nc_create(a_file.c_str(), NC_CLOBBER, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  int nCell = a_validData.size();
  int cellDimID;
  
  //define the netCDF dimensions etc
  if ( (rc = nc_def_dim(ncID, "cell" , nCell, &cellDimID)) != NC_NOERR)
    {
      MayDay::Error("failed to define cell dimension");
    }

  std::string xname[SpaceDim] = {D_DECL("x","y","z")};
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      if ( (rc = nc_def_var(ncID, xname[dir].c_str(), NC_DOUBLE,
			    1, &cellDimID, &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to define space variable");
	}
      
      std::string s = "projection_"  + xname[dir] + "_coordinate";
      
      if ( (rc =  nc_put_att_text (ncID, varID, "standard_name", 
				   s.size(), s.c_str())) != NC_NOERR)
	{
	  MayDay::Error("failed to add standard_name attribute");
	}
      
      if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				   1, "m")) != NC_NOERR)
	{
	  MayDay::Error("failed to add units attribute to x");
	}
      
    }
  
  for (int ic = 0; ic < a_validData.nComp(); ic++)
    {
      
      const std::string name = a_names[ic];
      size_t find = name.find("/");
      if (find == string::npos )
	{
	  
	  if ( (rc = nc_def_var(ncID, a_names[ic].c_str(), NC_DOUBLE,
			    1, &cellDimID, &varID)) != NC_NOERR)
	    {
	      MayDay::Error("failed to define field variable");
	    }

	  std::string s = xname[0]; 
	  for (int dir = 1; dir < SpaceDim; dir++)
	    {
	      s+= " " + xname[dir];
	    }

	  if ( (rc =  nc_put_att_text (ncID, varID, "coordinates", 
				       s.size(), s.c_str())) != NC_NOERR)
	    {
	      MayDay::Error("failed to add field attribute");
	    }

	}

    } 



  if ( (rc = nc_enddef(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to define netcdf file");
    }

 
  //write data
  //spatial fields
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      int varID;
      if ( (rc = nc_inq_varid(ncID, xname[dir].c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
      size_t start = 0;
      size_t count = nCell;
      if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, a_validData.x(dir))) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}


    }

  //scalar fields
  for (int ic = 0; ic < a_validData.nComp(); ic++)
    {
      
      const std::string name = a_names[ic];
      size_t find = name.find("/");
      if (find == string::npos )
	{
	  if ( (rc = nc_inq_varid(ncID, name.c_str(), &varID)) != NC_NOERR)
	    {
	      MayDay::Error("failed to find variable id");
	}
	  size_t start = 0;
	  size_t count = nCell;
	  if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, a_validData.field(ic))) != NC_NOERR)
	    {
	      MayDay::Error("failed to write data");
	    } 
	}
    }
  
  //close file
  if ( (rc = nc_close(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to close netcdf file");
    }

}

