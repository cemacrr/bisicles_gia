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
// ValidIO.cpp
// Operations related to reading and writing data from/to valid regions
// of AMR Hierarchies
//===========================================================================
#include "ValidIO.H"
#include "FieldNames.H"
#include "NamespaceHeader.H"


ValidData::ValidData(int a_nComp,  const RealVect& a_crseDx, const Vector<int>& a_ratio, const RealVect& a_x0)
{
  define(a_nComp, a_crseDx, a_ratio,  a_x0);
}

void ValidData::define(int a_nComp,  const RealVect& a_crseDx, const Vector<int>& a_ratio, const RealVect& a_x0)
{
  m_nComp = a_nComp;
  m_ratio = a_ratio;
  m_x0 = a_x0;
  m_nLevel = m_ratio.size();
  m_dx.resize(m_nLevel);
  m_dx[0] = a_crseDx;
  for (int lev = 1; lev < m_nLevel; lev++)
    m_dx[lev] = m_dx[lev-1] / Real(m_ratio[lev-1]);
    
  m_field.resize(m_nComp);
  m_x.resize(SpaceDim);
  m_iv.resize(SpaceDim);
  m_isDefined = true;
}


void ValidData::append
(int a_lev, int a_boxID, const IntVect& a_iv,  const Vector<Real>& a_data)
{
  CH_TIME("ValidData::append");

  CH_assert(m_isDefined);
  CH_assert(a_lev < m_nLevel);
  CH_assert(a_data.size() == m_field.size());

    m_level.push_back(a_lev);
    m_boxID.push_back(a_boxID);

    for (int dir = 0; dir < SpaceDim; dir++)
      {
	m_iv[dir].push_back(a_iv[dir]);  
      }
   
    for (int ic = 0; ic < m_nComp; ic++)
      {
	m_field[ic].push_back(a_data[ic]);
      }

    for (int dir = 0; dir < SpaceDim; dir++)
      {
	m_x[dir].push_back(m_x0[dir] + m_dx[a_lev][dir]*(Real(a_iv[dir]) + 0.5));  
      }
  }

///fill  a_nodeCoord with the a_dir co-ordinate of node a_node. 
/**
   \param Vector<Real>& a_nodeCoord will be resized to this.size(). 
   \param a_node[dir] = 0 implies lo side,   1 implies high side
                        
*/
void ValidData::computeNodeCoord(Vector<Real>& a_nodeCoord, int a_dir, const IntVect& a_node) const
{
  
  CH_TIME("ValidData::computeNodeCoord");
  CH_assert(m_isDefined);
  if (a_nodeCoord.size() != nCell())
    {
      a_nodeCoord.resize(nCell());
    }

  CH_assert( a_node[a_dir] == 0 || a_node[a_dir]  == 1);

  Real f =  (a_node[a_dir] == 1)?0.5:-0.5;
  for (int i =0; i <  a_nodeCoord.size(); i++)
    {
      a_nodeCoord[i] = m_x[a_dir][i] + f * m_dx[m_level[i]][a_dir];
    }

}

//given block structured data and a valid region mask, produce unstructured data 
void ValidIO::BStoValid 
( ValidData& a_validData , 
  const Vector<LevelData<FArrayBox>*>& a_bsData, 
  const Vector<int>& a_ratio)
{
  CH_TIME("ValidIO::BStoValid");
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
	  int boxID = dit.i().datInd();


	  for (BoxIterator bit(box);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	     
	      if (levelMask[dit](iv) == 1)
		{
		  Vector<Real> v(levelData.nComp());
		  for (int ic = 0; ic < v.size(); ic++)
		    v[ic] = levelData[dit](iv,ic);
		  a_validData.append(lev,boxID, iv,v);
		  
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

/// read valid data from a NetCDF-CF compliant file
/// along with the mesh data needed to reconstruct 
/// a Chombo AMR hierarchy  from it

void ValidIO::readCF ( ValidData& a_validData, 
		       Vector<std::string>& a_names, 
		       const std::string& a_file)
{
  
  CH_TIME("ValidIO::readCF");
  int rc; int ncID; int varID;
  std::string ivname[SpaceDim] = {D_DECL("i","j","k")};
 
  
  //open file
  if ( (rc = nc_open(a_file.c_str(), NC_NOWRITE, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  //read number of levels, coarse dx, and refinement ratios
  {  
    int levelDimID;
    if ( (rc = nc_inq_dimid(ncID, "level", &levelDimID) ) != NC_NOERR)
      {
	MayDay::Error("failed to detemine level dimension id ");
      }
    
    size_t nLevel;
    if ( (rc = nc_inq_dimlen(ncID, levelDimID, &nLevel) ) != NC_NOERR)
      {
	MayDay::Error("failed to detemine number of levels ");
      }

    int ratioID;
    if ( (rc = nc_inq_varid(ncID, "amr_ratio", &ratioID) ) != NC_NOERR)
      {
	MayDay::Error("failed to determine refinement ratios id");
      }
    
    Vector<int> ratio(nLevel);
    size_t start(0); size_t count(nLevel);
    if ( ( rc =  nc_get_vara_int(ncID, ratioID, &start, &count, &ratio[0]) ) != NC_NOERR)
      {
	MayDay::Error("nc_get_vara_int failed to read ratio");
      } 

    a_validData.define(1, RealVect::Unit, ratio);

  }

  if ( (rc = nc_close(ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to close netcdf file");
    }
}


/// write valid data to a NetCDF-CF compliant file, 
/// along with the mesh data needed to reconstruct 
/// a Chombo AMR hierarchy  from it
void ValidIO::writeCF ( const std::string& a_file,  
			const ValidData& a_validData , 
			const Vector<std::string>& a_names, 
			const Transformation& a_latlonTransformation)
{

  CH_TIME("ValidIO::writeCF");

  int rc; int ncID; int varID;

  std::string ivname[SpaceDim] = {D_DECL("i","j","k")};
  std::string xname[SpaceDim] = {D_DECL("x","y","z")};

  //create new file
  if ( (rc = nc_create(a_file.c_str(), NC_CLOBBER, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  size_t nCell = a_validData.nCell();
  int cellDimID;
  size_t nVertex = 4; //quads
  int vertexDimID;
  size_t nLevel = a_validData.nLevel();
  int levelDimID;
  
  //define the netCDF dimensions 
  if ( (rc = nc_def_dim(ncID, "cell" , nCell, &cellDimID)) != NC_NOERR)
    {
      MayDay::Error("failed to define cell dimension");
    }

   
  if ( (rc = nc_def_dim(ncID, "vertex" , nVertex, &vertexDimID)) != NC_NOERR)
    {
      MayDay::Error("failed to define vertex dimension");
    }

  
  if ( (rc = nc_def_dim(ncID, "level" , nLevel, &levelDimID)) != NC_NOERR)
    {
      MayDay::Error("failed to define level dimension");
    }

  //set global attributes
  {
    std::string s("CF-1.6");
    if ( (rc =  nc_put_att_text (ncID, NC_GLOBAL , "Conventions", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add Conventions attribute");
      }
  }

  //definition of data related to block structured mesh : 
  //(ratio, dx, level,boxID, IntVect components)
  {
    std::string s = "amr_ratio" ;
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_INT,
			  1, &levelDimID, &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define amr_ratio variable");
      }

     if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				 1, "1")) != NC_NOERR)
      {
	MayDay::Error("failed to add units attribute");
      }
    
    if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add long_name attribute");
      }

  }
   {
    std::string s = "amr_dx" ;
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  1, &levelDimID, &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define amr_dx variable");
      }
    if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				 1, "m")) != NC_NOERR)
      {
	MayDay::Error("failed to add units attribute");
      }
    
    if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add long_name attribute");
      }

  }
  

  {
    std::string s = "amr_level" ;
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_INT,
			  1, &cellDimID, &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define level variable");
      }
    
    if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				 1, "1")) != NC_NOERR)
      {
	MayDay::Error("failed to add units attribute");
      }
    
    if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add long_name attribute");
      }

  }

  {
    std::string s = "amr_boxid" ;
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_INT,
			  1, &cellDimID, &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define boxid variable");
      }
    
    
    
    if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add long_name attribute");
      }


    if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				 1, "1")) != NC_NOERR)
      {
	MayDay::Error("failed to add units attribute");
      }
    
  }
  {
    
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	std::string s = "amr_levelgrid_" + ivname[dir] + "_index"; 

	if ( (rc = nc_def_var(ncID, s.c_str(), NC_INT,
			      1, &cellDimID, &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to define IntVect component variable");
	  }
	
	
	if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
				     s.size(), s.c_str())) != NC_NOERR)
	  {
	    MayDay::Error("failed to add long_name attribute");
	  }
	
	if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				     1, "1")) != NC_NOERR)
	  {
	    MayDay::Error("failed to add units attribute");
	  }


      }
  }
  

  //definition of cell center x,y,[z]
  
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
      
      s = xname[dir] + "_vertices";
      if ( (rc =  nc_put_att_text (ncID, varID, "bounds", 
				   s.size(), s.c_str())) != NC_NOERR)
	{
	  MayDay::Error("failed to add bounds attribute to x");
	}


    }

  //definition of cell-centered lat and lon variables
  {
    std::string s("lat");
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  1, &cellDimID, &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define lat variable");
      }
    s = "latitude";
    if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add long_name attribute to lat");
      }
    s = "degrees_north";
    if ( (rc =  nc_put_att_text (ncID, varID, "units", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add units attribute to lat");
      }
    s = "lat_vertices";
    if ( (rc =  nc_put_att_text (ncID, varID, "bounds", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add bounds attribute to lat");
      }
    
  }
  {
    std::string s("lon");
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  1, &cellDimID, &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define lat variable");
      }  

    s = "longitude";
    if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add long_name attribute to lon");
      }
    s = "degrees_east";
    if ( (rc =  nc_put_att_text (ncID, varID, "units", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add units attribute to lon");
      }
    s = "lon_vertices";
    if ( (rc =  nc_put_att_text (ncID, varID, "bounds", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add bounds attribute to lon");
      }

  }
  
  //definition of node-centered x/y values
  {
    int dimID[2] = {cellDimID, vertexDimID};

    std::string s("x_vertices");
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  2, &dimID[0], &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define x_vertices variable");
      }

    s = "y_vertices";
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  2, &dimID[0], &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define y_vertices variable");
      }
  }

  
  //definition of node-centered lat/lon values
  {
    int dimID[2] = {cellDimID, vertexDimID};

    std::string s("lat_vertices");
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  2, &dimID[0], &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define lat_vertices variable");
      }

    s = "lon_vertices";
    if ( (rc = nc_def_var(ncID, s.c_str(), NC_DOUBLE,
			  2, &dimID[0], &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define lat_vertices variable");
      }
  }

  //convert the input names to CF names
  Vector<FieldNames::CFRecord> cfRecord;
  FieldNames::CFLookup(cfRecord,  a_names);
  

  //definition of field variables
  for (int ic = 0; ic < a_validData.nComp(); ic++)
    {
       
      if ( (rc = nc_def_var(ncID, cfRecord[ic].name().c_str(), NC_DOUBLE,
			    1, &cellDimID, &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to define field variable");
	}
      
      std::string s = "lat lon x y";
      
      if ( (rc =  nc_put_att_text (ncID, varID, "coordinates", 
				   s.size(), s.c_str())) != NC_NOERR)
	{
	  MayDay::Error("failed to add coordinates attribute");
	}
      
      s = cfRecord[ic].standardName();
      if (s.size() > 0)
	{
	  if ( (rc =  nc_put_att_text (ncID, varID, "standard_name", 
				       s.size(), s.c_str())) != NC_NOERR)
	    {
	      MayDay::Error("failed to add standard_name attribute");
	    }
	}

      s = cfRecord[ic].longName();
      if (s.size() > 0)
	{
	  if ( (rc =  nc_put_att_text (ncID, varID, "long_name", 
				       s.size(), s.c_str())) != NC_NOERR)
	    {
	      MayDay::Error("failed to add long_name attribute");
	    }
	}

      s = cfRecord[ic].units();
      if (s.size() > 0)
	{
	  if ( (rc =  nc_put_att_text (ncID, varID, "units", 
				       s.size(), s.c_str())) != NC_NOERR)
	    {
	      MayDay::Error("failed to add units attribute");
	    }
	}
    } 
  

  if ( (rc = nc_enddef(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to define netcdf file");
    }

 
  //write data

  //amr data / ratio
  {
    int varID;
    if ( (rc = nc_inq_varid(ncID, "amr_ratio", &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to find variable amr_ratio");
      }
    size_t start = 0;
    size_t count = nLevel;
    if ( (rc = nc_put_vara_int(ncID, varID, &start, &count, &(a_validData.ratio()[0]))) != NC_NOERR)
      {
	MayDay::Error("failed to write data");
      }
  }

  {
    //amr data / dx
    int varID;
    if ( (rc = nc_inq_varid(ncID, "amr_dx", &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to find variable amr_ratio");
      }
    size_t start = 0;
    size_t count = nLevel;

    Vector<Real> dx(nLevel);
    for (int lev = 0; lev < nLevel; lev++)
      {
	dx[lev] = a_validData.dx()[lev][0];
      }

    if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &dx[0])) != NC_NOERR)
      {
	MayDay::Error("failed to write data");
      }
  }

  //amr data / level
  {
    int varID;
    if ( (rc = nc_inq_varid(ncID, "amr_level", &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to find variable amr_level");
      }
    size_t start = 0;
    size_t count = nCell;
    if ( (rc = nc_put_vara_int(ncID, varID, &start, &count, &(a_validData.level()[0]))) != NC_NOERR)
      {
	MayDay::Error("failed to write data");
      }
  }

  // amr data / box id
  {
    int varID;
    if ( (rc = nc_inq_varid(ncID, "amr_boxid", &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to find variable amr_boxid");
      }
    size_t start = 0;
    size_t count = nCell;
    if ( (rc = nc_put_vara_int(ncID, varID, &start, &count, &(a_validData.boxID()[0]))) != NC_NOERR)
      {
	MayDay::Error("failed to write data");
      }
  }

  // amr data / grid vectors 
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      std::string s = "amr_levelgrid_" + ivname[dir] + "_index"; 
      int varID;
      if ( (rc = nc_inq_varid(ncID, s.c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable amr_boxid");
	}
      size_t start = 0;
      size_t count = nCell;
      if ( (rc = nc_put_vara_int(ncID, varID, &start, &count, &(a_validData.iv(dir)[0]))) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}
    }


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
      if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &(a_validData.x(dir)[0]))) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}


    }
  


  {
    // cell-center latitude and longitude data
    D_TERM(const Vector<Real>& x = a_validData.x(0);,
	   const Vector<Real>& y = a_validData.x(1);,
	   const Vector<Real>& z = a_validData.x(2););
      
    Vector<Real> lat(nCell);
    Vector<Real> lon(nCell);
    
    for (int i = 0; i < x.size(); i++)
      {
	RealVect X(D_DECL(x[i],y[i],z[i]));
	RealVect L = a_latlonTransformation.transform(X);
	lat[i] = L[0] / M_PI * 180.0;
	lon[i] = L[1] / M_PI * 180.0;
      }
    
    int varID;
    size_t start = 0;
    size_t count = nCell;

    if ( (rc = nc_inq_varid(ncID,"lat", &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to find variable id");
      }
    
    if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &(lat[0]))) != NC_NOERR)
      {
	MayDay::Error("failed to write lat data");
      }
    
      if ( (rc = nc_inq_varid(ncID,"lon", &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
      
      if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &(lon[0]))) != NC_NOERR)
	{
	  MayDay::Error("failed to write lon data");
	}	
	   
  }

  {
    // cell node x,y, latitude and longitude data

    Vector< Vector<Real> > nodeX(SpaceDim) ;
    Vector<Real> nodeLat(nCell);
    Vector<Real> nodeLon(nCell);

    Box unit(IntVect::Zero,IntVect::Unit);
    //for rectangular cells, iterating across every node
    //is a bit wasteful as e.g x(0,1) = x(0,0), but maybe
    //one day we might have mapped co-oords. Presumably the
    //I/O will dominate anyway
    size_t start[2] = {0,0};
    size_t count[2] = {nCell,1};
    for (BoxIterator bit(unit);bit.ok();++bit)
      {
	const IntVect& iv = bit();
	for (int dir = 0; dir < SpaceDim; dir++)
	  {
	    a_validData.computeNodeCoord( nodeX[dir], dir, iv);
	  }

	for (int i = 0; i < nCell; i++)
	  {
	    RealVect X(D_DECL(nodeX[0][i],nodeX[1][i],nodeX[2][i]));
	    RealVect L = a_latlonTransformation.transform(X);
	    nodeLat[i] = L[0] / M_PI * 180.0;
	    nodeLon[i] = L[1] / M_PI * 180.0;
	    


	  }

	if ( (rc = nc_inq_varid(ncID,"x_vertices", &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to find variable id");
	  }
	
	if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeX[0][0]))) != NC_NOERR)
	  {
	    MayDay::Error("failed to write x_vertices data");
	  }
	
	if ( (rc = nc_inq_varid(ncID,"y_vertices", &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to find variable id");
	  }
	
	if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeX[1][0]))) != NC_NOERR)
	  {
	    MayDay::Error("failed to write y_vertices data");
	  }
	


	if ( (rc = nc_inq_varid(ncID,"lat_vertices", &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to find variable id");
	  }
	
	if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeLat[0]))) != NC_NOERR)
	  {
	    MayDay::Error("failed to write lat_vertices data");
	  }
	
	
	    if ( (rc = nc_inq_varid(ncID,"lon_vertices", &varID)) != NC_NOERR)
	      {
		MayDay::Error("failed to find variable id");
	      }
	    
	    if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeLon[0]))) != NC_NOERR)
	      {
		MayDay::Error("failed to write lon_vertices data");
	      }
	  

	start[1] += 1;
      }
  
  }

  //scalar fields
  for (int ic = 0; ic < a_validData.nComp(); ic++)
    {
      
      
      if ( (rc = nc_inq_varid(ncID, cfRecord[ic].name().c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
      size_t start = 0;
      size_t count = nCell;
      if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &(a_validData.field(ic)[0]))) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	} 
    }
  
  
  //close file
  if ( (rc = nc_close(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to close netcdf file");
    }

}



#include "NamespaceFooter.H"
