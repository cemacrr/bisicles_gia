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

#include "NamespaceHeader.H"


ValidData::ValidData(int a_nComp,  const Vector<RealVect>& a_dx,  const RealVect& a_x0 )
  {
    m_nComp = a_nComp;
    m_dx = a_dx;
    m_x0 = a_x0;
    m_nLevel = m_dx.size();
    m_field.resize(m_nComp);
    m_x.resize(SpaceDim);
  }

void ValidData::append
(const int a_lev, const IntVect& a_iv,  const Vector<Real>& a_data)
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
  if (a_nodeCoord.size() != size())
    {
      a_nodeCoord.resize(size());
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
void ValidIO::writeCF ( const std::string& a_file,  
			const ValidData& a_validData , 
			const Vector<std::string>& a_names, 
			const Transformation& a_latlonTransformation)
{
  int rc; int ncID; int varID;
  //create new file
  if ( (rc = nc_create(a_file.c_str(), NC_CLOBBER, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  int nCell = a_validData.size();
  int cellDimID;
  int nVertex = 4; //quads
  int vertexDimID;
  //define the netCDF dimensions etc
  if ( (rc = nc_def_dim(ncID, "cell" , nCell, &cellDimID)) != NC_NOERR)
    {
      MayDay::Error("failed to define cell dimension");
    }

    //define the netCDF dimensions etc
  if ( (rc = nc_def_dim(ncID, "vertex" , nVertex, &vertexDimID)) != NC_NOERR)
    {
      MayDay::Error("failed to define vertex dimension");
    }


  {
    std::string s("CF-1.6");
    if ( (rc =  nc_put_att_text (ncID, NC_GLOBAL , "Conventions", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add Convetions attribute");
      }
  }

  //definition of cell center x,y,[z]

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
      
      // if ( (rc =  nc_put_att_text (ncID, varID, "axis", 
      // 				   1, xname[dir].c_str())) != NC_NOERR)
      // 	{
      // 	  MayDay::Error("failed to add units attribute to x");
      // 	}

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

  //definition of field variables
  
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

	  std::string s = "lat lon";

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
      if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &(a_validData.x(dir)[0]))) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}


    }
  


  {
    // cell-center latitude and longitude data
    const Vector<Real>& x = a_validData.x(0);
    const Vector<Real>& y = a_validData.x(1);
      Vector<Real> lat(a_validData.size());
  Vector<Real> lon(a_validData.size());

    for (int i = 0; i < x.size(); i++)
      {
	RealVect X(x[i],y[i]);
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
    Vector<Real> nodeLat(a_validData.size());
    Vector<Real> nodeLon(a_validData.size());

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

	for (int i = 0; i < a_validData.size(); i++)
	  {
	    RealVect X(nodeX[0][i],nodeX[1][i]);
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
	  if ( (rc = nc_put_vara_double(ncID, varID, &start, &count, &(a_validData.field(ic)[0]))) != NC_NOERR)
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



#include "NamespaceFooter.H"
