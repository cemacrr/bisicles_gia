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
 C-interface functions to read Chombo's amr files, read and modify their data, 
 and write them 
==========================================================================*/
#include "libamrfile.H"
#include "AMRIO.H"
#include "NamespaceHeader.H"

class AMRHierarchy
{
  Vector<LevelData<FArrayBox>* > m_data;
  Vector<DisjointBoxLayout> m_grids;
  Vector<std::string> m_names;
  Vector<int> m_ratio;
  Vector<Real> m_dx;
  Real m_dt;
  Real m_time;
  int m_nLevel;
  bool m_ok;
  
  AMRHierarchy();

public:

  AMRHierarchy(const std::string& file)
  {
  
    Real crseDx;
    Box crseBox;

    int status = ReadAMRHierarchyHDF5(file, m_grids, m_data, m_names, crseBox, crseDx, m_dt,
				      m_time, m_ratio, m_nLevel);

    m_dx.resize(m_nLevel,crseDx);
    for (int lev = 1; lev < m_nLevel; lev++)
      {
	m_dx[lev] = m_dx[lev-1]/Real(m_ratio[lev-1]);

      } 
    for (int lev = 0; lev < m_nLevel; lev++)
      {
	m_data[lev]->exchange();
      }
    m_ok = status == 0;
  }

  bool ok() const {return m_ok;}
  const Vector<LevelData<FArrayBox>* >& data() const {return m_data;}
  const Vector<DisjointBoxLayout>& grids() const {return m_grids;}
  const Vector<Real>& dx() const  {return m_dx;}
  int nLevel() const  {return m_nLevel;}

  ~AMRHierarchy()
  {
    for (int lev = 0; lev < m_data.size(); lev++)
      {
	if (m_data[lev] != NULL)
	  {
	    delete m_data[lev];
	    m_data[lev] = NULL;
	  }
      }
  }

};

namespace libamrfile
{
  std::map<int, AMRHierarchy*> g_store; 
}

void amr_read_file(int *status, int *amr_id, const char *file)
{

  if (!status)
    return;

  std::cout << file << std::endl;

  AMRHierarchy* h = new AMRHierarchy(file);
 
  if (h->ok())
    {
      if (libamrfile::g_store.size() == 0)
	*amr_id = 0;
      else
	*amr_id  = (--libamrfile::g_store.end())->first;
      libamrfile::g_store[*amr_id] = h;
      *status = 0;
    }
  else
    {
      *status = 1;
    }
}
void amr_read_file_R(int *status, int *amr_id, char **file)
{
  amr_read_file(status,amr_id,*file);
}


void amr_free_all()
{
  for (std::map<int, AMRHierarchy*>::iterator i = libamrfile::g_store.begin();
       i != libamrfile::g_store.end(); ++i)
    {
      if (i->second != NULL)
	{
	  delete i->second;
	  i->second = NULL;
	  libamrfile::g_store.erase(i->first);
	}
    }
  
}

void amr_free(int *status, int *amr_id)
{

  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  if (i->second != NULL)
	    {
	      delete i->second;
	      i->second = NULL;
	      libamrfile::g_store.erase(i->first);
	    }
	  *status = 0;
	}
      else
	{
	  std::cout << "amr_free failed" << std::endl;
	  *status = 1;
	}
    }
  else
    {
      *status = 1;
    }

}


void amr_query_n_level(int *status, int *n_level, const int *amr_id)
{
  
  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  *n_level = i->second->nLevel();
	  *status = 0;
	}
      else
	{
	  *status = 1;
	  std::cout << "amr_query_n_level failed: bad *amr_id" << *amr_id << std::endl;
	}
    }
  else
    {
      *status = 1;
       std::cout << "amr_query_n_level: bad amr_id" << std::endl;
    }

}

void amr_query_n_fab(int *status, int *n_fab, const int *amr_id, const int *level_id)
{

  if (amr_id && level_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();
	  if (*level_id < nLevel)
	    {
	      *n_fab = i->second->grids()[*level_id].size();
	      *status = 0;
	    }
	  else
	    {
	      *status =  1;
	    }
	}
      else
	{
	  *status = 1;
	}
    }
  else
    {
      *status = 1;
    }



}

#if CH_SPACEDIM == 2
void amr_query_fab_dimensions_2d(int *status, int *nx, int *ny, int *ncomp, 
				 const int *amr_id, const int *level_id, 
				 const int *fab_id)
{

  if (!status)
    return;

  if (nx && ny && ncomp && amr_id && level_id && fab_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();

	  if (*level_id < nLevel)
	    {
	      const LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
	      DataIterator dit(i->second->grids()[*level_id]);
	      for (int j=0; j < *fab_id; j++)
		{
		  ++dit;
		}
		if (dit.ok())
		  {
		    const Box& b = i->second->grids()[*level_id][dit];
		    *nx = b.bigEnd()[0] - b.smallEnd()[0] + 1;
		    *ny = b.bigEnd()[1] - b.smallEnd()[1] + 1;
		    *ncomp = ldf.nComp();
		    *status = 0;
		  }
		else
		  {
		    *status = 1;
		  }
	    }
	  else
	    {
	      *status =  1;
	    }
	}
      else
	{
	  *status = 1;
	}
    }
  else
    {
      *status = 1;
    }
}

void amr_read_fab_data_2d(int *status, double *fab_data, double *x_data, double *y_data, 
			  const int *amr_id, const int *level_id, 
			  const int* fab_id, const int *comp_id, const int* nghost)
{

  
  if (!status)
    return;

  if (fab_data && x_data && y_data && amr_id && level_id && fab_id && comp_id && nghost)
    {

     

      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();

	  if (*level_id < nLevel)
	    {
	      const LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
	      int maxnghost=std::min(ldf.ghostVect()[0],ldf.ghostVect()[1]);

		if ( *nghost < 0 || ( *nghost >  maxnghost ))
		{
		  *status = 1;
		  return;
		}
		  
	      if (*comp_id < 0 || *comp_id >= ldf.nComp())
		{
		  *status = 1;
		  return;
		}


	      Real dx = i->second->dx()[*level_id];
	      DataIterator dit(i->second->grids()[*level_id]);
	      for (int j=0; j < *fab_id; j++)
		{
		  ++dit;
		}
		if (dit.ok())
		  {
		    Box b = i->second->grids()[*level_id][dit];
		    b.grow(*nghost);


		    const FArrayBox& fab = ldf[dit];
		    fab.linearOut((void*)fab_data,b,Interval(*comp_id,*comp_id));

		    double *xptr = x_data;
		    for (int ix = b.smallEnd()[0] ; ix <= b.bigEnd()[0]; ix++)
		      {
			*xptr++ = (double(ix)+0.5)*dx;
		      }

		    double *yptr = y_data;
		    for (int ix = b.smallEnd()[1] ; ix <= b.bigEnd()[1]; ix++)
		      {
			*yptr++ = (double(ix)+0.5)*dx;
		      }

		    *status = 0;
		  }
		else
		  {
		    *status = 1;
		  }
	    }
	  else
	    {
	      *status =  1;
	    }
	}
      else
	{
	  *status = 1;
	}
    }
  else
    {
      *status = 1;
    }





}
#endif 

#include "NamespaceFooter.H"
