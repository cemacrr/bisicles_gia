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
  Box m_crseBox;

  AMRHierarchy();

public:

  AMRHierarchy(const std::string& file)
  {
  
    Real crseDx;

    int status = ReadAMRHierarchyHDF5(file, m_grids, m_data, m_names, 
				      m_crseBox, crseDx, m_dt,
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
  const Vector<std::string>& names() const  {return m_names;}
  int nLevel() const  {return m_nLevel;}

  int write(const std::string& file)
  {
    int status = LIBAMRFILE_ERR_WRITE_FAILED;
    if (m_ok)
      {
	WriteAMRHierarchyHDF5(file, m_grids, m_data, 
			      m_names, m_crseBox, m_dx[0], m_dt,
			      m_time, m_ratio, m_nLevel);
      }
    if (status != 0)
      status = LIBAMRFILE_ERR_WRITE_FAILED;

    return status;
  }

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
      *status = LIBAMRFILE_ERR_READ_FAILED;
    }
}
void amr_read_file_R(int *status, int *amr_id, char **file)
{
  amr_read_file(status,amr_id,*file);
}

void amr_write_file(int *status, int *amr_id, const char *file)
{

  if (!status)
    return;

  if (!amr_id)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }

  std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
  if (i == libamrfile::g_store.end())
    {
      *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
      return;
    }

  AMRHierarchy* h = i->second;
  if (!h || !(h->ok()))
    {
      *status = LIBAMRFILE_ERR_BAD_AMR_HIERARCHY;
      return;
    }

  *status = h->write(file);

}


void amr_write_file_R(int *status, int *amr_id, char **file)
{
  amr_write_file(status,amr_id,*file);
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
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
    }

}

void amr_query_comp_name(int *status, char *file, const int* amr_id, const int* comp, const int* buflen)
{

  if (!status)
    return;

  if (amr_id && comp && buflen && file)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  if (i->second)
	    {
	      const AMRHierarchy& h = *i->second;
	      if (*comp < h.names().size())
		{
		  int c = *comp;
		  int l = *buflen;
		  strncpy( file, h.names()[c].c_str(), l); 
		}
	    }
	  else
	    {
	      *status = LIBAMRFILE_ERR_NULL_POINTER;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }

}

void amr_query_comp_name_R(int *status, char **file, const int* amr_id, const int* comp, const int* buflen)
{
  amr_query_comp_name(status, *file, amr_id, comp, buflen);
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
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}

void amr_query_n_comp(int *status, int* n_comp, const int* amr_id)
{

  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  AMRHierarchy* h = i->second;
	  if (h)
	    {
	      *n_comp = h->data()[0]->nComp();
	      *status = 0;
	    }
	  else
	    {
	      *status = LIBAMRFILE_ERR_NULL_POINTER;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
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
	      *status =  LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }



}

void amr_query_fab_dimensions(int *status, 
			      D_DECL(int *nx, int *ny, int *nz),
			      int *ncomp, 
			      const int *amr_id, 
			      const int *level_id, 
			      const int* fab_id)
{
  
  if (!status)
    return;

  if (!nx)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }

#if CH_SPACEDIM > 1
  if (!ny)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }
#if CH_SPACEDIM > 2
  if (!nz)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }
#endif
#endif


  if (ncomp && amr_id && level_id && fab_id)
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
#if CH_SPACEDIM > 1
		    *ny = b.bigEnd()[1] - b.smallEnd()[1] + 1;
#if CH_SPACEDIM > 2
		    *nz = b.bigEnd()[2] - b.smallEnd()[2] + 1;
#endif
#endif
		    *ncomp = ldf.nComp();
		    *status = 0;
		  }
		else
		  {
		    *status = LIBAMRFILE_ERR_NO_SUCH_FAB;
		  }
	    }
	  else
	    {
	      *status =  LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
    }
}

void amr_read_fab_data(int *status, 
		       double *fab_data, 
		       D_DECL(double *x_data, double *y_data, double *z_data),
		       const int *amr_id, 
		       const int *level_id, 
		       const int* fab_id,
		       const int* comp_id,
		       const int* nghost)
{

  
  if (!status)
    return;

  if (!x_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER  ;
      return;
    }

#if CH_SPACEDIM > 1
  if (!y_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }
#if CH_SPACEDIM > 2
  if (!z_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }
#endif
#endif
  
  if (fab_data  && amr_id && level_id && fab_id && comp_id && nghost)
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
		  *status = LIBAMRFILE_ERR_BAD_NGHOST;
		  return;
		}
		  
	      if (*comp_id < 0 || *comp_id >= ldf.nComp())
		{
		  *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
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
#if CH_SPACEDIM > 1
		    double *yptr = y_data;
		    for (int ix = b.smallEnd()[1] ; ix <= b.bigEnd()[1]; ix++)
		      {
			*yptr++ = (double(ix)+0.5)*dx;
		      }
#if CH_SPACEDIM > 2
		    double *zptr = y_data;
		    for (int ix = b.smallEnd()[2] ; ix <= b.bigEnd()[2]; ix++)
		      {
			*zptr++ = (double(ix)+0.5)*dx;
		      }
#endif
#endif 


		    *status = 0;
		  }
		else
		  {
		    *status = LIBAMRFILE_ERR_NO_SUCH_FAB;
		  }
	    }
	  else
	    {
	      *status =   LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}

void amr_write_fab_data(int *status, 
			double *fab_data, 
			D_DECL(int *nx, int *ny, int *nz),
			const int *amr_id, 
			const int *level_id, 
			const int* fab_id,
			const int* comp_id,
			const int* nghost)
{
  if (!status)
    return;
  
   if (!nx)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }

#if CH_SPACEDIM > 1
   if (!ny)
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER ;
       return;
     }
#if CH_SPACEDIM > 2
   if (!nz)
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER ;
       return;
     }
#endif
#endif

   if ( !(fab_data && amr_id && level_id && fab_id && comp_id && nghost))
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER ;
       return;
     }
     
   std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
   if (i == libamrfile::g_store.end())
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
       return;
     }

   AMRHierarchy* h = i->second;
   if (!h)
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER;
       return;
     }

   int nLevel =  i->second->nLevel();
   if (*level_id >= nLevel)
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_LEVEL;
       return;
     }
   
   LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
   int maxnghost=std::min(ldf.ghostVect()[0],ldf.ghostVect()[1]);
   
   if (*nghost < 0 || ( *nghost >  maxnghost ))
     {
       *status = LIBAMRFILE_ERR_BAD_NGHOST;
       return;
     }
   
   if (*comp_id < 0 || *comp_id >= ldf.nComp())
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
       return;
     }

   DataIterator dit(i->second->grids()[*level_id]);

   for (int j=0; j < *fab_id; j++)
     {
       ++dit;
     }
   if (! dit.ok())
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_FAB;
       return;
     }

   Box b = i->second->grids()[*level_id][dit];
   b.grow(*nghost);
   FArrayBox& fab = ldf[dit];
   fab.linearIn((void*)fab_data,b,Interval(*comp_id,*comp_id));

   status = 0;
		    
}


#include "NamespaceFooter.H"
