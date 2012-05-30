#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// concrete class encapsulating surface fluxes determined  
// by copying, coarsening, or interpolating a LevelData<FArrayBox>
// covering an entire domain. Objects can be defined either
// by specifying a RefCountedPtr<LevelData<FArrayBox> > , or by specifying
// a std::map<Real,string> mapping time ti to an hdf5 file f. In the
// latter case, the flux at time t is found by linear interploation
// between time max(ti <= t) and min(ti > t) 

// \todo replace the std::map<Real,string> mechanism with
// a suitable abstraction

#include "LevelDataSurfaceFlux.H"
#include "ReadLevelData.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "NamespaceHeader.H"

void LevelDataSurfaceFlux::surfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{
    
  if (m_timeFileMap != NULL)
    {
      std::map<Real,std::string>::const_iterator start,end;
      if (m_timeFileMap->size() > 1)
	{
	  end = ++m_timeFileMap->begin();
	  while (end != m_timeFileMap->end() && end->first <= a_time)
	    {
	      end++;
	    }
	  start = end;
	  --start;
	  if (end == m_timeFileMap->end())
	    end = start;
	}
      else
	{
	  start = end = m_timeFileMap->begin();
	}
      
      Vector<std::string> name(1,m_name);
      
      if (start->first != m_startTime)
	{
	  //load m_startFlux
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_startFlux);
	  Real dx;
	  readLevelData(data,dx,start->second,name,1);
	  for (int dir=0;dir<SpaceDim;dir++)
	    m_dx[dir] = dx;
	  m_startTime = start->first;
	}
      if (end->first != m_endTime)
	{
	  if (start == end)
	    {
	      m_endFlux = m_startFlux;
	    }
	  else
	    {
	      //load m_endFlux
	      Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_endFlux);
	      Real dx;
	      readLevelData(data, dx , start->second,name,1);
	      for (int dir=0;dir<SpaceDim;dir++)
		CH_assert(m_dx[dir] = dx);
	      m_endTime = end->first;
	    }
	  m_endTime = end->first;
	}
    }
  pout() << " LevelDataSurfaceFlux::m_startTime = " <<   m_startTime;
  
  for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }

  FillFromReference(a_flux, *m_startFlux, a_coordSys.dx(),m_dx,m_verbose);
  if (a_time > m_startTime && m_startTime < m_endTime)
    {
      
      Real w = std::min(1.0 , (a_time - m_startTime) / (m_endTime - m_startTime)); 

      pout() << " LevelDataSurfaceFlux::m_endTime = " <<   m_endTime
	     << " w = " << w; 

      LevelData<FArrayBox> tmp; tmp.define(a_flux);
      for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
	{
	  tmp[dit].setVal(0.0);
	}

      FillFromReference(a_flux, *m_endFlux, a_coordSys.dx(),m_dx,m_verbose);
      for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
  	{
  	  tmp[dit] *=w;
  	  a_flux[dit] *= (1.0-w);
  	  a_flux[dit] += tmp[dit];
  	}
    }
  
  const ProblemDomain& domain = a_flux.disjointBoxLayout().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_flux, domain, dir, Side::Lo);
	ReflectGhostCells(a_flux, domain, dir, Side::Hi);
      }
    }
  a_flux.exchange();
  
  Real f = 0.0;
  for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
    {
      f = max(f, a_flux[dit].max());
    }
  pout() << " f = " << f << std::endl;

}



#include "NamespaceFooter.H"
