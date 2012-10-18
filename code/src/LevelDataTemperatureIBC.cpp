#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  LevelDataTemperatureIBC.cpp
// ============
//
// PhysIBC-derived class which stores initial temperature  data
// and imposes either periodic or reflection boundary conditions

#include "LevelDataTemperatureIBC.H"
#include "IceConstants.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "ReadLevelData.H"
#include "NamespaceHeader.H"


LevelDataTemperatureIBC* 
LevelDataTemperatureIBC::parse(ParmParse& a_pp)
{

  	std::string infile;
	a_pp.get("temperatureFile",infile);
	std::string temperatureName = "temp000000";
	a_pp.query("temperatureName",temperatureName);
	RefCountedPtr<LevelData<FArrayBox> > levelTemp
	  (new LevelData<FArrayBox>());
	Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
	vectData.push_back(levelTemp);
	Vector<std::string> names(1);
	names[0] = temperatureName;
	Real dx;
	ParmParse ppAmr ("amr");
	Vector<int> ancells(3); 
	ppAmr.getarr("num_cells", ancells, 0, ancells.size());
	readLevelData(vectData,dx,infile,names,ancells[2]);
	RealVect levelDx = RealVect::Unit * dx;

	RefCountedPtr<LevelData<FArrayBox> > levelSurfaceTemp(new LevelData<FArrayBox>());
	std::string surfaceTemperatureName = "";
	a_pp.query("surfaceTemperatureName",surfaceTemperatureName);
	if (surfaceTemperatureName == "")
	  {
	    //not ideal, but in this case copy the top layer temperature to the surface.
	    levelSurfaceTemp->define( levelTemp->disjointBoxLayout(), 1, levelTemp->ghostVect());
	    levelTemp->copyTo(Interval(0,0),*levelSurfaceTemp,Interval(0,0));
	  }
	else
	  {
	    names.resize(1);
	    names[0] = surfaceTemperatureName;
	    vectData[0] = levelSurfaceTemp;
	    readLevelData(vectData,dx,infile,names,1);
	    ; 
	    if (dx != levelDx[0])
	      {
		pout() << "surface temperature dx = " << dx << " but bulk temperature mesh dx = " << levelDx[0] << endl;
		CH_assert(dx == levelDx[0]);
		MayDay::Error("dx != levelDx[0]");
	      }
	  }

	RefCountedPtr<LevelData<FArrayBox> > levelBasalTemp(new LevelData<FArrayBox>());
	std::string basalTemperatureName = "";
	a_pp.query("basalTemperatureName",basalTemperatureName);
	if (basalTemperatureName == "")
	  {
	   
	    //not ideal, but in this case copy the bottome layer temperature to the base
	    levelBasalTemp->define( levelTemp->disjointBoxLayout(), 1, levelTemp->ghostVect());
	    const int n = levelTemp->nComp() - 1;
	    levelTemp->copyTo(Interval(n,n),*levelBasalTemp,Interval(0,0));
	  }
	else
	  {
	    names.resize(1);
	    names[0] = basalTemperatureName;
	    vectData[0] = levelBasalTemp;
	    readLevelData(vectData,dx,infile,names,1);
	    if (dx != levelDx[0])
	      {
		pout() << "basal temperature dx = " << dx << " but bulk temperature mesh dx = " << levelDx[0] << endl;
		CH_assert(dx == levelDx[0]);
		MayDay::Error("dx != levelDx[0]");
	      }
	  }
	return new LevelDataTemperatureIBC
	  (levelTemp,levelSurfaceTemp,levelBasalTemp,levelDx);

}

LevelDataTemperatureIBC::LevelDataTemperatureIBC
(RefCountedPtr<LevelData<FArrayBox> > a_temp, 
#if BISICLES_Z == BISICLES_LAYERED
 RefCountedPtr<LevelData<FArrayBox> > a_surfaceTemp,
 RefCountedPtr<LevelData<FArrayBox> > a_basalTemp,
#endif
 const RealVect& a_dx)
{
  m_temp = a_temp;
#if BISICLES_Z == BISICLES_LAYERED
  m_surfaceTemp = a_surfaceTemp;
  m_basalTemp = a_basalTemp;
#endif
  m_dx = a_dx;
  
  for (DataIterator dit( m_temp->disjointBoxLayout());dit.ok();++dit)
    {
      Real bulkMaxTemperature = (*m_temp)[dit].max();
      CH_assert(bulkMaxTemperature <= triplepoint);
   
      Real bulkMinTemperature = (*m_temp)[dit].min();
      CH_assert(bulkMinTemperature > 0);
    }
  for (DataIterator dit( m_surfaceTemp->disjointBoxLayout());dit.ok();++dit)
    {
      Real surfaceMaxTemperature = (*m_surfaceTemp)[dit].max();
      CH_assert(surfaceMaxTemperature <= triplepoint);
     
      Real surfaceMinTemperature = (*m_surfaceTemp)[dit].min();
      CH_assert(surfaceMinTemperature > 0);
    }
  for (DataIterator dit( m_basalTemp->disjointBoxLayout());dit.ok();++dit)
    {
      Real basalMaxTemperature = (*m_basalTemp)[dit].max();
      CH_assert(basalMaxTemperature <= triplepoint);
      
      Real basalMinTemperature = (*m_basalTemp)[dit].min();
      CH_assert(basalMinTemperature > 0);
    }
}

LevelDataTemperatureIBC::~LevelDataTemperatureIBC()
{
  
}

void LevelDataTemperatureIBC::define(const ProblemDomain& a_domain,
				     const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

#if BISICLES_Z == BISICLES_LAYERED
void LevelDataTemperatureIBC::initializeIceTemperature(LevelData<FArrayBox>& a_T, 
			      LevelData<FArrayBox>& a_surfaceT, 
			      LevelData<FArrayBox>& a_basalT,
			      const LevelSigmaCS& a_coordSys)
{

  if (true)
    {
      pout() << " LevelDataIBC::initializeIceTemperature" << endl;
    }

  FillFromReference(a_T,*m_temp,a_coordSys.dx(),m_dx,true);
  FillFromReference(a_surfaceT,*m_surfaceTemp,a_coordSys.dx(),m_dx,true);
  FillFromReference(a_basalT,*m_basalTemp,a_coordSys.dx(),m_dx,true);
  
  const ProblemDomain& domain = a_coordSys.grids().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_T, domain, dir, Side::Lo);
	ReflectGhostCells(a_T, domain, dir, Side::Hi);
	ReflectGhostCells(a_surfaceT, domain, dir, Side::Lo);
	ReflectGhostCells(a_surfaceT, domain, dir, Side::Hi);
	ReflectGhostCells(a_basalT, domain, dir, Side::Lo);
	ReflectGhostCells(a_basalT, domain, dir, Side::Hi);
      }
    }
}

#elif BISICLES_Z == BISICLES_FULLZ
void LevelDataTemperatureIBC::initializeIceTemperature(LevelData<FArrayBox>& a_T)
{
  FillFromReference(a_T,*m_temp,a_coordSys.dx(),m_dx,true);
}
#endif



LevelDataTemperatureIBC* 
LevelDataTemperatureIBC::new_temperatureIBC()
{
#if BISICLES_Z == BISICLES_LAYERED
  return new LevelDataTemperatureIBC(m_temp,m_surfaceTemp,m_basalTemp,m_dx);
#elif BISICLES_Z == BISICLES_FULLZ
  return new LevelDataTemperatureIBC(m_temp,m_dx);
#endif
}




#include "NamespaceFooter.H"
