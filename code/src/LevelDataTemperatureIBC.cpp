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
#include "IceThermodynamics.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "ReadLevelData.H"
#include "AmrIceBase.H"
#include "NamespaceHeader.H"


LevelDataTemperatureIBC* 
LevelDataTemperatureIBC::parse(ParmParse& a_pp)
{

  	std::string infile;
	a_pp.get("temperatureFile",infile);
	std::string temperatureName = "temp000000";
	a_pp.query("temperatureName",temperatureName);
	Real defaultTemperature = 258.0; // 
	a_pp.query("value", defaultTemperature);
	
	RefCountedPtr<LevelData<FArrayBox> > levelTemp
	  (new LevelData<FArrayBox>());
	Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
	vectData.push_back(levelTemp);
	Vector<std::string> names(1);
	names[0] = temperatureName;
	Real dx;
	ParmParse ppAmr ("amr");
	Vector<int> ancells(3,0); 
	ppAmr.queryarr("num_cells", ancells, 0, ancells.size());
	if (ancells[0] == 0)
	  {
	    ParmParse ppGeo ("geometry");
	    ppGeo.getarr("num_cells", ancells, 0, ancells.size());
	  }
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

	RefCountedPtr<LevelData<FArrayBox> > levelBasalHeatFlux(new LevelData<FArrayBox>());
	std::string basalHeatFluxName = "";
	a_pp.query("basalHeatFluxName",basalHeatFluxName);
	if (basalHeatFluxName == "")
	  {
	    //if no basal heat flux is given, assume zero flux
	    levelBasalHeatFlux->define( levelTemp->disjointBoxLayout(), 1, levelTemp->ghostVect());
	    for (DataIterator dit = levelBasalHeatFlux->dataIterator(); dit.ok(); ++dit)
	      {
		(*levelBasalHeatFlux)[dit].setVal(0.0);
	      }
	  }
	else
	  {
	    names.resize(1);
	    names[0] = basalHeatFluxName;
	    vectData[0] = levelBasalHeatFlux;
	    readLevelData(vectData,dx,infile,names,1);
	    if (dx != levelDx[0])
	      {
		pout() << "basal temperature dx = " << dx << " but bulk temperature mesh dx = " << levelDx[0] << endl;
		CH_assert(dx == levelDx[0]);
		MayDay::Error("dx != levelDx[0]");
	      }
	  }
	return new LevelDataTemperatureIBC
	  (levelTemp,levelSurfaceTemp,levelBasalHeatFlux,levelDx, defaultTemperature);

}

LevelDataTemperatureIBC::LevelDataTemperatureIBC
(RefCountedPtr<LevelData<FArrayBox> > a_temp, 
 RefCountedPtr<LevelData<FArrayBox> > a_surfaceTemp,
 RefCountedPtr<LevelData<FArrayBox> > a_basalHeatFlux,
 const RealVect& a_dx, const Real& a_defaultTemperature)
{
  m_temp = a_temp;
  m_surfaceTemp = a_surfaceTemp;
  m_basalHeatFlux = a_basalHeatFlux;
  m_defaultTemperature =  a_defaultTemperature;
  m_dx = a_dx;
  
  for (DataIterator dit( m_temp->disjointBoxLayout());dit.ok();++dit)
    {
      Real bulkMaxTemperature = (*m_temp)[dit].max();
      //CH_assert(bulkMaxTemperature <= triplepoint);
   
      Real bulkMinTemperature = (*m_temp)[dit].min();
      CH_assert(bulkMinTemperature > 0);
    }
  for (DataIterator dit( m_surfaceTemp->disjointBoxLayout());dit.ok();++dit)
    {
      Real surfaceMaxTemperature = (*m_surfaceTemp)[dit].max();
      //CH_assert(surfaceMaxTemperature <= triplepoint);
     
      Real surfaceMinTemperature = (*m_surfaceTemp)[dit].min();
      CH_assert(surfaceMinTemperature > 0);
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

void LevelDataTemperatureIBC::basalHeatFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  FillFromReference(a_flux,*m_basalHeatFlux,a_amrIce.dx(a_level),m_dx,true);
  const ProblemDomain& domain = a_amrIce.grids(a_level).physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_flux, domain, dir, Side::Lo);
	ReflectGhostCells(a_flux, domain, dir, Side::Hi);
      }
    }

}

void LevelDataTemperatureIBC::initializeIceInternalEnergy(LevelData<FArrayBox>& a_E,
							  LevelData<FArrayBox>& a_tillWaterDepth,
							  LevelData<FArrayBox>& a_surfaceE, 
							  LevelData<FArrayBox>& a_basalE,
							  const AmrIceBase& a_amrIce, 
							  int a_level, Real a_dt)
{

  if (true)
    {
      pout() << " LevelDataIBC::initializeIceTemperature" << endl;
    }

  
  const LevelSigmaCS& coordSys = *a_amrIce.geometry(a_level); 
  const DisjointBoxLayout dbl = coordSys.grids();
  LevelData<FArrayBox> T(dbl, a_E.nComp(), a_E.ghostVect());
  LevelData<FArrayBox> sT(dbl, 1, a_E.ghostVect());

  // set default tempeature - applies to regions of the domain not
  // covered by the data
  for (DataIterator dit(dbl);dit.ok();++dit)
    {
      T[dit].setVal(m_defaultTemperature);
      sT[dit].setVal(m_defaultTemperature);
    }
  
  FillFromReference(T,*m_temp,coordSys.dx(),m_dx,true);
  FillFromReference(sT,*m_surfaceTemp,coordSys.dx(),m_dx,true);
  {
  const ProblemDomain& domain = coordSys.grids().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(T, domain, dir, Side::Lo);
	ReflectGhostCells(T, domain, dir, Side::Hi);
	ReflectGhostCells(sT, domain, dir, Side::Lo);
	ReflectGhostCells(sT, domain, dir, Side::Hi);
      }
    }
  }

  for (DataIterator dit(dbl);dit.ok();++dit)
    {
      FArrayBox w( a_E[dit].box(), a_E[dit].nComp()); //water fraction, set to zero for now
      w.setVal(0.0);
      IceThermodynamics::composeInternalEnergy(a_E[dit],T[dit],w, a_E[dit].box() );
    }

 
   for (DataIterator dit(dbl);dit.ok();++dit)
    {
      FArrayBox w( a_surfaceE[dit].box(), 1); //water fraction, set to zero for now 
      w.setVal(0.0);
      IceThermodynamics::composeInternalEnergy(a_surfaceE[dit],sT[dit],w,a_surfaceE[dit].box() );
    }
   
  for (DataIterator dit(dbl);dit.ok();++dit)
    {
      a_tillWaterDepth[dit].setVal(0.0);
    }
   
  const ProblemDomain& domain = coordSys.grids().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_E, domain, dir, Side::Lo);
	ReflectGhostCells(a_E, domain, dir, Side::Hi);
	ReflectGhostCells(a_surfaceE, domain, dir, Side::Lo);
	ReflectGhostCells(a_surfaceE, domain, dir, Side::Hi);
      }
    }

  a_E.exchange();
  a_surfaceE.exchange();

}


LevelDataTemperatureIBC* 
LevelDataTemperatureIBC::new_internalEnergyIBC()
{
  return new LevelDataTemperatureIBC(m_temp,m_surfaceTemp,m_basalHeatFlux,m_dx,m_defaultTemperature);
}




#include "NamespaceFooter.H"
