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
#include "NamespaceHeader.H"

LevelDataTemperatureIBC::LevelDataTemperatureIBC
(RefCountedPtr<LevelData<FArrayBox> > a_temp, 
 const RealVect& a_dx)
{
  m_temp = a_temp;
  m_dx = a_dx;
  // Construction means nothing to me. It's a lot like Vienna.
  for (DataIterator dit( m_temp->disjointBoxLayout());dit.ok();++dit)
    {
      CH_assert( (*m_temp)[dit].min() <= triplepoint);
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
  
  const ProblemDomain& domain = a_coordSys.grids().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_T, domain, dir, Side::Lo);
	ReflectGhostCells(a_T, domain, dir, Side::Hi);
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
  return new LevelDataTemperatureIBC(m_temp,m_dx);
}




#include "NamespaceFooter.H"
