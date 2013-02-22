#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "MuCoefficient.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "NamespaceHeader.H"


MuCoefficient* UnitMuCoefficient::new_muCoefficient() const
{
  UnitMuCoefficient* ptr = new UnitMuCoefficient();
  return static_cast<MuCoefficient*>(ptr);

}

void UnitMuCoefficient::setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
					 LevelSigmaCS& a_coordSys,
					 Real a_time,
					 Real a_dt)
{
  for (DataIterator dit = a_cellMuCoef.dataIterator(); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].setVal(1.0);
    }
}

MuCoefficient* LevelDataMuCoefficient::new_muCoefficient() const
{
  LevelDataMuCoefficient* ptr = new LevelDataMuCoefficient(m_muCoef, m_dx);
  return static_cast<MuCoefficient*>(ptr);

}
void LevelDataMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{

  for (DataIterator dit = a_cellMuCoef.dataIterator(); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].setVal(1.0);
      
    }

  FillFromReference(a_cellMuCoef, *m_muCoef, a_coordSys.dx(),m_dx,m_verbose);
  
   for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const ProblemDomain domain = a_cellMuCoef.disjointBoxLayout().physDomain();
      if (!(domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_cellMuCoef, domain, dir, Side::Lo);
	  ReflectGhostCells(a_cellMuCoef, domain, dir, Side::Hi);
	}
    }
}

MuCoefficient* MultiLevelDataMuCoefficient::new_muCoefficient() const
{
  MultiLevelDataMuCoefficient* ptr = new MultiLevelDataMuCoefficient(m_muCoef, m_dxCrse, m_ratio);
  return static_cast<MuCoefficient*>(ptr);

}


void MultiLevelDataMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt) 
{
    
  
    setMuCoefficient(a_cellMuCoef,a_coordSys.dx(),a_time,a_dt); 
  
}

void MultiLevelDataMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 RealVect a_dx,
 Real a_time,
 Real a_dt)
{


  RealVect dx(m_dxCrse);
  for (int refDataLev = 0; refDataLev < m_muCoef.size(); refDataLev++)
    {
      
      if (a_dx[0] < dx[0] && refDataLev > 0)
	{
	  //We need a LevelData<FArrayBox> whose DisjointBoxLayout covers
	  //a_muCoef's, but m_muCoef[refDataLev] doesn't necessarily provide one. Since we secretly
	  //want to be lisp programmers, we turn to recursion to build one,
	  //throwing performance worries to the wind.
	  //\todo consider modifying the BasalFriction interface to avoid recompuing levels 1 to n-1
	  const int& nRef = int(dx[0]/a_dx[0]);
	  DisjointBoxLayout crseDBL;
	  coarsen(crseDBL, a_cellMuCoef.disjointBoxLayout(), nRef);
	  LevelData<FArrayBox> crseMuCoef(crseDBL, 1 , IntVect::Unit);
	  this->setMuCoefficient(crseMuCoef , dx, a_time, a_dt);
	  FillFromReference(a_cellMuCoef, crseMuCoef , a_dx, dx ,m_verbose);
	  
	}
      else
	{
	  FillFromReference(a_cellMuCoef, *m_muCoef[refDataLev], a_dx, dx ,m_verbose);
	}
      dx /= Real(m_ratio[refDataLev]);
    }

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const ProblemDomain domain = a_cellMuCoef.disjointBoxLayout().physDomain();
      if (!(domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_cellMuCoef, domain, dir, Side::Lo);
	  ReflectGhostCells(a_cellMuCoef, domain, dir, Side::Hi);
	}
    }

  
}


#include "NamespaceFooter.H"
