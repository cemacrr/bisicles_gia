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
#include "CellToEdge.H"
#include "FillFromReference.H"
#include "NamespaceHeader.H"


MuCoefficient* UnitMuCoefficient::new_muCoefficient() const
{
  UnitMuCoefficient* ptr = new UnitMuCoefficient();
  return static_cast<MuCoefficient*>(ptr);

}

void UnitMuCoefficient::setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
					 LevelData<FluxBox>& a_faceMuCoef,
					 LevelSigmaCS& a_coordSys,
					 Real a_time,
					 Real a_dt)
{
  for (DataIterator dit = a_cellMuCoef.dataIterator(); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].setVal(1.0);
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  a_faceMuCoef[dit][dir].setVal(1.0);
	}
    }
}

MuCoefficient* LevelDataMuCoefficient::new_muCoefficient() const
{
  LevelDataMuCoefficient* ptr = new LevelDataMuCoefficient(m_muCoef, m_dx);
  return static_cast<MuCoefficient*>(ptr);

}
void LevelDataMuCoefficient::setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
					 LevelData<FluxBox>& a_faceMuCoef,
					 LevelSigmaCS& a_coordSys,
					 Real a_time,
					 Real a_dt)
{

  for (DataIterator dit = a_cellMuCoef.dataIterator(); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].setVal(1.0);
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  a_faceMuCoef[dit][dir].setVal(1.0);
	}
    }

  FillFromReference(a_cellMuCoef, *m_muCoef, a_coordSys.dx(),m_dx,m_verbose);
  CellToEdge(a_cellMuCoef, a_faceMuCoef);
}

MuCoefficient* MultiLevelDataMuCoefficient::new_muCoefficient() const
{
  MultiLevelDataMuCoefficient* ptr = new MultiLevelDataMuCoefficient(m_muCoef, m_dxCrse, m_ratio);
  return static_cast<MuCoefficient*>(ptr);

}
void MultiLevelDataMuCoefficient::setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
					 LevelData<FluxBox>& a_faceMuCoef,
					 LevelSigmaCS& a_coordSys,
					 Real a_time,
					 Real a_dt)
{

  for (DataIterator dit = a_cellMuCoef.dataIterator(); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].setVal(1.0);
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  a_faceMuCoef[dit][dir].setVal(1.0);
	}
    }
  RealVect dx(m_dxCrse);
  for (int refDataLev = 0; refDataLev < m_muCoef.size(); refDataLev++)
    {
      FillFromReference(a_cellMuCoef, *m_muCoef[refDataLev], a_coordSys.dx(), dx ,m_verbose);
      dx /= Real(m_ratio[refDataLev]);
    }
  CellToEdge(a_cellMuCoef, a_faceMuCoef);
}


#include "NamespaceFooter.H"
