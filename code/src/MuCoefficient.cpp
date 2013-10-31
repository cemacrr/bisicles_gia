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
#include "ParmParse.H"
#include "ReadLevelData.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "NamespaceHeader.H"

///assemble a MuCoefficient  object from ParmParse input, return pointer
MuCoefficient* MuCoefficient::parseMuCoefficient(const char* a_prefix)
{
  MuCoefficient* ptr = NULL;
  ParmParse muPP(a_prefix);
  std::string muCoefType = "unit";
  muPP.query("type",muCoefType );
  if (muCoefType == "unit")
    {
      ptr = static_cast<MuCoefficient*>(new UnitMuCoefficient());
    }
  else if (muCoefType == "LevelData")
    {
      //read a one level muCoef from an AMR Hierarchy, and  store it in a LevelDataMuCoeffcient
      ParmParse ildPP("inputLevelData");
      std::string infile;
      ildPP.get("muCoefFile",infile);
      std::string frictionName = "muCoef";
      ildPP.query("muCoefName",frictionName);
      RefCountedPtr<LevelData<FArrayBox> > levelMuCoef (new LevelData<FArrayBox>());
      Vector<RefCountedPtr<LevelData<FArrayBox> > > vectMuCoef;
      vectMuCoef.push_back(levelMuCoef);
      Vector<std::string> names(1);
      names[0] = frictionName;
      Real dx;
      readLevelData(vectMuCoef,dx,infile,names,1);
      RealVect levelDx = RealVect::Unit * dx;
      ptr = static_cast<MuCoefficient*>
	(new LevelDataMuCoefficient(levelMuCoef,levelDx));
      
      }
    else if (muCoefType == "MultiLevelData")
      {
	//read a multi level muCoef from an AMR Hierarchy, and  store it in a MultiLevelDataMuCoeffcient
	 ParmParse ildPP("inputLevelData");
	 std::string infile;
	 ildPP.get("muCoefFile",infile);
	 std::string muCoefName = "muCoef";
	 ildPP.query("muCoefName",muCoefName);
	 
	 Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectMuCoef;
	 Vector<std::string> names(1);
	 names[0] = muCoefName;
	 Real dx;
	 Vector<int> ratio;
	 readMultiLevelData(vectMuCoef,dx,ratio,infile,names,1);
	 RealVect dxCrse = RealVect::Unit * dx;
	 ptr = static_cast<MuCoefficient*>
	   (new MultiLevelDataMuCoefficient(vectMuCoef[0],dxCrse,ratio));
	 
      }
#ifdef HAVE_PYTHON
    else if (muCoefType == "Python")
      {
	ParmParse pyPP("PythonMuCoefficient");
	std::string module;
	pyPP.get("module",module);
	std::string funcName = "muCoefficient";
	pyPP.query("function",funcName);
	ptr =  static_cast<MuCoefficient*>
	  (new PythonInterface::PythonMuCoefficient(module, funcName));
      }
#endif
    else
      {
	MayDay::Error("undefined MuCoefficient in inputs");
      }
  
  return ptr;
  
}


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


MuCoefficient* AxbyMuCoefficient::new_muCoefficient() const
{
  return static_cast<MuCoefficient*>
    (new AxbyMuCoefficient(m_a, m_x,m_b, m_y) );
}

AxbyMuCoefficient::AxbyMuCoefficient
(const Real& a_a, MuCoefficient* a_x, 
 const Real& a_b, MuCoefficient* a_y)
{

  m_a = a_a;
  m_b = a_b;
  
  CH_assert(a_x != NULL);
  CH_assert(a_y != NULL);
  m_x = a_x->new_muCoefficient();
  m_y = a_y->new_muCoefficient();
  CH_assert(m_x != NULL);
  CH_assert(m_y != NULL);

}

AxbyMuCoefficient::~AxbyMuCoefficient()
{
  if (m_x != NULL)
    {
      delete m_x; m_x = NULL;
    }
  if (m_y != NULL)
    {
      delete m_y; m_y = NULL;
    }
}

void AxbyMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{

  LevelData<FArrayBox> y_mu(a_cellMuCoef.disjointBoxLayout(),1,a_cellMuCoef.ghostVect());
  m_x->setMuCoefficient(a_cellMuCoef, a_coordSys, a_time, a_dt);
  m_y->setMuCoefficient(y_mu, a_coordSys, a_time, a_dt);
  for (DataIterator dit(a_cellMuCoef.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].axby(a_cellMuCoef[dit],y_mu[dit],m_a,m_b);
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
