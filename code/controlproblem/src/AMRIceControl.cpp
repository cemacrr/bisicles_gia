#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//================================================================
// AMRIceControl.cpp
// methods to contruct and solve ice sheet velocity problems
// that do not involve time-stepping
// 
// Primarily aimed at solving a  control problem (similar to 
// that of Vieli and Payne, Ann Glaciol, 2003. vol 26)
//
//=================================================================

#include "AMRIceControl.H"
#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "CGOptimize.H"
#include "FillFromReference.H"
#include "ControlF_F.H"
#include "BisiclesF_F.H"
#include "ReflectGhostCells.H"
#include "QuadCFInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "computeNorm.H"
#include "computeSum.H"
#include "IceVelocity.H"
#include "IceConstants.H"
#include "JFNKSolver.H"
#include <sstream>
#include "NamespaceHeader.H"

void AMRIceControl::define(IceThicknessIBC* a_ibcPtr,
			   IceTemperatureIBC* a_tempIBCPtr,
			   RateFactor* a_rateFactor,
			   ConstitutiveRelation* a_constRelPtr, 
			   BasalFrictionRelation* a_bfRelPtr,
			   const Vector<RealVect>& a_dataDx,
			   const Vector<RefCountedPtr<LevelData<FArrayBox > > >& a_referenceCOrigin,
			   const Vector<RefCountedPtr<LevelData<FArrayBox > > >& a_referenceMuCoefOrigin,
			   const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceVelObs,
			   const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceVelCoef,
			   const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceDivUHObs,
			   const Vector<RefCountedPtr<LevelData<FArrayBox > > >& a_referenceDivUHCoef)
{
  m_dataDx = a_dataDx;
  m_referenceCOrigin = a_referenceCOrigin;
  m_referenceMuCoefOrigin =  a_referenceMuCoefOrigin;
  m_referenceVelObs = a_referenceVelObs;
  m_referenceVelCoef = a_referenceVelCoef;
  m_referenceDivUHObs = a_referenceDivUHObs;
  m_referenceDivUHCoef =  a_referenceDivUHCoef;
  m_ibcPtr = a_ibcPtr;
  m_tempIBCPtr = a_tempIBCPtr;
  m_rateFactor = a_rateFactor;
  m_constRelPtr = a_constRelPtr;
  m_bfRelPtr = a_bfRelPtr;

  //base level geometry set up
  m_finestLevel = 0;
  
  ParmParse ppGeom("geometry");
  Vector<int> ncells(3); 
  ppGeom.getarr("num_cells",ncells,0,ncells.size());
  bool periodic[SpaceDim] = {D_DECL(false,false,false)};
  {
    Vector<int> is_periodic(SpaceDim);
    ppGeom.getarr("is_periodic", is_periodic, 0, SpaceDim);
    for (int dir=0; dir<SpaceDim; dir++) 
      {
	periodic[dir] = (is_periodic[dir] == 1);
      }
  }
  
  IntVect lo(IntVect::Zero);
  IntVect hi(D_DECL(ncells[0]-1, ncells[1]-1, ncells[2]-1));
  ProblemDomain domain(lo,hi,periodic);
  RealVect domainSize;
  Vector<Real> domSize(SpaceDim);
  ppGeom.getarr("domain_size", domSize, 0, SpaceDim);
  domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));


  
  m_faceSigma.resize(ncells[2]+1);
  Real dsigma = 1.0 / ncells[2];
  m_faceSigma[0] = 0.0;
  for (int j = 1; j < m_faceSigma.size()-1; ++j)
    m_faceSigma[j] = m_faceSigma[j-1] + dsigma;
  m_faceSigma[m_faceSigma.size()-1] = 1.0;

  ppGeom.queryarr("sigma",m_faceSigma,0,m_faceSigma.size());

  CH_assert(abs(m_faceSigma[0] - 0.0) < 1.0e-6);
  CH_assert(abs(m_faceSigma[m_faceSigma.size()-1] - 1.0) < 1.0e-6);

  
  ParmParse ppAMR("amr");
  m_blockFactor = 16;
  ppAMR.query("block_factor",m_blockFactor);
  m_maxBoxSize = 64;
  ppAMR.query("max_box_size",m_maxBoxSize);
  m_nestingRadius = 1;
  ppAMR.query("nesting_radius",m_nestingRadius);
  m_fillRatio = 0.85;
  ppAMR.query("fill_ratio",m_fillRatio);
  m_tagsGrow = 0;
  ppAMR.query("tags_grow",m_tagsGrow);

  m_maxFinestLevel = 0;
  ppAMR.query("max_level",m_maxFinestLevel);
  m_refRatio.resize(m_maxFinestLevel);
  m_maxFinestLevelFloating = m_maxFinestLevel;
  m_maxFinestLevelGrounded = m_maxFinestLevel;

  ppAMR.query("max_level_floating",m_maxFinestLevelFloating);
  ppAMR.query("max_level_grounded",m_maxFinestLevelGrounded);
   

  if (m_maxFinestLevel > 0)
    {
      ppAMR.queryarr("refinement_ratios",m_refRatio,0,m_maxFinestLevel);
    }

  m_dx.resize(m_maxFinestLevel + 1);
  m_dx[0] = domainSize;
  for (int dir = 0; dir < SpaceDim; ++dir)
    m_dx[0][dir] /= Real(ncells[dir]);

  m_domain.resize(m_maxFinestLevel + 1);
  m_domain[0] = domain;
  for (int lev=1;lev <= m_maxFinestLevel;lev++)
    {
      m_dx[lev] = m_dx[lev-1]/Real(m_refRatio[lev-1]);
      m_domain[lev] = refine(m_domain[lev-1],m_refRatio[lev-1]);
    }

  Vector<Box> baseBoxes;
  domainSplit(domain, baseBoxes,m_maxBoxSize,m_blockFactor);
  Vector<int> procAssign(baseBoxes.size());
  LoadBalance(procAssign,baseBoxes);
  
  m_grids.resize( m_maxFinestLevel + 1 );
  m_grids[0].define(baseBoxes, procAssign, domain);

  //resize all the AMR Hierarchy vectors
  m_coordSys.resize(m_maxFinestLevel + 1);
  //m_velObs.resize(m_maxFinestLevel + 1,NULL);
  //m_velCoef.resize(m_maxFinestLevel + 1,NULL);
  //m_divUHObs.resize(m_maxFinestLevel + 1,NULL);
  //m_divUHCoef.resize(m_maxFinestLevel + 1,NULL);
  m_maxVelDx = 3.0e6;
  ppAMR.query("max_vel_dx",m_maxVelDx);
  //base level data setup
  levelSetup(0);
  

  //mesh refinement
  for (int fl = 0; fl < m_maxFinestLevel; fl++)
    {
      Vector<IntVectSet> tagVect(m_finestLevel+1);
      
      for (int lev = 0; lev <= m_finestLevel; ++lev)
	{
	  levelTag(lev,tagVect[lev]);
	}
      
      BRMeshRefine meshrefine(m_domain[0], m_refRatio,
			      m_fillRatio, m_blockFactor, 
			      m_nestingRadius, m_maxBoxSize);
      
      //turn current grids into a Vector<Box> 
      Vector<Vector<Box> > oldGrids(m_finestLevel+1);
      Vector<Vector<Box> > newGrids(m_finestLevel+1);
      for (int lev=0; lev<= m_finestLevel; lev++)
	{
	  const DisjointBoxLayout& levelDBL = m_grids[lev];
	  oldGrids[lev].resize(levelDBL.size());
	  LayoutIterator lit = levelDBL.layoutIterator();
	  int boxIndex = 0;
	  for (lit.begin(); lit.ok(); ++lit, ++boxIndex) 
	    {
	      oldGrids[lev][boxIndex] = levelDBL[lit()];
	    }
	}
      m_finestLevel = meshrefine.regrid
	(newGrids, tagVect,  0 , m_finestLevel, oldGrids);
      

      for (int lev=0; lev<= m_finestLevel; lev++)
	{
	  int numGridsNew = newGrids[lev].size();
	  Vector<int> procIDs(numGridsNew);
	  LoadBalance(procIDs, newGrids[lev]);
	  const DisjointBoxLayout newDBL(newGrids[lev], procIDs,
					 m_domain[lev]);
	  
	  m_grids[lev] = newDBL;
	  levelSetup(lev);
	}
    } 

  
  create(m_velb,SpaceDim,IntVect::Unit);
  setToZero(m_velb);
  create(m_vels,SpaceDim,IntVect::Unit);
  setToZero(m_vels);

  create(m_adjVel,SpaceDim,IntVect::Unit);

  create(m_rhs,SpaceDim,IntVect::Unit);
  IceVelocity::defineRHS(m_rhs, m_coordSys, m_grids, m_dx);
  
  create(m_adjRhs,SpaceDim,IntVect::Unit);

  create(m_C,1,IntVect::Unit);
  create(m_Ccopy,1,IntVect::Unit);
  create(m_C0,1,IntVect::Zero);
  setToZero(m_C0);
  create(m_lapC,1,IntVect::Zero);
  create(m_gradCSq,1,IntVect::Zero);
  create(m_muCoef,1,IntVect::Unit);
  create(m_lapMuCoef,1,IntVect::Zero);
  create(m_gradMuCoefSq,1,IntVect::Zero);
  create(m_muCoef,1,IntVect::Unit);
  create(m_faceMuCoef,1,IntVect::Zero);

  create(m_divUH,1,IntVect::Unit);
  create(m_faceThckFlux,1,IntVect::Unit);
  
  create(m_temperature,m_faceSigma.size()-1,IntVect::Unit);
  create(m_sTemperature,1,IntVect::Unit);
  create(m_bTemperature,1,IntVect::Unit);
  create(m_A,m_faceSigma.size()-1,IntVect::Unit);
  create(m_faceA,m_faceSigma.size()-1,IntVect::Unit);

  for (int lev=0;lev<=m_finestLevel;lev++)
    {
       m_tempIBCPtr->initializeIceTemperature
	 (*m_temperature[lev], *m_sTemperature[lev], *m_bTemperature[lev],*m_coordSys[lev]);

      IceVelocity::computeA(*m_A[lev],m_coordSys[lev]->getSigma(), *m_coordSys[lev], m_rateFactor, *m_temperature[lev]);
      CellToEdge(*m_A[lev], *m_faceA[lev]);
    }
  
#if 0
  {
    Vector<std::string> names;
    for (int i =0; i < m_A[0]->nComp(); i++)
      {
	char s[64]; sprintf(s,"T%i",i);
	names.push_back(s);
      }
    WriteAMRHierarchyHDF5("T.2d.hdf5",m_grids,m_temperature,names, m_domain[0].domainBox(),
			  m_dx[0][0], 0.0, 0.0 , m_refRatio, m_A.size());
  }

  {
    Vector<std::string> names;
    for (int i =0; i < m_A[0]->nComp(); i++)
      {
	char s[64]; sprintf(s,"A%i",i);
	names.push_back(s);
      }
    WriteAMRHierarchyHDF5("A.2d.hdf5",m_grids,m_A,names, m_domain[0].domainBox(),
			  m_dx[0][0], 0.0, 0.0 , m_refRatio, m_A.size());
  }
#endif

  Real Amin = computeMin(m_A,  m_refRatio, Interval(0,m_A[0]->nComp()-1));
  Real Amax = computeMax(m_A,  m_refRatio, Interval(0,m_A[0]->nComp()-1));
  pout() << Amin << " <= A(x,y,sigma) <= " << Amax << std::endl;

  create(m_vtFace, SpaceDim, IntVect::Unit);
 
  create(m_velocityMisfit,1,IntVect::Unit);
  create(m_massBalanceMisfit,1,IntVect::Unit);
  create(m_TikhonovPenalty,2,IntVect::Unit);
  create(m_barrierPenalty,2,IntVect::Unit);


  m_velocityInitialised = false;

  m_referenceCOrigin = NULL;
  m_referenceMuCoefOrigin = NULL; 
  m_referenceVelObs = NULL;
  m_referenceVelCoef = NULL;
  m_referenceDivUHObs =NULL;
  m_referenceDivUHCoef =  NULL;



  pout() << " AMIceControl::define completed" << std::endl;
  
}

void AMRIceControl::levelTag(int a_lev, IntVectSet& a_tagset)
{

  int lev = a_lev;
  const DisjointBoxLayout& grids = m_grids[lev];
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      
      const FArrayBox& velObs =  (*m_velObs[lev])[dit];
      const FArrayBox& velCoef =  (*m_velCoef[lev])[dit];
      const BaseFab<int>& mask = m_coordSys[lev]->getFloatingMask()[dit];
      Box box = grids[dit];

      //tag cells where |vel| * dx > m_maxVelDx 
      if (m_maxVelDx > 0.0)
	{
	  for (BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      const Real& vc = velCoef(iv);
	      Real tol = -0.001;
	      if ( (vc > tol) && (
		   (mask(iv) == GROUNDEDMASKVAL && lev < m_maxFinestLevelGrounded) ||
		   (mask(iv) == FLOATINGMASKVAL && lev < m_maxFinestLevelFloating) ) )
		{
		  Real v = 0.0;
		  for (int dir = 0; dir < SpaceDim; ++dir)
		    v += std::pow(velObs(iv,dir),2);
		  if (sqrt(v)*m_dx[lev][0] > m_maxVelDx)
		    a_tagset |= iv;
		}   
	    }
	}

      


    }

  a_tagset.grow(m_tagsGrow);
  
  // if ( a_lev == 0)
  //   {
  //     Box box(m_domain[0].domainBox());
  //     box.grow(-m_blockFactor);
  //     a_tagset &= box;
  //   }


}




void AMRIceControl::levelSetup(int a_lev)
{
  m_vectOps.resize(m_finestLevel+1);

  //co-ordinate system
  m_coordSys[a_lev] = RefCountedPtr<LevelSigmaCS>
    (new LevelSigmaCS(m_grids[a_lev], m_dx[a_lev], IntVect::Unit));
 
  int nRef = (a_lev > 0)?m_refRatio[a_lev-1]:0;
  LevelSigmaCS* crseCoordsPtr = (a_lev > 0)?(&(*m_coordSys[a_lev-1])):NULL;

  ParmParse ppCon("constants");
  Real iceDensity = 910;
  ppCon.query("ice_density",iceDensity); 
  m_coordSys[a_lev]->setIceDensity(iceDensity);
    
  Real seaWaterDensity=1028;
  ppCon.query("sea_water_density",seaWaterDensity);
  m_coordSys[a_lev]->setWaterDensity(seaWaterDensity);

  Real gravity=9.81;
  ppCon.query("gravity",gravity);
 m_coordSys[a_lev]->setGravity(gravity);

  m_ibcPtr->initializeIceGeometry
	(*m_coordSys[a_lev],m_dx[a_lev],RealVect::Zero, 0.0, crseCoordsPtr, nRef);

  m_ibcPtr->setGeometryBCs
    (*m_coordSys[a_lev], m_domain[a_lev],m_dx[a_lev],0.0,0.0);

  m_coordSys[a_lev]->recomputeGeometry(crseCoordsPtr, nRef);
  m_coordSys[a_lev]->setFaceSigma(m_faceSigma);

  //observations
  fill(m_velObs, m_referenceVelObs, m_dataDx, a_lev, IntVect::Unit);
  fill(m_velCoef, m_referenceVelCoef, m_dataDx, a_lev, IntVect::Unit);
  fill(m_divUHObs, m_referenceDivUHObs, m_dataDx, a_lev, IntVect::Unit);
  fill(m_divUHCoef, m_referenceDivUHCoef, m_dataDx, a_lev, IntVect::Unit);

  //initial C
  fill(m_COrigin, m_referenceCOrigin, m_dataDx , a_lev, IntVect::Unit);
  //initial muCoef
  fill(m_muCoefOrigin, m_referenceMuCoefOrigin, m_dataDx , a_lev, IntVect::Unit);
}


void AMRIceControl::fill
(Vector<LevelData<FArrayBox>*>& a_data, const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_refData, 
 const Vector<RealVect>& a_refDx, const int a_lev, const IntVect& a_ghost)
{

  if (a_data.size() <= a_lev)
    a_data.resize(a_lev+1);

  if (a_data[a_lev] != NULL)
    a_data[a_lev]->define(m_grids[a_lev], a_refData[0]->nComp(), a_ghost); 
  else
    a_data[a_lev] = new LevelData<FArrayBox>(m_grids[a_lev], a_refData[0]->nComp(), a_ghost);

  for (int refDataLev = 0; refDataLev < a_refData.size(); refDataLev++)
    {
      FillFromReference(*a_data[a_lev], *a_refData[refDataLev], m_dx[a_lev], a_refDx[refDataLev] ,true); 
    }
 
}


void AMRIceControl::fill
(Vector<LevelData<FArrayBox>*>& a_data, const LevelData<FArrayBox>& a_refData, 
 const int a_lev, const IntVect& a_ghost)
{

  MayDay::Error("this AMRIceControl::fill is to be removed");
  if (a_data.size() <= a_lev)
    a_data.resize(a_lev+1);

  if (a_data[a_lev] != NULL)
    a_data[a_lev]->define(m_grids[a_lev], a_refData.nComp(), a_ghost); 
  else
    a_data[a_lev] = new LevelData<FArrayBox>(m_grids[a_lev], a_refData.nComp(), a_ghost);

  FillFromReference(*a_data[a_lev], a_refData, m_dx[a_lev], m_dataDx[0] ,true);  
}

void AMRIceControl::solveControl()
{
  
  ParmParse pp("control");
  int CGmaxIter = 500;
  pp.query("CGmaxIter",CGmaxIter);
  Real CGtol = 1.0e-6;
  pp.query("CGtol",CGtol);
  Real CGsecantParameter = 1.0e-8;
  pp.query("CGsecantParameter",CGsecantParameter);
 
  
  Real CGsecantStepMaxGrow = -1.0;
  pp.query("CGsecantStepMaxGrow",CGsecantStepMaxGrow);

  int CGsecantMaxIter = 10;
  pp.query("CGsecantMaxIter",CGsecantMaxIter);

  Real CGsecantTol = 1.0e-6;
  pp.query("CGsecantTol",CGsecantTol);
    
  m_velMisfitCoefficient = 1.0;
  pp.query("velMisfitCoefficient",m_velMisfitCoefficient);

  m_massImbalanceCoefficient = 0.0;
  pp.query("massImbalanceCoefficient",m_massImbalanceCoefficient);

  m_gradCsqRegularization = 0.0;
  pp.query("gradCsqRegularization",m_gradCsqRegularization);

  m_gradMuCoefsqRegularization = m_gradCsqRegularization;
  pp.query("gradMuCoefsqRegularization",m_gradMuCoefsqRegularization);

  m_X0Regularization = 0.0;
  pp.query("X0Regularization",m_X0Regularization);

  m_X1Regularization = 0.0;
  pp.query("X1Regularization",m_X1Regularization);

  m_muCoefLEQOne = false;
  pp.query("muCoefLEQOne",m_muCoefLEQOne);
  m_muCoefLEQOneGround = false;
  pp.query("muCoefLEQOneGround",m_muCoefLEQOneGround);
  m_muCoefGEQOne = false;
  pp.query("muCoefGEQOne",m_muCoefGEQOne);
  m_muCoefGEQOneGround = false;
  pp.query("muCoefGEQOneGround",m_muCoefGEQOneGround);

  

  m_boundArgX0 = 3.0;
  m_boundArgX1 = 2.0;
  pp.query("boundArgX0",m_boundArgX0);
  pp.query("boundArgX1",m_boundArgX1);

  m_lowerX0 = - m_boundArgX0;
  m_upperX0 = + m_boundArgX0;
  pp.query("lowerX0",m_lowerX0);
  pp.query("upperX0",m_upperX0);

  m_lowerX1 = - m_boundArgX1;
  m_upperX1 = + m_boundArgX1;
  pp.query("lowerX1",m_lowerX1);
  pp.query("upperX1",m_upperX1);

  {
    std::string s = "None";
    pp.query("boundMethod",s);
    if (s == "None")
      {
	m_boundMethod = None;
      }
    else if (s == "Projection")
      {
	m_boundMethod = Projection;
      }
    else
      {
	MayDay::Error("unknown control.boundMethod");
      }
  }


  Vector<LevelData<FArrayBox>* > X;
  create(X,2,IntVect::Zero);
  setToZero(X);

  m_writeInnerSteps = false;
  pp.query("writeInnerSteps",m_writeInnerSteps);

  m_innerStepFileNameBase = "ControlInner";
  pp.query("innerStepFileNameBase",m_innerStepFileNameBase);
  
  m_outerStepFileNameBase = "ControlOuter";
  pp.query("outerStepFileNameBase",m_outerStepFileNameBase);

  m_outerCounter = 0;
  pp.query("restart",m_outerCounter);


  
  if (m_outerCounter > 0)
    {
      
      std::string restartfile = outerStateFile();
      pout() << "restarting from " << restartfile << std::endl;
      readState(restartfile, m_innerCounter, X);
      m_outerCounter++;
      m_velocityInitialised = true;      
    }
  else
    {
      m_innerCounter = 0;
    }


  //ensure X sits within constraints
  checkBounds(X);

  CGOptimize(*this ,  X , CGmaxIter , CGtol , 
	     CGsecantParameter, CGsecantStepMaxGrow, 
	     CGsecantMaxIter , CGsecantTol);

  free(X);
  

}

void AMRIceControl::solveForward()
{
  MayDay::Error("forward problem not yet implemented");
}

AMRIceControl::~AMRIceControl()
{
  free(m_velb);
  free(m_vels);
  free(m_velObs);
  free(m_velCoef);
  free(m_divUHObs);
  free(m_divUHCoef);
  free(m_adjVel);
  free(m_rhs);
  free(m_C);
  free(m_Ccopy);
  free(m_C0);
  free(m_lapC);
  free(m_gradCSq);
  free(m_muCoef);
  free(m_lapMuCoef);
  free(m_gradMuCoefSq);
  free(m_faceMuCoef);
  free(m_divUHObs);
  free(m_faceThckFlux);
  free(m_temperature);
  free(m_sTemperature);
  free(m_bTemperature);
  free(m_A);
  free(m_vtFace);
  free(m_velocityMisfit);
  free(m_massBalanceMisfit);
  free(m_TikhonovPenalty);
  free(m_barrierPenalty);
}

//apply preconditioner s = M^{-1}r
void AMRIceControl::preCond
(Vector<LevelData<FArrayBox>* >& a_s, 
 const Vector<LevelData<FArrayBox>* >& a_r)
{
  assign(a_s,a_r);
}

// duplicate storage of a_b in a_a
void AMRIceControl::create
(Vector<LevelData<FArrayBox>* >& a_a, 
 const  Vector<LevelData<FArrayBox>* >& a_b)
{
  a_a.resize(a_b.size());
  for (int lev = 0; lev < a_b.size(); lev++)
    {
      if (a_a[lev] == NULL)
	a_a[lev] = new LevelData<FArrayBox>();
      m_vectOps[lev].create(*a_a[lev],*a_b[lev]);
    }
} 



  // set a_x = s * a_x
void AMRIceControl::setToZero
(Vector<LevelData<FArrayBox>* >& a_x)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].setToZero(*a_x[lev]);
}

  // set a_x = s * a_x
void AMRIceControl::scale
(Vector<LevelData<FArrayBox>* >& a_x, 
 const  Real a_s)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].scale(*a_x[lev],a_s);
}

  // set a_y = a_x
void AMRIceControl::assign
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x)
{
  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].assign(*a_y[lev],*a_x[lev]);
}
  
// set a_y = a_y + a_s * a_x
void AMRIceControl::incr
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x, Real a_s)
{
for (int lev = 0; lev < a_x.size(); lev++)
  m_vectOps[lev].incr(*a_y[lev],*a_x[lev],a_s);
}
  
  // return a_y.a_x
Real AMRIceControl::dotProduct
(Vector<LevelData<FArrayBox>* >& a_1, 
 const Vector<LevelData<FArrayBox>* >& a_2)
{
  
  //AMR calculation copied and pasted from MultiLevelLinearOp

  CH_TIME("AMRIceControl::dotProduct");

  // want to do this in an AMR way, so need to generate a temporary
  Vector<LevelData<FArrayBox>* > temp1, temp2;
  create(temp1, a_1);
  create (temp2, a_2);

  // first set to zero, then call assign, since assign only sets
  // valid regions (that way ghost cells are set to zero)

  setToZero(temp1);
  setToZero(temp2);

  assign(temp1, a_1);
  assign(temp2, a_2);

  // now set covered regions to zero
  for (int level =0 ; level<temp1.size()-1; level++)
    {
      LevelData<FArrayBox>& temp1Level = *temp1[level];
      LevelData<FArrayBox>& temp2Level = *temp2[level];

      CH_assert(temp1[level]->getBoxes() == temp2[level]->getBoxes());
      CH_assert(temp1[level+1] != NULL);

      int nRefFine = m_refRatio[level];
      const DisjointBoxLayout& finerGrids = temp1[level+1]->getBoxes();
      const DisjointBoxLayout& levelGrids = temp1[level]->getBoxes();
      DataIterator levelDit = levelGrids.dataIterator();
      LayoutIterator finerLit = finerGrids.layoutIterator();

      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          const Box& thisBox = levelGrids[levelDit];
          for (finerLit.begin(); finerLit.ok(); ++finerLit)
            {
              Box testBox = finerGrids[finerLit];
              testBox.coarsen(nRefFine);
              testBox &= thisBox;
              if (!testBox.isEmpty())
                {
                  temp1Level[levelDit].setVal(0.0, testBox,
                                              0, temp1Level.nComp());

                  temp2Level[levelDit].setVal(0.0, testBox,
                                              0, temp2Level.nComp());
                }
            } // end loop over finer boxes
        } // end loop over boxes on this level
    } // end loop over levels for setting covered regions to zero;

  // now loop over levels and call AMRLevelOp dotProduct
  Real prod = 0.0;
  for (int level = 0; level<temp1.size(); level++)
    {
      Real levelProd = m_vectOps[level].dotProduct(*temp1[level],
						   *temp2[level]);

      // incorporate scaling for AMR dot products
      RealVect dxLev = m_dx[level];
      Real scale = D_TERM(dxLev[0], *dxLev[1], *dxLev[2]);
      levelProd *= scale;
      prod += levelProd;
    }

  free(temp1);
  free(temp2);

  return prod;

}

// how many degrees of freedom are there in x?
int AMRIceControl::nDoF(const Vector<LevelData<FArrayBox>* >& x)
{
  return 10;
}
void AMRIceControl::checkBounds(Vector<LevelData<FArrayBox>* >& a_X)
{
  for (int lev =0; lev < a_X.size(); lev++)
    {
      for (DataIterator dit = a_X[lev]->dataIterator();dit.ok();++dit)
	{
	  FArrayBox& x = (*a_X[lev])[dit];
	  FORT_BOUNDCTRL(CHF_FRA1(x,0),
	  		 CHF_CONST_REAL(m_lowerX0),
			 CHF_CONST_REAL(m_upperX0),
	  		 CHF_BOX(x.box()));
	  FORT_BOUNDCTRL(CHF_FRA1(x,1),
	  		 CHF_CONST_REAL(m_lowerX1),
			 CHF_CONST_REAL(m_upperX1),
	  		 CHF_BOX(x.box()));

	}
    }

}

void AMRIceControl::computeObjectiveAndGradient
(Real& a_fm, Real& a_fp,
 Vector<LevelData<FArrayBox>* >& a_g, 
 const  Vector<LevelData<FArrayBox>* >& a_x, 
 bool a_inner)
{

  //compute C and muCoef from a_x : 
  //C = exp(a_x[0]) * COrigin
  // muCoef = exp(a_x[1]) * muCoefOrigin

#define CCOMP 0
#define MUCOMP 1

  setToZero(a_g);


  // probably excessive, but ...
  for (int lev=m_finestLevel; lev > 0 ;lev--)
    {
      CoarseAverage avg(m_grids[lev],a_x[lev]->nComp(),m_refRatio[lev-1]);
      avg.averageToCoarse(*a_x[lev-1],*a_x[lev]);
    }

  for (int lev=0; lev <= m_finestLevel;lev++)
    {

      const LevelData<FArrayBox>& levelX =  *a_x[lev];
      const LevelSigmaCS& levelCS =  *m_coordSys[lev];
      LevelData<FArrayBox>& levelC =  *m_C[lev];
      LevelData<FArrayBox>& levelCcopy =  *m_Ccopy[lev];
      LevelData<FArrayBox>& levelCOrigin =  *m_COrigin[lev];
      LevelData<FArrayBox>& levelMuCoef =  *m_muCoef[lev];
      LevelData<FArrayBox>& levelMuCoefOrigin =  *m_muCoefOrigin[lev];
      const DisjointBoxLayout levelGrids =  m_grids[lev];
      
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  const Box& box = levelGrids[dit];

	  
	  FArrayBox& thisC = levelC[dit];
	  thisC.copy(levelCOrigin[dit]);
	  FORT_BOUNDEXPCTRL(CHF_FRA1(thisC,0),
	  		    CHF_CONST_FRA1(levelX[dit],CCOMP),
	  		    CHF_CONST_REAL(m_lowerX0),
			    CHF_CONST_REAL(m_upperX0),
	  		    CHF_BOX(box));


	  FArrayBox& thisCcopy = levelCcopy[dit];
	  thisCcopy.copy(thisC);

	  FArrayBox& thisMuCoef = levelMuCoef[dit];
	  //Real initialMuCoef = 1.0;
	  //thisMuCoef.setVal(initialMuCoef);
	  thisMuCoef.copy(levelMuCoefOrigin[dit]);
	  FORT_BOUNDEXPCTRL(CHF_FRA1(thisMuCoef,0),
	  		    CHF_CONST_FRA1(levelX[dit],MUCOMP),
	  		    CHF_CONST_REAL(m_lowerX1),
			    CHF_CONST_REAL(m_upperX1),
	  		    CHF_BOX(box));
	  	  

	  if (m_muCoefLEQOne || m_muCoefLEQOneGround)
	    {

	      
	      const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	      for (BoxIterator bit(thisMuCoef.box());bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if (m_muCoefLEQOne || (m_muCoefLEQOneGround && mask(iv) == GROUNDEDMASKVAL))
		    {
		      thisMuCoef(iv) = std::min(1.0, thisMuCoef(iv));
		    }
		}
	      
	    }
	  if (m_muCoefGEQOne || m_muCoefGEQOneGround)
	    {

	      const LevelSigmaCS& levelCS =  *m_coordSys[lev];
	      const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	      for (BoxIterator bit(thisMuCoef.box());bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if (m_muCoefGEQOne || (m_muCoefGEQOneGround && mask(iv) == GROUNDEDMASKVAL))
		    {
		      thisMuCoef(iv) = std::max(1.0, thisMuCoef(iv));
		    }
		}
	      
	    }
	  
	}
    
      //boundary values for C and muCoef: extrapolation or reflection should be sufficient
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  if (! m_domain[lev].isPeriodic(dir))
	    {
	      ReflectGhostCells(levelC, m_domain[lev], dir, Side::Lo);
	      ReflectGhostCells(levelC, m_domain[lev], dir, Side::Hi);
	      ReflectGhostCells(levelMuCoef, m_domain[lev], dir, Side::Lo);
	      ReflectGhostCells(levelMuCoef, m_domain[lev], dir, Side::Hi);
	    }
	}

      if (lev > 0)
	{

	  CoarseAverage avg(m_grids[lev],1,m_refRatio[lev-1]);
	  avg.averageToCoarse(*m_C[lev-1],*m_C[lev]);
	  avg.averageToCoarse(*m_muCoef[lev-1],*m_muCoef[lev]);
		  
	  //  QuadCFInterp qcfi(levelGrids, &m_grids[lev-1], m_dx[lev][0], 
	  //		    m_refRatio[lev-1],1, m_domain[lev]);
	  //
	  //qcfi.coarseFineInterp(levelC, *m_C[lev-1]);
	  //qcfi.coarseFineInterp(levelMuCoef, *m_muCoef[lev-1]);
	  int nGhost = 1;
	  PiecewiseLinearFillPatch li(levelGrids, 
				      m_grids[lev-1],
				      1, 
				      m_domain[lev-1],
				      m_refRatio[lev-1],
				      nGhost);
              
	  // since we're not subcycling, don't need to interpolate in time
	  Real time_interp_coeff = 0.0;
	  li.fillInterp(levelMuCoef,
			*m_muCoef[lev-1],
			*m_muCoef[lev-1],
			time_interp_coeff,
			0, 0, 1);
	  li.fillInterp(levelC,
			*m_C[lev-1],
			*m_C[lev-1],
			time_interp_coeff,
			0, 0, 1);

             

        }

    
      levelC.exchange();
      levelMuCoef.exchange();
    
      CellToEdge(levelMuCoef, *m_faceMuCoef[lev]);
      IceVelocity::applyHelmOp(*m_lapC[lev], levelC, 0.0, 1.0, 
			       m_grids[lev], m_dx[lev]);
      IceVelocity::applyGradSq(*m_gradCSq[lev], levelC,m_grids[lev], m_dx[lev]);
      IceVelocity::applyHelmOp(*m_lapMuCoef[lev], levelMuCoef, 0.0, 1.0,
			       m_grids[lev], m_dx[lev]);
      IceVelocity::applyGradSq(*m_gradMuCoefSq[lev], levelMuCoef,  m_grids[lev], m_dx[lev]);
    }

  {
    Real maxC = computeMax(m_C,m_refRatio);
    Real minC = computeMin(m_C,m_refRatio);
    Real maxMuCoef = computeMax(m_muCoef,m_refRatio);
    Real minMuCoef = computeMin(m_muCoef,m_refRatio);
    Real maxX0 = computeMax(a_x,m_refRatio,Interval(0,0));
    Real minX0 = computeMin(a_x,m_refRatio,Interval(0,0));
    Real maxX1 = computeMax(a_x,m_refRatio,Interval(1,1));
    Real minX1 = computeMin(a_x,m_refRatio,Interval(1,1));
    
    pout() << " max X[0] = " << maxX0
	   << " min X[0] = " << minX0
	   << " max X[1] = " << maxX1
	   << " min X[1] = " << minX1 
	   << std::endl
	   << " max C = " << maxC 
	   << " min C = " << minC
	   << " max muCoef = " << maxMuCoef
	 << " min muCoef = " << minMuCoef
	   << std::endl;
  }
  //solve the forward problem (to update mu)



  IceVelocity::setFloatingC(m_C, m_coordSys, m_grids, 0.0);

  // add drag due to ice in contact with ice-free rocky walls
	 
  ParmParse ppamr("amr");
  bool wallDrag = true; 
  ppamr.query("wallDrag", wallDrag);
  if (wallDrag)
    {
      Real wallDragExtra = 0.0;
      ppamr.query("wallDragExtra", wallDragExtra);

      for (int lev=0; lev <= m_finestLevel;lev++)
	{
	  const LevelSigmaCS& levelCS =  *m_coordSys[lev];
	  LevelData<FArrayBox>& levelC =  *m_C[lev];
	  const LevelData<FArrayBox>& levelCCopy =  *m_Ccopy[lev]; //not set to zero in shelves
	  const DisjointBoxLayout levelGrids =  m_grids[lev];
	  for (DataIterator dit(levelGrids);dit.ok();++dit)
	    {

	      FArrayBox wallC(levelGrids[dit],1);
	      wallC.setVal(0.0);
	      IceVelocity::addWallDrag(wallC, 
				       levelCS.getFloatingMask()[dit], levelCS.getSurfaceHeight()[dit],
				       levelCS.getH()[dit], levelCS.getTopography()[dit], 
				       levelCCopy[dit], wallDragExtra,m_dx[lev],levelGrids[dit]);
	      
	      levelC[dit] += wallC;
	    }
	}
    }


  

  if (!m_velocityInitialised)
    {
      ParmParse pp("control");
      Real initialMu = 1.0e+6;
      pp.query("initialMu",initialMu);
      constMuRelation cmr; cmr.setConstVal(initialMu);
      ConstitutiveRelation *tmpPtr = m_constRelPtr;
      m_constRelPtr = &cmr;
      solveForwardProblem(m_velb,true,m_rhs,m_C,m_C0,m_A,m_faceMuCoef);
      m_constRelPtr = tmpPtr;
      m_velocityInitialised = true;
    }
  solveForwardProblem(m_velb,false,m_rhs,m_C,m_C0,m_A,m_faceMuCoef);
#if BISCICLES_Z == BISICLES_LAYERED
  //compute the surface velocity

  LevelData<FArrayBox>* cellDiffusivity = NULL;

  for (int lev = 0; lev <= m_finestLevel ;lev++)
    {
      int nLayer = m_A[lev]->nComp();
      LevelData<FArrayBox> layerSFaceXYVel(m_grids[lev],SpaceDim*(nLayer+1), IntVect::Unit);
      LevelData<FluxBox> layerXYFaceXYVel(m_grids[lev],nLayer, IntVect::Unit);
      LevelData<FluxBox> faceVelAdvection(m_grids[lev],1, IntVect::Unit);
      LevelData<FluxBox> faceVelTotal(m_grids[lev],1, IntVect::Unit);
      LevelData<FluxBox> faceDiffusivity(m_grids[lev],1, IntVect::Unit);
       
      LevelData<FArrayBox>* crseCellDiffusivityPtr = 
	(lev > 0)?cellDiffusivity:NULL;
      
      cellDiffusivity = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
      CH_assert(cellDiffusivity != NULL);
      LevelData<FArrayBox> sA(m_grids[lev],1, IntVect::Unit);
      //m_A[lev]->copyTo(Interval(0,0),sA,Interval(0,0));
      LevelData<FArrayBox> bA(m_grids[lev],1, IntVect::Unit);
      //m_A[lev]->copyTo(Interval(nLayer-1,nLayer-1),bA,Interval(0,0));

      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  const FArrayBox& A = (*m_A[lev])[dit];
	  sA[dit].copy(A,0,0);
	  bA[dit].copy(A,nLayer-1,0);
	}

      IceVelocity::computeFaceVelocity
	(faceVelAdvection, faceVelTotal, faceDiffusivity, *cellDiffusivity, 
	 layerXYFaceXYVel,  layerSFaceXYVel,
	 *m_velb[lev], *m_coordSys[lev],*m_A[lev], sA, bA,
	 (lev > 0)?m_velb[lev-1]:NULL,
	 crseCellDiffusivityPtr,
	 (lev > 0)?m_refRatio[lev-1]:1,
	 m_constRelPtr, true);
      
      if (crseCellDiffusivityPtr != NULL)
	delete crseCellDiffusivityPtr;

      layerSFaceXYVel.copyTo(Interval(0,1),*m_vels[lev],Interval(0,1));
    }
  if (cellDiffusivity != NULL)
    delete cellDiffusivity;

#endif

  //this might not be necessary, but it makes the plotted fields look better
  for (int lev=m_finestLevel; lev > 0 ;lev--)
    {
      CoarseAverage avg(m_grids[lev],SpaceDim,m_refRatio[lev-1]);
      avg.averageToCoarse(*m_velb[lev-1],*m_velb[lev]);
      avg.averageToCoarse(*m_vels[lev-1],*m_vels[lev]);
    }

  



  //compute cell-centered div(UH), needed to compute part of the adjoint equation's rhs
  for (int lev=0; lev <= m_finestLevel ;lev++)
    {
      IceVelocity::computeFaceFlux
	(*m_faceThckFlux[lev],*m_velb[lev],m_coordSys[lev]->getH(),m_grids[lev]);
    }
  for (int lev=m_finestLevel; lev > 0 ;lev--)
    {
      CoarseAverageFace avg(m_grids[lev],1,m_refRatio[lev-1]);
      avg.averageToCoarse(*m_faceThckFlux[lev-1],*m_faceThckFlux[lev]);
    }
  for (int lev=0; lev <= m_finestLevel ;lev++)
    {
      IceVelocity::applyDiv(*m_divUH[lev],*m_faceThckFlux[lev],m_grids[lev],m_dx[lev]);
      if (lev > 0)
	{
	  // fill ghost cells 
	  QuadCFInterp qcfi(m_grids[lev], &m_grids[lev-1], m_dx[lev][0], 
			    m_refRatio[lev-1],1, m_domain[lev]);
	  qcfi.coarseFineInterp(*m_divUH[lev], *m_divUH[lev-1]);
	  qcfi.coarseFineInterp(*m_divUHObs[lev], *m_divUHObs[lev-1]);
	}
    }
  //compute RHS for adjoint problem, and set m_adjVel = m_vel 
  //so that the viscosities are computed correctly. 
  //\todo see about pulling the viscosities from the forward problem

  for (int lev=0; lev <= m_finestLevel ;lev++)
    {
      const LevelData<FArrayBox>& levelVelObs =  *m_velObs[lev];
      LevelData<FArrayBox>& levelVelCoef =  *m_velCoef[lev];
      const LevelData<FArrayBox>& levelVel =  *m_velb[lev];
      LevelData<FArrayBox>& levelAdjRhs = *m_adjRhs[lev];
      LevelData<FArrayBox>& levelAdjVel = *m_adjVel[lev];
      
      

      //contribution due to speed misfit
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{

	  levelAdjVel[dit].copy(levelVel[dit]);

	  FArrayBox& adjRhs = levelAdjRhs[dit];
	  const FArrayBox& um = levelVel[dit];
	  const FArrayBox& uo = levelVelObs[dit];
	  const Box& box = (m_grids[lev])[dit];
	  
	  FArrayBox& misfit = (*m_velocityMisfit[lev])[dit];

	  FORT_ADJRHSSPEEDCTRL(CHF_FRA1(adjRhs,0), CHF_FRA1(adjRhs,1),
			       CHF_CONST_FRA1(misfit,0),
			       CHF_CONST_FRA1(um,0), CHF_CONST_FRA1(um,1),
			       CHF_CONST_FRA1(uo,0), CHF_CONST_FRA1(uo,1),
			       CHF_BOX(box));    
			     
	  

	  for (BoxIterator bit(levelVelCoef[dit].box());bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (levelVelCoef[dit](iv) < 0.975)
		levelVelCoef[dit](iv) = 0.0;
	    }

	  misfit.mult(levelVelCoef[dit]);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      levelAdjRhs[dit].mult(levelVelCoef[dit],0,dir);
	    }

	  levelAdjRhs[dit] *= m_velMisfitCoefficient;
	  misfit *= m_velMisfitCoefficient;
	  
	}

      const LevelData<FArrayBox>& levelDivUHObs =  *m_divUHObs[lev];
      LevelData<FArrayBox>& levelDivUHCoef =  *m_divUHCoef[lev];
      const LevelData<FArrayBox>& levelDivUH =  *m_divUH[lev];
      const LevelSigmaCS& levelCS =  *m_coordSys[lev];


      //contribution due to mass imbalance
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  
	  FArrayBox& adjRhs = levelAdjRhs[dit];
	  const FArrayBox& divuhm = levelDivUH[dit];
	  const FArrayBox& divuho = levelDivUHObs[dit];
	  const Box& box = m_grids[lev][dit];
	  //Box gbox(box);
	  //gbox.grow(1);
	  FArrayBox& misfit = (*m_massBalanceMisfit[lev])[dit];
	  misfit.copy(divuhm);misfit.minus(divuho);
	  Real dx = m_dx[0][0];
	  
	  FArrayBox tmp(box,SpaceDim);
	  Box sbox(box);
	  Box sdom(m_domain[lev].domainBox());
	  sdom.grow(-1);
	  sbox &= sdom;
	  tmp.setVal(0.0);
	  FORT_ADJRHSMASSCTRL(CHF_FRA(tmp),
			      CHF_CONST_FRA1(misfit,0),
			      CHF_CONST_FRA1(levelCS.getH()[dit],0),
			      CHF_CONST_REAL(dx),
			      CHF_BOX(sbox)); 

	  for (BoxIterator bit(levelVelCoef[dit].box());bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (levelDivUHCoef[dit](iv) < 0.975)
		levelDivUHCoef[dit](iv) = 0.0;
	    }

	  misfit.mult(levelDivUHCoef[dit]);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      tmp.mult(levelDivUHCoef[dit],0,dir);
	    }
	  
	  tmp *= m_massImbalanceCoefficient;
	  misfit *= m_massImbalanceCoefficient;
	  adjRhs += tmp;
	}

    }
  
  //compute the objective function 
  Real vobj = computeSum(m_velocityMisfit, m_refRatio, m_dx[0][0]);
  Real mobj = computeSum(m_massBalanceMisfit, m_refRatio, m_dx[0][0]);
  Real sumGradCSq = computeSum(m_gradCSq,m_refRatio, m_dx[0][0]);
  Real sumGradMuSq = computeSum(m_gradMuCoefSq,m_refRatio, m_dx[0][0]);
  Real normX0 = computeNorm(a_x,m_refRatio, m_dx[0][0], Interval(0,0));
  Real normX1 = computeNorm(a_x,m_refRatio, m_dx[0][0], Interval(1,1));

  a_fm = vobj + mobj ;
  a_fp =  m_gradCsqRegularization * sumGradCSq
    + m_gradMuCoefsqRegularization * sumGradMuSq
    + m_X0Regularization * normX0*normX0
    + m_X1Regularization * normX1*normX1;

  pout() << " ||velocity misfit||^2 = " << vobj
	 << " ||mass balance misfit||^2 = " << mobj
	 << std::endl;
  //solve the adjoint problem
  solveForwardProblem(m_adjVel,true,m_adjRhs,m_C,m_C0,m_A,m_faceMuCoef);
  //this might not be necessary, but it makes the plotted fields look better
  for (int lev=m_finestLevel; lev > 0 ;lev--)
    {
      CoarseAverage avg(m_grids[lev],SpaceDim,m_refRatio[lev-1]);
      avg.averageToCoarse(*m_adjVel[lev-1],*m_adjVel[lev]);
    }

  //initialize  the gradient vector with the regularization terms R = (- a C lap(C)) 
 

  
  for (int lev = 0; lev <= m_finestLevel; lev++)
    {
      LevelData<FArrayBox>& levelG =  *a_g[lev];
      const LevelData<FArrayBox>& levelX =  *a_x[lev];
      const LevelData<FArrayBox>& levelC =  *m_C[lev];
      const LevelData<FArrayBox>& levelMuCoef  =  *m_muCoef[lev];
      const LevelData<FArrayBox>& levelLapC  =  *m_lapC[lev];
      const LevelData<FArrayBox>& levelLapMuCoef  =  *m_lapMuCoef[lev];

      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  FArrayBox& thisG =  levelG[dit];
	  const FArrayBox& thisX =  levelX[dit];
	  thisG.setVal(0.0,CCOMP);
	  thisG.setVal(0.0,MUCOMP);

	  FArrayBox t(thisG.box(),1);
	  t.copy(levelLapC[dit]);t*= levelC[dit]; t*= -m_gradCsqRegularization;
	  thisG.plus(t,0,CCOMP);

	
	  t.copy(levelLapMuCoef[dit]);t*= levelMuCoef[dit]; t*= -m_gradMuCoefsqRegularization;
	  thisG.plus(t,0,MUCOMP);

	  
	  t.copy(thisX,CCOMP,0);
	  t *= m_X0Regularization;
	  thisG.plus(t,0,CCOMP);

	  t.copy(thisX,MUCOMP,0);
	  t *= m_X1Regularization;
	  thisG.plus(t,0,MUCOMP);
	}

    }
  
  //compute unregularized part of gradient with respect to a_x component 0 (basal friction coefficient C)
  // = - adjVel * vel * C 
  for (int lev = 0; lev <= m_finestLevel; lev++)
    {
      const LevelData<FArrayBox>& levelLambda =  *m_adjVel[lev];
      const LevelData<FArrayBox>& levelU =  *m_velb[lev];
      const LevelData<FArrayBox>& levelC =  *m_C[lev];
      LevelData<FArrayBox>& levelG =  *a_g[lev];
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  const FArrayBox& thisU = levelU[dit];
	  const FArrayBox& thisC = levelC[dit];
	  const FArrayBox& thisLambda = levelLambda[dit];
	  
	 
	  FArrayBox& thisG = levelG[dit];
	  //FArrayBox t(thisG.box(),1);
	  //t.copy(thisU,0,0);
	  //t.mult(thisLambda,0,0);
	  //t.mult(thisC,0,0);
	  //thisG.minus(t,0,CCOMP);
	  //t.copy(thisU,1,0);
	  //t.mult(thisLambda,1,0);
	  //t.mult(thisC,0,0);
	  //thisG.minus(t,0,CCOMP);
	  
	  FArrayBox t(thisG.box(),1);
	  FORT_CADOTBCTRL(CHF_FRA1(t,0),
			  CHF_CONST_FRA1(thisC,0),
			  CHF_CONST_FRA(thisU),
			  CHF_CONST_FRA(thisLambda),
			  CHF_BOX(t.box()));

	  t *= -1.0;
	  if (m_boundMethod == Projection)
	    {
	      const LevelData<FArrayBox>& levelX =  *a_x[lev];
	      const FArrayBox& thisX =  levelX[dit];
	      //FORT_MULTHATCTRL
	      FORT_HARDPOINTINCTRL(CHF_FRA1(t,0),
			       CHF_CONST_FRA1(thisX,CCOMP),
			       CHF_CONST_REAL(m_lowerX0),
			       CHF_CONST_REAL(m_upperX0),
			       CHF_BOX(t.box()));
	    }

	  thisG.plus(t,0,CCOMP);
	  
	  
	}
    }
  {
    //compute unregularized part of gradient with respect to a_x component 1 
    //(viscosity multiplier muCoef)
    IceJFNKstate state(m_grids, m_refRatio, m_domain, m_dx,
		       m_coordSys, m_velb, m_C, m_C0, m_finestLevel ,
		       *m_constRelPtr, *m_bfRelPtr,
		       *m_ibcPtr, m_A, m_faceA, 0.0, 0.0, 0, 0.0);
    state.setState(m_velb);
    state.computeViscousTensorFace(m_vtFace);

    for (int lev = 0; lev <= m_finestLevel; lev++)
      {
	LevelData<FArrayBox>& levelG =  *a_g[lev];
	const LevelData<FluxBox>& levelVT =  *m_vtFace[lev];
	const LevelData<FArrayBox>& levelMuCoef  =  *m_muCoef[lev];
	const LevelData<FArrayBox>& levelLambda =  *m_adjVel[lev];
	for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	  {
	    const FArrayBox& muCoef = levelMuCoef[dit];
	    const FArrayBox& lambda = levelLambda[dit];
	    const Box& box = m_grids[lev][dit]; 
	    const FluxBox& vt = levelVT[dit];
	    FArrayBox t(box,1);
	    t.setVal(0.0);
	    for (int facedir = 0; facedir < SpaceDim; facedir++)
	      {
		const IntVect e = BASISV(facedir);
		for (BoxIterator bit(box);bit.ok();++bit)
		  {
		    for (int dir = 0; dir < SpaceDim; dir++)
		      {
			const IntVect& iv = bit();
			
			t(iv) += (lambda(iv+e,dir)-lambda(iv,dir))* vt[facedir](iv+e,dir)
			  + (lambda(iv,dir)-lambda(iv-e,dir))* vt[facedir](iv,dir);
		      }
		  }
	      }

	    RealVect oneOnTwoDx = 1.0/ (2.0 * m_dx[lev]);
	    t.mult(-oneOnTwoDx[0]);
	    t.mult(muCoef);
	    
	    if (m_boundMethod == Projection)
	      {
		const LevelData<FArrayBox>& levelX =  *a_x[lev];
		const FArrayBox& thisX =  levelX[dit];
		//FORT_MULTHATCTRL
		
		FORT_HARDPOINTINCTRL(CHF_FRA1(t,0),
				     CHF_CONST_FRA1(thisX,MUCOMP),
				     CHF_CONST_REAL(m_lowerX1),
				     CHF_CONST_REAL(m_upperX1),
				     CHF_BOX(t.box()));
	    }


	    FArrayBox& thisG =  levelG[dit];
	    thisG.plus(t,0,MUCOMP);

	  }
      }

  }

  //dump data
  if (a_inner)
    {
      if (m_writeInnerSteps)
	writeState(innerStateFile(), m_innerCounter, a_x, a_g);
      m_innerCounter++;
    }
  else
    {
      writeState(outerStateFile(), m_innerCounter, a_x, a_g);
      m_innerCounter++;
      m_outerCounter++;
    }
}
 
std::string AMRIceControl::outerStateFile() const
{
  std::stringstream ss;
  ss << m_outerStepFileNameBase;
  ss.width(6);ss.fill('0');ss << m_outerCounter;
  ss.width(0);ss << ".2d.hdf5";
  return ss.str(); 
}    

std::string AMRIceControl::innerStateFile() const
{
  std::stringstream ss;
  ss << m_innerStepFileNameBase;
  ss.width(6);ss.fill('0');ss << m_innerCounter;
  ss.width(0);ss << ".2d.hdf5";
  return ss.str(); 
}


void AMRIceControl::readState(const std::string& a_file, int& a_counter,
			      Vector<LevelData<FArrayBox>* >& a_x) 
{
  pout() << "reading state from " << a_file << std::endl;
  Vector<std::string> restartNames;
  Vector<LevelData<FArrayBox>* > restartData;
  Vector<DisjointBoxLayout> restartGrids;
  Vector<int> restartRatio;
  Real restartDx = 0.0, restartDt = 0.0, restartTime = 0.0;
  Box restartDomBox;
  int restartNumLevels;
  int status = ReadAMRHierarchyHDF5
    (a_file,restartGrids,restartData,restartNames,
     restartDomBox,restartDx,restartDt,restartTime,
     restartRatio,restartNumLevels);
  CH_assert(restartDx == m_dx[0][0]);
  CH_assert(status == 0);
  CH_assert(restartNumLevels == m_finestLevel + 1);
    
  for (int lev = restartNumLevels-1;lev > 0 ; lev--)
    {
      LevelData<FArrayBox>& fine = *restartData[lev];
      LevelData<FArrayBox>& crse = *restartData[lev-1];
      CoarseAverage ca(restartGrids[lev],fine.nComp(),restartRatio[lev-1]);
      ca.averageToCoarse(crse,fine);
    }

  a_counter = int(restartTime) + 1;
  for (int lev = 0; lev < restartNumLevels; lev++)
    {
      int j = 0;
      restartData[lev]->copyTo(Interval(j,j+1),*(a_x[lev]),Interval(0,1));
      j+=2;
      //restartData[lev]->copyTo(Interval(j,j),*(m_Ccopy[lev]),Interval(0,0));
      j++;
      //restartData[lev]->copyTo(Interval(j,j),*(m_C[lev]),Interval(0,0));
      j++;
      //restartData[lev]->copyTo(Interval(j,j),*(m_muCoef[lev]),Interval(0,0));
      j++;
      restartData[lev]->copyTo(Interval(j,j+1),*(m_velb[lev]),Interval(0,1));j+=2;
      
      //restartData[lev]->copyTo(Interval(j,j+1),*(m_vels[lev]),Interval(0,1));
      j+=2;

      //restartData[lev]->copyTo(Interval(j,j+1),*(m_velObs[lev]),Interval(0,1));
      j+=2;
      //restartData[lev]->copyTo(Interval(j,j),*(m_divUH[lev]),Interval(0,0));
      j++;
      //restartData[lev]->copyTo(Interval(j,j),*(m_divUHObs[lev]),Interval(0,0));j++;
      restartData[lev]->copyTo(Interval(j,j+1),*(m_adjVel[lev]),Interval(0,1));
      j+=2;
      restartData[lev]->copyTo(Interval(j,j+1),*(m_adjRhs[lev]),Interval(0,1));
      j+=2;
      //restartData->copyTo(Interval(j,j+1),a_g[lev],Interval(0,1));
      j+=2;
      //restartData[lev]->copyTo(Interval(j,j),*(m_velCoef[lev]),Interval(0,0));
      j++;
      //restartData[lev]->copyTo(Interval(j,j),*(m_divUHCoef[lev]),Interval(0,0));
      j++;
    }
    
  

}


void AMRIceControl::writeState(const std::string& a_file, int a_counter,
			       const Vector<LevelData<FArrayBox>* >& a_x,
			       const Vector<LevelData<FArrayBox>* >& a_g) const
{
  pout() << "writing state to " << a_file << std::endl;
  Vector<std::string> names;
  names.resize(0);
  names.push_back("X0");
  names.push_back("X1");
  names.push_back("C");
  names.push_back("Cwshelf");
  names.push_back("muCoef");
  names.push_back("xVelb");
  names.push_back("yVelb");
  names.push_back("xVels");
  names.push_back("yVels");
  names.push_back("xVelo");
  names.push_back("yVelo");
  names.push_back("divuh");
  names.push_back("divuho");
  names.push_back("xAdjVel");
  names.push_back("yAdjVel");
  names.push_back("xAdjRhs");
  names.push_back("yAdjRhs");
  names.push_back("gradJC");
  names.push_back("gradJMuCoef");
  names.push_back("velc");
  names.push_back("divuhc");
  names.push_back("topg");
  names.push_back("thck");
  names.push_back("usrf");

  Vector<LevelData<FArrayBox>*> vdata(m_finestLevel+1);
  for (int lev = 0; lev <= m_finestLevel;lev++)
    {
      vdata[lev] = new LevelData<FArrayBox>(m_grids[lev],names.size(),IntVect::Zero);
      LevelData<FArrayBox>& data = *vdata[lev];
      int j = 0;
      a_x[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      a_x[lev]->copyTo(Interval(1,1),data,Interval(j,j));j++;
      //use m_Ccopy rather than m_C because m_C is set to zero in shelves etc
      m_Ccopy[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_C[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;  
      m_muCoef[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++; 
      m_velb[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_vels[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_velObs[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_divUH[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_divUHObs[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_adjVel[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_adjRhs[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      a_g[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      a_g[lev]->copyTo(Interval(1,1),data,Interval(j,j));j++;
      m_velCoef[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_divUHCoef[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_coordSys[lev]->getTopography().copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_coordSys[lev]->getH().copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_coordSys[lev]->getSurfaceHeight().copyTo(Interval(0,0),data,Interval(j,j));j++;
    }
  const Real dt = 1.0;
  const Real time = Real(a_counter);
  
  for (int lev = vdata.size()-1;lev > 0 ; lev--)
    {
      LevelData<FArrayBox>& fine = *vdata[lev];
      LevelData<FArrayBox>& crse = *vdata[lev-1];
      CoarseAverage ca(m_grids[lev],fine.nComp(),m_refRatio[lev-1]);
      ca.averageToCoarse(crse,fine);
    }

  WriteAMRHierarchyHDF5(a_file,m_grids,vdata,names, m_domain[0].domainBox(),
			m_dx[0][0], dt, time , m_refRatio, vdata.size());
  
  
  
  for (int lev = 0; lev <= m_finestLevel;lev++)
    {
      delete vdata[lev];
    }
  
}





void AMRIceControl::solveForwardProblem(Vector<LevelData<FArrayBox>* >& a_u,
					const bool a_linear,
					const Vector<LevelData<FArrayBox>* >& a_rhs,
					const Vector<LevelData<FArrayBox>* >& a_C,
					const Vector<LevelData<FArrayBox>* >& a_C0,
					const Vector<LevelData<FArrayBox>* >& a_A,
					const Vector<LevelData<FluxBox>* >& a_muCoef)
{

  JFNKSolver jfnkSolver;
  jfnkSolver.define(m_domain[0], m_constRelPtr , m_bfRelPtr,
		    m_grids, m_refRatio, m_dx[0], m_ibcPtr, m_finestLevel+1);
    
  Real initialNorm = 1.0; Real finalNorm = 1.0; Real convergenceMetric=-1.0;

  if (a_linear)
    {
      //jfnkSolver.m_RelaxRelTol = 1.0e-10;
      jfnkSolver.m_maxRelaxIter = jfnkSolver.m_maxRelaxIter * 2;
      jfnkSolver.m_RelaxHang = 0.0; // don't allow hangs
      //jfnkSolver.m_writeResiduals = true;
    }

  jfnkSolver.solve(a_u, initialNorm, finalNorm, convergenceMetric, a_linear, 
		   a_rhs, a_C, a_C0, a_A, a_muCoef, m_coordSys, 0.0 , 0, m_finestLevel);


}

void AMRIceControl::CviscoustoCplastic(LevelData<FArrayBox>& a_Cplastic, 
				       const LevelData<FArrayBox>& a_Cviscous, 
				       const LevelData<FArrayBox>&  a_u,
				       const LevelData<BaseFab<int> >&  a_mask) const
{

  CH_assert(SpaceDim < 3); //not done for spacedim = 3...

  const DisjointBoxLayout& grids = a_Cplastic.disjointBoxLayout();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      FArrayBox& Cp =  a_Cplastic[dit];
      const FArrayBox& Cv = a_Cviscous[dit];
      const FArrayBox& u = a_u[dit];
      const BaseFab<int>& mask = a_mask[dit];
      //Cp.setVal(0.0);
      Cp.copy(Cv);
      for (BoxIterator bit(Cp.box());bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) == GROUNDEDMASKVAL)
	    {
	      Real usq = (D_TERM(u(iv,0)*u(iv,0), + u(iv,1)*u(iv,1), +0.0));
	      usq = std::min(usq,1e+8);
	      Cp(iv) = Cv(iv) * std::pow(usq,1.0/3.0);
	    }
	}

    }

}
  
#include "NamespaceFooter.H"
