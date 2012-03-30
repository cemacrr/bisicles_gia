#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// SolveControlProblem.cpp
// contruct and solve a control problem
// by a method related to Vieli and Payne, Ann Glaciol, 2003. vol 26
// assumes a 2D SSA model or L1L2 model
//===========================================================================

#include "SolveControlProblem.H"
#include "JFNKSolver.H"
#include "CGOptimize.H"
#include "IceVelocityControlProblem.H"
#include "IceVelocity.H"
#include "FillFromReference.H"
#include "LevelDataIBC.H"
#include "ParmParse.H"


void ParseConstants(LevelSigmaCS& a_coords)
  {
    
    ParmParse ppCon("constants");
    Real iceDensity = 910;
    ppCon.query("ice_density",iceDensity); 
    a_coords.setIceDensity(iceDensity);
    
    Real seaWaterDensity=1028;
    ppCon.query("sea_water_density",seaWaterDensity);
    a_coords.setWaterDensity(seaWaterDensity);

    Real gravity=9.81;
    ppCon.query("gravity",gravity);
    a_coords.setGravity(gravity);

  }

/// define and solve a control problem .
/// 
/// miminise L = (u - u_obs)^2 + q (fwd(u)) with respect to C
///
/// the foward problem is solved for a velcity field u 
/// 
/// div(H*mu(u)*(grad(u)+...)) - C nu(u)*u = rho * g * H * grad(s) 
///
/// the adjoint problem is solved for a largange multiplier q
/// adj(q) = 0
/// div(H*mu(u)*(grad(q)+...)) - C nu(u)*q = rhs_adj
///
/// then dL/dC = q * u 
/// a_constRelPtr : constitutuve relation for mu calculation
/// a_basalFrictionRealPtr : constitutuve relation for nu calculation
/// a_levelThck : LevelData<FArrayBox>& of ice thickness, spans the domain
/// a_levelTopg : topography data
/// a_levelVelObs : velocity observations
/// a_levelC : on entry, initial guess for C. on succesful exit, 
///           a C which reduces (u - u_obs)^2.           
void solveControl(LevelData<FArrayBox>& a_levelC,
		  const LevelData<FArrayBox>& a_levelVelObs,
		  const LevelData<FArrayBox>& a_levelVelCoef,
		  const LevelData<FArrayBox>& a_levelDivUHObs,
		  const ProblemDomain& a_problemDomain,
		  const DisjointBoxLayout& a_grids,
#if BISICLES_Z == BISICLES_LAYERED
		  const Vector<Real>& a_faceSigma,
#endif
		  IceThicknessIBC* a_bc,
		  ConstitutiveRelation* a_constRelPtr, 
		  BasalFrictionRelation* a_basalFrictionRelPtr,
		  const RealVect& a_baseDx, const RealVect& a_dataDx)
{

  
 
  
  ParmParse pp("control");

  Real velMisfitCoefficient = 1.0;
  pp.query("velMisfitCoefficient",velMisfitCoefficient);

  Real massImbalanceCoefficient = 0.0;
  pp.query("massImbalanceCoefficient",massImbalanceCoefficient);

  Real gradCsqRegularization = 0.0;
  pp.query("gradCsqRegularization",gradCsqRegularization);

  Real gradMuCoefsqRegularization = gradCsqRegularization;
  pp.query("gradMuCoefsqRegularization",gradCsqRegularization);
  

  Vector<DisjointBoxLayout> vectGrids(1); vectGrids[0] = a_grids;
  Vector<int> refRatios(1); refRatios[0] = 2;
  Vector<RefCountedPtr<LevelSigmaCS > > vectCoordSys(1);
  vectCoordSys[0] = RefCountedPtr<LevelSigmaCS>
    (new LevelSigmaCS(vectGrids[0], a_baseDx, IntVect::Unit));
  a_bc->initializeIceGeometry(*vectCoordSys[0],a_baseDx,RealVect::Zero,0.0,NULL,2);
  a_bc->setGeometryBCs(*vectCoordSys[0],a_problemDomain,a_baseDx,0.0,0.0);
  ParseConstants(*vectCoordSys[0]);
  vectCoordSys[0]->recomputeGeometry(NULL,0);
#if BISICLES_Z == BISICLES_LAYERED
  vectCoordSys[0]->setFaceSigma(a_faceSigma);
#endif
 

  Vector<LevelData<FArrayBox>*> vectC(1,NULL); 
  vectC[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);

  Vector<LevelData<FArrayBox>*> vectC0(1,NULL); 
  vectC0[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);

  Vector<LevelData<FArrayBox>*> vectRhs(1,NULL); 
  vectRhs[0] =  new LevelData<FArrayBox>(vectGrids[0],SpaceDim,IntVect::Zero);

  Vector<LevelData<FArrayBox>*> vectVel(1,NULL); 
  vectVel[0] =  new LevelData<FArrayBox>(vectGrids[0],SpaceDim,IntVect::Unit);

  Vector<LevelData<FArrayBox>*> vectVelObs(1,NULL); 
  vectVelObs[0] =  new LevelData<FArrayBox>(vectGrids[0],SpaceDim,IntVect::Unit);

  Vector<LevelData<FArrayBox>*> vectVelCoef(1,NULL); 
  vectVelCoef[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);
 
  Vector<LevelData<FArrayBox>*> vectDivUHObs(1,NULL); 
  vectDivUHObs[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);

  Vector<LevelData<FArrayBox>*> vectA(1,NULL); 
#if BISICLES_Z == BISICLES_LAYERED
  vectA[0] =  new LevelData<FArrayBox>
    (vectGrids[0],a_faceSigma.size()-1,IntVect::Unit);   
#else
  vectA[0] =  new LevelData<FArrayBox>
    (vectGrids[0],1,IntVect::Unit); 
#endif
  Vector<LevelData<FArrayBox>*> vectCOrigin(1,NULL); 
  vectCOrigin[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);
  
  

  //parameters to be found : C and coefficient of mu
  Vector<LevelData<FArrayBox>*> vectParameter(1,NULL); 
  vectParameter[0] =  new LevelData<FArrayBox>(vectGrids[0],2,IntVect::Zero);

  Vector<RealVect> vectDx(1);
  vectDx[0] = a_baseDx;

  for (DataIterator dit(vectGrids[0]);dit.ok();++dit)
    {
      (*vectA[0])[dit].setVal(4.0e-17);
      //(*vectC[0])[dit].copy(a_levelC[dit]);
      (*vectC0[0])[dit].setVal(0.0);				     
      (*vectVel[0])[dit].setVal(0.0);
      (*vectParameter[0])[dit].setVal(0.0);
      //(*vectCOrigin[0])[dit].copy(a_levelC[dit]);
      //(*vectVelObs[0])[dit].copy(a_levelVelObs[dit]);
      //(*vectVelCoef[0])[dit].copy(a_levelVelCoef[dit]);
      //(*vectDivUHObs[0])[dit].copy(a_levelDivUHObs[dit]);
    }
  FillFromReference(*vectC[0],a_levelC,a_baseDx,a_dataDx,true);
  FillFromReference(*vectCOrigin[0],a_levelC,a_baseDx,a_dataDx,true);
  std::string restartfile("");
  pp.query("restartfile",restartfile);
  
  if (restartfile != "")
    {
      pout() << "restarting from " << restartfile << std::endl;
      Vector<std::string> restartNames;
      Vector<LevelData<FArrayBox>* > restartData;
      Vector<DisjointBoxLayout> restartGrids;
      Vector<int> restartRatio;
      Real restartDx = 0.0, restartDt = 0.0, restartTime = 0.0;
      Box restartDomBox;
      int restartNumLevels;
      int status = ReadAMRHierarchyHDF5
	(restartfile,restartGrids,restartData,restartNames,
	 restartDomBox,restartDx,restartDt,restartTime,
	 restartRatio,restartNumLevels);
      CH_assert(restartDx == a_baseDx[0]);
      CH_assert(status == 0);
     
      for (int i = 0; i < restartNames.size(); i++)
	{
	  if (restartNames[i] == "X0")
	    {
	      restartData[0]->copyTo(Interval(i,i),*vectParameter[0],Interval(0,0));
	    }
	  if (restartNames[i] == "X1")
	      {
		restartData[0]->copyTo(Interval(i,i),*vectParameter[0],Interval(1,1));
	      }
	}
    }
  
  
  
  FillFromReference(*vectVelObs[0],a_levelVelObs,a_baseDx,a_dataDx,true);
  FillFromReference(*vectVelCoef[0],a_levelVelCoef,a_baseDx,a_dataDx,true);
  FillFromReference(*vectDivUHObs[0],a_levelDivUHObs,a_baseDx,a_dataDx,true);

  IceVelocity::defineRHS(vectRhs, vectCoordSys, vectGrids, vectDx);
  {
    //give the first nonlinear solve in 
    //IceVelocityControl::computeObjectiveAngGradient 
    //a fighting chance by solving a linear problem first.
    //maybe move this into IceVelocityControl::computeObjectiveAngGradient
    IceVelocity::setFloatingC(vectC,vectCoordSys, vectGrids, 0.0);
    //linear solve
    solveForward(vectVel, true, vectC, vectC0, vectA, 
	       vectRhs, vectCoordSys, a_problemDomain, 
		 vectGrids, refRatios,a_bc, a_constRelPtr, 
		 a_basalFrictionRelPtr, a_baseDx);
  }
  
  IceVelocityControlProblem cp(vectVelObs, vectVelCoef, vectDivUHObs, 
			       vectVel, vectRhs, vectCOrigin,
			       vectC0, vectA, vectCoordSys, vectGrids, 
			       a_problemDomain, refRatios, a_baseDx,
			       a_bc, a_constRelPtr, a_basalFrictionRelPtr,
			       velMisfitCoefficient,
			       massImbalanceCoefficient,
			       gradCsqRegularization,
			       gradMuCoefsqRegularization);

  int CGmaxIter = 500;
  pp.query("CGmaxIter",CGmaxIter);

  Real CGtol = 1.0e-6;
  pp.query("CGtol",CGtol);

  Real CGsecantParameter = 1.0e-8;
  pp.query("CGsecantParameter",CGsecantParameter);

  int CGsecantMaxIter = 10;
  pp.query("CGsecantMaxIter",CGsecantMaxIter);

  Real CGsecantTol = 1.0e-6;
  pp.query("CGsecantTol",CGsecantTol);

  

  CGOptimize(cp , vectParameter , CGmaxIter , CGtol , 
	     CGsecantParameter, CGsecantMaxIter , CGsecantTol);

 

}

//construct and solve, for u, the forward problem
// div(H*mu(u)*(grad(u)+...)) - C nu(u)*u = rho * g * H * grad(s)
void solveForward(Vector<LevelData<FArrayBox>* >& a_u,
		  bool a_linear,
		  const Vector<LevelData<FArrayBox>* >& a_C,
		  const Vector<LevelData<FArrayBox>* >& a_C0,
		  const Vector<LevelData<FArrayBox>* >& a_A,
		  const Vector<LevelData<FArrayBox>* >& a_rhs,
		  Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		  const ProblemDomain& a_crseDomain,
		  const Vector<DisjointBoxLayout>& a_grids,
		  const Vector<int> a_refRatios,
		  IceThicknessIBC* a_bc,
		  ConstitutiveRelation* a_constRelPtr, 
		  BasalFrictionRelation* a_basalFrictionRelPtr,
		  const RealVect& a_crseDx)
{

  int maxLevel = a_u.size()-1;
  
  constMuRelation* constMuPtr = NULL;
  ConstitutiveRelation* constRelPtr;
  if (a_linear)
    {
       constMuPtr = new constMuRelation;
       constMuPtr->setConstVal(1.0e+7);
       //jfnkSolver.m_maxIter = 1;
       //jfnkSolver.m_RelaxRelTol = 1.0e-10;
       //jfnkSolver.m_maxRelaxIter = 100;
       constRelPtr = constMuPtr;
    }
  else
    {
      constRelPtr = a_constRelPtr;
    }

  JFNKSolver jfnkSolver;
  jfnkSolver.define(a_crseDomain, constRelPtr , a_basalFrictionRelPtr,
		    a_grids, a_refRatios, a_crseDx, a_bc, maxLevel+1);
    
  Real initialNorm = 10; Real finalNorm = 1.0; Real convergenceMetric=-1.0;

  Vector<LevelData<FluxBox>* > muCoef(maxLevel + 1,NULL);

  jfnkSolver.solve(a_u, initialNorm, finalNorm, convergenceMetric, a_linear, 
		   a_rhs, a_C, a_C0, a_A, muCoef, a_coordSys, 0.0 , 0, maxLevel);


  if (a_linear)
    {
      delete constMuPtr;
    }

}



//construct and solve, for u, the forward problem
//div(H*mu(u)*(grad(u)+...)) - C nu(u)*u = rho * g * H * grad(s)
void solveForward(LevelData<FArrayBox>& a_levelVel,
		  LevelData<FArrayBox>& a_levelDivUH,
		  const LevelData<FArrayBox>& a_levelC,
		  const ProblemDomain& a_problemDomain,
		  const DisjointBoxLayout& a_grids,
#if BISICLES_Z == BISICLES_LAYERED
		  const Vector<Real>& a_faceSigma,
#endif
		  IceThicknessIBC* a_bc,
		  ConstitutiveRelation* a_constRelPtr, 
		  BasalFrictionRelation* a_basalFrictionRelPtr,
		  const RealVect& a_baseDx, const RealVect& a_dataDx)
{

  Vector<DisjointBoxLayout> vectGrids(1); vectGrids[0] = a_grids;
  Vector<int> refRatios(1); refRatios[0] = 2;
  Vector<RefCountedPtr<LevelSigmaCS > > vectCoordSys(1);
  vectCoordSys[0] = RefCountedPtr<LevelSigmaCS>
    (new LevelSigmaCS(vectGrids[0], a_baseDx, IntVect::Unit));
  a_bc->initializeIceGeometry(*vectCoordSys[0],a_baseDx,RealVect::Zero,0.0,NULL,2);
  a_bc->setGeometryBCs(*vectCoordSys[0],a_problemDomain,a_baseDx,0.0,0.0);
  ParseConstants(*vectCoordSys[0]);
  vectCoordSys[0]->recomputeGeometry(NULL,0);
#if BISICLES_Z == BISICLES_LAYERED
  vectCoordSys[0]->setFaceSigma(a_faceSigma);
#endif
  

  Vector<LevelData<FArrayBox>*> vectC(1,NULL); 
  vectC[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);
  Vector<LevelData<FArrayBox>*> vectC0(1,NULL); 
  vectC0[0] =  new LevelData<FArrayBox>(vectGrids[0],1,IntVect::Zero);
  Vector<LevelData<FArrayBox>*> vectRhs(1,NULL); 
  vectRhs[0] =  new LevelData<FArrayBox>(vectGrids[0],SpaceDim,IntVect::Zero);
  Vector<LevelData<FArrayBox>*> vectA(1,NULL); 
  Vector<ProblemDomain> vectDomain(1,a_problemDomain);

#if BISICLES_Z == BISICLES_LAYERED
  vectA[0] =  new LevelData<FArrayBox>
    (vectGrids[0],a_faceSigma.size()-1,IntVect::Unit);   
#else
  vectA[0] =  new LevelData<FArrayBox>
    (vectGrids[0],1,IntVect::Unit); 
#endif
  Vector<RealVect> vectDx(1);
  vectDx[0] = a_baseDx;

  Vector<LevelData<FArrayBox>*> vectVel(1,NULL); 
  vectVel[0] = new LevelData<FArrayBox>(vectGrids[0],SpaceDim,IntVect::Unit);
  
  FillFromReference(*vectC[0],a_levelC,a_baseDx,a_dataDx,true);


  for (DataIterator dit(vectGrids[0]);dit.ok();++dit)
    {
      (*vectA[0])[dit].setVal(4.0e-17);
      (*vectC0[0])[dit].setVal(0.0);				     
      (*vectVel[0])[dit].setVal(0.0);
      
    }

  //set up
  IceVelocity::defineRHS(vectRhs, vectCoordSys, vectGrids, vectDx);
  IceVelocity::setFloatingC(vectC,vectCoordSys, vectGrids, 0.0);
  //linear solve
  solveForward(vectVel, true, vectC, vectC0, vectA, 
	       vectRhs, vectCoordSys, a_problemDomain, 
	       vectGrids, refRatios,a_bc, a_constRelPtr, 
	       a_basalFrictionRelPtr, a_baseDx);
  
  //nonlinear solve
  solveForward(vectVel, false, vectC, vectC0, vectA, 
	       vectRhs, vectCoordSys, a_problemDomain, 
	       vectGrids, refRatios,a_bc, a_constRelPtr, 
	       a_basalFrictionRelPtr, a_baseDx);

  LevelData<FluxBox> flux(a_grids,1,IntVect::Zero);

  FillFromReference(a_levelVel,*vectVel[0],a_dataDx,a_baseDx,true);

  //IceVelocity::computeFaceFlux
  //  (flux,a_levelVel,vectCoordSys[0]->getH(),a_grids);
  //IceVelocity::applyDiv(a_levelDivUH,flux,a_grids,a_dx);
  pout() << "solveForwardProblem done";
}


//fill an AMR hierarchy from a single level of data
//and a vector of refinement ratios.
//The level data resides on level a_level of the hierarchy. 
void initData(Vector<RefCountedPtr<LevelSigmaCS > >& a_coords,
	      Vector<LevelData<FArrayBox>*>& a_C,
	      Vector<LevelData<FArrayBox>*>& a_A,
	      Vector<LevelData<FArrayBox>*>& a_velObs,
	      Vector<LevelData<FArrayBox>*>& a_vel,
	      Vector<LevelData<FArrayBox>*>& a_DivUHObs,
	      Vector<DisjointBoxLayout>& a_grids,
	      Vector<ProblemDomain>& a_domain,
	      Vector<RealVect>& a_dx,
	      const int& a_maxLevel,
	      const Vector<int>& a_refRatio,
	      const int& a_level,
	      LevelDataIBC* a_bc,
	      LevelData<FArrayBox>& a_levelC,
	      const LevelData<FArrayBox>& a_levelVelObs,
	      const LevelData<FArrayBox>& a_levelVelCoef,
	      const LevelData<FArrayBox>& a_levelDivUHObs,
#if BISICLES_Z == BISICLES_LAYERED
	      const Vector<Real>& a_faceSigma,
#endif	     
	      const RealVect& a_crseDx,
	      const ProblemDomain& a_crseDomain,
	      const DisjointBoxLayout& a_crseGrids)
{



}
