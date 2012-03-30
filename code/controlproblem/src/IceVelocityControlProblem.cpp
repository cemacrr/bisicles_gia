#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IceVelocityControlProblem.H"
#include "BoxIterator.H"
#include "JFNKSolver.H"
#include "ReflectGhostCells.H"
#include "IceVelocity.H"
#include "IceConstants.H"
#include "ControlF_F.H"
//#include "computeNorm.H"
#include <sstream>
#include "NamespaceHeader.H"



IceVelocityControlProblem::IceVelocityControlProblem
(const Vector<LevelData<FArrayBox>* >& a_velObs,
 const Vector<LevelData<FArrayBox>* >& a_velCoef,
 const Vector<LevelData<FArrayBox>* >& a_divUHObs,
 const Vector<LevelData<FArrayBox>* >& a_vel,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 const Vector<LevelData<FArrayBox>* >& a_C_origin,
 const Vector<LevelData<FArrayBox>* >& a_C0,
 const Vector<LevelData<FArrayBox>* >& a_A,
 Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<DisjointBoxLayout>& a_grids,
 const ProblemDomain& a_crseDomain,
 const Vector<int>& a_refRatios,
 const RealVect& a_crseDx,
 IceThicknessIBC* a_bcPtr,
 ConstitutiveRelation* a_constRelPtr,
 BasalFrictionRelation* a_basalFrictionRelPtr,
 const Real& a_velMisfitCoefficient,
 const Real& a_massImbalanceCoefficient,
 const Real& a_gradCsqRegularization,
 const Real& a_gradMuCoefsqRegularization
)
  :m_velObs(a_velObs), m_velCoef(a_velCoef), 
   m_divUHObs(a_divUHObs), 
   m_vel(a_vel), m_rhs(a_rhs),
   m_C_origin(a_C_origin),m_C0(a_C0),m_A(a_A),
   m_bcPtr(a_bcPtr),m_constRelPtr(a_constRelPtr),
   m_basalFrictionRelPtr(a_basalFrictionRelPtr),
   m_coordSys(a_coordSys),m_grids(a_grids),
   m_crseDomain(a_crseDomain),m_crseDx(a_crseDx),
   m_refRatios(a_refRatios), m_nComputeGradientCalls(0),
   m_velMisfitCoefficient(a_velMisfitCoefficient),
   m_massImbalanceCoefficient(a_massImbalanceCoefficient),
   m_gradCsqRegularization(a_gradCsqRegularization),
   m_gradMuCoefsqRegularization(a_gradMuCoefsqRegularization)
   
{
  int numLevels = m_velObs.size();
  m_vectOps.resize(numLevels);
  //create(m_C,m_C_origin);
  // C needs a ghost cell, becaue we want to compute its Laplacian
  m_C.resize(numLevels,NULL);
  for (int lev = 0; lev < numLevels; ++lev)
    {
      m_C[lev] = new LevelData<FArrayBox>( a_grids[lev],1,IntVect::Unit);
    }
 
  
  m_muCoef.resize(numLevels,NULL);
  m_faceMuCoef.resize(numLevels,NULL);
  m_faceA.resize(numLevels,NULL);
 
  for (int lev = 0; lev < numLevels; ++lev)
    {
      m_muCoef[lev] = new LevelData<FArrayBox>( a_grids[lev],1,IntVect::Unit);
      m_faceMuCoef[lev] = new LevelData<FluxBox>( a_grids[lev],1,IntVect::Zero);
      m_faceA[lev] = new LevelData<FluxBox>( a_grids[lev],m_A[lev]->nComp(),IntVect::Zero);
      CellToEdge(*m_A[lev],*m_faceA[lev]);
    }

  create(m_zero,m_C_origin);
  setToZero(m_zero);

  create(m_lapC,m_C_origin);
  create(m_lapMuCoef,m_C_origin);
  create(m_divUH,m_C_origin);
  create(m_adjVel,m_vel);
  create(m_adjRhs,m_rhs);
  create(m_wk2compA,m_rhs);

  CH_assert(m_velMisfitCoefficient >= 0.0);
  CH_assert(m_massImbalanceCoefficient >= 0.0);
  CH_assert(m_gradCsqRegularization >= 0.0);
  CH_assert(m_gradMuCoefsqRegularization >= 0.0);

}

IceVelocityControlProblem::~IceVelocityControlProblem()
{
  free(m_C);
  free(m_adjVel);
  free(m_adjRhs);
  free(m_wk2compA);
  free(m_zero);
  free(m_muCoef);
  free(m_faceA);
  free(m_faceMuCoef);
    
}

void IceVelocityControlProblem::solveForwardProblem
(Vector<LevelData<FArrayBox>* >& a_u,
 bool a_linear,
 const Vector<LevelData<FArrayBox>* >& a_C,
 const Vector<LevelData<FArrayBox>* >& a_C0,
 const Vector<LevelData<FArrayBox>* >& a_A,
 const Vector<LevelData<FluxBox>* >& a_muCoef,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const ProblemDomain& a_crseDomain,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatios,
 IceThicknessIBC* a_bc,
 ConstitutiveRelation* a_constRelPtr, 
 BasalFrictionRelation* a_basalFrictionRelPtr,
 const RealVect& a_crseDx)
{
  
  int maxLevel = a_u.size()-1;
  JFNKSolver jfnkSolver;
  jfnkSolver.define(a_crseDomain, a_constRelPtr , a_basalFrictionRelPtr,
		    a_grids, a_refRatios, a_crseDx, a_bc, maxLevel+1);
    
  Real initialNorm = 1.0; Real finalNorm = 1.0; Real convergenceMetric=-1.0;

  jfnkSolver.solve(a_u, initialNorm, finalNorm, convergenceMetric, a_linear, 
		   a_rhs, a_C, a_C0, a_A, a_muCoef, a_coordSys, 0.0 , 0, maxLevel);

}

void IceVelocityControlProblem::solveAdjointProblem
(Vector<LevelData<FArrayBox>* >& a_u,
 const Vector<LevelData<FArrayBox>* >& a_C,
 const Vector<LevelData<FArrayBox>* >& a_C0,
 const Vector<LevelData<FArrayBox>* >& a_A,
 const Vector<LevelData<FluxBox>* >& a_muCoef,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const ProblemDomain& a_crseDomain,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatios,
 IceThicknessIBC* a_bc,
 ConstitutiveRelation* a_constRelPtr, 
 BasalFrictionRelation* a_basalFrictionRelPtr,
 const RealVect& a_crseDx)
{
  int maxLevel = a_u.size()-1;
  JFNKSolver jfnkSolver;
  jfnkSolver.define(a_crseDomain, a_constRelPtr , a_basalFrictionRelPtr,
		    a_grids, a_refRatios, a_crseDx, a_bc, maxLevel+1);
    
  Real initialNorm = 1.0; Real finalNorm = 1.0; Real convergenceMetric=1.0;
  jfnkSolver.m_RelaxRelTol = 1.0e-10;
  jfnkSolver.solve(a_u, initialNorm, finalNorm, convergenceMetric, true, 
		   a_rhs, a_C, a_C0, a_A, a_muCoef, a_coordSys, 0.0 , 0, maxLevel);

}

//members required by CGOptimize

void IceVelocityControlProblem::computeObjectiveAndGradient
(Real& a_f, Vector<LevelData<FArrayBox>* >& a_g, 
 const  Vector<LevelData<FArrayBox>* >& a_x, bool a_inner)
{
  m_nComputeGradientCalls++;
  a_f = 0;
  //compute C = C0 * exp(a_x[0]) and muCoef = exp(a_x[1])
  for (int lev=0; lev < a_x.size();lev++)
    {
      const LevelData<FArrayBox>& levelX =  *a_x[lev];
      LevelData<FArrayBox>& levelC =  *m_C[lev];
      LevelData<FArrayBox>& levelC_origin =  *m_C_origin[lev];
      LevelData<FArrayBox>& levelMuCoef =  *m_muCoef[lev];
      const DisjointBoxLayout levelGrids =  m_grids[lev];
      levelC_origin.copyTo(levelC);
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  FArrayBox& thisMuCoef = levelMuCoef[dit];
	  thisMuCoef.setVal(1.0);
	}

      //boundary values for C : extrapolation or reflection should be sufficient
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  if (! m_crseDomain.isPeriodic(dir))
	    {
	      ReflectGhostCells(levelC, m_crseDomain, dir, Side::Lo);
	      ReflectGhostCells(levelC, m_crseDomain, dir, Side::Hi);
	      ReflectGhostCells(levelMuCoef, m_crseDomain, dir, Side::Lo);
	      ReflectGhostCells(levelMuCoef, m_crseDomain, dir, Side::Hi);
	    }
	}
      
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  FArrayBox& thisC = levelC[dit];
	  FArrayBox& thisMuCoef = levelMuCoef[dit];
	  Real boundArg = 2.0;
	  for (BoxIterator bit(levelGrids[dit]);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      Real arg = std::min(boundArg,std::max(levelX[dit](iv,0),-boundArg));
	      thisC(iv) *= exp(arg);

	      arg = std::min(boundArg,std::max(levelX[dit](iv,1),-boundArg));
	      thisMuCoef(iv) *= exp(arg);
	    }
	  //thisC *= levelC_origin[dit];
	  pout() << " min(C) = " << thisC.min(levelGrids[dit]) 
		 << " max(C) = " << thisC.max(levelGrids[dit]) << endl; 
	  
	  pout() << " min(muCoef) = " << thisMuCoef.min(levelGrids[dit]) 
		 << " max(muCoef) = " << thisMuCoef.max(levelGrids[dit]) << endl; 

	}
      levelC.exchange();
      levelMuCoef.exchange();

      CellToEdge( levelMuCoef, *m_faceMuCoef[lev]);

      //\todo need to make calculation of Laplacian(C) and Laplacian(MuCoef) 
      // work for multiple levels
      // mainly need to QuadCFInterp ghost cells of levelC to do so
      // also set level dx correctky
      RealVect levelDx = m_crseDx;
      LevelData<FArrayBox>& levelLapC =  *m_lapC[lev];
      IceVelocity::applyHelmOp(levelLapC, levelC, 0.0, 1.0, m_grids[lev],levelDx);
      LevelData<FArrayBox>& levelLapMuCoef =  *m_lapMuCoef[lev];
      IceVelocity::applyHelmOp(levelLapMuCoef, levelMuCoef, 0.0, 1.0, m_grids[lev],levelDx);

    }
  
  //solve forward problem
  IceVelocity::setFloatingC(m_C,m_coordSys, m_grids, 0.0);
  solveForwardProblem
    (m_vel,false,m_C,m_C0,m_A,m_faceMuCoef,m_rhs,m_coordSys,
     m_crseDomain, m_grids, m_refRatios, m_bcPtr,
     m_constRelPtr, m_basalFrictionRelPtr, m_crseDx);
  
  for (int lev=0; lev < a_x.size();lev++)
   {
     //\todo need to make calcualtion of div(uH) work for multiple levels
     //IceVelocity::applyDiv(levelDivUH,levelVel);
     LevelData<FArrayBox>& levelDivUH =  *m_divUH[lev];
     const LevelSigmaCS& levelCS = *m_coordSys[lev];
     const LevelData<FArrayBox>& levelVel =  *m_vel[lev];
     RealVect levelDx = m_crseDx;
     LevelData<FluxBox> flux(m_grids[lev],1,IntVect::Zero);
     IceVelocity::computeFaceFlux
       (flux,levelVel,levelCS.getH(),m_grids[lev]);
     IceVelocity::applyDiv(levelDivUH,flux,m_grids[lev],levelDx);

   }
  //compute RHS for adjoint problem, copy vel into adjVel
  for (int lev=0; lev < a_x.size();lev++)
    {
      const LevelData<FArrayBox>& levelVelObs =  *m_velObs[lev];
     LevelData<FArrayBox>& levelVelCoef =  *m_velCoef[lev];
      const LevelData<FArrayBox>& levelVel =  *m_vel[lev];
      LevelData<FArrayBox>& levelDivUHObs =  *m_divUHObs[lev];
      LevelData<FArrayBox>& levelDivUH =  *m_divUH[lev];
      const LevelSigmaCS& levelCS = *m_coordSys[lev];
      LevelData<FArrayBox>& levelAdjRhs = *m_adjRhs[lev];
      LevelData<FArrayBox>& levelAdjVel = *m_adjVel[lev];
    
#ifdef FITVELCOMP  
      // contribution due to velocity misfit
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{

	  FArrayBox& adjRhs = levelAdjRhs[dit];
	  
	  //contribution due to velocity misfit
	  levelAdjRhs[dit].copy(levelVelObs[dit]);
	  levelAdjRhs[dit].minus(levelVel[dit]);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      for (BoxIterator bit(levelVelCoef[dit].box());bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if (levelVelCoef[dit](iv) < 0.975)
		    levelVelCoef[dit](iv) = 0.0;
		}


	      levelAdjRhs[dit].mult(levelVelCoef[dit],0,dir);
	    }
	  pout() << " min(velocity misfit) = " << adjRhs.min() 
		 << " max(velocity misfit) = " << adjRhs.max() << endl;
	  levelAdjRhs[dit] *= m_velMisfitCoefficient;

	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      a_f += levelAdjRhs[dit].norm(dir,box,2,dir,1)
	    }
	}
#else     
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  
	  FArrayBox& adjRhs = levelAdjRhs[dit];
	  const FArrayBox& um = levelVel[dit];
	  const FArrayBox& uo = levelVelObs[dit];
	  const Box& box = (m_grids[lev])[dit];
	  FArrayBox misfit(box,1);

	  FORT_ADJRHSSPEEDCTRL(CHF_FRA1(adjRhs,0), CHF_FRA1(adjRhs,1),
			       CHF_CONST_FRA1(misfit,0),
			       CHF_CONST_FRA1(um,0), CHF_CONST_FRA1(um,1),
			       CHF_CONST_FRA1(uo,0), CHF_CONST_FRA1(uo,1),
			       CHF_BOX(box));    
			     
	  misfit.mult(levelVelCoef[dit]);

	  for (BoxIterator bit(levelVelCoef[dit].box());bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (levelVelCoef[dit](iv) < 0.975)
		levelVelCoef[dit](iv) = 0.0;
	    }


	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      levelAdjRhs[dit].mult(levelVelCoef[dit],0,dir);
	    }

	 
	  pout() << " max(velocity misfit) = " << std::sqrt(misfit.max()) << endl;
	  levelAdjRhs[dit] *= m_velMisfitCoefficient;
	  misfit *= m_velMisfitCoefficient;
	  a_f += misfit.norm(box,2,0,1);
	  
	}  


#endif

      //contribution due to mass imbalance
      {
	for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  // \todo : make this work for mutliple levels/FABs
	  const Box& box = m_grids[lev][dit];
	  FArrayBox t(box,1);
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];

	  t.copy(levelDivUHObs[dit]);t.minus(levelDivUH[dit]);

	  a_f += m_massImbalanceCoefficient * t.norm(box,2,0,1);

	  RealVect oneOnTwoDx = 1.0/ (2.0 * m_crseDx);

	  Box sbox(box);
	  sbox.grow(-1);
	  FArrayBox grad(sbox,2);
	  for (int dir =0; dir < SpaceDim; ++dir)
	  grad.setVal(0.0,dir);
	  for (BoxIterator bit(sbox);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (mask(iv) != GROUNDEDMASKVAL)
		{
		  t(iv) = 0.0;
		}
	      else 
		{
		  for (int dir =0; dir < SpaceDim; ++dir)
		    {
		      IntVect ivp = iv + BASISV(dir);
		      IntVect ivm = iv - BASISV(dir);
		      grad(iv,dir) = ( t(ivp) - t(ivm) ) * oneOnTwoDx[dir];
		    }

		}
	    }
	  pout() << " min(mass imbalance) = " << t.min() 
		 << " max(mass imbalance) = " << t.max() << endl; 

	  const FArrayBox& H = levelCS.getH()[dit];
	  grad.mult(H,0,0); grad.mult(H,0,1);
	  grad *= m_massImbalanceCoefficient;
	  for (int dir =0; dir < SpaceDim; ++dir)
	    grad.mult(levelVelCoef[dit],0,dir);
	  levelAdjRhs[dit] -= grad;
	}

	

	for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  // initialize adjvel to vel, tm ensure that the
	  // correct coeffiecients are computed
	  levelAdjVel[dit].copy(levelVel[dit]);
	}
      }
    }

  //solve adjoint problem
  solveAdjointProblem
    (m_adjVel,m_C,m_C0,m_A,m_faceMuCoef,m_adjRhs,m_coordSys,
      m_crseDomain, m_grids, m_refRatios, m_bcPtr,
      m_constRelPtr, m_basalFrictionRelPtr, m_crseDx);

  //compute gradient with respct to a_x component 0 
  // = - lambda * u * C - b^2 * lap(C)
  for (int lev=0; lev < a_x.size();lev++)
    {
      const LevelData<FArrayBox>& levelLambda =  *m_adjVel[lev];
      const LevelData<FArrayBox>& levelU =  *m_vel[lev];
      //const LevelData<FArrayBox>& levelX =  *a_x[lev];
      const LevelData<FArrayBox>& levelC =  *m_C[lev];
      LevelData<FArrayBox>& levelLapC =  *m_lapC[lev];
      LevelData<FArrayBox>& levelG =  *a_g[lev];
      for (DataIterator dit(levelU.dataIterator());dit.ok();++dit)
	{
	  const FArrayBox& thisU = levelU[dit];
	  const FArrayBox& thisC = levelC[dit];
	  //const FArrayBox& thisX = levelX[dit];
	  const FArrayBox& thisLambda = levelLambda[dit];
	  
	  FArrayBox& thisG =  levelG[dit];
	  thisG.setVal(0.0,0);
	  
	  //gradient with respect to a_x component 0,
          //where C = Corigin * exp(a_x)
	  FArrayBox t(thisG.box(),1);
	  t.copy(thisU,0,0);
	  t.mult(thisLambda,0,0);
	  thisG.minus(t,0,0);
	  t.copy(thisU,1,0);
	  t.mult(thisLambda,1,0);
	  thisG.minus(t,0,0);
	  thisG.mult(thisC,0,0);

	

	  //regularization : penalise large {grad(C)}^2
	  FArrayBox& thisLapC = levelLapC[dit];
	  t.copy(thisLapC); t*= thisC; t*= -m_gradCsqRegularization;
	  thisG.plus(t,0,0);

	}
    }

  //compute gradient with respect to a_x component 1 
  {
    int nlev = m_vel.size();
    Vector<ProblemDomain> domain(nlev);
    Vector<RealVect> dx(nlev);
    domain[0] = m_crseDomain;
    dx[0] = m_crseDx;
    for (int lev = 1; lev < nlev; lev++)
      {
	domain[lev] = domain[lev-1]; domain[lev].refine(m_refRatios[lev-1]);
	dx[lev] = dx[lev-1] / Real( m_refRatios[lev-1]);
      }

    IceJFNKstate state(m_grids, m_refRatios, domain, dx,
		       m_coordSys, m_vel, m_zero, m_zero, nlev-1,
		       *m_constRelPtr, *m_basalFrictionRelPtr,
		       *m_bcPtr, m_A, m_faceA, 0.0, 0.0, 0, 0.0);
    state.setState(m_vel);

    Vector<LevelData<FluxBox>*> vtFace(nlev,NULL);
    for (int lev = 0; lev < nlev; lev++)
      {
	vtFace[lev] = new LevelData<FluxBox>(m_grids[lev], SpaceDim, IntVect::Unit);
      }
    state.computeViscousTensorFace(vtFace);

    for (int lev = 0; lev < nlev; lev++)
      {
	LevelData<FArrayBox>& levelG =  *a_g[lev];
	const LevelData<FluxBox>& levelVT =  *vtFace[lev];
	LevelData<FArrayBox>& levelMuCoef  =  *m_muCoef[lev];
	LevelData<FArrayBox>& levelLapMuCoef  =  *m_lapMuCoef[lev];
	const LevelData<FArrayBox>& levelX =  *a_x[lev];
	const LevelData<FArrayBox>& levelLambda =  *m_adjVel[lev];
	

	for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	  {
	    const FArrayBox& thisMuCoef = levelMuCoef[dit];
	   
	    const Box& box = m_grids[lev][dit]; 
	    FArrayBox t(box,1);
	    t.setVal(0.0);
	    const FArrayBox& lambda = levelLambda[dit];
	    const FluxBox& vt = levelVT[dit];
	    
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

	    RealVect oneOnTwoDx = 1.0/ (2.0 * m_crseDx);
	    t.mult(-oneOnTwoDx[0]);
	    t.mult(thisMuCoef);

	    FArrayBox& thisG =  levelG[dit];
	    thisG.copy(t,0,1);
	    

	    //regularization : penalise large X[1]
	    //const FArrayBox& thisX = levelX[dit];
	    //t.copy(thisX,1,0); t *= 1e+6 * m_gradMuCoefsqRegularization;
	    //thisG.plus(t,0,1);
	    //regularization : penalise large {grad(MuCoef)}^2
	    const FArrayBox& thisLapMuCoef = levelLapMuCoef[dit];
	    t.copy(thisLapMuCoef); t*= thisMuCoef; t*= -m_gradMuCoefsqRegularization;
	    thisG.plus(t,0,1);
	  }
      }

    for (int lev = 0; lev < nlev; lev++)
      {
	delete vtFace[lev]; vtFace[lev] = NULL;
      }

    
  }
  


   {
     //output current state
     
     //recompute C (to see what is happening under the shelf
     Vector<LevelData<FArrayBox>*> Cplus(a_x.size(),NULL);
     for (int lev=0; lev < a_x.size();lev++)
       {
	 const LevelData<FArrayBox>& levelX =  *a_x[lev];
	 LevelData<FArrayBox>& levelC_origin =  *m_C_origin[lev];
	 const DisjointBoxLayout levelGrids =  m_grids[lev];

	 Cplus[lev] = new LevelData<FArrayBox>(levelGrids,1,IntVect::Zero);
	 LevelData<FArrayBox>& levelC = *Cplus[lev];
	 for (DataIterator dit(levelGrids);dit.ok();++dit)
	   {
	     FArrayBox& thisC = levelC[dit];
	     for (BoxIterator bit(levelGrids[dit]);bit.ok();++bit)
	       {
		 const IntVect& iv = bit();
		 thisC(iv) = exp(levelX[dit](iv));
	    }
	     thisC *= levelC_origin[dit];
	   }
	 levelC.exchange();
       }

    std::stringstream ss;
    ss << "out";
    ss.width(6);ss.fill('0');ss << m_nComputeGradientCalls - 1;
    ss.width(0);ss << ".2d.hdf5";
    Vector<std::string> names;
    names.push_back("X0");
    names.push_back("X1");
    names.push_back("C");
    names.push_back("muCoef");
    names.push_back("Cplus");
    names.push_back("xVel");
    names.push_back("yVel");
    names.push_back("xVelo");
    names.push_back("yVelo");
    names.push_back("divuh");
    names.push_back("divuho");
    names.push_back("xLaMu");
    names.push_back("yLaMu");
    names.push_back("gradfC");
    names.push_back("gradfMuCoef");
    names.push_back("velc");
    LevelData<FArrayBox> data(m_grids[0],names.size(),IntVect::Zero);

    int j = 0;
    a_x[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;
    a_x[0]->copyTo(Interval(1,1),data,Interval(j,j));j++;
    m_C[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;
    m_muCoef[0]->copyTo(Interval(0,0),data,Interval(j,j));j++; 
    Cplus[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;
    m_vel[0]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
    m_velObs[0]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
    m_divUH[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;
    m_divUHObs[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;
    m_adjVel[0]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
    a_g[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;
    a_g[0]->copyTo(Interval(1,1),data,Interval(j,j));j++;
    m_velCoef[0]->copyTo(Interval(0,0),data,Interval(j,j));j++;

    delete Cplus[0];

    const Real dt = 1.0;
    const Real time = Real(m_nComputeGradientCalls) - 1.0;
    Vector<LevelData<FArrayBox>* > vectData(1,&data);
    WriteAMRHierarchyHDF5(ss.str(),m_grids,vectData,names, m_crseDomain.domainBox(),
			  m_crseDx[0], dt, time, m_refRatios, vectData.size());
			  
  }
 



  
}



//apply preconditioner s = M^{-1}r
void IceVelocityControlProblem::preCond
(Vector<LevelData<FArrayBox>* >& a_s, 
 const Vector<LevelData<FArrayBox>* >& a_r)
{
  assign(a_s,a_r);
}

// duplicate storage of a_b in a_a
void IceVelocityControlProblem::create
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

void IceVelocityControlProblem::free
(Vector<LevelData<FArrayBox>* >& a_a)
{
 
  for (int lev = 0; lev < a_a.size(); lev++)
    {
      if (a_a[lev] != NULL)
	{
	  delete a_a[lev];
	  a_a[lev] = NULL;
	}
    } 
}

void IceVelocityControlProblem::free
(Vector<LevelData<FluxBox>* >& a_a)
{
 
  for (int lev = 0; lev < a_a.size(); lev++)
    {
      if (a_a[lev] != NULL)
	{
	  delete a_a[lev];
	  a_a[lev] = NULL;
	}
    } 
}

  // set a_x = s * a_x
void IceVelocityControlProblem::setToZero
(Vector<LevelData<FArrayBox>* >& a_x)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].setToZero(*a_x[lev]);
}

  // set a_x = s * a_x
void IceVelocityControlProblem::scale
(Vector<LevelData<FArrayBox>* >& a_x, 
 const  Real a_s)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].scale(*a_x[lev],a_s);
}

  // set a_y = a_x
void IceVelocityControlProblem::assign
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x)
{
  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].assign(*a_y[lev],*a_x[lev]);
}
  
// set a_y = a_y + a_s * a_x
void IceVelocityControlProblem::incr
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x, Real a_s)
{
for (int lev = 0; lev < a_x.size(); lev++)
  m_vectOps[lev].incr(*a_y[lev],*a_x[lev],a_s);
}
  
  // return a_y.a_x
Real IceVelocityControlProblem::dotProduct
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x)
{
  CH_assert((a_x.size() == 1) && (a_y.size() == a_x.size()));
  Real p = m_vectOps[0].dotProduct(*a_y[0],*a_x[0]);
  return p;
}

// how many degrees of freedom are there in x?
int IceVelocityControlProblem::nDoF(const Vector<LevelData<FArrayBox>* >& x)
{
  return 100;
}
#include "NamespaceFooter.H"
