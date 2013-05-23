#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "JFNKSolver.H"
#include "QuadCFInterp.H"
#include "CornerCopier.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"
#include "ParmParse.H"
#include "BisiclesF_F.H"
#include "LinearSolver.H"
#include "RelaxSolver.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "CGSolver.H"
#include <sstream>
#include "NamespaceHeader.H"

// static member initialization
int JFNKOp::m_bottom_solver_type = bicg;
int JFNKOp::m_MG_solver_depth = -1;

void JFNKSolver::imposeMaxAbs(Vector<LevelData<FArrayBox>*>& a_u,
		 Real a_limit)
{

  for (int lev = 0; lev < a_u.size(); ++lev)
    {
      LevelData<FArrayBox>& levelU = *a_u[lev];
      for (DataIterator dit (levelU.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  FArrayBox& thisU = levelU[dit];
	  FORT_ABSLIMITFAB(CHF_FRA(thisU), 
			   CHF_CONST_REAL(a_limit), 
			   CHF_BOX(thisU.box()));
	}
    }
  

}

int JFNKOp::m_residualID(0);

JFNKOp::JFNKOp(JFNKState* a_currentState , 
	       JFNKState* a_peturbedState, 
	       Vector<LevelData<FArrayBox>*>& a_u,
	       Real a_h,  Real a_err, Real a_umin, bool a_hAdaptive, 
	       Vector<DisjointBoxLayout>& a_grids,
	       Vector<int>& a_refRatio,
	       Vector<ProblemDomain>& a_domains,
	       Vector<RealVect>& a_dxs,
	       int a_lBase, 
	       int a_numMGSmooth,
	       int a_numMGIter,
	       SolverMode a_mode) 
  : m_u(a_currentState), m_uPerturbed(a_peturbedState),
    m_h(a_h), m_err(a_err), m_umin(a_umin), m_hAdaptive(a_hAdaptive),
    m_grids(a_grids), m_refRatio(a_refRatio),
    m_domains(a_domains), m_dxs(a_dxs), m_lBase(a_lBase),
    m_mode(a_mode),m_writeResiduals(false)
    
{
  CH_assert(a_currentState != a_peturbedState);
  
  Vector<LevelData<FArrayBox>*> localU( m_grids.size());
  for (int lev = 0; lev < m_grids.size(); ++lev)
    {
      localU[lev] = a_u[lev];
    }

  RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
    opFactoryPtr = m_u->opFactoryPtr();
  
  m_mlOp.m_bottom_solver_type = m_bottom_solver_type;
  m_mlOp.m_preCondSolverDepth = m_MG_solver_depth;
  m_mlOp.m_num_mg_smooth =  a_numMGSmooth;
  m_mlOp.m_num_mg_iterations = a_numMGIter;
  m_mlOp.define(m_grids , m_refRatio, m_domains, m_dxs, 
		opFactoryPtr ,a_lBase);
  
 
  RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
    pOpFactoryPtr = m_uPerturbed->opFactoryPtr();
  m_perturbedMlOp.define(m_grids , m_refRatio, m_domains, m_dxs, 
			 pOpFactoryPtr , a_lBase);
  
  
  create(m_fu,localU);
  if (m_mode == JFNK_SOLVER_MODE)
    create(m_uplushv,localU);
  setU(localU);
    
  }

void  JFNKOp::outerResidual(Vector<LevelData<FArrayBox>*>& a_lhs, 
			     const Vector<LevelData<FArrayBox>*>& a_u, 
			     const Vector<LevelData<FArrayBox>*>& a_rhs, 
			     bool a_homogeneous )
  {
    
    m_mlOp.applyOp(a_lhs, a_u, a_homogeneous);
    incr(a_lhs, a_rhs, -1);
    scale(a_lhs, -1.0);
    if (m_writeResiduals)
    {
      writeResidual(a_u,a_lhs);
    }
  }
  

void  JFNKOp::residual(Vector<LevelData<FArrayBox>*>& a_lhs, 
		       const Vector<LevelData<FArrayBox>*>& a_v, 
		       const Vector<LevelData<FArrayBox>*>& a_rhs, 
		       bool a_homogeneous)
{    
  applyOp(a_lhs, a_v, a_homogeneous);
  incr(a_lhs, a_rhs, -1.0);
  scale(a_lhs, -1.0);
  if (m_writeResiduals)
    {
      if (m_mode == JFNK_SOLVER_MODE)
	writeResidual(m_uplushv,a_lhs);
      else
	writeResidual(a_v,a_lhs);
    }

}


void JFNKOp::setU(Vector<LevelData<FArrayBox>*>& a_u)
{
  CH_TIME("JFNKOp::setU");
  m_u->setState(a_u);
  m_mlOp.applyOp(m_fu, a_u);
  
}
Real JFNKOp::finiteh(const Vector<LevelData<FArrayBox>*>& a_v)
{
  Real h = m_h;

  if (m_hAdaptive)
    {
      
      //calculation of h from PETSc SNES (see license notice at the end of this file)
      // h = err*u'v/||v||^2 if  |u'v| > umin*||v||_{1}
      //  = err*umin*sign(u'v)*||v||_{1}/||v||^2   otherwise
      Real utv = dotProduct(m_u->getState(), a_v); 
      Real scale = D_TERM(m_dxs[0][0], *m_dxs[0][1], *m_dxs[0][2]);
      utv /= scale;
      Real vL1norm = norm(a_v, 1);
      Real vL2norm = norm(a_v, 2);
      if (vL1norm > 0.0)
	{	  
	  if (utv >= 0.0 && utv <  m_umin*vL1norm) utv = m_umin*vL1norm;
	  else if (utv < 0.0 && utv > - m_umin*vL1norm) utv = -m_umin*vL1norm;
	  h = m_err * utv / (vL2norm *  vL2norm);
	}
      else
	{
	  h =  m_err;
	}
    }
  return h;

}


void JFNKOp::applyOp(Vector<LevelData<FArrayBox>*>& a_lhs, 
		       const Vector<LevelData<FArrayBox>*>& a_v, 
		       bool a_homogeneous )
{

  CH_TIME("JFNKOp::applyOp");
  
  if (m_mode == JFNK_SOLVER_MODE)
    {
      //"Newton mode"
      setToZero(m_uplushv);
      Real h = finiteh(a_v);
      axby(m_uplushv, m_u->getState(), a_v, 1.0, h);
      CH_assert(norm(m_uplushv,0) < HUGE_NORM);
      m_uPerturbed->setState(m_uplushv);
      m_perturbedMlOp.applyOp(a_lhs, m_uplushv , a_homogeneous);

      incr(a_lhs, m_fu, -1.0);
      scale(a_lhs, 1.0 / h);

    }
  else if (m_mode == PICARD_SOLVER_MODE)
    {
      //"Picard mode"
      m_mlOp.applyOp(a_lhs, a_v, a_homogeneous);
    }
  else 
    {
      MayDay::Error("Unknown SolverMode in JFNKOp::applyOp");
    }

}


void JFNKOp::writeResidual  
(const Vector<LevelData<FArrayBox> *>& a_u,
 const Vector<LevelData<FArrayBox> *>& a_residual)
{

  Vector<std::string> names;
  names.push_back("xU");
  names.push_back("yU");
  names.push_back("xRes");
  names.push_back("yRes");
  names.push_back("C");
  names.push_back("muSum");
  

  Vector<LevelData<FArrayBox>* > data(a_u.size(),NULL);
  for (int lev = 0; lev < a_u.size(); lev++)
    {
      data[lev] = new LevelData<FArrayBox>(m_grids[lev],int(names.size()),IntVect::Zero);
      int j = 0;
      a_u[lev]->copyTo(Interval(0,1),*data[lev],Interval(j,j+1)); j+=2;
      a_residual[lev]->copyTo(Interval(0,1),*data[lev],Interval(j,j+1)); j+=2;
      m_u->getDragCoef()[lev]->copyTo(Interval(0,0),*data[lev],Interval(j,j)); j+=1;

      const LevelData<FluxBox>& mu =  *m_u->getViscosityCoef()[lev];
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  FArrayBox muSum; 
	  muSum.define(  Interval(j,j),(*data[lev])[dit]);
	  for (BoxIterator bit( m_grids[lev][dit]);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      const IntVect ive = iv + BASISV(0);
	      const IntVect ivn = iv + BASISV(1); 
	      muSum(iv,0) = mu[dit][0](iv) + mu[dit][0](ive) 
	      	+ mu[dit][1](iv) +  mu[dit][1](ivn);
	    }
	}
    }

  // for (int lev = a_u.size() -1; lev > 0; lev--)
  //   {
  //     CoarseAverage avg(m_grids[lev],m_grids[lev-1],data[lev]->nComp(),m_refRatio[lev-1]);
  //     avg.averageToCoarse(*data[lev-1], *data[lev]);
  //   }
  
  //decide on the file name
  char file[32];
  sprintf(file,"jfnkopres.%06d.2d.hdf5",m_residualID);
  Real dt(1.0); Real time(m_residualID);
  pout() << "writing " << file << std::endl;

  WriteAMRHierarchyHDF5(file ,m_grids, data ,names, m_domains[0].domainBox(),
			m_dxs[0][0], dt, time, m_refRatio, data.size());


  for (int lev = 0; lev < a_u.size(); lev++)
    {
      delete data[lev];
    }
  m_residualID++;
}

void JFNKOp::preCond(Vector<LevelData<FArrayBox>*>& a_cor,
		     const Vector<LevelData<FArrayBox>*>& a_residual)
{
  //CH_assert(norm(a_residual,0) < HUGE_NORM);
  CH_TIME("JFNKOp::precond");
  m_mlOp.preCond( a_cor, a_residual);
  //CH_assert(norm(a_cor,0) < HUGE_NORM);
}

IceJFNKstate::~IceJFNKstate()
{
  delete m_bcPtr;
}
IceJFNKstate::IceJFNKstate
(const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatio,
 const Vector<ProblemDomain>& a_domains,
 const Vector<RealVect>& a_dxs,
 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<LevelData<FArrayBox>*>& a_u,
 const Vector<LevelData<FArrayBox>*>& a_C,
 const Vector<LevelData<FArrayBox>*>& a_C0,
 const int a_finestLevel,
 const ConstitutiveRelation& a_constRel,
 const BasalFrictionRelation& a_basalFrictionRel,
 IceThicknessIBC& a_bc,
 const Vector<LevelData<FArrayBox>*>& a_A,
 const Vector<LevelData<FluxBox>*>& a_faceA,
 const Real a_time,
 const Real a_vtopSafety,
 const int a_vtopRelaxMinIter,
 const Real a_vtopRelaxTol,
 const Real a_muMin ,
 const Real a_muMax)
  :m_u(a_u), m_C(a_C), m_C0(a_C0), m_grids(a_grids), m_refRatio(a_refRatio),
   m_domains(a_domains), m_dxs(a_dxs),m_finestLevel(a_finestLevel), 
   m_coordSys(a_coordSys), m_constRelPtr(&a_constRel), 
   m_basalFrictionRelPtr(&a_basalFrictionRel), m_bcPtr(a_bc.new_thicknessIBC()),
   m_A(a_A), m_faceA(a_faceA) , m_time(a_time), m_vtopSafety(a_vtopSafety),
   m_vtopRelaxMinIter(a_vtopRelaxMinIter),m_vtopRelaxTol(a_vtopRelaxTol),
   m_muMin(a_muMin),m_muMax(a_muMax)
{

  m_mu.resize(m_finestLevel + 1);
  m_lambda.resize(m_finestLevel + 1);
  m_alpha.resize(m_finestLevel + 1);
  m_muCoef.resize(m_finestLevel + 1);

  for (int lev =0; lev <= m_finestLevel; ++lev)
    {
      DisjointBoxLayout levelGrids = m_grids[lev];

      m_mu[lev] = RefCountedPtr<LevelData<FluxBox> >
	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Unit));

      m_muCoef[lev] = NULL;
      
      m_lambda[lev] = RefCountedPtr<LevelData<FluxBox> >
	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Unit));	
      
      m_alpha[lev] = RefCountedPtr<LevelData<FArrayBox> >
	(new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit));	
    }

   Real alpha = -1.0;
   Real beta = 1.0;
   // for the moment, at least, this only works for dx = dy:
   if (SpaceDim > 1) CH_assert(m_dxs[0][0] == m_dxs[0][1]);

   BCHolder velSolveBC = m_bcPtr->velocitySolveBC();
   m_opFactoryPtr = 
     RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
     (new ViscousTensorOpFactory(m_grids, m_mu, m_lambda, m_alpha, alpha, 
				 beta, m_refRatio , m_domains[0], m_dxs[0][0], 
				 velSolveBC, m_vtopSafety, m_vtopRelaxTol, 
				 m_vtopRelaxMinIter));


}


void IceJFNKstate::computeViscousTensorFace(const Vector<LevelData<FluxBox>*>& a_viscousTensor)
{
  if (m_opFactoryPtr == NULL)
    MayDay::Error("IceJFNKstate::computeViscousTensorFace missing m_opFactoryPtr");
  
  ViscousTensorOpFactory* vtopf = static_cast<ViscousTensorOpFactory*>(&(*m_opFactoryPtr));

  //compute cell centered grad(vel)
  Vector<LevelData<FArrayBox>* > gradVel(m_finestLevel + 1, NULL);
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      const DisjointBoxLayout& grids = m_grids[lev];
      gradVel[lev] = new LevelData<FArrayBox>(grids , SpaceDim*SpaceDim, IntVect::Unit);
      LevelData<FArrayBox>& vel = *m_u[lev];
      LevelData<FArrayBox>& grad = *gradVel[lev];
      
      if (lev > 0)
	{
	  const DisjointBoxLayout& crseGrids = m_grids[lev-1];
	  const ProblemDomain& pd = grids.physDomain();
	  TensorCFInterp cf(grids,&crseGrids,m_dxs[lev][0],m_refRatio[lev-1],SpaceDim,pd);
	  cf.coarseFineInterp(vel,grad,*m_u[lev-1]);
	}
      ViscousTensorOp* vtop = vtopf->AMRnewOp(m_domains[lev]);
      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	{
	  vtop->cellGrad(grad[dit],vel[dit],grids[dit]);
	}
      grad.exchange();
      delete vtop;
    }
  
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      ViscousTensorOp* vtop = vtopf->AMRnewOp(m_domains[lev]);
      const DisjointBoxLayout& grids = m_grids[lev];
      LevelData<FluxBox>& flux = *a_viscousTensor[lev];
      LevelData<FArrayBox>& vel = *m_u[lev];
      LevelData<FArrayBox>& grad = *gradVel[lev];
      
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& muFace  =  (*m_mu[lev])[dit][dir];
	      const FArrayBox& lambdaFace  = (*m_lambda[lev])[dit][dir];
	      Box faceBox =  surroundingNodes(grids[dit],dir);
	      Box domFaceBox = surroundingNodes(m_domains[lev].domainBox(), dir);
	      faceBox &= domFaceBox;
	      FArrayBox tmpFlux(faceBox,SpaceDim); //  vtop->getFlux tinkers with the FAB box
	      vtop->getFlux(tmpFlux, vel[dit], grad[dit] , muFace, lambdaFace, faceBox, dir, 1);
	      flux[dit][dir].setVal(0.0);
	      flux[dit][dir].copy(tmpFlux,faceBox);
	    }
	}
      delete vtop;
    }
  
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      if (gradVel[lev] != NULL)
	{
	  delete gradVel[lev];
	  gradVel[lev] = NULL; 
	}
    }
}



//store a_u and set the coeffients mu(u), lambda(u) = 2*mu(u) and alpha(u)
//in L[u] = div( mu * (grad(u) + grad(u)^T) + lambda * div(u)*I) - alpha*u
void IceJFNKstate::setState(const Vector<LevelData<FArrayBox>*>& a_u)
{
  m_u = a_u;
  
  CH_TIME("IceJFNKstate::setState");
  
  for (int lev=0; lev <= m_finestLevel ; lev++)
    {
      ProblemDomain levelDomain = m_domains[lev];
      const DisjointBoxLayout& levelGrids = m_grids[lev];
      LevelData<FArrayBox>& levelVel = *a_u[lev];
      LevelData<FArrayBox>* crseVelPtr = NULL;
      int nRefCrse = -1;
      if (lev > 0) 
        {
          crseVelPtr = a_u[lev-1];
          nRefCrse = m_refRatio[lev-1];
        }
      LevelData<FluxBox>& levelMu = *m_mu[lev];
     
      LevelSigmaCS& levelCoords = *m_coordSys[lev];
      const LevelData<FArrayBox>& levelC = *m_C[lev];
      LevelData<FArrayBox>& levelAlpha = *m_alpha[lev];

      // // first thing, if there is a finer level, average-down
      // // the current velocity field
      if (lev < m_finestLevel )
        {
          LevelData<FArrayBox>& finerVel = *a_u[lev+1];

          CoarseAverage averager(finerVel.getBoxes(),
                                 levelGrids,
                                 finerVel.nComp(),
                                 m_refRatio[lev]);

          averager.averageToCoarse(levelVel, finerVel);
        }

      // just in case, add an exchange here
      levelVel.exchange();
      
      // first set BC's on vel
      m_bcPtr->velocityGhostBC(levelVel,
                               levelCoords,
                               levelDomain, m_time);
      

      DataIterator dit = levelGrids.dataIterator();
      dit.reset();
#define TENSORCF
#ifndef TENSORCF 
      if (lev > 0) 
        {
          QuadCFInterp qcfi(levelGrids, &m_grids[lev-1],
                            m_dxs[lev][0], m_refRatio[lev-1], 
                            2, levelDomain);
          qcfi.coarseFineInterp(levelVel, *a_u[lev-1]);
        }

      //slc : qcfi.coarseFineInterp fills the edges of lev > 0 cells
      //but not the corners. We need them filled to compute the
      //rate-of-strain invariant, so here is a bodge for now
      if (SpaceDim == 2)
      	{
      	  for (dit.begin(); dit.ok(); ++dit)
      	    {
      	      Box sbox = levelVel[dit].box();
      	      sbox.grow(-1);
      	      FORT_EXTRAPCORNER2D(CHF_FRA(levelVel[dit]),
      				  CHF_BOX(sbox));
      	    }
	  
      	}

#endif
      // actually need to use a cornerCopier, too...
      CornerCopier cornerCopier(levelGrids, levelGrids, 
                                levelDomain,levelVel.ghostVect(),
                                true);
      levelVel.exchange(cornerCopier);
           
      (*m_constRelPtr).computeFaceMu(levelMu,
                                     levelVel,
                                     crseVelPtr,
                                     nRefCrse,
                                     *m_faceA[lev],
                                     levelCoords,
				     levelDomain,
                                     IntVect::Zero);

      // now limit and  multiply by ice thickness H
      const LevelData<FluxBox>& faceH = levelCoords.getFaceH();
      LevelData<FluxBox>* muCoefPtr = m_muCoef[lev];
      for (dit.begin(); dit.ok(); ++dit)
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
	      FArrayBox& thisMu = levelMu[dit][dir];
	      const Box& box = thisMu.box();
	       
	      FORT_MAXFAB1(CHF_FRA(thisMu),
	       		  CHF_CONST_REAL(m_muMin),
	       		  CHF_BOX(box));

	      if (muCoefPtr!=NULL)
		{
		  levelMu[dit][dir].mult( (*muCoefPtr)[dit][dir] );
		}

	      Box facebox = levelGrids[dit].surroundingNodes(dir);
	      CH_assert(levelMu[dit][dir].min(facebox) >= 0.0);

              levelMu[dit][dir].mult(faceH[dit][dir],
                                    levelMu[dit][dir].box(),0,0,1);
	      //levelMu[dit][dir]*=1000.0;
	      FORT_MINFAB1(CHF_FRA(thisMu),
	       		  CHF_CONST_REAL(m_muMax),
	       		  CHF_BOX(box));
	      
	    }
        
	  // also update alpha
          const Box& gridBox = levelGrids[dit];
	  
	  m_basalFrictionRelPtr->computeAlpha
	    (levelAlpha[dit], levelVel[dit], levelCoords.getThicknessOverFlotation()[dit] 
	     , levelC[dit] , 	levelCoords.getFloatingMask()[dit], gridBox);
	  levelAlpha[dit] += (*m_C0[lev])[dit];
	  CH_assert(levelAlpha[dit].min() >= 0.0);

	  
	} // end loop over grids		
    }

  //coarse average mu and alpha
  for (int lev = m_finestLevel; lev > 0 ; lev--)
  {
      CoarseAverage averageOp(m_grids[lev],1,m_refRatio[lev-1]);
      averageOp.averageToCoarse(*m_alpha[lev-1], *m_alpha[lev]);

      CoarseAverageFace averageOpFace(m_grids[lev],1,m_refRatio[lev-1]);
      averageOpFace.averageToCoarse(*m_mu[lev-1], *m_mu[lev] );
  }
  // lambda = 2*mu
  for (int lev=0; lev <= m_finestLevel ; lev++)
    {
      const DisjointBoxLayout& levelGrids = m_grids[lev];
      LevelData<FluxBox>& levelLambda = *m_lambda[lev];
      LevelData<FluxBox>& levelMu = *m_mu[lev];
      LevelData<FArrayBox>& levelAlpha = *m_alpha[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{

#if CH_SPACEDIM==2
	  {
	    Real mu0 = 1.0;
	    Real C0 = 1.0;
	    
	    FORT_ENFORCEWELLPOSEDCELL
	      (CHF_FRA1(levelAlpha[dit],0),
	       CHF_FRA1(levelMu[dit][0],0),
	       CHF_FRA1(levelMu[dit][1],0),
	       CHF_CONST_REAL(mu0),
	       CHF_CONST_REAL(C0),
	       CHF_BOX(levelGrids[dit]));

	  }
#endif


	  FluxBox& lambda = levelLambda[dit];
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      lambda[dir].copy(levelMu[dit][dir]);
	      lambda[dir] *= 2.0;
	    }
	} // end loop over grids
    } // end loop over levels

}





void JFNKSolver::setDefaultParameters()
{
//sensible defaults
  m_residualOnly = false;
  //m_linearSolverType=Relax; 
  m_linearSolverType=BiCGStab;
  m_maxIter = 15;
  m_absTol = 1.0e-10;
  m_relTol = 1.0e-10;
  m_BiCGStabRelTol = 1.0e-3;
  m_maxBiCGStabIter = 10;
  m_CGRelTol = 1.0e-3;
  m_maxCGIter = 10;
  m_RelaxRelTol = 1.0e-3;
  m_maxRelaxIter = 10;
  m_RelaxHang = 0.25;
  m_GMRESRelTol = 1.0e-3;
  m_maxGMRESIter = 10;
  m_normType = 0;
  m_verbosity = 5;
  m_vtopSafety = 0.5;
  m_vtopRelaxMinIter = 8;
  m_vtopRelaxTol = 1.0e-2;
  m_numMGSmooth = 8;
  m_numMGIter = 1;
  m_h = 0.025;
  m_hAdaptive=false;
  m_err = 1.0e-6;
  m_umin = 1.0e-6;
  m_switchRate = 1.5;
  m_minPicardIter = 0;
  m_uMaxAbs = 1.23456789e+300;
  m_muMax = 1.23456789e+300;
  m_muMin = 0.0;
  m_writeResiduals = false;
  m_minStepFactor = 1.0;
  // these ones don't need to be stored (at least for now), but should be set
  int mgAverageType  = CoarseAverageFace::arithmetic;
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  ViscousTensorOp::s_prolongType = mgProlongType;

}

void JFNKSolver::define(const ProblemDomain& a_coarseDomain,
			ConstitutiveRelation* a_constRel,
			BasalFrictionRelation* a_basalFrictionRel,
			const Vector<DisjointBoxLayout>& a_grids,
			const Vector<int>& a_refRatios,
			const RealVect& a_dxCrse,
			IceThicknessIBC* a_bc,
			int a_numLevels)
{

  // set parameters based on parmParse
  ParmParse pp("JFNKSolver");
  pp.query("maxIter",m_maxIter);
  pp.query("absTol", m_absTol);
  pp.query("relTol",m_relTol);
  pp.query("normType",m_normType);
  pp.query("verbosity",m_verbosity);
  pp.query("vtopSafety",m_vtopSafety);
  pp.query("vtopRelaxTol",m_vtopRelaxTol);
  pp.query("vtopRelaxMinIter",m_vtopRelaxMinIter);
  pp.query("numMGSmooth",m_numMGSmooth);
  pp.query("numMGIter",m_numMGIter);
  pp.query("h",m_h);
  pp.query("hAdaptive",m_hAdaptive);
  pp.query("err",m_err);
  pp.query("umin",m_umin);
  pp.query("switchRate",m_switchRate);
  pp.query("minPicardIterations", m_minPicardIter);
  pp.query("muMax", m_muMax);
  pp.query("muMin", m_muMin);
  pp.query("uMaxAbs", m_uMaxAbs);
  pp.query("writeResiduals", m_writeResiduals);
  pp.query("minStepFactor", m_minStepFactor);
  if (pp.contains("solverType") )
    {
      int solverIntType = m_linearSolverType;
      pp.query("solverType", solverIntType);
      if (solverIntType == BiCGStab)
        {
          m_linearSolverType = BiCGStab;
	  pp.query("BiCGStabRelTol",m_BiCGStabRelTol);
	  pp.query("maxBiCGStabIter",m_maxBiCGStabIter);
        }
      else if (solverIntType == Relax)
        {
          m_linearSolverType = Relax;
	  pp.query("RelaxTol", m_RelaxRelTol);
	  pp.query("RelaxRelTol", m_RelaxRelTol);
	  pp.query("maxRelaxIter", m_maxRelaxIter);
	  pp.query("RelaxHang", m_RelaxHang);
	  m_numMGIter = 1; // m_numMGIter > 1 doesn't help
        }
      else if (solverIntType == GMRES)
        {
          m_linearSolverType = GMRES;
	  MayDay::Error("JFNKSolver -- GMRES is on the blink");
        }
      else if (solverIntType == CG)
        {
	  m_linearSolverType = CG;
	  pp.query("CGRelTol",m_CGRelTol);
	  pp.query("maxCGIter",m_maxCGIter);
	}
      else if (solverIntType == petsc)
        {
#ifdef CH_USE_PETSC
          m_linearSolverType = petsc;
#else
          // default back to relax if petsc isn't compiled in
          // (this is just to simplify comparisons btwn petsc and MG)
          m_linearSolverType = Relax;
	  pp.query("RelaxTol", m_RelaxRelTol);
	  pp.query("RelaxRelTol", m_RelaxRelTol);
	  pp.query("maxRelaxIter", m_maxRelaxIter);
	  pp.query("RelaxHang", m_RelaxHang);
	  m_numMGIter = 1; // m_numMGIter > 1 doesn't help
#endif          
        }      
      else 
        {
          MayDay::Error("JFNKSolver -- bad linear solver type");
        }
    }
  
  int bs_type = JFNKOp::m_bottom_solver_type;
  pp.query("bottom_solver_type", bs_type);
  JFNKOp::m_bottom_solver_type = bs_type;

  int mg_depth = JFNKOp::m_MG_solver_depth;
  pp.query("mg_solver_depth", mg_depth);
  JFNKOp::m_MG_solver_depth = mg_depth;

  int mgAverageType  = CoarseAverageFace::arithmetic;
  pp.query("mgCoefficientAverageType", mgAverageType);
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  pp.query("mgProlongType", mgProlongType);

  ViscousTensorOp::s_prolongType = mgProlongType;
   
 

  m_constRelPtr = a_constRel;
  m_basalFrictionRelPtr = a_basalFrictionRel;
  m_bcPtr = a_bc;
  
  m_grids.resize(a_numLevels); 
  m_grids[0] = a_grids[0];

  m_refRatios.resize(a_numLevels);
  m_refRatios = a_refRatios;

  m_dxs.resize(a_numLevels);
  m_dxs[0] = a_dxCrse;

  m_domains.resize(a_numLevels);
  m_domains[0] = a_coarseDomain;
  
  
  for (int lev = 1; lev < a_numLevels; ++lev)
    {
      m_dxs[lev] = m_dxs[lev-1] / m_refRatios[lev-1];
      m_domains[lev] = m_domains[lev-1];
      m_domains[lev].refine(m_refRatios[lev-1]);
      m_grids[lev] = a_grids[lev];
      
    }

}

//IceVelocitySolver full solve
inline
int JFNKSolver::solve(Vector<LevelData<FArrayBox>* >& a_u,
		      Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		      const Real a_convergenceMetric,
		      const Vector<LevelData<FArrayBox>* >& a_rhs,
		      const Vector<LevelData<FArrayBox>* >& a_C,
		      const Vector<LevelData<FArrayBox>* >& a_C0,
		      const Vector<LevelData<FArrayBox>* >& a_A,
		      const Vector<LevelData<FluxBox>* >& a_muCoef,
		      Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		      Real a_time,  int a_lbase, int a_maxLevel)
{

  //Vector<LevelData<FluxBox>* > muCoef(a_maxLevel,NULL);

  return solve(a_u, a_initialResidualNorm, a_finalResidualNorm,
	       a_convergenceMetric, false, a_rhs, a_C, a_C0, a_A, 
	       a_muCoef, a_coordSys,
	       a_time, a_lbase,  a_maxLevel);

}			    

//general solve, linear or non-linear
int JFNKSolver::solve(Vector<LevelData<FArrayBox>* >& a_u,
		      Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		      const Real a_convergenceMetric,
		      const bool a_linear,
		      const Vector<LevelData<FArrayBox>* >& a_rhs,
		      const Vector<LevelData<FArrayBox>* >& a_C,
		      const Vector<LevelData<FArrayBox>* >& a_C0,
		      const Vector<LevelData<FArrayBox>* >& a_A,
		      const Vector<LevelData<FluxBox>* >& a_muCoef,
		      Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		      Real a_time,  int a_lbase, int a_maxLevel)
{
  CH_TIME("JFNKSolver::solve");

  int returnCode = 0;
  
  Vector<LevelData<FArrayBox>* > localU(a_maxLevel+1);
  Vector<LevelData<FArrayBox>* > localRhs(a_maxLevel+1);
  Vector<LevelData<FArrayBox>* > localC(a_maxLevel+1);
  Vector<LevelData<FArrayBox>* > localC0(a_maxLevel+1);
 
 for (int lev = 0; lev < a_maxLevel + 1; ++lev)
   {
     localU[lev] = a_u[lev];
     localRhs[lev] = a_rhs[lev];
     localC[lev] = a_C[lev]; 
     localC0[lev] = a_C0[lev]; 
   }
 
 
 //cell face A
 Vector<LevelData<FluxBox>* > faceA(a_maxLevel+1, NULL);

 for (int lev = 0; lev < a_maxLevel + 1; ++lev)
   {
     faceA[lev] = new LevelData<FluxBox>
       (m_grids[lev], a_A[lev]->nComp(), IntVect::Zero);
     CellToEdge(*a_A[lev] , *faceA[lev]);
   }
 

 //define state managers
 IceJFNKstate current
   (m_grids , m_refRatios, m_domains,m_dxs , a_coordSys, localU ,
    localC, localC0, a_maxLevel, *m_constRelPtr, *m_basalFrictionRelPtr, *m_bcPtr, 
    a_A, faceA, a_time, m_vtopSafety, m_vtopRelaxMinIter, m_vtopRelaxTol, 
    m_muMin, m_muMax);
 
 IceJFNKstate perturbed
   (m_grids , m_refRatios, m_domains, m_dxs , a_coordSys, localU ,
    localC, localC0, a_maxLevel, *m_constRelPtr,*m_basalFrictionRelPtr, *m_bcPtr, 
    a_A, faceA, a_time, m_vtopSafety, m_vtopRelaxMinIter, m_vtopRelaxTol,
    m_muMin, m_muMax);
 
 current.setFaceViscCoef(a_muCoef);perturbed.setFaceViscCoef(a_muCoef);

 //define a JFNKOp
 CH_assert(a_lbase == 0); //\todo : support lBase != 0
 JFNKOp jfnkOp
   (&current, &perturbed, localU, m_h, m_err, m_umin, m_hAdaptive,  m_grids, m_refRatios, 
    m_domains, m_dxs, a_lbase, m_numMGSmooth, m_numMGIter, PICARD_SOLVER_MODE);
 
  Vector<LevelData<FArrayBox>*> residual;
  jfnkOp.create(residual,localU);
  
  Vector<LevelData<FArrayBox>*> du;
  jfnkOp.create(du,localU);
  
  

  

  int iter = 0;
  Real oldResNorm;
  Real& resNorm = a_finalResidualNorm;
  jfnkOp.m_writeResiduals = m_writeResiduals;
  jfnkOp.outerResidual(residual,localU,localRhs);
  resNorm = a_initialResidualNorm = jfnkOp.norm(residual, m_normType);
  if (m_verbosity > 0)
    {
      pout() << "JFNK initial residual norm = " << resNorm << std::endl;
    }
  Real convergenceMetric = (a_convergenceMetric < 0.0)?a_initialResidualNorm:a_convergenceMetric;
  if (m_verbosity > 0)
    {
      pout() << "JFNK convergence metric = " <<  convergenceMetric << std::endl;
    }
  if (m_residualOnly)
    return 0;

  bool done =  (resNorm  < m_absTol) |  (resNorm < a_convergenceMetric * m_relTol);
  SolverMode mode = PICARD_SOLVER_MODE;
  if (m_minPicardIter < 0)
    {
      mode = JFNK_SOLVER_MODE;
    }


  while (!done && iter < m_maxIter)
    {
      oldResNorm = resNorm;
      
      // if (m_writeIterations)
      // 	writeIteration(iter,localU,localRhs,residual);

      if (iter == 0 || !a_linear)
	current.setState(a_u);  
      JFNKOp newOp
	(&current, &perturbed, localU, m_h, m_err, m_umin, m_hAdaptive, m_grids, 
	 m_refRatios, m_domains, m_dxs, a_lbase, m_numMGSmooth, m_numMGIter, mode);

      if (a_linear)
	 newOp.setToZero(localU);

      newOp.m_writeResiduals = m_writeResiduals;

      LinearSolver<Vector<LevelData<FArrayBox>* > >* krylovSolver = NULL;
      if (m_linearSolverType == BiCGStab)
	{
	  BiCGStabSolver<Vector<LevelData<FArrayBox>* > >* biCGStabSolver 
	    = new BiCGStabSolver<Vector<LevelData<FArrayBox>* > >;
	  
	  biCGStabSolver->define(&newOp , false);
	  biCGStabSolver->m_verbosity = m_verbosity - 1;
	  biCGStabSolver->m_reps = m_BiCGStabRelTol;
	  biCGStabSolver->m_imax = m_maxBiCGStabIter; 
	  biCGStabSolver->m_normType = m_normType;
	  //JFNK mode will only work well if 
	  //the linear system is solved quickly.
	  //By being intolerant here, we revert to 
          //the cheaper Picard mode sooner rather than later
	  if (mode == JFNK_SOLVER_MODE)
	    {
	      biCGStabSolver->m_numRestarts = 0;
	      biCGStabSolver->m_hang = m_RelaxHang;
	    }
	  krylovSolver = biCGStabSolver;
	}
      else if (m_linearSolverType == CG)
	{
	  CGSolver<Vector<LevelData<FArrayBox>* > >* cGSolver 
	    = new CGSolver<Vector<LevelData<FArrayBox>* > >;
	  
	  cGSolver->define(&newOp , false);
	  cGSolver->m_verbosity = m_verbosity - 1;
	  cGSolver->m_eps = m_CGRelTol;
	  cGSolver->m_imax = m_maxCGIter;
	  cGSolver->m_normType = m_normType;

	  cGSolver->m_numRestarts = 0;
	  cGSolver->m_hang = 1.0;
	  //JFNK mode will only work well if 
	  //the linear system is solved quickly.
	  //By being intolerant here, we revert to 
          //the cheaper Picard mode sooner rather than later
	  if (mode == JFNK_SOLVER_MODE)
	    {
	      cGSolver->m_numRestarts = 0;
	      cGSolver->m_hang = 1.0;
	    }
	  krylovSolver = cGSolver;
	}
      else if (m_linearSolverType == GMRES)
	{
	  GMRESSolver<Vector<LevelData<FArrayBox>* > >* gmresSolver 
	    = new GMRESSolver<Vector<LevelData<FArrayBox>* > >();
	  
	  gmresSolver->define(&newOp , false);
	  gmresSolver->m_verbosity = m_verbosity - 1;
	  gmresSolver->m_reps = m_GMRESRelTol;
	  gmresSolver->m_imax = m_maxGMRESIter;
	  gmresSolver->m_normType = m_normType;
	  //JFNK mode will only work well if 
	  //the linear system is solved quickly.
	  //By being intolerant here, we revert to 
          //the cheaper Picard mode sooner rather than later
	  //krylovSolver.m_hang = 1.0;
	  if (mode == JFNK_SOLVER_MODE)
	    {
	      //gmresSolver->m_hang = 1.0;
	    }
	  krylovSolver = gmresSolver;
	}
      else if (m_linearSolverType == Relax)
	{
	  
	  RelaxSolver<Vector<LevelData<FArrayBox>* > >* relaxSolver
	    = new RelaxSolver<Vector<LevelData<FArrayBox>* > >();
	  
	  relaxSolver->define(&newOp,false);
	  relaxSolver->m_verbosity = m_verbosity - 1;
	  relaxSolver->m_normType = m_normType;
	  relaxSolver->m_eps = m_RelaxRelTol;
	  relaxSolver->m_imax = m_maxRelaxIter;
	  relaxSolver->m_hang = m_RelaxHang;
	  krylovSolver = relaxSolver;
	}
#ifdef CH_USE_PETSC
      else if (m_linearSolverType == petsc)
        {
          // handle petsc solver a bit differently
          // since setup is more expensive _and_ it's easy 
          // to swap the coefficients, keep a pre-defined 
          // solver around and just reset the coefficients
          Real opAlpha, opBeta;
          opAlpha = -1.0;
          opBeta = 1.0;
          if (m_petscSolver == NULL)
            {
              m_petscSolver = new PetscAMRSolver;
	      m_petscSolver->m_petscCompMat.setDiri(true);
            }
	  RefCountedPtr<ConstDiriBC> bcfunc = RefCountedPtr<ConstDiriBC>(new ConstDiriBC(1,m_petscSolver->m_petscCompMat.getGhostVect()));
	  BCHolder bc(bcfunc);
	  m_petscSolver->m_petscCompMat.define(m_domains[0],m_grids,m_refRatios,bc,m_dxs[0]); // generic AMR setup
	  m_petscSolver->m_petscCompMat.defineCoefs(opAlpha,opBeta,current.m_mu,current.m_lambda,current.m_alpha);
        }
#endif
      else 
	{
	  MayDay::Error("JFNKSolver::solve unknown linear solver type");
	}

      
      if (mode == JFNK_SOLVER_MODE && !a_linear)
	{

	  newOp.setToZero(du);
          if (m_linearSolverType != petsc)
            {
              krylovSolver->solve(du, residual);
            }
#ifdef CH_USE_PETSC
          else if (m_linearSolverType == petsc)
            {
	      m_petscSolver->solve_mfree(du,residual,&newOp);
            }
#endif
          else 
            {
              MayDay::Error("linearSolverType is petsc, but code compiled w/o PETSC=TRUE");
            }

	  //take the Newton step du, and if the residual is not reduced try a sequence
	  //of smaller steps w*du, halving w until either the residual is reduced or
	  //w < m_minStepFactor
	  Real w = 1.0; 
	  newOp.incr(localU,du,w);
	  do 
            {
	      current.setState(localU);      
	      JFNKOp testOp (&current, &perturbed, localU, m_h, m_err, m_umin, m_hAdaptive, m_grids, 
			     m_refRatios, m_domains, m_dxs, a_lbase, 
			     m_numMGSmooth, m_numMGIter, mode);
	      testOp.outerResidual(residual,localU,localRhs);
	      resNorm = testOp.norm(residual, m_normType);
	      if (resNorm >= oldResNorm)
		{
		  w *= 0.5;
		  if (w >= m_minStepFactor)
		    {
		      //step halfway back to the last U
		      newOp.incr(localU,du,-w);
		      if (m_verbosity > 0)
			pout() << "JFNK iteration " << iter  << " halving step length, w =   " << w << std::endl; 	
		    }
		  else
		    {
		      //give up, return to the last U
		      newOp.incr(localU,du,-2.0*w);
		    } 
		}
	    }
          while ( resNorm >= oldResNorm && w >= m_minStepFactor);
	  
	  if (resNorm > oldResNorm)
	    {
	      if (m_verbosity > 0){
		pout() << "JFNK iteration " << iter  << 
		  " did not reduce residual" << std::endl;
                if (m_verbosity > 1)
                  {
                    pout() << " ---- old residual = " << oldResNorm 
                           << " -- new residual = " << resNorm
                           << std::endl;
                  }
              }
	      
	      //don't take the step and 
	      //switch back to Picard mode
	      //newOp.incr(localU,du,-w); 
	      //imposeMaxAbs(localU,m_uMaxAbs); 
	      mode = PICARD_SOLVER_MODE;
	      resNorm = oldResNorm;
	    }
	
	  
	  if (m_verbosity > 0)
	    {
	      pout() << " JFNK iteration " << iter 
		     << " residual norm = " << resNorm 
		     << " rate = " <<  oldResNorm / resNorm 
		     << std::endl;
	      
	    }
	  
	} // end if (mode == JFNK_SOLVER_MODE)
      else if ( (mode == PICARD_SOLVER_MODE) | a_linear)
	{
          if (m_linearSolverType != petsc)
            {
              krylovSolver->solve(localU, localRhs);
            }
#ifdef CH_USE_PETSC
          else if (m_linearSolverType == petsc)          
            {
              m_petscSolver->solve(localU,localRhs);
            }
#endif
          else
            {
              MayDay::Error("linearSolverType is petsc, but code compiled w/o PETSC=TRUE");
            }
          
          
	  if (a_linear)
	    {
	      //if we are solving  a linear equation ,there is no
              //need to update the state, and no point in switching to 
	      //JFNK mode
	     
	      newOp.outerResidual(residual,localU,localRhs);
	      resNorm = newOp.norm(residual, m_normType);
	    }
	  else
	    {
	      
	     
	     
	      imposeMaxAbs(localU,m_uMaxAbs);
	      current.setState(localU); 
	      newOp.outerResidual(residual,localU,localRhs);
	      resNorm = newOp.norm(residual, m_normType);
	      Real rate = oldResNorm / resNorm;
	      if (m_verbosity > 0)
		{
		  pout() << " Picard iteration " << iter 
			 << " residual norm = " << resNorm 
			 << " rate = " <<  rate
			 << std::endl;
		}
	  
	      if (iter >= m_minPicardIter && rate > 1.0 && rate < m_switchRate){
		//since the iterations are progressing slowly
		//but positively, try JFNK mode
		mode = JFNK_SOLVER_MODE;
	      }
	    }

	}// end if (mode == PICARD_SOLVER_MODE)
	
      if (m_linearSolverType != petsc)
        {
          delete krylovSolver;
        }
      done = a_linear || resNorm < m_absTol || resNorm < convergenceMetric * m_relTol;
      
      ++iter;
    }
  // if (m_writeIterations)
  //   writeIteration(iter,localU,localRhs,residual);
  // check to see if we actually solved the problem
  pout() << "JFNK final residual norm = " << resNorm << std::endl;

  

  if (done)
    {
      returnCode = 0;
      if (m_verbosity > 0)
        {
          pout() << "JFNKSolver converged -- final norm(resid) = "
                 << resNorm << " after " << iter << " iterations"
                 << endl;
        }
    }
  else 
    {
      returnCode = 1;
      if (m_verbosity > 0)
        {
          pout() << "JFNKSolver NOT CONVERGED -- final norm(resid) = "
                 << resNorm << " after " << iter << " iterations"
                 << endl;          
        }
    }
  
  // clean up storage
  for (int lev=0; lev<residual.size(); lev++)
    {

 
      if (residual[lev] != NULL)
	{
	  delete residual[lev];
	  residual[lev] = NULL;
	}
      
      if (du[lev] != NULL)
	{
	  delete du[lev];
	  du[lev] = NULL;
	}

      if (faceA[lev] != NULL)
	{
	  delete faceA[lev];
	  faceA[lev] = NULL;
	}
    }
  
  return returnCode;

} 

#include "NamespaceFooter.H"


/*
  PETSc Licensing Notification
  
  Permission to use, reproduce, prepare derivative works, and to redistribute to others this software, derivatives of this software, and future versions of this software as well as its documentation is hereby granted, provided that this notice is retained thereon and on all copies or modifications. This permission is perpetual, world-wide, and provided on a royalty-free basis. UChicago Argonne, LLC and all other contributors make no representations as to the suitability and operability of this software for any purpose. It is provided "as is" without express or implied warranty.
  Principal Software authors
  
  
  
  Mathematics and Computer Science Division
  Argonne National Laboratory,
  Argonne IL 60439
  Any questions or comments on the software may be directed to petsc-maint@mcs.anl.gov.
  
  Portions of this software are copyright by UChicago Argonne, LLC. Argonne National Laboratory with facilities in the state of Illinois, is owned by The United States Government, and operated by UChicago Argonne, LLC under provision of a contract with the Department of Energy.
  DISCLAIMER
  
  PORTIONS OF THIS SOFTWARE WERE PREPARED AS AN ACCOUNT OF WORK SPONSORED BY AN AGENCY OF THE UNITED STATES GOVERNMENT. NEITHER THE UNITED STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, NOR ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. REFERENCE HEREIN TO ANY SPECIFIC COMMERCIAL PRODUCT, PROCESS, OR SERVICE BY TRADE NAME, TRADEMARK, MANUFACTURER, OR OTHERWISE, DOES NOT NECESSARILY CONSTITUTE OR IMPLY ITS ENDORSEMENT, RECOMMENDATION, OR FAVORING BY THE UNITED STATES GOVERNMENT OR ANY AGENCY THEREOF. THE VIEW AND OPINIONS OF AUTHORS EXPRESSED HEREIN DO NOT NECESSARILY STATE OR REFLECT THOSE OF THE UNITED STATES GOVERNMENT OR ANY AGENCY THEREOF. */
