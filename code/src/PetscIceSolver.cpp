#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "PetscIceSolver.H"
#include "ViscousTensorOp.H"
#include "ParmParse.H"
#include "CornerCopier.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"
#include "BisiclesF_F.H"
#include "TensorCFInterp.H"
#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::PetscIceSolver()
////////////////////////////////////////////////////////////////////////
PetscIceSolver::PetscIceSolver()
{
  // default constructor leaves things in an undefined state
  setDefaultValues();
  m_bc = NULL;
  m_constRelPtr = NULL;
  m_basalFrictionRelPtr = NULL;
  m_opFactoryPtr = NULL;

  m_isOpDefined = false;
  m_vtopSafety = VTOP_DEFAULT_SAFETY;

  //m_refRatio = 0;
  //m_OpPtr = NULL;

  m_max_its = 20;
  m_rtol = 1.e-6;
  m_atol = 1.e-30;
}
////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::~PetscIceSolver() 
////////////////////////////////////////////////////////////////////////
PetscIceSolver::~PetscIceSolver() 
{
  if(m_opFactoryPtr != NULL)
    {
      delete m_opFactoryPtr;
      m_opFactoryPtr = NULL;
    }
  
  if (m_bc != NULL)
    {
      delete m_bc;
      m_bc = NULL;
    }
}

#ifdef CH_USE_PETSC
/* ------------------------------------------------------------------- */
/* 
   FormFunction

   Input Parameters:
   user - user-defined application context
.  x - vector
.  ctx - user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  f - vector
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction( SNES snes, Vec x, Vec f, void *ctx )
{
  CH_TIME("PetscIceSolver::FormFunction");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;
  PetscInt *pilev;

  PetscFunctionBegin;

  ierr = SNESGetApplicationContext(snes,(void**)&pilev); CHKERRQ(ierr);
pout() << "FormFunction ilev= " << *pilev << endl;
  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  ierr = solver->putPetscInChombo( *tthis->m_tphi2, x );     CHKERRQ(ierr);
  
  tthis->updateCoefs( *tthis->m_tphi2, (int)*pilev ); // needed because called before FormJacobian

  tthis->m_op[*pilev]->applyOp( *tthis->m_tphi, *tthis->m_tphi2 ); 

  ierr = solver->putChomboInPetsc( f, *tthis->m_tphi );  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian - this does not depend on PetscIceSolver so it could be moved into Mat solver.

   This can go into the PetscSolver class -- it is not BICICLES specific!!!

   Input Parameters:
.  snes - the SNES context
.  dummy - input vector
.  ctx - optional user-defined context, as set by SNESSetJacobian()

   Output Parameters:
.  jac - Jacobian matrix
.  prejac - different preconditioning matrix
.  flag - flag indicating matrix structure
*/
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian( SNES snes,Vec x,Mat *jac,Mat *prejac,MatStructure *flag, void *ctx )
{
  CH_TIME("PetscIceSolver::FormJacobian");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;
  PetscInt *pilev; // not used

  PetscFunctionBegin;

  ierr = SNESGetApplicationContext(snes,(void**)&pilev); CHKERRQ(ierr);
pout() << "\t\tFormJacobian ilev= " << *pilev << endl;
  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  ierr = solver->putPetscInChombo( *tthis->m_tphi2, x );     CHKERRQ(ierr);

  //tthis->updateCoefs(*tthis->m_tphi2, (int)*pilev); // needed because called after FormJacobian
  ierr = solver->formMatrix( *prejac, *tthis->m_tphi2 ); CHKERRQ(ierr);

  // 
  ierr = MatAssemblyBegin(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if(prejac!=jac)
    {
      ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

  *flag = SAME_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}
#endif // petsc

////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::define
////////////////////////////////////////////////////////////////////////
void 
PetscIceSolver::define(const ProblemDomain& a_coarseDomain,
		       ConstitutiveRelation* a_constRelPtr,
		       BasalFrictionRelation* a_basalFrictionRelPtr,
		       const Vector<DisjointBoxLayout>& a_vectGrids,
		       const Vector<int>& a_vectRefRatio,
		       const RealVect& a_dxCrse,
		       IceThicknessIBC* a_bc,
		       int a_numLevels)
{
  Real opAlpha, opBeta;

  getOperatorScaleFactors( opAlpha, opBeta );

  m_constRelPtr = a_constRelPtr;
  m_basalFrictionRelPtr = a_basalFrictionRelPtr;
  m_bc = a_bc->new_thicknessIBC();

  m_op.resize(a_numLevels);
  m_grid.resize(a_numLevels);
  m_domain.resize(a_numLevels);
  m_refRatio.resize(a_numLevels-1);
  m_Mu.resize(a_numLevels);
  m_Lambda.resize(a_numLevels);
  m_Beta.resize(a_numLevels);
  m_Beta0.resize(a_numLevels);
  m_C.resize(a_numLevels);
  m_petscSolver.resize(a_numLevels);
  
  m_domain[0] = a_coarseDomain;
  for(int ilev=0;ilev<a_numLevels;ilev++)
    {
      m_grid[ilev] = a_vectGrids[ilev];
      if( ilev < a_numLevels-1 )
	{
	  m_refRatio[ilev] = a_vectRefRatio[ilev];
	}

      m_Mu[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grid[ilev], 
									     1, 
									     IntVect::Zero) );
      
      m_Lambda[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grid[ilev],
										 1, 
										 IntVect::Zero) );
      
      // beta only has one component...
      // m_Beta[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grid[ilev],
      // 										   1, 
      // 										   IntVect::Zero));
      // m_Beta0[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grid[ilev],
      // 										    1, 
      // 										    IntVect::Zero));
      
      // C only has one component...
      m_C[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grid[ilev], 
										1,
										IntVect::Zero));
      // make domain
      if(ilev>0)
	{
	  m_domain[ilev] = m_domain[ilev-1];
	  m_domain[ilev].refine(m_refRatio[ilev-1]);
pout() << "PetscIceSolver::define WOW!!! new fine grid refrat:" << m_refRatio[ilev-1] << endl;
	}
    }
  
  // create op factory
  defineOpFactory( a_dxCrse, a_coarseDomain, a_numLevels );
  
  for(int ilev=0;ilev<a_numLevels;ilev++)
    {
      // this copies the unset data above, just needed here for dx &crdx.
      m_op[ilev] = RefCountedPtr<ViscousTensorOp>(m_opFactoryPtr->MGnewOp(m_domain[ilev],0)); 
#ifdef CH_USE_PETSC
      m_petscSolver[ilev] = RefCountedPtr<PetscSolverViscousTensor<LevelData<FArrayBox> > >(new PetscSolverViscousTensor<LevelData<FArrayBox> >);
      m_petscSolver[ilev]->define( &(*m_op[ilev]), false ); // dx & crdx
      m_petscSolver[ilev]->setVTParams( opAlpha, opBeta, 
					&(*m_C[ilev]), 
					&(*m_Mu[ilev]), 
					&(*m_Lambda[ilev]) );
      
      m_petscSolver[ilev]->setFunctionAndJacobian( FormFunction, FormJacobian );
#endif
    }
}

void
PetscIceSolver::computeAMRLevelsResidual( RefCountedPtr<LevelData<FArrayBox> > a_resid, 
					  int a_ilev, int a_lbase, int a_maxLevel, 
					  const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
					  const Vector<LevelData<FArrayBox>* >& a_rhs )
{
  if (a_ilev==a_lbase) 
    {
      if(a_ilev==a_maxLevel-1) // one level
	{
	  m_op[a_ilev]->residual(*a_resid, 
				 *a_horizontalVel[a_ilev], 
				 *a_rhs[a_ilev], true ); 
	}
      else // coarse grid
	{
	  m_op[a_ilev]->AMRResidualNC( *a_resid, 
				       *a_horizontalVel[a_ilev+1], 
				       *a_horizontalVel[a_ilev], 
				       *a_rhs[a_ilev], true, 
				       (a_ilev == a_maxLevel-1) ? 0 : &(*m_op[a_ilev+1]));
	}
    }
  else if (a_ilev == a_maxLevel-1) // fine grid
    {
      m_op[a_ilev]->AMRResidualNF( *a_resid, 
				   *a_horizontalVel[a_ilev],
				   *a_horizontalVel[a_ilev-1],
				   *a_rhs[a_ilev], true );
    }
  else
    {
      m_op[a_ilev]->AMRResidual( *a_resid, 
				 *a_horizontalVel[a_ilev+1], 
				 *a_horizontalVel[a_ilev], 
				 *a_horizontalVel[a_ilev-1], 
				 *a_rhs[a_ilev], true, 
				 &(*m_op[a_ilev+1]));
    }

  // zero covered
  if(a_ilev != a_maxLevel-1) // ! fine grid
    {
      const int nc = a_horizontalVel[a_ilev]->nComp();
      const DisjointBoxLayout& finedbl = a_horizontalVel[a_ilev+1]->disjointBoxLayout();
      DisjointBoxLayout dblCoarsenedFine;
      coarsen( dblCoarsenedFine, finedbl, m_refRatio[a_ilev] );
      LevelData<FArrayBox> coverLDF(dblCoarsenedFine,nc,IntVect::Unit); // use this for prolong!
      Copier copier(dblCoarsenedFine,a_horizontalVel[a_ilev]->disjointBoxLayout(),IntVect::Zero);
      m_op[a_ilev]->zeroCovered(*a_resid,coverLDF,copier); // copier needs to be cached!!!
    }
}

/// solve for isothermal ice
/** beta scales sliding coefficient C -- acoef in terms of the ViscousTensorOp
 */
int
PetscIceSolver::solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		       Real& a_initialResidualNorm, 
		       Real& a_finalResidualNorm,
		       const Real a_convergenceMetric,
		       const Vector<LevelData<FArrayBox>* >& a_rhs,
		       const Vector<LevelData<FArrayBox>* >& a_beta,
		       const Vector<LevelData<FArrayBox>* >& a_beta0, // not used
		       const Vector<LevelData<FArrayBox>* >& a_A,
		       const Vector<LevelData<FluxBox>* >& a_muCoef,
		       Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		       Real a_time,
		       int a_lbase, int a_maxLevel )
{
  CH_assert(a_lbase==0);
  CH_assert(m_isOpDefined);
  int returnCode = 0, ilev, ierr;
  const int nc = a_horizontalVel[0]->nComp();
  Real residNorm0,residNorm;
  IntVect ghostVect = a_horizontalVel[0]->ghostVect();

  // copy betas and form faceA
  Vector<RefCountedPtr<LevelData<FluxBox> > > faceAs(a_maxLevel+1);
  for (ilev=0;ilev<=a_maxLevel;ilev++)
    {
      m_Beta[ilev] = a_beta[ilev];
      m_Beta0[ilev] = a_beta0[ilev];

      // initial implementation -- redefine solver for every iteration. 
      RefCountedPtr<LevelData<FluxBox> > faceA(new LevelData<FluxBox>(m_grid[ilev], 
								      a_A[ilev]->nComp(), 
								      IntVect::Zero));
      CellToEdge(*a_A[ilev], *faceA);
      faceAs[ilev] = faceA;
    }
  
  Vector<RefCountedPtr<LevelData<FArrayBox> > > resid(a_maxLevel+1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > tempv(a_maxLevel+1);
  for (ilev=a_lbase;ilev<=a_maxLevel;ilev++)
    {
      RefCountedPtr<LevelData<FArrayBox> > vect(new LevelData<FArrayBox>(m_grid[ilev],nc,ghostVect));
      resid[ilev] = vect;
      RefCountedPtr<LevelData<FArrayBox> > vect2(new LevelData<FArrayBox>(m_grid[ilev],nc,ghostVect));
      tempv[ilev] = vect2;
    }

  if (m_verbosity >= 0)
    {
      // residual
      for (ilev=a_lbase,residNorm0=0.;ilev<=a_maxLevel;ilev++)
	{
	  computeAMRLevelsResidual(resid[ilev], ilev,  a_lbase, a_maxLevel+1, a_horizontalVel, a_rhs);
	  residNorm0 += m_op[ilev]->dotProduct(*resid[ilev],*resid[ilev]);
	  // residNorm0 = max(residNorm0,m_op[ilev]->localMaxNorm(*resid[ilev]);
	}
      residNorm0 = sqrt(residNorm0);
      pout() << "PetscIceSolver::solve: Initial residual |r|_2 = " << residNorm0 << endl;
    }
  
  // FMG solve - coarse grid
  ilev = 0;
  // call solver will create matrix, setup PC
  m_tphi = resid[ilev];
  m_tphi2 = tempv[ilev];
  m_tfaceA = faceAs[ilev];
  m_tcoordSys = a_coordSys[ilev];
  m_ttime = a_time;
#ifdef CH_USE_PETSC
  m_petscSolver[ilev]->m_ctx = (void*)this;
  // does c-f interp for finer grids, computes c-f interp if crs_vel provided
  computeMu( *a_horizontalVel[ilev], *faceAs[ilev], a_coordSys[ilev], 0, ilev, a_time ); 
  // creates Mat and Vecs, creates SNES
  m_petscSolver[ilev]->setup_solver( *a_horizontalVel[ilev] );  
  // set the level to index into this
  ierr = SNESSetApplicationContext( m_petscSolver[ilev]->m_snes, (void*)&ilev );CHKERRQ(ierr);
  m_petscSolver[ilev]->solve( *a_horizontalVel[ilev], *a_rhs[ilev] );
  // this signals for next solve that its new (nonlinear)
  m_petscSolver[ilev]->resetOperator();
#endif

  // go up grid hierarchy
  for (ilev=a_lbase+1;ilev<=a_maxLevel;ilev++)
    {
      pout() << "PetscIceSolver::solve: in FMG level " << ilev << endl;
      // prolongate


      // level solve, just like above


    }

  if (m_verbosity >= 0)
    {
      for (ilev=a_lbase,residNorm=0.;ilev<=a_maxLevel;ilev++)
	{
	  computeAMRLevelsResidual( resid[ilev], ilev,  a_lbase, a_maxLevel+1, a_horizontalVel, a_rhs);
	  residNorm += m_op[ilev]->dotProduct(*resid[ilev],*resid[ilev]);
	}
      residNorm = sqrt(residNorm);
      pout() << "PetscIceSolver::solve: PETSc nonlinear solve done |r|_2 = " << residNorm << ", rate=" << residNorm/residNorm0 << endl;
    }
  //cout.precision(5);
  //pout() << "\t" << it << ") |r|_2 = " << residNorm << ", rate=" << residNorm/lastR << endl;
  //lastR = residNorm;
  
  return returnCode;
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::setDefaultValues()
////////////////////////////////////////////////////////////////////////
void
PetscIceSolver::setDefaultValues()
{
  // use isothermal ice temp from Pattyn(2003)
  m_constThetaVal = 238.15;
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::defineOpFactory()
////////////////////////////////////////////////////////////////////////
void PetscIceSolver::defineOpFactory( RealVect a_Crsdx,
				      const ProblemDomain &a_domainCoar,
				      int a_numLevels )
{
  if( !m_isOpDefined )
    {
      if (SpaceDim == 2)
	{
	  CH_assert(m_opFactoryPtr==NULL);	  
	  Real alpha, beta;
	  BCHolder velSolveBC = m_bc->velocitySolveBC();
  	  
	  for(int ilev=0;ilev<a_numLevels;ilev++)
	    {
	      // so needs to be set to avoid assert failures when setting diag too early
	      DataIterator dit = m_C[ilev]->dataIterator();
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  (*m_C[ilev])[dit].setVal(1.0);
		  (*m_Lambda[ilev])[dit].setVal(1.0);
		  (*m_Mu[ilev])[dit].setVal(1.0);
		}
	    }

	  // OpFactory grabs a pointer to mu,lambda,acoef,etc.
	  getOperatorScaleFactors( alpha, beta );
	  m_opFactoryPtr = new ViscousTensorOpFactory( m_grid, m_Mu, m_Lambda, m_C, alpha, 
						       beta, m_refRatio, a_domainCoar, a_Crsdx[0], 
						       velSolveBC, m_vtopSafety);
	}
      else 
	{
	  MayDay::Error("PetscIceSolver::defineOpFactory not implemented for dim = SpaceDim");
	}
    }
  else 
    {
      MayDay::Error("PetscIceSolver::defineOpFactory called twice???");
    }

  m_isOpDefined = true;
}
////////////////////////////////////////////////////////////////////////
// PetscIceSolver::getOperatorScaleFactors()
////////////////////////////////////////////////////////////////////////
void
PetscIceSolver::getOperatorScaleFactors(Real& a_alpha, Real& a_beta) const
{
  a_alpha = -1.0; // sort of wrong signs, but that's what B does 
  a_beta = 1.0; 
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::updateCoefs()
////////////////////////////////////////////////////////////////////////
void
PetscIceSolver::updateCoefs( LevelData<FArrayBox> &a_horizontalVel, int a_ilev )
{
  computeMu( a_horizontalVel, 
	     *m_tfaceA, 
	     m_tcoordSys, 
	     0, // don't need to update ghosts with coarse grid, this must be done above!
	     a_ilev,
	     m_ttime );
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::computeMu()
//   side effect: sets m_Mu & m_Lambda  & m_C
////////////////////////////////////////////////////////////////////////
// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
void 
PetscIceSolver::computeMu( LevelData<FArrayBox> &a_horizontalVel,
			   const LevelData<FluxBox> &a_faceA, 
			   const RefCountedPtr<LevelSigmaCS> &a_coordSys,
			   LevelData<FArrayBox>* crseVelPtr,
			   int a_ilev,
			   Real a_time)
{
  CH_TIME("PetscIceSolver::computeMu");
  ProblemDomain levelDomain = m_domain[a_ilev];
  const DisjointBoxLayout& levelGrids = m_grid[a_ilev];
  const LevelSigmaCS& levelCS = *a_coordSys;
  LevelData<FArrayBox>& levelVel = a_horizontalVel;
  LevelData<FluxBox>& levelMu = *m_Mu[a_ilev];
  const LevelData<FluxBox>& levelA = a_faceA;
  LevelData<FluxBox>& levelLambda = *m_Lambda[a_ilev];
  const LevelData<FArrayBox>& levelBeta = *m_Beta[a_ilev];
  const LevelData<FArrayBox>& levelBeta0 = *m_Beta0[a_ilev];
  LevelData<FArrayBox>& levelC = *m_C[a_ilev];
  DataIterator dit = levelGrids.dataIterator();

  // this is needed
  levelVel.exchange();

  // first set BC's on vel
  m_bc->velocityGhostBC(levelVel,
			levelCS,
			levelDomain, a_time);

  //slc : qcfi.coarseFineInterp fills the edges of lev > 0 cells
  //but not the corners. We need them filled to compute the
  //rate-of-strain invariant, so here is a bodge for now
  // if (SpaceDim == 2)
  //   {
  //     for (dit.begin(); dit.ok(); ++dit)
  // 	{
  // 	  Box sbox = levelVel[dit].box();
  // 	  sbox.grow(-1);
  // 	  FORT_EXTRAPCORNER2D(CHF_FRA(levelVel[dit]),
  // 			      CHF_BOX(sbox));
  // 	}

  //   }
  
  // // actually need to use a cornerCopier, too...
  // CornerCopier cornerCopier(levelGrids, levelGrids, 
  // 			    levelDomain,levelVel.ghostVect(),
  // 			    true);
  // levelVel.exchange(cornerCopier);

  int nRefCrse = -1;  
  IntVect muGhost = IntVect::Zero;
  m_constRelPtr->computeFaceMu( levelMu,
				levelVel,
				crseVelPtr,
				nRefCrse,
				levelA,
				levelCS,
				levelDomain,
				muGhost);

  // now multiply by ice thickness H
  const LevelData<FluxBox>& faceH = levelCS.getFaceH();
  // Real muMax = 1.23456789e+300;
  // Real muMin = 0.0;
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  FArrayBox& thisMu = levelMu[dit][dir];
	  const Box& box = thisMu.box();
	  	  
	  // FORT_MAXFAB1(CHF_FRA(thisMu),
	  // 	       CHF_CONST_REAL(muMin),
	  // 	       CHF_BOX(box));
	  
	  thisMu.mult(faceH[dit][dir],box,0,0,1);

	  // FORT_MINFAB1(CHF_FRA(thisMu),
	  // 	       CHF_CONST_REAL(muMax),
	  // 	       CHF_BOX(box));
	}    
      
      // also update alpha (or C)
      const Box& gridBox = levelGrids[dit];
      m_basalFrictionRelPtr->computeAlpha
	(levelC[dit], levelVel[dit], levelCS.getThicknessOverFlotation()[dit], 
	 levelBeta[dit], levelCS.getFloatingMask()[dit] ,gridBox);
      
      levelC[dit] += levelBeta0[dit];

// #if CH_SPACEDIM==2
//       {
// 	Real mu0 = 1.0;
// 	Real C0 = 1.0;
	
// 	FORT_ENFORCEWELLPOSEDCELL
// 	  (CHF_FRA1(levelC[dit],0),
// 	   CHF_FRA1(levelMu[dit][0],0),
// 	   CHF_FRA1(levelMu[dit][1],0),
// 	   CHF_CONST_REAL(mu0),
// 	   CHF_CONST_REAL(C0),
// 	   CHF_BOX(levelGrids[dit]));
	
//       }
// #endif

      // lambda = 2*mu
      FluxBox& lambda = levelLambda[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  lambda[dir].copy(levelMu[dit][dir]);
	  lambda[dir] *= 2.0;
	}

    }
}
#include "NamespaceFooter.H"
