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
#include "QuadCFInterp.H"
#include "CornerCopier.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"
//#include "AMRIO.H"
#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::PetscIceSolver()
////////////////////////////////////////////////////////////////////////
PetscIceSolver::PetscIceSolver()
{
  // default constructor leaves things in an undefined state
  setDefaultValues();
  m_bc = NULL;
  m_OpPtr = NULL;
  m_constRelPtr = NULL;
  m_basalFrictionRelPtr = NULL;
  m_opFactoryPtr = NULL;

  m_isOpDefined = false;
  m_vtopSafety = VTOP_DEFAULT_SAFETY;

  m_refRatio = 0;

  m_max_its = 20;
  m_rtol = 1.e-6;
  m_atol = 1.e-30;
#ifdef CH_USE_PETSC
  m_petscSolver = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
#else
  MayDay::Error("PetscIceSolver::PetscIceSolver called w/o PETSc");
#endif
}
////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::~PetscIceSolver() 
////////////////////////////////////////////////////////////////////////
PetscIceSolver::~PetscIceSolver() 
{
#ifdef CH_USE_PETSC
  if (m_petscSolver != NULL)
    {
      delete m_petscSolver;
      m_petscSolver = NULL;
    }
#endif

  if(m_opFactoryPtr != NULL)
    {
      delete m_opFactoryPtr;
      m_opFactoryPtr = NULL;
    }
  if( m_OpPtr )
    {
      delete m_OpPtr;
      m_OpPtr = 0;
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
.  T - optional user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  f - vector
 */
PetscErrorCode FormFunction( SNES snes, Vec x, Vec f, void *dummy )
{
  CH_TIME("PetscIceSolver::FormFunction");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;
 
  PetscFunctionBegin;

  ierr = SNESGetApplicationContext(snes,(void**)&solver); CHKERRQ(ierr);
  tthis = (PetscIceSolver*)solver->m_ctx;

  ierr = solver->putPetscInChombo( *tthis->m_tphi2, x );     CHKERRQ(ierr);
  tthis->m_tphi2->exchange();
  
  tthis->updateCoefs( *tthis->m_tphi2 ); // needed because called before FormJacobian

  tthis->applyOp( *tthis->m_tphi, *tthis->m_tphi2 ); 

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
.  T - optional user-defined context, as set by SNESSetJacobian()

   Output Parameters:
.  jac - Jacobian matrix
.  prejac - different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode FormJacobian( SNES snes,Vec x,Mat *jac,Mat *prejac,MatStructure *flag, void *dummy )
{
  CH_TIME("PetscIceSolver::FormJacobian");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;

  PetscFunctionBegin;

  ierr = SNESGetApplicationContext(snes,(void**)&solver);CHKERRQ(ierr);
  tthis = (PetscIceSolver*)solver->m_ctx;

  ierr = solver->putPetscInChombo( *tthis->m_tphi2, x );     CHKERRQ(ierr);
  tthis->m_tphi2->exchange(); 

  ierr = solver->formMatrix( *prejac, *tthis->m_tphi2 ); CHKERRQ(ierr);
  
  // not sure why this needs to be called, we don't touch it
  ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

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
  CH_assert(a_numLevels==1); // no m_isDefined ???
  m_domain = a_coarseDomain;
  
  m_constRelPtr = a_constRelPtr;
  m_basalFrictionRelPtr = a_basalFrictionRelPtr;
  m_grid = a_vectGrids[0];
  m_bc = a_bc->new_thicknessIBC();

  m_refRatio = a_vectRefRatio[0];

  m_Mu = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grid, 
								   1, 
								   IntVect::Zero) );

  m_Lambda = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grid, 
								       1, 
								       IntVect::Zero) );
  
  // beta only has one component...
  m_Beta = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grid,
									 1, 
									 IntVect::Zero));

  // C only has one component...
  m_C = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grid, 
								      1,
								      IntVect::Zero));

  ParmParse petscPP("petscSolver");
  //petscPP.query("vtopSafety", m_vtopSafety);

  if ( !m_isOpDefined ) {
    defineOpFactory( a_dxCrse );    CH_assert(!m_OpPtr);
    // this copies the unset data above, just needed here for dx &crdx.
    m_OpPtr = m_opFactoryPtr->MGnewOp( m_domain, 0, true ); 
#ifdef CH_USE_PETSC
    Real opAlpha, opBeta;
    m_petscSolver->define( m_OpPtr, false ); // dx & crdx
    getOperatorScaleFactors( opAlpha, opBeta );
    m_petscSolver->setVTParams( opAlpha, opBeta, &(*m_C), &(*m_Mu), &(*m_Lambda) );
    if(true){
      m_petscSolver->setFunctionAndJacobian( FormFunction, FormJacobian );
    }
#endif
  }
}

/// solve for isothermal ice
/** beta scales sliding coefficient C -- acoef in terms of the ViscousTensorOp
 */
int
PetscIceSolver::solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
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
		      int a_lbase, int a_maxLevel)
{
  CH_assert(a_maxLevel==0);
  CH_assert(m_isOpDefined);
  int returnCode = 0;
  
  // initial implementation -- redefine solver for every 
  // iteration. not ideal, but a good place to start. Eventually grab
  // coefficients from operators and reset in existing solver
  //cell face A, needed to compute mu
  LevelData<FluxBox>* faceA = new LevelData<FluxBox>( m_grid, a_A[0]->nComp(), IntVect::Zero );
  CellToEdge(*a_A[0], *faceA);

  // copy beta into local storage (also not terribly efficient; we'll
  // fix that later as well)
  {
    LevelData<FArrayBox>& localBeta = *m_Beta;
    LevelData<FArrayBox>& argBeta = *a_beta[0];
    CH_assert(localBeta.nComp() == argBeta.nComp());

    DataIterator dit = localBeta.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
	localBeta[dit].copy(argBeta[dit]);
      }
  }

  LevelData<FArrayBox>* vect = new LevelData<FArrayBox>( m_grid, a_horizontalVel[0]->nComp(), IntVect::Zero );
  if(true)
    {
      Real residNorm0;
      // setup matrix stuff
      computeMu( *a_horizontalVel[0], *faceA, a_coordSys[0], a_time ); 
      CH_assert( m_OpPtr ); 
      //delete m_OpPtr;
      // MGnewOp copies mu,lambda,acoef copies stuff 
      //m_OpPtr = m_opFactoryPtr->MGnewOp( m_domain, 0, true );
      if (m_verbosity >= 0)
	{
	  m_OpPtr->residual( *vect, *a_horizontalVel[0], *a_rhs[0] );  
	  residNorm0 = sqrt(m_OpPtr->dotProduct(*vect,*vect));	  
	  pout() << "PetscIceSolver::solve: Initial residual |r|_2 = " << residNorm0 << endl;
	}
      // call solver will create matrix, setup PC
      m_tphi = vect;
      LevelData<FArrayBox>* vect2 = new LevelData<FArrayBox>(m_grid,a_horizontalVel[0]->nComp(),a_horizontalVel[0]->ghostVect());
      m_tphi2 = vect2;
      m_tfaceA = faceA;
      m_tcoordSys = a_coordSys[0];
      m_ttime = a_time;
#ifdef CH_USE_PETSC
      m_petscSolver->m_ctx = (void*)this;
      m_petscSolver->solve( *a_horizontalVel[0], *a_rhs[0] );	
#endif
      delete vect2;
      if (m_verbosity >= 0)
	{
	  m_OpPtr->residual( *vect, *a_horizontalVel[0], *a_rhs[0] );  
	  Real residNorm = sqrt(m_OpPtr->dotProduct(*vect,*vect));	  
	  pout() << "PetscIceSolver::solve: PETSc nonlinear solve done |r|_2 = " << residNorm << ", rate=" << residNorm/residNorm0 << endl;
	}
    }
  else{ // Picard/Newton
    for( int it=0 ; it < m_max_its ; it++ ) // temparary
      {
	computeMu( *a_horizontalVel[0], *faceA, a_coordSys[0], a_time ); 
	CH_assert( m_OpPtr ); 
	//delete m_OpPtr;
	// MGnewOp copies mu,lambda,acoef copies stuff
	//m_OpPtr = m_opFactoryPtr->MGnewOp( m_domain, 0, true );
	m_OpPtr->residual( *vect, *a_horizontalVel[0], *a_rhs[0] );  
	Real residNorm = sqrt(m_OpPtr->dotProduct(*vect,*vect)),lastR;
	if(it==0) 
	  {
	    a_initialResidualNorm = residNorm;
	    if (m_verbosity >= 0)
	      {
		pout() << "PetscIceSolver::solve: iteration 0 |r|_2 = " << residNorm << endl;
	      }
	  }
	else if (m_verbosity >= 0)
	  {
	    cout.precision(5);
	    pout() << "\t" << it << ") |r|_2 = " << residNorm << ", rate=" << residNorm/lastR << endl;
	  }
	lastR = residNorm;
	
	if( residNorm/a_initialResidualNorm < m_rtol || residNorm < m_atol ) break; 
#ifdef CH_USE_PETSC	  
	m_petscSolver->setInitialGuessNonzero( true );
	m_petscSolver->solve( *a_horizontalVel[0], *a_rhs[0] );	
	// this signals for next solve that its new (nonlinear)
	m_petscSolver->resetOperator(); 
#endif
      }
  }

  // clean up temporary storage
  delete vect;
  delete faceA;
  
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
void PetscIceSolver::defineOpFactory( RealVect a_dx )
{
  if( !m_isOpDefined )
    {
      if (SpaceDim == 2)
	{
	  CH_assert(m_opFactoryPtr == NULL);
	  
	  Real alpha, beta;
	  getOperatorScaleFactors( alpha, beta );
	  
	  BCHolder velSolveBC = m_bc->velocitySolveBC();
	  
	  Vector< DisjointBoxLayout > grids(1);                       grids[0] = m_grid;
	  Vector< RefCountedPtr< LevelData< FluxBox > > > eta(1);     eta[0] = m_Mu;
	  Vector< RefCountedPtr< LevelData< FluxBox > > > lambda(1);  lambda[0] = m_Lambda;
	  Vector< RefCountedPtr< LevelData< FArrayBox > > > acoef(1); acoef[0] = m_C;
	  Vector< int > refRatios(1);                                 refRatios[0] = m_refRatio;
	  // OpFactory grabs a pointer to mu,lambda,acoef,etc.
	  m_opFactoryPtr = new ViscousTensorOpFactory( grids, eta, lambda, acoef, alpha, 
						       beta, refRatios, m_domain, a_dx[0], velSolveBC, 
						       m_vtopSafety);
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
  a_alpha = -1.0;
  a_beta = 1.0; 
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::computeMu()
//   side effect: sets m_Mu & m_Lambda  & m_C
////////////////////////////////////////////////////////////////////////
// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
void 
PetscIceSolver::computeMu(LevelData<FArrayBox> &a_horizontalVel,
			  const LevelData<FluxBox> &a_A, 
			  const RefCountedPtr<LevelSigmaCS> &a_coordSys,
			  Real a_time)
{
  CH_TIME("PetscIceSolver::computeMu");
  ProblemDomain levelDomain = m_domain;
  const DisjointBoxLayout& levelGrids = m_grid;
  const LevelSigmaCS& levelCS = *a_coordSys;
  LevelData<FArrayBox>& levelVel = a_horizontalVel;
  LevelData<FluxBox>& levelMu = *m_Mu;
  const LevelData<FluxBox>& levelA = a_A;
  LevelData<FluxBox>& levelLambda = *m_Lambda;
  const LevelData<FArrayBox>& levelBeta = *m_Beta;
  LevelData<FArrayBox>& levelC = *m_C;

  // just in case, add an exchange here
  levelVel.exchange();
  //pout() << "PetscIceSolver::computeMu" << endl;
  // first set BC's on vel
  m_bc->velocityGhostBC(levelVel,
			levelCS,
			levelDomain, a_time);

  //slc : qcfi.coarseFineInterp fills the edges of lev > 0 cells
  //but not the corners. We need them filled to compute the
  //rate-of-strain invariant, so here is a bodge for now
  DataIterator dit = levelGrids.dataIterator();
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
  
  // actually need to use a cornerCopier, too...
  CornerCopier cornerCopier(levelGrids, levelGrids, 
			    levelDomain,levelVel.ghostVect(),
			    true);
  levelVel.exchange(cornerCopier);      

  // constant temperature...
  LevelData<FluxBox> theta(levelGrids, 1, IntVect::Unit);   
  LevelData<FArrayBox> cellTheta(levelGrids, 1, IntVect::Unit);         
  for (dit.begin(); dit.ok(); ++dit)
    {
      theta[dit].setVal(m_constThetaVal);
      cellTheta[dit].setVal(m_constThetaVal);          
    }
  
  LevelData<FArrayBox>* crseVelPtr = NULL;
  int nRefCrse = -1;  
  IntVect muGhost = IntVect::Zero;
  m_constRelPtr->computeFaceMu(levelMu,
			       levelVel,
			       crseVelPtr,
			       nRefCrse,
			       levelA,
			       levelCS,
			       m_domain,
			       muGhost);
  
  // now multiply by ice thickness H
  const LevelData<FluxBox>& faceH = levelCS.getFaceH();

  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  levelMu[dit][dir].mult(faceH[dit][dir],
				 levelMu[dit][dir].box(),0,0,1);
	}    

      // lambda = 2*mu
      FluxBox& lambda = levelLambda[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  lambda[dir].copy(levelMu[dit][dir]);
	  lambda[dir] *= 2.0;
	}
      
      // also update alpha (or C)
      const Box& gridBox = levelGrids[dit];
      m_basalFrictionRelPtr->computeAlpha
	(levelC[dit], levelVel[dit], levelBeta[dit] ,gridBox);
    }
}
#include "NamespaceFooter.H"
