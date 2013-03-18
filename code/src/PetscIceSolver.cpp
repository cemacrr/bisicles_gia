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
#include "CoarseAverageFace.H"
#include "IntInterpF_F.H"
#include "IceConstants.H"
#include "BisiclesF_F.H"
#include "TensorCFInterp.H"
#include "AMRIO.H"

#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::PetscIceSolver()
////////////////////////////////////////////////////////////////////////
PetscIceSolver::PetscIceSolver()
{
  // default constructor leaves things in an undefined state


  m_bc = NULL;
  m_constRelPtr = NULL;
  m_basalFrictionRelPtr = NULL;
  m_opFactoryPtr = NULL;

  // use isothermal ice temp from Pattyn(2003)
  m_constThetaVal = 238.15;

  m_isOpDefined = false;
  m_vtopSafety = VTOP_DEFAULT_SAFETY;

  m_max_its = 20;
  m_rtol = 1.e-6;
  m_atol = 1.e-30;
  m_minPicardIterations = 3;

  ParmParse pp("petsc");
  pp.query("maxIter",m_max_its);
  pp.query("absNLTol",m_atol);
  pp.query("relNLTol",m_rtol);
  pp.query("minPicardIterations",m_minPicardIterations);

  // these ones don't need to be stored (at least for now), but should be set
  int mgAverageType  = CoarseAverageFace::arithmetic;
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  ViscousTensorOp::s_prolongType = mgProlongType;
}
////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::~PetscIceSolver() 
////////////////////////////////////////////////////////////////////////
PetscIceSolver::~PetscIceSolver() 
{
  if (m_opFactoryPtr != NULL)
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

  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  ierr = solver->putPetscInChombo( *tthis->m_twork2, x );     CHKERRQ(ierr);

  if ( tthis->m_tphi0 ) tthis->m_op[*pilev]->incr( *tthis->m_twork2, *tthis->m_tphi0, 1.);
  tthis->updateCoefs( *tthis->m_twork2, (int)*pilev ); 
  if ( tthis->m_tphi0 ) tthis->m_op[*pilev]->incr( *tthis->m_twork2, *tthis->m_tphi0, -1.);

  tthis->m_op[*pilev]->applyOp( *tthis->m_twork1, *tthis->m_twork2 ); 

  ierr = solver->putChomboInPetsc( f, *tthis->m_twork1 );  CHKERRQ(ierr);

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
  
  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  // form Function was just called so do not need to update coefs
  ierr = solver->formMatrix( *prejac ); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (prejac!=jac)
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
  m_fineCover.resize(a_numLevels-1);
  m_projCopier.resize(a_numLevels-1);
  m_restCopier.resize(a_numLevels-1);
  m_Mu.resize(a_numLevels);
  m_Lambda.resize(a_numLevels);
  m_Beta.resize(a_numLevels);
  m_Beta0.resize(a_numLevels);
  m_C.resize(a_numLevels);

  m_domain[0] = a_coarseDomain;
  for (int ilev=0;ilev<a_numLevels;ilev++)
    {
      m_grid[ilev] = a_vectGrids[ilev];
      if ( ilev < a_numLevels-1 )
	{
	  m_refRatio[ilev] = a_vectRefRatio[ilev];
	}

      m_Mu[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grid[ilev], 
									     1, 
									     IntVect::Zero) );
      
      m_Lambda[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grid[ilev],
										 1, 
										 IntVect::Zero) );      
      // C only has one component...
      m_C[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grid[ilev], 
										1,
										IntVect::Zero));
      // not coarsest grid
      if (ilev>0)
	{
	  // make domain
	  m_domain[ilev] = m_domain[ilev-1];
	  m_domain[ilev].refine(m_refRatio[ilev-1]);

	  // make prologation stuff (index to coarse grid)
	  const int nc = 2;
	  const DisjointBoxLayout& finedbl = m_grid[ilev];
	  DisjointBoxLayout dblCoarsenedFine;
	  coarsen( dblCoarsenedFine, finedbl, m_refRatio[ilev-1]);	  
	  m_fineCover[ilev-1]=RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(dblCoarsenedFine,nc,IntVect::Unit));
	  // prolongator copier, from data to cover
	  m_projCopier[ilev-1].define(m_grid[ilev-1],dblCoarsenedFine,IntVect::Unit);
	  // restrict copier used for zero cover
	  m_restCopier[ilev-1].define(dblCoarsenedFine,m_grid[ilev-1], IntVect::Zero);
	}
    }
  
  // create op factory
  defineOpFactory( a_dxCrse, a_coarseDomain, a_numLevels );
  
  // create ops
  for (int ilev=0;ilev<a_numLevels;ilev++)
    {
      // this copies the unset data above, just needed here for dx &crdx.
      m_op[ilev] = RefCountedPtr<ViscousTensorOp>(m_opFactoryPtr->AMRnewOp(m_domain[ilev])); 
    }
}

// compute residaul on all levels, sets C-F ghosts and coarsens fine solutions
//
//
void
PetscIceSolver::computeAMRResidual( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_resid, 
				    const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
				    const Vector<LevelData<FArrayBox>* >& a_rhs,
				    int a_lbase, int a_maxLevel
				    )
{
  for (int ilev=a_lbase;ilev<=a_maxLevel;ilev++)
    {
      if (ilev==a_lbase) 
	{
	  if (ilev==a_maxLevel) // one level
	    {
	      m_op[ilev]->residual(*a_resid[ilev], 
				   *a_horizontalVel[ilev], 
				   *a_rhs[ilev], true ); 
	    }
	  else // coarse grid
	    {
	      m_op[ilev]->AMRResidualNC( *a_resid[ilev], 
					 *a_horizontalVel[ilev+1], 
					 *a_horizontalVel[ilev], 
					 *a_rhs[ilev], true, 
					 &(*m_op[ilev+1]));
	    }
	}
      else if (ilev == a_maxLevel) // fine grid
	{
	  m_op[ilev]->AMRResidualNF( *a_resid[ilev], 
				     *a_horizontalVel[ilev],
				     *a_horizontalVel[ilev-1],
				     *a_rhs[ilev], true );
	}
      else
	{
	  m_op[ilev]->AMRResidual( *a_resid[ilev], 
				   *a_horizontalVel[ilev+1], 
				   *a_horizontalVel[ilev], 
				   *a_horizontalVel[ilev-1], 
				   *a_rhs[ilev], true, 
				   &(*m_op[ilev+1]));
	}
      
      // zero covered
      // if (ilev != a_maxLevel) // not finest grid
      //   {
      //     m_op[ilev]->zeroCovered(*a_resid,*m_fineCover[ilev],m_restCopier[ilev]); 
      //   }
    }
}

void PetscIceSolver::AMRProlong( LevelData<FArrayBox>&       a_fu,
				 const LevelData<FArrayBox>& a_cu,
				 LevelData<FArrayBox>&       a_CrsCover,
				 Copier a_copier,
				 int a_refRatio
				 )
{
CH_TIME("PetscIceSolver::AMRProlong");
  
  DisjointBoxLayout dbl = a_fu.disjointBoxLayout();
  DisjointBoxLayout cdbl = a_CrsCover.disjointBoxLayout();
  
  a_cu.copyTo(a_CrsCover.interval(), a_CrsCover, a_CrsCover.interval(), a_copier);
  
  for ( DataIterator dit = a_fu.dataIterator(); dit.ok(); ++dit )
    {
      FArrayBox& phi =  a_fu[dit];
      FArrayBox& coarse = a_CrsCover[dit];
      Box region = dbl[dit];

      FORT_PROLONGQUAD_ICE(CHF_FRA(phi),
			   CHF_CONST_FRA(coarse),
			   CHF_BOX(region),
			   CHF_CONST_INT(a_refRatio));
    }
}

// Picard solve in residual correction form
void PetscIceSolver::picardSolve_private( int a_ilev,
					  LevelData<FArrayBox> &a_horizontalVel,
					  const LevelData<FArrayBox> &a_rhs,
					  int a_numIts, Real a_norm0, int &a_it )
{
#ifdef CH_USE_PETSC
  for (int ii=0;a_it<a_numIts;a_it++,ii++)
    {
      // keep in loop to refresh solver completely
      Real opAlpha, opBeta;
      PetscSolverViscousTensor<LevelData<FArrayBox> > *solver = new PetscSolverViscousTensor<LevelData<FArrayBox> >();
      getOperatorScaleFactors( opAlpha, opBeta );
      solver->define( &(*m_op[a_ilev]), false ); // dx & crdx
      solver->setVTParams( opAlpha, opBeta,
			   &(*m_C[a_ilev]),
			   &(*m_Mu[a_ilev]),
			   &(*m_Lambda[a_ilev]) );
      solver->m_ctx = (void*)this;

      // update coeficients 
      if ( m_tphi0 ) m_op[a_ilev]->incr( a_horizontalVel, *m_tphi0, 1.);
      updateCoefs( a_horizontalVel, a_ilev ); // needed because called before FormJacobian      
      if ( m_tphi0 ) m_op[a_ilev]->incr( a_horizontalVel, *m_tphi0, -1.);

      // in residual correction form (m_tphi0) and first iteration - no init guess.
      solver->setInitialGuessNonzero(true);
      // linear KSP solve
      solver->solve(a_horizontalVel,a_rhs);

      if (m_verbosity>0)
	{
	  pout() << a_it+1 << "/" << m_max_its <<  ") Picard iteration" << endl;	  
	}
      delete solver;
    }
#endif
}

void PetscIceSolver::jfnkSolve_private( int a_ilev,
					LevelData<FArrayBox> &a_horizontalVel,
					const LevelData<FArrayBox> &a_rhs,
					int a_numIts, Real a_norm0, int &a_it
					)
{  
#ifdef CH_USE_PETSC
  
  for (/* void */;a_it<a_numIts;a_it++)
    {
      Real opAlpha, opBeta, norm;
      getOperatorScaleFactors( opAlpha, opBeta );
      PetscSolverViscousTensor<LevelData<FArrayBox> > *solver = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
      solver->define( &(*m_op[a_ilev]), false ); // dx & crdx
      solver->setVTParams( opAlpha, opBeta,
			   &(*m_C[a_ilev]),
			   &(*m_Mu[a_ilev]),
			   &(*m_Lambda[a_ilev]) );
      solver->setFunctionAndJacobian( FormFunction, FormJacobian ); // NL solve
      solver->m_ctx = (void*)this;

      // update coeficients 
      if ( m_tphi0 ) m_op[a_ilev]->incr( a_horizontalVel, *m_tphi0, 1.);
      updateCoefs( a_horizontalVel, a_ilev ); // needed because called before FormJacobian      
      if ( m_tphi0 ) m_op[a_ilev]->incr( a_horizontalVel, *m_tphi0, -1.);
      
      // creates Mat and Vecs, creates SNES
      solver->setup_solver(a_horizontalVel); // creates SNES
      // set the level to index into this
      SNESSetTolerances(solver->m_snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);
      SNESSetApplicationContext( solver->m_snes,(void*)&a_ilev );
      solver->solve(a_horizontalVel,a_rhs);
      SNESGetFunctionNorm(solver->m_snes,&norm);
      delete solver;
      if (m_verbosity>0)
	{
	  pout() << a_it+1 << "/" << m_max_its <<  ") SNES |r|_2 = " << norm << ", rel norm = " << norm/a_norm0 << "/" << m_rtol << endl;	  
	}
      if (norm/a_norm0 < m_rtol){a_it++; break;}
    }

#endif
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
		       Real a_time, int a_lbase, int a_maxLevel )
{
  CH_assert(a_lbase==0); // todo?
  CH_assert(m_isOpDefined);
  int returnCode = 0, ilev;
  const int nc = a_horizontalVel[0]->nComp();
  Real residNorm0,levelNorm;
  IntVect ghostVect = a_horizontalVel[0]->ghostVect(); // can this be Zero???

  // copy betas and form faceA
  Vector<RefCountedPtr<LevelData<FluxBox> > > faceAs(a_maxLevel+1);
  for (ilev=a_maxLevel;ilev<=a_maxLevel;ilev++)
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
  Vector<RefCountedPtr<LevelData<FArrayBox> > > tempv2(a_maxLevel+1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > tempv3(a_maxLevel+1);
  for (ilev=a_maxLevel;ilev<=a_maxLevel;ilev++)
    {
      RefCountedPtr<LevelData<FArrayBox> > v1(new LevelData<FArrayBox>(m_grid[ilev],nc,IntVect::Zero));
      resid[ilev] = v1;
      RefCountedPtr<LevelData<FArrayBox> > v2(new LevelData<FArrayBox>(m_grid[ilev],nc,ghostVect));
      tempv[ilev] = v2;
      RefCountedPtr<LevelData<FArrayBox> > v3(new LevelData<FArrayBox>(m_grid[ilev],nc,IntVect::Zero));
      tempv2[ilev] = v3;
      RefCountedPtr<LevelData<FArrayBox> > v4(new LevelData<FArrayBox>(m_grid[ilev],nc,ghostVect));
      tempv3[ilev] = v4;
    }
  
  // update coeficients to get correct residuals, sets c-f values
  for (ilev=a_maxLevel;ilev<=a_maxLevel;ilev++)
    {
      computeMu( *a_horizontalVel[ilev],
      		 *faceAs[ilev],
      		 a_muCoef[ilev],
      		 a_coordSys[ilev], 
      		 ilev==a_lbase ? 0 : a_horizontalVel[ilev-1],
		 ilev==a_maxLevel ? 0 : a_horizontalVel[ilev+1], // 0
      		 ilev,
      		 a_time );
    }

  // residual
  ilev = a_maxLevel;
  m_op[ilev]->residual( *resid[ilev], 
			*a_horizontalVel[ilev], 
			*a_rhs[ilev], true ); 
  residNorm0 = levelNorm = m_op[ilev]->norm( *resid[ilev], 2 );

  pout() << "\tPetscIceSolver::solve: Initial AMR residual |r|_2 = " << residNorm0 << endl;
  
  // cache stuff for nonlinear solver
  m_twork1 = tempv2[ilev];
  m_twork2 = tempv3[ilev];
  m_tfaceA = faceAs[ilev];
  m_tmuCoef = a_muCoef[ilev];
  m_tcoordSys = a_coordSys[ilev];
  m_ttime = a_time;
  m_tfineVel = ilev==a_maxLevel ? 0 : a_horizontalVel[ilev+1]; // always 0
  
  int it = 0;
  if (a_lbase==a_maxLevel) // coursest grid
    {      
      m_tphi0 = 0;    // not defect correction form
      m_tcrseVel = 0;

      levelNorm = m_op[ilev]->norm(*a_rhs[ilev], 2); // use |b| to compare against

      // start with Picard solve
      picardSolve_private(ilev,*a_horizontalVel[ilev],*a_rhs[ilev],m_minPicardIterations,levelNorm,it);

      // finish wit jfnk
      jfnkSolve_private(ilev,*a_horizontalVel[ilev],*a_rhs[ilev],m_max_its,levelNorm,it);

      if (m_verbosity>0)pout() << "\t\tLevel 0: " << it << " total nonlinear iterations with " << 
			  m_minPicardIterations << " Picard iterations" << endl;
    }
  else
    {
      m_tphi0 = a_horizontalVel[ilev]; // cache phi_0 to get correct linearization
      m_tcrseVel = a_lbase==ilev ? 0 : a_horizontalVel[ilev-1]; 

      pout() << "PetscIceSolver::solve: in FMG level " << ilev << endl;
      
      // prolongate to ilev - done in AmrIce.cpp, use higher order?
      m_op[ilev]->setToZero( *a_horizontalVel[ilev] );
      //AMRProlong( *a_horizontalVel[ilev], *a_horizontalVel[ilev-1], *m_fineCover[ilev-1], 
      //	  m_projCopier[ilev-1], m_refRatio[ilev-1] );
            
      // start with Picard solve
      levelNorm = m_op[ilev]->norm(*a_rhs[ilev], 2); // use |b| to compare against
      it = 0;
      m_op[ilev]->setToZero( *tempv[ilev] );
      picardSolve_private( ilev, *tempv[ilev], *resid[ilev], m_minPicardIterations, levelNorm, it );
      m_op[ilev]->incr( *a_horizontalVel[ilev], *tempv[ilev], 1.); 

      // put in residual correction form, finish with jfnk
      m_op[ilev]->residual( *resid[ilev], *a_horizontalVel[ilev], *a_rhs[ilev], true );
      pout() << "\t\tLevel "<< ilev <<" |r|_2 = " << m_op[ilev]->norm(*resid[ilev], 2) << endl;

      m_op[ilev]->setToZero(*tempv[ilev]);
      jfnkSolve_private( ilev, *tempv[ilev], *resid[ilev], m_max_its, levelNorm, it );      
      m_op[ilev]->incr( *a_horizontalVel[ilev], *tempv[ilev], 1.);

      if (m_verbosity>0)pout() << "\t\tLevel "<< ilev <<" done with " << it << " nonlinear iterations " << endl;      
    }
  
  if (m_verbosity>0)
    {
      m_op[ilev]->residual( *resid[ilev], *a_horizontalVel[ilev], *a_rhs[ilev], true );
      levelNorm = m_op[ilev]->norm(*resid[ilev], 2);
      pout() << "\t\t|r|_2 on level " << ilev << " = " << levelNorm << endl;
    }
  
  // clean up with full JFNK solve

  
  return returnCode;
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::defineOpFactory()
////////////////////////////////////////////////////////////////////////
void PetscIceSolver::defineOpFactory( RealVect a_Crsdx,
				      const ProblemDomain &a_domainCoar,
				      int a_numLevels )
{
  if ( !m_isOpDefined )
    {
      if (SpaceDim == 2)
	{
	  CH_assert(m_opFactoryPtr==NULL);	  
	  Real alpha, beta;
	  BCHolder velSolveBC = m_bc->velocitySolveBC();
  	  
	  for (int ilev=0;ilev<a_numLevels;ilev++)
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
	     m_tmuCoef,
	     m_tcoordSys, 
	     m_tcrseVel, 
	     m_tfineVel,
	     a_ilev,
	     m_ttime );
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::computeMu()
//   side effect: sets m_Mu & m_Lambda & m_C
////////////////////////////////////////////////////////////////////////
// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
void 
PetscIceSolver::computeMu( LevelData<FArrayBox> &a_horizontalVel,
			   const LevelData<FluxBox> &a_faceA, 
			   const LevelData<FluxBox> *a_muCoef,
			   const RefCountedPtr<LevelSigmaCS> &a_coordSys,
			   LevelData<FArrayBox>* crseVelPtr,
			   LevelData<FArrayBox>* fineVelPtr,
			   int a_ilev,
			   Real a_time)
{
  CH_TIME("PetscIceSolver::computeMu");
  ProblemDomain levelDomain = m_domain[a_ilev];
  const DisjointBoxLayout& levelGrids = m_grid[a_ilev];
  const LevelSigmaCS& levelCS = *a_coordSys;
  LevelData<FArrayBox>& levelVel = a_horizontalVel;
  LevelData<FluxBox>& levelMu = *m_Mu[a_ilev];
pout() << "\tPetscIceSolver::computeMu: crseVelPtr = " << crseVelPtr << endl;
  const LevelData<FluxBox>& levelA = a_faceA;
  LevelData<FluxBox>& levelLambda = *m_Lambda[a_ilev];
  const LevelData<FArrayBox>& levelBeta = *m_Beta[a_ilev];
  const LevelData<FArrayBox>& levelBeta0 = *m_Beta0[a_ilev];
  LevelData<FArrayBox>& levelC = *m_C[a_ilev];
  DataIterator dit = levelGrids.dataIterator();

  // // first thing, if there is a finer level, average-down
  // // the current velocity field
  if ( fineVelPtr )
    {
      CoarseAverage averager(fineVelPtr->getBoxes(),
			     levelGrids,
			     fineVelPtr->nComp(),
			     m_refRatio[a_ilev]);
      
      averager.averageToCoarse(levelVel, *fineVelPtr);
    }

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

  int refToCrs = crseVelPtr ? m_refRatio[a_ilev-1] : -1;
  IntVect muGhost = IntVect::Zero;
  m_constRelPtr->computeFaceMu( levelMu,
				levelVel,
				crseVelPtr,
				refToCrs,
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
	  if (a_muCoef != NULL)
	    {
	      thisMu.mult((*a_muCoef)[dit][dir],box,0,0,1);
	    }
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
