#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "FASIceSolver.H" 
#include "ParmParse.H"
#include "QuadCFInterp.H"
#include "CornerCopier.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"
#include "CoarseAverageFace.H"
#include "AMRIO.H"
#include "NamespaceHeader.H"

static bool s_recompute_mu = false; // this does help convergance but it is much slower

///
/**
   FASIceViscouseTensorOp: FAS IceViscouseTensor Op class
*/
FASIceViscouseTensorOp::FASIceViscouseTensorOp( int a_o, 
						const DisjointBoxLayout &a_grid,
						const ConstitutiveRelation*  a_constRelPtr,
						const BasalFrictionRelation* a_basalFrictionRelPtr,
						IceThicknessIBC* a_bc
						) :
  AMRFAS_LDFOp( a_o, a_grid ),
  m_constThetaVal(238.15), // use isothermal ice temp from Pattyn(2003)
  m_vtopSafety(VTOP_DEFAULT_SAFETY),
  m_constRelPtr(a_constRelPtr),
  m_basalFrictionRelPtr(a_basalFrictionRelPtr),
  m_VTO(0),
  m_bc(a_bc)
{
  // beta only has one component...
  m_Beta = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_grid,
									 1, 
									 IntVect::Zero));
}

// ---------------------------------------------------------
// applyLevel - apply operator on one level - do BC's but no exchange or C-F
// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::applyLevel( LevelData<FArrayBox>& a_lhs,
				    const LevelData<FArrayBox>& a_phi )
{
  CH_TIME("FASIceViscouseTensorOp::applyLevel");

  // this is only called from base class apply()
  if( s_recompute_mu )
    {
      LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi; // need to force
      
      computeMu( phi, 0 ); // lets do this a lot for now
      // trick to get VTO to recompute D^-1 (m_relaxCoef)
      m_VTO->setAlphaAndBeta( m_VTO->getAlpha(), m_VTO->getBeta() ); 
      // this is used for G-S and Rich
      scale( m_VTO->m_relaxCoef, m_smoothing_damping_factor );
    }

  // this applies BCs
  m_VTO->applyOp( a_lhs, a_phi, true ); // homogenous ? 
}

// ---------------------------------------------------------
// levelGSRB
// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::levelGSRB( RefCountedPtr<LevelData<FArrayBox> > a_phi,
				   const RefCountedPtr<LevelData<FArrayBox> > a_rhs
				   )
{
  CH_TIME("FASIceViscouseTensorOp::levelGSRB");

  if( s_recompute_mu )
    {
      LevelData<FArrayBox>& phi = *a_phi; 
      
      computeMu( phi, 0 ); // lets do this a lot for now  
      // trick to get VTO to recompute D^-1 (m_relaxCoef)
      m_VTO->setAlphaAndBeta( m_VTO->getAlpha(), m_VTO->getBeta() ); 
      scale( m_VTO->m_relaxCoef, m_smoothing_damping_factor );
    }

  m_VTO->relax( *a_phi, *a_rhs, 1 );
}

// ---------------------------------------------------------
// levelRich
// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::levelRich( RefCountedPtr<LevelData<FArrayBox> > a_phi,
				   const  RefCountedPtr<LevelData<FArrayBox> > a_rhs
				   )
{
  CH_TIME("FASIceViscouseTensorOp::levelRich");
  
  CH_assert(a_phi->isDefined());
  CH_assert(a_rhs->isDefined());
  CH_assert(a_phi->ghostVect() >= IntVect::Unit);
  CH_assert(a_phi->nComp() == a_rhs->nComp());

  {
    Real omega = 1.0; // this is folded into m_VTO->m_relaxCoef
    LevelData<FArrayBox> tmp;
    create( tmp, *a_rhs );
    
    // apply has BCs, exchange and C-F interp done before (eg, in compute state)
    apply( tmp, *a_phi, 0, 0, 0, false ); // applies BCs, calls applyLevel

    // x = x + omega D^-1 (b - Ax) ; D == alaph - 2*beta*Dim/h^2
    axby( tmp, tmp, *a_rhs, -1.0, 1.0 );    // (b - Ax) 

    mult( tmp, m_VTO->m_relaxCoef );       // D^-1(b - Ax) 

    axby( *a_phi, *a_phi, tmp, 1.0, omega ); // X = X + omega D^-1 (b - Ax) 
  }
}

// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::restrictState( RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_fOp, // fine op
				       Copier &a_copier )    // copier from buffer to real distribution
{
  CH_TIME("FASIceViscouseTensorOp::restrictState");
  const FASIceViscouseTensorOp * const fop = dynamic_cast<FASIceViscouseTensorOp*>( &(*a_fOp) );
  CH_assert(fop);
  int nRef = fop->refToCoarser();

  // average material parameters to coarser (this)
  m_time = fop->m_time; // this just gets passed up
  m_coordSys = RefCountedPtr<LevelSigmaCS>(new LevelSigmaCS( m_grid, 
							     m_dx, 
							     *fop->m_coordSys,
							     nRef
							     ));
  const LevelData<FluxBox> &f_faceA = *fop->m_faceA;
  const DisjointBoxLayout& fdbl = fop->m_grid;  
  CH_assert( fdbl.coarsenable(nRef) );
  const int aNC = fop->m_faceA->nComp();
  if( !m_faceA )
    {
      m_faceA = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( m_grid,
									   aNC,
									   IntVect::Unit) );
    }
  
  // restrict stuff that was copied in solver: m_faceA, m_Beta
  CoarseAverageFace face_crs( fdbl, aNC, nRef );
  face_crs.averageToCoarse( *m_faceA, f_faceA );
  
  // restrict m_Beta
  {
    const LevelData<FArrayBox>& f_beta = *fop->m_Beta;
    LevelData<FArrayBox> crs_beta, &c_beta = *m_Beta;
    
    //fop->createCoarsened( crs_beta, f_beta, nRef );
    DisjointBoxLayout dblCoarsenedFine;
    coarsen( dblCoarsenedFine, f_beta.disjointBoxLayout(), nRef );
    crs_beta.define( dblCoarsenedFine, f_beta.nComp(), f_beta.ghostVect() );
    //fop->AMRRestrict( crs_beta, f_beta );
    //crs_beta.copyTo( crs_beta.interval(), c_beta, c_beta.interval(), a_copier );
    fop->AMRRestrict( c_beta, f_beta, crs_beta, a_copier );
  }
}

// ---------------------------------------------------------
bool
FASIceViscouseTensorOp::computeState( RefCountedPtr<LevelData<FArrayBox> > a_phi )
{
  CH_TIME("FASIceViscouseTensorOp::recomputeState");
  if( !s_recompute_mu )
    {
      computeMu( *a_phi, 0 );
      // trick to get VTO to recompute D^-1 (m_relaxCoef)
      m_VTO->setAlphaAndBeta( m_VTO->getAlpha(), m_VTO->getBeta() ); 
      // this is used for G-S and Rich
      scale( m_VTO->m_relaxCoef, m_smoothing_damping_factor );
      return true;
    }
  else return false;
}


// ---------------------------------------------------------
void
FASIceViscouseTensorOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
			       const LevelData<FArrayBox>&        a_phi,
			       LevelData<FArrayBox>&              a_residual,
			       AMRFASOp<LevelData<FArrayBox> >*   a_finerOp )
{
  CH_TIME("FASIceViscouseTensorOp::reflux");
  const FASIceViscouseTensorOp * const fop = dynamic_cast<FASIceViscouseTensorOp*>( a_finerOp );
  CH_assert(fop);

  // this constructs a LevelFluxRegister each time - yikes 
  m_VTO->reflux( a_phiFine,
		 a_phi,
		 a_residual,
		 fop->m_VTO );
}


///
/**
   FASIceViscouseTensorOpFactory: ViscouseTensor derived class for FAS-VTO operator factory
*/
  //! Constructor. 1 dof, 2nd order.
FASIceViscouseTensorOpFactory::FASIceViscouseTensorOpFactory(ConstitutiveRelation* a_constRelPtr,
							     BasalFrictionRelation*a_basalFrictionRelPtr,
							     IceThicknessIBC* a_bc
							     ) : 
  AMRFAS_LDFOpFactory( 1, 2 ),
  m_constRelPtr(a_constRelPtr),
  m_basalFrictionRelPtr(a_basalFrictionRelPtr),
  m_VTOFactory(0), // do not like weak constructors so need to wait on this (need a factory factory)
  m_bc(a_bc) // , m_sizeA(a_sizeA)
{
  AMRFASOpFactory<LevelData<FArrayBox> >::define( a_bc->velocitySolveBC() ); // base class define
}

// ---------------------------------------------------------
//  AMR Factory define 
//    - override base define
//
RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > 
FASIceViscouseTensorOpFactory::AMRNewOp( int a_ilev, const DisjointBoxLayout& a_grid, bool a_isSR_dummy )
{
  RefCountedPtr<FASIceViscouseTensorOp> newOp = 
    RefCountedPtr<FASIceViscouseTensorOp>( new FASIceViscouseTensorOp( 2, 
								       a_grid, 
								       m_constRelPtr,
								       m_basalFrictionRelPtr, 
								       m_bc
								       ) );

  // this is a defines the newOp's base class stuff
  AMRFAS_LDFOpFactory::AMRNewOp( a_ilev, newOp, false );

  // define the VTO part
  CH_assert(newOp->m_VTO == NULL);
  CH_assert(m_VTOFactory != NULL);
  newOp->m_VTO = m_VTOFactory->AMRnewOp( m_domains[a_ilev] ); // base factory will search for index

  return newOp;
}

///
/**
   FASIceSolver (Ice) class
*/
//
// main define method - start of define.  pure virtual from Ice base class so need this.
//
void 
FASIceSolver::define( const ProblemDomain& a_coarseDomain,
		      ConstitutiveRelation* a_constRelPtr,
		      BasalFrictionRelation* a_basalFrictionRelPtr,
		      const Vector<DisjointBoxLayout>& a_vectGrids,
		      const Vector<int>& a_vectRefRatio,
		      const RealVect& a_dxCrse,
		      IceThicknessIBC* a_bc,
		      int a_numLevels 
		      )
{
  CH_TIME("FASIceSolver::define");

  // define factory
  CH_assert(m_opFactoryPtr == NULL);
    
  // for the moment, at least, this only works for dx = dy:
  CH_assert(a_dxCrse[0] == a_dxCrse[1]);

  // I don't have all the levels yet, so defer complete define
  m_opFactoryPtr = RefCountedPtr<FASIceViscouseTensorOpFactory>
    (new FASIceViscouseTensorOpFactory( a_constRelPtr, a_basalFrictionRelPtr, a_bc ) );

  // define base solver class, makes coarse levels as needed, calls FASIceViscouseTensorOpFactory::define
  AMRFAS<LevelData<FArrayBox> >::define( a_coarseDomain, 
					 a_dxCrse, 
					 a_vectGrids, 
					 a_vectRefRatio, 
					 *m_opFactoryPtr 
					 );
}

// ---------------------------------------------------------
//  AMR Factory define 
//    - override base define
//
void FASIceViscouseTensorOpFactory::define( const ProblemDomain& a_coarseDomain, 
					    const RealVect&      a_crsDx,
					    const Vector<DisjointBoxLayout>& a_grids,
					    const Vector<int>&   a_refRatios,
					    int a_nSRGrids
					    )
{
  CH_TIME("FASIceViscouseTensorOpFactory::define");
  
  // call base implimentation
  AMRFASOpFactory<LevelData<FArrayBox> >::define( a_coarseDomain,
						  a_crsDx,
						  a_grids,
						  a_refRatios,
						  a_nSRGrids
						  );

  // set alpha & beta here
  const Real alpha = -1.0;  // sort of wrong signs, but that's what B does 
  const Real beta = 1.0;    

  int num_levels = a_grids.size();
  Vector<RefCountedPtr<LevelData<FluxBox> > >   eta(num_levels);
  Vector<RefCountedPtr<LevelData<FluxBox> > >   lambda(num_levels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > c(num_levels);

  for (int i = 0; i < num_levels; i++)
    {
      eta[i] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( a_grids[i], 
									  1,
									  IntVect::Zero));
 
      lambda[i] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( a_grids[i], 
									     1,
									     IntVect::Zero));

      c[i] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>( a_grids[i], 
									    1,
									    IntVect::Zero));

      // this is used to get relaxation values too early so needs to be set to avoid assert failures!
      DataIterator dit = c[i]->dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  (*c[i])[dit].setVal(1.0);
	  (*lambda[i])[dit].setVal(1.0);
	  (*eta[i])[dit].setVal(1.0);
	}
    }
  // complete define for VTO factory
  CH_assert(m_VTOFactory == NULL);
  BCHolder bc = m_bc->velocitySolveBC();

  m_VTOFactory = new ViscousTensorOpFactory( a_grids,
					     eta,
					     lambda,
					     c,
					     alpha,
					     beta,
					     a_refRatios,
					     a_coarseDomain,
					     a_crsDx[0],
					     bc
					     );
}

/// solve for isothermal ice
/** beta scales sliding coefficient C == acoef in terms of the ViscousTensorOp
 */
int
FASIceSolver::solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		     Real& a_initialResidualNorm, 
		     Real& a_finalResidualNorm,
		     const Real a_convergenceMetric,
		     const Vector<LevelData<FArrayBox>* >& a_rhs,
		     const Vector<LevelData<FArrayBox>* >& a_beta,  // Basal C (a_beta)
		     const Vector<LevelData<FArrayBox>* >& a_beta0, // not used
		     const Vector<LevelData<FArrayBox>* >& a_A,
		     const Vector<LevelData<FluxBox>* >& a_muCoef,  // not used, this is computed
		     Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		     Real a_time,
		     int a_lbase, 
		     int a_maxLevel
		     )
{
  CH_assert(a_maxLevel+1 == a_horizontalVel.size());
  CH_TIME("FASIceSolver::solve");

  // copy data into ops, acting like a factory ...
  FASIceViscouseTensorOp *crs_op = 0, *op;
  for (int lev = 0; lev < a_maxLevel + 1; ++lev, crs_op = op )
    {
      op = dynamic_cast<FASIceViscouseTensorOp*>( &(*getOp(lev)) );
      CH_assert(op);

      if( !op->m_faceA )
	{
	  op->m_faceA = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( op->m_grid,
	  									   a_A[lev]->nComp(),
	  									   IntVect::Zero) );
	}
      //
      LevelData<FArrayBox>& argA = *a_A[lev];
      CellToEdge( argA, *op->m_faceA );

      // copy beta into local storage, m_faceA computed for coarse grids in callback in solver
      LevelData<FArrayBox>& localBeta = *op->m_Beta;
      LevelData<FArrayBox>& argBeta =   *a_beta[lev];
      CH_assert(localBeta.nComp() == argBeta.nComp());

      DataIterator dit = localBeta.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          localBeta[dit].copy(argBeta[dit]);
        }
      
      // cache args
      op->m_time = a_time;
      op->m_coordSys = a_coordSys[lev];
    }
  // solver will copy internal levels

  // solver wants RefCountedPtr, allocate and copy in
  Vector<RefCountedPtr<LevelData<FArrayBox> > > phi(a_maxLevel + 1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > rhs(a_maxLevel + 1);
  for (int lev = 0; lev < a_maxLevel + 1; ++lev)
    { 
      phi[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
      rhs[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);

      op = dynamic_cast<FASIceViscouseTensorOp*>( &(*getOp(lev)) );
      CH_assert(op);

      phi[lev]->define( op->m_grid, a_horizontalVel[0]->nComp(), a_horizontalVel[0]->ghostVect() ); 
      op->assignLocal( *phi[lev], *a_horizontalVel[lev] );
      rhs[lev]->define( op->m_grid, a_rhs[0]->nComp(), IntVect::Zero );
      op->assignLocal( *rhs[lev], *a_rhs[lev] );
    }

  // solve -- should we zero out vel?
  int returnCode;
  a_finalResidualNorm = AMRFAS<LevelData<FArrayBox> >::solve( phi, rhs, &returnCode );

  a_initialResidualNorm = 1.0; // don't have 

  // copy out
  for (int lev = 0; lev < a_maxLevel + 1; ++lev)
    { 
      op = dynamic_cast<FASIceViscouseTensorOp*>( &(*getOp(lev)) );
      CH_assert(op);
      op->assignLocal( *a_horizontalVel[lev], *phi[lev] );
    }
  
  if( !m_avoid_norms && m_plot_residual )
    {
#ifdef CH_USE_HDF5
      string fname = "residual.";
      char suffix[30];
      sprintf(suffix, "%dd.hdf5",SpaceDim);
      fname += suffix;
      
      Vector<string> varNames(SpaceDim);
      varNames[0] = "residual-x";
      varNames[1] = "residual-y";

      Real bogusVal = 1.0;
      // writes all grids
      if( 0 )
	{
	  op = dynamic_cast<FASIceViscouseTensorOp*>( &(*getOp(a_maxLevel)) );
	  Real n1 = op->norm( *phi[a_maxLevel], 0, 0);
	  Real n2 = op->norm( *phi[a_maxLevel], 0, 1);
	  pout() << "FASIceSolver::solve plot residual " << m_residual.size() << " levels, |phi_x|= " << setprecision(15) << n1 << " |phi_y|= " << n2 << endl;
	}

      Vector<LevelData<FArrayBox>*> res( m_residual.size() );
      for (int lev = 0; lev < m_residual.size() ; ++lev)
	{ 
	  res[lev] = &(*m_residual[lev]); 
	}

      WriteAMRHierarchyHDF5(fname,
			    m_opFactoryPtr->m_grids,
			    res,
			    varNames,
			    m_op[0]->m_domain.domainBox(),
			    m_op[0]->dx(),
			    bogusVal,
			    bogusVal,
			    m_opFactoryPtr->m_refRatios,
			    m_residual.size());
#endif // end if HDF5
    }

  return returnCode;
}

// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
//
// side effect: does exchange on a_horizontalVel (otherwise it should be const)
//
void 
FASIceViscouseTensorOp::computeMu( LevelData<FArrayBox>& a_horizontalVel,
				   const LevelData<FArrayBox>* a_crsVel
				   )
{ 
  const DisjointBoxLayout& levelGrids =    m_grid;
  const LevelSigmaCS&      levelCS =      *m_coordSys;
  LevelData<FArrayBox>&    levelVel =      a_horizontalVel;
  LevelData<FluxBox>&      levelMu =      *m_VTO->getEta();
  LevelData<FluxBox>&      levelLambda =  *m_VTO->getLambda();
  const LevelData<FluxBox>&levelA =       *m_faceA;
  const LevelData<FArrayBox>& levelBeta = *m_Beta;
  LevelData<FArrayBox>&    levelC =       *m_VTO->getACoef(); // not intuitive: VTO acoef == solver C

  
  // side effect, does exchange
  levelVel.exchange( m_exchangeCopier );

  //levelVel.exchange();

  // first set BC's on vel
  m_bc->velocityGhostBC( levelVel,
			 levelCS,
			 m_domain, 
			 m_time);

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
  
  // actually need to use a cornerCopier, too ... really ??
  CornerCopier cornerCopier(levelGrids, levelGrids, 
			    m_domain,levelVel.ghostVect(),
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

  IntVect muGhost = IntVect::Zero;
  (*m_constRelPtr).computeFaceMu( levelMu,
				  levelVel,
				  a_crsVel,
				  refToCoarser(),
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
	(levelC[dit], levelVel[dit], levelCS.getThicknessOverFlotation()[dit], levelBeta[dit] ,
	 levelCS.getFloatingMask()[dit],gridBox);

    }
}

#include "NamespaceFooter.H"
