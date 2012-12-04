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

#include "AMRIO.H"
#include "NamespaceHeader.H"

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
  
  m_faceA =  RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( a_grid, 
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
  
  m_VTO->applyOp( a_lhs, a_phi, true ); // homogenous ? 
}

// ---------------------------------------------------------
// levelGSRB
// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::levelGSRB( LevelData<FArrayBox>&       a_phi,
				   const LevelData<FArrayBox>& a_rhs,
				   const LevelData<FArrayBox> *a_phiCoarse,
				   const LevelData<FArrayBox> *a_phiFine,
				   AMRFASOp<LevelData<FArrayBox> > *a_finerOp )
{
  CH_TIME("FASIceViscouseTensorOp::levelGSRB");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi; // fortran chokes on const
  computeMu( phi, a_phiCoarse ); // lets do this a lot for now  
  
  m_VTO->relax( a_phi, a_rhs, 1);
}

// ---------------------------------------------------------
// levelRich
// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::levelRich( LevelData<FArrayBox>&    a_phi,
				const LevelData<FArrayBox>& a_rhs,
				const LevelData<FArrayBox> *a_phiCoarse,
				const LevelData<FArrayBox> *a_phiFine,
				AMRFASOp<LevelData<FArrayBox> > *a_finerOp )
{
  CH_TIME("FASIceViscouseTensorOp::levelRich");
  
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  
  {
    CH_TIME("levelRich::exchange");
    a_phi.exchange( m_exchangeCopier );
  }

  {
    // need to add max(beta(x)) ????
    CH_assert(0);
    Real diag = m_VTO->getAlpha() - m_VTO->getBeta()*2.*SpaceDim/(m_dx[0]*m_dx[0]); 
    Real omega = .666; // ??
    LevelData<FArrayBox> tmp;
    create( tmp, a_rhs );
    
    computeMu( a_phi, a_phiCoarse ); // lets do this a lot for now  

    apply( tmp, a_phi, a_phiCoarse, a_phiFine, a_finerOp, false ); // applies BCs, calls applyLevel

    // x = x + omega D^-1 (b - Ax) -- D == alaph - 2*beta*Dim/h^2
    axby( tmp, tmp, a_rhs, -1.0, 1.0 );
    axby( a_phi, a_phi, tmp, 1.0, omega/diag );
  }

}

// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::restrictState( RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_fOp, // fine op
				       RefCountedPtr<LevelData<FArrayBox> > a_CrsCover,       // buffer
				       Copier &a_copier )    // copier from buffer to real distribution
{
  const FASIceViscouseTensorOp * const fop = dynamic_cast<FASIceViscouseTensorOp*>( &(*a_fOp) );
  CH_assert(fop);

  // average material parameters to coarser (this)
  m_time = fop->m_time;
  m_coordSys = RefCountedPtr<LevelSigmaCS>(new LevelSigmaCS());
  m_coordSys->define( *fop->m_coordSys, m_refToFiner );

  // restrict stuff that was copied in solver: m_faceA, m_Beta, m_A
  fop->AMRRestrict( *a_CrsCover,         // output - this cell A
		    *fop->m_VTO->getACoef() );  // input
  assignCopier( *m_VTO->getACoef(), *a_CrsCover, a_copier ); // move to me,

  // form m_faceA
  CellToEdge( *m_VTO->getACoef(), *m_faceA );

  // restrict m_Beta
  fop->AMRRestrict( *a_CrsCover,     // output
		    *fop->m_Beta );   // input
  assignCopier( *m_Beta, *a_CrsCover, a_copier );
}

// ---------------------------------------------------------
void 
FASIceViscouseTensorOp::recomputeState( RefCountedPtr<LevelData<FArrayBox> > a_phi,
					const LevelData<FArrayBox> *a_phiCoarse )
{
  CH_TIME("FASIceViscouseTensorOp::recomputeState");
  
  computeMu( *a_phi, a_phiCoarse ); // -- we are doing this all the time now
}


///
/**
   FASIceViscouseTensorOpFactory: ViscouseTensor derived class for FAS-VTO operator factory
*/
  //! Constructor. 1 dof, 2nd order.
FASIceViscouseTensorOpFactory::FASIceViscouseTensorOpFactory(ConstitutiveRelation* a_constRelPtr,
							     BasalFrictionRelation*a_basalFrictionRelPtr,
							     IceThicknessIBC* a_bc ) : 
  AMRFAS_LDFOpFactory( 1, 2 ),
  m_constRelPtr(a_constRelPtr),
  m_basalFrictionRelPtr(a_basalFrictionRelPtr),
  m_VTOFactory(0), // do not like weak constructors so need to wait on this (need a factory factory)
  m_bc(a_bc)
{
  AMRFASOpFactory<LevelData<FArrayBox> >::define( a_bc->velocitySolveBC() ); // base class define
}
// ---------------------------------------------------------
//  AMR Factory define 
//    - override base define
//
void FASIceViscouseTensorOpFactory::define( const ProblemDomain& a_coarseDomain, 
					    const RealVect&      a_crsDx,
					    const Vector<DisjointBoxLayout>& a_grids,
					    const Vector<int>&   a_refRatios
					    )
{
  CH_TIME("FASIceViscouseTensorOpFactory::define");
  
  // call base implimentation
  AMRFASOpFactory<LevelData<FArrayBox> >::define(a_coarseDomain,a_crsDx,a_grids,a_refRatios);

  // set alpha & beta here
  const Real alpha = -1.0;  // sort of wrong signs, but that's what B does 
  const Real beta = 1.0;    

  int num_levels = a_grids.size();
  Vector<RefCountedPtr<LevelData<FluxBox> > >   eta(num_levels);
  Vector<RefCountedPtr<LevelData<FluxBox> > >   lambda(num_levels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > acoef(num_levels);

  for (int i = 0; i < num_levels; i++)
    {
      eta[i] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( a_grids[i], 
									  1,
									  IntVect::Zero));

      lambda[i] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>( a_grids[i], 
									     1,
									     IntVect::Zero));

      acoef[i] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>( a_grids[i], 
										1,
										IntVect::Zero));
    }

  // complete define for VTO factory
  CH_assert(m_VTOFactory == NULL);
  BCHolder bc = m_bc->velocitySolveBC();
  m_VTOFactory = new ViscousTensorOpFactory( a_grids,
					     eta,
					     lambda,
					     acoef,
					     alpha,
					     beta,
					     a_refRatios,
					     a_coarseDomain,
					     a_crsDx[0],
					     bc
					     );
}

// ---------------------------------------------------------
//  AMR Factory define 
//    - override base define
//
RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > 
FASIceViscouseTensorOpFactory::AMRNewOp( int a_ilev, const DisjointBoxLayout& a_grid )
{
  RefCountedPtr<FASIceViscouseTensorOp> newOp = 
    RefCountedPtr<FASIceViscouseTensorOp>( new FASIceViscouseTensorOp( 2, 
								       a_grid, 
								       m_constRelPtr,
								       m_basalFrictionRelPtr, 
								       m_bc
								       ) );

  // this is a defines the newOp's base class stuff
  AMRFAS_LDFOpFactory::AMRNewOp( a_ilev, newOp );

  // define the VTO part
  CH_assert(newOp->m_VTO == NULL);
  CH_assert(m_VTOFactory != NULL);
  newOp->m_VTO = m_VTOFactory->MGnewOp( m_domains[a_ilev], // base factory will search for index
					0, // flag to use factories eta,lambda,acoef
					true ); // homogenous ? 
}

///
/**
   FASIceSolver (Ice) class
*/
//
// main define method - start of define.  Virtual from Ice base class so need this.
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

  // I don't have all the levels yet so 
  m_opFactoryPtr = RefCountedPtr<FASIceViscouseTensorOpFactory>
    (new FASIceViscouseTensorOpFactory( a_constRelPtr, a_basalFrictionRelPtr, a_bc ) );
    
  // define base solver class, makes coarse levels as needed, calls factory define with grids
  AMRFAS<LevelData<FArrayBox> >::define( a_coarseDomain, 
					 a_dxCrse, 
					 a_vectGrids, 
					 a_vectRefRatio, 
					 *m_opFactoryPtr 
					 );
}

/// solve for isothermal ice
/** beta scales sliding coefficient C -- acoef in terms of the ViscousTensorOp
 */
int
FASIceSolver::solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		  Real& a_initialResidualNorm, 
		  Real& a_finalResidualNorm,
		  const Real a_convergenceMetric,
		  const Vector<LevelData<FArrayBox>* >& a_rhs,
		  const Vector<LevelData<FArrayBox>* >& a_beta,  // Basal C (a_C)
		  const Vector<LevelData<FArrayBox>* >& a_beta0, // not used
		  const Vector<LevelData<FArrayBox>* >& a_A,
		  const Vector<LevelData<FluxBox>* >& a_muCoef,  // not used, this is computed
		  Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		  Real a_time,
		  int a_lbase, 
		  int a_maxLevel
		  )
{
  CH_assert(a_A[0]->nComp()==1);
  CH_TIME("FASIceSolver::solve");
  
  // copy data into ops, acting like a factory ...
  FASIceViscouseTensorOp *crs_op = 0, *op;
  for (int lev = 0; lev < a_maxLevel + 1; ++lev, crs_op = op )
    {
      op = dynamic_cast<FASIceViscouseTensorOp*>( &(*getOp(lev)) );
      CH_assert(op);

      //
      LevelData<FArrayBox>& argA = *a_A[lev];
      CellToEdge( argA, *op->m_faceA );

      // copy A & beta into local storage, m_faceA computed for coarse grids in callback in solver
      LevelData<FArrayBox>& localBeta = *op->m_Beta;
      LevelData<FArrayBox>& argBeta =   *a_beta[lev];
      CH_assert(localBeta.nComp() == argBeta.nComp());
      LevelData<FArrayBox>& localA = *op->m_VTO->getACoef();

      DataIterator dit = localBeta.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          localBeta[dit].copy(argBeta[dit]);
	  localA[dit].copy(argA[dit]);
        }
    }
  // solver will copy internal levels
  
  Real& residNorm = a_finalResidualNorm;
  Real& initialResid = a_initialResidualNorm;

  // solver wants RefCountedPtr
  Vector<RefCountedPtr<LevelData<FArrayBox> > > phi(a_maxLevel + 1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > rhs(a_maxLevel + 1);
  for (int lev = 0; lev < a_maxLevel + 1; ++lev)
    { 
      phi[lev] = RefCountedPtr<LevelData<FArrayBox> >( a_horizontalVel[lev] );
      rhs[lev] = RefCountedPtr<LevelData<FArrayBox> >( a_rhs[lev] );
    }
  
  int returnCode;
  Real rnorm = AMRFAS<LevelData<FArrayBox> >::solve( phi, rhs, &returnCode, true );

  return returnCode;
}

// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
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
  LevelData<FluxBox>&      levelA =       *m_faceA;
  const LevelData<FArrayBox>& levelBeta = *m_Beta;
  LevelData<FArrayBox>&    levelC =       *m_VTO->getACoef(); // not intuitive: acoef <- C

  // just in case, add an exchange here
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
  
  // actually need to use a cornerCopier, too...
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
