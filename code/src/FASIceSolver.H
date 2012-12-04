#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _FASICESOLVER_H_
#define _FASICESOLVER_H_

#include "IceVelocitySolver.H"
#include "ViscousTensorOp.H"
#include "AMRFAS.H"

#include "NamespaceHeader.H"

///
/**
   Factory to create AMR FASIceViscouseTensorOp(s) -- these two classes should be cloned to add operators
*/
class FASIceViscouseTensorOpFactory: public AMRFAS_LDFOpFactory
{
public:
  //! Constructor.
  FASIceViscouseTensorOpFactory( ConstitutiveRelation* a_constRelPtr,
				 BasalFrictionRelation* a_basalFrictionRelPtr,
				 IceThicknessIBC* a_bc );

  virtual ~FASIceViscouseTensorOpFactory()    
  {
    if(m_VTOFactory) delete m_VTOFactory;
  }

  // overide base class' define to add VTO stuff, create m_VTOFactory
  virtual void define( const ProblemDomain& a_coarseDomain, 
		       const RealVect&      a_crsDx,
		       const Vector<DisjointBoxLayout>& a_grids,
		       const Vector<int>&   a_refRatios
		       );
  
  /// override base
  virtual RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > 
  AMRNewOp(int a_ilev, const DisjointBoxLayout& a_grid );

  // Ice op data, copy into op
  const ConstitutiveRelation* m_constRelPtr;
  const BasalFrictionRelation* m_basalFrictionRelPtr;
  
  ViscousTensorOpFactory* m_VTOFactory;
  IceThicknessIBC* m_bc; // this masks m_bc in base class
};

/**
   ViscouseTensor derived class for FAS operator
*/
class FASIceViscouseTensorOp : public AMRFAS_LDFOp
{
  friend class FASIceViscouseTensorOpFactory;
public:
  //! Constructor.
  FASIceViscouseTensorOp( int a_o, 
			  const DisjointBoxLayout &a_grid,
			  const ConstitutiveRelation*  a_constRelPtr,
			  const BasalFrictionRelation* a_basalFrictionRelPtr,
			  IceThicknessIBC* a_bc
			  );
  
  virtual ~FASIceViscouseTensorOp()
  {
    if(m_VTO) delete m_VTO;
  }

  virtual void restrictState( RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > >,
			      RefCountedPtr<LevelData<FArrayBox> > a_CrsCover, 
			      Copier &a_copier );

  virtual void recomputeState( RefCountedPtr<LevelData<FArrayBox> > a_phi,
			       const LevelData<FArrayBox> *a_phiCoarse );

protected:
  virtual void applyLevel( LevelData<FArrayBox>& a_LofPhi,
			   const LevelData<FArrayBox>& a_phi
			   );

  virtual void reflux(const LevelData<FArrayBox>&        a_phiFine,
		      const LevelData<FArrayBox>&        a_phi,
                      LevelData<FArrayBox>&              a_residual,
		      AMRFASOp<LevelData<FArrayBox> >*   a_finerOp )
  {
    const FASIceViscouseTensorOp * const fop = dynamic_cast<FASIceViscouseTensorOp*>( a_finerOp );
    CH_assert(fop);
    m_VTO->reflux( a_phiFine,
		   a_phi,
		   a_residual,
		   fop->m_VTO );
  }

  virtual void levelGSRB( LevelData<FArrayBox>&       a_phi,
			  const LevelData<FArrayBox>& a_rhs,
			  const LevelData<FArrayBox> *a_phiCoarse,
			  const LevelData<FArrayBox> *a_phiFine,
			  AMRFASOp<LevelData<FArrayBox> > *a_finerOp );

  virtual void levelRich( LevelData<FArrayBox>&       a_phi,
			  const LevelData<FArrayBox>& a_rhs,
			  const LevelData<FArrayBox> *a_phiCoarse,
			  const LevelData<FArrayBox> *a_phiFine,
			  AMRFASOp<LevelData<FArrayBox> > *a_finerOp );

  // compute face-centered coefficients for tensor solver (really
  // winds up being H*mu) -- non-isothermal version...
  void computeMu( LevelData<FArrayBox> &a_vel,
		  const LevelData<FArrayBox>* a_crseVel
		  );

public:  
  // constant op data
  const Real m_constThetaVal;
  const Real m_vtopSafety;

  // passed from solver to factory to op
  const ConstitutiveRelation*  m_constRelPtr;
  const BasalFrictionRelation* m_basalFrictionRelPtr;

  // cached in solver
  Real m_time;
  RefCountedPtr<LevelSigmaCS> m_coordSys;

  // beta is used to scale C below
  RefCountedPtr<LevelData<FArrayBox> > m_Beta;     // set in solver from a_beta
  RefCountedPtr<LevelData<FluxBox> >   m_faceA;    // derived in solver from m_A
  // Viscous tensor operator
  ViscousTensorOp* m_VTO;
  IceThicknessIBC* m_bc; // this masks m_bc in base class
};

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
class FASIceSolver : public IceVelocitySolver, public AMRFAS<LevelData<FArrayBox> >
{

public:

  FASIceSolver()
  {
  }

  virtual ~FASIceSolver()
  {
  }

  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_basalFrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_vectRefRatio,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels);
  
  virtual void setTolerance( Real a_tolerance ) { m_rtol = a_tolerance; }

  /// set "absolute tolerance"
  /** If norm(residual) drops below the absolute tolerance (regardless 
      of the initial residual), we declare ourselves "done" -- this is a 
      simple way to avoid solver hangs.
  */
  //virtual void setAbsoluteTolerance(Real a_tolerance) {m_atol = a_tolerance;}

  /// sets maximum number of FAS iterations
  virtual void setMaxIterations(int a_max_iter) { m_max_iter = a_max_iter; }
  
  virtual void setVerbosity(int a_verbosity) {
    IceVelocitySolver::m_verbosity = a_verbosity;
    AMRFAS<LevelData<FArrayBox> >::m_verbosity = a_verbosity;
  }
  
  /// full solve for non-isothermal ice
  virtual int solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		     Real& a_initialResidualNorm, 
		     Real& a_finalResidualNorm,
		     const Real a_convergenceMetric,
		     const Vector<LevelData<FArrayBox>* >& a_rhs,
		     const Vector<LevelData<FArrayBox>* >& a_C,
		     const Vector<LevelData<FArrayBox>* >& a_C0,
		     const Vector<LevelData<FArrayBox>* >& a_A,
		     const Vector<LevelData<FluxBox>* >& a_muCoef,
		     Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		     Real a_time,
		     int a_lbase, 
		     int a_maxLevel
		     );

protected:
  RefCountedPtr<FASIceViscouseTensorOpFactory> m_opFactoryPtr;
  bool m_writeResidToPlotFile;
};

#include "NamespaceFooter.H"

#endif
