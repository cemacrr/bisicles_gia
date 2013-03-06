#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _PETSCICESOLVER_H_
#define _PETSCICESOLVER_H_

#include "IceVelocitySolver.H"
#include "ViscousTensorOp.H"
#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#endif

#include "NamespaceHeader.H"

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
/** This is a virtual base class to facilitate the implementation of
    different solution algorithms.
*/
class PetscIceSolver : public IceVelocitySolver
{

public:

  PetscIceSolver(); 

  virtual ~PetscIceSolver();

  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_basalFrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_vectRefRatio,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels);

  void setSolverType(int a_solver_type) 
  {
    pout() << "PetscIceSolver::setSolverType=" << a_solver_type << endl;
  }

  void setTolerance(Real a_tolerance) { m_rtol = a_tolerance; }

  // set "absolute tolerance"
  /** If norm(residual) drops below the absolute tolerance (regardless 
      of the initial residual), we declare ourselves "done" -- this is a 
      simple way to avoid solver hangs.
  */
  void setAbsoluteTolerance(Real a_tolerance){ m_rtol = a_tolerance; }

  // sets maximum number of iterations
  void setMaxIterations(int a_max_iter) { m_max_its = a_max_iter; }

  void setVerbosity(int a_verbosity) {m_verbosity = a_verbosity;}

  void computeAMRLevelsResidual( RefCountedPtr<LevelData<FArrayBox> > a_resid, 
				 int a_ilev, int a_lbase, int a_maxLevel, 
				 const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
				 const Vector<LevelData<FArrayBox>* >& a_rhs );

  // solve for isothermal ice -- returns 0 if converged, 1 if not
  virtual int solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		     Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		     const Real a_convergenceMetric,
		     const Vector<LevelData<FArrayBox>* >& a_rhs,
		     const Vector<LevelData<FArrayBox>* >& a_C,
		     const Vector<LevelData<FArrayBox>* >& a_C0,
		     const Vector<LevelData<FArrayBox>* >& a_A,
		     const Vector<LevelData<FluxBox>* >& a_muCoef,
		     Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		     Real a_time,
		     int a_lbase, int a_maxLevel);

  Vector<ProblemDomain> m_domain;
  ViscousTensorOpFactory *m_opFactoryPtr;
  Vector<RefCountedPtr<ViscousTensorOp> > m_op;
protected:
  ConstitutiveRelation* m_constRelPtr;
  BasalFrictionRelation* m_basalFrictionRelPtr;
  Vector<DisjointBoxLayout> m_grid;
  Real m_vtopSafety;
  int m_verbosity;
  bool m_isOpDefined;
  bool m_isSolverDefined;

  IceThicknessIBC* m_bc;
  int m_max_its;
  Real m_rtol;
  Real m_atol;

  Vector<int> m_refRatio;
#ifdef CH_USE_PETSC
public:
  Vector<RefCountedPtr<PetscSolverViscousTensor<LevelData<FArrayBox> > > > m_petscSolver;
protected:
#endif

  // materials
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_Mu;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_Lambda;
  // beta is used to scale C below
  Vector<const LevelData<FArrayBox> *> m_Beta;
  Vector<const LevelData<FArrayBox> *> m_Beta0;
  // C here is the sliding coefficient -- the "acoeff" in terms of the viscousTensorOp
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_C;

  Real m_constThetaVal; // not used?

  void setDefaultValues();
  
  /// define operator factory. In 2D, this is a ViscousTensorOp.
  void defineOpFactory(RealVect,const ProblemDomain &,int);

  /// get ViscousTensorOp scaling factors (alpha and beta)
  /** put this in a function because we need them in two separate places */
  void getOperatorScaleFactors( Real& a_alpha, Real& a_beta ) const;
public:
  // compute face-centered coefficients for tensor solver (really
  // winds up being H*mu) -- non-isothermal version...
  void computeMu( LevelData<FArrayBox> &a_horizontalVel, 
		  const LevelData<FluxBox> &a_A, 
		  const RefCountedPtr<LevelSigmaCS> &a_coordSys,
		  LevelData<FArrayBox>* crseVelPtr,
		  int a_ilev,
		  Real a_time);
  // abstract methods for (PETSc) nonlinear solvers
  virtual void updateCoefs(LevelData<FArrayBox> &a_horizontalVel, int);

  // cache 
  RefCountedPtr<LevelData<FArrayBox> >m_tphi;  // place to cache a vector
  RefCountedPtr<LevelData<FArrayBox> >m_tphi2; // place to cache a vector
  RefCountedPtr<LevelData<FluxBox> >  m_tfaceA;
  RefCountedPtr<LevelSigmaCS> m_tcoordSys;
  Real m_ttime;
};
#include "NamespaceFooter.H"
#endif
