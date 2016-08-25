#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _PICARDSOLVER_H_
#define _PICARDSOLVER_H_

#include "IceVelocitySolver.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "AMRMultiGrid.H"
#include "LinearSolver.H"
#include "MultilevelIceVelOp.H"
#ifdef CH_USE_PETSC
#include "PetscAMRSolver.H"
#endif

#include "NamespaceHeader.H"

enum linearSolverTypes {multigrid = 0,
                        BiCGStab,
                        GMRES,
                        relaxSolver,
                        noOpSolver,
                        petsc,
                        num_solver_types};

enum operatorType {viscousTensorOp = 0,
                   num_operator_types};

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
/** This is a virtual base class to facilitate the implementation of
    different solution algorithms.
*/
class PicardSolver : public IceVelocitySolver
{

public:

  PicardSolver(); 

  virtual ~PicardSolver();

  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_basalFrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_vectRefRatio,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels);

  void setSolverType(int a_solver_type) {m_solver_type = a_solver_type;}

  void setTolerance(Real a_tolerance) {m_solver_tolerance = a_tolerance;}

  /// set "absolute tolerance"
  /** If norm(residual) drops below the absolute tolerance (regardless 
      of the initial residual), we declare ourselves "done" -- this is a 
      simple way to avoid solver hangs.
  */
  void setAbsoluteTolerance(Real a_tolerance) 
  {m_absolute_tolerance = a_tolerance;}
  
  /// sets maximum number of Picard iterations
  void setMaxIterations(int a_max_iter) {m_max_iter = a_max_iter;}
  
  void setVerbosity(int a_verbosity) {m_verbosity = a_verbosity;}
  
  /// full solve for non-isothermal ice
  virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		    Vector<LevelData<FArrayBox>* >& a_calvedIce,
		    Vector<LevelData<FArrayBox>* >& a_addedIce,
		    Vector<LevelData<FArrayBox>* >& a_removedIce,
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

  
protected:

  ProblemDomain m_coarseDomain;
  Vector<ProblemDomain> m_domains;
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > m_opFactoryPtr;
  ConstitutiveRelation* m_constRelPtr;
  BasalFrictionRelation* m_basalFrictionRelPtr;

  Vector<DisjointBoxLayout> m_vectGrids;
  Vector<int> m_vectRefRatio;
  RealVect m_dxCrse;
  int m_numLevels;

  Real m_solver_tolerance;
  /// if residual drops below this number, call ourselves "done"
  Real m_absolute_tolerance;
  int m_max_iter;
  int m_solver_type;
  // bottom solver type for multigrid
  int m_bottom_solver_type;

  Real m_vtopSafety;
  Real m_linearsolver_tolerance;
  int m_max_linearsolver_iter;

  int m_verbosity;
  bool m_isOpDefined;
  bool m_isSolverDefined;

  bool m_writeResidToPlotFile;

  IceThicknessIBC* m_bc;

  MultilevelIceVelOp* m_mlOpPtr;
  BiCGStabSolver<Vector<LevelData<FArrayBox>* > >* m_biCGStabSolverPtr;
  GMRESSolver<Vector<LevelData<FArrayBox>* > >* m_GMRESSolverPtr;
  LinearSolver<LevelData<FArrayBox> >* m_bottomSolverPtr;

  AMRMultiGrid<LevelData<FArrayBox> >* m_MGSolverPtr;
  
#ifdef CH_USE_PETSC
  PetscAMRSolver *m_petscSolver;
#endif

  Vector<RefCountedPtr<LevelData<FluxBox> > > m_vectMu;

  Vector<RefCountedPtr<LevelData<FluxBox> > > m_vectLambda;
  
  // beta is used to scale C below
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_vectBeta;
  
  // C here is the sliding coefficient -- the "acoeff" in terms of the viscousTensorOp
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_vectC;

  Real m_constThetaVal;

  void setDefaultValues();

  void defineLinearSolver();

  void clearLinearSolvers();
  
  /// define operator factory. In 2D, this is a ViscousTensorOp.
  /// We'll get to 3D someday...
  void defineOpFactory();

  ///define multigrid solver
  /** MG solver is needed regardless whether we're doing straight MG
      or MG-preconditioned BiCGStab
  */
  void defineMGSolver(AMRMultiGrid<LevelData<FArrayBox> >& a_mgSolver,
                      LinearSolver<LevelData<FArrayBox> >* a_bottomSolverPtr);
  
  /// define BiCGStab solver
  /** returns pointer to the MultilevelLinearOp used by the solver. 
      As things stand now, we'll need to manage that one ourselves.
  */
  MultilevelIceVelOp*  defineBiCGStabSolver(BiCGStabSolver<Vector<LevelData<FArrayBox>* > >& a_bicgstab,
                                                       AMRMultiGrid<LevelData<FArrayBox> >& a_mgSolver);

  /// define GMRES solver
  /** returns pointer to the MultilevelLinearOp used by the solver. 
      As things stand now, we'll need to manage that one ourselves.
  */
  MultilevelIceVelOp*  defineGMRESSolver(GMRESSolver<Vector<LevelData<FArrayBox>* > >& a_gmres,
                                                    AMRMultiGrid<LevelData<FArrayBox> >& a_mgSolver);


  /// get ViscousTensorOp scaling factors (alpha and beta)
  /** put this in a function because we need them in two separate places */
  void getOperatorScaleFactors(Real& a_alpha, Real& a_beta) const;

  // // compute face-centered coefficients for tensor solver (really
  // // winds up being H*mu) -- isothermal version...
  // void computeMu(Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
  //                Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys, 
  //                Real a_time);

  // compute face-centered coefficients for tensor solver (really
  // winds up being H*mu) -- non-isothermal version...
  void computeMu(Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
                 Vector<LevelData<FluxBox>* >& a_A, 
                 Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                 Real a_time);

};

#include "NamespaceFooter.H"

#endif
