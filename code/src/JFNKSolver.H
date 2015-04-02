#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _JFNKSOLVER_H_
#define _JFNKSOLVER_H_


#include "IceVelocitySolver.H"
#include "NonlinearViscousTensor.H"
#include "AMRMultiGrid.H"
#include "MultilevelIceVelOp.H"
#include "ViscousTensorOp.H"
#ifdef CH_USE_PETSC
#include "PetscAMRSolver.H"
#endif

#include "NamespaceHeader.H"

/// Jacobian-Free Newton Krylov solver for the nonlinear ice-sheet/shelf momentum equations
enum LinearizationMode {PICARD_LINEARIZATION_MODE, JFNK_LINEARIZATION_MODE, MAX_LINEARIZATION_MODE};

/// LinearOp implementation which uses a finite difference approximation
/// to evaluate L v = 1/h * ( f(u + h*v) - f(u)) at the current iterate u, where
/// f is a nonlinear function 
/// f:<Vector<LevelData<FArrayBox>*> -> <Vector<LevelData<FArrayBox>*>
/// which is defined in terms of Chombo's MultilevelLinearOp interface.
/// If we say that an operator L[u] is a linear operator defined
/// at the state u, then f(u) = L[u]u, and so 
/// ( f(u + h*v) - f(u)) = L[u+hv](u+hv) - L[u](u)
/// JFNK methods need a preconditioner: in this case we employ 
/// the multigrid method of  MultilevelLinearOp<FArrayBox>,
/// with its operator factory defined at u. 
class JFNKSolver;


class LinearizedVTOp :  public LinearOp<Vector<LevelData<FArrayBox>*> >
{  
  friend class JFNKSolver;
#ifdef CH_USE_PETSC
  friend class PetscAMRSolver;
#endif
  
  NonlinearViscousTensor* m_u; // current state
  NonlinearViscousTensor* m_uPerturbed; // perturbed state

  Vector<LevelData<FArrayBox>*> m_fu; // current value of f(u)
  Vector<LevelData<FArrayBox>*> m_uplushv; // current value of u + v

  //constant perturbation coefficient h
  Real m_h;
  //parameters that appear in the adaptive calculation of h
  Real m_err, m_umin; bool m_hAdaptive;
  //compute a perturbation coefficient h given a direction a_v
  Real finiteh(const Vector<LevelData<FArrayBox>*>& a_v);
  

  // m_mlOp provides the multigrid precondioner
  // and is used to evaluate f(u)
  MultilevelIceVelOp m_mlOp;
  
  //m_perturbedMlOp is used to evaluate f(u + h*v)
  MultilevelIceVelOp m_perturbedMlOp;
  
  /// bottom solver type for MG solver in MultiLevelIceVelOps
  static int m_bottom_solver_type;

  /// max depth for multigrid solvers
  static int m_MG_solver_depth;

  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_refRatio;
  Vector<ProblemDomain> m_domains;
  Vector<RealVect> m_dxs;
  int m_lBase;

  LinearizationMode m_mode;

  //options related to writing out the residual
  bool m_writeResiduals;
  static int m_residualID;

  void writeResidual
  (const Vector<LevelData<FArrayBox> *>& a_u,
   const Vector<LevelData<FArrayBox> *>& a_residual);

public:

  NonlinearViscousTensor& current() const {return *m_u;}

  virtual ~LinearizedVTOp()
  {

    if (m_uPerturbed != NULL)
      {
	delete m_uPerturbed;
	m_uPerturbed = NULL;
      }

    for (int lev = 0; lev < m_fu.size(); ++lev)
      {
	if (m_fu[lev] != NULL)
	  delete m_fu[lev];
      }
    for (int lev = 0; lev < m_uplushv.size(); ++lev)
      {
	if (m_uplushv[lev] != NULL)
	  delete m_uplushv[lev];	
      }
  }
  LinearizedVTOp(NonlinearViscousTensor* a_currentState , 
	 Vector<LevelData<FArrayBox>*>& a_u,
	 Real a_h,  Real m_err, Real m_umin, bool m_hAdaptive,
	 Vector<DisjointBoxLayout>& a_grids,
	 Vector<int>& a_refRatio,
	 Vector<ProblemDomain>& a_domains,
	 Vector<RealVect>& a_dxs,
	 int a_lBase, 
	 int a_numMGSmooth, 
	 int a_numMGIter,
	 LinearizationMode a_mode
	 );

  // set the current solution u (so that J * v \approx (f(u+hv)-f(u))/h
  virtual void setU(Vector<LevelData<FArrayBox>*>& a_u);
  
  // a_lhs = (f(u+ h * (a_v) ) -f(u))/ h (in JFNK mode)
  // or L[u](a_v) (in Picard mode)
  virtual void applyOp(Vector<LevelData<FArrayBox>*>& a_lhs, 
		       const Vector<LevelData<FArrayBox>*>& a_v, 
		       bool a_homogeneous = false);

  // a_lhs = (f(u+ h * (a_v) )-f(u))/h - a_rhs (in JFNK mode)
  // or L[u](a_v) (in Picard mode)
  virtual void residual(Vector<LevelData<FArrayBox>*>& a_lhs, 
			const Vector<LevelData<FArrayBox>*>& a_v, 
			const Vector<LevelData<FArrayBox>*>& a_rhs, 
			bool a_homogeneous = false);
   
  virtual void outerResidual(Vector<LevelData<FArrayBox>*>& a_lhs, 
			     const Vector<LevelData<FArrayBox>*>& a_u, 
			     const Vector<LevelData<FArrayBox>*>& a_rhs, 
			     bool a_homogeneous = false);

  
  virtual void outerLHS(Vector<LevelData<FArrayBox>*>& a_lhs, 
			     const Vector<LevelData<FArrayBox>*>& a_u,  
			     bool a_homogeneous = false)
  {
    
    m_mlOp.applyOp(a_lhs, a_u, a_homogeneous);
  }
  


  //the rest of LinearOp can be palmed off on m_mlOp.
  virtual void preCond(Vector<LevelData<FArrayBox>*>& a_cor,
		       const Vector<LevelData<FArrayBox>*>& a_residual);

    
  virtual void create(Vector<LevelData<FArrayBox>*>& a_lhs,
		      const Vector<LevelData<FArrayBox>*>& a_rhs)
  { 
    m_mlOp.create(a_lhs , a_rhs);
  }

  virtual void clear(Vector<LevelData<FArrayBox>*>& a_lhs)
  {
    m_mlOp.clear(a_lhs);
  }

  // helpful function which clears memory allocated by create
  virtual void clearStorage(Vector<LevelData<FArrayBox> * >& a_lhs)
  {
    for (int lev=0; lev<a_lhs.size(); lev++)
      {
	if (a_lhs[lev] != NULL)
	  {
	    delete a_lhs[lev];
	    a_lhs[lev] = NULL;
	  }
      }
  }

  virtual void assign(Vector<LevelData<FArrayBox>*>& a_lhs,
		      const Vector<LevelData<FArrayBox>*>& a_rhs)
  { 
    m_mlOp.assign(a_lhs, a_rhs);
  }
    
    
  virtual Real dotProduct(const Vector<LevelData<FArrayBox>*>& a_1,
			  const Vector<LevelData<FArrayBox>*>& a_2)
  {
    return m_mlOp.dotProduct(a_1,a_2);
  }

  virtual void incr(Vector<LevelData<FArrayBox>*>& a_lhs,
		    const Vector<LevelData<FArrayBox>*>& a_x,
		    Real a_scale)
  {
    m_mlOp.incr(a_lhs, a_x, a_scale);
  }

    
  virtual void axby(Vector<LevelData<FArrayBox>*>& a_lhs,
		    const Vector<LevelData<FArrayBox>*>& a_x,
		    const Vector<LevelData<FArrayBox>*>& a_y,
		    Real a_a,
		    Real a_b)
  {
    m_mlOp.axby(a_lhs, a_x, a_y, a_a, a_b);
  }

  virtual void scale(Vector<LevelData<FArrayBox>*>& a_lhs,
		     const Real& a_scale)
  {
    m_mlOp.scale(a_lhs , a_scale);
  }

  virtual Real norm(const Vector<LevelData<FArrayBox>*>& a_rhs,
		    int a_ord)
  {
    return m_mlOp.norm(a_rhs, a_ord);
  }

  virtual void setToZero(Vector<LevelData<FArrayBox>*>& a_lhs)
  {
    m_mlOp.setToZero(a_lhs);
  }
  

};


/// JFNK implemenation of IceVelocitySolver
class JFNKSolver : public IceVelocitySolver
{

  Vector<ProblemDomain> m_domains;
  Vector<DisjointBoxLayout> m_grids;
  Vector<RealVect> m_dxs;
  Vector<int> m_refRatios;

  ConstitutiveRelation* m_constRelPtr;
  BasalFrictionRelation* m_basalFrictionRelPtr;
  IceThicknessIBC* m_bcPtr;

  enum LinearSolverType {Relax,BiCGStab,GMRES,CG,petsc};

  /// petsc solvers are more expensive to set up, so keep one around
#ifdef CH_USE_PETSC
  PetscAMRSolver *m_petscSolver;
#endif

  LinearSolver<Vector<LevelData<FArrayBox>* > >* newLinearSolver(LinearizedVTOp& a_op, 
								 LinearizationMode a_mode);

  /// find the largest w <= 1 such that  the residual norm ||f(a_u + w * a_du) - a_rhs|| is not increased
  /// return residual norm
  Real lineSearch(Vector<LevelData<FArrayBox>* >& a_u, 
		  Vector<LevelData<FArrayBox>* >& a_residual,
		  const Vector<LevelData<FArrayBox>* >& a_rhs,
		  const Vector<LevelData<FArrayBox>* >& a_du, 
		  NonlinearViscousTensor& a_nvt, Real a_minW, 
		  bool a_resetOnFail, bool a_linear);

  ///set up and solve the linear system a_op * a_u = a_rhs
  void linearSolve(LinearizedVTOp& a_op,
		   Vector<LevelData<FArrayBox>* >& a_u,
		   const Vector<LevelData<FArrayBox>* >& a_rhs, 
		   LinearizationMode a_mode);

public:
  LinearSolverType m_linearSolverType;
  //maximum number of Newton steps;
  int m_maxIter;
  // absolute tolerance for convergence test;
  Real m_absTol; 
  // relative tolerance for convergence test; 
  Real m_relTol;
  // BiCGstab solver relative tolerance
  Real m_BiCGStabRelTol;
  //maximum number of BiCGStab iterations
  int m_maxBiCGStabIter;
  // CG solver relative tolerance
  Real m_CGRelTol;
  //maximum number of CG iterations
  int m_maxCGIter;
  // GMRES solver relative tolerance
  Real m_GMRESRelTol;
  //maximum number of GMRES iterations
  int m_maxGMRESIter;
  // Relax solver relative tolerance
  Real m_RelaxRelTol;
  //maximum number of Relax iterations
  int m_maxRelaxIter;
  //hang parameter for RelaxSolver;
  Real m_RelaxHang;

  // type of norm to test convergence against
  int m_normType; 
  // verbosity
  int m_verbosity;
  // safety number for viscous tensor op
  Real m_vtopSafety;
   // relax tolerance for  viscous tensor op
  Real m_vtopRelaxTol;
  //  minimum number of smnooths in viscous tensor op mg relax
  int m_vtopRelaxMinIter;
  // (maximum) number of smooths in mg relax
  int m_numMGSmooth;
  // number of MG v-cycles per BiCGStab iteration
  int m_numMGIter;
  // size of finite differenc
  Real m_h;
  //parameters that appear in the adaptive calculation of h
  Real m_err, m_umin; bool m_hAdaptive;
  // rate at which to switch from Picard to JFNK mode
  Real m_switchRate;
  // minumum number of Picard Iterations before switching to JFNK
  /** This is useful when the first iteration or two of Picard don't 
      get very far, but don't get close enough to the solution to warrant 
      switching to JFNK mode. Default is zero.
  */
  int m_minPicardIter;
  //smallest acceptable factor to try when reducing the newton step length
  Real m_minStepFactor;
  //options to do with writing out residual data
  bool m_writeResiduals;
  
  ///throw away fast ice at the end of each iteration? 
  bool m_eliminateFastIce;
  Real m_eliminateFastIceSpeed;
  Real m_eliminateRemoteIceTol;
  Real m_eliminateRemoteIceMaxIter;
  bool m_eliminateFastIceEdgeOnly;

  Real m_uMaxAbs, m_muMin, m_muMax;
private:

  Real m_convergenceMetric;


public:
  virtual ~JFNKSolver() {
 #ifdef CH_USE_PETSC
    if (m_petscSolver) delete m_petscSolver;
#endif
  }

  JFNKSolver() 
#ifdef CH_USE_PETSC
    : m_petscSolver(NULL)
#endif
  {
    setDefaultParameters();
  }


  virtual void define(const ProblemDomain& a_coarseDomain,
		      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_basalFrictionRel,
		      const Vector<DisjointBoxLayout>& a_vectGrids,
		      const Vector<int>& a_vectRefRatio,
		      const RealVect& a_dxCrse,
		      IceThicknessIBC* a_bc,
		      int a_numLevels);

  // set internal parameters to sensible defaults
  void setDefaultParameters();

  // //isothermal solve -- returns 0 if solved, 1 if not
  // virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
  // 		    const Vector<LevelData<FArrayBox>* >& a_rhs,
  // 		    const Vector<LevelData<FArrayBox>* >& a_C,
  // 		    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
  // 		    Real a_time,
  // 		    int a_lbase, int a_maxLevel);

  bool m_residualOnly;
  
  // void writeIteration(const int a_iter, 
  // 		      const Vector<LevelData<FArrayBox> *> a_u,
  // 		      const Vector<LevelData<FArrayBox> *> a_rhs,
  // 		      const Vector<LevelData<FArrayBox> *> a_residual);

  //IceVelocitySolver full solve
  virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		    Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		    const Real a_convergenceMetric,
                    const Vector<LevelData<FArrayBox>* >& a_rhs,
                    const Vector<LevelData<FArrayBox>* >& a_C,
		    const Vector<LevelData<FArrayBox>* >& a_C0,
                    const Vector<LevelData<FArrayBox>* >& a_A,
		    const Vector<LevelData<FluxBox>* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel) ; 

  //general solve, linear or non-linear
  virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		    Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		    const Real a_convergenceMetric,
		    const bool a_linear,
                    const Vector<LevelData<FArrayBox>* >& a_rhs,
                    const Vector<LevelData<FArrayBox>* >& a_C,
		    const Vector<LevelData<FArrayBox>* >& a_C0,
                    const Vector<LevelData<FArrayBox>* >& a_A,
		    const Vector<LevelData<FluxBox>* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel);

  
protected:

  void imposeMaxAbs(Vector<LevelData<FArrayBox>*>& a_u,
		  Real a_limit);
  
};







#include "NamespaceFooter.H"
#endif
