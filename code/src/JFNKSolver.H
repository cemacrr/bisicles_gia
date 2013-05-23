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
#include "AMRMultiGrid.H"
#include "MultilevelIceVelOp.H"
#include "ViscousTensorOp.H"
#ifdef CH_USE_PETSC
#include "PetscAMRSolver.H"
#endif

#include "NamespaceHeader.H"

/// Jacobian-Free Newton Krylov solver for the nonlinear ice-sheet/shelf momentum equations
enum SolverMode {PICARD_SOLVER_MODE, JFNK_SOLVER_MODE, MAX_SOLVER_MODE};

// virtual base for classes which encpasulate a state,
class JFNKState
{
  friend class JFNKSolver;
public:
  virtual ~JFNKState(){};
  
  // set the state 
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_state)=0;

  // retreive the state 
  virtual const Vector<LevelData<FArrayBox>*>& getState() const =0;

  //retrieve the coefficient C
  virtual const Vector<RefCountedPtr< LevelData<FArrayBox> > >& getDragCoef() const =0;
  //retrieve the coefficient mu 
  virtual const Vector<RefCountedPtr< LevelData<FluxBox> > >& getViscosityCoef() const =0;
  // provide a pointer to an AMRLevelOpFactory<LevelData<FArrayBox> >
  // defined at the current state
  virtual RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > opFactoryPtr()=0;
};

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



class JFNKOp :  public LinearOp<Vector<LevelData<FArrayBox>*> >
{  
  friend class JFNKSolver;
#ifdef CH_USE_PETSC
  friend class PetscAMRSolver;
#endif
  //CoefficientProvider* m_coeffProvider;
  
  JFNKState* m_u; // current state
  JFNKState* m_uPerturbed; // perturbed state

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

  SolverMode m_mode;

  //options related to writing out the residual
  bool m_writeResiduals;
  static int m_residualID;

  void writeResidual
  (const Vector<LevelData<FArrayBox> *>& a_u,
   const Vector<LevelData<FArrayBox> *>& a_residual);


public:
  virtual ~JFNKOp()
  {
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
  JFNKOp(JFNKState* a_currentState , 
	 JFNKState* a_peturbedState, 
	 Vector<LevelData<FArrayBox>*>& a_u,
	 Real a_h,  Real m_err, Real m_umin, bool m_hAdaptive,
	 Vector<DisjointBoxLayout>& a_grids,
	 Vector<int>& a_refRatio,
	 Vector<ProblemDomain>& a_domains,
	 Vector<RealVect>& a_dxs,
	 int a_lBase, 
	 int a_numMGSmooth, 
	 int a_numMGIter,
	 SolverMode a_mode
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

//manage coeffients mu(u), lambda(u) = 2*mu(u) and alpha(u) 
//in L[u] = div( mu(theta) * (grad(u) + grad(u)^T) + lambda * div(u)*I)  - alpha*u
//and provide the appropriate ViscousTensorOpFactory
//alpha(u) and mu(u) are determined by BasalFrictionRelation and 
//ConstitutiveRelation members.
class IceJFNKstate : public JFNKState
{
  friend class JFNKSolver;

  Vector<LevelData<FArrayBox>*> m_u;
  Vector<LevelData<FArrayBox>*> m_C;
  Vector<LevelData<FArrayBox>*> m_C0;
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > m_opFactoryPtr;
  
 

  //mesh parameters
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_refRatio;
  Vector<ProblemDomain> m_domains;
  Vector<RealVect> m_dxs;
  int m_finestLevel;
  Vector<RefCountedPtr<LevelSigmaCS > > m_coordSys;
  
  
  //constitutive relations
  const ConstitutiveRelation* m_constRelPtr;
  const BasalFrictionRelation* m_basalFrictionRelPtr;
  
  //boundary conditions
  IceThicknessIBC* m_bcPtr;

  //pointers to temperature fields
  Vector<LevelData<FArrayBox>*> m_A;
  Vector<LevelData<FluxBox>*> m_faceA;

  //additional face-centred coefficient of H*mu, usually NULL
  Vector<LevelData<FluxBox>*> m_muCoef;

  // time : this class doesn't care, but m_bcPtr does
  Real m_time;

  // safety factor for VicousTensorOp
  Real m_vtopSafety;
  //  minimum number of smnooths in viscous tensor op mg relax
  int m_vtopRelaxMinIter;
   // relax tolerance for  viscous tensor op
  Real m_vtopRelaxTol;

  //lower and upper bounds on the viscosity
  Real m_muMin, m_muMax;

  
protected:
  //coefficient storage. 
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_mu;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_lambda;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_alpha;
  

public:

  ~IceJFNKstate();
  

  IceJFNKstate(const Vector<DisjointBoxLayout>& a_grids,
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
	       const Vector<LevelData<FluxBox>*>& a_faceFlowLawCoefficient,
	       const Real a_time, 
	       const Real a_vtopSafety,
	       const int a_vtopRelaxMinIter,
	       const Real a_vtopRelaxTol,
	       const Real a_muMin = 0.0,
	       const Real a_muMax = 1.234567e+300);

 

  //implementation of JFNKState interface 
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_u);
  virtual const Vector<LevelData<FArrayBox>*>& getState() const
  {
    return m_u;
  }
  virtual const Vector<LevelData<FArrayBox>*>& getC() const
  {
    return m_C;
  }
  virtual const Vector<RefCountedPtr<LevelData<FluxBox> > >& getViscosityCoef() const
  {
    return m_mu;
  }
  virtual const Vector<RefCountedPtr<LevelData<FArrayBox> > >& getDragCoef() const
  {
    return m_alpha;
  }

  virtual RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
  opFactoryPtr() {
    return m_opFactoryPtr;
  } 

  //an overloaded setState, because beta changes with time
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_u, 
			const Vector<LevelData<FArrayBox>*>& a_basalFrictionCoefficent)
  {
    m_C = a_basalFrictionCoefficent;
    setState(a_u);
  }

  //allows an additional coefficient f of the face-centered viscosity to be specified,
  // so that H*mu*(grad(u) ... ) is replaced with  f*H*mu*(grad(u) ... ).
  virtual void setFaceViscCoef(const Vector<LevelData<FluxBox>*>& a_muCoef)
  {
    int nlev = min(a_muCoef.size(),m_muCoef.size());

    for (int lev=0; lev < nlev; lev++)
      {
	m_muCoef[lev] = a_muCoef[lev];
      }
  }
  virtual void unsetFaceViscCoef()
  {
    for (int lev=0; lev < m_muCoef.size(); lev++)
      {
	m_muCoef[lev] = NULL;
      }
  }

  //provide the components of the Visocus Tensor G[ij] at cell faces, such that
  //at j-faces we have G[ij] for all i and so on
  virtual void computeViscousTensorFace(const Vector<LevelData<FluxBox>*>& a_viscousTensor);
    


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
  
  Real m_uMaxAbs, m_muMin, m_muMax;
private:

  Real m_convergenceMetric;


public:
  virtual ~JFNKSolver() {;}

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
