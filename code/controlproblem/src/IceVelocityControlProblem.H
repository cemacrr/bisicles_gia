#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _ICEVELOCITYCONTROLPROBLEM_H_
#define _ICEVELOCITYCONTROLPROBLEM_H_

#include "LevelDataOps.H"
#include "Vector.H"
#include "LevelSigmaCS.H"
#include "JFNKSolver.H"
#include "NamespaceHeader.H"
/// class IceVelocityControlProblem defines members sufficent
/// for CGOptimize(IceVelcocityProblem& F, LevelData<FArrayBox>& x) 
class IceVelocityControlProblem
{
  
  Vector<LevelDataOps<FArrayBox> > m_vectOps;
  
  //observed velocities
  Vector<LevelData<FArrayBox>* > m_velObs;

  // 1 / (2 sigma^2) for observed velocities
  Vector<LevelData<FArrayBox>* > m_velCoef;

  //observed div(u H)
  Vector<LevelData<FArrayBox>* > m_divUHObs;

  //solution to the forward problem
  Vector<LevelData<FArrayBox>* > m_vel;

  //rhs for the forward problem
  Vector<LevelData<FArrayBox>* > m_rhs;

  //initial value of basal friction coeffcient
  Vector<LevelData<FArrayBox>* > m_C_origin;

  //basal friction coeffcient C_0
  Vector<LevelData<FArrayBox>* > m_C0;

  //rate factor
  Vector<LevelData<FArrayBox>* > m_A;
  Vector<LevelData<FluxBox>* > m_faceA;
  
  //a bunch of workspaces
  Vector<LevelData<FArrayBox>* > m_wk2compA;
  
  //boundary conditions for the forward problem
  IceThicknessIBC* m_bcPtr;

  //constitutive relation for the forward problem
  ConstitutiveRelation* m_constRelPtr;

  //basal friction relation for the forward problem
  BasalFrictionRelation* m_basalFrictionRelPtr;

  //solution to the adjoint problem;
  Vector<LevelData<FArrayBox>* > m_adjVel;
 
  //rhs for the adjoint problem
  Vector<LevelData<FArrayBox>* > m_adjRhs;

  //basal friction coeffcient C
  Vector<LevelData<FArrayBox>* > m_C;
  //laplacian of basal friction coeffcient C
  Vector<LevelData<FArrayBox>* > m_lapC;
  //cell-centred divergence of uH (model u)
  Vector<LevelData<FArrayBox>* > m_divUH;

  //coefficient of mu, cell- and face-centered
  Vector<LevelData<FArrayBox>* > m_muCoef;
  Vector<LevelData<FluxBox>* > m_faceMuCoef;
  //laplacian of basal friction coeffcient C
  Vector<LevelData<FArrayBox>* > m_lapMuCoef;
  Vector<LevelData<FArrayBox>* > m_zero;
  

  //geometry
  Vector<RefCountedPtr<LevelSigmaCS > > m_coordSys;
  
  //meshes
  Vector<DisjointBoxLayout> m_grids;

  //problem domain on the coarsest level
  ProblemDomain m_crseDomain;

  //mesh spacing on the coarsest level
  RealVect m_crseDx;

  //refinement ratios
  Vector<int> m_refRatios;

  // number of calls to computeGradient
  int m_nComputeGradientCalls;

  // Various coefficients that appear in the objective function
  // 2J = sum[ p * (u - u_obs)^2  + q * ( div (uH) - a)^2 + r grad(C)^2]
  // eventually these may become Vector<LevelData<FArrayBox> *>
  // coeffcient of velocity misfit in objective function J.
  Real m_velMisfitCoefficient;
  // coeffcient of mass imbalance in objective function J.
  Real m_massImbalanceCoefficient;
  // coeffcient of grad(C)^2 (r) in objective function J.
  Real m_gradCsqRegularization;
  // coeffcient of grad(muCoef)^2 (r) in objective function J.
  Real m_gradMuCoefsqRegularization;

public:

  IceVelocityControlProblem
  (const Vector<LevelData<FArrayBox>* >& a_velObs,
   const Vector<LevelData<FArrayBox>* >& a_velCoef,
   const Vector<LevelData<FArrayBox>* >& a_divUHObs,
   const Vector<LevelData<FArrayBox>* >& a_vel,
   const Vector<LevelData<FArrayBox>* >& a_rhs,
   const Vector<LevelData<FArrayBox>* >& a_C_origin,
   const Vector<LevelData<FArrayBox>* >& a_C0,
   const Vector<LevelData<FArrayBox>* >& a_A,
   Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const Vector<DisjointBoxLayout>& a_grids,
   const ProblemDomain& a_crseDomain,
   const Vector<int>& a_refRatios,
   const RealVect& a_crseDx,
   IceThicknessIBC* a_bcPtr,
   ConstitutiveRelation* a_constRelPtr,
   BasalFrictionRelation* a_basalFrictionRelPtr,
   const Real& a_velMisfitCoefficient,
   const Real& a_massImbalanceCoefficient,
   const Real& a_gradCsqRegularization,
   const Real& a_gradMuCoefsqRegularization
   );

  ~IceVelocityControlProblem();


  static void solveForwardProblem
  (Vector<LevelData<FArrayBox>* >& a_u,
   bool a_linear,
   const Vector<LevelData<FArrayBox>* >& a_C,
   const Vector<LevelData<FArrayBox>* >& a_C0,
   const Vector<LevelData<FArrayBox>* >& a_A,
   const Vector<LevelData<FluxBox>* >& a_muCoef,
   const Vector<LevelData<FArrayBox>* >& a_rhs,
   Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const ProblemDomain& a_crseDomain,
   const Vector<DisjointBoxLayout>& a_grids,
   const Vector<int>& a_refRatios,
   IceThicknessIBC* a_bc,
   ConstitutiveRelation* a_constRelPtr, 
   BasalFrictionRelation* a_basalFrictionRelPtr,
   const RealVect& a_crseDx);

  static void solveAdjointProblem
  (Vector<LevelData<FArrayBox>* >& a_u,
   const Vector<LevelData<FArrayBox>* >& a_C,
   const Vector<LevelData<FArrayBox>* >& a_C0,
   const Vector<LevelData<FArrayBox>* >& a_A,
   const Vector<LevelData<FluxBox>* >& a_muCoef,
   const Vector<LevelData<FArrayBox>* >& a_rhs,
   Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const ProblemDomain& a_crseDomain,
   const Vector<DisjointBoxLayout>& a_grids,
   const Vector<int>& a_refRatios,
   IceThicknessIBC* a_bc,
   ConstitutiveRelation* a_constRelPtr, 
   BasalFrictionRelation* a_basalFrictionRelPtr,
   const RealVect& a_crseDx);

  //members required by CGOptimize

  //set a_f = objective a_g = gradient of the objective at a_x
  void computeObjectiveAndGradient
  (Real& a_f,
   Vector<LevelData<FArrayBox>* >& a_g, 
   const  Vector<LevelData<FArrayBox>* >& a_x,
   bool a_inner);
  
  //apply preconditioner s = M^{-1}r
  void preCond(Vector<LevelData<FArrayBox>* >& a_s, 
	       const Vector<LevelData<FArrayBox>* >& a_r);

  // duplicate storage of a_b in a_a
  void create(Vector<LevelData<FArrayBox>* >& a_a, 
	      const  Vector<LevelData<FArrayBox>* >& a_b);
  
  // free (any) storage in a_a
  void free(Vector<LevelData<FArrayBox>* >& a_a);
  void free(Vector<LevelData<FluxBox>* >& a_a);

  void setToZero(Vector<LevelData<FArrayBox>* >& a_a);

  // set a_x = s * a_x
  void scale(Vector<LevelData<FArrayBox>* >& a_x, 
	     const  Real a_s);

  // set a_y = a_x
  void assign(Vector<LevelData<FArrayBox>* >& a_y, 
	      const Vector<LevelData<FArrayBox>* >& a_x);
  
  // set a_y = a_y + a_s * a_x
  void incr(Vector<LevelData<FArrayBox>* >& y, 
	    const Vector<LevelData<FArrayBox>* >& x, Real s);
  
  // return a_y.a_x
  Real dotProduct(Vector<LevelData<FArrayBox>* >& a_y, 
		  const Vector<LevelData<FArrayBox>* >& a_x);

  // how many degrees of freedom are there in x?
  int nDoF(const Vector<LevelData<FArrayBox>* >& x);
};
#include "NamespaceFooter.H"
#endif /*_ICEVELOCITYCONTROLPROBLEM_H_*/
