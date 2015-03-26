#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef AMRICECONTROL_H
#define AMRICECONTROL_H
//================================================================
// AMRIceVelocity.H
// methods to contruct and solve ice sheet velocity problems
// that do not involve time-stepping
// 
// Primarily aimed at solving a  control problem (similar to 
// that of Vieli and Payne, Ann Glaciol, 2003. vol 26)
//
//=================================================================

#include "ConstitutiveRelation.H"
#include "BasalFrictionRelation.H"
#include "LevelDataIBC.H"
#include "IceTemperatureIBC.H"
#include "Vector.H"
#include "AMRIO.H"
#include "LevelDataOps.H"
#include "SurfaceFlux.H"
#include <map>
#include <iostream>
#include <fstream>
#include "NamespaceHeader.H"

#if BISICLES_Z == BISICLES_LAYERED



class AMRIceControl
{

#define CCOMP 0
#define MUCOMP 1
#define NXCOMP MUCOMP + 1
  // finest level permitted
  int m_maxFinestLevel;

  //finest levels permitted on floating and grounded regions
  int m_maxFinestLevelFloating;
  int m_maxFinestLevelGrounded;

  // current finest level;
  int m_finestLevel;

  //various AMR parameters
  int m_blockFactor,  m_maxBoxSize,  m_nestingRadius, m_tagsGrow;
  Real m_fillRatio, m_maxVelDx, m_maxLapVelDx;

  //grids for each level
  Vector<DisjointBoxLayout> m_grids;

  //refinement ratios
  Vector<int> m_refRatio;
  
  //mesh spacings
  Vector<RealVect> m_dx;
  //mesh spacing on the data levels
  Vector<RealVect> m_dataDx;

  //problem domains
  Vector<ProblemDomain> m_domain;

  // provides ice thickness and topography, and velocity bcs
  IceThicknessIBC* m_ibcPtr;

  // provide ice temperature
  IceTemperatureIBC* m_tempIBCPtr;

  // provide viscosity and basal friction calculation
  RateFactor* m_rateFactor; 
  ConstitutiveRelation* m_constRelPtr; 
  BasalFrictionRelation* m_bfRelPtr;

  // co-ordinate system, stores thickness, topography, surfaces, mask
  Vector<RefCountedPtr<LevelSigmaCS > > m_coordSys;

  //Storage for lap(X) and grad(X)^2
  Vector<LevelData<FArrayBox>*> m_lapX;
  Vector<LevelData<FArrayBox>*> m_gradXSq;

  // coefficients of basal friction.
  // m_C = exp(X[0]) * m_COrigin;
  Vector<LevelData<FArrayBox>*> m_COrigin;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_referenceCOrigin;
  Vector<LevelData<FArrayBox>*> m_C0;
  Vector<LevelData<FArrayBox>*> m_C;
  Vector<LevelData<FArrayBox>*> m_Ccopy;
  Vector<LevelData<FArrayBox>*> m_lapC;
  Vector<LevelData<FArrayBox>*> m_gradCSq;
 
  // coefficients of mu.
  // m_muCoef = exp(X[1]) * m_referenceMuCoefOrigin;
  Vector<LevelData<FArrayBox>*> m_muCoefOrigin;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_referenceMuCoefOrigin;
  Vector<LevelData<FArrayBox>*> m_muCoef;
  Vector<LevelData<FArrayBox>*> m_lapMuCoef;
  Vector<LevelData<FArrayBox>*> m_gradMuCoefSq;
  Vector<LevelData<FluxBox>*> m_faceMuCoef;

  //thickness
  //thickness = exp(X[2]) * m_thicknessOrigin;
  Vector<LevelData<FArrayBox>*> m_thicknessOrigin;
  Vector<LevelData<FArrayBox>*> m_pointThicknessData;

  //right-hand side of the velocity equation
  Vector<LevelData<FArrayBox>*> m_rhs;
  //right-hand side of the adjoint velocity equation
  Vector<LevelData<FArrayBox>*> m_adjRhs;
  //velocity (at base of ice)
  Vector<LevelData<FArrayBox>*> m_velb;
  //velocity (at surface of ice)
  Vector<LevelData<FArrayBox>*> m_vels;

  bool m_additionalVelocity;

  bool m_velocityInitialised;
  //lagrange multipler (at base of ice)
  Vector<LevelData<FArrayBox>*> m_adjVel;

  //observations of velocity data (at surface)
  Vector<LevelData<FArrayBox>*> m_velObs; // AMR Hierarchy of observations
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_referenceVelObs; // reference (input) observations

  //coefficient of velocity data in objective fucntion (ie 1/sigma^2)
  Vector<LevelData<FArrayBox>*> m_velCoef; 
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_referenceVelCoef; 

  //coefficient of thickness data in objective fucntion (ie 1/sigma^2)
  Vector<LevelData<FArrayBox>*> m_thkCoef; 
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_referenceThkCoef; 

  //observations of div(uH) 
  Vector<LevelData<FArrayBox>*> m_divUHObs;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_referenceDivUHObs;

  //coefficient of div(UH) in objective funvtion (ie 1/sigma^2)
  Vector<LevelData<FArrayBox>*> m_divUHCoef; 
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  m_referenceDivUHCoef;
  //calculations of div(uH)  
  Vector<LevelData<FArrayBox>*> m_divUH;
  //storage for thickness flux
  Vector<LevelData<FluxBox>*> m_faceThckFlux;

  //components of the (local, unregularized) misift function:
  //the global objective function is the L2-norm of the local.
  Vector<LevelData<FArrayBox>*> m_velocityMisfit;
  Vector<LevelData<FArrayBox>*> m_massBalanceMisfit;
  //components of the (local) penalty functions
  Vector<LevelData<FArrayBox>*> m_TikhonovPenalty;
  Vector<LevelData<FArrayBox>*> m_barrierPenalty;
  
  // temperature
  Vector<LevelData<FArrayBox>*> m_temperature; //cell / layer centers
  Vector<LevelData<FArrayBox>*> m_bTemperature; // cell centered surface temperature
  Vector<LevelData<FArrayBox>*> m_sTemperature; // cell centered basal temperature
  // A in glen's flow law
  Vector<LevelData<FArrayBox>*> m_A;
  Vector<LevelData<FluxBox>*> m_faceA;
  //storage for the viscous tensor at cell faces
  Vector<LevelData<FluxBox>*> m_vtFace;

  //sigma levels;
  Vector<Real> m_faceSigma;

  //standard operations for LevelData<FArrayBox> objects
  Vector<LevelDataOps<FArrayBox> > m_vectOps;

  //control problem parameters
  Real m_boundArgX0, 
    m_boundArgX1, 
    m_velMisfitCoefficient, 
    m_massImbalanceCoefficient,
    m_gradCsqRegularization,
    m_gradMuCoefsqRegularization, 
    m_X0Regularization, 
    m_X1Regularization,
    m_gradX0sqRegularization, 
    m_gradX1sqRegularization;
 
  enum VelocityMisfitType {Speed, Velocity, LogSpeed, MAX_VELOCITY_MISFIT_TYPE};
  VelocityMisfitType m_velMisfitType;

  enum ThicknessUpdateType {GroundedSurfaceConstant, BedrockDownSurfaceDown, MAX_THICKNESS_UPDATE_TYPE};
  ThicknessUpdateType m_evolveThicknessUpdateType;


  Real m_lowerX0, m_upperX0, m_lowerX1, m_upperX1;

  enum BoundMethod {None, Projection, MAX_BOUND_METHOD};
  BoundMethod m_boundMethod;

  int m_outerCounter;
  int m_innerCounter;

  int m_restartInterval;


  //options for writing data
  bool m_writeInnerSteps;
  std::string m_innerStepFileNameBase;
  std::string m_outerStepFileNameBase;

  


public:
  ~AMRIceControl();
  // set up  geometry, define observations etc
  void define(IceThicknessIBC* a_ibcPtr,
	      IceTemperatureIBC* a_tempIBCPtr,
	      RateFactor* a_rateFactor,
	      ConstitutiveRelation* a_constRelPtr, 
	      BasalFrictionRelation* a_bfRelPtr,
	      const Vector<RealVect>& a_dataDx,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >& a_referenceCOrigin,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >& a_referenceMuCoefOrigin,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceVelObs,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceVelCoef,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceThkCoef,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >&a_referenceDivUHObs,
	      const Vector<RefCountedPtr<LevelData<FArrayBox > > >& a_referenceDivUHCoef);


  void solveControl();

  void solveForward();

  

  //members needed by CGOptimize
  void computeObjectiveAndGradient
  (Real& a_fm, Real& a_fp,
   Vector<LevelData<FArrayBox>* >& a_g, 
   const  Vector<LevelData<FArrayBox>* >& a_x, 
   bool a_inner);
  
  void restart();



  //apply preconditioner s = M^{-1}r
  void preCond(Vector<LevelData<FArrayBox>* >& a_s, 
	       const Vector<LevelData<FArrayBox>* >& a_r);

  // duplicate storage of a_b in a_a
  void create(Vector<LevelData<FArrayBox>* >& a_a, 
	      const  Vector<LevelData<FArrayBox>* >& a_b);
  
  template <class T>
  void free(Vector<LevelData<T>*>& a_data)
  {
    for (int lev =0; lev < a_data.size(); lev++)
      {
	if (a_data[lev] != NULL)
	  {
	    delete a_data[lev];
	    a_data[lev] = NULL;
	  }
      }
  }

  // free (any) storage in a_a
  //void free(Vector<LevelData<FluxBox>* >& a_a);
  //void free(Vector<LevelData<FArrayBox>*>& a_a);

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

  // no of degrees of freedom
  int nDoF(const Vector<LevelData<FArrayBox>* >& x);

private:

  template <class T>
  void create(Vector<LevelData<T>*>& a_data, 
	      const int a_ncomp, const IntVect& a_ghost)
  {
  if (a_data.size() <= m_finestLevel)
    a_data.resize(m_finestLevel+1);

  for (int lev =0; lev <= m_finestLevel; lev++)
    {
      if (a_data[lev] != NULL)
	delete a_data[lev];
      a_data[lev] = new LevelData<T>(m_grids[lev], a_ncomp, a_ghost);
    }
    
  }
  
  void checkBounds(Vector<LevelData<FArrayBox>* >& a_X);

  
  //void create(Vector<LevelData<FluxBox>*>& a_data, 
  //	      const int a_ncomp, const IntVect& a_ghost);

  void levelSetup(int a_level);
  
  void levelTag(int a_lev, IntVectSet& a_tagset);
  
  void fill(Vector<LevelData<FArrayBox>*>& a_data, const LevelData<FArrayBox>& a_refData, 
	    const int a_lev, const IntVect& a_ghost = IntVect::Zero);

  void fill(Vector<LevelData<FArrayBox>*>& a_data, 
	    const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_refData, 
	    const Vector<RealVect>& a_refDx, const int a_lev, 
	    const IntVect& a_ghost);

  void solveForwardProblem(Vector<LevelData<FArrayBox>* >& a_u,
			   const bool a_linear,
			   const Vector<LevelData<FArrayBox>* >& a_rhs,
			   const Vector<LevelData<FArrayBox>* >& a_C,
			   const Vector<LevelData<FArrayBox>* >& a_C0,
			   const Vector<LevelData<FArrayBox>* >& a_A,
			   const Vector<LevelData<FluxBox>* >& a_muCoef);
  
  std::string outerStateFile() const;
  std::string innerStateFile() const;
  void writeState(const std::string& a_file, const int a_counter,
		  const Vector<LevelData<FArrayBox>* >& a_x,
		  const Vector<LevelData<FArrayBox>* >& a_g) const;

  void readState(const std::string& a_file, int& a_counter,
		 Vector<LevelData<FArrayBox>* >& a_x);

  
  void helmholtzSolve(Vector<LevelData<FArrayBox>* >& a_phi,
				     const Vector<LevelData<FArrayBox>* >& a_rhs,
				     Real a_alpha, Real a_beta) const;
  
  void helmholtzSolve(Vector<LevelData<FArrayBox>* >& a_phi, 
				     Real a_alpha, Real a_beta) const;



  /// flag whether an evolution is required or not at the start of the 
  /// next call to computeObjectiveAndGradient
  bool m_evolve;
  /// time step parameters for evolution
  int m_evolveSteps; // no of outer iterations (CG then evolve)
  Real m_evolveDt; // timestep unless exceeds CFL
  Real m_evolveCFL; // CFL number : dt < n_CFL * dx / max|u|
  Real m_evolveTime; // time to evolve on each iteration
  Real m_evolveRegularizationCoefficient;
  Real m_evolveLimitDH; // maximum change in thickness
  Real m_evolveMeltRateLengthScale;
  Real m_evolveThicknessLengthScale;
  Real m_evolveObservedVelocityMinVelc; // use observed velocity in the evolve step where velc > m_evolveObservedVelocityMinVelc
  /// evolve the geoemtry on CG restarts ?
  bool m_evolveOnRestart;
  /// evolve geoemtry according to the thickness transport equation
  void evolveGeometry(Real a_dt, Real a_time);
  ///compute a smoothed melt-rate such that m_divUHObs \approx m_divUH in the ice shelf
  void evolveMeltRate();
  ///eliminate remote ice after evolving ?
  bool m_eliminateRemoteIce;
  ///maximum number of elimination iterations
  int m_eliminateRemoteIceMaxIter;
  /// thin ice tolerance
  Real m_eliminateRemoteIceTol;

};

#else
#error AMRIceControl not implemented for  BISICLES_Z != BISICLES_LAYERED

#endif
#include "NamespaceFooter.H"
#endif
