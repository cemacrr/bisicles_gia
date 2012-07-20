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
#include "NamespaceHeader.H"
#include <map>
#if BISICLES_Z == BISICLES_LAYERED

class AMRIceControl
{

  // finest level permitted
  int m_maxFinestLevel;

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
  //mesh spacing on the data level
  RealVect m_dataDx;

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

  
 

  // coefficients of basal friction.
  // m_C = exp(X[0]) * m_COrigin;
  Vector<LevelData<FArrayBox>*> m_COrigin;
  RefCountedPtr<LevelData<FArrayBox> > m_dataLevelCOrigin;
  Vector<LevelData<FArrayBox>*> m_C0;
  Vector<LevelData<FArrayBox>*> m_C;
  Vector<LevelData<FArrayBox>*> m_Ccopy;
  Vector<LevelData<FArrayBox>*> m_lapC;

  // coefficients of mu.
  Vector<LevelData<FArrayBox>*> m_muCoef;
  Vector<LevelData<FArrayBox>*> m_lapMuCoef;
  Vector<LevelData<FluxBox>*> m_faceMuCoef;

  //right-hand side of the velocity equation
  Vector<LevelData<FArrayBox>*> m_rhs;
  //right-hand side of the adjoint velocity equation
  Vector<LevelData<FArrayBox>*> m_adjRhs;
  //velocity (at base of ice)
  Vector<LevelData<FArrayBox>*> m_velb;
  //velocity (at surface of ice)
  Vector<LevelData<FArrayBox>*> m_vels;

  bool m_velocityInitialised;
  //lagrange multipler (at base of ice)
  Vector<LevelData<FArrayBox>*> m_adjVel;

  //observations of velocity data (at surface)
  Vector<LevelData<FArrayBox>*> m_velObs; // AMR Hierarchy of observations
  RefCountedPtr<LevelData<FArrayBox> > m_dataLevelVelObs; // data level observations

  //coefficient of velocity data in objective fucntion (ie 1/sigma^2)
  Vector<LevelData<FArrayBox>*> m_velCoef; 
  RefCountedPtr<LevelData<FArrayBox> > m_dataLevelVelCoef;

  //observations of div(uH) 
  Vector<LevelData<FArrayBox>*> m_divUHObs;
  RefCountedPtr<LevelData<FArrayBox> > m_dataLevelDivUHObs;

  //coefficient of div(UH) in objective funvtion (ie 1/sigma^2)
  Vector<LevelData<FArrayBox>*> m_divUHCoef; 
  RefCountedPtr<LevelData<FArrayBox> > m_dataLevelDivUHCoef;
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
  Real m_boundArgX0, m_boundArgX1, 
    m_velMisfitCoefficient, m_massImbalanceCoefficient,
    m_gradCsqRegularization,m_gradMuCoefsqRegularization,
    m_X0Regularization, m_X1Regularization;
  bool m_muCoefLEQOne, m_muCoefLEQOneGround;
  bool m_muCoefGEQOne, m_muCoefGEQOneGround;

  int m_outerCounter;
  int m_innerCounter;

  

public:
  ~AMRIceControl();
  // set up  geoemtry, define observations etc
  void define(IceThicknessIBC* a_ibcPtr,
	      IceTemperatureIBC* a_tempIBCPtr,
	      RateFactor* a_rateFactor,
	      ConstitutiveRelation* a_constRelPtr, 
	      BasalFrictionRelation* a_bfRelPtr,
	      const RealVect& a_dataDx,
	      RefCountedPtr<LevelData<FArrayBox> > a_dataLevelCOrigin,
	      RefCountedPtr<LevelData<FArrayBox> > a_dataLevelVelObs,
	      RefCountedPtr<LevelData<FArrayBox> > a_dataLevelVelCoef,
	      RefCountedPtr<LevelData<FArrayBox> > a_dataLevelDivUHObs,
	      RefCountedPtr<LevelData<FArrayBox> > a_dataLevelDivUHCoef);


  void solveControl();

  void solveForward();

  

  //members needed by CGOptimize
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

  // maximum number of CG iterations before restart?
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

  void CviscoustoCplastic(LevelData<FArrayBox>& a_Cplastic, 
			  const LevelData<FArrayBox>& a_Cviscous, 
			  const LevelData<FArrayBox>&  a_u,
			  const LevelData<BaseFab<int> >&  a_mask) const;
};

#else
#error AMRIceControl not implemented for  BISICLES_Z != BISICLES_LAYERED

#endif
#include "NamespaceFooter.H"
#endif
