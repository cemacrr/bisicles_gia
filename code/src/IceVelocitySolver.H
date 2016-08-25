#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _ICEVELOCITYSOLVER_H_
#define _ICEVELOCITYSOLVER_H_

#include "ConstitutiveRelation.H"
#include "BasalFrictionRelation.H"
#include "RealVect.H"
#include "AMRMultiGrid.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "IceThicknessIBC.H"


#include "NamespaceHeader.H"

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
/** This is a virtual base class to facilitate the implementation of
    different solution algorithms.
*/
class IceVelocitySolver
{
public:
  
  IceVelocitySolver() {m_isDefined = false;}
  
  virtual ~IceVelocitySolver() {;}

  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_FrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_vectRefRatio,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels) = 0;
  
  /// solve for isothermal ice
  /** return value is 0 if solver converges successfully
   */
  //virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
  //                  const Vector<LevelData<FArrayBox>* >& a_rhs,
  //                  const Vector<LevelData<FArrayBox>* >& a_C,
  //                  Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
  //                  Real a_time,
  //                  int a_lbase, int a_maxLevel) = 0;


  /// full solve for non-isothermal ice
  /** return value is 0 if solver converges successfully
   */
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
                    int a_lbase, int a_maxLevel) = 0;

  virtual void setVerbosity(int a_verbosity) {m_verbosity = a_verbosity;}

  /// if relevant, set max number of solver iterations
  /** default implementation does nothing -- it's assumed that 
      derived classes will implement this in a meaningful way 
  */
  virtual void setMaxIterations(int a_max_iter) {;}

  /// if relevant, set solver tolerance
  /** default implementation does nothing -- it's assumed that 
      derived classes will implement this in a meaningful way 
  */
  virtual void setTolerance(Real a_tolerance) {;}

  virtual bool isDefined() {return m_isDefined;}


protected:

  bool m_isDefined;
  int m_verbosity;

};

#include "NamespaceFooter.H"

#endif
