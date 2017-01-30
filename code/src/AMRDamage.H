#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRDAMAGE_H_
#define _AMRDAMAGE_H_

#include "AmrIce.H"
#include "NamespaceHeader.H"

/// Continuum damage model. 
class AMRDamage
{

  int m_time_step;
  Real m_time;

  /// cell-centered damage 
  Vector<LevelData<FArrayBox>*> m_damage;
  /// AMR mesh hierarchy
  int m_finestLevel;
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_ratio;
  Vector<RealVect> m_dx;
  

public:
  AMRDamage();

  void define(const Vector<DisjointBoxLayout>& a_grids, 
	      const Vector<int>& a_ratio,
	      int a_finestLevel,
	      const RealVect& a_crseDx);

  // advance by timestep a_dt. a_dt must satisfy CFL
  void timestep(Real a_dt, const Vector<LevelData<FArrayBox>* >& a_cellVel, 
		const Vector<LevelData<FluxBox>* >& a_faceVel);

private:
  void computeSource(Real a_dt);

  void computeFlux( Vector<LevelData<FluxBox>* >& a_faceFlux,
		    const Vector<LevelData<FArrayBox>* >& a_damage,
		    const Vector<LevelData<FArrayBox>* >& a_cellVel, 
		    const Vector<LevelData<FluxBox>* >& a_faceVel,
		    Real a_dt);

  void updateDamage(Vector<LevelData<FArrayBox>* >& a_damage,
		    const Vector<LevelData<FluxBox>* >& a_faceFlux,
		    Real a_dt);
  
  
};

/// Contains an AMRDamage, responds to changes in the ice sheet state
class DamageIceObserver : public AmrIce::Observer 
{
  AMRDamage* m_damagePtr;

public:

  DamageIceObserver();
  
  ~DamageIceObserver();

  const AMRDamage& damage() const;

  void notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce);
};



class DamagePhysIBC : public PhysIBC
{

public:

  DamagePhysIBC()
  {
  }

  void define(const ProblemDomain& a_domain,
	      const Real&          a_dx);
  
  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  PhysIBC* new_physIBC();

  /// Set up initial conditions
  void initialize(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values
  void primBC(FArrayBox&            a_WGdnv,
	      const FArrayBox&      a_Wextrap,
	      const FArrayBox&      a_W,
	      const int&            a_dir,
	      const Side::LoHiSide& a_side,
	      const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one-sided difference
     approximations.  If this function doesn't change them, they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
		     const Real&      a_time);
  
  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  void artViscBC(FArrayBox&       a_F,
		 const FArrayBox& a_U,
		 const FArrayBox& a_divVel,
		 const int&       a_dir,
		 const Real&      a_time);   
};



#include "NamespaceFooter.H"

#endif
