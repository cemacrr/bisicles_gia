
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BASALFRICTIONRELATION_H_
#define _BASALFRICTIONRELATION_H_

///virtual base class for basal Friction relations 
/** sublclasses must, at least, compute a positive C
   given basal velocity u_b such that basal traction 
   \f T_b  = C(u_b, theta_b) u_b \f
*/

#include "FArrayBox.H"
#include "NamespaceHeader.H"

class LevelSigmaCS;
class DataIterator;

class BasalFrictionRelation
{
public:

  virtual ~BasalFrictionRelation() {;}


  /// computes cell-centered \f$ \alpha \f$ based on the cell-centered velocity
  /// and a coefficient field C such that  \f T_b  = \alpha (u_b, C) u_b \f
  /** 
      \param a_alpha -- \f$C\f$ based on the local velocity field. 
      \param a_basalVel -- Cell-centered basal velocity field.
      \param a_C -- Cell-centered coefficient field
      \param a_coords -- whole level geometry data
      \param a_dit -- index needed to retrieve data from a_coords
      \param a_box: cell-centered box over which to do this computation
  */

  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_C,
			    const LevelSigmaCS& a_coords, 
			    const DataIterator& a_dit,
			    const Box& a_box) const = 0;

  /// computes cell-centered heat source \f$ s = T_b . u_b  / k_ice \f$
  /// based on the cell-centered velocity
  /// and a coefficient field beta such that \f T_b  = \alpha (u_b, C) u_b \f
  /** 
      \param a_dissipation -- \f$C\f$ based on the local velocity field. 
      \param a_basalVel -- Cell-centered basal velocity field.
      \param a_C: Cell-centered coefficient field
      \param a_coords -- whole level geometry data
      \param a_dit -- index needed to retrieve data from a_coords
      \param a_box: cell-centered box over which to do this computation
  */
  virtual void computeDissipation(FArrayBox& a_dissipation,  
				  const FArrayBox& a_basalVel,
				  const FArrayBox& a_C,
				  const LevelSigmaCS& a_coords, 
				  const DataIterator& a_dit,
				  const Box& a_box) const;

  
  /// creates a new copy of subclass objects.
  virtual BasalFrictionRelation* getNewBasalFrictionRelation() const = 0;

  virtual Real power() const = 0;


  static BasalFrictionRelation* parseBasalFrictionRelation(const char* a_prefix, int a_recursion );

};


/// implements BasalFrictionRelation to provide a power law
/**
   \f T_b  = - C |u_b|^{m-1} u_b \f
   
   when m = 1, we have the usual linear law.
        m = 0  gives Coulomb Friction
        m = 1/3 is often used (e.g in the MISMIP tests)
	

 */
class BasalFrictionPowerLaw : public BasalFrictionRelation
{
public:

  BasalFrictionPowerLaw(Real a_m, bool a_includeEffectivePressure) 
    : m_m(a_m), m_includeEffectivePressure(a_includeEffectivePressure) {;}
  
  virtual ~BasalFrictionPowerLaw() {;}

  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_C,
			    const LevelSigmaCS& a_coords, 
			    const DataIterator& a_dit,
			    const Box& a_box) const;

  BasalFrictionRelation* getNewBasalFrictionRelation() const{
    return static_cast<BasalFrictionRelation*> (new BasalFrictionPowerLaw(m_m, m_includeEffectivePressure));
  }

 
  virtual Real power() const {return m_m;} 

private:
  Real m_m;
  bool m_includeEffectivePressure;
  BasalFrictionPowerLaw();

};

/// Basal friction relation of the form Tb = min(Tb', a*p*u/|u|)
/**
   Tb' is any BasalFrictionRelation
   a is a constant (for now) coefficient, p is the local (hydrostatic) pressure

   for a ~ 1, Tb = Tb' is applied most places, with the Coulomb law 
   Tb = a*p*u/|u| imposed as the ice approaches flotation

*/ 
class PressureLimitedBasalFrictionRelation : public BasalFrictionRelation
{

public:
  enum Model {Tsai,Leguy,MAX_MODEL};

private:
  Real m_a, m_p;
  Model m_model;
  BasalFrictionRelation* m_bfr; 
  
public:

  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_C,
			    const LevelSigmaCS& a_coords, 
			    const DataIterator& a_dit,
			    const Box& a_box) const;


  PressureLimitedBasalFrictionRelation(Real a_a, Real a_p, Model a_model, 
				       BasalFrictionRelation* a_bfr)
    : m_a(a_a), m_p(a_p), m_model(a_model),  
      m_bfr(a_bfr->getNewBasalFrictionRelation())
  {

  }


  virtual BasalFrictionRelation* getNewBasalFrictionRelation() const
  {
    PressureLimitedBasalFrictionRelation* ptr = 
      new PressureLimitedBasalFrictionRelation(m_a, m_p, m_model, m_bfr);

    return static_cast<BasalFrictionRelation*>(ptr);
  }
    

  virtual Real power() const
  {
    return m_bfr->power();
  }

};



#include "NamespaceFooter.H"

#endif
