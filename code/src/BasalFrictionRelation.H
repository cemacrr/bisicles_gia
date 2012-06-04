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
class BasalFrictionRelation
{
public:

  virtual ~BasalFrictionRelation() {;}


  /// computes cell-centered \f$ \alpha \f$ based on the cell-centered velocity
  /// and a coefficient field C such that  \f T_b  = \alpha (u_b, C) u_b \f
  /** 
      \param a_alpha -- \f$C\f$ based on the local velocity field. 
      \param a_basalVel -- Cell-centered basal velocity field.
      \param a_thckOverFlotation -- Cell-centered thckness over flotation field.
      \param a_C: Cell-centered coefficient field
      \param a_box: cell-centered box over which to do this computation
  */

  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_thckOverFlotation,
			    const FArrayBox& a_C,
			    const BaseFab<int>& a_mask,
			    const Box& a_box) const = 0;

  /// computes cell-centered heat source \f$ s = T_b . u_b  / k_ice \f$
  /// based on the cell-centered velocity
  /// and a coefficient field beta such that \f T_b  = \alpha (u_b, C) u_b \f
  /** 
      \param a_dissipation -- \f$C\f$ based on the local velocity field. 
      \param a_basalVel -- Cell-centered basal velocity field.
      \param a_C: Cell-centered coefficient field
      \param a_box: cell-centered box over which to do this computation
  */
  virtual void computeDissipation(FArrayBox& a_dissipation,  
				  const FArrayBox& a_basalVel,
				  const FArrayBox& a_thckOverFlotation,
				  const FArrayBox& a_C,
				  const BaseFab<int>& a_mask,
				  const Box& a_box) const;

  
  /// creates a new copy of subclass objects.
  virtual BasalFrictionRelation* getNewBasalFrictionRelation() const = 0;

  virtual Real power() const = 0;


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
			    const FArrayBox& a_thckOverFlotation,
			    const FArrayBox& a_C,
			    const BaseFab<int>& a_mask,
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



#include "NamespaceFooter.H"

#endif
