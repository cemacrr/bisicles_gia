#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "BasalFrictionRelation.H"
#include "BasalFrictionRelationF_F.H"
#include "IceConstants.H"
#include "LevelSigmaCS.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"


/// computes cell-centered \f$ s = - T_b . u_b \f$
/// based on the cell-centered velocity
/// and a coefficient field beta such that \f T_b  = - \alpha (u, C) u_b \f
/** 
    \param a_C -- \f$C\f$ based on the local velocity field. 
    \param a_basalVel -- Cell-centered basal velocity field.
    \param a_box: cell-centered box over which to do this computation 
*/
void 
BasalFrictionRelation::computeDissipation(FArrayBox& a_dissipation,  
					  const FArrayBox& a_basalVel,
					  const FArrayBox& a_C,
					  const LevelSigmaCS& a_coords, 
					  const DataIterator& a_dit,
					  const Box& a_box) const
{
    
  computeAlpha(a_dissipation, a_basalVel, a_C, a_coords, a_dit, a_box);
  
  FORT_BFRICTIONAUU(CHF_FRA1(a_dissipation,0),
		    CHF_CONST_FRA(a_basalVel),
		    CHF_BOX(a_box));

}


void
BasalFrictionPowerLaw::computeAlpha(FArrayBox& a_alpha,
				    const FArrayBox& a_basalVel,
				    const FArrayBox& a_C,
				    const LevelSigmaCS& a_coords, 
				    const DataIterator& a_dit,
				    const Box& a_box) const
{
  const Real mm1 = m_m - 1.0;
  const FArrayBox& thckOverFlotation = a_coords.getThicknessOverFlotation()[a_dit];
  const BaseFab<int>& mask = a_coords.getFloatingMask()[a_dit];
  if (std::abs(mm1) > 1.0e-3)
    {
      const Real pexp = (m_includeEffectivePressure)?m_m:0.0;
      FORT_BFRICTIONPOWER(CHF_FRA1(a_alpha,0),
			  CHF_CONST_FRA(a_basalVel),
			  CHF_CONST_FRA1(a_C,0),
			  CHF_CONST_FRA1(thckOverFlotation,0),
			  CHF_CONST_FIA1(mask,0),
			  CHF_CONST_REAL(mm1),
			  CHF_CONST_REAL(pexp),
			  CHF_BOX(a_box));	
    }
  else
    {
      a_alpha.copy(a_C,a_box);
      if (m_includeEffectivePressure)
	a_alpha *= thckOverFlotation;
    }
  
}


void
PressureLimitedBasalFrictionRelation::computeAlpha
(FArrayBox& a_alpha,
 const FArrayBox& a_basalVel,
 const FArrayBox& a_C,
 const LevelSigmaCS& a_coords, 
 const DataIterator& a_dit,
 const Box& a_box) const
{
  m_bfr->computeAlpha(a_alpha, a_basalVel,  a_C, a_coords,  a_dit, a_box);

  //a * hydrostatic pressure - generalize?
  FArrayBox limit(a_box,1);
  const FArrayBox& thckOverFlotation = a_coords.getThicknessOverFlotation()[a_dit];
  limit.copy( thckOverFlotation );
  limit *= (m_a * a_coords.iceDensity() * a_coords.gravity());

  FORT_BFRICTIONPLIMIT(CHF_FRA1(a_alpha,0),
		       CHF_CONST_FRA(a_basalVel),
		       CHF_CONST_FRA1(limit,0),
		       CHF_BOX(a_box));
    
}

BasalFrictionRelation* 
BasalFrictionRelation::parseBasalFrictionRelation(const char* a_prefix, int a_recursion)
{
  //can't do more than one recursion without changing the input file syntax
  CH_assert(a_recursion < 2);

  BasalFrictionRelation* basalFrictionRelationPtr = NULL;
  std::string type = "powerLaw";
  ParmParse pp(a_prefix);
  pp.query("basalFrictionRelation",type);

  if (type == "powerLaw")
      {
	ParmParse plPP("BasalFrictionPowerLaw");

	Real m = 1.0;
	plPP.query("m",m);
	bool includeEffectivePressure = false;
	plPP.query("includeEffectivePressure",includeEffectivePressure);
	BasalFrictionPowerLaw*  pl = new BasalFrictionPowerLaw(m,includeEffectivePressure);
	basalFrictionRelationPtr = static_cast<BasalFrictionRelation*>(pl);
      }
  else if (type == "pressureLimitedLaw")
    {
     
      std::string prefix("BasalFrictionPressureLimitedLaw");
      ParmParse ppPLL(prefix.c_str());
      Real a = 1.0;
      ppPLL.query("coefficient",a);

      BasalFrictionRelation* innerLaw = parseBasalFrictionRelation(prefix.c_str(), a_recursion + 1);
      PressureLimitedBasalFrictionRelation* p = 
	new PressureLimitedBasalFrictionRelation(a, innerLaw);
      basalFrictionRelationPtr = static_cast<BasalFrictionRelation*>(p);
    }
  else
    {
      MayDay::Error("undefined basalFrictionRelation in inputs");
    }
  return basalFrictionRelationPtr;
}


#include "NamespaceFooter.H"
