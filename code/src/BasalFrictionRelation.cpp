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
#include "NamespaceHeader.H"


/// computes cell-centered \f$ s = - T_b . u_b \f$
/// based on the cell-centered velocity
/// and a coefficient field beta such that \f T_b  = - \alpha (u, C) u_b \f
/** 
    \param a_C -- \f$C\f$ based on the local velocity field. 
    \param a_basalVel -- Cell-centered basal velocity field.
    \param a_basalTheta: Cell-centered temperature field
    \param a_box: cell-centered box over which to do this computation 
*/
void 
BasalFrictionRelation::computeDissipation(FArrayBox& a_dissipation, 
					  const FArrayBox& a_basalVel, 
					  const FArrayBox& a_thckOverFlotation,
					  const FArrayBox& a_C,
					  const BaseFab<int>& a_mask,
					  const Box& a_box) const
{
    
  computeAlpha(a_dissipation,a_basalVel,a_thckOverFlotation,a_C,a_mask,a_box);
  
  FORT_BFRICTIONAUU(CHF_FRA1(a_dissipation,0),
		    CHF_CONST_FRA(a_basalVel),
		    CHF_BOX(a_box));

}


void
BasalFrictionPowerLaw::computeAlpha(FArrayBox& a_alpha,
				    const FArrayBox& a_basalVel, 
				    const FArrayBox& a_thckOverFlotation,
				    const FArrayBox& a_C,
				    const BaseFab<int>& a_mask,
				    const Box& a_box) const
{
  Real mm1 = m_m - 1.0;
  FArrayBox C(a_C.box(),1);
  C.copy(a_C);

  if (m_includeEffectivePressure)
    {
      C *= a_thckOverFlotation;
      //C *= 300.0;
    }

  if (std::abs(mm1) > 1.0e-3)
    {
      FORT_BFRICTIONPOWER(CHF_FRA1(a_alpha,0),
			  CHF_CONST_FRA( a_basalVel),
			  CHF_CONST_FRA1(C,0),
			  CHF_CONST_FIA1(a_mask,0),
			  CHF_CONST_REAL(mm1),
			  CHF_BOX(a_box));	
    }
  else
    {
      a_alpha.copy(C,a_box);
    } 
  
}

#include "NamespaceFooter.H"
