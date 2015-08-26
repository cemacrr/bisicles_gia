#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IceThermodynamics.H"
#include "IceThermodynamicsF_F.H"
#include "IceConstants.H"
#include "NamespaceHeader.H"

//default values for static data members
Real IceThermodynamics::s_ice_conductivity = ICECONDUCTIVITY;
Real IceThermodynamics::s_ice_heat_capacity = ICEHEATCAPACITY;
Real IceThermodynamics::s_latent_heat_fusion = ICELATENTHEAT;
Real IceThermodynamics::s_water_conductivity = 1.0e-2 * ICECONDUCTIVITY;
Real IceThermodynamics::s_pressure_melt_factor = ICEPMELTFACTOR;
Real IceThermodynamics::s_triple_point = TRIPLEPOINT;

///compose internal energy F(T,w) from temperature T and water fraction  
void IceThermodynamics::composeInternalEnergy
(FArrayBox& a_F, const FArrayBox& a_T, const FArrayBox& a_W, 
 const Box& a_box, bool a_test)
{
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_T.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_W.nComp() == a_T.nComp());
  CH_assert(a_F.nComp() == a_T.nComp());

  FORT_COMPOSEINTERNALENERGYICE(CHF_FRA(a_F), 
				CHF_CONST_FRA(a_T), 
				CHF_CONST_FRA(a_W), 
				CHF_CONST_REAL(s_ice_heat_capacity),
				CHF_CONST_REAL(s_latent_heat_fusion),
				CHF_BOX(a_box));

  if (a_test)
    {
      CH_assert(a_F.norm(a_box, 0, 0, a_F.nComp() ) <= 2.0* triplepoint * iceheatcapacity);
    } 

}

///decompose internal energy F into temperarure T and water fraction W, given pressure P
void IceThermodynamics::decomposeInternalEnergy
(FArrayBox& a_T, FArrayBox& a_W, 
 const FArrayBox& a_F, const FArrayBox& a_P, const Box& a_box)
{
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_T.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_P.box().contains(a_box));
  CH_assert(a_W.nComp() == a_T.nComp());
  CH_assert(a_F.nComp() == a_T.nComp());
  CH_assert(a_P.nComp() == a_T.nComp());
  
  FORT_DECOMPOSEINTERNALENERGYICE(CHF_FRA(a_T), 
				  CHF_FRA(a_W), 
				  CHF_CONST_FRA(a_F),
				  CHF_CONST_FRA(a_P),
				  CHF_CONST_REAL(s_ice_heat_capacity),
				  CHF_CONST_REAL(s_latent_heat_fusion),
				  CHF_CONST_REAL(s_pressure_melt_factor),
				  CHF_CONST_REAL(s_triple_point),
				  CHF_BOX(a_box));


}
#include "NamespaceFooter.H"
