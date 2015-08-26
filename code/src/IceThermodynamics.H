#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef ICE_THERMODYNAMICS_H
#define ICE_THERMODYNAMICS_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "NamespaceHeader.H"

/// common functions for ice thermodyamics. All public static for now so could
/// just be a namespace.
class IceThermodynamics
{
public:
  static Real s_ice_conductivity, 
    s_ice_heat_capacity, 
    s_latent_heat_fusion,
    s_water_conductivity, 
    s_pressure_melt_factor,
    s_triple_point;

  ///compose internal energy F(T,w) from temperature T and water fraction  
  static void composeInternalEnergy(FArrayBox& a_F, const FArrayBox& a_T, const FArrayBox& a_w, 
				    const Box& a_box, bool a_test = true);
    
  ///compose internal energy F(T,w) from temperature T assuming water fraction is zero  
  static void composeInternalEnergy(FArrayBox& a_F, const FArrayBox& a_T,  
				    const Box& a_box, bool a_test = true)
  {
    FArrayBox w(a_box,a_T.nComp()); w.setVal(0.0);
    composeInternalEnergy(a_F, a_T, w, a_box, a_test);
  }


  ///decompose internal energy F into temperarure T and water fraction W, given pressure P
  static void decomposeInternalEnergy(FArrayBox& a_T,  FArrayBox& a_w, 
			       const FArrayBox& a_F, const FArrayBox& a_P, 
			       const Box& a_box);
  
  

};


#include "NamespaceFooter.H"
#endif 
