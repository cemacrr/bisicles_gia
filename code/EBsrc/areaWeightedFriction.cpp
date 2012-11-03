#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "areaWeightedFriction.H"
#include "IceConstants.H"
#include "CONSTANTS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"


areaWeightedFriction::areaWeightedFriction() 
: m_isValSet(false)
{
}

areaWeightedFriction::areaWeightedFriction(Real a_value)
: m_frictionVal(a_value), m_isValSet(true)
{

}

BasalFriction* areaWeightedFriction::new_basalFriction() const
{
  areaWeightedFriction* newPtr = new areaWeightedFriction;
  
  newPtr->m_frictionVal = m_frictionVal;
  newPtr->m_isValSet    = m_isValSet;
  
  return static_cast<BasalFriction*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void areaWeightedFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                        LevelSigmaCS        & a_coordSys,
                                        Real                  a_time,
                                        Real                  a_dt)
{
  CH_assert(m_isValSet);
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_betaSqr[dit].setVal(m_frictionVal);
    }
}

///
void areaWeightedFriction::setFrictionVal(const Real& a_betaSqr) 
{
  m_frictionVal = a_betaSqr; 
  m_isValSet    = true;
}


#include "NamespaceFooter.H"
