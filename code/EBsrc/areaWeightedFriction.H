#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _AREAWEIGHTEDFRICTION_H_
#define _AREAWEIGHTEDFRICTION_H_

#include "BasalFriction.H"
#ifdef CH_USE_EB
#include "EBIndexSpace.H"
#endif
#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

class areaWeightedFriction : public BasalFriction
{
public:
  /// constructor
  areaWeightedFriction();

  areaWeightedFriction(BasalFriction* frictionPtr);

  /// destructor
  virtual ~areaWeightedFriction() {}

  /// factory method
  /** return a pointer to a new BasalFriction object
   */
  virtual BasalFriction* new_basalFriction() const;

  /// define basal friction coefficient beta^2 and place in a_betaSqr
  /** time and dt are included in case this is time-dependent. Units 
      should be Pa*a/m (any conversion to mks units is internal to the
      AmrIce code)       
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                LevelSigmaCS        & a_coordSys,
                                Real                  a_time,
                                Real                  a_dt);

protected:
  BasalFriction* m_frictionPtr;
  
public:
  // Disallowed for all the usual reasons
  void operator=(const areaWeightedFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  areaWeightedFriction(const areaWeightedFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

#ifndef CH_USE_EB
#if CH_SPACEDIM == 1
  void setGroundingLineData(const IntVect                & a_groundingLineIv,
                            const IndexTM<Real,SpaceDim> & a_physCoordGroundingPt,
                            const Real                   & a_lengthFraction);

  IntVect                m_groundingLineIv;
  IndexTM<Real,SpaceDim> m_physCoordGroundingPt;
  Real                   m_lengthFraction;
#endif
#endif
};

#include "NamespaceFooter.H"

#endif 
