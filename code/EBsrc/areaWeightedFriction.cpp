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


#ifdef CH_USE_EB
#include "EBISBox.H"
#include "EBIndexSpace.H"
#include "VoFIterator.H"
#endif

#include "NamespaceHeader.H"


areaWeightedFriction::areaWeightedFriction() 
  :m_frictionPtr(NULL)
{
}

areaWeightedFriction::areaWeightedFriction(BasalFriction* a_frictionPtr)
  : m_frictionPtr(a_frictionPtr)
{
}

BasalFriction* areaWeightedFriction::new_basalFriction() const
{
  areaWeightedFriction* newPtr = new areaWeightedFriction;
  
  newPtr->m_frictionPtr = m_frictionPtr;
   
  return static_cast<BasalFriction*>(newPtr);
}

  /// set the friction coefficient as grounded-area-Fraction*grounded value.
  /** set the friction coefficient as grounded-area-Fraction*grounded value.
  */
void areaWeightedFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                            LevelSigmaCS        & a_coordSys,
                                            Real                  a_time,
                                            Real                  a_dt)
{
  CH_assert(m_frictionPtr != NULL);
  m_frictionPtr->setBasalFriction(a_betaSqr,a_coordSys,a_time,a_dt);

#ifdef CH_USE_EB

  // ebisPtr
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  // if there are multiple levels, we need the problemDomain for the finest level
  int validForSingleLevelCalculationOnly = 0;
  ProblemDomain domain = ebisPtr->getBox(validForSingleLevelCalculationOnly);
 
  // two extra ghost cells because Dan G. recommends extra ghost cells in the ebisl
  EBISLayout ebisl;
  int nGhost = 2;

  // dbl
  DisjointBoxLayout dbl = a_betaSqr.disjointBoxLayout();

  // fill ebisl
  ebisPtr->fillEBISLayout(ebisl, dbl, domain, nGhost);

  //area weight the friction value
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // geometry data and (trivial) connectivity for this box
      EBISBox ebisBox = ebisl[dit()];
      
      // modify the value in the irreg ivs
      Box box = dbl[dit()];
          
      // iterate over all irregular volumes of fluid
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          IntVect iv        = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          if (vofs.size() == 0)
            {
              a_betaSqr[dit()](iv,0) = 0.0;
            }
            else if (vofs.size() == 1)
            {
              Real areaFraction = ebisBox.volFrac(vofs[0]);
              a_betaSqr[dit()](iv,0) *= areaFraction;
            }
            else
            {
              MayDay::Error("Not coded for multi-cells");
            }
        }

    }

#elif CH_SPACEDIM == 1
  {
    // dbl
    DisjointBoxLayout dbl = a_betaSqr.disjointBoxLayout();

    DataIterator dit = a_betaSqr.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        Box box = dbl[dit()];
        
        // iterate over all irregular volumes of fluid
        for (BoxIterator bit(box); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            if (iv[0] == m_groundingLineIv[0])
              {
                a_betaSqr[dit()](iv,0) *= m_lengthFraction;
              }
            else if (iv[0] > m_groundingLineIv[0])
              {
                a_betaSqr[dit()](iv,0) = 0.0;
              }
          }
      }
  }
#endif
    
}

#ifndef CH_USE_EB
#if CH_SPACEDIM == 1
void areaWeightedFriction::setGroundingLineData(const IntVect                & a_groundingLineIv,
                                                const IndexTM<Real,SpaceDim> & a_physCoordGroundingPt,
                                                const Real                   & a_lengthFraction)

{
  m_groundingLineIv      = a_groundingLineIv;
  m_physCoordGroundingPt = a_physCoordGroundingPt;
  m_lengthFraction       = a_lengthFraction;   
}        
#endif
#endif

#include "NamespaceFooter.H"
