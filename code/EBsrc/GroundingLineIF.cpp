#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )

#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>

#include "BaseIF.H"
#include "DataFileIF.H"
#include "GroundingLineIF.H"
#include "CH_Attach.H"
#include "NamespaceHeader.H"

GroundingLineIF::GroundingLineIF(const Real                 & a_rhoIce  ,
                                 const Real                 & a_rhoWater,
                                 const Real                 & a_seaLevel,
                                 const RefCountedPtr<BaseIF>& a_topographyIF)
  : m_seaLevel    (a_seaLevel),
    m_rhoWater    (a_rhoWater),
    m_rhoIce      (a_rhoIce),
    m_topographyIF(a_topographyIF),
    m_depthIFSet  (false)
{
}

GroundingLineIF::GroundingLineIF(const GroundingLineIF& a_inputIF)

  // Remember the parameters
  : m_seaLevel    (a_inputIF.m_seaLevel),
    m_rhoWater    (a_inputIF.m_rhoWater),
    m_rhoIce      (a_inputIF.m_rhoIce),
    m_topographyIF(a_inputIF.m_topographyIF),
    m_depthIFSet  (a_inputIF.m_depthIFSet),
    m_depthIF     (a_inputIF.m_depthIF)
{
}

GroundingLineIF::~GroundingLineIF()
{
}

Real GroundingLineIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                            const IndexTM<Real,GLOBALDIM>& a_point) const
{
  CH_assert(m_depthIFSet);
  CH_assert(m_topographyIF!= NULL);
  
  Real thicknessPartial = m_depthIF     ->value(a_partialDerivative,a_point);
  Real zBottomPartial   = m_topographyIF->value(a_partialDerivative,a_point);
  
  Real retval =  m_rhoWater*(m_seaLevel - zBottomPartial) - m_rhoIce*thicknessPartial;
  
  //pout()<<"point = "           <<a_point         <<endl;
  //pout()<<"thicknessPartial = "<<thicknessPartial<<endl;
  //pout()<<"zBottomPartial = "  <<zBottomPartial  <<endl;
  //pout()<<"retval = "          <<retval          <<endl;
  //pout()<<""                   <<""              <<endl;
  
  return retval;
}

Real GroundingLineIF::value(const RealVect & a_point) const
{

  CH_assert(m_depthIFSet);
  CH_assert(m_topographyIF != NULL);
  IndexTM<Real,GLOBALDIM> pt;

  //check does SpaceDim = GLOBALDIM
  if (GLOBALDIM==3 && SpaceDim==2)
    {
      MayDay::Abort("GroundingLineIF should be wrapped in ReferenceHeightIF when GLOBALDIM==3 and SpaceDim==2");
    }
  else
    {
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          pt[idir] = a_point[idir];
        }
    }

  return value(pt);
}

Real GroundingLineIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  CH_assert(m_depthIFSet);
  CH_assert(m_topographyIF != NULL);
  
  IndexTM<int,GLOBALDIM> noDerivative = IndexTM<int,GLOBALDIM>::Zero;
 
  return value(noDerivative,a_point);
}

void GroundingLineIF::setDepth(const RefCountedPtr<DataFileIF>  & a_depth)
{
  m_depthIF    = a_depth;
  m_depthIFSet = true;
}

BaseIF* GroundingLineIF::newImplicitFunction() const
{
  GroundingLineIF* groundingLinePtr = new GroundingLineIF(m_rhoIce,
                                                          m_rhoWater,
                                                          m_seaLevel,
                                                          m_topographyIF);
  
  groundingLinePtr->setDepth(m_depthIF);

  return static_cast<BaseIF*>(groundingLinePtr);
}

#include "NamespaceFooter.H"
