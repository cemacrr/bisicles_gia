#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BasalFriction.H"
#include "IceConstants.H"
#include "CONSTANTS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

// -----------------------------------------------------------------
//   zeroFriction
// -----------------------------------------------------------------

/// factory method
/** return a pointer to a new BasalFriction object
 */
BasalFriction* 
zeroFriction::new_basalFriction() const
{
  zeroFriction* newPtr = new zeroFriction;
  return static_cast<BasalFriction*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
zeroFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                               LevelSigmaCS& a_coordSys,
                               Real a_time,
                               Real a_dt)
{
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_betaSqr[dit].setVal(0.0);
    }
}

// -----------------------------------------------------------------
//   constantFriction
// -----------------------------------------------------------------

constantFriction::constantFriction() : m_isValSet(false)
{
}

constantFriction::constantFriction(Real a_value)
  : m_frictionVal(a_value), m_isValSet(true)
{

}

BasalFriction* 
constantFriction::new_basalFriction() const
{
  constantFriction* newPtr = new constantFriction;
  newPtr->m_frictionVal = m_frictionVal;
  newPtr->m_isValSet = m_isValSet;
  return static_cast<BasalFriction*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
constantFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
				   LevelSigmaCS& a_coordSys,
				   Real a_time,
				   Real a_dt)
{
  CH_assert(m_isValSet);
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_betaSqr[dit].setVal(m_frictionVal);
    }
}


///
void
constantFriction::setFrictionVal(const Real& a_betaSqr) 
{
  m_frictionVal = a_betaSqr; 
  m_isValSet = true;
}


// -----------------------------------------------------------------
//   sinusoidalFriction
// -----------------------------------------------------------------



/// constructor
sinusoidalFriction::sinusoidalFriction()
{
}

sinusoidalFriction::sinusoidalFriction(const Real& a_betaVal,
                                       const RealVect& a_omega,
                                       const Real& a_eps,
                                       const RealVect& a_domainSize) 
  : m_betaVal(a_betaVal), m_omega(a_omega), m_eps(a_eps), m_domainSize(a_domainSize)
{

}
                                                 

/// destructor
sinusoidalFriction::~sinusoidalFriction()
{
}

/// factory method
/** return a pointer to a new BasalFriction object
 */
BasalFriction* 
sinusoidalFriction::new_basalFriction() const
{
  sinusoidalFriction* newPtr = new sinusoidalFriction;
  newPtr->m_betaVal = m_betaVal;
  newPtr->m_omega = m_omega;
  newPtr->m_eps = m_eps;
  newPtr->m_domainSize = m_domainSize;
  return static_cast<BasalFriction*>(newPtr);
}

/// define source term for thickness evolution and place it in friction
/** dt is included in case one needs integrals or averages over a
    timestep. friction should be defined in meters/second in the current 
    implementation. 
*/
void
sinusoidalFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
				     LevelSigmaCS& a_coordSys,
				     Real a_time,
				     Real a_dt)
{
  RealVect twoPiOmega(m_omega);
  twoPiOmega /= m_domainSize;
  twoPiOmega *= 2*Pi;

  RealVect dx = a_coordSys.dx();

  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisBeta = a_betaSqr[dit];

      BoxIterator bit(thisBeta.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;
          
          Real betaVal = m_betaVal*(1.0+m_eps+D_TERM(cos(twoPiOmega[0]*x[0]),
                                                     *cos(twoPiOmega[1]*x[1]),
                                                     *cos(twoPiOmega[2]*x[2])));
          thisBeta(iv,0) = betaVal;
        } // end loop over cells
    } // end loop over grids

}

/// set friction value in Pa*a/m)
void
sinusoidalFriction::setSinParameters(const Real& a_betaVal, 
                                     const RealVect& a_omega,
                                     const Real& a_eps,
                                     const RealVect& a_domainSize)
{
  m_betaVal = a_betaVal;
  m_omega = a_omega;
  m_eps = a_eps;
  m_domainSize = a_domainSize;
}

#include "NamespaceFooter.H"
