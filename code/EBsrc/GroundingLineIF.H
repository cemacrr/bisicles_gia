#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GROUNDINGLINEIF_H_
#define _GROUNDINGLINEIF_H_

#include "MayDay.H"
#include "IndexTM.H"
#include "REAL.H"
#include "LevelData.H"
#include "BaseIF.H"
#include "DataFileIF.H"
#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a plane templated in dimension.
 */
class GroundingLineIF: public BaseIF
{
public:
  ///
  /**
     Constructor specifying the coefficients of the plane normal, a_normal,
     and a point on the plane, a_point.
    */
  GroundingLineIF(const Real                  & a_rhoIce  ,
                  const Real                  & a_rhoWater,
                  const Real                  & a_sealevel,
                  const RefCountedPtr <BaseIF>& a_topographyIF);

  /// Copy constructor
  GroundingLineIF(GroundingLineIF& a_inputIF);

  /// Destructor
  virtual ~GroundingLineIF();
  /**
      Return the value of the function at a_point.
   */

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual IndexTM<Real,GLOBALDIM> normal(const IndexTM<Real,GLOBALDIM>& a_point) const ;

  virtual Vector<IndexTM<Real,GLOBALDIM> > gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual void setDepth(const RefCountedPtr<DataFileIF>& m_depth);

  virtual BaseIF* newImplicitFunction() const;

  //protected:

  // sea level
  Real m_seaLevel;
  Real m_rhoWater;
  Real m_rhoIce  ;
  
  const RefCountedPtr<BaseIF> m_topographyIF;
  RefCountedPtr  <DataFileIF> m_depthIF;

private:
  GroundingLineIF()
  {
    MayDay::Abort("GroundingLineIF uses strong construction");
  }

  void operator=(const GroundingLineIF& a_inputIF)
  {
    MayDay::Abort("GroundingLineIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
