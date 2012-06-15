#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BisiclesF_F.H"
#include "PiecewiseLinearFlux.H"

PiecewiseLinearFlux::PiecewiseLinearFlux(const Vector<Real>& a_abscissae, 
					 const Vector<Real>& a_ordinates)
  :m_abscissae(a_abscissae),m_ordinates(a_ordinates)
{
  CH_assert(m_abscissae.size() == m_ordinates.size());
}


SurfaceFlux* PiecewiseLinearFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new PiecewiseLinearFlux(m_abscissae,m_ordinates));
}

void PiecewiseLinearFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					       const LevelSigmaCS& a_coordSys,
					       Real a_time,
					       Real a_dt)
{
  Vector<Real> dx(m_abscissae.size());
  Vector<Real> db(m_abscissae.size());

  for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
    {
      FORT_PWLFILL(CHF_FRA1(a_flux[dit],0),
		   CHF_CONST_FRA1(a_coordSys.getH()[dit],0),
		   CHF_CONST_VR(m_abscissae),
		   CHF_CONST_VR(m_ordinates),
		   CHF_VR(dx),CHF_VR(db),
		   CHF_BOX(a_flux[dit].box()));
    }

}
