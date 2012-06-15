#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _PWLFLUX_H_
#define _PWLFLUX_H_
#include "SurfaceFlux.H"

#include "NamespaceHeader.H"

/// surface flux which has a piecwise linear dependence
/// on thickness H. 
class PiecewiseLinearFlux : public SurfaceFlux
{
public:
  //constructor
  PiecewiseLinearFlux(const Vector<Real>& a_abscissae, 
		      const Vector<Real>& a_ordinates);

  /// destructor
  ~PiecewiseLinearFlux(){};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep. flux should be defined in meters/second in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                     const LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt);
  


protected:
  Vector<Real> m_abscissae;
  Vector<Real> m_ordinates;

private:
   // Disallowed for all the usual reasons
  void operator=(const PiecewiseLinearFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PiecewiseLinearFlux(const PiecewiseLinearFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  PiecewiseLinearFlux()
  {
    MayDay::Error("invalid operator");
  }

};
#include "NamespaceFooter.H"
#endif
