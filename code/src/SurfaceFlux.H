#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Virtual base class encapsulating surface fluxes for ice-sheet model

#ifndef _SURFACEFLUX_H_
#define _SURFACEFLUX_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

/// interface class encapsulating surface fluxes for ice-sheet model
/** This is a virtual base class -- it is expected that specific cases
    will be derived from this
*/
class SurfaceFlux
{
public:
  /// constructor
  SurfaceFlux() {};

  /// destructor
  virtual ~SurfaceFlux() {};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux() = 0;

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep.  flux should be defined in meters/second in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                    LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt) = 0;
  
  //assemble a SurfaceFlux* object from ParmParse input
  static SurfaceFlux* parseSurfaceFlux(const char* a_prefix);


private:
  // Disallowed for all the usual reasons
  void operator=(const SurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  SurfaceFlux(const SurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


// simple version which sets surface flux to zero
class zeroFlux : public SurfaceFlux
{
public:
  /// constructor
  zeroFlux() {};

  /// destructor
  virtual ~zeroFlux() {};

  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep.  flux should be defined in meters/second in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                    LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt);


private:
  // Disallowed for all the usual reasons
  void operator=(const zeroFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  zeroFlux(const zeroFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


// simple version which sets surface flux to a constant in space and time
class constantFlux : public SurfaceFlux
{
public:
  /// constructor
  constantFlux();

  /// destructor
  virtual ~constantFlux() {}

  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep. flux should be defined in meters/second in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                    LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt);

  /// set flux value in meters/year
  void setFluxVal(const Real& a_fluxVal);

protected:
  Real m_fluxVal;

  bool m_isValSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const constantFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  constantFlux(const constantFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


/// class which takes an input fortran array 
/** averages or interpolates as necessary to fill the flux
 */
class fortranInterfaceFlux : public SurfaceFlux
{
public:
  /// constructor
  fortranInterfaceFlux();

  /// destructor
  virtual ~fortranInterfaceFlux() {}

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
                                    LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt);

  /// set fortran array-valued surface flux
  void setFluxVal(Real* a_data_ptr,
                  const int* a_dimInfo,
                  const Real* a_dew, const Real* a_dns);

protected:
  FArrayBox m_fluxVal;
  RealVect m_fluxDx;

  bool m_isValSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const fortranInterfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  fortranInterfaceFlux(const fortranInterfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};



// class which applies a different surface flux 
// depending on the mask
class MaskedFlux : public SurfaceFlux
{

public:
  /// constructor
  MaskedFlux(SurfaceFlux* a_groundedIceFlux, SurfaceFlux* a_floatingIceFlux);

  /// destructor
  virtual ~MaskedFlux(){};

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
                                    LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt);

protected:
  SurfaceFlux* m_groundedIceFlux;
  SurfaceFlux* m_floatingIceFlux;

private:
  // Disallowed for all the usual reasons
  void operator=(const MaskedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MaskedFlux(const MaskedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  MaskedFlux()
  {
    MayDay::Error("invalid operator");
  }

};
  



#include "NamespaceFooter.H"


#endif // multiple-include preventer
