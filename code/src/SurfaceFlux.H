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
#ifdef HAVE_PYTHON
#include "Python.h"
#endif

#include "NamespaceHeader.H"

/// interface class encapsulating surface fluxes for ice-sheet model
/** This is a virtual base class -- it is expected that specific cases
    will be derived from this
*/

class AmrIce;

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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt) = 0;
  
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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);


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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);

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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);

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
  MaskedFlux(SurfaceFlux* a_groundedIceFlux, SurfaceFlux* a_floatingIceFlux,
	     SurfaceFlux* a_openSeaFlux, SurfaceFlux* a_openLandFlux);

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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);


protected:
  SurfaceFlux* m_groundedIceFlux;
  SurfaceFlux* m_floatingIceFlux;
  SurfaceFlux* m_openSeaFlux;
  SurfaceFlux* m_openLandFlux;
  
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
  

//class which applies a flux inside a box with Real boundaries 
//(and optionally time boundaries) and zero outside
class BoxBoundedFlux : public SurfaceFlux
{
  RealVect m_lo,m_hi;
  Real m_startTime, m_endTime;
  SurfaceFlux* m_fluxPtr;
public:
  BoxBoundedFlux
(const RealVect& a_lo,
 const RealVect& a_hi,
 const Real& a_startTime,
 const Real& a_endTime,
 SurfaceFlux* a_fluxPtr)
    : m_lo(a_lo),m_hi(a_hi),m_startTime(a_startTime),m_endTime(a_endTime)
  {
    m_fluxPtr = a_fluxPtr->new_surfaceFlux();
    if (!(m_hi > m_lo))
      {
	pout() <<  "!(m_hi > m_lo), m_lo = " << m_lo << " m_hi = " << m_hi << endl;
	CH_assert( (m_hi > m_lo) );
	MayDay::Error("BoxBoundedFlux::  !(m_hi > m_lo)");
      }
  }
  

   /// destructor
  virtual ~BoxBoundedFlux(){};

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
				  const AmrIce& a_amrIce, 
				  int a_level, Real a_dt);

private:
  // Disallowed for all the usual reasons
  void operator=(const BoxBoundedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  BoxBoundedFlux(const BoxBoundedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  BoxBoundedFlux()
  {
    MayDay::Error("invalid operator");
  }

};

//class which contains two fluxes (x any y)
//and two scalars (a and b) and provide the
//flux a x + b y
class AxbyFlux : public SurfaceFlux
{
  SurfaceFlux* m_x;
  SurfaceFlux* m_y;
  Real m_a,m_b;
  
public:
  AxbyFlux(const Real& a_a, SurfaceFlux* a_x, 
	   const Real& a_b, SurfaceFlux* a_y);

  virtual ~AxbyFlux();

  virtual SurfaceFlux* new_surfaceFlux();
  
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);

private:
  
  // Disallowed for all the usual reasons
  void operator=(const AxbyFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AxbyFlux(const AxbyFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  AxbyFlux()
  {
    MayDay::Error("invalid operator");
  }

};


//class which adds together the results of a list of fluxes
class CompositeFlux : public SurfaceFlux
{
  
  Vector<SurfaceFlux*> m_fluxes;

public:
  CompositeFlux(const Vector<SurfaceFlux*>& a_fluxes);

  virtual ~CompositeFlux();
    
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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);

  
private:
  // Disallowed for all the usual reasons
  void operator=(const CompositeFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  CompositeFlux(const CompositeFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  CompositeFlux()
  {
    MayDay::Error("invalid operator");
  }


};

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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);
  


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




#ifdef HAVE_PYTHON

///Use a python function to evaluate the surface flux
/** the user needs to provide a python function that
    takes x,y,thickness,topography as args
    returns a single scalar (the flux)
    def <fluxf>(x,y,thickness,topography):
    ...
    return flux

 */
class PythonSurfaceFlux : public SurfaceFlux
{
  PyObject*  m_pModule;
  PyObject*  m_pFunc;
 
public:
  PythonSurfaceFlux(const std::string& a_pyModule,
		    const std::string& a_pyFunc);

  virtual ~PythonSurfaceFlux();
    
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
				    const AmrIce& a_amrIce, 
				    int a_level, Real a_dt);

  
private:
  // Disallowed for all the usual reasons
  void operator=(const CompositeFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PythonSurfaceFlux(const CompositeFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  PythonSurfaceFlux()
  {
    MayDay::Error("invalid operator");
  }

  PythonSurfaceFlux(PyObject*  a_pModule,  PyObject*  a_pFunc)
  {
    m_pModule = a_pModule;
    m_pFunc = a_pFunc;
    Py_XINCREF(m_pModule);
    Py_XINCREF(m_pFunc);
  }
  

};

#endif



#include "NamespaceFooter.H"


#endif // multiple-include preventer
