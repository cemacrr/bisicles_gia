#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifdef HAVE_PYTHON
#ifndef _PYTHONINTERFACE_H_
#define _PYTHONINTERFACE_H_

#include "Python.h"
#include "IceThicknessIBC.H"
#include "IceTemperatureIBC.H"
#include "SurfaceFlux.H"
#include "BasalFriction.H"
#include "MuCoefficient.H"
#include "VelocityBC.H"
#include "NamespaceHeader.H"
//
// PythonInterface.H
// ============
// BISICLES embedded Python interface. Includes
//
// PythonInterface, a namespace for common parts of the interface
//
// PythonInterface::PythonIBC, a PhysIBC-derived class which allows the bedrock and topography to be set by
// a user-defined python function f(x,y) . Reflection boundary conditions are imposed on
// each edge. 
//
// PythonInterface::PythonTemperatureIBC, much as PythonIBC, but provides initial / boundary conditions for the temperature
// field
//
// PythonInterface::PythonSurfaceFlux , a SurfaceFlux derived class which allows a thickness source to be set
// by a python function f(x,y,thickness,topography)

class AmrIce;

namespace PythonInterface
{
  void InitializePythonModule(PyObject **a_pModule,
			      const std::string& a_pyModuleName);
  
  void InitializePythonFunction(PyObject **a_pFunc,
				PyObject *a_pModule,
				const std::string& a_pyFuncName);
  
  void PythonEval(PyObject* a_pyFunc, Vector<Real>& a_value, Vector<Real>& a_arg); 

  /// PythonInterface::PythonIBC, a PhysIBC-derived class which allows the bedrock and topography to be set by  a user-defined python function f(x,y) . 
  /** Reflection boundary conditions are imposed on
      each edge. 
  */
  class PythonIBC : public IceThicknessIBC
  {
    PyObject*  m_pModule;
    PyObject*  m_pFuncThck;
    PyObject*  m_pFuncTopg;
    PyObject*  m_pFuncRhs;
    PyObject*  m_pFuncFaceVel;
    
    enum BCType {IceDivide, NoSlip, MaxBCType};
    
    RefCountedPtr<BCFunction> BCFactory(BCType a_type, 
					int a_dir, 
					Side::LoHiSide a_side)
    {
      if (a_type == IceDivide)
	{
	  return RefCountedPtr<BCFunction>(new IceDivideBCFunction(a_dir, a_side));
	}
      else if (a_type == NoSlip)
	{
	  return RefCountedPtr<BCFunction>(new NoSlipBCFunction(a_dir, a_side));
	}
      
      CH_assert(a_type < MaxBCType && a_type > 0);

    }

  public:
    
    PythonIBC(const std::string& a_pyModuleName , 
	      const std::string& a_pyFuncThckName, 
	      const std::string& a_pyFuncTopgName, 
	      const std::string& a_pyFuncRhsName,
	      const std::string& a_pyFuncFaceVelName);
    
    PythonIBC();
    
    virtual ~PythonIBC();
    
    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object.
    */
    virtual void define(const ProblemDomain& a_domain,
			const Real&          a_dx);
    
    /// Factory method - this object is its own factory
    /**
       Return a pointer to a new IceThicknessIBC object with 
       m_isDefined = false (i.e., its define() must be called before it
       is used). 
    */
    virtual IceThicknessIBC* new_thicknessIBC();
    
    
    /// Set up initial conditions
    /**
     */
    virtual void initialize(LevelData<FArrayBox>& a_U);
    
    void regridIceGeometry(LevelSigmaCS& a_coords,
			   const RealVect& a_dx,
			   const RealVect& a_domainSize,
			   const Real& a_time,
			   const LevelSigmaCS* a_crseCoords,
			   const int a_refRatio);
    
    /// set up initial ice state
    /** reads info from ParmParse and sets up ice sheet geometry
     */
    virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
				       const RealVect& a_dx,
				       const RealVect& a_domainSize,
				       const Real& a_time,
				       const LevelSigmaCS* a_crseCoords,
				       const int a_refRatio);

  /// if appropriate, modify velocity solve RHS in a problem-dependent way. 
  /** default is to do nothing
   */  
    virtual void modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
				   LevelSigmaCS& a_coords,
				   const ProblemDomain& a_domain,
				   Real a_time, Real a_dt);
    /// if appropriate, modify face velocities in a problem-dependent way. 
    /** default is to do nothing
     */
    virtual void modifyFaceVelocity(LevelData<FluxBox>& a_faceVel,
				    const LevelSigmaCS& a_coords,
				    const ProblemDomain& a_domain) const ;

    /// Set boundary fluxes
    /**
     */
    virtual void primBC(FArrayBox&            a_WGdnv,
			const FArrayBox&      a_Wextrap,
			const FArrayBox&      a_W,
			const int&            a_dir,
			const Side::LoHiSide& a_side,
			const Real&           a_time);

    /// Set boundary slopes
    /**
       The boundary slopes in a_dW are already set to one sided difference
       approximations.  If this function doesn't change them they will be
       used for the slopes at the boundaries.
    */
    virtual
    void setBdrySlopes(FArrayBox&       a_dW,
		       const FArrayBox& a_W,
		       const int&       a_dir,
		       const Real&      a_time);

    /// Adjust boundary fluxes to account for artificial viscosity
    /**
     */
    virtual
    void artViscBC(FArrayBox&       a_F,
		   const FArrayBox& a_U,
		   const FArrayBox& a_divVel,
		   const int&       a_dir,
		   const Real&      a_time);

    /// return boundary condition for Ice velocity solve
    /** 
     */
    virtual BCHolder velocitySolveBC();

    


    /// set non-periodic ghost cells for surface height z_s. 
    virtual void setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
				     LevelSigmaCS& a_coords,
				     const ProblemDomain& a_domain,
				     const RealVect& a_dx, 
				     Real a_time, Real a_dt);
    /// set non-periodic geometry (thickness, topography) ghost cells
    virtual void setGeometryBCs(LevelSigmaCS& a_coords,
				const ProblemDomain& a_domain,
				const RealVect& a_dx,
				Real a_time, Real a_dt);

  private:

    RealVect m_dx;
    BCFunc m_velBCs;
    RefCountedPtr<BCFunction> m_velBC;
    bool m_isBCsetUp;
    bool m_verbose;
    void setupBCs();

  private:
    // Disallowed for all the usual reasons
    void operator=(const PythonIBC& a_input)
    {
      MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    PythonIBC(const PythonIBC& a_input)
    {
      MayDay::Error("invalid operator");
    }

    PythonIBC(PyObject*  a_pModule,  PyObject*  a_pFuncThck, 
	      PyObject*  a_pFuncTopg, PyObject*  a_pFuncRhs,
	      PyObject*  a_pFuncFaceVel)
    {
      m_isBCsetUp = false;
      m_pModule = a_pModule;
      m_pFuncThck = a_pFuncThck;
      m_pFuncTopg = a_pFuncTopg;
      m_pFuncRhs = a_pFuncRhs;
      m_pFuncFaceVel = a_pFuncFaceVel;
      Py_XINCREF(m_pModule);
      Py_XINCREF(m_pFuncThck);
      Py_XINCREF(m_pFuncTopg);
      if (m_pFuncRhs != NULL)
	{
	  Py_XINCREF(m_pFuncRhs);
	}
      if (m_pFuncFaceVel != NULL)
	{
	  Py_XINCREF(m_pFuncFaceVel);
	}
    }
    
  };
  /// Use a python function to evaluate the surface flux
  /** the user needs to provide a python function foo that
      takes x,y,thickness,topography as args
      and returns a single scalar (the flux)
      def foo(x,y,thickness,topography):
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
  /// PhysIBC-derived class that allows the temperature to set a a user-defined python function f(x,y,sigma,thickness,topography,surface) . 
  /** Reflection boundary conditions are imposed on
      each edge. 
  */
  class PythonIceTemperatureIBC : public IceTemperatureIBC
  {
  public:
    virtual ~PythonIceTemperatureIBC();
    PythonIceTemperatureIBC(const std::string& a_pyModuleName,
			    const std::string& a_pyFuncName) ;

    virtual PythonIceTemperatureIBC* new_temperatureIBC();
    
    /// set basal heat flux
    virtual void basalHeatFlux(LevelData<FArrayBox>& a_flux,
			       const AmrIce& a_amrIce, 
			       int a_level, Real a_dt);

#if BISICLES_Z == BISICLES_LAYERED
    virtual void initializeIceTemperature(LevelData<FArrayBox>& a_T, 
					  LevelData<FArrayBox>& a_surfaceT, 
					  LevelData<FArrayBox>& a_basalT,
					  const LevelSigmaCS& a_coordSys);

    virtual void setIceTemperatureBC(LevelData<FArrayBox>& a_T, 
				    LevelData<FArrayBox>& a_surfaceT, 
				    LevelData<FArrayBox>& a_basalT,
				    const LevelSigmaCS& a_coordSys);

#elif BISICLES_Z == BISICLES_FULLZ
    virtual void initializeIceTemperature(LevelData<FArrayBox>& a_T,
					  const LevelSigmaCS& a_coordSys);
    virtual void setIceTemperatureBC(LevelData<FArrayBox>& a_T,
					  const LevelSigmaCS& a_coordSys);
#endif

private:
    PyObject*  m_pModule;
    PyObject*  m_pFunc;
    PythonIceTemperatureIBC() {;} 

    PythonIceTemperatureIBC(PyObject*  a_pModule,  PyObject*  a_pFunc)
    {
      m_pModule = a_pModule;
      m_pFunc = a_pFunc;
      Py_XINCREF(m_pModule);
      Py_XINCREF(m_pFunc);
    }

  };



  class PythonBasalFriction : public BasalFriction
  {

    PyObject*  m_pModule;
    PyObject*  m_pFunc;
 
  public:
    PythonBasalFriction(const std::string& a_pyModule,
			const std::string& a_pyFunc);

    virtual ~PythonBasalFriction();

    virtual BasalFriction* new_basalFriction() const;
    
    virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
				  LevelSigmaCS& a_coordSys,
				  Real a_time,
				  Real a_dt);
  private:
    void operator=(const PythonBasalFriction& a_input)
    {
      MayDay::Error("invalid operator");
    }
    
    PythonBasalFriction(const PythonBasalFriction& a_input)
    {
      MayDay::Error("invalid operator");
    }

    PythonBasalFriction()
    {

    }

    
    PythonBasalFriction(PyObject*  a_pModule,  PyObject*  a_pFunc)
    {
      m_pModule = a_pModule;
      m_pFunc = a_pFunc;
      Py_XINCREF(m_pModule);
      Py_XINCREF(m_pFunc);
    }

  };


  class PythonMuCoefficient : public MuCoefficient
  {

    PyObject*  m_pModule;
    PyObject*  m_pFunc;
 
  public:
    PythonMuCoefficient(const std::string& a_pyModule,
			const std::string& a_pyFunc);
    
    virtual ~PythonMuCoefficient();

    virtual MuCoefficient* new_muCoefficient() const;
    
    virtual void setMuCoefficient(LevelData<FArrayBox>& a_C,
				  LevelSigmaCS& a_coordSys,
				  Real a_time,
				  Real a_dt);
  private:
    void operator=(const PythonMuCoefficient& a_input)
    {
      MayDay::Error("invalid operator");
    }
    
    PythonMuCoefficient(const PythonMuCoefficient& a_input)
    {
      MayDay::Error("invalid operator");
    }

    PythonMuCoefficient()
    {

    }

    
    PythonMuCoefficient(PyObject*  a_pModule,  PyObject*  a_pFunc)
    {
      m_pModule = a_pModule;
      m_pFunc = a_pFunc;
      Py_XINCREF(m_pModule);
      Py_XINCREF(m_pFunc);
    }

  };

}

#include "NamespaceFooter.H"
#endif
#endif
