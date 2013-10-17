#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  LevelDataTemperatureIBC.H
// ============
//
// PhysIBC-derived class which stores initial temperature data
// and imposes either periodic or reflection boundary conditions

#ifndef _LEVELDATATEMPERATUREIBC_H_
#define _LEVELDATATEMPERATUREIBC_H_

#include "IceTemperatureIBC.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   
 */
class LevelDataTemperatureIBC : public IceTemperatureIBC
{  
  RefCountedPtr<LevelData<FArrayBox> > m_temp;
#if BISICLES_Z == BISICLES_LAYERED
  RefCountedPtr<LevelData<FArrayBox> > m_surfaceTemp ;
  RefCountedPtr<LevelData<FArrayBox> > m_basalHeatFlux ;
#endif
  RealVect m_dx;
public:
  /// Constructor
  /**
   */
  LevelDataTemperatureIBC(RefCountedPtr<LevelData<FArrayBox> > a_temp, 
#if BISICLES_Z == BISICLES_LAYERED
			  RefCountedPtr<LevelData<FArrayBox> > a_surfaceT,
			  RefCountedPtr<LevelData<FArrayBox> > a_basalT,
#endif
			  const RealVect& a_dx);

  /// Destructor
  /**
   */
  virtual ~LevelDataTemperatureIBC();

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
  virtual LevelDataTemperatureIBC* new_temperatureIBC();

  /// set a basal heat flux. units are Joules / Year
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


  static LevelDataTemperatureIBC* parse(ParmParse& a_pp);

private:
  


  // Disallowed for all the usual reasons
  void operator=(const LevelDataTemperatureIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelDataTemperatureIBC(const LevelDataTemperatureIBC& a_input)
    :m_temp(a_input.m_temp)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
