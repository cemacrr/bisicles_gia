#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef AMR_ICE_BASE_H
#define AMR_ICE_BASE_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"

/// abstract base class for amr ice sheet models (AmrIce, AMRIceControl)
class AmrIceBase
{

public:

  virtual ~AmrIceBase(){}

  virtual Real time() const = 0;

  //acess the ice sheet geometry
  //virtual int finestLevel() const = 0;
  virtual const RealVect dx(int a_level) const = 0;
  virtual const RefCountedPtr<LevelSigmaCS> geometry(int a_level) const = 0;
  virtual const DisjointBoxLayout& grids(int a_level) const = 0;
  //virtual const Vector<DisjointBoxLayout>& grids() const = 0;
  //virtual const Vector<RefCountedPtr<LevelSigmaCS> >& amrGeometry() const = 0;
  //virtual const Vector<int>& refRatios() const  = 0;
  //virtual const Vector<Real>& amrDx() const = 0;

  //access to the velocity
  //virtual const LevelData<FArrayBox>* velocity(int a_level) const = 0;
  //virtual const Vector<LevelData<FArrayBox>* >& amrVelocity() const = 0;

  //access to the upper surface heat boundary data
  virtual SurfaceFlux& surfaceHeatBoundaryData() const = 0;
  virtual bool surfaceHeatBoundaryDirichlett() const = 0;
  virtual bool surfaceHeatBoundaryTemperature() const = 0;
  
  //access to the lower surface heat boundary data
  virtual SurfaceFlux& basalHeatBoundaryData() const =0;
   
  //grouding line proximity
  virtual const LevelData<FArrayBox>* groundingLineProximity(int a_level) const = 0;
  virtual Real groundingLineProximityScale() const = 0;
  
  // viscous tensor
  enum viscousTensorComponent{xx,xy,yx,yy,max};
  //virtual const LevelData<FArrayBox>* viscousTensor(int a_level) const = 0;
  //virtual const LevelData<FArrayBox>* dragCoefficient(int a_level) const = 0;
  // virtual const LevelData<FArrayBox>* viscosityCoefficient(int a_level) const = 0;

  /// access function for layered face-centered velocities
  //virtual const Vector<LevelData<FluxBox>* >& faceVelocitiesLayered() const = 0;
  /// access function for calved ice thickness
  //virtual const  Vector<LevelData<FArrayBox>* >& calvedIceThickness() const = 0;
  /// access function for surface temperature
  //virtual const  Vector<LevelData<FArrayBox>* >& surfaceInternalEnergy() const  = 0;
  /// access function for surface heat flux
  //virtual const  Vector<LevelData<FArrayBox>* >& surfaceHeatFlux() const = 0;

};

#include "NamespaceFooter.H"
#endif
