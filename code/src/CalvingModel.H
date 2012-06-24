#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CALVINGMODEL_H_
#define _CALVINGMODEL_H_

#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

///virtual base for calving models
class CalvingModel
{

public:
  //alter the geometry and velocity fields after 
  //thickness has ben updated but before 
  //a_coordSys.recomputeGeometry() is called.
  virtual void postUpdateThickness(LevelSigmaCS& a_coordSys,
				   LevelData<FArrayBox>& a_vel,
				   const Real& a_time) const=0;

  //modify a thickness flux
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const LevelSigmaCS& a_coordSys,
				    const LevelData<FArrayBox>& a_vel,
				    Real a_time, Real a_dt) const=0;


  virtual ~CalvingModel(){};


};

///the simplest calving model : don't do anything
class NoCalvingModel : public CalvingModel
{

public:

  virtual void postUpdateThickness(LevelSigmaCS& a_coordSys,
				   LevelData<FArrayBox>& a_vel,
				   const Real& a_time) const
  {

  }

  
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const LevelSigmaCS& a_coordSys,
				    const LevelData<FArrayBox>& a_vel,
				    Real a_time, Real a_dt) const
  {

  }


};


// A Deglacation-only calving model which 
// i. Makes calving irreversible : if a cell mask is open sea / land, set thickness
//    there to zero. No need to worry about the velocity fields in cells
//    which were previously calved.
// ii. Sets any floating ice below a given thickness and in a given depth of
//     ocean to zero thickness
// iii. prevents thickness from dropping below a given thickness
class DeglaciationCalvingModelA : public CalvingModel
{

public:
  virtual void postUpdateThickness(LevelSigmaCS& a_coordSys,
				   LevelData<FArrayBox>& a_vel,
				   const Real& a_time) const;
 

  
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const LevelSigmaCS& a_coordSys,
				    const LevelData<FArrayBox>& a_vel,
				    Real a_time, Real a_dt) const;
  

  DeglaciationCalvingModelA(const Real& a_calvingThickness,
			    const Real& a_calvingOceanDepth,
			    const Real& a_minThickness, 
			    const Real& a_startTime, 
			    const Real& a_endTime)
    : m_calvingThickness(a_calvingThickness),
      m_calvingOceanDepth(a_calvingOceanDepth),
      m_minThickness(a_minThickness),
      m_startTime(a_startTime),
      m_endTime(a_endTime)
  {
    ;
  }
 
protected:
  Real m_calvingThickness, m_calvingOceanDepth, m_minThickness, m_startTime, m_endTime;
  
private:
  DeglaciationCalvingModelA();
};

/// Ensure open ocean on one or more domain edges
class DomainEdgeCalvingModel :  public CalvingModel
{

public:
  virtual void postUpdateThickness(LevelSigmaCS& a_coordSys,
				   LevelData<FArrayBox>& a_vel,
				   const Real& a_time) const;

  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const LevelSigmaCS& a_coordSys,
				    const LevelData<FArrayBox>& a_vel,
				    Real a_time, Real a_dt) const;


  DomainEdgeCalvingModel(const Vector<int>& a_frontLo,
			 const Vector<int>& a_frontHi,
			 bool a_preserveSea, 
			 bool a_preserveLand)
    
    :m_frontLo(a_frontLo),m_frontHi(a_frontHi),
     m_preserveSea(a_preserveSea),m_preserveLand(a_preserveLand)
  {
    ;
  }

protected:
  Vector<int> m_frontLo;
  Vector<int> m_frontHi;
  bool m_preserveSea, m_preserveLand;
private:
  DomainEdgeCalvingModel();
};


#include "NamespaceFooter.H"

#endif
