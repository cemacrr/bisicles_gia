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

  DomainEdgeCalvingModel(const Vector<int>& a_frontLo,
			 const Vector<int>& a_frontHi)
    
    :m_frontLo(a_frontLo),m_frontHi(a_frontHi)
  {
    ;
  }

protected:
  Vector<int> m_frontLo;
  Vector<int> m_frontHi;
private:
  DomainEdgeCalvingModel();
};


#include "NamespaceFooter.H"

#endif
