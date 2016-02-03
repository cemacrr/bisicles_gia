
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

class AmrIce;

///virtual base for calving models
class CalvingModel
{

public:

  //alter the thickness field prior to the first time step. 
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness,  
				  const AmrIce& a_amrIce,
				  int a_level) = 0;


  //alter the thickness field at the end of a time step
  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level) = 0;
  
//alter the thickness field after a regid
  virtual void regridModifyState(LevelData<FArrayBox>& a_thickness, 
				 const AmrIce& a_amrIce,
				 int a_level) = 0;

  //modify a thickness flux
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  LevelData<FArrayBox>& a_calvingMelt,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt) =0;


  virtual ~CalvingModel(){};

  static CalvingModel* parseCalvingModel(const char* a_prefix);


};

///the simplest calving model : don't do anything
class NoCalvingModel : public CalvingModel
{

public:

  
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness, 
				  const AmrIce& a_amrIce,
				  int a_level)
  {

  }

  
  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level)
  {

  }

  virtual void regridModifyState(LevelData<FArrayBox>& a_thickness, 
				 const AmrIce& a_amrIce,
				 int a_level) 
  {

  }
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  LevelData<FArrayBox>& a_calvingMelt,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt)
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

  
  //alter the thickness field prior to the first time step. Same as endTimeStepModifyState
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness, 
				  const AmrIce& a_amrIce,
				  int a_level)
  {
    endTimeStepModifyState(a_thickness, a_amrIce, a_level);
  }


  //alter the thickness field at the end of a time step
  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level);

  //alter the thickness field after a regrid
  virtual void regridModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level)
  {
    //do nothing... for backward compatibility
  }


  //modify a thickness flux
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  LevelData<FArrayBox>& a_calvingMelt,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt)
  {

  }


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

  
  //alter the thickness field prior to the first time step. 
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness, 
				  const AmrIce& a_amrIce,
				  int a_level)
  {
    endTimeStepModifyState(a_thickness, a_amrIce, a_level);
  }


  //alter the thickness field at the end of a time step
  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level);

  //alter the thickness field after a regrid
  virtual void regridModifyState(LevelData<FArrayBox>& a_thickness, 
				 const AmrIce& a_amrIce,
				 int a_level)
  {
    endTimeStepModifyState(a_thickness, a_amrIce, a_level);
  }

  //modify a thickness flux
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  LevelData<FArrayBox>& a_calvingMelt,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt)
  {

  }


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


///Calve ice based on proximity to the grounding line
class ProximityCalvingModel : public CalvingModel
{

  Real m_proximity, m_velocity,  m_startTime, m_endTime;

public:
  ProximityCalvingModel(const Real& a_proximity, 
			const Real& a_velocity,
			const Real& a_startTime, 
			const Real& a_endTime)
    : m_proximity(a_proximity), 
      m_velocity(a_velocity),
      m_startTime(a_startTime), 
      m_endTime(a_endTime)
  {
  }
  
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness, 
				  const AmrIce& a_amrIce,
				  int a_level)
  {
    endTimeStepModifyState( a_thickness,  a_amrIce, a_level);
  }

  //alter the thickness field after a regrid
  virtual void regridModifyState(LevelData<FArrayBox>& a_thickness, 
				 const AmrIce& a_amrIce,
				 int a_level)
  {
    //nothing
  }
  
  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level);


  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  LevelData<FArrayBox>& a_calvingMelt,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt);


};



// A Deglacation-only calving model which 
// i. Makes calving irreversible : if a cell mask is open sea / land, set thickness
//    there to zero. No need to worry about the velocity fields in cells
//    which were previously calved.
// ii. Sets any floating ice below a given thickness to the min thickness in (iii) (this is where it differs from ModelA, which doesn't actually use the CalvingThickness)
// iii. prevents thickness from dropping below a given thickness
class DeglaciationCalvingModelB : public CalvingModel
{

public:

  
  //alter the thickness field prior to the first time step. Same as endTimeStepModifyState
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness, 
				  const AmrIce& a_amrIce,
				  int a_level)
  {
    endTimeStepModifyState(a_thickness, a_amrIce, a_level);
  }


  //alter the thickness field at the end of a time step
  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level);


  //alter the thickness field after a regrid
  virtual void regridModifyState(LevelData<FArrayBox>& a_thickness, 
                                 const AmrIce& a_amrIce,
                                 int a_level)
  {
    //do nothing... for backward compatibility
  }
  
  //modify a thickness flux
  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  LevelData<FArrayBox>& a_calvingMelt,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt)
  {

  }


  DeglaciationCalvingModelB(const Real& a_calvingThickness,
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
  DeglaciationCalvingModelB();
};




#include "NamespaceFooter.H"

#endif
