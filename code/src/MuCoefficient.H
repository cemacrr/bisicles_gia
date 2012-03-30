#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Virtual base class encapsulating coefficient of effective viscosity 
// for ice-sheet model. NB : this is an additional coefficient a, so
// that the viscous tensor is (a * mu * rate-of-strain)

#ifndef _MUCOEFFICIENT_H_
#define _MUCOEFFICIENT_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "NamespaceHeader.H"


class MuCoefficient
{
public:
  /// constructor
  MuCoefficient() {};

  /// destructor
  virtual ~MuCoefficient() {};

  /// factory method
  /** return a pointer to a new MuCoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const = 0;
 /// define cell and face-centered coefficient
  virtual void setMuCoefficient(LevelData<FArrayBox>& a_muCoef,
				LevelData<FluxBox>& a_faceMuCoef,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt) = 0;


private:
  // Disallowed for all the usual reasons
  void operator=(const MuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MuCoefficient(const MuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


// simple version which sets mu coefficient to one
class UnitMuCoefficient : public MuCoefficient
{
public:
  /// constructor
  UnitMuCoefficient() {};

  /// destructor
  virtual ~UnitMuCoefficient() {};

  /// factory method
  /** return a pointer to a new Mucoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const;

  /// define cell and face-centered coefficient
  /** time and dt are included in case this is time-dependent.
  */
  virtual void setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
				LevelData<FluxBox>& a_faceMuCoef,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);


private:
  // Disallowed for all the usual reasons
  void operator=(const UnitMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  UnitMuCoefficient(const UnitMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

class LevelDataMuCoefficient : public MuCoefficient
{

public:
  //constructor
  LevelDataMuCoefficient
  (RefCountedPtr<LevelData<FArrayBox> > a_muCoef,
   const RealVect& a_dx)
    :m_muCoef(a_muCoef),m_dx(a_dx),m_verbose(true)
  {
  }
    /// factory method
  /** return a pointer to a new Mucoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const;

  virtual void setMuCoefficient
  (LevelData<FArrayBox>& a_muCoef,
   LevelData<FluxBox>& a_faceMuCoef,
   LevelSigmaCS& a_coordSys,
   Real a_time,
   Real a_dt) ;

protected:
  RefCountedPtr< LevelData<FArrayBox> > m_muCoef;
  RealVect m_dx;
  bool m_verbose;

private:
  void operator=(const LevelDataMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  LevelDataMuCoefficient(const LevelDataMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"


#endif // multiple-include preventer
