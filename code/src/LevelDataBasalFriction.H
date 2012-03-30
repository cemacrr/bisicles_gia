#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// coefficient of basal friction stored in a LevelDataIBC
#ifndef _LEVELDATABASALFRICTION_H_
#define _LEVELDATABASALFRICTION_H_

#include "BasalFriction.H"
#include "FillFromReference.H"
#include "NamespaceHeader.H"


class LevelDataBasalFriction : public BasalFriction
{
public:
  /// constructor
  LevelDataBasalFriction(RefCountedPtr<LevelData<FArrayBox> > a_C,
			 const RealVect& a_dx)
    :m_C(a_C),m_dx(a_dx),m_verbose(true)
  {
  }
  
  /// destructor
  virtual ~LevelDataBasalFriction()
  {;}
  
  /// factory method
  /** return a pointer to a new BasalFriction object
   */
  virtual BasalFriction* new_basalFriction() const
  {
    LevelDataBasalFriction* ptr = new LevelDataBasalFriction(m_C,m_dx);
    return static_cast<BasalFriction*>(ptr);
  }

  /// define basal friction coefficient C and place in a_betaSqr
  /** time and dt are included in case this is time-dependent. Units 
      should be Pa*a/m 
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
				Real a_dt)
  {
    FillFromReference(a_C, *m_C, a_coordSys.dx(),m_dx,m_verbose);
  }

protected:
  RefCountedPtr< LevelData<FArrayBox> > m_C;
  RealVect m_dx;
  bool m_verbose;

private:
  void operator=(const LevelDataBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  LevelDataBasalFriction(const LevelDataBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

 

};

#include "NamespaceFooter.H"

#endif
