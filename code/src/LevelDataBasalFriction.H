#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// coefficient of basal friction stored in a LevelData 
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

class MultiLevelDataBasalFriction : public BasalFriction
{
public:
  /// constructor
  MultiLevelDataBasalFriction
  (const Vector<RefCountedPtr<LevelData<FArrayBox> > >&a_C,
   const RealVect& a_dxCrse, const Vector<int> a_ratio)
    :m_C(a_C),m_dxCrse(a_dxCrse),m_ratio(a_ratio),m_verbose(true)
  {
  }
  
  /// destructor
  virtual ~MultiLevelDataBasalFriction()
  {;}
  
  /// factory method
  /** return a pointer to a new BasalFriction object
   */
  virtual BasalFriction* new_basalFriction() const
  {
    MultiLevelDataBasalFriction* ptr = new MultiLevelDataBasalFriction(m_C,m_dxCrse,m_ratio);
    return static_cast<BasalFriction*>(ptr);
  }

  void setBasalFriction(LevelData<FArrayBox>& a_C, 
			RealVect a_dx,
			Real a_time,
			Real a_dt)
  {
    RealVect dx(m_dxCrse);
    for (int refDataLev = 0; refDataLev < m_C.size(); refDataLev++)
    {
     
      if (a_dx[0] < dx[0] && refDataLev > 0)
	{
	  //We need a LevelData<FArrayBox> whose DisjointBoxLayout covers
	  //a_C's, but m_C[refDataLev] doesn't necessarily provide one. Since we secretly
	  //want to be lisp programmers, we turn to recursion to build one,
	  //throwing performance worries to the wind.
	  //\todo consider modifying the BasalFriction interface to avoid recompuing levels 1 to n-1
	  const int& nRef = int(dx[0]/a_dx[0]);
	  DisjointBoxLayout crseDBL;
	  coarsen(crseDBL, a_C.disjointBoxLayout(), nRef);
	  LevelData<FArrayBox> crseC(crseDBL, 1 , IntVect::Zero);
	  this->setBasalFriction(crseC , dx, a_time, a_dt);
	  FillFromReference(a_C, crseC , a_dx, dx ,m_verbose);
	}
      else
	{
	  FillFromReference(a_C, *m_C[refDataLev], a_dx, dx ,m_verbose);
	}
      dx /= Real(m_ratio[refDataLev]);
    }
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
    setBasalFriction(a_C,a_coordSys.dx(),a_time,a_dt);
  }

protected:
  Vector<RefCountedPtr< LevelData<FArrayBox> > > m_C;
  RealVect m_dxCrse;
  Vector<int> m_ratio;
  bool m_verbose;

private:
  void operator=(const MultiLevelDataBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  MultiLevelDataBasalFriction(const MultiLevelDataBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"

#endif
