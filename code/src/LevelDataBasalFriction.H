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
#include "NamespaceHeader.H"


class LevelDataBasalFriction : public BasalFriction
{
public:
  /// constructor
  LevelDataBasalFriction(RefCountedPtr<LevelData<FArrayBox> > a_C,
			 const RealVect& a_dx)
    :m_timeFileMap(NULL),m_startTime(-1.2345678e+300),m_endTime(1.2345678e+300),
     m_C(a_C),m_endC(a_C),m_dx(a_dx),m_verbose(true),m_name("")
  {
  }
  
  
  LevelDataBasalFriction(RefCountedPtr<std::map<Real, std::string> >& a_timeFileMap, 
		       const std::string& a_name)
    :m_timeFileMap(a_timeFileMap),m_startTime(1.2345678e+299),m_endTime(-1.2345678e+299),
     m_C(new LevelData<FArrayBox>()),m_endC(new LevelData<FArrayBox>())
    ,m_verbose(true),m_name(a_name)
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
    LevelDataBasalFriction* ptr;
    if (m_timeFileMap == NULL)
      {
	ptr = new LevelDataBasalFriction(m_C,m_dx);
      }
    else
      {
	RefCountedPtr<std::map<double, std::basic_string<char> > > r = m_timeFileMap;
	  
	ptr = new LevelDataBasalFriction(r,m_name);
      }
    return static_cast<BasalFriction*>(ptr);
  }

  /// define basal friction coefficient C and place in a_betaSqr
  /** Either copy data from a single stored level, or load files as needed
      depending on a_time
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);

  
protected:
  RefCountedPtr<std::map<Real,std::string> > m_timeFileMap;
  Real m_startTime, m_endTime;
  RefCountedPtr< LevelData<FArrayBox> > m_C; // C at t = m_startTime
  RefCountedPtr< LevelData<FArrayBox> > m_endC; // C at t = m_endTime
  RealVect m_dx;
  bool m_verbose;
  std::string m_name;

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
  /// define basal friction coefficient C and place in a_betaSqr
  /** Copy data from a single stored AMR HIerarchy
  */
  void setBasalFriction(LevelData<FArrayBox>& a_C, 
			RealVect a_dx,
			Real a_time,
			Real a_dt);
 
  /// define basal friction coefficient C and place in a_betaSqr
  virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt)
  {
    setBasalFriction( a_C,a_coordSys.dx(),a_time,a_dt);
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
