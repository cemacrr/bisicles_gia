#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BENNCALVINGMODEL_H_
#define _BENNCALVINGMODEL_H_


#include "LevelSigmaCS.H"
#include "CalvingModel.H"

#include "NamespaceHeader.H"


/// Calving model based on viscous stress 
/// see (for 1D) Benn (2007), Earth-Science Reviews vol 82 p 143 doi:10.1016/j.earscirev.2007.02.002
/// Nick et al (2011) J. Glaciol vol 56 p 781 doi:10.3189/002214310794457344
/// author : Andrew Taylor, University of Bristol
class BennCalvingModel : public CalvingModel
{

public:

  
  virtual void initialModifyState(LevelData<FArrayBox>& a_thickness, 
				  const AmrIce& a_amrIce,
				  int a_level);


  virtual void endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
				      const AmrIce& a_amrIce,
				      int a_level);


  virtual void modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIce& a_amrIce,
					  int a_level,
					  Real a_dt);

  BennCalvingModel (const Real& a_waterDepth,
		    const Vector<int>& a_frontLo,
		    const Vector<int>& a_frontHi,
		    bool a_preserveSea, 
		    bool a_preserveLand,
		    bool a_inclBasalCrev,
		    bool a_oneDimCrev,
		    bool a_RedFreqCalv,
		    bool a_setZeroThck,
		    bool a_oldMeltRate,
		    const Real& a_NoiseScale,
		    const Real& a_calvingFreq,
		    const Real& a_decay,
		    const Real& a_timescale)
    : m_waterDepth(a_waterDepth),
      m_inclBasalCrev(a_inclBasalCrev),
      m_oneDimCrev(a_oneDimCrev),
      m_RedFreqCalv(a_RedFreqCalv),
      m_setZeroThck(a_setZeroThck),
      m_oldMeltRate(a_oldMeltRate),
      m_NoiseScale(a_NoiseScale),
      m_calvingFreq(a_calvingFreq),
      m_decay(a_decay),
      m_timescale(a_timescale),
      m_domainEdgeCalvingModel(a_frontLo,a_frontHi, 
			       a_preserveSea, a_preserveLand)
  {
    
  }

  ~BennCalvingModel();

private:
  Real m_waterDepth;
  bool m_inclBasalCrev;
  bool m_oneDimCrev;
  bool m_RedFreqCalv;
  bool m_setZeroThck;
  bool m_oldMeltRate;
  Real m_NoiseScale;
  Real m_calvingFreq;
  Real m_decay;
  Real m_timescale;
  DomainEdgeCalvingModel(m_domainEdgeCalvingModel);
  Vector<LevelData<BaseFab<int> >* > m_selected;

  //BennCalvingModel() : 
};

#include "NamespaceFooter.H"

#endif
