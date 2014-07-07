#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

/// Calving model based on viscous stress 
/// see (for 1D) Benn (2007), Earth-Science Reviews vol 82 p 143 doi:10.1016/j.earscirev.2007.02.002
/// Nick et al (2011) J. Glaciol vol 56 p 781 doi:10.3189/002214310794457344
/// author : Andrew Taylor, University of Bristol

#include "CalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "BennCalvingModel.H"
#include "NamespaceHeader.H"
#include <ctime>

void BennCalvingModel::modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
						  const AmrIce& a_amrIce,
						  int a_level,
						  Real a_dt)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const Real& rhoi = levelCoords.iceDensity();
  const Real& rhoo = levelCoords.waterDensity();
  const Real& grav = levelCoords.gravity();
  const LevelData<FArrayBox>& vt  = *a_amrIce.viscousTensor(a_level);
  const LevelData<FArrayBox>& levelVel = *(a_amrIce.velocity(a_level));
  srand(time(NULL)+ a_level + procID());
  const Real runTime = a_amrIce.time();
  Real roundedTime = std::floor(runTime + (1-m_calvingFreq));
  bool calvingActive = (roundedTime <= floor(runTime));
  pout() << " calvingActive = " << calvingActive
         << std::endl;
  if (calvingActive)
    {
      
      // expand m_selected if < a_level
      if (m_selected.size() < a_level + 1)
	{
	  m_selected.resize(a_level + 1,NULL);
	}
      //clear up any old memory
      if (m_selected[a_level] != NULL)
	{
	  delete m_selected[a_level];
	}
      //create one element of m_selected
      m_selected[a_level] = new LevelData<BaseFab<int> > (levelCoords.grids(),1,IntVect::Zero);
  
      for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
	{
	  FArrayBox& source = a_flux[dit];
	  const FArrayBox& thck = levelCoords.getH()[dit];
	  const FArrayBox& VT = vt[dit];
	  const FArrayBox& usrf = levelCoords.getSurfaceHeight()[dit];
	  const FArrayBox& vel = levelVel[dit];
	  const FArrayBox& Hab = levelCoords.getThicknessOverFlotation()[dit];
	  BaseFab<int>& selected = (*m_selected[a_level])[dit];
	  const Box& b = levelCoords.grids()[dit];

	  Box remnantBox = b; remnantBox.grow(1); //we need one layer of ghost cells to do upwind calculations
	  FArrayBox remnant(remnantBox,1);
	  //compute the total remnant size (surface + basal) [needs to go in a fortran kernel]
	  
	  //set selected to zero everywhere
	  selected.setVal(0);
	  for (BoxIterator bit(remnantBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      const Real& sxx = VT(iv,0);
	      const Real& syy = VT(iv,3);
	      const Real& sxy = 0.5 *(VT(iv,2) + VT(iv,1));

	      //vertically integrated first principal stress
	      Real s1 = 0.5 * (sxx + syy)  
		+ std::sqrt (  std::pow ( 0.5*(sxx - syy), 2) + std::pow(sxy,2));
	      
	      //vertically averaged first principal stress
	      s1 *= thck(iv) / (1.0e-6 + thck(iv)*thck(iv));

	      //surface crevasse depth
	      Real Ds = std::max(s1,0.0) / (grav*rhoi) + rhoi/rhoo * m_waterDepth;
	      
	      Real random1 = ((Real)(rand())+0.0001)/((Real)(RAND_MAX)+0.0001);
	      Real random2 = ((Real)(rand())+0.0001)/((Real)(RAND_MAX)+0.0001);
	      
	      Real normalRandom = std::cos(8.*std::atan(1.)*random2)*std::sqrt(-2.*std::log(random1));
	      
	      Real noise = normalRandom * m_NoiseScale;

	      if (m_inclBasalCrev == true)
		{
		  //explicit basal crevasse depth calculation
		  Real Db = (rhoi/(rhoo-rhoi)) * ((s1/(grav*rhoi)) - Hab(iv));
		  remnant(iv) = std::min(thck(iv),thck(iv) - (Db + Ds + noise));
		}
	      else
		{
		  //assume basal crevasse reaches water level if surface crevasses do 
		  remnant(iv) = std::min(thck(iv),usrf(iv)-(Ds + noise));
		}
	      remnant(iv) = std::max(-0.0, remnant(iv));
	    }
	  
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      
	      //compute a melt-rate which should be dominated by the contribution
	      //of thick upwind cells: it will only be large if all upwind cells
	      //are close to fractured and the current cell
	 
	      Real upwRemnant = remnant(iv)*thck(iv) ;
	      Real umod = std::sqrt(vel(iv,0)*vel(iv,0)+vel(iv,1)*vel(iv,1)) + 1.0e-10;
	      Real upwThck = thck(iv) + 1.0e-10; 
	      
	      for (int dir = 0; dir < SpaceDim; dir++)
		{
		  for (int side = -1; side <=1; side += 2)
		    {
		      const IntVect ivu = iv + side*BASISV(dir);
		      Real fu = (side < 0)?
			(std::max(vel(ivu,dir),0.0)):
			(std::max(-vel(ivu,dir),0.0));
		      fu /= umod;
		      
		      upwThck += thck(ivu)*fu;
		      upwRemnant += thck(ivu)*remnant(ivu)*fu;
		      
		    }
		}
	      // set mask to 1 to indicate calving event
	      if ((m_setZeroThck == true) && (upwRemnant < 5.0))
		{
		  selected(iv) = 1;
		}
	      else if (m_oldMeltRate == true) // use a_dt melt rate
		{
		  const Real decay = 3.0e+2;
		  source(iv) -= 10.0 * thck(iv)/a_dt * std::min(upwThck,1.0) 
		    * std::min(1.0, std::exp( - decay * upwRemnant/(upwThck)));
		}
	      else // use decay and ts melt rate
		{
		  source(iv) -=  thck(iv)/m_timescale * std::min(upwThck,1.0) 
		    * std::min(1.0, std::exp( - m_decay * upwRemnant/(upwThck)));
		  
		}
	    }  
	}
    }
}


void BennCalvingModel::initialModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{
  //endTimeStepModifyState(a_thickness, a_amrIce,a_level);
m_domainEdgeCalvingModel.endTimeStepModifyState(a_thickness, a_amrIce, a_level);
}
	      		 
void BennCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)

{
  if ((m_setZeroThck == true) && (m_selected.size() >= a_level + 1) 
      && (m_selected[a_level] != NULL))
    {
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
       {
	 FArrayBox& thickness = a_thickness[dit];
	 const BaseFab<int>& selected = (*m_selected[a_level])[dit];
	 const Box& b = levelCoords.grids()[dit];      
	 for (BoxIterator bit(b); bit.ok(); ++bit)
	   {
	     const IntVect& iv = bit();
	     if (selected(iv) == 1)
	       {
		 thickness(iv) = 0.0;
	       }
	   }
       }
   }  
  m_domainEdgeCalvingModel.endTimeStepModifyState(a_thickness, a_amrIce, a_level);

}

BennCalvingModel::~BennCalvingModel()
{
  for (int i = 0; i < m_selected.size(); i++)
    {
      if (m_selected[i] != NULL)
	{
	  delete m_selected[i];
	  m_selected[i] = NULL;
	}
    }
}
#include "NamespaceFooter.H"	      
	      
	  

	  
