#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

/// Calving model based on crevasse opening
/// see Van der Veen (1998) , Cold Regions Sci. Tech. vol 27 p 231

#include "CalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "VdVCalvingModel.H"
#include "CrevasseF_F.H"
#include "NamespaceHeader.H"
#include <ctime>

void VdVCalvingModel::modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
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
	  m_selected[a_level] = NULL;
	}
      //create one element of m_selected
      m_selected[a_level] = new LevelData<BaseFab<int> > (levelCoords.grids(),1,IntVect::Zero);
  
      
      // expand m_stressIntensity if < a_level
      if (m_stressIntensity.size() < a_level + 1)
	{
	  m_stressIntensity.resize(a_level + 1,NULL);
	}
      //clear up any old memory
      if (m_stressIntensity[a_level] != NULL)
	{
	  delete m_stressIntensity[a_level];
	  m_stressIntensity[a_level] = NULL;
	}
      //create one element of _stressIntensity 
      m_stressIntensity[a_level] = new LevelData<FArrayBox> (levelCoords.grids(),20,IntVect::Unit);


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
	  
	  FArrayBox s1(remnantBox,1); // first principal stress
	  //compute the first principal stress [needs to go in a fortran kernel]
	  for (BoxIterator bit(remnantBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      const Real& sxx = VT(iv,0);
	      const Real& syy = VT(iv,3);
	      const Real& sxy = 0.5 *(VT(iv,2) + VT(iv,1));
	      //vertically integrated first principal stress
	      s1(iv) = 0.5 * (sxx + syy)  
		+ std::sqrt (  std::pow ( 0.5*(sxx - syy), 2) + std::pow(sxy,2));
	      //vertically averaged first principal stress
	      s1(iv) *= thck(iv) / (1.0e-6 + thck(iv)*thck(iv));
	    }
	  
	  {
	    // Van Der Veen's net stress intensity for surface crevasses that 
	    // reach sea level or full thickness 
	    FArrayBox depth(remnantBox,1);
	    FArrayBox layerdepth(remnantBox,1);
	    FArrayBox& vdvK = (*m_stressIntensity[a_level])[dit];
	    
	    for (BoxIterator bit(remnantBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      depth(iv) = std::min(usrf(iv),.999*thck(iv));
	    }
	    
	    FArrayBox waterDepth(remnantBox,1);waterDepth.setVal(m_waterDepth);


	    for (int l =0; l < 10; l++)
	      {
		layerdepth.copy(depth);
		layerdepth*=Real(l+1)/10.0;
		FORT_VDVSTRESSS( CHF_FRA1(vdvK,l),
				 CHF_CONST_FRA1(thck,0),
				 CHF_CONST_FRA1(waterDepth,0),
				 CHF_CONST_FRA1(layerdepth,0),
				 CHF_CONST_FRA1(s1,0),
				 CHF_CONST_REAL(rhoi),
				 CHF_CONST_REAL(rhoo),
				 CHF_CONST_REAL(grav),
				 CHF_BOX(remnantBox));
	      }
	   

	    // Van Der Veen's net stress intensity for basal crevasses that 
	    // reach sea level 

	    FArrayBox thckp(remnantBox,1); thckp.setVal(0.0);
	    for (BoxIterator bit(remnantBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      thckp(iv) =  0.9*std::max(thck(iv)-usrf(iv),0.0);
	      depth(iv) = thck(iv)-depth(iv) ; //are there other possible relations?
	    }

	    for (int l =0; l < 10; l++)
	      {
		layerdepth.copy(depth);
		layerdepth*=Real(l+1)/10.0;
		
		FORT_VDVSTRESSB( CHF_FRA1(vdvK,l+10),
				 CHF_CONST_FRA1(thck,0),
				 CHF_CONST_FRA1(thckp,0),
				 CHF_CONST_FRA1(layerdepth,0),
				 CHF_CONST_FRA1(s1,0),
				 CHF_CONST_REAL(rhoi),
				 CHF_CONST_REAL(rhoo),
				 CHF_CONST_REAL(grav),
				 CHF_BOX(remnantBox));
	      }

	  }
			   
	  //compute the total remnant size (surface + basal) [needs to go in a fortran kernel]
	  FArrayBox remnant(remnantBox,1);
	  FArrayBox& vdvK = (*m_stressIntensity[a_level])[dit];
	  //set selected to zero everywhere
	  selected.setVal(0);
	  for (BoxIterator bit(remnantBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();

	      //surface crevasse depth
	      Real Ds = std::max(s1(iv),0.0) / (grav*rhoi) + rhoi/rhoo * m_waterDepth;
	      
	      Real random1 = ((Real)(rand())+0.0001)/((Real)(RAND_MAX)+0.0001);
	      Real random2 = ((Real)(rand())+0.0001)/((Real)(RAND_MAX)+0.0001);
	      
	      Real normalRandom = std::cos(8.*std::atan(1.)*random2)*std::sqrt(-2.*std::log(random1));
	      
	      Real noise = normalRandom * m_NoiseScale;

	      if (m_inclBasalCrev == true)
		{
		  //explicit basal crevasse depth calculation
		  //Real Db = (rhoi/(rhoo-rhoi)) * ((s1(iv)/(grav*rhoi)) - Hab(iv));
		  //remnant(iv) = std::min(thck(iv),thck(iv) - (Db + Ds + noise));
		   
		   if ( vdvK(iv,0) + vdvK(iv,1) > 0.1e6)
		     {
		       remnant(iv) = 0.0;
		     }
		}
	      else
		{
		  //assume basal crevasse reaches water level if surface crevasses do 
		  //remnant(iv) = std::min(thck(iv),usrf(iv)-(Ds + noise));
		  if ( vdvK(iv,0)  > 0.1e6)
		     {
		       remnant(iv) = 0.0;
		     }
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
  if (a_level == a_amrIce.finestLevel())
    {
      std::string file("vdvK.2d.hdf5");
      Real dt = 0.0; 
      Real time = 0.0;
      Vector<std::string> name(20);
      for (int l =0; l < 10; l++)
	{
	  char str[4];
	  sprintf(str,"ks%d",l);
	  name[l] = str;
	  sprintf(str,"kb%d",l);
	  name[l+10] = str;
	}

      name[0] = "ks";name[1]="kb";
      WriteAMRHierarchyHDF5(file ,a_amrIce.grids(), m_stressIntensity ,name, 
			    m_stressIntensity[0]->disjointBoxLayout().physDomain().domainBox(),
			    a_amrIce.dx(0)[0], dt, time, a_amrIce.refRatios(), m_stressIntensity.size());
      
    }
}
//re-apply internal domainEdgeCalvingModel after a regrid
void VdVCalvingModel::regridModifyState(LevelData<FArrayBox>& a_thickness, 
			       const AmrIce& a_amrIce,
			       int a_level)
{
  m_domainEdgeCalvingModel.regridModifyState(a_thickness, a_amrIce, a_level);
}

void VdVCalvingModel::initialModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{
  //endTimeStepModifyState(a_thickness, a_amrIce,a_level);
m_domainEdgeCalvingModel.endTimeStepModifyState(a_thickness, a_amrIce, a_level);
}
	      		 
void VdVCalvingModel::endTimeStepModifyState
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

VdVCalvingModel::~VdVCalvingModel()
{
  for (int i = 0; i < m_selected.size(); i++)
    {
      if (m_selected[i] != NULL)
	{
	  delete m_selected[i];
	  m_selected[i] = NULL;
	}
    }

  for (int i = 0; i < m_stressIntensity.size(); i++)
    {
      if (m_stressIntensity[i] != NULL)
	{
	  delete m_stressIntensity[i];
	  m_stressIntensity[i] = NULL;
	}
    }


}
#include "NamespaceFooter.H"	      
	      
	  

	  
