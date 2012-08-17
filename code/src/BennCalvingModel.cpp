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
//#define XX 0
//#define YX 1
//#define XY 2
//#define YY 3

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

  for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& source = a_flux[dit];
      const FArrayBox& VT = vt[dit];
      const FArrayBox& thck = levelCoords.getH()[dit];
      const FArrayBox& topg = levelCoords.getTopography()[dit];
      const FArrayBox& usrf = levelCoords.getSurfaceHeight()[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) == FLOATINGMASKVAL)
	    {

	      Real sxx = VT(iv,0);
	      Real syy = VT(iv,3);
	      Real sxy = 0.5 *(VT(iv,2) + VT(iv,1));

	      Real s1 = 0.5 * (sxx + syy)  + std::sqrt ( std::pow ( 0.5*(sxx - syy), 2)
							 + std::pow(sxy,2));

	      s1 *= thck(iv) / (1.0e-6 + thck(iv)*thck(iv));
	      Real Ds = s1 / (grav*rhoi) + rhoi/rhoo * m_waterDepth;
	      if (usrf(iv) - Ds < 1.0e-10)
		{
		  for (int dir =0; dir < SpaceDim; dir++)
		    {
		      for (int sign = -1; sign <= 1; sign+=2)
			{
			  IntVect ivp = iv + sign*BASISV(dir);			  
			  if (usrf(ivp) < 1.0)
			    {
			      source(iv) -= 10.0 * thck(iv);
			    }
			}
		    }

		}
	      //single floating pixels have a stupid viscous tensor at the mo
	      bool single = true;
	      for (int dir =0; dir < SpaceDim; dir++)
		{
		  for (int sign = -1; sign <= 1; sign+=2)
		    {
		      IntVect ivp = iv + sign*BASISV(dir);
		      single &= thck(ivp) < 1.0;
		    }
		}
	      if (single)
		source(iv) -=  10.0 *thck(iv);

	    }

	}
    }
  
}

void BennCalvingModel::initialModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{
  endTimeStepModifyState(a_thickness, a_amrIce,a_level);
}

					 
void BennCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{

  m_domainEdgeCalvingModel.endTimeStepModifyState(a_thickness, a_amrIce, a_level);
  
  {
    const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
    for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
      {
	const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	FArrayBox& thck = a_thickness[dit];
	const Box& b = levelCoords.grids()[dit];
	for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) == OPENSEAMASKVAL || mask(iv) == FLOATINGMASKVAL)
	    {
	      bool single = true;
	      for (int dir =0; dir < SpaceDim; dir++)
		{
		  for (int sign = -1; sign <= 1; sign+=2)
		    {
		      IntVect ivp = iv + sign*BASISV(dir);
		      
		      //single &= ((mask(ivp) == OPENSEAMASKVAL) || (mask(ivp) == OPENLANDMASKVAL));
		      single &= thck(ivp) < 1.0;
		    }
		}
	      if (single || thck(iv) < 50.0)
		thck(iv) = 0.0;
	    }

	  



	}
      }
  }
  pout() << "Calving... BennCalvingModel::endTimeStepModifyState" << endl;
}

#include "NamespaceFooter.H"	      
	      
	  

	  
