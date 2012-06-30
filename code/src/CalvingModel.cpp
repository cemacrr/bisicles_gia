#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "NamespaceHeader.H"

void 
DeglaciationCalvingModelA::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{

  Real time = a_amrIce.time();
  bool calvingActive = (time >= m_startTime && time < m_endTime);
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
 
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      const FArrayBox& topg = levelCoords.getTopography()[dit];
      Box b = thck.box();
      b &= topg.box();
      FArrayBox water(b,1); //ocean depth beneath ice shelf
      water.copy(levelCoords.getSurfaceHeight()[dit]);
      water -= thck;
      water -= topg;

      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();

	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (calvingActive && mask(iv) == FLOATINGMASKVAL)
	    {
	      if (thck(iv) < m_calvingThickness && water(iv) > m_calvingOceanDepth) 
		thck(iv) =  0.0;
	    }
	  
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	}


    }

}

void DomainEdgeCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
  const IntVect ghost = a_thickness.ghostVect();
  
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      //const Box& gridBox = grids[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  if (!domain.isPeriodic(dir))
	    {

	      if (m_frontLo[dir] > 0)
		{
		  Box loBox = adjCellLo(domain,dir,ghost[dir]);
		  loBox &= a_thickness[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv + BASISV(dir);
		      if (levelMask[dit](ip) != GROUNDEDMASKVAL)
			a_thickness[dit](iv) = 0.0;
		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
		  hiBox &= a_thickness[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv - BASISV(dir);
		      if (levelMask[dit](ip) != GROUNDEDMASKVAL)
			a_thickness[dit](iv) = 0.0;
		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	}

    } // end loop over boxes

}



#include "NamespaceFooter.H"
