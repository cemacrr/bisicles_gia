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
#include "NamespaceHeader.H"

void 
DeglaciationCalvingModelA::postUpdateThickness
(LevelSigmaCS& a_coordSys, 
 LevelData<FArrayBox>& a_vel) const
{

  const DisjointBoxLayout& levelGrids = a_coordSys.grids();
  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
      FArrayBox& thck = a_coordSys.getH()[dit];
      const FArrayBox& topg = a_coordSys.getTopography()[dit];
      Box b = thck.box();
      b &= thck.box();
      b &= topg.box();
      FArrayBox water(b,1); //ocean depth beneath ice shelf
      water.copy(a_coordSys.getSurfaceHeight()[dit]);
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
	  else if (mask(iv) == FLOATINGMASKVAL)
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


void DomainEdgeCalvingModel::postUpdateThickness
(LevelSigmaCS& a_coordSys, 
 LevelData<FArrayBox>& a_vel) const
{
  const DisjointBoxLayout& grids = a_coordSys.grids();
  const ProblemDomain domain = grids.physDomain();
  LevelData<FArrayBox>& levelH = a_coordSys.getH();
  const IntVect ghost = levelH.ghostVect();
  
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
		  loBox &= levelH[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      levelH[dit](bit()) = 0.0;
		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
		  hiBox &= levelH[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      levelH[dit](bit()) = 0.0;
		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
    } // end loop over boxes

}


#include "NamespaceFooter.H"
