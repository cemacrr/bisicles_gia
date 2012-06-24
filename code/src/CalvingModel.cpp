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
 LevelData<FArrayBox>& a_vel, 
 const Real& a_time) const
{
  // pout() << "DeglaciationCalvingModelA::postUpdateThickness"
  // 	 << "time = " << a_time 
  // 	 << " start time = " << m_startTime 
  // 	 << " end time = " << m_endTime 
  // 	 << std::endl;
  bool calvingActive = (a_time >= m_startTime && a_time < m_endTime);
  //if (calvingActive)
    // pout() << " calving active " << std::endl;

  const DisjointBoxLayout& levelGrids = a_coordSys.grids();
  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
      FArrayBox& thck = a_coordSys.getH()[dit];
      FArrayBox& u = a_vel[dit];
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

      
      // if (a_time < 0.26)
      // 	{
      // 	  for (BoxIterator bit(u.box()); bit.ok(); ++bit)
      // 	    {
      // 	      const IntVect& iv = bit();
      // 	      if( mask(iv) == FLOATINGMASKVAL)
      // 		{
      // 		  Real umod = std::sqrt( u(iv,0)*u(iv,0) + u(iv,1)*u(iv,1));
      // 		  if (umod > 2.0e+3)
      // 		    {
      // 		      thck(iv) = 5.0;
      // 		      u(iv,0) = 0.0;
      // 		      u(iv,1) = 0.0;
      // 		    }
      // 		}
      // 	    }
      // 	}


    }

}

void DeglaciationCalvingModelA::modifySurfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const LevelSigmaCS& a_coordSys,
 const LevelData<FArrayBox>& a_vel,
 Real a_time, Real a_dt) const
{

  // const DisjointBoxLayout& levelGrids = a_coordSys.grids();
  // for (DataIterator dit(levelGrids); dit.ok(); ++dit)
  //   {
  //     const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
  
  //     if (a_time < -1.25)
  // 	{
  // 	  const FArrayBox& u = a_vel[dit];
  // 	  const FArrayBox& thck = a_coordSys.getH()[dit];
  // 	  FArrayBox& flux = a_flux[dit];
  // 	  for (BoxIterator bit(u.box()); bit.ok(); ++bit)
  // 	    {
  // 	      const IntVect& iv = bit();
  // 	      if (mask(iv) == FLOATINGMASKVAL)
  // 		{
  // 		  Real umod = std::sqrt(u(iv,0)*u(iv,0) + u(iv,1)*u(iv,1));
  // 		  if (umod > 1.5e+3)
  // 		    {
  // 		      flux(iv) -= 4 * thck(iv);
  // 		    }
  // 		}
  // 	    }
  // 	}
  //   }
}


void DomainEdgeCalvingModel::postUpdateThickness
(LevelSigmaCS& a_coordSys, 
 LevelData<FArrayBox>& a_vel, 
 const Real& a_time) const
{
  const DisjointBoxLayout& grids = a_coordSys.grids();
  const ProblemDomain domain = grids.physDomain();
  LevelData<FArrayBox>& levelH = a_coordSys.getH();
  const LevelData<BaseFab<int> >& levelMask = a_coordSys.getFloatingMask();
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
		      const IntVect& iv = bit();
		      const IntVect ip = iv + BASISV(dir);
		      if (levelMask[dit](ip) != GROUNDEDMASKVAL)
			levelH[dit](iv) = 0.0;
		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
		  hiBox &= levelH[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv - BASISV(dir);
		      if (levelMask[dit](ip) != GROUNDEDMASKVAL)
			levelH[dit](iv) = 0.0;
		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
      FArrayBox& thck = a_coordSys.getH()[dit];
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

void DomainEdgeCalvingModel::modifySurfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const LevelSigmaCS& a_coordSys,
 const LevelData<FArrayBox>& a_vel,
 Real a_time, Real a_dt) const
{

}


#include "NamespaceFooter.H"
