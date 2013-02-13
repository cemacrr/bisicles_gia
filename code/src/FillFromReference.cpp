#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  FillFromReference.H
// ============
//
// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference FArrayBox
//

#include "FillFromReference.H"
#include "QuadCFInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "IceConstants.H"

#include "NamespaceHeader.H"

/// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference FArrayBox
/**
   
 */

void FillFromReference(LevelData<FArrayBox>& a_destData,
		       const FArrayBox& a_srcData,
		       const RealVect& a_destDx,
		       const RealVect& a_srcDx,
		       const IntVect& a_srcGhost,
		       const bool a_verbose)
{
  // tolerance used when converting refinement ratio from Real -> integer
  Real tolerance = 1.0e-6;

  if (a_verbose)
    {
      pout() << "in FillFromReference" << endl;
    }
  
  // refinement ratio
  Real refRatio = a_srcDx[0]/a_destDx[0];

  // reality check
  Real testDx = a_destDx[0]*refRatio;
  if (abs(testDx - a_srcDx[0])/a_srcDx[0] > TINY_NORM)
    {
      MayDay::Error("FillFromReference::incompatible src and dest dx");
    }

  // passed-in data on a LevelData
  Box srcBox(a_srcData.box());
  srcBox.grow(-a_srcGhost);
  Vector<Box> srcBoxes(1,srcBox);
  
  // grab periodicity info from a_coords
  const DisjointBoxLayout& destGrids = a_destData.getBoxes();
  const ProblemDomain& destDomain = destGrids.physDomain();

  ProblemDomain srcDomain(destDomain);
  // now need to either coarsen or refine srcDomain 
  if (refRatio > 1+tolerance) 
    {
      // srcDomain is coarser than destDomain
      int nRef = (int)(refRatio + tolerance);
      srcDomain.coarsen(nRef);
    }
  else if (refRatio < 1-tolerance)
    {
      // srcDomain is finer than destDomain
      Real refRatioInv = 1.0/refRatio;
      int nRef = (int) (refRatioInv + tolerance);
      srcDomain.refine(nRef);

      // for this case, in order for CoarseAverage to work, we
      // need to redefine the source to be a refinement of the 
      // destGrids
      Vector<Box> destBoxes = destGrids.boxArray();
      srcBoxes.resize(destBoxes.size());
      for (int i=0; i<srcBoxes.size(); i++)
        {
          srcBoxes[i] = destBoxes[i];
          srcBoxes[i].refine(nRef);
        }
    }
  else
    {
      // same size -- do nothing
    }

  Vector<int> procAssign(srcBoxes.size(),0);
  
  DisjointBoxLayout srcDBL(srcBoxes, procAssign, srcDomain);
  LevelData<FArrayBox> srcLD(srcDBL, 1);
  
  DataIterator dit = srcDBL.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      srcLD[dit].copy(a_srcData);
    }

  // now go through and either interpolate, average, or copy srcData->destData

  if (refRatio > 1+tolerance) 
    {
      // need to interpolate data
      int nRef = (int)(refRatio + tolerance); 
      
      if (a_verbose)
	{
	  pout() << " ...interpolating data with refinement ratio = " 
		 << nRef << endl;
	}

      FineInterp interpolator(destGrids, 1, nRef, destDomain);
      interpolator.interpToFine(a_destData, srcLD);

      // now fill in ghost cells, using PiecewiseLinearFillPatch
      IntVect ghostVect = a_destData.ghostVect();

      if (ghostVect[0] != 0)
        {
          PiecewiseLinearFillPatch pwl(destGrids, srcDBL,
                                       1, srcDomain, nRef, 
                                       ghostVect[0]);
          
          Real time_interp_coeff = 0.0;
          pwl.fillInterp(a_destData, srcLD, srcLD,
                         time_interp_coeff, 0, 0, 1);
          
        }
    }
  else if (refRatio < 1-tolerance)
    {
      // need to average data
      refRatio = 1.0/refRatio;
      int nRef = (int) (refRatio + tolerance);
      if (a_verbose)
	{
	  pout() << " ...averaging data with refinement ratio = " 
		 << nRef << endl;
	}
      CoarseAverage averager(srcDBL,destGrids,
			     1, nRef, a_destData.ghostVect());
      averager.averageToCoarse(a_destData, srcLD);
    }
  else
    {
      // same size
      if (a_verbose)
	{
	  pout() << " ...same-size copy of data" << endl;
	}

      srcLD.copyTo(a_destData);
    }


}

/// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference  LevelData<FArrayBox>
/**
   
 */
void FillFromReference(LevelData<FArrayBox>& a_destData,
		       const LevelData<FArrayBox>&  a_srcData,
		       const RealVect&  a_destDx,
		       const RealVect&  a_srcDx,
		       const bool a_verbose)
{

  if (a_verbose)
    {
      pout() << "in LDF FillFromReference" << endl;
    }

  // tolerance used when converting refinement ratio from Real -> integer 
  Real tolerance = 1.0e-6;
  // refinement ratio
  Real refRatio = a_srcDx[0]/a_destDx[0];
  // reality check
  Real testDx = a_destDx[0]*refRatio;
  if (abs(testDx - a_srcDx[0])/a_srcDx[0] > TINY_NORM)
    {
      MayDay::Error("FillFromReference(LevelData<FarrayBox>& ,const LevelData<FarrayBox>&,... ::incompatible src and dest dx");
    }

  const DisjointBoxLayout& srcGrids = a_srcData.disjointBoxLayout();
  const DisjointBoxLayout& destGrids = a_destData.disjointBoxLayout();

  if (a_verbose)
    {
      pout() << "refRatio = " << refRatio << ", srcDx = " << a_srcDx[0]
             << ", destDx = " << a_destDx[0] << endl;
    }

  if (refRatio > 1+tolerance) 
     {
       //interpolate data
       if (a_verbose)
         {
           pout() << "Interpolating data... " << endl;
         }

       int nRef = (int)(refRatio + tolerance); 
       bool coarsenable = destGrids.coarsenable(nRef);
       //avoid attempting to refine with a ratio greater than block_factor 
       // by interpolating coarse data in stages until nRef == block_factor
       if (coarsenable)
	 {
	   if (a_verbose)
	     {
	       pout() << " ...interpolating data with refinement ratio = " 
		      << nRef << endl;
	     }
	   
	   const ProblemDomain& fineDomain = destGrids.physDomain();
	   FineInterp interpolator(destGrids, a_destData.nComp(), nRef, fineDomain);
	   interpolator.interpToFine(a_destData, a_srcData);
	   
	   // now fill in ghost cells, using PiecewiseLinearFillPatch
	   IntVect ghostVect = a_destData.ghostVect();
	   if (ghostVect[0] != 0)
	     {
	       if (a_verbose)
		 {
		   pout() << " ...interpolating " << ghostVect[0] 
			  << " cells of ghost data along coarse-fine interfaces refinement ratio =  "
			  << nRef << endl;
		 }
	       ProblemDomain coarseDomain = fineDomain;
	       coarseDomain.coarsen(nRef);
	       //const ProblemDomain coarseDomain = srcGrids.physDomain();
	       PiecewiseLinearFillPatch pwl(destGrids, srcGrids, a_destData.nComp(), coarseDomain, nRef,  ghostVect[0]);
	       Real time_interp_coeff = 0.0;
	       pwl.fillInterp(a_destData, a_srcData, a_srcData,time_interp_coeff, 0, 0, a_destData.nComp());
	     }
	 }
       else if (nRef%2 == 0)
	 {
	   if (a_verbose)
	     {
	       pout() << " ...interpolating data with refinement ratio = " 
		      << nRef << endl
		      << " ...recursively (potentially memory intensive as this refines entire coarse levels -  consider increasing amr.block_factor) " << endl;
	     }
 
	   DisjointBoxLayout stepGrids;
	   refine(stepGrids,srcGrids,2);
	   LevelData<FArrayBox> stepData(stepGrids,a_srcData.nComp(), a_srcData.ghostVect());
	   FillFromReference(stepData,  a_srcData, a_srcDx * 0.5, a_srcDx, a_verbose);
	   FillFromReference(a_destData, stepData, a_destDx, a_srcDx*0.5 , a_verbose);


	 }
       else
	 {
	   MayDay::Error("FillFromReference(LevelData<FarrayBox>& ,const LevelData<FarrayBox>&,... ) odd refinement ratio");
	 }
     }
  else if (refRatio < 1-tolerance)
    {
      if (a_verbose)
        {
          pout() << "averaging data" << endl;
        }
      // need to average data
      refRatio = 1.0/refRatio;
      int nRef = (int) (refRatio + tolerance);
      	{
          if (a_verbose)
            {
              pout() << " ...averaging data with refinement ratio = " 
                     << nRef << endl;
            }
	}
	CoarseAverage averager(srcGrids,destGrids,
			       a_destData.nComp(), nRef, a_destData.ghostVect());
	averager.averageToCoarse(a_destData, a_srcData);
    }
  else
    {
      // same size
      if (a_verbose)
	{
	  pout() << " ...same-size copy of data" << endl;
	}
      a_srcData.copyTo(a_destData);
    }

  a_destData.exchange();

  // 
  if (a_verbose)
    {
      pout() << "Leaving FillFromReference" << endl;
    }

}



#include "NamespaceFooter.H"
