#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//===========================================================================
// stats.cpp
// read in a bisicles plotfile (which must include thickness and Z_base) 
// and an optional mask, and write out a bunch of stats about the ice sheet. 
// These are
// 0. time
// 1. Volume of ice
// 2. Volume of ice above flotation
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LevelSigmaCS.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "FillFromReference.H"

void createMaskedDEM( Vector<LevelData<FArrayBox>* >& topography,  
		      Vector<LevelData<FArrayBox>* >& thickness, 
		      Vector<LevelData<FArrayBox>* >& mdata, 
		      Vector<std::string>& name, 
		      Vector<LevelData<FArrayBox>* >& data,
		      Vector<Real>& dx,
		      Real mcrseDx,
		      int maskNo)
{
  CH_TIME("createMaskedDEM");
  int numLevels = data.size();
  for (int lev = 0; lev < numLevels; lev++)
    {
      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "Z_base")
	    {
	      data[lev]->copyTo(Interval(j,j),*topography[lev],Interval(0,0));
	      topography[lev] -> exchange();
	    }
	  else if (name[j] == "thickness")
	    {
	      data[lev]->copyTo(Interval(j,j),*thickness[lev],Interval(0,0));
	      thickness[lev] -> exchange();
	    }	  
	}
    }
  if (mdata.size() != 0)
    {
      // there is mask data, so set thickness to zero outside the selected region
      for (int lev =0; lev < numLevels; lev++)
	{
	  const DisjointBoxLayout& levelGrids = topography[lev]->disjointBoxLayout();
	  LevelData<FArrayBox> levelMask(levelGrids,1,IntVect::Zero);
	  FillFromReference(levelMask, *mdata[0], RealVect::Unit*dx[lev], RealVect::Unit*mcrseDx,true);
	  LevelData<FArrayBox>& levelData = *thickness[lev];
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      const Box& b = levelGrids[dit];
	      FArrayBox coef(b,1); coef.setVal(0.0);
	      for (BoxIterator bit(b); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		    
		  if ( std::abs ( levelMask[dit](iv) - maskNo) < 1.0e-6)
		    {
		      coef(iv) = 1.0;
		    }
		    
		}
	      levelData[dit].mult(coef,0,0,1);
		
	    }
	    
	}
    }
}

void computeStats(Vector<LevelData<FArrayBox>* >& topography, 
		  Vector<LevelData<FArrayBox>* >& thickness, Vector<Real>& dx, Vector<int>& ratio,
		  Real iceDensity, Real waterDensity, Real gravity)
{ 

  CH_TIMERS("computeStats");
  CH_TIMER("createSigmaCS",t1);
  CH_TIMER("integrateH",t2);
  CH_TIMER("integrateHab",t3);
  CH_TIMER("integrateGA",t4);
  
  int numLevels = topography.size();

  CH_START(t2);
   
  //Compute the total thickness
  // for (int lev=0; lev< numLevels; lev++)
  //  {
  //    coords[lev]->getH().copyTo(Interval(0,0),*tmp[lev],Interval(0,0));
  //  }
  Real iceVolumeAll = computeSum(thickness, ratio, dx[0], Interval(0,0), 0);
 
  CH_STOP(t2);
  Real iceVolumeAbove = 0.0;
  Real groundedArea = 0.0;

  //Creating a LevelSigmaCS is expensive, so only do it if there is some ice
  if (iceVolumeAll > 1.0e-10)
    {
      CH_START(t1);
      
      Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
      IntVect sigmaCSGhost(2*IntVect::Unit);
       
      for (int lev = 0; lev < numLevels; lev++)
	{
	  const DisjointBoxLayout& levelGrids = topography[lev]->disjointBoxLayout();
	   
	  coords[lev] = RefCountedPtr<LevelSigmaCS> 
	    (new LevelSigmaCS(levelGrids, RealVect::Unit*dx[lev], sigmaCSGhost));
	  coords[lev]->setIceDensity(iceDensity);
	  coords[lev]->setWaterDensity(waterDensity);
	  coords[lev]->setGravity(gravity);
	  topography[lev]->copyTo(Interval(0,0),coords[lev]->getTopography(),Interval(0,0));   
	  thickness[lev]->copyTo(Interval(0,0),coords[lev]->getH(),Interval(0,0));
	  coords[lev]->recomputeGeometry(NULL,0);
	   
	   
	}
       
      Vector<LevelData<FArrayBox>* > tmp(numLevels, NULL);
      for (int lev=0; lev< numLevels; lev++)
	{
	  const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
	  tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

	}
      CH_STOP(t1);
       
      CH_START(t3);
       
      {
	//Compute the total thickness above flotation;
	for (int lev=0; lev< numLevels; lev++)
	  {
	    coords[lev]->getThicknessOverFlotation().copyTo(Interval(0,0),*tmp[lev],Interval(0,0));
	  }
	iceVolumeAbove = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
	 
      }
      CH_STOP(t3);
      CH_START(t4);
       
      {
	//grounded area
	 
	for (int lev=0; lev< numLevels; lev++)
	  {
	     
	    const DisjointBoxLayout& grids = coords[lev]->grids();
	    tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
	    for (DataIterator dit(grids);dit.ok();++dit)
	      {
		const BaseFab<int>& mask =  coords[lev]->getFloatingMask()[dit];
		const Box& b = grids[dit];
		FArrayBox& a = (*tmp[lev])[dit];
		a.setVal(0.0);
		for (BoxIterator bit(b);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();
		    if (mask(iv) == GROUNDEDMASKVAL)
		      {
			a(iv) = 1.0;
		      }
		
		  }
	      }
	  }
	groundedArea = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
   
    
      }
      CH_STOP(t4);
      for (int lev=0; lev< numLevels; lev++)
	{
	  if (tmp[lev] != NULL)
	    {
	      delete tmp[lev]; tmp[lev];
	    } 
	}
    }
    
  pout() << " iceVolumeAll = " << iceVolumeAll << " ";
  pout() << " iceVolumeAbove = " << iceVolumeAbove << " ";
  pout() << " groundedArea = " << groundedArea << " ";


 
  
}


int main(int argc, char* argv[]) {

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif
   
    CH_TIMERS("stats");
    CH_TIMER("loadplot",tp);
    if(argc < 5) 
      { 
	std::cerr << " usage: " << argv[0] << " <plot file> <ice_density> <water_density> <gravity> [mask_file] [mask_no_start = 0] [mask_no_end = mask_no_start] " << std::endl; 
	exit(0); 
      }
    char* plotFile = argv[1];
    Real iceDensity = atof(argv[2]);
    Real waterDensity = atof(argv[3]);
    Real gravity = atof(argv[4]);
    char* maskFile = (argc > 5)?argv[5]:NULL;
    int maskNoStart = 0;
    if (maskFile && argc > 6)
      {
    	maskNoStart = atoi(argv[6]);
      }
    int maskNoEnd = maskNoStart;
    if (maskFile && argc > 7)
      {
    	maskNoEnd = atoi(argv[7]);
      }

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    
    Real dt ,crseDx, time;
    CH_START(tp);
    ReadAMRHierarchyHDF5(std::string(plotFile), grids, data, name , 
			 domainBox, crseDx, dt, time, ratio, numLevels);

    Vector<ProblemDomain> domain(numLevels,domainBox);
    Vector<RealVect> vdx(numLevels,RealVect::Unit*crseDx);
    Vector<Real> dx(numLevels,crseDx);
    for (int lev=1;lev<numLevels;++lev)
      {
	dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
	vdx[lev] = vdx[lev-1] / Real(ratio[lev-1]);
	domain[lev] = domain[lev-1];
	domain[lev].refine(ratio[lev-1]);
      }

    //load the sector mask, if it exists
    //Vector<RefCountedPtr<LevelData<FArrayBox> > > sectorMask;
    
    	
    Box mdomainBox;
    Vector<std::string> mname;
    Vector<LevelData<FArrayBox>* > mdata;
    Vector<DisjointBoxLayout > mgrids;
    Vector<int > mratio;
    int mnumLevels;
    Real mdt ,mcrseDx, mtime;

    if (maskFile)
      {
	ReadAMRHierarchyHDF5(std::string(maskFile), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);
      }
    
    CH_STOP(tp);

    Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
    for (int lev=0;lev<numLevels;++lev)
      {
	thickness[lev] = new LevelData<FArrayBox>(grids[lev],1,2*IntVect::Unit);
	topography[lev] = new LevelData<FArrayBox>(grids[lev],1,2*IntVect::Unit);
      }


    pout().setf(ios_base::scientific,ios_base::floatfield); 
    pout().precision(12);

   
    for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
      {
	createMaskedDEM(topography, thickness, mdata, name, data,dx, mcrseDx, maskNo);
	pout() << " time = " << time  ;
	computeStats(topography, thickness, dx, ratio, iceDensity, waterDensity,gravity);
	if (maskFile)

	  pout() << " sector = " << maskNo; 
	pout() << endl;
      }
    
    for (int lev=0;lev<numLevels;++lev)
      {
	if (thickness[lev] != NULL) delete thickness[lev];
	if (topography[lev] != NULL) delete topography[lev];
      }

		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
