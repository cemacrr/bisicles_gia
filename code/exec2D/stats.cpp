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
   
    if(argc < 5) 
      { 
	std::cerr << " usage: " << argv[0] << " <plot file> <ice_density> <water_density> <gravity> [mask_file] [mask_no = 0] " << std::endl; 
	exit(0); 
      }
    char* plotFile = argv[1];
    Real iceDensity = atof(argv[2]);
    Real waterDensity = atof(argv[3]);
    Real gravity = atof(argv[4]);
    char* maskFile = (argc > 5)?argv[5]:NULL;
    int maskNo = 0;
    if (maskFile && argc > 6)
      {
	maskNo = atoi(argv[6]);
      }

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    Real dt ,crseDx, time;
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
    if (maskFile)
      {
	Box mdomainBox;
	Vector<std::string> mname;
	Vector<LevelData<FArrayBox>* > mdata;
	Vector<DisjointBoxLayout > mgrids;
	Vector<int > mratio;
	int mnumLevels;
	Real mdt ,mcrseDx, mtime;
	ReadAMRHierarchyHDF5(std::string(maskFile), mgrids, mdata, mname , 
			 mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);

	CH_assert(mnumLevels = 1);
	for (int lev =0; lev < numLevels; lev++)
	  {
	    LevelData<FArrayBox> levelMask(grids[lev],1,IntVect::Zero);
	    FillFromReference(levelMask, *mdata[0], vdx[lev], RealVect::Unit*mcrseDx,true);
	    LevelData<FArrayBox>& levelData = *data[lev];
	    for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	      {
		const Box& b = grids[lev][dit];
		FArrayBox coef(b,1); coef.setVal(0.0);
		for (BoxIterator bit(b); bit.ok(); ++bit)
		  {
		    const IntVect& iv = bit();
		    
		    if ( std::abs ( levelMask[dit](iv) - maskNo) < 1.0e-6)
		      {
			coef(iv) = 1.0;
		      }
		    
		  }
		
		for (int i =0; i < levelData.nComp(); ++i)
		  levelData[dit].mult(coef,0,i,1);

	      }

	  }
       
      }


    //extract the topography and thickness fields
    Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
    IntVect sigmaCSGhost(2*IntVect::Unit);

    for (int lev = 0; lev < numLevels; lev++)
      {
	coords[lev] = RefCountedPtr<LevelSigmaCS> 
	  (new LevelSigmaCS(grids[lev], vdx[lev], sigmaCSGhost));
	coords[lev]->setIceDensity(iceDensity);
	coords[lev]->setWaterDensity(waterDensity);
	coords[lev]->setGravity(gravity);
	//coords[lev]->setFaceSigma(faceSigma);

	if (lev > 0)
	  {
	    coords[lev]->interpFromCoarse(*coords[lev-1],ratio[lev-1]);
	  }
	data[lev]->exchange();
	for (int j = 0; j < name.size(); j++)
	  {
	    if (name[j] == "Z_base")
	      {
		data[lev]->copyTo(Interval(j,j),coords[lev]->getTopography(),Interval(0,0));
	      }
	    else if (name[j] == "thickness")
	      {
		data[lev]->copyTo(Interval(j,j),coords[lev]->getH(),Interval(0,0));
	      }
	    
	  }
	if (lev > 0)
	  {
	    coords[lev]->recomputeGeometry(coords[lev-1],ratio[lev-1]);
	  }
	else
	  {
	    coords[lev]->recomputeGeometry(NULL,0);
	  }

      }

    



    pout() << " time = " << time << endl;

    {
      //Compute the total thickness
      Vector<LevelData<FArrayBox>* > thk(numLevels, NULL);
      for (int lev=0; lev< numLevels; lev++)
	{
	  thk[lev] = const_cast<LevelData<FArrayBox>*>(&coords[lev]->getH());
	}
      Real iceVolumeAll = computeSum(thk, ratio, dx[0], Interval(0,0), 0);
      pout() << " iceVolumeAll = " << iceVolumeAll << endl;
    }
    
    
    {
      //Compute the total thickness above flotation
      Vector<LevelData<FArrayBox>* > thk(numLevels, NULL);
      for (int lev=0; lev< numLevels; lev++)
	{
	  thk[lev] = const_cast<LevelData<FArrayBox>*>(&coords[lev]->getThicknessOverFlotation());
	}
      Real iceVolumeAbove = computeSum(thk, ratio, dx[0], Interval(0,0), 0);
      pout() << " iceVolumeAbove = " << iceVolumeAbove << endl;
    }


		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
