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
// flatten.cpp
// read a multi-level AMR hierarchy and write out a single level AMR hierachy
// as either a Chombo hdf5 file or a netcdf file
//===========================================================================

#include <iostream>
#include "AMRIO.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "fabncio.H"
#include "NamespaceHeader.H"

bool verbose = true;

enum out_file_type_enum {hdf5,nc};

int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

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

    if(argc < 4) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> level [x0 [y0 [z0]]]" << std::endl; 
	exit(0); 
      }

    char* in_file = argv[1];
    char* out_file = argv[2];

 
    out_file_type_enum out_file_type;
    std::string ofile(out_file);
    if (ofile.find(".hdf5") == ofile.size()-5)
      {
	out_file_type = hdf5;
      }
    else if (ofile.find(".nc") == ofile.size()-3)
      {
	out_file_type = nc;
      }
    else
      {
	pout() << "unknown output file type [" <<  ofile << "]" << std::endl;
	MayDay::Error("unknown output file type");
      }

    int flatLevel = atoi(argv[3]);

    RealVect x0 = RealVect::Zero;
    
    if (argc == 4+SpaceDim)
      {
	for (int dir=0;dir < SpaceDim; dir++)
	  x0[dir] = atof(argv[4+dir]);
      }

    
    //pout() << "flattening AMR file " << in_file << " on to level " 
    //	   << flatLevel << " and writing to " << out_file << std::endl;

    Vector<std::string> names;

    if (verbose)
      {
        pout ();
        pout() << "reading AMR file..." << endl;
      }
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout> grids;
    Vector<int> ratio;
    Real crseDx = 0.0, dt = 0.0, time = 0.0;
    Box crseBox;
    int numLevels;
    int status = ReadAMRHierarchyHDF5
      (in_file,grids,data,names,crseBox,crseDx,dt,time,
       ratio,numLevels);
    if (status != 0)
      {
        MayDay::Error("failed to read AMR hierarchy");
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
      }

    Box flatBox(crseBox);
    Real flatDx = crseDx;
    if (flatLevel >= 0)
      {
        for (int lev=0; lev < flatLevel;lev++)
          {
            flatBox.refine(ratio[lev]);
            flatDx /= Real(ratio[lev]);
          }
      }
    else
      {
        // if flatlevel < 0, then coarsen level 0 domain by a factor of 2 
        // for each level
        for (int lev=0; lev>flatLevel; lev--)
          {
            flatBox.coarsen(2);
            flatDx *= 2.0;
          }
      }

    ProblemDomain pd(flatBox);
    Vector<Box> boxes(1,flatBox);
    Vector<int> procAssign(1,uniqueProc(SerialTask::compute));
    DisjointBoxLayout flatDBL(boxes, procAssign, pd);
    LevelData<FArrayBox> flatLevelData(flatDBL,names.size(),IntVect::Unit);
    LevelData<FArrayBox> fiData;
    int nRef;
    int nComp = names.size();
    Interval ivl(0,nComp-1);
    
    //interpolate from coarser levels
    for (int lev=0; lev < flatLevel;lev++)
      {
	nRef = 1;
	for (int l=lev; l < flatLevel;l++)
	  {
	    nRef *= ratio[l];
	  }

        const DisjointBoxLayout& fineGrids = flatLevelData.getBoxes();
        FineInterp fi(fineGrids,nComp,nRef,fineGrids.physDomain());
        fi.interpToFine(flatLevelData,*data[lev], true);
      }

    //direct copy on same level
    if (flatLevel >= 0)
      {
        data[flatLevel]->copyTo(ivl,flatLevelData,ivl);
      }
    
    //average from finer levels
    nRef = 1;
    int startLevel = max(0,flatLevel +1);
    if (flatLevel < 0)
      {
        for (int lev=flatLevel; lev<startLevel; lev++)
          {
            nRef *= 2;
          }
      }
    for (int lev = startLevel ; lev < data.size() ; lev++)
      {
        if (lev > 0)
          {
            nRef *= ratio[lev-1];
          }
	CoarseAverage av(grids[lev], nComp, nRef);
	av.averageToCoarse(flatLevelData,*data[lev]);		 
      }

    if (out_file_type == hdf5)
      {
	Vector<int> one(1,2);
	Vector<DisjointBoxLayout> fdbl(1,flatDBL);
	Vector<LevelData<FArrayBox>* > fdata(1,&flatLevelData);
	WriteAMRHierarchyHDF5(ofile, fdbl, fdata, names, 
			      flatBox, flatDx, dt, time, one, 1);
      }
    else if (out_file_type == nc)
      {
	if (procID() == uniqueProc(SerialTask::compute))
	  {
	    DataIterator dit(flatDBL); 
	    const FArrayBox& fab = flatLevelData[dit];

	    // Looking for illegal / character in variables names
	    for (int iname=0; iname < names.size();iname++)
	      {
		int locslash = names[iname].find("/"); 
		if (locslash <= names[iname].size())
		  {
		    pout() << "About to write netCDF file..." << endl;
		    pout() << "Changing variable name " <<  names[iname] << std::endl;
		    names[iname].replace(locslash,1,"_");
		    pout() << "... to " <<  names[iname] << std::endl;	    
		  }
	      }
            
            // remove ghost cells from netcdf data, since they're 
            // going to be unset anyway
            Box gridBox = flatLevelData.getBoxes()[dit];
            FArrayBox validFab(gridBox, fab.nComp());
            validFab.copy(fab);
            
	    writeNetCDF(out_file, names, validFab, flatDx, x0);
	  }
      }
   
  
    // clean up memory
    for (int lev=0; lev<data.size(); lev++)
      {
        if (data[lev] != NULL)
          {
            delete data[lev];
            data[lev] = NULL;
          }
      }

  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
#include "NamespaceFooter.H"
