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
// nctoamr.cpp
// read in data from a netcdf file, write out a single level AMR hierarchy
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LoadBalance.H"
#include "fabncio.H"
#include "FineInterp.H"
#include "CoarseAverage.H"


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
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> <var 1> [<var 2>, ...] " 
		  << std::endl; 
	exit(0); 
      }
    
    char* in_file = argv[1];
    char* out_file = argv[2];
    
    Vector<std::string> var;
    for (int i = 3; i < argc; i++)
      {
	var.push_back(std::string(argv[i]));
      }  
    
    //pout() << "converting netcdf file " << in_file << " to AMR file " << out_file << std::endl;
    
    FArrayBox fab;
    Box box;
    Real dx;
    if (procID() == uniqueProc(SerialTask::compute))
      {
#ifdef HAVE_NETCDF
	NCIO::readFAB(in_file,var,fab,dx);
#else
	MayDay::Error("netcdf input requested but netcdf support not built")
#endif
	box = fab.box();
      }// end if serial compute
    broadcast(box,uniqueProc(SerialTask::compute));
    ProblemDomain pd(box);
    Vector<Box> boxes(1,pd.domainBox());
    Vector<int> procAssign(1,uniqueProc(SerialTask::compute));
    DisjointBoxLayout grids(boxes, procAssign, pd);
    LevelData<FArrayBox> levelData(grids,var.size(),IntVect::Zero);

    for (DataIterator dit(grids); dit.ok(); ++dit)
      {
	levelData[dit].copy(fab,0,0,fab.nComp());
      }


    Vector<LevelData<FArrayBox>* > vectData(1,&levelData);
    Vector<DisjointBoxLayout > vectGrids(1,grids);
    Vector<int > vectRatio;
    
    const Real dt = 1.0;
    const Real time = 0.0;
    WriteAMRHierarchyHDF5(std::string(out_file), vectGrids, vectData, var , 
			  pd.domainBox(), dx, dt, time, vectRatio, 1);
    
    
  
    
  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
