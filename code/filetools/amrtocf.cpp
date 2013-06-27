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
// amrtocf.cpp
// Read data from hdf5 files containing Chombo block-structured AMR hierachies
// Write unstructurd data *plus* the grid data needed to reconstruct
// the block structured data to a CF compliant netcdf file
// Or the other way round
//===========================================================================

#include <iostream>
#include "AMRIO.H"
#include "ValidIO.H"


void AMRtoCF(const std::string& ifile, const std::string& ofile);
void CFtoAMR(const std::string& ifile, const std::string& ofile);




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
  

    if(argc < 3) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> " 
		  << std::endl; 
	exit(0); 
      }
    
    char* in_file = argv[1];
    char* out_file = argv[2];
    
  
    std::string ifile(in_file);
    std::string ofile(out_file);
 
    if ( (ifile.size() >= 5) && (ifile.find(".hdf5") == ifile.size()-5))
      {
	AMRtoCF(ifile,ofile);
      }
    else if ((ifile.size() >= 3) && (ifile.find(".nc") == ifile.size()-3))
      {
        CFtoAMR(ifile,ofile);
      }
    else
      {
	pout() << "unknown input file type [" <<  ifile << "]" << std::endl;
	MayDay::Error("unknown input file type");
      }


   

   
    
  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
void AMRtoCF(const std::string& ifile, const std::string& ofile)
{
  

  //read the AMR data
  Vector<std::string> names;
  Vector<LevelData<FArrayBox>* > data;
  Vector<DisjointBoxLayout> grids;
  Vector<int> ratio;
  Real crseDx = 0.0, dt = 0.0, time = 0.0;
  Box domBox;
  int numLevels;
  int status = ReadAMRHierarchyHDF5
    (ifile,grids,data,names,domBox,crseDx,dt,time,
     ratio,numLevels);
  if (status != 0)
    {
      MayDay::Error("failed to read AMR hierarchy");
    }

  Vector<RealVect> dx(numLevels);
  dx[0] = crseDx * IntVect::Unit;
  for (int lev = 1; lev < numLevels; lev++)
    {
      dx[lev] = dx[lev-1]/Real(ratio[lev-1]);
    }

  RealVect x0 = RealVect::Unit * -434.123e+3; //\todo  NEEDS to be a sensible parameter

  //extract valid data
  ValidData validData(names.size(), dx, x0);
  ValidIO::BStoValid(validData, data, ratio);

  
  PolarStereographicCartesianToPolarTransformation transformation
    (0.08181922,  6.3781370e+6, 1.0 , 0.0, RealVect::Zero);
  //eccentricity, equatorial radius, 

  //write to CF file
  ValidIO::writeCF ( ofile,  validData , names, transformation );



  for (int lev = 0; lev < numLevels; lev++)
    {
      if (data[lev] != NULL)
	{
	  delete data[lev];data[lev]=NULL;
	}
    }
}
void CFtoAMR(const std::string& ifile, const std::string& ofile)
{
 MayDay::Error("CFtoAMR not implemented");
}


