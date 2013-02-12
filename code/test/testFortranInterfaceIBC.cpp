#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;


#include "parstream.H"

#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "FABView.H"
#include "CONSTANTS.H"

#include "FortranInterfaceIBC.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testFortranInterfaceIBC" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
//static bool verbose = false ;
static bool verbose = true ;

#ifdef CH_USE_DOUBLE
//static Real precision = 1.0e-15;
#else
//static Real precision = 1.0e-7;
#endif

//static Real rateTolerance = 0.02;

// probtype enum
enum basalEnum{sinusoidalZb = 0,
               num_probtype};

enum thicknessTypeEnum{sinusoidalH = 0,
                       num_thicknesstype};

int basalType = sinusoidalZb;

int thicknessType = sinusoidalH;


///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for( int i = 1 ; i < argc ; ++i )
    {
      if( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              //argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return;
}


void initData(LevelData<FArrayBox>& a_thickness,
              LevelData<FArrayBox>& a_topography,
              const RealVect& a_dx,
              const RealVect& a_domainSize)
{
  DataIterator dit = a_thickness.dataIterator(); 
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& topoFab = a_topography[dit];
      FArrayBox& thickFab = a_thickness[dit];
      BoxIterator ccBit(thickFab.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (ccBit.begin(); ccBit.ok(); ++ccBit)
        {
          IntVect iv = ccBit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= a_dx;
          
          if (thicknessType == sinusoidalH)
            {
              Real thisH = 0.5*D_TERM(sin(2*Pi*mappedLoc[0]/a_domainSize[0]),
                                      *sin(2*Pi*mappedLoc[1]/a_domainSize[1]),
                                      NOT_IMPLEMENTED_FOR_3D_YET);
              thickFab(iv,0) = thisH;
            }
          
          
          if (basalType == sinusoidalZb)
            {
              Real thisZb = D_TERM(cos(2*Pi*mappedLoc[0]/a_domainSize[0]),
                                   *cos(2*Pi*mappedLoc[1]/a_domainSize[1]),
                                   NOT_IMPLEMENTED_FOR_3D_YET);
              topoFab(iv,0) = thisZb;
            }
        }
    }
  
}


int
testFortranInterfaceIBC();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testFortranInterfaceIBC();

  if( status == 0 )
    pout() << indent << "FortranInterfaceIBC test" << " passed." << endl ;
  else
    pout() << indent << "FortranInterfaceIBC test" << " failed with return code "
         << status << endl ;



  if( status == 0 )
    pout() << indent << pgmname << " passed." << endl ;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testFortranInterfaceIBC()
{
  int status = 0;

  // this is the resolution of the "fortran" mesh...  
  IntVect numCells(32*IntVect::Unit);
  if (CH_SPACEDIM == 3)
    {
      numCells[0] = 8;
    }
  IntVect hiVect = numCells - IntVect::Unit;
  Box domainBox(IntVect::Zero, hiVect);
  ProblemDomain entireDomain(domainBox);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      entireDomain.setPeriodic(dir,true);
    }

  RealVect domainSize(10.0*RealVect::Unit);
  if (CH_SPACEDIM == 3)
    {
      domainSize[0] = 1.0;
    }
  RealVect dx = domainSize;
  dx /= numCells;

  //int maxBoxSize = 8;
  int numSplit = int(sqrt(float(numProc())));

  int maxBoxSize = 32/numSplit;

  pout() << "numSplit = " << numSplit
         << ", maxBoxSize = " << maxBoxSize << endl;

  Vector<Box> gridBoxes;
  domainSplit(domainBox, gridBoxes, maxBoxSize); 
  
  Vector<int> procAssign(gridBoxes.size(), 0);
  LoadBalance(procAssign, gridBoxes);

  // first generate data
  IntVect baseGhostVect = IntVect::Zero;
  DisjointBoxLayout sourceGrids(gridBoxes, procAssign, entireDomain);
  LevelData<FArrayBox> baseThickness(sourceGrids, 1, baseGhostVect );
  LevelData<FArrayBox> baseZb(sourceGrids, 1, baseGhostVect);

  initData(baseThickness, baseZb,
           dx, domainSize);


  // now construct a FortranInterfaceIBC
  FortranInterfaceIBC baseIBC;

  int diminfo[4];
  diminfo[0] = SpaceDim;
  D_TERM(diminfo[2] = numCells[0];,
         diminfo[3] = numCells[1];,
         diminfo[1] = numCells[2]);


  // second dx is a bit of a kluge so that this will work in 1D
  // assume that there is only one box per processor.
  DataIterator dit = baseThickness.dataIterator();
  dit.begin();
  // grab the box on this processor
  const Box& gridBox = sourceGrids[dit];
  int ewlb = gridBox.loVect()[0];
  int ewub = gridBox.hiVect()[0];
  int nslb = gridBox.loVect()[1];
  int nsub = gridBox.hiVect()[1];

  int lb[SpaceDim];
  int ub[SpaceDim];
  D_TERM(lb[0] = ewlb;,
         lb[1] = nslb;,
         lb[2] = numCells[2];)

  D_TERM(ub[0] = ewub;,
         ub[1] = nsub;,
         ub[2] = numCells[2];)

  baseIBC.define(entireDomain, dx[0]);

  IntVect offset = IntVect::Zero;

  baseIBC.setThickness(baseThickness[dit].dataPtr(),
                       diminfo,
                       lb, ub,
                       &dx[0], D_SELECT(&dx[0],&dx[1],&dx[1]),
                       offset);

  baseIBC.setTopography(baseZb[dit].dataPtr(),
                        diminfo,
                        lb, ub,
                        &dx[0], D_SELECT(&dx[0],&dx[1],&dx[1]),
                        offset);


 
  // now set up grid hiearchy, starting with level 0 one level coarser 
  int numLevels = 4;
  
  ProblemDomain levelDomain = entireDomain;
  levelDomain.coarsen(2);


  RealVect dxLevel(dx);
  dxLevel *= 2;

  ProblemDomain level0domain = levelDomain;
  RealVect dx0 = dxLevel;
  
  Vector<DisjointBoxLayout> vectGrids(numLevels);
  Vector<int> vectNref(numLevels, 2);
  Vector<LevelSigmaCS* > vectCS(numLevels, NULL);
  for (int lev=0; lev<numLevels; lev++)
    {
      IntVect domLo = levelDomain.domainBox().smallEnd();
      IntVect domHi = levelDomain.domainBox().bigEnd();
      IntVect domSize = levelDomain.domainBox().size();
      
      Vector<Box> levelBoxes;
      if (lev == 0)
        {      
          // level 0 spans domain 
          levelBoxes.resize(1);
          levelBoxes[0] = levelDomain.domainBox();
        }
      else if (lev == 1)
        {
          levelBoxes.resize(1);

          IntVect boxlo = domLo + domSize/8;
          IntVect boxhi = ((domHi+1) - domSize/8)-1;
          Box gridBox(boxlo, boxhi);
          levelBoxes[0] = gridBox;
        }
      else if (lev == 2)
        {
          levelBoxes.resize(1);
          IntVect boxlo = domLo + domSize/4;
          IntVect boxhi = ((domHi+1) - domSize/4)-1;
          Box gridBox(boxlo, boxhi);
          levelBoxes[0] = gridBox;
        }
      else if (lev == 3)
        {
          levelBoxes.resize(1);
          IntVect boxlo = domLo + 3*domSize/8;
          IntVect boxhi = ((domHi+1) - 3*domSize/8)-1;
          Box gridBox(boxlo, boxhi);
          levelBoxes[0] = gridBox;
        }

      Vector<int> procAssign(levelBoxes.size(),0);
      LoadBalance(procAssign, levelBoxes);

      DisjointBoxLayout levelGrids(levelBoxes, procAssign, levelDomain);

      vectGrids[lev] = levelGrids;
      IntVect ghostVect = IntVect::Zero;
      vectCS[lev] = new LevelSigmaCS(levelGrids, dxLevel, ghostVect);
      
      IceThicknessIBC* levelIBC = baseIBC.new_thicknessIBC();
      levelIBC->define(levelDomain, dxLevel[0]);

      Real time = 0.0;
      LevelSigmaCS* crseCoords = NULL;
      int nRefCrse = -1;
      levelIBC->initializeIceGeometry(*vectCS[lev], dxLevel, domainSize,
                                      time, crseCoords, nRefCrse );

      delete levelIBC;

      // set up for next level
      levelDomain.refine(vectNref[lev]);
      dxLevel /= vectNref[lev];

    } // end loop over levels
          

  
  bool writePlotFiles = true;
  if (writePlotFiles)
        {
          int nPlotComp = 2;
          Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);
          for (int lev=0; lev<numLevels; lev++)
            {
              const DisjointBoxLayout& levelGrids = vectGrids[lev];
              plotData[lev] = new LevelData<FArrayBox>(levelGrids, 
                                                       nPlotComp,
                                                       IntVect::Zero);

              LevelSigmaCS& levelCS = *vectCS[lev];
              const LevelData<FArrayBox>& levelH = levelCS.getH();
              const LevelData<FArrayBox>& levelZb = levelCS.getBaseHeight();


              // fill plot data
              DataIterator dit = levelGrids.dataIterator();
              for (dit.begin(); dit.ok(); ++dit)
                {

                  FArrayBox& thisData = (*plotData[lev])[dit];
                  const FArrayBox& thisH = levelH[dit];
                  const FArrayBox& thisZb = levelZb[dit];

                  thisData.copy(thisH,0,0,1);
                  thisData.copy(thisZb,0,1,1);
                }
            } // end loop over levels

          Vector<string> plotNames(nPlotComp);
          plotNames[0] = "thickness";
          plotNames[1] = "topography";

          string fname;
          if (SpaceDim == 1 )
            {
              fname = "fortranInterface.1d.hdf5";
            }
          else if (SpaceDim == 2 )
            {
              fname = "fortranInterface.2d.hdf5";
            }
          else if (SpaceDim == 3)
            {
              fname = "fortranInterface.3d.hdf5";
            }
          else
            {
              MayDay::Error("unknown dimensionality");
            }
        
          Real fakeDt = 1.0;
          Real fakeTime = 1.0;
          WriteAMRHierarchyHDF5(fname, vectGrids, 
                                plotData,
                                plotNames,
                                level0domain.domainBox(),
                                dx0[0], fakeDt, fakeTime, 
                                vectNref, 
                                numLevels);

          

          // clean up memory
          for (int lev=0; lev<plotData.size(); lev++)
            {
              delete plotData[lev];
              plotData[lev] = NULL;
            }

        } // end if write plot files
 
  // clean up memory
  for (int lev=0; lev<vectCS.size(); lev++)
    {
      if (vectCS[lev] != NULL)
        {
          delete vectCS[lev];
          vectCS[lev] = NULL;
        }
    }
  return status;
      
}
