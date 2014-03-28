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
// driver.cpp
//
//===========================================================================
#include <iostream>
#include <fstream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "AmrIce.H"
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
#include "BasalFrictionRelation.H"
#include "MuCoefficient.H"
#include "twistyStreamFriction.H"
#include "singularStreamFriction.H"
#include "GaussianBumpFriction.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "HumpIBC.H"
#include "LevelDataIBC.H"
#include "MultiLevelDataIBC.H"
#include "IceTemperatureIBC.H"
#include "LevelDataTemperatureIBC.H"
#include "LevelDataBasalFriction.H"
#include "PiecewiseLinearFlux.H"
#include "SurfaceFlux.H"
#include "IceConstants.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
//#include "LevelDataSurfaceFlux.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "ReadLevelData.H"
#ifdef CH_USE_PETSC
#include "petsc.h"
#endif 

/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
    guassianBump is used for the MISMIP3D perturbations tests.
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
			 gaussianBump,
			 singularStream,
                         NUM_BETA_TYPES};

//===========================================================================
// 2D shallow-shelf ice sheet model
//
//===========================================================================
int main(int argc, char* argv[]) {
  
  int ierr = 0;

#ifdef CH_USE_PETSC
  ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 
#endif // end petsc conditional


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
    
    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);
    ParmParse pp2("main");

    std::string poutBaseName = "pout";
    pp2.query("poutBaseName",poutBaseName);
    setPoutBaseName(poutBaseName);


    RealVect domainSize;
    Vector<Real> domSize(SpaceDim);
    pp2.getarr("domain_size", domSize, 0, SpaceDim);
    domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));


    AmrIce amrObject;
    // ---------------------------------------------
    // set constitutive relation & rate factor
    // ---------------------------------------------
   
    std::string rateFactorType = "constRate";
    pp2.query("rateFactor", rateFactorType);
    if (rateFactorType == "constRate")
      {
	ParmParse crPP("constRate");
	Real A = 9.2e-18;
	crPP.query("A", A);
	ConstantRateFactor rateFactor(A);
	amrObject.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "arrheniusRate")
      {
	ArrheniusRateFactor rateFactor;
	ParmParse arPP("ArrheniusRate");
	amrObject.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "patersonRate")
      {
	PatersonRateFactor rateFactor;
	ParmParse arPP("PatersonRate");
	amrObject.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "zwingerRate")
      {
	ZwingerRateFactor rateFactor;
	ParmParse arPP("ZwingerRate");
	amrObject.setRateFactor(&rateFactor);
      }

    ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");

    if (constRelPtr == NULL)
      {
	MayDay::Error("undefined constitutiveRelation in inputs");
      }

    amrObject.setConstitutiveRelation(constRelPtr);
 
    std::string basalRateFactorType = "";
    pp2.query("basalRateFactor", basalRateFactorType);
    
    if (basalRateFactorType == "patersonRate")
      {
	PatersonRateFactor rateFactor;
	rateFactor.setA0(1.0);
	amrObject.setBasalRateFactor(&rateFactor);
      }

    // ---------------------------------------------
    // set surface flux. 
    // ---------------------------------------------

    SurfaceFlux* surf_flux_ptr = SurfaceFlux::parseSurfaceFlux("surfaceFlux");
    if (surf_flux_ptr == NULL)
      {
	const std::string err("failed to parse surfaceFlux (maybe you have the old style surface_flux_type?");
	pout() << err << endl;
	MayDay::Error(err.c_str());
      }

    amrObject.setSurfaceFlux(surf_flux_ptr);
  
    // ---------------------------------------------
    // set basal (lower surface) flux. 
    // ---------------------------------------------
    
    SurfaceFlux* basal_flux_ptr = SurfaceFlux::parseSurfaceFlux("basalFlux");
    if (basal_flux_ptr == NULL)
      {
	const std::string err("failed to parse basalFlux (maybe you have the old style basal_flux_type?");
	pout() << err << endl;
	MayDay::Error(err.c_str());
      }

    amrObject.setBasalFlux(basal_flux_ptr); 



    // ---------------------------------------------
    // set mu coefficient
    // ---------------------------------------------
    {
      MuCoefficient* muCoefPtr =  MuCoefficient::parseMuCoefficient("muCoefficient");

      if (muCoefPtr == NULL)
	{
	  const std::string err("failed to parse muCoefficient");
	  pout() << err << endl;
	  MayDay::Error(err.c_str());
	}
     
      amrObject.setMuCoefficient(muCoefPtr);
      delete muCoefPtr;
    }

    // ---------------------------------------------
    // set basal friction coefficient and relation
    // ---------------------------------------------

    ParmParse geomPP("geometry");
    
    BasalFriction* basalFrictionPtr 
      = BasalFriction::parseBasalFriction("geometry", domainSize);
    
    if (basalFrictionPtr == NULL)
      {
	MayDay::Error("undefined  geometry.beta_type in inputs");
      }
    
    amrObject.setBasalFriction(basalFrictionPtr);
    
    BasalFrictionRelation* basalFrictionRelationPtr;
    std::string basalFrictionRelType = "powerLaw";
    pp2.query("basalFrictionRelation", basalFrictionRelType);
    
    if (basalFrictionRelType == "powerLaw")
      {
	ParmParse plPP("BasalFrictionPowerLaw");

	Real m = 1.0;
	plPP.query("m",m);
	bool includeEffectivePressure = false;
	plPP.query("includeEffectivePressure",includeEffectivePressure);
	BasalFrictionPowerLaw*  pl = new BasalFrictionPowerLaw(m,includeEffectivePressure);
	basalFrictionRelationPtr = static_cast<BasalFrictionRelation*>(pl);
      }
    else
      {
	MayDay::Error("undefined basalFrictionRelation in inputs");
      }

    amrObject.setBasalFrictionRelation(basalFrictionRelationPtr);

    // ---------------------------------------------
    // set IBC -- this includes initial ice thickness, 
    // and basal geometry
    // ---------------------------------------------

    
    IceThicknessIBC* thicknessIBC;

    std::string problem_type;
    geomPP.get("problem_type", problem_type);
    if (problem_type == "basic")
      {
        thicknessIBC = new BasicThicknessIBC;
      }
    else if (problem_type == "VieliPayne")
      {
        VieliPayneIBC* ibcPtr = new VieliPayneIBC;
        ParmParse pvPP("vieliPayne");

        Real thickness, seaLevel, originElevation;
        RealVect basalSlope;
        pvPP.get("thickness", thickness);
        seaLevel = 0.0;
        pvPP.query("seaLevel", seaLevel);
        
        Vector<Real> vect(SpaceDim);
        pvPP.getarr("basal_slope", vect, 0, SpaceDim);
        basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));

        pvPP.get("originElevation", originElevation);

        ibcPtr->setParameters(thickness, basalSlope, 
                              originElevation, seaLevel);       

        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "marineIceSheet")
      {
        MarineIBC* ibcPtr = new MarineIBC;
        ParmParse mPP("marineIceSheet");
	
	Real  seaLevel;
	seaLevel = 0.0;
	mPP.query("seaLevel", seaLevel);

	std::string thicknessType = "constant";
	mPP.query("thickness_type",thicknessType);
	RefCountedPtr<RealFunction<RealVect> > thicknessFunction;
	Vector<RefCountedPtr<RealFunction<RealVect> > > bedrockFunction(1);
	if (thicknessType == "constant")
	  {
	    Real thickness;
	    mPP.get("thickness", thickness);
	    RefCountedPtr<RealFunction<RealVect> > ptr(new ConstantRealFunction<RealVect>(thickness));
	    thicknessFunction =ptr;
	  }
	else if (thicknessType == "flowline")
	  {
	    Real dx; 
	    std::string file, set;
	    mPP.get("thickness_flowline_dx", dx);
	    mPP.get("thickness_flowline_file", file);
	    mPP.get("thickness_flowline_set", set);
	    RefCountedPtr<RealFunction<RealVect> > ptr(new ExtrudedPieceWiseLinearFlowline(file,set,dx));
	    thicknessFunction = ptr;
	  }
	else 
	  {
	    MayDay::Error("bad marineIceSheet.thicknessType");
	  }

	
	std::string geometry = "plane";
	mPP.query("geometry",geometry);

	if (geometry == "plane")
	  {
	    //inclined plane geometry
	    RealVect basalSlope;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("basal_slope", vect, 0, SpaceDim);
	    basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
	    Real originElevation;
	    mPP.get("originElevation", originElevation);

	   
	    RefCountedPtr<RealFunction<RealVect> > ptr(new InclinedPlaneFunction(originElevation, basalSlope));
	    bedrockFunction[0] =  ptr;
	  }
	else if (geometry == "symmetricPlane")
	  {
	    //inclined plane geometry, symmetric about origin
	    RealVect basalSlope;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("basal_slope", vect, 0, SpaceDim);
	    basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));

	    Real originElevation;
	    mPP.get("originElevation", originElevation);

            RealVect symmetryPoint(RealVect::Zero);
            mPP.getarr("symmetryPoint", vect, 0, SpaceDim);
            symmetryPoint = RealVect(D_DECL(vect[0], vect[1], vect[2]));
	   
	    RefCountedPtr<RealFunction<RealVect> > ptr(new SymmetricInclinedPlaneFunction(originElevation, basalSlope, symmetryPoint));
	    bedrockFunction[0] =  ptr;
	  }

        else if (geometry == "regroundingTest")
          {
	    //inclined plane geometry with a Gaussian bump
            bedrockFunction.resize(2);

	    RealVect basalSlope;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("basal_slope", vect, 0, SpaceDim);
	    basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
	    Real originElevation;
	    mPP.get("originElevation", originElevation);

            // compose flat plane with Gaussian hump
	    RefCountedPtr<RealFunction<RealVect> > ptr1(new InclinedPlaneFunction(originElevation, basalSlope));
	    bedrockFunction[0] =  ptr1;
            
            Real bumpCenter;
            Real bumpRad;
            Real bumpMag;
            
            mPP.get("bumpCenter", bumpCenter);
            mPP.get("bumpRad", bumpRad);
            mPP.get("bumpMag",bumpMag);
            RefCountedPtr<RealFunction<RealVect> > ptr2(new GaussianFunctionX(bumpCenter, bumpRad, bumpMag));

	    bedrockFunction[1] =  ptr2;
          }
	else if (geometry == "Schoof")
	  {
	    //geometry of Schoof, 2007
	    Real originElevation;
	    mPP.get("originElevation", originElevation);
	    
	    Real lengthScaleFactor = 1.0;
	    mPP.query("schoofLengthScaleFactor",lengthScaleFactor);

	   
	    Real schoofCoeff2, schoofCoeff4, schoofCoeff6;
	    mPP.get("schoofCoeff2", schoofCoeff2);
	    mPP.get("schoofCoeff4", schoofCoeff4);
	    mPP.get("schoofCoeff6", schoofCoeff6);
	    
	    //RefCountedPtr<RealFunction<RealVect> > schoofBedrock
	    RefCountedPtr<RealFunction<RealVect> > ptr(new SchoofBedrockElevation(domainSize[SpaceDim-2] * lengthScaleFactor,
										  originElevation,
										  schoofCoeff2, schoofCoeff4, 
										  schoofCoeff6));

	    bedrockFunction[0] = ptr;
	   

	  }

	else if (geometry == "Katz")
	  {
	    //geometry of Katz and Worster, 2010
	    Real originElevation;
	    mPP.get("originElevation", originElevation);
	    
	    Real lengthScaleFactor = 1.0;
	    mPP.query("schoofLengthScaleFactor",lengthScaleFactor);

	    Real katzAlpha, katzSigma;
	    mPP.get("katzAlpha", katzAlpha);
	    mPP.get("katzSigma", katzSigma);

	    Real schoofCoeff2, schoofCoeff4, schoofCoeff6;
	    mPP.get("schoofCoeff2", schoofCoeff2);
	    mPP.get("schoofCoeff4", schoofCoeff4);
	    mPP.get("schoofCoeff6", schoofCoeff6);
	    
	    //RefCountedPtr<RealFunction<RealVect> > katzBedrock
	    RefCountedPtr<RealFunction<RealVect> > ptr(new KatzBedrockElevation(domainSize[SpaceDim-2],
										domainSize[SpaceDim-1],
										originElevation,
										katzAlpha, katzSigma,
										lengthScaleFactor,
										schoofCoeff2, schoofCoeff4, 
										schoofCoeff6));

	    bedrockFunction[0] = ptr;
	  

	  }
	else
	  {
	    MayDay::Error("bad marineIceSheet.geometry");
	  }
	
	ibcPtr->setParameters(thicknessFunction, bedrockFunction ,  seaLevel);
        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "hump")
       {
         HumpIBC* ibcPtr = new HumpIBC;
         ParmParse humpPP("hump");

         Real maxThickness, radSqr, baseElevation, minThickness, seaLevel;
         RealVect center, widthScale;
         
         // default values to be equivalent to hump in Glimmer-CISM
         radSqr = 0.125*domainSize[0]*domainSize[1];
         maxThickness = 2000.0*pow(radSqr,0.5);
         baseElevation = 0.0;
         minThickness = 0.0;
         widthScale = RealVect::Unit;
         // this just lowers the sea level so that it's not relevant...
         seaLevel = -10.0;
         center = 0.5*domainSize;

         humpPP.query("radSqr", radSqr);
         humpPP.query("maxThickness", maxThickness);
         humpPP.query("baseElevation", baseElevation);
         humpPP.query("minThickness", minThickness);
         if (humpPP.contains("center"))
           {
             Vector<Real> centerArr(SpaceDim);
             humpPP.getarr("center", centerArr, 0, SpaceDim);
             center = RealVect(D_DECL(centerArr[0], centerArr[1], 
                                      centerArr[2]));
           }

         if (humpPP.contains("widthScale"))
           {
             Vector<Real> factorArr(SpaceDim);
             humpPP.getarr("widthScale", factorArr, 0, SpaceDim);
             widthScale = RealVect(D_DECL(factorArr[0], factorArr[1], 
                                           factorArr[2]));
           }

         ibcPtr->setParameters(maxThickness,
                               radSqr,
                               baseElevation,
                               minThickness,
                               center,
                               seaLevel, 
                               widthScale);
         
         thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
       }
     else if (problem_type == "LevelData")
       {
	 //read geometry from an AMR Hierarchy, store in LevelDataIBC
	 ParmParse ildPP("inputLevelData");
	 std::string infile;
	 ildPP.get("geometryFile",infile);
	 std::string thicknessName = "thck";
	 ildPP.query("thicknessName",thicknessName);
	 std::string topographyName = "topg";
	 ildPP.query("topographyName",topographyName);
	 
	 RefCountedPtr<LevelData<FArrayBox> > levelThck
	   (new LevelData<FArrayBox>());
	 RefCountedPtr<LevelData<FArrayBox> > levelTopg
	   (new LevelData<FArrayBox>());

	 Real dx;

	 Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
	 vectData.push_back(levelThck);
	 vectData.push_back(levelTopg);

	 Vector<std::string> names(2);
	 names[0] = thicknessName;
	 names[1] = topographyName;
	 readLevelData(vectData,dx,infile,names,1);

         // this is about removing ice from regions which
         // don't affect the dynamics of the region, but which 
         // can cause our solvers problems. 
         // for now, store these regions in a separate file. 
         // There's probably a better way to do this.
         
         // this will contain the boxes in the index space of the 
         // original LevelData in which the thickness will be cleared.
         
         if (ildPP.contains("clearThicknessRegionsFile"))
           {
             Vector<Box> clearBoxes;
             std::string clearFile;
             ildPP.get("clearThicknessRegionsFile", clearFile);
             
             if (procID() == uniqueProc(SerialTask::compute))
               {
                 ifstream is(clearFile.c_str(), ios::in);
                 if (is.fail())
                   {
                     MayDay::Error("Cannot open file with regions for thickness clearing");
                   }
                 // format of file: number of boxes, then list of boxes.
                 int numRegions;
                 is >> numRegions;
                 
                 // advance pointer in file
                 while (is.get() != '\n');
                 
                 clearBoxes.resize(numRegions);
                 
                 for (int i=0; i<numRegions; i++)
                   {
                     Box bx;
                     is >> bx;
                     while (is.get() != '\n');
                     
                     clearBoxes[i] = bx;
                   }
                 
               } // end if serial proc
             // broadcast results
             broadcast(clearBoxes, uniqueProc(SerialTask::compute));
             
             // now loop over the thickness levelData and set intersections
             // with boxes to zero
             
             DataIterator dit = levelThck->dataIterator();
             for (dit.begin(); dit.ok(); ++dit)
               {
                 FArrayBox& thickFab = levelThck->operator[](dit);
                 const Box& fabBox = thickFab.box();
                 for (int boxno=0; boxno<clearBoxes.size(); boxno++)
                   {
                     Box intersectBox(fabBox);
                     intersectBox &= clearBoxes[boxno];
                     if (!intersectBox.isEmpty())
                       {
                         thickFab.setVal(0.0,intersectBox,0);
                       } // end if there's an intersection
                   } // end loop over clearboxes
               } // end loop over grids in thickness levelData

           } // end if we're setting thickness to zero
       
	 RealVect levelDx = RealVect::Unit * dx;
	 LevelDataIBC* ptr = new LevelDataIBC(levelThck,levelTopg,levelDx);
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
     else if (problem_type == "MultiLevelData")
       {
	 //read geometry from an AMR Hierarchy, store in MultiLevelDataIBC
	 ParmParse ildPP("inputLevelData");
	 std::string infile;
	 ildPP.get("geometryFile",infile);
	 std::string thicknessName = "thck";
	 ildPP.query("thicknessName",thicknessName);
	 std::string topographyName = "topg";
	 ildPP.query("topographyName",topographyName);
	
	
	 Real dx;
	 Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectData;
	

	 Vector<std::string> names(2);
	 names[0] = thicknessName;
	 names[1] = topographyName;
	 Vector<int> refRatio;
	 readMultiLevelData(vectData,dx,refRatio,infile,names,1);
	
	 RealVect crseDx = RealVect::Unit * dx;
	 MultiLevelDataIBC* ptr = new MultiLevelDataIBC
	   (vectData[0],vectData[1],crseDx,refRatio);
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);

       }
#ifdef HAVE_PYTHON
     else if (problem_type == "Python")
       {
	 
	 ParmParse pyPP("PythonIBC");
	 std::string module;
	 pyPP.get("module",module);
	 std::string thckFuncName = "thickness";
	 pyPP.query("thicknessFunction",thckFuncName);
	 std::string topgFuncName = "topography";
	 pyPP.query("topographyFunction",topgFuncName);
	 std::string rhsFuncName = "";
	 pyPP.query("RHSFunction",rhsFuncName);
	 std::string faceVelFuncName = "";
	 pyPP.query("faceVelFunction",faceVelFuncName);
	 PythonInterface::PythonIBC* ptr = new PythonInterface::PythonIBC
	   (module, thckFuncName, topgFuncName, rhsFuncName,faceVelFuncName);
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
#endif
     else 
       {
         MayDay::Error("bad problem type");
       }

    amrObject.setThicknessBC(thicknessIBC);

    IceTemperatureIBC* temperatureIBC = NULL;
    ParmParse tempPP("temperature");
    std::string tempType("constant");
    tempPP.query("type",tempType);
    if (tempType == "constant")
      {
	Real T = 258.0;
	tempPP.query("value",T);
	ConstantIceTemperatureIBC* ptr = new ConstantIceTemperatureIBC(T);
	temperatureIBC  = static_cast<IceTemperatureIBC*>(ptr);
      }
    else if (tempType == "LevelData")
      {
	ParmParse ildPP("inputLevelData");
	LevelDataTemperatureIBC* ptr = NULL;
	ptr = LevelDataTemperatureIBC::parse(ildPP); CH_assert(ptr != NULL);
	temperatureIBC  = static_cast<IceTemperatureIBC*>(ptr);
      }
#ifdef HAVE_PYTHON
    else if (tempType == "Python")
      {
	ParmParse pyPP("PythonIceTemperatureIBC");
	std::string module;
	pyPP.get("module",module);
	std::string funcName = "temperature";
	pyPP.query("function",funcName);
	temperatureIBC  = static_cast<IceTemperatureIBC*>
	  (new PythonInterface::PythonIceTemperatureIBC(module, funcName));
      }
#endif
    else 
      {
	MayDay::Error("bad temperature type");
      }	
	
    amrObject.setTemperatureBC(temperatureIBC);
  
    {
      // ---------------------------------------------
      // set surface heat boundary data 
      // ---------------------------------------------
      
      SurfaceFlux* surf_heat_boundary_data_ptr = SurfaceFlux::parseSurfaceFlux("surfaceHeatBoundaryData");
      ParmParse pps("surfaceHeatBoundaryData");
      bool diri = true; //Dirichlett boundary data by default
      pps.query("Dirichlett",diri);
      if (surf_heat_boundary_data_ptr == NULL)
	{
	  if (!diri)
	    {
	      surf_heat_boundary_data_ptr = new zeroFlux();
	    }
	}
      
      amrObject.setSurfaceHeatBoundaryData(surf_heat_boundary_data_ptr, diri);
      if (surf_heat_boundary_data_ptr != NULL)
	{
	  delete surf_heat_boundary_data_ptr;
	  surf_heat_boundary_data_ptr=NULL;
	}
    
      // ---------------------------------------------
      // set basal (lower surface) heat boundary data. 
      // ---------------------------------------------
      
      SurfaceFlux* basal_heat_boundary_data_ptr = SurfaceFlux::parseSurfaceFlux("basalHeatBoundaryData");
      // if (basal_heat_boundary_data_ptr == NULL)
      // 	{
      // 	  basal_heat_boundary_data_ptr = new zeroFlux();
      // 	}
      
      amrObject.setBasalHeatBoundaryData(basal_heat_boundary_data_ptr);
      if (basal_heat_boundary_data_ptr != NULL)
	{
	  delete basal_heat_boundary_data_ptr;
	  basal_heat_boundary_data_ptr=NULL;
	}
      
    }
      
    amrObject.setDomainSize(domainSize);
    // set up initial grids, initialize data, etc.
    amrObject.initialize();
  

    int maxStep;
    Real maxTime;
    //Real startTime;
    pp2.get("maxTime", maxTime);
    pp2.get("maxStep", maxStep);
    
    amrObject.run(maxTime, maxStep);
    
    // clean up
    if (constRelPtr != NULL)
      {
        delete constRelPtr;
        constRelPtr = NULL;
      }

    if (surf_flux_ptr != NULL)
      {
        delete surf_flux_ptr;
        surf_flux_ptr = NULL;
      }

    if (basal_flux_ptr != NULL)
      {
        delete basal_flux_ptr;
        basal_flux_ptr = NULL;
      }

    if (basalFrictionPtr != NULL)
      {
	delete basalFrictionPtr;
	basalFrictionPtr = NULL;
      }

    if (basalFrictionRelationPtr != NULL)
      {
	delete basalFrictionRelationPtr;
	basalFrictionRelationPtr = NULL;
      }

    if (thicknessIBC != NULL)
      {
        delete thicknessIBC;
        thicknessIBC=NULL;
      }

    if (temperatureIBC != NULL)
      {  
	delete temperatureIBC;
        temperatureIBC=NULL;
      }

  }  // end nested scope
  

  CH_TIMER_REPORT();

#ifdef HAVE_PYTHON
  Py_Finalize();
#endif

#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional

  return ierr;
}


