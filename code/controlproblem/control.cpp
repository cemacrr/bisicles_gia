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
// control.cpp
// contruct and solve a control problem
// by a method related to Vieli and Payne, Ann Glaciol, 2003. vol 26
// assumes a 2D SSA model or L1L2 model
// provide main
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
#include "BasalFrictionRelation.H"
#include "FillFromReference.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "LevelDataIBC.H"
#include "MultiLevelDataIBC.H"
#include "LevelDataTemperatureIBC.H"
#include "AMRIceControl.H"
#include "ReadLevelData.H"
//#include "FlatHDF5.H"

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
    

    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);
    ParmParse ppMain("main");

    Vector<RealVect> dataDx;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > xVelObs;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > yVelObs;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > velObs;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > velCoef;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > originC;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > originMuCoef;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > divUHObs;
    Vector<RefCountedPtr<LevelData<FArrayBox> > > divUHCoef;
    
    //read input data
    std::string inputMethod = "";
    ppMain.get("inputMethod", inputMethod );

    if (inputMethod == "LevelData")
      {
	ParmParse ildPP("inputLevelData");
	std::string infile;
	ildPP.get("inputFile",infile);

	Vector<std::string> names;
	

	//read names of the input fields
	std::string frictionName = "beta";
	ildPP.query("frictionName",frictionName);
	names.push_back(frictionName);

	std::string muCoefName = "";
	ildPP.query("muCoefName",muCoefName);
	if (muCoefName != "")
	  {
	    //loading muCoef is optional
	    names.push_back(muCoefName);
	  }

	std::string xvelName = "xvel";
	ildPP.query("xvelName",xvelName);
	names.push_back(xvelName);
	std::string yvelName = "yvel";
	ildPP.query("yvelName",yvelName);
	names.push_back(yvelName);
	std::string velcoefName = "velcoef";
	ildPP.query("velcoefName",velcoefName);
	names.push_back(velcoefName);
	std::string divuhName = "divuh";
	ildPP.query("divuhName",divuhName);
	names.push_back(divuhName);
	std::string divuhcoefName = "divuhcoef";
	ildPP.query("divuhcoefName",divuhcoefName);
	names.push_back(divuhcoefName);

	Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectData;
	Real dx;
	Vector<int> refRatio;
	readMultiLevelData(vectData,dx,refRatio,infile,names,1);

	

	for (int lev =0 ; lev < vectData[0].size(); lev++)
	  {
	    int j = 0;
	    originC.push_back( vectData[j++][lev]);
	    if (muCoefName != "")
	      {
		originMuCoef.push_back( vectData[j++][lev]);
	      }
	    else
	      {
		//if there was no initial muCoef, then muCoef = 1.0 is appropriate
		//for backward compatibility and in general.
		originMuCoef.push_back
		  (RefCountedPtr<LevelData<FArrayBox> >
		   (new  LevelData<FArrayBox>
		    (originC[lev]->disjointBoxLayout(),originC[lev]->nComp(),originC[lev]->ghostVect())));
		for (DataIterator dit(originMuCoef[lev]->disjointBoxLayout());dit.ok();++dit)
		  {
		    (*originMuCoef[lev])[dit].setVal(1.0);
		  }
	      }
	    xVelObs.push_back( vectData[j++][lev]);
	    yVelObs.push_back( vectData[j++][lev]);
	    velCoef.push_back( vectData[j++][lev]);
	    divUHObs.push_back( vectData[j++][lev]);
	    divUHCoef.push_back( vectData[j++][lev]);
	    velObs.push_back 
	      (RefCountedPtr<LevelData<FArrayBox> > 
	       (new LevelData<FArrayBox>
		(xVelObs[lev]->disjointBoxLayout(),SpaceDim,xVelObs[lev]->ghostVect())));
	    xVelObs[lev]->copyTo(Interval(0,0),*velObs[lev],Interval(0,0));
	    yVelObs[lev]->copyTo(Interval(0,0),*velObs[lev],Interval(1,1));
	  }
	dataDx.resize(xVelObs.size());
	dataDx[0][0]=dx;dataDx[0][1]=dx;
	for (int lev = 1; lev < refRatio.size(); lev++)
	  {
	    dataDx[lev] = dataDx[lev-1] / Real(refRatio[lev-1]);
	  }

	pout() << "read in data from hdf5 file " << infile  << std::endl;

	

      }
    else
      {
	MayDay::Error("unknown inputMethod");
      }

    // ---------------------------------------------
    // set constitutive relation & rate factor
    // ---------------------------------------------
    std::string constRelType;
    ppMain.get("constitutiveRelation", constRelType);
    ConstitutiveRelation* constRelPtr = NULL;
    GlensFlowRelation* gfrPtr = NULL;
    if (constRelType == "constMu")
      {
        constMuRelation* newPtr = new constMuRelation;
        ParmParse crPP("constMu");
        Real muVal;
        crPP.get("mu", muVal);
        newPtr->setConstVal(muVal);
        constRelPtr = static_cast<ConstitutiveRelation*>(newPtr);
      }
    else if (constRelType == "GlensLaw")
      {
        constRelPtr = new GlensFlowRelation;
	gfrPtr = dynamic_cast<GlensFlowRelation*>(constRelPtr);
      }
    else if (constRelType == "L1L2")
      {
	L1L2ConstitutiveRelation* l1l2ptr = new L1L2ConstitutiveRelation;
	l1l2ptr->parseParameters();
	gfrPtr = l1l2ptr->getGlensFlowRelationPtr();
        constRelPtr = l1l2ptr;
      }
    else 
      {
        MayDay::Error("bad Constitutive relation type");
      }
    Real epsSqr0 = 1.0e-9;
    std::string rateFactorType = "constRate";
    RateFactor*  rateFactorPtr = NULL;
    ppMain.query("rateFactor", rateFactorType);
    if (rateFactorType == "constRate")
      {
	ParmParse crPP("constRate");
	Real A = 9.2e-18;
	crPP.query("A", A);
	ConstantRateFactor rateFactor(A);
	rateFactorPtr = rateFactor.getNewRateFactor();
	crPP.query("epsSqr0", epsSqr0);
	
	if (gfrPtr) 
	  {
	    gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	  }
      }
    else if (rateFactorType == "arrheniusRate")
      {
	ArrheniusRateFactor rateFactor;
	ParmParse arPP("ArrheniusRate");
	arPP.query("epsSqr0", epsSqr0);
	rateFactorPtr = rateFactor.getNewRateFactor();
	if (gfrPtr) 
	  {
	    gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	  }
      }
    else if (rateFactorType == "patersonRate")
      {
	PatersonRateFactor rateFactor;
	ParmParse arPP("PatersonRate");
	arPP.query("epsSqr0", epsSqr0);
	rateFactorPtr = rateFactor.getNewRateFactor();
	if (gfrPtr) 
	  {
	    gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	  }
      }
  
    // ---------------------------------------------------
    // set initial basal friction coefficient and relation
    // ---------------------------------------------------

    BasalFrictionRelation* basalFrictionRelationPtr = NULL;
    std::string basalFrictionRelType = "powerLaw";
    ppMain.query("basalFrictionRelation", basalFrictionRelType);
    
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
   
    

    IceThicknessIBC* thicknessIBC = NULL;
    ParmParse geomPP("geometry");
    std::string geomType("LevelData");
    geomPP.query("type",geomType);
    if (geomType == "LevelData" || geomType == "MultiLevelData")
      {
	// MultiLevelThicknessIBC can replace LevelDataThicknessIBC
	ParmParse ildPP("inputLevelData");
	
	std::string infile="";
	ildPP.query("inputFile",infile);
	ildPP.query("geometryFile",infile);
	if (infile == "")
	  {
	    MayDay::Error("inputLevelData.inputFile or inputLevelData.geometryFile must be specified");
	  }
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
	
	thicknessIBC  = static_cast<IceThicknessIBC*>(ptr);



      }
    else 
      {
	MayDay::Error("bad geometry type");
      }

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
	CH_assert( (ptr = LevelDataTemperatureIBC::parse(ildPP)) != NULL);
	temperatureIBC  = static_cast<IceTemperatureIBC*>(ptr);

      }
    else 
      {
	MayDay::Error("bad temperature type");
      }	


    AMRIceControl amrIceControl;
    amrIceControl.define(thicknessIBC, temperatureIBC,  rateFactorPtr, constRelPtr , 
			 basalFrictionRelationPtr,dataDx,originC,originMuCoef,velObs,
			 velCoef,divUHObs,divUHCoef);

    std::string problem = "control";
    ppMain.query("problem", problem);
    if (problem == "forward")
      {
	amrIceControl.solveForward();
      }
    else if (problem == "control")
      {
	amrIceControl.solveControl();
      }
    else
      {
	MayDay::Error("undefined problem in inputs");
      }

    if (constRelPtr != NULL)
      {
        delete constRelPtr;
        constRelPtr = NULL;
      }

    if (rateFactorPtr != NULL)
      {
        delete rateFactorPtr;
        rateFactorPtr = NULL;
      }

    if (thicknessIBC != NULL)
      {
	delete thicknessIBC;
	thicknessIBC = NULL;
      }
    
    if (temperatureIBC != NULL)
      {
	delete temperatureIBC;
	temperatureIBC = NULL;
      }

    if (basalFrictionRelationPtr != NULL)
      {
	delete basalFrictionRelationPtr;
	basalFrictionRelationPtr = NULL;
      }
 


  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
