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

//     //geometry set up
//     Vector<int> ncells(3); 
//     ppMain.getarr("num_cells",ncells,0,ncells.size());
//     bool periodic[SpaceDim] = {D_DECL(false,false,false)};
//     {
//       Vector<int> is_periodic(SpaceDim);
//       ppMain.getarr("is_periodic", is_periodic, 0, SpaceDim);
//       for (int dir=0; dir<SpaceDim; dir++) 
// 	{
// 	  periodic[dir] = (is_periodic[dir] == 1);
// 	}
//     }


//     IntVect lo(IntVect::Zero);
//     IntVect hi(D_DECL(ncells[0]-1, ncells[1]-1, ncells[2]-1));
//     ProblemDomain pd(lo,hi,periodic);
// #if BISICLES_Z == BISICLES_LAYERED
//     Vector<Real> faceSigma(ncells[2]+1);
//     Real dz = 1.0 / ncells[2];
//     faceSigma[0] = 0.0;
//     for (int j = 1; j < faceSigma.size()-1; ++j)
// 	faceSigma[j] = faceSigma[j-1] + dz;
//     faceSigma[faceSigma.size()-1] = 1.0;
// #endif

   

//     RealVect domainSize;
//     Vector<Real> domSize(SpaceDim);
//     ppMain.getarr("domain_size", domSize, 0, SpaceDim);
//     domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));
   
//     RealVect dx (domainSize);
//     for (int dir = 0; dir < SpaceDim; ++dir)
//       dx[dir] /= Real(ncells[dir]);
//     RealVect dataDx = dx;

//     int blockFactor = 8;
//     int maxBoxSize = 1000000;

//     Vector<Box> baseBoxes;
//     domainSplit(pd, baseBoxes,maxBoxSize,blockFactor);
//     Vector<int> procAssign(baseBoxes.size());
//     LoadBalance(procAssign,baseBoxes);
//     DisjointBoxLayout grids(baseBoxes, procAssign, pd);

    

    RealVect dataDx;
    RefCountedPtr<LevelData<FArrayBox> > levelThck;
    RefCountedPtr<LevelData<FArrayBox> > levelTopg;
    RefCountedPtr<LevelData<FArrayBox> > levelXVel;
    RefCountedPtr<LevelData<FArrayBox> > levelYVel;
    RefCountedPtr<LevelData<FArrayBox> > levelVel;
    RefCountedPtr<LevelData<FArrayBox> > levelVelCoef;
    RefCountedPtr<LevelData<FArrayBox> > levelC;
    RefCountedPtr<LevelData<FArrayBox> > levelDivUH;
    RefCountedPtr<LevelData<FArrayBox> > levelDivUHCoef;
    //read input data
    std::string inputMethod = "";
    ppMain.get("inputMethod", inputMethod );

    if (inputMethod == "LevelData")
      {
	ParmParse ildPP("inputLevelData");
	std::string infile;
	ildPP.get("inputFile",infile);

	Vector<std::string> names;
	 Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;

	std::string thicknessName = "thk";
	ildPP.query("thicknessName",thicknessName);
	names.push_back(thicknessName);
	levelThck =  RefCountedPtr<LevelData<FArrayBox> >
	 (new LevelData<FArrayBox>); 
	vectData.push_back(levelThck);

	std::string topographyName = "topg";
	ildPP.query("topographyName",topographyName);
	names.push_back(topographyName);
	levelTopg = RefCountedPtr<LevelData<FArrayBox> >
	 (new LevelData<FArrayBox>);
	vectData.push_back(levelTopg);

	std::string frictionName = "beta";
	ildPP.query("frictionName",frictionName);
	names.push_back(frictionName);
	levelC = RefCountedPtr<LevelData<FArrayBox> >
	  (new LevelData<FArrayBox>);
	vectData.push_back(levelC);

	std::string xvelName = "xvel";
	ildPP.query("xvelName",xvelName);
	names.push_back(xvelName);
	levelXVel = RefCountedPtr<LevelData<FArrayBox> >
	  (new LevelData<FArrayBox>);
	vectData.push_back(levelXVel);

	std::string yvelName = "yvel";
	ildPP.query("yvelName",yvelName);
	names.push_back(yvelName);
	levelYVel = RefCountedPtr<LevelData<FArrayBox> >
	  (new LevelData<FArrayBox>);
	vectData.push_back(levelYVel);
	
	std::string velcoefName = "velcoef";
	ildPP.query("velcoefName",velcoefName);
	names.push_back(velcoefName);
	levelVelCoef = RefCountedPtr<LevelData<FArrayBox> >
	   (new LevelData<FArrayBox>);
	vectData.push_back(levelVelCoef);
	

	std::string divuhName = "divuh";
	ildPP.query("divuhName",divuhName);
	names.push_back(divuhName);
	levelDivUH = RefCountedPtr<LevelData<FArrayBox> >
	   (new LevelData<FArrayBox>);
	vectData.push_back(levelDivUH);

	std::string divuhcoefName = "divuhcoef";
	ildPP.query("divuhcoefName",divuhcoefName);
	names.push_back(divuhcoefName);
	levelDivUHCoef = RefCountedPtr<LevelData<FArrayBox> >
	  (new LevelData<FArrayBox>);
	vectData.push_back(levelDivUHCoef);

	Real dx = 0.0;
	readLevelData(vectData,dx,infile,names,1);
	
	dataDx[0]=dx;dataDx[1]=dx;

	pout() << "read in data from hdf5 file " << infile  << std::endl;

	levelVel = RefCountedPtr<LevelData<FArrayBox> >
	  (new LevelData<FArrayBox>(levelXVel->disjointBoxLayout(),
				    SpaceDim,IntVect::Zero));
	levelXVel->copyTo(Interval(0,0),*levelVel,Interval(0,0));
	levelYVel->copyTo(Interval(0,0),*levelVel,Interval(1,1));

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
   
    LevelDataIBC ibc(levelThck,levelTopg,dataDx);


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
	std::string infile;
	ildPP.get("temperatureFile",infile);
	std::string temperatureName = "temp000000";
	ildPP.query("temperatureName",temperatureName);
	RefCountedPtr<LevelData<FArrayBox> > levelTemp
	  (new LevelData<FArrayBox>());
	Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
	vectData.push_back(levelTemp);
	Vector<std::string> names(1);
	names[0] = temperatureName;
	Real dx;
	ParmParse ppAmr ("geometry");
	Vector<int> ancells(3); 
	ppAmr.getarr("num_cells", ancells, 0, ancells.size());
	readLevelData(vectData,dx,infile,names,ancells[2]);
	RealVect levelDx = RealVect::Unit * dx;
	LevelDataTemperatureIBC* ptr = new LevelDataTemperatureIBC(levelTemp,levelDx);
	temperatureIBC  = static_cast<IceTemperatureIBC*>(ptr);
      }
    else 
      {
	MayDay::Error("bad temperature type");
      }	


    AMRIceControl amrIceControl;
    amrIceControl.define(&ibc, temperatureIBC,  rateFactorPtr, constRelPtr , 
			 basalFrictionRelationPtr,dataDx,levelC, levelVel,
			 levelVelCoef,levelDivUH,levelDivUHCoef);

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
