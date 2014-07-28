
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
#ifdef CH_USE_PETSC
#include "petsc.h"
#endif 

int main(int argc, char* argv[]) {

#ifdef CH_USE_PETSC
  int ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
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
    ParmParse ppMain("main");

    std::string poutBaseName = "pout";
    ppMain.query("poutBaseName",poutBaseName);
    setPoutBaseName(poutBaseName);

    AMRIceControl amrIceControl;
    IceThicknessIBC* thicknessIBC = NULL;
    IceTemperatureIBC* temperatureIBC = NULL;
    BasalFrictionRelation* basalFrictionRelationPtr = NULL;
    ConstitutiveRelation* constRelPtr = NULL;
    RateFactor*  rateFactorPtr = NULL;

    {
      Vector<RealVect> dataDx;
      Vector<RefCountedPtr<LevelData<FArrayBox> > > xVelObs;
      Vector<RefCountedPtr<LevelData<FArrayBox> > > yVelObs;
      Vector<RefCountedPtr<LevelData<FArrayBox> > > velObs;
      Vector<RefCountedPtr<LevelData<FArrayBox> > > velCoef;
      Vector<RefCountedPtr<LevelData<FArrayBox> > > thkCoef;
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
	  
	  std::string yvelName = "";
	  ildPP.query("yvelName",yvelName);
	  if (yvelName != "")
	    {
	      //loading yvel is optional
	      names.push_back(yvelName);
	    }
	  
	  std::string velcoefName = "velcoef";
	  ildPP.query("velcoefName",velcoefName);
	  names.push_back(velcoefName);
	  
	  std::string thkcoefName = "";
	  ildPP.query("thkcoefName",thkcoefName);
	  if (thkcoefName != "")
	    {
	      //loading thkcoef is optional
	      names.push_back(thkcoefName);
	    }

	  std::string divuhName = "";
	  ildPP.query("divuhName",divuhName);
	  if (divuhName != "")
	    {
	      //loading divuh is optional
	      names.push_back(divuhName);
	      //but divuhcoef is required as well
	      std::string divuhcoefName = "divuhcoef";
	      ildPP.query("divuhcoefName",divuhcoefName);
	      names.push_back(divuhcoefName);
	    }
	  
	  
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
	      if (yvelName != "")
		{
		  yVelObs.push_back( vectData[j++][lev]);
		}
	      else
		{
		  //yvelobs = 0 by default (in which case mod(u) = u_x)
		  yVelObs.push_back
		    (RefCountedPtr<LevelData<FArrayBox> >
		     (new  LevelData<FArrayBox>
		      (originC[lev]->disjointBoxLayout(),originC[lev]->nComp(),originC[lev]->ghostVect())));
		  for (DataIterator dit(yVelObs[lev]->disjointBoxLayout());dit.ok();++dit)
		    {
		      (*yVelObs[lev])[dit].setVal(0.0);
		    }
		}
	      
	      velCoef.push_back( vectData[j++][lev]);
	      
	      if (thkcoefName != "")
		{
		  thkCoef.push_back( vectData[j++][lev]);
		}
	      else
		{
		  //thkCoef = 1 by default
		  thkCoef.push_back
		    (RefCountedPtr<LevelData<FArrayBox> >
		     (new  LevelData<FArrayBox>
		      (originC[lev]->disjointBoxLayout(),1,originC[lev]->ghostVect())));
		  for (DataIterator dit(thkCoef[lev]->disjointBoxLayout());dit.ok();++dit)
		    {
		      (*thkCoef[lev])[dit].setVal(1.0);
		    }
		}
	      
	      if (divuhName != "")
		{
		  divUHObs.push_back( vectData[j++][lev]);
		  divUHCoef.push_back( vectData[j++][lev]);
		}
	      else
		{
		  //divuh,divuhcoef = 0,1 by default (in which case mod(u) = u_x)
		  divUHObs.push_back
		    (RefCountedPtr<LevelData<FArrayBox> >
		     (new  LevelData<FArrayBox>
		      (originC[lev]->disjointBoxLayout(),originC[lev]->nComp(),originC[lev]->ghostVect())));

		  divUHCoef.push_back
		    (RefCountedPtr<LevelData<FArrayBox> >
		     (new  LevelData<FArrayBox>
		      (originC[lev]->disjointBoxLayout(),originC[lev]->nComp(),originC[lev]->ghostVect())));
		  
		  for (DataIterator dit(divUHObs[lev]->disjointBoxLayout());dit.ok();++dit)
		    {
		      (*divUHObs[lev])[dit].setVal(0.0);
		      (*divUHCoef[lev])[dit].setVal(1.0);
		    }
		}
	      
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
      ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");
      
      if (constRelPtr == NULL)
	{
	  MayDay::Error("undefined constitutiveRelation in inputs");
	}

      Real epsSqr0 = 1.0e-9;
      std::string rateFactorType = "constRate";
     
      ppMain.query("rateFactor", rateFactorType);
      if (rateFactorType == "constRate")
	{
	  ParmParse crPP("constRate");
	  Real A = 9.2e-18;
	  crPP.query("A", A);
	  ConstantRateFactor rateFactor(A);
	  rateFactorPtr = rateFactor.getNewRateFactor();
	  crPP.query("epsSqr0", epsSqr0);
	
	 
	}
      else if (rateFactorType == "arrheniusRate")
	{
	  ArrheniusRateFactor rateFactor;
	  ParmParse arPP("ArrheniusRate");
	  arPP.query("epsSqr0", epsSqr0);
	  rateFactorPtr = rateFactor.getNewRateFactor();
	 
	}
      else if (rateFactorType == "patersonRate")
	{
	  PatersonRateFactor rateFactor;
	  ParmParse arPP("PatersonRate");
	  arPP.query("epsSqr0", epsSqr0);
	  rateFactorPtr = rateFactor.getNewRateFactor();
	 
	}
    else if (rateFactorType == "zwingerRate")
	{
	  ZwingerRateFactor rateFactor;
	  ParmParse arPP("ZwingerRate");
	  arPP.query("epsSqr0", epsSqr0);
	  rateFactorPtr = rateFactor.getNewRateFactor();
	 
	}
      // ---------------------------------------------------
      // set initial basal friction coefficient and relation
      // ---------------------------------------------------

     
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

	  if (geomType == "LevelData")
	    {
	      //TODO : MultiLevelThicknessIBC can replace LevelDataThicknessIBC
	      //when it is finished 
	      CH_assert(vectData[0].size() == 1);
	      CH_assert(vectData[1].size() == 1);
	      LevelDataIBC* ptr = new LevelDataIBC
		(vectData[0][0],vectData[1][0],crseDx);
	      thicknessIBC  = static_cast<IceThicknessIBC*>(ptr);
	    }
	  else
	    {
	      //MultiLevelThicknessIBC can replace LevelDataThicknessIBC
	      MultiLevelDataIBC* ptr = new MultiLevelDataIBC
		(vectData[0],vectData[1],crseDx,refRatio);
	      thicknessIBC  = static_cast<IceThicknessIBC*>(ptr);
	    }
	}
      else 
	{
	  MayDay::Error("bad geometry type");
	}

     
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


    
      amrIceControl.define(thicknessIBC, temperatureIBC,  rateFactorPtr, constRelPtr , 
			   basalFrictionRelationPtr,dataDx,originC,originMuCoef,velObs,
			   velCoef,thkCoef, divUHObs,divUHCoef);

    }

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

#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif// mpi conditional
#endif // petsc conditional
  return 0;
}
