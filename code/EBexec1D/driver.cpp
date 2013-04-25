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
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "AmrIce.H"
#include "EBAmrIce.H"
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
#include "areaWeightedFriction.H"
#include "BasalFrictionRelation.H"
#include "MuCoefficient.H"
#include "twistyStreamFriction.H"
#include "GaussianBumpFriction.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "HumpIBC.H"
#include "LevelDataIBC.H"
#include "IceTemperatureIBC.H"
#include "LevelDataTemperatureIBC.H"
#include "LevelDataBasalFriction.H"
#include "PiecewiseLinearFlux.H"
#include "SurfaceFlux.H"
#include "IceConstants.H"
#include "CH_Attach.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "ReadLevelData.H"
#include "PetscSolver.H"

extern "C" {
    // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  
  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}

/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
    guassianBump is used for the MISMIP3D perturbations tests.
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
			 gaussianBump,
                         areaWeighted,
                         NUM_BETA_TYPES};

RateFactor           * rateFactor          (Real& a_epsSqr0);
ConstitutiveRelation * constitutiveRelation(RateFactor* a_rateFactor, const Real& a_epsSqr0);

SurfaceFlux          * surfaceFlux          ();
SurfaceFlux          * basalSurfaceFlux     ();
MuCoefficient        * muCoef               ();

RealVect               readDomainSize       ();
BasalFrictionRelation* basalFrictionRelation();

BasalFriction        * basalFriction        (const RealVect& a_domainSize);
IceThicknessIBC      * thicknessIBC         (const RealVect& a_domainSize);

IceTemperatureIBC    * temperatureIBC       ();




int main(int argc, char* argv[]) {
  
  int ierr = 0;
 
  
  // petsc ifdef
#ifdef CH_USE_PETSC
  ierr = PetscInitialize(&argc, &argv,"./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  registerDebugger();
#endif 
#endif 

  { 
    
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
      { 
        std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); 
      }
    
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);
        
    // computes evolution of ice thickness, velocity, and grounding line
    EBAmrIce amrObject;
            
    // todo: need comment for epsSqr0
    Real epsSqr0;
    
    // helper classes constructed from input file parameters
    RateFactor           *rateFactorPtr = rateFactor          (              epsSqr0);
    ConstitutiveRelation *constRelPtr   = constitutiveRelation(rateFactorPtr,epsSqr0);
    
    SurfaceFlux          *surfFluxPtr               = surfaceFlux          ();
    SurfaceFlux          *basalFluxPtr              = basalSurfaceFlux     ();
    MuCoefficient        *muCoefficientPtr          = muCoef               ();
    RealVect              domainSize                = readDomainSize       ();
    BasalFrictionRelation*basalFrictionRelationPtr  = basalFrictionRelation();
    BasalFriction        *basalFrictionPtr          = basalFriction        (domainSize);
    IceThicknessIBC      *thicknessIBCPtr           = thicknessIBC         (domainSize);
    IceTemperatureIBC    *temperatureIBCPtr         = temperatureIBC       ();
    
    ParmParse pp2("main");
    bool useAreaWeightedFriction;
    pp2.get("useAreaWeightedFriction",useAreaWeightedFriction);
        
    // define an EBAmrIce using areaWeighted friction
    areaWeightedFriction* areaWeightedFrictionPtr = new areaWeightedFriction(basalFrictionPtr);
      
    if (useAreaWeightedFriction)
      {
        // set member data pointers (the set functions call "new")
        amrObject.define(constRelPtr,
                         rateFactorPtr,
                         surfFluxPtr,
                         basalFluxPtr,
                         muCoefficientPtr,
                         basalFrictionRelationPtr,
                         areaWeightedFrictionPtr,
                         thicknessIBCPtr,
                         temperatureIBCPtr,
                         domainSize);
      }
    else
      {
         // set member data pointers (the set functions call "new")
        amrObject.define(constRelPtr,
                         rateFactorPtr,
                         surfFluxPtr,
                         basalFluxPtr,
                         muCoefficientPtr,
                         basalFrictionRelationPtr,
                         basalFrictionPtr,
                         thicknessIBCPtr,
                         temperatureIBCPtr,
                         domainSize);
      }
    
    // num steps
    int maxStep;
    pp2.get("maxStep", maxStep);
    
    // time
    Real maxTime;
    pp2.get("maxTime", maxTime);
        
    
    // run the simulation
    amrObject.run(maxTime, maxStep);
    
    // clean up
    if (constRelPtr != NULL)
      {
        delete constRelPtr;
        constRelPtr = NULL;
      }

    if (surfFluxPtr != NULL)
      {
        delete surfFluxPtr;
        surfFluxPtr = NULL;
      }

    if (basalFluxPtr != NULL)
      {
        delete basalFluxPtr;
        basalFluxPtr = NULL;
      }

    if (basalFrictionPtr != NULL)
      {
	delete basalFrictionPtr;
	basalFrictionPtr = NULL;
      }
    
    if (areaWeightedFrictionPtr != NULL)
      {
        delete areaWeightedFrictionPtr;
        areaWeightedFrictionPtr = NULL;
      }
      
    if (basalFrictionRelationPtr != NULL)
      {
	delete basalFrictionRelationPtr;
	basalFrictionRelationPtr = NULL;
      }

    if (thicknessIBCPtr != NULL)
      {
        delete thicknessIBCPtr;
        thicknessIBCPtr=NULL;
      }

    if (temperatureIBCPtr != NULL)
      {  
	delete temperatureIBCPtr;
        temperatureIBCPtr=NULL;
      }

  } 

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

ConstitutiveRelation* constitutiveRelation(RateFactor* a_rateFactor, const Real& a_epsSqr0)
{
  
  ConstitutiveRelation* constRelPtr = NULL;

  // read constitutive relation type
  ParmParse pp2("main");
  std::string constRelType;
  pp2.get("constitutiveRelation", constRelType);
  
  // set constitutive relation 
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
      GlensFlowRelation* gfrPtr = new GlensFlowRelation;
      gfrPtr->setParameters(3.0 , a_rateFactor, a_epsSqr0);
      constRelPtr = gfrPtr;
    }
  else if (constRelType == "L1L2")
    {
      L1L2ConstitutiveRelation* l1l2Ptr = new L1L2ConstitutiveRelation;
      l1l2Ptr->parseParameters();
      l1l2Ptr->getGlensFlowRelationPtr()->setParameters(3.0 , a_rateFactor, a_epsSqr0);
      constRelPtr = l1l2Ptr;
    }
  else 
    {
      MayDay::Error("bad Constitutive relation type");
    }

  return constRelPtr;
}

RateFactor* rateFactor(Real& a_epsSqr0) 
{ 
  RateFactor* rateFactor;

  // set rate factor. Also, Glen's flow relation parameters, if appropriate
  ParmParse pp2("main");

  std::string rateFactorType = "constRate";
  pp2.query("rateFactor", rateFactorType);
  if (rateFactorType == "constRate")
    {
      ParmParse crPP("constRate");

      Real A = 9.2e-18;
      crPP.query("A", A);
      rateFactor = new ConstantRateFactor(A);
  
      crPP.query("epsSqr0", a_epsSqr0);
    }
  else if (rateFactorType == "arrheniusRate")
    {
      ParmParse arPP("ArrheniusRate");
      
      rateFactor = new ArrheniusRateFactor;
      
      arPP.query("epsSqr0", a_epsSqr0);
    }

  return rateFactor;
}

SurfaceFlux* surfaceFlux()
{
  SurfaceFlux* surfFluxPtr(SurfaceFlux::parseSurfaceFlux("surfaceFlux"));
  if (surfFluxPtr == NULL)
    {
      MayDay::Error("invalid surface flux type");
    }

  return surfFluxPtr;
}

SurfaceFlux* basalSurfaceFlux()
{
  SurfaceFlux* basalFluxPtr = SurfaceFlux::parseSurfaceFlux("basalFlux");
     
  if (basalFluxPtr == NULL)
    {
      MayDay::Error("invalid basal flux type");
    }

  return basalFluxPtr;
}

MuCoefficient* muCoef()
{
  MuCoefficient* muCoefficientPtr;
  // set mu coefficient
  ParmParse muPP("muCoefficient");
  std::string muCoefType = "unit";
  muPP.query("type",muCoefType );
  if (muCoefType == "unit")
    {
      muCoefficientPtr = static_cast<MuCoefficient*>(new UnitMuCoefficient());
    }
  else if (muCoefType == "LevelData")
    {
      //read a one level muCoef from an AMR Hierarchy, and  store it in a LevelDataMuCoeffcient
      ParmParse ildPP("inputLevelData");
      std::string infile;
      ildPP.get("muCoefFile",infile);
      std::string frictionName = "muCoef";
      ildPP.query("muCoefName",frictionName);
      RefCountedPtr<LevelData<FArrayBox> > levelMuCoef (new LevelData<FArrayBox>());
      Vector<RefCountedPtr<LevelData<FArrayBox> > > vectMuCoef;
      vectMuCoef.push_back(levelMuCoef);
      Vector<std::string> names(1);
      names[0] = frictionName;
      Real dx;
      readLevelData(vectMuCoef,dx,infile,names,1);
      RealVect levelDx = RealVect::Unit * dx;
      muCoefficientPtr = static_cast<MuCoefficient*>(new LevelDataMuCoefficient(levelMuCoef,levelDx));
    }
  else
    {
      MayDay::Error("undefined MuCoefficient in inputs");
    }

 return muCoefficientPtr;
}

BasalFrictionRelation* basalFrictionRelation()
{
  BasalFrictionRelation* basalFrictionRelationPtr;
  ParmParse pp2("main");

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

  return basalFrictionRelationPtr;
}

BasalFriction* basalFriction(const RealVect& a_domainSize)
{
  ParmParse geomPP("geometry");
  
  BasalFriction* basalFrictionPtr = NULL;
  
  std::string beta_type;
  geomPP.get("beta_type", beta_type);
  
  // read in type of beta^2 distribution
  if (beta_type == "constantBeta")
    {
      Real betaVal;
      geomPP.get("betaValue", betaVal);
      basalFrictionPtr = static_cast<BasalFriction*>(new constantFriction(betaVal));
    }
  else if (beta_type == "sinusoidalBeta")
    {
      Real betaVal, eps;
      RealVect omega(RealVect::Unit);
      Vector<Real> omegaVect(SpaceDim);
      geomPP.get("betaValue", betaVal);
      if (geomPP.contains("omega"))
        {
          geomPP.getarr("omega", omegaVect, 0, SpaceDim);
          omega = RealVect(D_DECL(omegaVect[0], omegaVect[1], omegaVect[2]));
        }
      geomPP.get("betaEps", eps);
      basalFrictionPtr = static_cast<BasalFriction*>(new sinusoidalFriction(betaVal, 
                                                                            omega, 
                                                                            eps,
                                                                            a_domainSize));
    }
  
  // special case of sinusoidalBeta
  else if (beta_type == "sinusoidalBetay")
    {
      Real betaVal, eps, omegaVal;
      RealVect omega(RealVect::Zero);
      omega[1] = 1;
      
      geomPP.get("betaValue", betaVal);
      if (geomPP.contains("omega"))
        {
          geomPP.get("omega", omegaVal);
          omega[1] = omegaVal;
        }
      geomPP.get("betaEps", eps);
      basalFrictionPtr = static_cast<BasalFriction*>(new sinusoidalFriction(betaVal, 
                                                                            omega, 
                                                                            eps,
                                                                            a_domainSize));
      
    }
    else if (beta_type == "twistyStreamx")
      {
        Real betaVal, eps, magOffset;
        magOffset = 0.25;
        RealVect omega(RealVect::Unit);
        Vector<Real> omegaVect(SpaceDim);
        geomPP.get("betaValue", betaVal);
        if (geomPP.contains("omega"))
          {
            geomPP.getarr("omega", omegaVect, 0, SpaceDim);
            omega = RealVect(D_DECL(omegaVect[0], omegaVect[1], omegaVect[2]));
          }
        geomPP.query("magOffset", magOffset);
        geomPP.get("betaEps", eps);
        basalFrictionPtr = static_cast<BasalFriction*>(new twistyStreamFriction(betaVal, 
                                                                                omega, 
                                                                                magOffset, 
                                                                                eps,
                                                                                a_domainSize));          
      }
    else if (beta_type == "gaussianBump")
      {
	int nt;
	geomPP.get("gaussianBump_nt", nt);
	Vector<Real> t(nt-1);
	Vector<Real> C0(nt),a(nt);
	Vector<RealVect> b(nt),c(nt);
        
	geomPP.getarr("gaussianBump_t", t, 0, nt-1);
	geomPP.getarr("gaussianBump_C", C0, 0, nt);
	geomPP.getarr("gaussianBump_a", a, 0, nt);

#if CH_SPACEDIM == 1
	Vector<Real> xb(nt),xc(nt);
	geomPP.getarr("gaussianBump_xb", xb, 0, nt);
	geomPP.getarr("gaussianBump_xc", xc, 0, nt);
	for (int i = 0; i < nt; ++i)
	  {
	    b[i][0] = xb[i];
	    c[i][0] = xc[i];
	  }     
#elif CH_SPACEDIM == 2
	Vector<Real> xb(nt),yb(nt),xc(nt),yc(nt);
	geomPP.getarr("gaussianBump_xb", xb, 0, nt);
	geomPP.getarr("gaussianBump_xc", xc, 0, nt);
	geomPP.getarr("gaussianBump_yb", yb, 0, nt);
	geomPP.getarr("gaussianBump_yc", yc, 0, nt);
	for (int i = 0; i < nt; ++i)
	  {
	    b[i][0] = xb[i];
	    b[i][1] = yb[i];
	    c[i][0] = xc[i];
	    c[i][1] = yc[i];
	  }
#else
        MayDay::Error("beta_type = gaussianBump not implemented for CH_SPACEDIM > 2")
#endif
          basalFrictionPtr = static_cast<BasalFriction*>(new GaussianBumpFriction(t, C0, a, b, c));
      }
    else if (beta_type == "LevelData")
      {
        //read a one level beta^2 from an AMR Hierarchy, and  store it in a LevelDataBasalFriction
        ParmParse ildPP("inputLevelData");
        std::string infile;
        ildPP.get("frictionFile",infile);
        std::string frictionName = "btrc";
        ildPP.query("frictionName",frictionName);
        
        RefCountedPtr<LevelData<FArrayBox> > levelC (new LevelData<FArrayBox>());
        
        Real dx;
        
        Vector<RefCountedPtr<LevelData<FArrayBox> > > vectC;
        vectC.push_back(levelC);
        
        Vector<std::string> names(1);
        names[0] = frictionName;
        
        readLevelData(vectC,dx,infile,names,1);
	
        RealVect levelDx = RealVect::Unit * dx;
        basalFrictionPtr = static_cast<BasalFriction*>(new LevelDataBasalFriction(levelC,levelDx));
      }
#ifdef HAVE_PYTHON
    else if (beta_type == "Python")
      {
	ParmParse pyPP("PythonBasalFriction");
	std::string module;
	pyPP.get("module",module);
	std::string funcName = "friction";
	pyPP.query("function",funcName);
	basalFrictionPtr = static_cast<BasalFriction*>(new PythonInterface::PythonBasalFriction(module, funcName));
        
      }
#endif
    else 
      {
        MayDay::Error("undefined beta_type in inputs");
      }
  return basalFrictionPtr;
}

IceThicknessIBC* thicknessIBC(const RealVect& a_domainSize)
{
  // set IBC -- this includes initial ice thickness and basal geometry
  IceThicknessIBC* thicknessIBCPtr;

  ParmParse geomPP("geometry");
  std::string problem_type;
  geomPP.get("problem_type", problem_type);

  if (problem_type == "basic")
    {
      thicknessIBCPtr = new BasicThicknessIBC;
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
      
      thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
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
          RefCountedPtr<RealFunction<RealVect> > ptr(new SchoofBedrockElevation(a_domainSize[SpaceDim-2] * lengthScaleFactor,
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
	    RefCountedPtr<RealFunction<RealVect> > ptr(new KatzBedrockElevation(a_domainSize[SpaceDim-2],
										a_domainSize[SpaceDim-1],
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
      
      ibcPtr->setParameters(thicknessFunction, bedrockFunction, seaLevel);
      thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
    }
  else if (problem_type == "hump")
    {
      HumpIBC* ibcPtr = new HumpIBC;
      ParmParse humpPP("hump");
      
      Real maxThickness, radSqr, baseElevation, minThickness, seaLevel;
      RealVect center, widthScale;
      
      // default values to be equivalent to hump in Glimmer-CISM
      radSqr = 0.125*a_domainSize[0]*a_domainSize[1];
      maxThickness = 2000.0*pow(radSqr,0.5);
      baseElevation = 0.0;
      minThickness = 0.0;
      widthScale = RealVect::Unit;
      // this just lowers the sea level so that it's not relevant...
      seaLevel = -10.0;
      center = 0.5*a_domainSize;
      
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
      
      thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
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
      
      RefCountedPtr<LevelData<FArrayBox> > levelThck(new LevelData<FArrayBox>());
      RefCountedPtr<LevelData<FArrayBox> > levelTopg(new LevelData<FArrayBox>());
	 
      Real dx;
      
      Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
      vectData.push_back(levelThck);
      vectData.push_back(levelTopg);
      
      Vector<std::string> names(2);
      names[0] = thicknessName;
      names[1] = topographyName;
      readLevelData(vectData,dx,infile,names,1);
      
      RealVect levelDx = RealVect::Unit * dx;
      LevelDataIBC* ptr = new LevelDataIBC(levelThck,levelTopg,levelDx);
      thicknessIBCPtr = static_cast<IceThicknessIBC*>( ptr);
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
      
      PythonInterface::PythonIBC* ptr = new PythonInterface::PythonIBC(module, thckFuncName, topgFuncName);
      thicknessIBCPtr = static_cast<IceThicknessIBC*>( ptr);
    }
#endif
  else 
    {
      MayDay::Error("bad problem type");
    }

  return thicknessIBCPtr;
}

IceTemperatureIBC* temperatureIBC()
{
  IceTemperatureIBC* temperatureIBCPtr = NULL;
  ParmParse tempPP("temperature");
  std::string tempType("constant");
  tempPP.query("type",tempType);
  if (tempType == "constant")
    {
      Real T = 258.0;
      tempPP.query("value",T);
      ConstantIceTemperatureIBC* ptr = new ConstantIceTemperatureIBC(T);
      temperatureIBCPtr                 = static_cast<IceTemperatureIBC*>(ptr);
    }
  else if (tempType == "LevelData")
    {
      ParmParse ildPP("inputLevelData");
      LevelDataTemperatureIBC* ptr = NULL;
      CH_assert( (ptr = LevelDataTemperatureIBC::parse(ildPP)) != NULL);
      temperatureIBCPtr  = static_cast<IceTemperatureIBC*>(ptr);
    }
  else 
    {
      MayDay::Error("bad temperature type");
    }

  return temperatureIBCPtr;
}

RealVect readDomainSize()
{
  ParmParse pp2("main");
  
  Vector<Real> domSize(SpaceDim);
  pp2.getarr("domain_size", domSize, 0, SpaceDim);
  
  return RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));
}
