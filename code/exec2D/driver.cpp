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
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
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
#include "PetscSolver.H"


/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
    guassianBump is used for the MISMIP3D perturbations tests.
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
			 gaussianBump,
                         NUM_BETA_TYPES};

//===========================================================================
// 2D shallow-shelf ice sheet model
//
//===========================================================================
int main(int argc, char* argv[]) {
  
  int ierr = 0;

#ifdef CH_USE_PETSC
  ierr = PetscInitialize(&argc, &argv,"./.petscrc",PETSC_NULL); CHKERRQ(ierr);
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

    RealVect domainSize;
    Vector<Real> domSize(SpaceDim);
    pp2.getarr("domain_size", domSize, 0, SpaceDim);
    domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));


    AmrIce amrObject;
    // ---------------------------------------------
    // set constitutive relation & rate factor
    // ---------------------------------------------
    std::string constRelType;
    pp2.get("constitutiveRelation", constRelType);
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
        L1L2ConstitutiveRelation* l1l2Ptr = new L1L2ConstitutiveRelation;
	l1l2Ptr->parseParameters();
	gfrPtr = l1l2Ptr->getGlensFlowRelationPtr();
        constRelPtr = l1l2Ptr;
      }
    else 
      {
        MayDay::Error("bad Constitutive relation type");
      }
    Real epsSqr0 = 1.0e-9;
    std::string rateFactorType = "constRate";
    pp2.query("rateFactor", rateFactorType);
    if (rateFactorType == "constRate")
      {
	ParmParse crPP("constRate");
	Real A = 9.2e-18;
	crPP.query("A", A);
	ConstantRateFactor rateFactor(A);

	crPP.query("epsSqr0", epsSqr0);
	amrObject.setRateFactor(&rateFactor);
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
	amrObject.setRateFactor(&rateFactor);
	if (gfrPtr) 
	  {
	    gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	  }
      }

    amrObject.setConstitutiveRelation(constRelPtr);
 
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

#if 0
    if (surf_flux_ptr == NULL)
      {
	// chunk for compatiblity with older input files
	MayDay::Warning("trying to parse old style surface_flux_type");
	surf_flux_ptr = NULL;
	std::string surfaceFluxType = "zeroFlux";
	pp2.query("surface_flux_type", surfaceFluxType);
	
	if (surfaceFluxType == "zeroFlux")
	  {
            surf_flux_ptr = new zeroFlux;
	  }
	else if (surfaceFluxType == "constantFlux")
	  {
	    constantFlux* constFluxPtr = new constantFlux;
	    Real fluxVal;
	    ParmParse ppFlux("constFlux");
	    ppFlux.get("flux_value", fluxVal);
	    constFluxPtr->setFluxVal(fluxVal);
	    
	    surf_flux_ptr = static_cast<SurfaceFlux*>(constFluxPtr);
	  }
      }
    
    if (surf_flux_ptr == NULL)
      {
	MayDay::Error("invalid surface flux type");
      }
#endif
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


#if 0    
    if (basal_flux_ptr == NULL)
      {

	//chunk for compatiblity with older input files
	MayDay::Warning("trying to parse old style basal_flux_type");
	std::string basalFluxType = "zeroFlux";
	pp2.query("basal_flux_type", basalFluxType);
	if (basalFluxType == "zeroFlux")
	  {
	    basal_flux_ptr = new zeroFlux;
	  }
	else if (basalFluxType == "constantFlux")
	  {
	    constantFlux* constFluxPtr = new constantFlux;
	    Real fluxVal;
	    ParmParse ppFlux("basalConstFlux");
	    ppFlux.get("flux_value", fluxVal);
	    constFluxPtr->setFluxVal(fluxVal);
	    basal_flux_ptr = static_cast<SurfaceFlux*>(constFluxPtr);
	    
	  }
	else if (basalFluxType == "maskedFlux")
	  {
	    SurfaceFlux* grounded_basal_flux_ptr = NULL;
	    std::string groundedBasalFluxType = "zeroFlux";
	    ParmParse ppbmFlux("basalMaskedFlux");
	    ppbmFlux.query("grounded_flux_type",groundedBasalFluxType);
	    if (groundedBasalFluxType == "zeroFlux")
	      {
		grounded_basal_flux_ptr = new zeroFlux;
	      }
	    else if (groundedBasalFluxType == "constantFlux")
	      {
		constantFlux* constFluxPtr = new constantFlux;
		Real fluxVal;
		ParmParse ppgFlux("groundedBasalConstFlux");
		ppgFlux.get("flux_value", fluxVal);
		constFluxPtr->setFluxVal(fluxVal);
		grounded_basal_flux_ptr = static_cast<SurfaceFlux*>(constFluxPtr);
	      }
	    else
	      {
		MayDay::Error("invalid grounded basal flux type");
	      }
	    
	    SurfaceFlux* floating_basal_flux_ptr = NULL;
	    std::string floatingBasalFluxType = "zeroFlux";
	    
	    ppbmFlux.query("floating_flux_type",floatingBasalFluxType);
	    if (floatingBasalFluxType == "zeroFlux")
	      {
		floating_basal_flux_ptr = new zeroFlux;
	      }
	    else if (floatingBasalFluxType == "constantFlux")
	      {
		constantFlux* constFluxPtr = new constantFlux;
		Real fluxVal;
		ParmParse ppfFlux("floatingBasalConstFlux");
		ppfFlux.get("flux_value", fluxVal);
		constFluxPtr->setFluxVal(fluxVal);
		floating_basal_flux_ptr = static_cast<SurfaceFlux*>(constFluxPtr);
	      }
	    else if (floatingBasalFluxType == "piecewiseLinearFlux")
	      {
		ParmParse pwlFlux("floatingBasalPWLFlux");
		int n = 1;  
		pwlFlux.query("n",n);
		Vector<Real> vabs(n,0.0);
		Vector<Real> vord(n,0.0);
		pwlFlux.getarr("abscissae",vabs,0,n);
		pwlFlux.getarr("ordinates",vord,0,n);
		PiecewiseLinearFlux* ptr = new PiecewiseLinearFlux(vabs,vord);
		floating_basal_flux_ptr = static_cast<SurfaceFlux*>(ptr);
	      }
	    else
	      {
		MayDay::Error("invalid floating basal flux type");
	      }
	    
	    SurfaceFlux* openland_basal_flux_ptr = new zeroFlux;
	    SurfaceFlux* opensea_basal_flux_ptr = new zeroFlux;
	    
	    basal_flux_ptr = static_cast<SurfaceFlux*>
	      (new MaskedFlux(grounded_basal_flux_ptr->new_surfaceFlux(), 
			      floating_basal_flux_ptr->new_surfaceFlux(),
			      opensea_basal_flux_ptr->new_surfaceFlux(),
			      openland_basal_flux_ptr->new_surfaceFlux()));
	    
	    delete grounded_basal_flux_ptr;
	    delete floating_basal_flux_ptr;
	    delete opensea_basal_flux_ptr;
	    delete openland_basal_flux_ptr;
	  }

      }
#endif    
    if (basal_flux_ptr == NULL)
      {
	MayDay::Error("invalid basal flux type");
      }
    
    amrObject.setBasalFlux(basal_flux_ptr); 

    // ---------------------------------------------
    // set mu coefficient
    // ---------------------------------------------
    ParmParse muPP("muCoefficient");
    std::string muCoefType = "unit";
    muPP.query("type",muCoefType );
    if (muCoefType == "unit")
      {
	MuCoefficient* ptr = static_cast<MuCoefficient*>(new UnitMuCoefficient());
	amrObject.setMuCoefficient(ptr);
	delete ptr;
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
	 MuCoefficient* ptr = static_cast<MuCoefficient*>
	   (new LevelDataMuCoefficient(levelMuCoef,levelDx));
	 amrObject.setMuCoefficient(ptr);
	 delete ptr;
      }
    else
      {
	MayDay::Error("undefined MuCoefficient in inputs");
      }


    // ---------------------------------------------
    // set basal friction coefficient and relation
    // ---------------------------------------------

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
                                                                              domainSize));
      }
    // keep this one around for backward compatibility, even if it
    // is a special case of sinusoidalBeta
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
                                                                              domainSize));

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
                                                                                domainSize));          
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
        basalFrictionPtr = static_cast<BasalFriction*>
	  (new GaussianBumpFriction(t, C0, a, b, c));
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
	 basalFrictionPtr = static_cast<BasalFriction*>
	   (new LevelDataBasalFriction(levelC,levelDx));
       }
#ifdef HAVE_PYTHON
    else if (beta_type == "Python")
      {
	ParmParse pyPP("PythonBasalFriction");
	std::string module;
	pyPP.get("module",module);
	std::string funcName = "friction";
	pyPP.query("function",funcName);
	basalFrictionPtr = static_cast<BasalFriction*>
	  (new PythonInterface::PythonBasalFriction(module, funcName));

      }
#endif
    else 
      {
        MayDay::Error("undefined beta_type in inputs");
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

	 PythonInterface::PythonIBC* ptr = new PythonInterface::PythonIBC(module, thckFuncName, topgFuncName);
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
	CH_assert( (ptr = LevelDataTemperatureIBC::parse(ildPP)) != NULL);
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


