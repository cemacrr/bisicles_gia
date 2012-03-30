
#include <iostream>
#include <fstream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "AmrIce.H"
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
#include "twistyStreamFriction.H"
#include "IceThicknessIBC.H"
#include "FortranInterfaceIBC.H"
#include "FortranInterfaceBasalFriction.H"
#include "SurfaceFlux.H"
#include "PiecewiseLinearFlux.H"
#include "bike_driver.H"


static int maxStep;
static Real maxTime;
static ParmParse* ppPtr;

//int bike_store(int obj_index, AmrIce * amr_object, int mode);

extern "C" void bike_driver_();

int
bike_store(int obj_index, Bike ** bike_object, int mode)
{
  static Bike * bike_store_ptr_arr[DYCORE_MODEL_COUNT];

  switch (mode) {
  case 0: bike_store_ptr_arr[obj_index] = *bike_object;
    cout << "In bike_store, mode = 0 -- Storing Bike Object # " 
	 << obj_index << ", Address = " << *bike_object << endl;
    break;
  case 1: *bike_object = bike_store_ptr_arr[obj_index];
    cout << "In bike_store, mode = 1 -- Retrieving Bike Object # " 
	 << obj_index << ", Address = " << *bike_object << endl;
    break;
  default: ;
  }
  return 0;
}

void bike_driver_run(int argc, int exec_mode);
 
/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
                         NUM_BETA_TYPES};

void bike_driver_init(int argc, int exec_mode,BisiclesToGlimmer * btg_ptr, const char * input_fname)
{ 

  char *argv[3];
  //char argv1[] = "/home/ranken/util/BISICLES/code/interface/inputs32.glimmer";
  // this one assumed we were running in gc1/parallel/src/fortran
  //char argv1[] = "../../../..//BISICLES/code/interface/inputs.glimmer";
  // this one assumes we're running in gc1/parallel/bin
  char argv1[] = "../../..//BISICLES/code/interface/inputs.glimmer";
  //char argv1[] = "inputs.glimmer";
  char argv0[] = "bike_driver";

  //  argv[0] = 
  argv[0] = argv0;
  argv[1] = argv1;


#ifdef MPI
  MPI_Init(&argc, &argv);
#endif
  
  cout << "In bike_driver..." << endl;
  
  
  { // Begin nested scope
    
     

#ifdef MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    //static AmrIce amrObject;
    //    AmrIce amrObject;

    Bike* bikePtr = new Bike();

    bikePtr->amrIce = new AmrIce;
    AmrIce* amrObjectPtr = bikePtr->amrIce;
    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    //    char* in_file = argv[1];
    const char* in_file = input_fname;
    cout << "Parsing: " << in_file << endl;
    //ParmParse pp(argc-2,argv+2,NULL,in_file);
    bikePtr->parmParse = new ParmParse(argc-2,argv+2,NULL,in_file);

    ParmParse pp2("main");

    ParmParse interfacePP("glimmerInterface");

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
        constRelPtr = new L1L2ConstitutiveRelation;
	gfrPtr = (dynamic_cast<L1L2ConstitutiveRelation*>(constRelPtr))->getGlensFlowRelationPtr();
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
	amrObjectPtr->setRateFactor(&rateFactor);
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
	amrObjectPtr->setRateFactor(&rateFactor);
	if (gfrPtr) 
	  {
	    gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	  }
      }


    //amrObject.setConstitutiveRelation(constRelPtr);
    amrObjectPtr -> setConstitutiveRelation(constRelPtr);  
    

    // ---------------------------------------------
    // set (upper) surface flux. 
    // ---------------------------------------------

    // for now, set it to be zero
    SurfaceFlux* surf_flux_ptr = NULL;
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
    else
      {
        MayDay::Error("invalid surface flux type");
      }
    
     amrObjectPtr -> setSurfaceFlux(surf_flux_ptr);

    // ---------------------------------------------
    // set basal (lower surface) flux. 
    // ---------------------------------------------
    
    SurfaceFlux* basal_flux_ptr = NULL;
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

	basal_flux_ptr = static_cast<SurfaceFlux*>
	  (new MaskedFlux(grounded_basal_flux_ptr->new_surfaceFlux(), 
			  floating_basal_flux_ptr->new_surfaceFlux()));

	delete grounded_basal_flux_ptr;
	delete floating_basal_flux_ptr;

	  
      }
    else
      {
	MayDay::Error("invalid basal flux type");
      }

    amrObjectPtr->setBasalFlux(basal_flux_ptr); 
   
    // ---------------------------------------------
    // set basal friction coefficient
    // ---------------------------------------------

    ParmParse geomPP("geometry");

    // ---------------------------------------------
    // set IBC -- this includes initial ice thickness, 
    // and basal geometry
    // ---------------------------------------------
   
    
    IceThicknessIBC* thicknessIBC = NULL;
    BasalFriction* basalFrictionPtr = NULL;
    std::string problem_type;
 
 
    RealVect domainSize;
    
    geomPP.get("problem_type", problem_type);
    if (problem_type =="fortran")
      {
        FortranInterfaceIBC* ibcPtr = new FortranInterfaceIBC;
        // need to set thickness and topography

        Vector<int> nGhost(SpaceDim, 0);
        IntVect ghostVect;
        interfacePP.queryarr("numGhost", nGhost, 0, SpaceDim);
        {
          ghostVect = IntVect(D_DECL(nGhost[0], nGhost[1], nGhost[2]));
        }

        bool nodal = true;
        interfacePP.query("nodalInitialData", nodal);
        
        // this is about removing ice from regions which
        // don't affect the dynamics of the region, but which 
        // can cause our solvers problems. Siple Island comes to mind here
        // for now, store these regions in a separate file. There's probably 
        // a better way to do this.
        
        // this will contain the boxes in the index space of the 
        // original data in which the thickness will be cleared.
        Vector<Box> clearBoxes;

        bool clearThicknessRegions = false;
        if (interfacePP.contains("clearThicknessRegionsFile"))
          {
            clearThicknessRegions = true;
            std::string clearFile;
            interfacePP.get("clearThicknessRegionsFile", clearFile);

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
            
            ibcPtr->setThicknessClearRegions(clearBoxes);
          }
        
        

        long * dimInfo;
        
        double * thicknessDataPtr, *topographyDataPtr;
        int i, reg_index;      

        dimInfo = btg_ptr -> getLongVar("dimInfo","geometry");
        thicknessDataPtr = btg_ptr -> getDoubleVar("thck","geometry");
        topographyDataPtr = btg_ptr -> getDoubleVar("topg","geometry");


        Real dew, dns;

        //dew = 1000.;
        //dns = 1000.;

        dew = *(btg_ptr -> getDoubleVar("dew","numerics"));
        dns = *(btg_ptr -> getDoubleVar("dns","numerics"));
	cout << "In bike_driver: dew, dns = " << dew << "  " << dns << endl;

        int * dimInfoGeom = new int[dimInfo[0]+1];    

        for (i=0;i<=dimInfo[0];i++) dimInfoGeom[i] = dimInfo[i];   
        cout << "DimInfoGeom  in bike_driver: " << endl;
        for (i=0;i<=dimInfoGeom[0];i++) cout << dimInfoGeom[i] << " ";
        cout << endl;

        CH_assert(SpaceDim == 2);
#if 1
        if (nodal)
          {
            // slc : thck and topg are defined at cell nodes on the glimmer grid,
            //       while bisicles needs them at cell centers. To get round this,
            //       decrement the grid dimensions by 1, then interpolate to the
            //       cell centers. This means  extra work, but is required
            //       if we think the domain boundaries should be  in the same place.
            dimInfoGeom[2] -= 1; dimInfoGeom[3] -=1;
          }
        domainSize[0] = dew*(dimInfoGeom[2]); 
        domainSize[1] = dns*(dimInfoGeom[3]);     
                 
        ibcPtr->setThickness(thicknessDataPtr, dimInfoGeom, &dew, &dns, 
                             ghostVect, nodal);                             
        ibcPtr->setTopography(topographyDataPtr, dimInfoGeom, &dew, &dns, 
                              ghostVect, nodal);
#else
	domainSize[0] = dew*(dimInfoGeom[2]); 
	domainSize[1] = dns*(dimInfoGeom[3]);	
        ibcPtr->setThickness(thicknessDataPtr, dimInfoGeom, &dew, &dns, ghostVect);
        ibcPtr->setTopography(topographyDataPtr, dimInfoGeom, &dew, &dns, ghostVect);
#endif
        // if desired, smooth out topography to fill in holes
        bool fillTopographyHoles = false;
        geomPP.query("fill_topography_holes", fillTopographyHoles);
        if (fillTopographyHoles)
          {
            Real holeVal = 0.0;
            geomPP.query("holeFillValue", holeVal);
            int numPasses = 1;
            geomPP.query("num_fill_passes", numPasses);

            for (int pass=0; pass<numPasses; pass++)
              {
                ibcPtr->fillTopographyHoles(holeVal);
              }
          }

        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);

	dimInfo = btg_ptr -> getLongVar("dimInfo","velocity");
        
        int * dimInfoVelo = new int[dimInfo[0]+1];
    
        for (i=0;i<=dimInfo[0];i++) dimInfoVelo[i] = dimInfo[i];      
        
        cout << "DimInfoVelo in bike_driver: " << endl;
        for (i=0;i<=dimInfoVelo[0];i++) cout << dimInfoVelo[i] << " ";
        cout << endl; 
	
	string beta_type;
	geomPP.get("beta_type", beta_type);
	if (beta_type == "fortran")
	  {
	    cout << "setting up basal friction coefficient " << endl;
	    double * basalTractionCoeffPtr = btg_ptr -> getDoubleVar("btrc","velocity"); 
	    FortranInterfaceBasalFriction* fibfPtr = new FortranInterfaceBasalFriction();
	    RealVect dx;
	    dx[0] = dew;
	    dx[1] = dns;

	    bool bfnodal = false;
	    interfacePP.query("nodalBasalFrictionData", bfnodal);
        
	    if (!bfnodal)
	      {
		// (DFM -- 6/19/11) -- friction is always cell-centered, 
		// regardless of thickness centering
		cout << "basal friction data is assumed to be cell-centered" << endl;
		fibfPtr->setReferenceFAB(basalTractionCoeffPtr, dimInfoVelo, dx, 
					 ghostVect, false);
	      }
	    else
	      {
		// (SLC -- 11/25/11) -- if we are taking glimmer's thickness to
		// be cell-centered, then its friction must be node centered 
		cout << "basal friction data will be averaged to cell-centers" << endl;
		fibfPtr->setReferenceFAB(basalTractionCoeffPtr, dimInfoVelo, dx, 
					 ghostVect, !nodal);
	      }
	    basalFrictionPtr = static_cast<BasalFriction*>(fibfPtr);
	  }
        
        // get Glimmer surface mass balance
        double * surfMassBalPtr;

        dimInfo = btg_ptr -> getLongVar("dimInfo","climate"); 

        int * dimInfoClim = new int[dimInfo[0]+1];

        for (i=0;i<=dimInfo[0];i++) dimInfoClim[i] = dimInfo[i];      
        surfMassBalPtr = btg_ptr -> getDoubleVar("acab","climate");
        cout << "DimInfoClim in bike_driver: " << endl;
        for (i=0;i<=dimInfoClim[0];i++) cout << dimInfoClim[i] << " ";
        cout << endl;      
        
      }
    else 
      {
        MayDay::Error("bad problem type");
      }


    // this lets us over-ride the Glimmer domain size
    // for the case where we're not using the entire glimmer domain
    // (normally due to the fact that they tend to choose really 
    //  odd (and un-coarsenable) domain sizes, so we often want 
    // to throw away a row or two of cells in each direction in 
    // order to make Multigrid have a chance of converging
    
    if (pp2.contains("domain_size"))
      {
        Vector<Real> domSize(SpaceDim);
        pp2.getarr("domain_size", domSize, 0, SpaceDim);
        domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));
      }

    amrObjectPtr -> setDomainSize(domainSize);
            
    // amrObject.setDomainSize(domainSize);
    amrObjectPtr -> setThicknessBC(thicknessIBC);

    
    

    std::string beta_type;
    geomPP.get("beta_type", beta_type);
    // read in type of beta^2 distribution
    
    if (beta_type == "constantBeta")
      {
        Real betaVal;
        geomPP.get("betaValue", betaVal);
        basalFrictionPtr = static_cast<BasalFriction*>(new constantFriction(betaVal));
      }
    else if (beta_type == "fortran")
      {
	// use the basal traction coefficient  btrc(x,y) supplied by glimmer
	// in this case, we already set up basalFrictionPtr
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
    
    else 
      {
        MayDay::Error("undefined beta_type in inputs");
      }
  
    // amrObject.setBasalFriction(basalFrictionPtr);
    amrObjectPtr -> setBasalFriction(basalFrictionPtr);   


    BasalFrictionRelation* basalFrictionRelationPtr = NULL;
    std::string basalFrictionRelType = "powerLaw";
    pp2.query("basalFrictionRelation", basalFrictionRelType);
    
    if (basalFrictionRelType == "powerLaw")
      {
	ParmParse plPP("BasalFrictionPowerLaw");

	Real m = 1.0;
	plPP.query("m",m);

	BasalFrictionPowerLaw*  pl = new BasalFrictionPowerLaw(m);
	basalFrictionRelationPtr = static_cast<BasalFrictionRelation*>(pl);
      }
    else
      {
	MayDay::Error("undefined basalFrictionRelation in inputs");
      }

    amrObjectPtr->setBasalFrictionRelation(basalFrictionRelationPtr);



    bike_store(btg_ptr -> getDyCoreIndex(), &bikePtr,0);

    // set up initial grids, initialize data, etc.
    cout << "Calling initialize..." << endl;  
    amrObjectPtr -> initialize();
    cout << "AMR object initialized." << endl;

    //Real startTime;

    pp2.get("maxTime", maxTime);
    pp2.get("maxStep", maxStep);

    // clean up
    pout () << "exec mode = " << exec_mode << endl;
    if (exec_mode == 0) {
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
       

    if (thicknessIBC != NULL)
      {
	delete thicknessIBC;
        thicknessIBC=NULL;
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
    } // if exec_mode == 2)
  }  

 

}

void bike_driver_run(BisiclesToGlimmer * btg_ptr, float cur_time_yr, float time_inc_yr)
{
  Bike *bikePtr;
  
  cout << "In bike_driver_run, cur_time, time_inc = " 
       << cur_time_yr << "   " << time_inc_yr << endl;
 
  bike_store(btg_ptr -> getDyCoreIndex(), &bikePtr ,1);

  AmrIce * amrObject_ptr = bikePtr->amrIce;
  cout << "Calling Amr.run..." << endl;

  amrObject_ptr -> run(maxTime, maxStep);

  // btg_ptr -> copyDoubleVar(thicknessPtr,"thck","geometry");

  cout << "Amr.run returned." << endl; 

#if 0
  //  reg_index = dycore_registry(0,1,&model_index,&dtg_ptr,-1);
  amrObject.getIceThickness(thicknessPtr, dim_info,
			    dew, dns);
#endif

}
  

void bike_driver_finalize(int amr_obj_index)
{
  Bike* bikePtr;

  pout() << "In bike_driver_finalize..." << endl;

  bike_store(amr_obj_index, &bikePtr, 1);
  
  if (bikePtr != NULL)
    {
      delete bikePtr; 
      bikePtr = NULL;
    }
  pout() << "Bike Object deleted." << endl << endl; 
#ifdef MPI
  MPI_Finalize();
#endif
  
  //  return 0;
}

Bike::~Bike()
{    
  if (amrIce != NULL)
    {
      delete amrIce;
      amrIce = NULL;
    }
  if (parmParse != NULL)
    {
      delete parmParse;
      parmParse = NULL;
    }
}
