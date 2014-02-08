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
// cwrapper.cpp
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
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "ReadLevelData.H"
#include "PetscSolver.H"
#include "cwrapper.H"
#include "cdriverconstants.h"
#include "LevelDataSurfaceFlux.H"

struct BisiclesWrapper
{
  AmrIce m_amrIce;
  std::string m_input_fname;
  LevelDataSurfaceFlux* m_surface_flux;
  LevelDataSurfaceFlux* m_basal_flux;
  LevelDataSurfaceFlux* m_floating_ice_basal_flux;
  LevelDataSurfaceFlux* m_grounded_ice_basal_flux;
  BisiclesWrapper()
  {
    m_surface_flux = NULL;
    m_basal_flux = NULL;
    m_floating_ice_basal_flux = NULL;
    m_grounded_ice_basal_flux = NULL;
  }

  ~BisiclesWrapper()
  {

  }
};

namespace bisicles_c_wrapper
{
  std::map<int, BisiclesWrapper*> instances;
}

//fortran wrappers
void bisicles_new_instance_(int *instance_id,  char *input_fname, const int *len_fname)
  {
    input_fname[*len_fname - 1] = 0; // null terminate the string
    bisicles_new_instance(instance_id, input_fname);
  }
   void bisicles_free_instance_(int *instance_id)
  {
    bisicles_free_instance(instance_id);
  }
 void bisicles_set_2d_data_(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi)
  {
    bisicles_set_2d_data(instance_id,  data_ptr, field, 
			      dx, dims, boxlo, boxhi);
  }


  void bisicles_get_2d_data_(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			    const int *boxlo, const int *boxhi)
  {
    bisicles_get_2d_data(intance_id, data_ptr, field,
			      dx, dims, boxlo, boxhi);
    
  }
 void bisicles_init_instance_(int *instance_id)
  {
    bisicles_init_instance(instance_id);
  }
void bisicles_advance_(int *instance_id, double *max_time, int *max_step)
  {
    bisicles_advance(instance_id, max_time, max_step);
  }

// initialize the AmrIce object in a wrapper. Any surface
// fluxes specified in the wrapper will be added to whatever
// is specified in the input file (a_innputfile) 
// computes the initial solve / load from checkpoint etc
// such that the object is ready to be used in timestepping
void init_bisicles_instance( int argc, char *argv[], const char *a_inputfile, BisiclesWrapper& a_wrapper)
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
    { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
  
  ParmParse pp(argc-2,argv+2,NULL,a_inputfile);
  ParmParse pp2("main");

  //\todo Check how pout will work with multiple instances
  std::string poutBaseName = "pout";
  pp2.query("poutBaseName",poutBaseName);
  setPoutBaseName(poutBaseName);

  RealVect domainSize;
  Vector<Real> domSize(SpaceDim);
  pp2.getarr("domain_size", domSize, 0, SpaceDim);
  domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));

  AmrIce& amrObject = a_wrapper.m_amrIce;

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
  else if (rateFactorType == "patersonRate")
    {
      PatersonRateFactor rateFactor;
      ParmParse arPP("PatersonRate");
      arPP.query("epsSqr0", epsSqr0);
      amrObject.setRateFactor(&rateFactor);
      if (gfrPtr) 
	{
	  gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	}
    }
    else if (rateFactorType == "zwingerRate")
    {
      ZwingerRateFactor rateFactor;
      ParmParse arPP("ZwingerRate");
      arPP.query("epsSqr0", epsSqr0);
      amrObject.setRateFactor(&rateFactor);
      if (gfrPtr) 
	{
	  gfrPtr->setParameters(3.0 , &rateFactor, epsSqr0);
	}
    }
  amrObject.setConstitutiveRelation(constRelPtr);
 
  // -------------------------------------------------
  // set surface flux.
  // --------------------------------------------------

  SurfaceFlux* surf_flux_ptr = SurfaceFlux::parseSurfaceFlux("surfaceFlux");
  if ( surf_flux_ptr == NULL )
    {
      const std::string err("failed to parse surfaceFlux (maybe you have the old style surface_flux_type?");
      pout() << err << endl;
      MayDay::Error(err.c_str());
    }

  if (a_wrapper.m_surface_flux)
    {
      AxbyFlux* ptr = new AxbyFlux(1.0, a_wrapper.m_surface_flux, 1.0, surf_flux_ptr);	
      amrObject.setSurfaceFlux(ptr);
    }
  else
    {
      amrObject.setSurfaceFlux(surf_flux_ptr);
    }
  
 

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

  if ( a_wrapper.m_floating_ice_basal_flux || a_wrapper.m_grounded_ice_basal_flux)
    {
      //need to specify floating and grounded ice fluxes etc
      SurfaceFlux* f;
      if (a_wrapper.m_floating_ice_basal_flux == NULL)
	{
	  f = new zeroFlux();
	}
      else
	{
	  f = a_wrapper.m_floating_ice_basal_flux;
	}

      SurfaceFlux* g;
      if (a_wrapper.m_grounded_ice_basal_flux == NULL)
	{
	  g = new zeroFlux();
	}
      else
	{
	  g = a_wrapper.m_grounded_ice_basal_flux;
	}


      MaskedFlux* m = new MaskedFlux(g,f,f,g);
      AxbyFlux* ptr = new AxbyFlux(1.0, m , 1.0, basal_flux_ptr);
      amrObject.setBasalFlux(ptr); 

    }
  else if (a_wrapper.m_basal_flux)
    {
      AxbyFlux* ptr = new AxbyFlux(1.0, a_wrapper.m_basal_flux, 1.0, basal_flux_ptr);
      amrObject.setBasalFlux(ptr);
    }
  else
    {
      amrObject.setBasalFlux(basal_flux_ptr);
    }




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
  else if (beta_type == "MultiLevelData")
    {
      //read a multi level beta^2 from an AMR Hierarchy, and  store it in a MultiLevelDataBasalFriction
      ParmParse ildPP("inputLevelData");
      std::string infile;
      ildPP.get("frictionFile",infile);
      std::string frictionName = "btrc";
      ildPP.query("frictionName",frictionName);
      Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectC;
      Real dx;
      Vector<int> ratio;
      Vector<std::string> names(1);
      names[0] = frictionName;
      readMultiLevelData(vectC,dx,ratio,infile,names,1);
      RealVect dxCrse = RealVect::Unit * dx;
      basalFrictionPtr = static_cast<BasalFriction*>
	(new MultiLevelDataBasalFriction(vectC[0],dxCrse,ratio));
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
    
     
  amrObject.setDomainSize(domainSize);

  // set up initial grids, initialize data, etc. 
  amrObject.initialize();

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

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif


}  
void advance_bisicles_instance(BisiclesWrapper* wrapper_ptr, double max_time, int max_step  )
{
  if (wrapper_ptr != NULL)
    {
      wrapper_ptr->m_amrIce.run(max_time,max_step);
    }
}

namespace cwrapper
{

  /// On-the-fly definition of a DisjointBoxLayout which assumes level data arranged as one  rectangular box per processor
  /** assume that each process contributes one rectangular block of data and do it will call this function once per data field.
      If a processor has no data, then specify a box outside of the problem domain 
      we will need to think again if we want anything flash
  */
  void defineDBL(DisjointBoxLayout& a_dbl,  const int *dims, const int *boxlo, const int *boxhi)
  {
    IntVect lo; D_TERM(lo[0] = boxlo[0];, lo[1] = boxlo[1];, lo[2] = boxlo[2]);
    IntVect hi; D_TERM(hi[0] = boxhi[0];, hi[1] = boxhi[1];, hi[2] = boxhi[2]);

    Box domBox;
    if (procID() == uniqueProc(SerialTask::compute))
      {
	IntVect dlo; D_TERM(dlo[0] = 0;, dlo[1] = 0;, dlo[2] = 0);
	IntVect dhi; D_TERM(dhi[0] = dims[0] - 1;, dhi[1] = dims[1] - 1;, dhi[2] = dims[2]-1);
	domBox.define(Box(dlo,dhi));
      }
    broadcast(domBox,uniqueProc(SerialTask::compute));
   
    //first off, accumulate the boxes from all the processors
    Vector<Box> tboxes;
    gather(tboxes, Box(lo,hi),  uniqueProc(SerialTask::compute));
    Vector<int> tprocIDs;
    gather(tprocIDs, procID(), uniqueProc(SerialTask::compute));
    

    //now discard boxes outside the domain 
    Vector<Box> boxes;
    Vector<int> procIDs;
    if (procID() == uniqueProc(SerialTask::compute))
      {
	for (int i = 0; i < tboxes.size(); i++)
	  {
	    if (tboxes[i].intersects(domBox))
	      {
		boxes.push_back(tboxes[i]);
		procIDs.push_back(tprocIDs[i]);
	      }
	  }
      }

   
      
    broadcast(boxes,uniqueProc(SerialTask::compute));
    broadcast(procIDs,uniqueProc(SerialTask::compute));
    a_dbl.define(boxes, procIDs);

    
    pout() << "proc: " << procID() << " " << domBox << " ";
    for (int i = 0; i < procIDs.size(); i++)
      {
	pout() <<  procIDs[i] << " " << boxes[i] << endl;
      }

  }
}

void bisicles_set_2d_data
(BisiclesWrapper* wrapper_ptr,   double *data_ptr, const int *field, 
 const double *dx, const int *dims, const int *boxlo, const int *boxhi)
{
  if (wrapper_ptr != NULL)
    {
      DisjointBoxLayout dbl;
      cwrapper::defineDBL(dbl,dims,boxlo,boxhi);
     
      RefCountedPtr<LevelData<FArrayBox> > ptr(new LevelData<FArrayBox> (dbl, 1, IntVect::Zero));
      DataIterator dit(dbl);
      dit.reset();

      if (dit.ok())
	{
	  (*ptr)[dit].define(dbl[dit], 1, data_ptr) ;
	}

      RealVect dxv; D_TERM(dxv[0] = dx[0];, dxv[1] = dx[1];, dxv[2] = dx[2]);

      switch (*field)
	{
	case BISICLES_FIELD_SURFACE_FLUX:
	  wrapper_ptr->m_surface_flux = new LevelDataSurfaceFlux(ptr, dxv);
	  break;
	case BISICLES_FIELD_BASAL_FLUX:
	  wrapper_ptr->m_basal_flux = new LevelDataSurfaceFlux(ptr, dxv);
	  break;
	case BISICLES_FIELD_FLOATING_ICE_BASAL_FLUX:
	  wrapper_ptr->m_floating_ice_basal_flux = new LevelDataSurfaceFlux(ptr, dxv);
	  break;
	case BISICLES_FIELD_GROUNDED_ICE_BASAL_FLUX:
	  wrapper_ptr->m_grounded_ice_basal_flux = new LevelDataSurfaceFlux(ptr, dxv);
	  break;
	default: 
	  MayDay::Error("bisicles_get_2d_data: unknown (or unimplemented) field");
	}
      

    }
}




void bisicles_get_2d_data
(BisiclesWrapper* wrapper_ptr,   double *data_ptr, const int *field, 
 const double *dx, const int *dims, const int *boxlo, const int *boxhi)
{
  if (wrapper_ptr != NULL)
    {
      DisjointBoxLayout dbl;
      cwrapper::defineDBL(dbl,dims,boxlo,boxhi);
      RefCountedPtr<LevelData<FArrayBox> > ptr(new LevelData<FArrayBox> (dbl, 1, IntVect::Zero));

      DataIterator dit(dbl);
      dit.reset();
      if (dit.ok())
	{
	  (*ptr)[dit].define(dbl[dit], 1, data_ptr) ;
	}

      RealVect dxv; D_TERM(dxv[0] = dx[0];, dxv[1] = dx[1];, dxv[2] = dx[2]);
      AmrIce& amrIce = wrapper_ptr->m_amrIce;
      int n = amrIce.finestLevel() + 1;
      Vector<LevelData<FArrayBox>* > data(n);
      Vector<RealVect> amrDx(n);

      switch (*field)
	{
	case BISICLES_FIELD_SURFACE_ELEVATION:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getSurfaceHeight()));
	      amrDx[lev] = amrIce.dx(lev);
	    }

	  flattenCellData(*ptr,dxv,data,amrDx,true);
	  break;
	  
	case BISICLES_FIELD_BEDROCK_ELEVATION:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getTopography()));
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true);	
	  
	  break;

	case BISICLES_FIELD_ICE_THICKNESS:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getH()));
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true);	
	  
	  break;

	  

	default: 
	  MayDay::Error("bisicles_get_2d_data: unknown (or unimplemented) field");
	}
     

    }
}





void bisicles_new_instance(int *instance_id, const char *input_fname)
{
  //fake argc
#define NARGS 2

  BisiclesWrapper* ptr = new BisiclesWrapper;
  ptr->m_input_fname = input_fname;
  CH_assert(ptr != NULL);

  if  (bisicles_c_wrapper::instances.size() == 0)
    {
      *instance_id = 0;
    }
  else
    {
      *instance_id  = (--bisicles_c_wrapper::instances.end())->first + 1;
    }
  bisicles_c_wrapper::instances[*instance_id] = ptr;

  CH_assert(bisicles_c_wrapper::instances.size() == 1); //just one instance for now
  
}


void bisicles_init_instance(int *instance_id)
{
  if (instance_id)
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_id) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL)
	    {
#define NARGS 2
	      int argc = NARGS;
	      char *argv[NARGS];
	      char argv0[] = "cwrapper";
	      argv[0] = argv0;
	      char argv1[] = "drivel";
	      argv[1] = argv1;
	      init_bisicles_instance( argc, argv, i->second->m_input_fname.c_str(), *(i->second));
	    }
	}
    } 
}

void bisicles_free_instance(int *instance_id)
{
  if (instance_id)
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_id) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL)
	    {
	      delete i->second;
	      i->second = NULL;
	      bisicles_c_wrapper::instances.erase(i->first);
	    }
	}
    } 
}


void bisicles_advance(int *instance_id, double *max_time, int *max_step)
{
  if (instance_id && max_time && max_step)
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_id) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  advance_bisicles_instance(i->second, *max_time, *max_step);
	}
    }
}

void bisicles_set_2d_data(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi)
{
  if (instance_id) //\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_id) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  bisicles_set_2d_data(i->second, data_ptr, field, dx, dims, boxlo, boxhi);
	}
    }
}

void bisicles_get_2d_data(int *instance_id,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi)
{
  if (instance_id) //\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_id) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  bisicles_get_2d_data(i->second, data_ptr, field, dx, dims, boxlo, boxhi);
	}
    }
}
