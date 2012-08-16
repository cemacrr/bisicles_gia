#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

using std::ifstream; 
using std::ios;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;


#include "Box.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "parstream.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "FineInterp.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "AmrIce.H"
#include "computeNorm.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "PiecewiseLinearFillPatch.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "DivergenceF_F.H"
#include "computeSum.H"
#include "CONSTANTS.H"
#include "IceConstants.H"
#include "ExtrapBCF_F.H"
#include "amrIceF_F.H"
#include "BisiclesF_F.H"
#include "PicardSolver.H"
#include "JFNKSolver.H"
#include "PetscIceSolver.H"
#include "RelaxSolver.H"
#include "KnownVelocitySolver.H"
#include "VCAMRPoissonOp2.H"
#include "AMRPoissonOpF_F.H"
#include "CH_HDF5.H"
#include "IceVelocity.H"
#include "LevelMappedDerivatives.H"
#include "NamespaceHeader.H"

// small parameter defining when times are equal
#define TIME_EPS 1.0e-12

int AmrIce::s_verbosity = 1;

#if 0
void zeroBCValue(Real* pos,
                 int* dir,
                 Side::LoHiSide* side,
                 Real* a_values)
{
  a_values[0]=0.0;
}


void iceNeumannBC(FArrayBox& a_state,
                  const Box& a_valid,
                  const ProblemDomain& a_domain,
                  Real a_dx,
                  bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}


void iceDirichletBC(FArrayBox& a_state,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}

#endif


/// diagnostic function -- integrates thickness over domain
Real
AmrIce::computeTotalIce() const
{
  Vector<LevelData<FArrayBox>* > thickness(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      // need a const_cast to make things all line up right
      // (but still essentially const)
      thickness[lev] = const_cast<LevelData<FArrayBox>* >(&levelCoords.getH());
    }

  Interval thicknessInt(0,0);
  Real totalIce = computeSum(thickness, m_refinement_ratios,
                             m_amrDx[0], thicknessInt, 0);


  return totalIce;

}



/// fill flattened Fortran array of data with ice thickness
void
AmrIce::getIceThickness(Real* a_data_ptr, int* a_dim_info, 
			Real* a_dew, Real* a_dns) const
{
    // dimInfo is (SPACEDIM, nz, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering
  IntVect hiVect(D_DECL(a_dim_info[2]-1,a_dim_info[3]-1, a_dim_info[1]-1));
  Box fabBox(IntVect::Zero, hiVect);

  FArrayBox exportHfab(fabBox, 1, a_data_ptr); 

  // now pack this into a LevelData  -- this will need to be expanded in MPI
  Vector<Box> exportBoxes(1,fabBox);
  Vector<int> procAssign(1,0);
  // ignore periodicity, since we don't have ghost cells
  DisjointBoxLayout exportGrids(exportBoxes, procAssign);
  LevelData<FArrayBox> exportLDF(exportGrids, 1);

  // this isn't correct in 3d...
  CH_assert(SpaceDim != 3);
  RealVect exportDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // assume that dx = dy, at least for now
  CH_assert (exportDx[0] == exportDx[1]);

  // start at level 0, then work our way up to finest level, 
  // over-writing as we go. An optimzation would be to check 
  // to see if finer levels cover the entire domain...
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      const LevelSigmaCS& levelCS = *(m_vect_coordSys[lev]);
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      const LevelData<FArrayBox>& levelThickness = levelCS.getH();

      // refinement ratio
      Real refRatio = exportDx[0]/m_amrDx[lev];

      Real tolerance = 1.0e-6;

      if (refRatio > 1.0 + tolerance)
        {
          // current level finer than export level -- average solution
          // onto output
          int nRef = (int)(refRatio + tolerance);

          CoarseAverage averager(levelGrids,exportGrids,
                                 1, nRef);
          averager.averageToCoarse(exportLDF, levelThickness);
        }
      else if (refRatio < 1.0-tolerance)
        {
          // current level coarser than export level -- interpolate solution
          int nRef = (int)(1.0/refRatio + tolerance);
          
          // FineInterp needs a problem domain
          ProblemDomain exportDomain(m_amrDomains[lev]);
          exportDomain.refine(nRef);
                                     

          FineInterp interpolator(exportGrids, 1, nRef, exportDomain);
          interpolator.interpToFine(exportLDF, levelThickness);
        }
      else
        {
          // same resolution -- just do a copyTo
          levelThickness.copyTo(exportLDF);
        }

    } // end loop over levels

  // now copy to input fab
  DataIterator exportDit = exportLDF.dataIterator();
  for (exportDit.begin(); exportDit.ok(); ++exportDit)
    {
      // in parallel, this should only be on proc 0
      exportHfab.copy(exportLDF[exportDit]);
    }
  
}

bool 
AmrIce::isDefined() const
{
  return m_is_defined;
}


AmrIce::AmrIce() : m_velSolver(NULL),
                   m_constitutiveRelation(NULL),
		   m_basalFrictionRelation(NULL),
                   m_thicknessPhysPtr(NULL),
                   m_thicknessIBCPtr(NULL), 
                   m_surfaceFluxPtr(NULL),
		   m_basalFluxPtr(NULL),
                   m_basalFrictionPtr(NULL)
{
  setDefaults();
}

void
AmrIce::setDefaults()
{
  m_sigmaSet = false;
  // set some bogus values as defaults 
  m_is_defined = false;
  m_max_level = -1;
  m_finest_level = -1;
  m_finest_timestep_level = -1;
  m_tag_cap = 100;
  m_block_factor = -1;
  m_fill_ratio = -1;
  m_do_restart = false;
  m_restart_step = -1;
  //  m_constitutiveRelation = NULL;
  m_solverType = Picard;
  // at the moment, 1 is the only one which works
  m_temporalAccuracy = 1;
  m_num_thickness_ghost = 4;
  // default is -1, which means use the solver's own defaults
  m_maxSolverIterations = -1;
  
  m_velocity_solver_tolerance = 1e-10;

  //m_velSolver = NULL;
  m_domainSize = -1*RealVect::Unit;

  //m_beta_type = constantBeta;
  m_betaVal = 1000.0;
  m_betaEps = 0.0;
  m_basalSlope = RealVect::Zero;
  
  m_interpolate_zb = true;

  // set the rest of these to reasonable defaults
  m_nesting_radius = 1;
  m_tagOnGradVel = false;
  m_tagging_val = 0.1;
  m_tagOnLapVel = false;
  m_tagOnGroundedLapVel = false;
  m_laplacian_tagging_val = 1.0;
  m_laplacian_tagging_max_basal_friction_coef = 1.2345678e+300;
  m_tagOnEpsSqr = false;  
  m_epsSqr_tagVal =0.1;
  m_tagOnVelRHS = false;
  m_velRHS_tagVal = 1.0;
  m_tagOndivHgradVel = false;
  m_divHGradVel_tagVal = 1.0;
  m_tagGroundingLine  = false;
  m_tagMargin  = false;
  m_tagAllIce  = false;
  m_groundingLineTaggingMinVel = 200.0;
  m_groundingLineTaggingMaxBasalFrictionCoef = 1.2345678e+300;
  
  m_tags_grow = 1;
  m_cfl = 0.25;
  m_max_dt_grow = 1.5;
  m_dt = 1.0e20;
  m_max_box_size = 64;
  m_max_base_grid_size = -1;
  m_isothermal = true;
  m_isothermalTemperature = 238.5;
  m_iceDensity = 910.0; 
  m_seaWaterDensity = 1028.0;
  m_gravity = 9.81;

  m_report_grounded_ice = false;

  m_plot_prefix = "plot";
  m_plot_interval = 10000000;
  m_plot_time_interval = 1.0e+12;
  m_write_presolve_plotfiles = false;
  m_write_solver_rhs = false;
  m_write_dHDt = true;
  m_write_fluxVel = true;
  m_write_viscousTensor = false;
  m_write_baseVel = true;
  m_write_solver_rhs = false;
  m_write_temperature = false;
  m_write_map_file = false;
  m_write_thickness_sources = false;
  m_write_layer_velocities = false;
 
  m_check_prefix = "chk";
  m_check_interval = -1;
  m_check_overwrite = true;
  m_check_exit = false;

  m_diffusionTreatment = NONE;
  m_additionalDiffusivity = 0.0;
  m_additionalVelocity = false;
  m_timeStepTicks = false;
  m_fixed_dt  = zero;

  m_limitVelRHS = false;
  m_gradLimitRadius = 10;
  
  m_limitFluxSpeed = -1.0;

  m_velocitySolveInitialResidualNorm = -1.0;
  m_doInitialVelSolve = true; 
  m_doInitialVelGuess = false; 
  m_initialGuessType = SlidingLaw;
  m_initialGuessConstMu = 1.0e+7; // a number that seems plausible, only needed
                                  // if m_initialGuessType = ConstMu
  m_initialGuessSolverType = Picard; 
  m_initialGuessConstVel = RealVect::Zero; // only needed if m_initialGuessType = ConstMu *and* the basal traction relation is nonlinear

  m_basalFrictionDecay = -1.0 ; // set basal friction to zero wherever ioce is floating
  m_basalLengthScale = 0.0; // don't mess about with the basal friction / rhs by default
  
  m_wallDrag = false; //by default, don't compute  additional drag due to contact with rocky walls
  m_wallDragExtra = 0.0; // by default, assume wall drag proportional to basal drag

  m_groundingLineRegularization = 0.0;
  m_groundingLineCorrection = true;
  m_evolve_thickness = true;
  m_grounded_ice_stable = false;
  m_floating_ice_stable = false;
  
  m_groundingLineProximityScale = 1.0e+4;

  //cache validity flags
  m_A_valid = false;
  m_groundingLineProximity_valid = false;
  m_viscousTensor_valid = false;

  constantFlux* cfptr = new constantFlux;
  cfptr->setFluxVal(0.0);
  m_basalFluxPtr = cfptr;
  
  m_calvingModelPtr = NULL;//new NoCalvingModel; 

  
  m_offsetTime = 0.0;

}

AmrIce::~AmrIce()
{

  // clean up memory
  for (int lev=0; lev<m_velocity.size(); lev++)
    {
      if (m_velocity[lev] != NULL)
        {
          delete m_velocity[lev];
          m_velocity[lev] = NULL;
        }
    }
 

  // clean up velocityRHS storage if appropriate
  for (int lev=0; lev<m_velRHS.size(); lev++)
    {
      if (m_velRHS[lev] != NULL)
        {
          delete m_velRHS[lev];
          m_velRHS[lev] = NULL;
        }
    }
 
  // clean up basal C storage if appropriate
  for (int lev=0; lev < m_velBasalC.size(); lev++)
    {
      if (m_velBasalC[lev] != NULL)
        {
          delete m_velBasalC[lev];
          m_velBasalC[lev] = NULL;
        }
    }

  // clean up cell centered mu coefficient storage if appropriate
  for (int lev=0; lev < m_cellMuCoef.size(); lev++)
    {
      if (m_cellMuCoef[lev] != NULL)
        {
          delete m_cellMuCoef[lev];
          m_cellMuCoef[lev] = NULL;
        }
    }
  // clean up face-centered mu coef storage if appropriate
  for (int lev=0; lev < m_faceMuCoef.size(); lev++)
    {
      if (m_faceMuCoef[lev] != NULL)
        {
          delete m_faceMuCoef[lev];
          m_faceMuCoef[lev] = NULL;
        }
    }


  // clean up face velocity storage if appropriate
  for (int lev=0; lev < m_faceVel.size(); lev++)
    {
      if (m_faceVel[lev] != NULL)
        {
          delete m_faceVel[lev];
          m_faceVel[lev] = NULL;
        }
    }
  // clean up surface thickness storage if appropriate
  for (int lev=0; lev < m_surfaceThicknessSource.size(); lev++)
    {
      if (m_surfaceThicknessSource[lev] != NULL)
	{
	  delete m_surfaceThicknessSource[lev];
	  m_surfaceThicknessSource[lev] = NULL;
	}
      if (m_basalThicknessSource[lev] != NULL)
	{
	  delete m_basalThicknessSource[lev];
	  m_basalThicknessSource[lev] = NULL;
	}
    }
  
  for (int lev=0; lev < m_balance.size(); lev++)
    {
      if (m_balance[lev] != NULL)
	{
	  delete m_balance[lev];
	  m_balance[lev] = NULL;
	}
    }

#if BISCICLES_Z == BISICLES_LAYERED
  for (int lev=0; lev < m_layerSFaceXYVel.size(); lev++)
    {
      if (m_layerSFaceXYVel[lev] != NULL)
        {
          delete m_layerSFaceXYVel[lev];
          m_layerSFaceXYVel[lev] = NULL;
        }
    }

  for (int lev=0; lev < m_layerXYFaceXYVel.size(); lev++)
    {
      if (m_layerXYFaceXYVel[lev] != NULL)
        {
          delete m_layerXYFaceXYVel[lev];
          m_layerXYFaceXYVel[lev] = NULL;
        }
    }

#endif
  

  // clean up temperature  storage if appropriate
  for (int lev=0; lev < m_temperature.size(); lev++)
    {
      if (m_temperature[lev] != NULL)
        {
          delete m_temperature[lev];
          m_temperature[lev] = NULL;
        }
    }
  for (int lev=0; lev < m_A.size(); lev++)
    {
      if (m_A[lev] != NULL)
        {
          delete m_A[lev];
          m_A[lev] = NULL;
        }
      
      if (m_sA[lev] != NULL)
        {
          delete m_sA[lev];
          m_sA[lev] = NULL;
        }

    }
#if BISCICLES_Z == BISICLES_LAYERED

  for (int lev=0; lev < m_sA.size(); lev++)
    {
      if (m_sA[lev] != NULL)
        {
          delete m_sA[lev];
          m_sA[lev] = NULL;
        }
    }
  for (int lev=0; lev < m_sTemperature.size(); lev++)
    {
      if (m_sTemperature[lev] != NULL)
	{
	  delete m_sTemperature[lev];
	  m_sTemperature[lev] = NULL;
	}
    }
  
  for (int lev=0; lev < m_bTemperature.size(); lev++)
    {
      if (m_bTemperature[lev] != NULL)
        {
          delete m_bTemperature[lev];
          m_bTemperature[lev] = NULL;
        }
    }
  
  for (int lev=0; lev < m_bA.size(); lev++)
    {
      if (m_bA[lev] != NULL)
	{
	  delete m_bA[lev];
	  m_bA[lev] = NULL;
	}
    }
#endif

  for (int lev = 0; lev < m_old_thickness.size(); lev++) 
    {
      if (m_old_thickness[lev] != NULL) 
        {
          delete m_old_thickness[lev];
          m_old_thickness[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_groundingLineProximity.size(); lev++) 
    {
      if (m_groundingLineProximity[lev] != NULL) 
        {
          delete m_groundingLineProximity[lev];
          m_groundingLineProximity[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_dragCoef.size(); lev++) 
    {
      if (m_dragCoef[lev] != NULL) 
        {
          delete m_dragCoef[lev];
          m_dragCoef[lev] = NULL;
        }
    }
  
  for (int lev = 0; lev < m_viscosityCoefCell.size(); lev++) 
    {
      if (m_viscosityCoefCell[lev] != NULL) 
        {
          delete m_viscosityCoefCell[lev];
          m_viscosityCoefCell[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_viscousTensorCell.size(); lev++) 
    {
      if (m_viscousTensorCell[lev] != NULL) 
        {
          delete m_viscousTensorCell[lev];
          m_viscousTensorCell[lev] = NULL;
        }
    }
  for (int lev = 0; lev < m_viscousTensorFace.size(); lev++) 
    {
      if (m_viscousTensorFace[lev] != NULL) 
        {
          delete m_viscousTensorFace[lev];
          m_viscousTensorFace[lev] = NULL;
        }
    }



  if (m_velSolver != NULL)
    {
      // code crashes here. comment out until I figure out the problem...
      delete m_velSolver;
      m_velSolver = NULL;
    }

  if (m_constitutiveRelation != NULL)
    {
      delete m_constitutiveRelation;
      m_constitutiveRelation = NULL;
    }
  if (m_basalFrictionRelation != NULL)
    {
      delete m_basalFrictionRelation;
      m_basalFrictionRelation = NULL;
    }
  
  for (int lev=0; lev<m_thicknessPatchGodVect.size(); lev++)
    {
      if (m_thicknessPatchGodVect[lev] != NULL)
	{
	  delete m_thicknessPatchGodVect[lev];
	  m_thicknessPatchGodVect[lev] = NULL;
	}
    }


  if (m_thicknessPhysPtr != NULL)
    {
      delete m_thicknessPhysPtr;
      m_thicknessPhysPtr = NULL;
    }

  if (m_thicknessIBCPtr != NULL)
    {
      delete m_thicknessIBCPtr;
      m_thicknessIBCPtr = NULL;
    }

  for (int lev=0; lev<m_temperaturePatchGodVect.size(); lev++)
    {
      if (m_temperaturePatchGodVect[lev] != NULL)
	{
	  delete m_temperaturePatchGodVect[lev];
	  m_temperaturePatchGodVect[lev] = NULL;
	}
    }

  if (m_temperaturePhysPtr != NULL)
    {
      // this is segfaulting...
      delete m_temperaturePhysPtr;
      m_temperaturePhysPtr = NULL;
    }

  if (m_temperatureIBCPtr != NULL)
    {
      delete m_temperatureIBCPtr;
      m_temperatureIBCPtr = NULL;
    }

  

  if (m_surfaceFluxPtr != NULL)
    {
      delete m_surfaceFluxPtr;
      m_surfaceFluxPtr = NULL;
    }
  if (m_basalFluxPtr != NULL)
    {
      delete m_basalFluxPtr;
      m_basalFluxPtr = NULL;
    }
  if (m_calvingModelPtr != NULL)
    {
      delete m_calvingModelPtr;
      m_calvingModelPtr = NULL;
    }
  if (m_basalFrictionPtr != NULL)
    {
      delete m_basalFrictionPtr;
      m_basalFrictionPtr = NULL;
    }

  // that should be it!

}


void
AmrIce::initialize()
{

  if (s_verbosity > 3) 
    {
      pout() << "AmrIce::initialize" << endl;
    }

  // set time to be 0 -- do this now in case it needs to be changed later
  m_time = 0.0;
  m_cur_step = 0;
 
  // first, read in info from parmParse file
  ParmParse ppCon("constants");
  ppCon.query("ice_density",m_iceDensity);
  ppCon.query("sea_water_density",m_seaWaterDensity);
  ppCon.query("gravity",m_gravity);
  ParmParse ppAmr("amr");
  Vector<int> ancells(3); 
  // slc : SpaceDim == 2 implies poor-mans multidum mode, in which we still
  // care about the number of vertical layers. 
  Vector<Real> domsize(SpaceDim);

  // assumption is that domains are not periodic
  bool is_periodic[SpaceDim];
  for (int dir=0; dir<SpaceDim; dir++)
    is_periodic[dir] = false;
  Vector<int> is_periodic_int(SpaceDim, 0);

  ppAmr.get("maxLevel", m_max_level);
  ppAmr.query("tagCap",m_tag_cap);
  
  ppAmr.getarr("num_cells", ancells, 0, ancells.size());

# if 0
  // this is now handled in main and passed in using the
  // setDomainSize function
  if (ppAmr.contains("domain_size"))
    {
      ppAmr.getarr("domain_size", domsize, 0, SpaceDim);
      m_domainSize = RealVect(D_DECL(domsize[0], domsize[1], domsize[2]));
    }
  
#endif
  
 

  ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++) 
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }

  ppAmr.query("cfl", m_cfl);

  m_initial_cfl = m_cfl;
  ppAmr.query("initial_cfl", m_initial_cfl);

  ppAmr.query("max_dt_grow_factor", m_max_dt_grow);

  ppAmr.query("time_step_ticks", m_timeStepTicks);
  // max_dt_grow must be at least two in this case - or it will never grow
  if (m_timeStepTicks)
    m_max_dt_grow = std::max(m_max_dt_grow, two);

  ppAmr.query("fixed_dt", m_fixed_dt);
  ppAmr.query("offsetTime", m_offsetTime);


  ppAmr.query("isothermal",m_isothermal);
  ppAmr.query("isothermalTemperature",m_isothermalTemperature);

  ppAmr.query("plot_interval", m_plot_interval);
  ppAmr.query("plot_time_interval", m_plot_time_interval);
  ppAmr.query("plot_prefix", m_plot_prefix);

  ppAmr.query("write_map_file", m_write_map_file);

  ppAmr.query("write_preSolve_plotfiles", m_write_presolve_plotfiles);

  ppAmr.query("write_flux_velocities", m_write_fluxVel);
  ppAmr.query("write_viscous_tensor", m_write_viscousTensor);
  ppAmr.query("write_base_velocities", m_write_baseVel);
  ppAmr.query("write_temperature", m_write_temperature);
  ppAmr.query("write_thickness_sources", m_write_thickness_sources);
  ppAmr.query("write_layer_velocities", m_write_layer_velocities);

  ppAmr.query("evolve_thickness", m_evolve_thickness);
  ppAmr.query("grounded_ice_stable", m_grounded_ice_stable);
  ppAmr.query("floating_ice_stable", m_floating_ice_stable);

  ppAmr.query("check_interval", m_check_interval);
  ppAmr.query("check_prefix", m_check_prefix);
  ppAmr.query("check_overwrite", m_check_overwrite);
  ppAmr.query("check_exit", m_check_exit);

  bool tempBool = m_write_dHDt;
  ppAmr.query("write_dHDt", tempBool);
  m_write_dHDt = (tempBool == 1);

  ppAmr.query("write_solver_rhs", m_write_solver_rhs);

  if (m_max_level > 0)
    {
      //m_refinement_ratios.resize(m_max_level+1,-1);
      ppAmr.getarr("ref_ratio", m_refinement_ratios, 0, m_max_level);
    }
  else
    {
      m_refinement_ratios.resize(1);
      m_refinement_ratios[0] = -1;
    }

  ppAmr.query("verbosity", s_verbosity);

  ppAmr.get("regrid_interval", m_regrid_interval);
  m_n_regrids = 0;
  
  ppAmr.query("interpolate_zb", m_interpolate_zb);

  ppAmr.get("blockFactor", m_block_factor);

  // int n_tag_subset_boxes = 0;
  // m_tag_subset.define();
  // ppAmr.query("n_tag_subset_boxes",n_tag_subset_boxes);
  // if ( n_tag_subset_boxes > 0)
  //   {
      
  //     Vector<int> corners(2*SpaceDim*n_tag_subset_boxes);
  //     ppAmr.getarr("tag_subset_boxes", corners, 0, corners.size());
  //     int j = 0;
  //     for (int i =0; i < n_tag_subset_boxes; i++)
  // 	{
  // 	  IntVect small;
  // 	  for (int dir = 0; dir < SpaceDim; ++dir)
  // 	    {
  // 	      small[dir] = corners[j++];
  // 	    }
  // 	  IntVect big;
  // 	  for (int dir = 0; dir < SpaceDim; ++dir)
  // 	    {
  // 	      big[dir] = corners[j++];
  // 	    }
  // 	  m_tag_subset |= Box(small,big);
  // 	}
  //   }

 
 


  ppAmr.get("fill_ratio", m_fill_ratio);

  ppAmr.query("nestingRadius", m_nesting_radius);
  
  bool isThereATaggingCriterion = false;
  ppAmr.query("tag_on_grad_velocity", m_tagOnGradVel);
  isThereATaggingCriterion |= m_tagOnGradVel;

  ppAmr.query("tagging_val", m_tagging_val);

  ppAmr.query("tag_on_laplacian_velocity", m_tagOnLapVel);
  isThereATaggingCriterion |= m_tagOnLapVel;

  ppAmr.query("tag_on_grounded_laplacian_velocity", m_tagOnGroundedLapVel);
  isThereATaggingCriterion |= m_tagOnGroundedLapVel;

  // if we set either of these to be true, require that we also provide the threshold
  if (m_tagOnLapVel | m_tagOnGroundedLapVel)
    {
      ppAmr.get("lap_vel_tagging_val", m_laplacian_tagging_val);
      ppAmr.query("lap_vel_tagging_max_basal_friction_coef", m_laplacian_tagging_max_basal_friction_coef);
    }


  ppAmr.query("tag_on_strain_rate_invariant",m_tagOnEpsSqr);
  isThereATaggingCriterion |= m_tagOnEpsSqr;
  // if we set this to be true, require that we also provide the threshold
  if (m_tagOnEpsSqr)
    {
      ppAmr.get("strain_rate_invariant_tagging_val", m_epsSqr_tagVal);
    }


  ppAmr.query("tag_on_velocity_rhs",m_tagOnVelRHS);
  isThereATaggingCriterion |= m_tagOnVelRHS;

  // if we set this to be true, require that we also provide the threshold
  if (m_tagOnVelRHS)
    {
      ppAmr.get("velocity_rhs_tagging_val", m_velRHS_tagVal);
    }
  
  ppAmr.query("tag_grounding_line", m_tagGroundingLine);
  isThereATaggingCriterion |= m_tagGroundingLine;
  // if we set this to be true, require that we also provide the threshold
  if (m_tagGroundingLine)
    {
      ppAmr.get("grounding_line_tagging_min_vel",m_groundingLineTaggingMinVel);
      ppAmr.query("grounding_line_tagging_max_basal_friction_coef", m_groundingLineTaggingMaxBasalFrictionCoef);
    }
  
  
  ppAmr.query("tag_ice_margin", m_tagMargin);
  isThereATaggingCriterion |= m_tagMargin;

  ppAmr.query("tag_all_ice", m_tagAllIce);
  isThereATaggingCriterion |= m_tagAllIce;


  ppAmr.query("tag_on_div_H_grad_vel",m_tagOndivHgradVel);
  isThereATaggingCriterion |= m_tagOndivHgradVel;

  // if we set this to be true, require that we also provide the threshold
  if (m_tagOndivHgradVel)
    {
      ppAmr.get("div_H_grad_vel_tagging_val", m_divHGradVel_tagVal);
    }


  // here is a good place to set default to grad(vel)
  //if ((!m_tagOnGradVel) && (!m_tagOnLapVel) && (!m_tagOnEpsSqr))
  if (!isThereATaggingCriterion)
    {
      m_tagOnGradVel = true;
    }
  
  ppAmr.query("tags_grow", m_tags_grow);

  ppAmr.query("max_box_size", m_max_box_size);

  if (ppAmr.contains("max_base_grid_size") )
    {
      ppAmr.get("max_base_grid_size", m_max_base_grid_size);
    }
  else 
    {
      m_max_base_grid_size = m_max_box_size;
    }

  ppAmr.query("report_sum_grounded_ice",   m_report_grounded_ice);
  
  // get temporal accuracy
  ppAmr.query("temporal_accuracy", m_temporalAccuracy);

  // number of ghost cells depends on what scheme we're using
  if (m_temporalAccuracy < 3)
    {
      m_num_thickness_ghost = 4;
    }
  else 
    {
      m_num_thickness_ghost = 1;      
    }

  // get solver type
  ppAmr.query("velocity_solver_type", m_solverType);

  ppAmr.query("max_solver_iterations",m_maxSolverIterations);

  ppAmr.query("velocity_solver_tolerance", m_velocity_solver_tolerance);

  ppAmr.query("limit_velocity_rhs", m_limitVelRHS);
  ppAmr.query("limit_rhs_radius", m_gradLimitRadius);

  ppAmr.query("limit_flux_speed",m_limitFluxSpeed);

  ppAmr.query("do_initial_velocity_solve", m_doInitialVelSolve);
  ppAmr.query("do_initial_velocity_guess", m_doInitialVelGuess);
  ppAmr.query("initial_velocity_guess_type", m_initialGuessType);
  ppAmr.query("initial_velocity_guess_const_mu", m_initialGuessConstMu);
  ppAmr.query("initial_velocity_guess_solver_type", m_initialGuessSolverType);

  {
    Vector<Real> t(SpaceDim,0.0);
    ppAmr.queryarr("initial_velocity_guess_const_vel", t, 0, SpaceDim);
    m_initialGuessConstVel = RealVect(D_DECL(t[0], t[1], t[2]));
  }

  ppAmr.query("additional_velocity",m_additionalVelocity);


  //thickness diffusion options
  std::string diffusionTreatment = "none";
  ppAmr.query("diffusion_treatment", diffusionTreatment);
  if (diffusionTreatment == "implicit")
    {
      MayDay::Error("implicit diffusion invalid for now");
      m_diffusionTreatment = IMPLICIT;
    }
  else if (diffusionTreatment == "explicit")
    m_diffusionTreatment = EXPLICIT;
  ppAmr.query("additional_diffusivity",m_additionalDiffusivity);


  //option to advance thickness/temperature only on coarser levels
  ppAmr.query("finest_timestep_level",m_finest_timestep_level);

  ppAmr.query("basal_length_scale", m_basalLengthScale);
  ppAmr.query("basal_friction_decay",m_basalFrictionDecay);

  ppAmr.query("wallDrag",m_wallDrag);
  ppAmr.query("wallDragExtra",m_wallDragExtra);

  ppAmr.query("grounding_line_regularization",m_groundingLineRegularization);
  ppAmr.query("grounding_line_correction",m_groundingLineCorrection);

  //calving model options
  m_calvingModelPtr = CalvingModel::parseCalvingModel("CalvingModel");
  if (m_calvingModelPtr == NULL)
    {
      MayDay::Warning("trying to parse old style amr.calving_model_type");

      std::string calvingModelType = "NoCalvingModel";
      ppAmr.query("calving_model_type",calvingModelType);
      if (calvingModelType == "NoCalvingModel")
	{
	  m_calvingModelPtr = new NoCalvingModel;
	}
      else if (calvingModelType == "DeglaciationCalvingModelA")
	{
	  ParmParse ppc("DeglaciationCalvingModelA");
	  Real minThickness = 0.0;
	  ppc.get("min_thickness", minThickness );
	  Real calvingThickness = 0.0;
	  ppc.get("calving_thickness", calvingThickness );
	  Real calvingDepth = 0.0;
	  ppc.get("calving_depth", calvingDepth );
	  Real startTime = -1.2345678e+300;
	  ppc.query("start_time",  startTime);
	  Real endTime = 1.2345678e+300;
	  ppc.query("end_time",  endTime);
	  DeglaciationCalvingModelA* ptr = new DeglaciationCalvingModelA
	    (calvingThickness,  calvingDepth, minThickness, startTime, endTime);
	  m_calvingModelPtr = ptr;
	  
	}
      else if (calvingModelType == "DomainEdgeCalvingModel")
	{
	  ParmParse ppc("DomainEdgeCalvingModel");
	  Vector<int> frontLo(2,false); 
	  ppc.getarr("front_lo",frontLo,0,frontLo.size());
	  Vector<int> frontHi(2,false);
	  ppc.getarr("front_hi",frontHi,0,frontHi.size());
	  bool preserveSea = false;
	  ppc.query("preserveSea",preserveSea);
	  bool preserveLand = false;
	  ppc.query("preserveLand",preserveLand);
	  DomainEdgeCalvingModel* ptr = new DomainEdgeCalvingModel
	    (frontLo, frontHi,preserveSea,preserveLand);
	  m_calvingModelPtr = ptr;
	}
      else
	{
	  MayDay::Error("Unknown calving model");
	}
    }

  // now set up problem domains
  {
    IntVect loVect = IntVect::Zero;
    IntVect hiVect(D_DECL(ancells[0]-1, ancells[1]-1, ancells[2]-1));
#if BISICLES_Z == BISICLES_LAYERED
    {
      int nLayers = ancells[2];
      Vector<Real> faceSigma(nLayers+1);
      Real dsigma = 1.0 / Real(nLayers);
      for (unsigned int l = 0; l < faceSigma.size(); ++l)
	faceSigma[l] = dsigma * (Real(l));
      ppAmr.queryarr("sigma",faceSigma,0,faceSigma.size());
      setLayers(faceSigma);
    }
#endif
    ProblemDomain baseDomain(loVect, hiVect);
    // now set periodicity
    for (int dir=0; dir<SpaceDim; dir++) 
      baseDomain.setPeriodic(dir, is_periodic[dir]);

    // now set up vector of domains
    m_amrDomains.resize(m_max_level+1);
    m_amrDx.resize(m_max_level+1);

    m_amrDomains[0] = baseDomain;
    m_amrDx[0] = m_domainSize[0]/baseDomain.domainBox().size(0);

    for (int lev=1; lev<= m_max_level; lev++)
      {
        m_amrDomains[lev] = refine(m_amrDomains[lev-1],
                                   m_refinement_ratios[lev-1]);
        m_amrDx[lev] = m_amrDx[lev-1]/m_refinement_ratios[lev-1];
      }
  } // leaving problem domain setup scope
  
  std::string tagSubsetBoxesFile = "";
  m_vectTagSubset.resize(m_max_level);
  
  ppAmr.query("tagSubsetBoxesFile",tagSubsetBoxesFile);
  
  if (tagSubsetBoxesFile != "")
    {
      if (procID() == uniqueProc(SerialTask::compute))
  	{
	 
  	  ifstream is(tagSubsetBoxesFile.c_str(), ios::in);
  	  int lineno = 1;
  	  if (is.fail())
  	    {
  	      pout() << "Can't open " << tagSubsetBoxesFile << std::endl;
  	      MayDay::Error("Cannot open refine boxes file");
  	    }

      // keep this around until changes propagate through inputs...
#define BASE_SUBSET_ON_LEVEL true
#if BASE_SUBSET_ON_LEVEL    
          
  	  for (int lev = 0; lev < m_max_level; lev++)
  	    {
	     
  	      const char level[6] = "level";
  	      char s[6];
  	      is >> s;
  	      if (std::string(level) != std::string(s))
  		{
  		  pout() << "expected '" << level << "' at line " << lineno << ", got " << s << std::endl;
  		  MayDay::Error("bad input file");
  		}
  	      int inlev;
  	      is >> inlev;
  	      if (inlev != lev)
  		{
  		  pout() << "expected ' " << lev << "' at line " << lineno << std::endl;
  		  MayDay::Error("bad input file");
  		}
  	      //advance to next line
  	      while (is.get() != '\n');lineno++;
  	      int nboxes;
  	      is >> nboxes;
  	      if (nboxes > 0)
  		{
  		  for (int i = 0; i < nboxes; ++i)
  		    {
  		      Box box;
  		      is >> box;while (is.get() != '\n');lineno++;
  		      m_vectTagSubset[lev] |= box;
  		      pout() << " level " << lev << " refine box : " << box << std::endl;
  		    }
  		}
  	      //advance to next line
  	      while (is.get() != '\n');lineno++;
	     
	      if (lev > 0)
		{
		  //add lower level's subset to this subset
		  IntVectSet crseSet (m_vectTagSubset[lev-1]);
		  if (!crseSet.isEmpty())
		    {
		      crseSet.refine(m_refinement_ratios[lev-1]);
		      // crseSet.nestingRegion(m_block_factor,m_amrDomains[lev]);
		      if (m_vectTagSubset[lev].isEmpty())
			{
			  m_vectTagSubset[lev] = crseSet;
			} 
		      else
			{
			  m_vectTagSubset[lev] &= crseSet;
			} 
		    }
		 
		}
	      
  	    } // end loop over levels
#else //not BASE_SUBSET_ON_LEVEL    
  
  	  for (int lev = 0; lev < m_max_level; lev++)
  	    {
              // first set empty domains for all levels
              m_vectTagSubset[lev] = IntVectSet();
            }


  	      const char  levelChar[6] = "level";
              const char domainChar[7] = "Domain";
  	      char s[6];
  	      is >> s;
  	      if (std::string(level) != std::string(s))
  		{
  		  pout() << "expected '" << level << "' at line " << lineno << ", got " << s << std::endl;
  		  MayDay::Error("bad input file");
  		}
  	      int inlev;
  	      is >> inlev;
  	      if (inlev != lev)
  		{
  		  pout() << "expected ' " << lev << "' at line " << lineno << std::endl;
  		  MayDay::Error("bad input file");
  		}
  	      //advance to next line
  	      while (is.get() != '\n');lineno++;
  	      int nboxes;
  	      is >> nboxes;
  	      if (nboxes > 0)
  		{
  		  for (int i = 0; i < nboxes; ++i)
  		    {
  		      Box box;
  		      is >> box;while (is.get() != '\n');lineno++;
  		      m_vectTagSubset[lev] |= box;
  		      pout() << " level " << lev << " refine box : " << box << std::endl;
  		    }
  		}
  	      //advance to next line
  	      while (is.get() != '\n');lineno++;
	     
	      if (lev > 0)
		{
		  //add lower level's subset to this subset
		  IntVectSet crseSet (m_vectTagSubset[lev-1]);
		  if (!crseSet.isEmpty())
		    {
		      crseSet.refine(m_refinement_ratios[lev-1]);
		      // crseSet.nestingRegion(m_block_factor,m_amrDomains[lev]);
		      if (m_vectTagSubset[lev].isEmpty())
			{
			  m_vectTagSubset[lev] = crseSet;
			} 
		      else
			{
			  m_vectTagSubset[lev] &= crseSet;
			} 
		    }
		 
		}
	      
  	    } // end loop over levels        
#endif // not BASE_SUBSET_ON_LEVEL    
	} // end if serial compute
      for (int lev = 0; lev < m_max_level; lev++)
	broadcast(m_vectTagSubset[lev], uniqueProc(SerialTask::compute));
    }
  /// PatchGodunov used for thickness advection
  if (m_temporalAccuracy < 3)
    {
      // get PatchGodunov options -- first set reasonable defaults.
      // can over-ride from ParmParse
      int normalPredOrder = 2;
      bool useFourthOrderSlopes = true;
      bool usePrimLimiting = true;
      bool useCharLimiting = false;
      bool useFlattening = false;
      bool useArtificialViscosity = false;
      Real artificialViscosity = 0.0;
      
      // define AdvectionPhysics ptr
      // (does this need to be a special Ice-advection pointer due
      // to H factors?      
      
      m_thicknessPhysPtr = new AdvectPhysics;
      m_thicknessPhysPtr->setPhysIBC(m_thicknessIBCPtr);
    
      

      m_thicknessPatchGodVect.resize(m_max_level+1, NULL);

      for (int lev=0; lev<=m_max_level; lev++)
        {


          m_thicknessPatchGodVect[lev] = new PatchGodunov;
          m_thicknessPatchGodVect[lev]->define(m_amrDomains[lev],
                                               m_amrDx[lev],
                                               m_thicknessPhysPtr,
                                               normalPredOrder,
                                               useFourthOrderSlopes,
                                               usePrimLimiting,
                                               useCharLimiting,
                                               useFlattening,
                                               useArtificialViscosity,
                                               artificialViscosity);
        }
      

      //m_temperatureIBCPtr = new IceTemperatureIBC;
      m_temperaturePhysPtr = new AdvectPhysics;
      m_temperaturePhysPtr->setPhysIBC(m_temperatureIBCPtr);
      m_temperaturePatchGodVect.resize(m_max_level+1, NULL);
      for (int lev=0; lev<=m_max_level; lev++)
	{
	  m_temperaturePatchGodVect[lev] = new PatchGodunov;
	  
	  int normalPredOrder = 2;
	  bool useFourthOrderSlopes = true;
	  bool usePrimLimiting = true;
	  bool useCharLimiting = false;
	  bool useFlattening = false;
	  bool useArtificialViscosity = false;
	  Real artificialViscosity = 0.0;
	  m_temperaturePatchGodVect[lev]->define(m_amrDomains[lev],
						 m_amrDx[lev],
						 m_temperaturePhysPtr,
						 normalPredOrder,
						 useFourthOrderSlopes,
						 usePrimLimiting,
						 useCharLimiting,
						 useFlattening,
						 useArtificialViscosity,
						 artificialViscosity);
	}

    } // end if temporal accuracy < 3
  
 

  // check to see if we're using predefined grids
  bool usePredefinedGrids = false;
  std::string gridFile;
  if (ppAmr.contains("gridsFile"))
    {
      usePredefinedGrids = true;
      ppAmr.get("gridsFile",gridFile);
    }

  
  ParmParse geomPP("geometry");
  if (geomPP.contains("basalSlope") )
    {
      Vector<Real> basalSlope(SpaceDim, 0.0);
      geomPP.getarr("basalSlope", basalSlope, 0, SpaceDim);
      D_TERM(
             m_basalSlope[0] = basalSlope[0];,
             m_basalSlope[1] = basalSlope[1];,
             m_basalSlope[2] = basalSlope[2];);
    }


  // check to see if we're restarting from a checkpoint file
  if (!ppAmr.contains("restart_file"))
    {
      // if we're not restarting
      
      // now set up data holders
      m_old_thickness.resize(m_max_level+1, NULL);
      m_velocity.resize(m_max_level+1, NULL);
      m_faceVel.resize(m_max_level+1, NULL);
      m_diffusivity.resize(m_max_level+1);
      m_velBasalC.resize(m_max_level+1,NULL);
      m_cellMuCoef.resize(m_max_level+1,NULL);
      m_faceMuCoef.resize(m_max_level+1,NULL);
      m_velRHS.resize(m_max_level+1, NULL);
      m_surfaceThicknessSource.resize(m_max_level+1, NULL);
      m_basalThicknessSource.resize(m_max_level+1, NULL);
      m_balance.resize(m_max_level+1, NULL);
      m_temperature.resize(m_max_level+1, NULL);
#if BISICLES_Z == BISICLES_LAYERED
      m_layerXYFaceXYVel.resize(m_max_level+1, NULL);
      m_layerSFaceXYVel.resize(m_max_level+1, NULL);
      m_sTemperature.resize(m_max_level+1, NULL);
      m_bTemperature.resize(m_max_level+1, NULL);
#endif
      // allocate storage for m_old_thickness and m_velocity
      for (int lev=0; lev<m_velocity.size(); lev++)
        {
          m_old_thickness[lev] = new LevelData<FArrayBox>;
          m_velocity[lev] = new LevelData<FArrayBox>;
	  m_faceVel[lev] = new LevelData<FluxBox>;
	  m_velBasalC[lev] = new LevelData<FArrayBox>;
	  m_cellMuCoef[lev] = new LevelData<FArrayBox>;
	  m_faceMuCoef[lev] = new LevelData<FluxBox>;
	  m_velRHS[lev] = new LevelData<FArrayBox>;
	  m_diffusivity[lev] = 
	    RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>);
	  m_temperature[lev] = new LevelData<FArrayBox>;

#if BISICLES_Z == BISICLES_LAYERED
	  m_sTemperature[lev] = new LevelData<FArrayBox>;
	  m_bTemperature[lev] = new LevelData<FArrayBox>;
	  m_layerXYFaceXYVel[lev] = new LevelData<FluxBox>;
	  m_layerSFaceXYVel[lev] = new LevelData<FArrayBox>;
#endif

        }

      int finest_level = -1;
      if (usePredefinedGrids)
        {
          setupFixedGrids(gridFile);
        } 
      else
        {
          // now create  grids
          initGrids(finest_level);
        }
      
      // last thing to do is to set this to true from here on out...
      m_doInitialVelSolve = true;

      // that should be it
    }
  else
    {
      // we're restarting from a checkpoint file
      string restart_file;
      ppAmr.get("restart_file", restart_file);
      m_do_restart = true;
      restart(restart_file);
      
    }


  // set up counter of number of cells
  m_num_cells.resize(m_max_level+1, 0);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LayoutIterator lit = levelGrids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& thisBox = levelGrids.get(lit());
          m_num_cells[lev] += thisBox.numPts();
        }
    }


  // finally, set up covered_level flags
  m_covered_level.resize(m_max_level+1, 0);

  // note that finest level can't be covered.
  for (int lev=m_finest_level-1; lev>=0; lev--)
    {

      // if the next finer level is covered, then this one is too.
      if (m_covered_level[lev+1] == 1)
        {
          m_covered_level[lev] = 1;
        }
      else
        {
          // see if the grids finer than this level completely cover it
          IntVectSet fineUncovered(m_amrDomains[lev+1].domainBox());
          const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];

          LayoutIterator lit = fineGrids.layoutIterator();
          for (lit.begin(); lit.ok(); ++lit)
            {
              const Box& thisBox = fineGrids.get(lit());
              fineUncovered.minus_box(thisBox);
            }

          if (fineUncovered.isEmpty()) 
            {
              m_covered_level[lev] = 1;
            }
        }
    } // end loop over levels to determine covered levels

  m_initialSumIce = computeTotalIce();
  m_lastSumIce = m_initialSumIce;
  if (m_report_grounded_ice)
    {
      m_initialSumGroundedIce = computeTotalGroundedIce();
      m_lastSumGroundedIce = m_initialSumGroundedIce;
    }
}  
  
/// set BC for thickness advection
void
AmrIce::setThicknessBC( IceThicknessIBC* a_thicknessIBC)
{
  m_thicknessIBCPtr = a_thicknessIBC->new_thicknessIBC(); 
}

/// set BC for temperature advection
void 
AmrIce::setTemperatureBC( IceTemperatureIBC* a_temperatureIBC)
{
  m_temperatureIBCPtr = a_temperatureIBC->new_temperatureIBC();
}

void 
AmrIce::defineSolver()
{
  if (m_solverType == Picard)
    {
      // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }

      m_velSolver = new PicardSolver;
      
      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;

      int numLevels = m_finest_level +1;

      // make sure that the IBC has the correct grid hierarchy info
      m_thicknessIBCPtr->setGridHierarchy(m_vect_coordSys, m_amrDomains);

      m_velSolver->define(m_amrDomains[0],
                          m_constitutiveRelation,
			  m_basalFrictionRelation,
                          m_amrGrids,
                          m_refinement_ratios,
                          dxCrse,
                          m_thicknessIBCPtr,
                          numLevels);
      m_velSolver->setVerbosity(s_verbosity);

      m_velSolver->setTolerance(m_velocity_solver_tolerance);

      if (m_maxSolverIterations > 0)
        {
          m_velSolver->setMaxIterations(m_maxSolverIterations);
	  
        }
    }
  else  if (m_solverType == JFNK)
    {
      // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
     
      JFNKSolver* jfnkSolver;

      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;
      int numLevels = m_finest_level +1;
      
      // make sure that the IBC has the correct grid hierarchy info
      m_thicknessIBCPtr->setGridHierarchy(m_vect_coordSys, m_amrDomains);
      
 

      if (m_velSolver != NULL)
	{
	  // assume that any extant solver is also a JFNKSolver
	  jfnkSolver = static_cast<JFNKSolver*>(m_velSolver);	  
	}
      else {
	jfnkSolver = new JFNKSolver();
      }
      
      jfnkSolver->define(m_amrDomains[0],
			 m_constitutiveRelation,
			 m_basalFrictionRelation,
			 m_amrGrids,
			 m_refinement_ratios,
			 dxCrse,
			 m_thicknessIBCPtr,
			 numLevels);

      m_velSolver = jfnkSolver;

    }
  else  if (m_solverType == KnownVelocity)
    {
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }
      m_velSolver = new KnownVelocitySolver;
      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;
      int numLevels = m_finest_level +1;
      m_velSolver->define(m_amrDomains[0],
                          m_constitutiveRelation,
			  m_basalFrictionRelation,
                          m_amrGrids,
                          m_refinement_ratios,
                          dxCrse,
                          m_thicknessIBCPtr,
                          numLevels);
    }
#ifdef CH_USE_PETSC
  else if (m_solverType == PetscNLSolver)
    {
     // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }

      m_velSolver = new PetscIceSolver;
      
      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;

      int numLevels = m_finest_level +1;

      // make sure that the IBC has the correct grid hierarchy info
      m_thicknessIBCPtr->setGridHierarchy(m_vect_coordSys, m_amrDomains);

      m_velSolver->define(m_amrDomains[0],
                          m_constitutiveRelation,
			  m_basalFrictionRelation,
                          m_amrGrids,
                          m_refinement_ratios,
                          dxCrse,
                          m_thicknessIBCPtr,
                          numLevels);
      m_velSolver->setVerbosity(s_verbosity);

      m_velSolver->setTolerance(m_velocity_solver_tolerance);

      if (m_maxSolverIterations > 0)
        {
          m_velSolver->setMaxIterations(m_maxSolverIterations);	  
        }
    }
#endif
  else
    {
      MayDay::Error("unsupported velocity solver type");
    }
 
}

//inline 
//Real remainder(Real a, Real b)
//{
//  Real p = a/b; int i(p);
//  return std::min( p - i, p - 1 - i);
//}

void
AmrIce::run(Real a_max_time, int a_max_step)
{

  if (s_verbosity > 3) 
    {
      pout() << "AmrIce::run -- max_time= " << a_max_time 
             << ", max_step = " << a_max_step << endl;
    }


  Real dt;
  // only call computeInitialDt if we're not doing restart
  if (!m_do_restart)
    {
      dt = computeInitialDt();
    } 
  else
    {
      dt = computeDt();
    }

  
  // advance solution until done
  if ( !(m_plot_time_interval > TIME_EPS) || m_plot_time_interval > a_max_time) m_plot_time_interval = a_max_time;

  while ( a_max_time > m_time && (m_cur_step < a_max_step))
    {
      Real next_plot_time = std::min(m_plot_time_interval
				     *(1.0 + Real(int((m_time/m_plot_time_interval)))),
				     a_max_time);
      while ( (next_plot_time > m_time) && (m_cur_step < a_max_step)
	      && (dt > TIME_EPS))
	{
	  
	  // dump plotfile before regridding
	  if ( (m_cur_step%m_plot_interval == 0) && m_plot_interval > 0)
	    {
	      writePlotFile();
	    }
	  
	  if ((m_cur_step != 0) && (m_cur_step%m_regrid_interval ==0))
	    {
	      regrid();
	    }
	  
	  if (m_cur_step != 0)
	    {
	      // compute dt after regridding in case number of levels has changed
	      dt = computeDt();           
	    }
	  
	  if (next_plot_time - m_time + TIME_EPS < dt) 
	   dt =  std::max(2 * TIME_EPS, next_plot_time - m_time);
	  
	  
	  if ((m_cur_step%m_check_interval == 0) && (m_check_interval > 0)
	      && (m_cur_step != m_restart_step))
	    {
	      writeCheckpointFile();
	      if (m_cur_step > 0 && m_check_exit)
		{
		  if (s_verbosity > 2)
		    {
		      pout() << "AmrIce::exit on checkpoint" << endl;
		      return;
		    }
		}


	    }
	  
	  
	  timeStep(dt);
	  
	} // end of plot_time_interval
      if (m_plot_interval >= 0)
	writePlotFile();

    } // end timestepping loop

  // dump out final plotfile, if appropriate
  if (m_plot_interval >= 0)
    {
      writePlotFile();
    }
    
	  
  if (s_verbosity > 2)
    {
      pout() << "AmrIce::run finished" << endl;
    }
}


void
AmrIce::timeStep(Real a_dt)
{

  if (s_verbosity >=2) 
    {
      pout() << "Timestep " << m_cur_step 
             << " Advancing solution from time " 
             << m_time << " ( " << time() << ")" " with dt = " << a_dt << endl;
    }

  

  // first copy thickness into old thickness   
  for (int lev=0; lev <= m_finest_level ; lev++)
    {
      
      LevelData<FArrayBox>& oldThickness = *m_old_thickness[lev];
      LevelData<FArrayBox>& currentThickness = m_vect_coordSys[lev]->getH();

      // this way we avoid communication and maintain ghost cells...
      DataIterator dit = oldThickness.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          oldThickness[dit].copy(currentThickness[dit],0, 0, 1);
        }

    }
        

  // assumption here is that we've already computed the current velocity 
  // field, most likely at initialization or at the end of the last timestep...
  // so, we don't need to recompute the velocity at the start.

  if (m_temporalAccuracy < 3)
    {
      // use PatchGodunov hyperbolic solver

      // need a grown velocity field
      IntVect grownVelGhost(2*IntVect::Unit);
      Vector<LevelData<FArrayBox>* > grownVel(m_finest_level+1, NULL);
      Vector<LevelData<FluxBox>* > H_half(m_finest_level+1,NULL);
     

      for (int lev = finestTimestepLevel() ; lev>=0 ; lev--)
        {
	  
          const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	
          H_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, 
                                               IntVect::Unit);
          LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];

         
          
          // ensure that ghost cells for thickness  are filled in
          if (lev > 0)
            {          
              int nGhost = levelOldThickness.ghostVect()[0];
              PiecewiseLinearFillPatch thicknessFiller(levelGrids, 
                                                       m_amrGrids[lev-1],
                                                       1, 
                                                       m_amrDomains[lev-1],
                                                       m_refinement_ratios[lev-1],
                                                       nGhost);
              
              // since we're not subcycling, don't need to interpolate in time
              Real time_interp_coeff = 0.0;
              thicknessFiller.fillInterp(levelOldThickness,
                                         *m_old_thickness[lev-1],
                                         *m_old_thickness[lev-1],
                                         time_interp_coeff,
                                         0, 0, 1);

             
              
            }
          // these are probably unnecessary...
          levelOldThickness.exchange();
          	     
	 
          // do we need to also do a coarseAverage for the vel here?
        }
      
      
      for (int lev=0; lev<= finestTimestepLevel();  lev++)
        {
          
          // get AdvectPhysics object from PatchGodunov object
          PatchGodunov* patchGod = m_thicknessPatchGodVect[lev];
          AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>(patchGod->getGodunovPhysicsPtr());
          if (advectPhysPtr == NULL)
            {
              MayDay::Error("AmrIce::timestep -- unable to upcast GodunovPhysics to AdvectPhysics");
            }
          
          patchGod->setCurrentTime(m_time);
          
          // loop over grids on this level and compute H-Half
          const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
          LevelSigmaCS& levelCoordSys = *(m_vect_coordSys[lev]);
          //LevelData<FArrayBox>& levelVel = *m_velocity[lev];
          LevelData<FluxBox>& levelFaceVel = *m_faceVel[lev];
          LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];
          LevelData<FluxBox>& levelHhalf = *H_half[lev];

	  LevelData<FArrayBox>& levelSTS = *m_surfaceThicknessSource[lev];
	  LevelData<FArrayBox>& levelBTS = *m_basalThicknessSource[lev];
          CH_assert(m_surfaceFluxPtr != NULL);

          // set surface thickness source
          m_surfaceFluxPtr->surfaceThicknessFlux(levelSTS, *this, lev, a_dt);
	  
	  // set basal thickness source
	  m_basalFluxPtr->surfaceThicknessFlux(levelBTS, *this, lev, a_dt);

	  m_calvingModelPtr->modifySurfaceThicknessFlux(levelBTS, *this, lev, a_dt);


	  // modify face velocity, removing the diffusive component,
          // -D/H grad(H) and compute a source term div(D grad H))
    
	  // LevelData<FArrayBox> levelDiffusionSrc;
	  // if (m_diffusionTreatment == IMPLICIT)
	  //   { 
	  //     levelDiffusionSrc.define(levelGrids,1,IntVect::Zero);
	  //     LevelData<FluxBox> levelDiffusionVel(levelGrids,levelFaceVel.nComp(),IntVect::Unit);
	  //     computeDiffusionTerms(levelDiffusionVel,levelDiffusionSrc,lev);
	  //     levelDiffusionVel.exchange();

	  //     for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	  // 	{
	  // 	  for (int dir = 0; dir < SpaceDim; ++dir)
	  // 	    {
	  // 	      levelFaceVel[dit][dir] -= levelDiffusionVel[dit][dir];
	  // 	    }
	  // 	}
	  //   }
	  LevelData<FArrayBox> levelCCVel(levelGrids, SpaceDim, IntVect::Unit);
	  EdgeToCell( levelFaceVel, levelCCVel);
          
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              patchGod->setCurrentBox(levelGrids[dit]);
              advectPhysPtr->setVelocities(&(levelCCVel[dit]), 
                                           &(levelFaceVel[dit]));
              
	    
	      
	      FArrayBox advectiveSource(levelSTS[dit].box(),1);
	      advectiveSource.copy(levelSTS[dit]);
	      advectiveSource.plus(levelBTS[dit]);
	      // if (m_diffusionTreatment == IMPLICIT)
	      // 	{
	      // 	  advectiveSource.minus(levelDiffusionSrc[dit]);
	      // 	}
              patchGod->computeWHalf(levelHhalf[dit],
                                     levelOldThickness[dit],
                                     advectiveSource,
                                     a_dt,
                                     levelGrids[dit]);
	     

            } //end loop over grids

	 
	} // end loop over levels for computing Whalf
      
      // coarse average new H-Half to covered regions
      for (int lev= finestTimestepLevel(); lev>0; lev--)
        {
          CoarseAverageFace faceAverager(m_amrGrids[lev],
                                         1, m_refinement_ratios[lev-1]);
          faceAverager.averageToCoarse(*H_half[lev-1], *H_half[lev]);
        }
       
      Vector<LevelData<FluxBox>* > layerTH_half(m_finest_level+1,NULL);
      Vector<LevelData<FluxBox>* > layerH_half(m_finest_level+1,NULL);
      if (!m_isothermal)
	computeTHalf(layerTH_half, layerH_half, m_layerXYFaceXYVel, a_dt, m_time);
     
      // slc : having found H_half we can define temporary LevelSigmaCS at t + dt / 2
      // We want this for the metric terms in  temperature advection, and also when
      // m_temporalAccuracy == 2 to compute a new velocity field 
      Vector<RefCountedPtr<LevelSigmaCS> > vectCoords_half (m_finest_level+1);
     
      if (m_temporalAccuracy == 1)
	{
	  // just use the old time LevelSigmaCS
	  for (int lev=0; lev<= m_finest_level; lev++)
	    vectCoords_half[lev] = m_vect_coordSys[lev];
	}
      else if (m_temporalAccuracy == 2)
	{
	  for (int lev=0; lev<= finestTimestepLevel(); lev++)
	    {
	      IntVect sigmaCSGhost = IntVect::Unit;
	      RealVect dx = m_amrDx[lev]*RealVect::Unit;
	      vectCoords_half[lev] = RefCountedPtr<LevelSigmaCS> 
		(new LevelSigmaCS(m_amrGrids[lev], dx, sigmaCSGhost));
	      LevelSigmaCS& levelCoords_half = *vectCoords_half[lev];
	      LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
	      
	      ///todo : Here, assume that the base height doesn't change during the
	      ///timestep, which is not strictly true. Instead, we should perform 
	      ///an isostasy calculation at this point.
	      levelCoords_half.setBaseHeight(levelCoords.getBaseHeight());
	      levelCoords_half.setFaceSigma(levelCoords.getFaceSigma());
	      levelCoords_half.setIceDensity(levelCoords.iceDensity());
	      levelCoords_half.setGravity(levelCoords.gravity());
	      levelCoords_half.setWaterDensity(levelCoords.waterDensity());
	      //now set the thickness from H_half
	      LevelData<FluxBox>& levelH = *H_half[lev];
	      LevelData<FluxBox>& levelFaceH = levelCoords_half.getFaceH();
	      for (DataIterator dit( m_amrGrids[lev]); dit.ok(); ++dit)
		{
		  FluxBox& faceH = levelFaceH[dit];
		  faceH.copy(levelH[dit], levelH[dit].box());
		}
	      {
		LevelSigmaCS* crseCoords = (lev > 0)?&(*vectCoords_half[lev-1]):NULL;
		int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
		levelCoords_half.recomputeGeometryFace(crseCoords, refRatio);
	      }
	    }
	  // updateSurfaceGradient(vectCoords_half, m_amrGrids, 
	  // 			m_amrDomains, m_refinement_ratios, 
	  // 			m_amrDx, m_basalSlope, m_time, m_dt,
	  // 			(m_limitVelRHS)?m_gradLimitRadius:0,
	  // 			0, m_finest_level, m_thicknessIBCPtr);
	}

      // do velocity solve for half-time velocity field
      if (m_temporalAccuracy == 2)
        {

          // first, reset H in coordSys using H_half 
	  // (slc :: calculation was already done above and we will need the old time
          // also, so change solveVelcoityField so we can just swap LevelSigmaCSPointers)
	  MayDay::Error("m_temporalAccuracy ==  doesn't work yet");
          for (int lev=0; lev<= m_finest_level; lev++)
            {
              LevelData<FluxBox>& levelH = *H_half[lev];
              LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
              LevelData<FluxBox>& levelFaceH = levelCoords.getFaceH();
              DataIterator dit = levelH.dataIterator();
              for (dit.begin(); dit.ok(); ++dit)
                {
                  FluxBox& faceH = levelFaceH[dit];
                  faceH.copy(levelH[dit], levelH[dit].box());
                }
	      {
		LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
		int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
		levelCoords.recomputeGeometryFace(crseCoords, refRatio);
              }
                  
            } // end loop over levels
          
          // compute new ice velocity field
          solveVelocityField(m_velocity);
          
          // average cell-centered velocity field to faces just like before
          
        }
      
      // construct fluxes (could also store in face Velocity holder
      // instead of creating new data holder, but do it this way for now
      Vector<LevelData<FluxBox>* > vectFluxes(m_finest_level+1, NULL);
      
      for (int lev=0; lev<=finestTimestepLevel(); lev++)
        {
          LevelData<FluxBox>& levelFaceVel = *m_faceVel[lev];
          LevelData<FluxBox>& levelFaceH = *H_half[lev];
          // if we're doing AMR, we'll need to average these fluxes
          // down to coarser levels. As things stand now, 
          // CoarseAverageFace requires that the coarse LevelData<FluxBox>
          // have a ghost cell. 
          IntVect ghostVect = IntVect::Unit;
          vectFluxes[lev] = new LevelData<FluxBox>(m_amrGrids[lev],1, ghostVect);
          const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              FluxBox& faceVel = levelFaceVel[dit];
              FluxBox& faceH = levelFaceH[dit];
              FluxBox& flux = (*vectFluxes[lev])[dit];
	      
              const Box& gridBox = levelGrids[dit];
              
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  Box faceBox(gridBox);
                  faceBox.surroundingNodes(dir);
                  flux[dir].copy(faceH[dir], faceBox);
                  flux[dir].mult(faceVel[dir], faceBox, 0, 0, 1);
                }
#ifdef FO_UPWIND
	      const FArrayBox& H = m_vect_coordSys[lev]->getH()[dit];
              for (int dir=0; dir < SpaceDim; dir++)
                {
                  Box faceBox(gridBox);
                  faceBox.surroundingNodes(dir);
                  for (BoxIterator bit(faceBox); bit.ok(); ++bit)
                    {
                      const IntVect& iv = bit();
                      if (faceVel[dir](iv) > 0.0)
                        {
                          flux[dir](iv) = H(iv-BASISV(dir));
                        }
                      else
                        {
                          flux[dir](iv) = H(iv);
                        }
                    }
                  flux[dir].mult(faceVel[dir], faceBox, 0, 0, 1);
		  CH_assert( flux[dir].norm(0) < 1.0e+8);
                }
#endif      

            }
        } // end loop over levels

      if (m_diffusionTreatment == IMPLICIT)
      	{
      	  //remove diffusive flux, if any, from vectFluxes
      	  subtractDiffusiveFlux(vectFluxes);
      	}

      // average fine fluxes down to coarse levels
      for (int lev=finestTimestepLevel(); lev>0; lev--)
        {
          CoarseAverageFace faceAverager(m_amrGrids[lev],
                                         1, m_refinement_ratios[lev-1]);
          faceAverager.averageToCoarse(*vectFluxes[lev-1], *vectFluxes[lev]);
        }
      
      // make a copy of m_vect_coordSys before it is overwritten
     
      Vector<RefCountedPtr<LevelSigmaCS> > vectCoords_old (m_finest_level+1);
      for (int lev=0; lev<= m_finest_level; lev++)
	 {
	   IntVect sigmaCSGhost = IntVect::Unit;
	   RealVect dx = m_amrDx[lev]*RealVect::Unit;
	   vectCoords_old[lev] = RefCountedPtr<LevelSigmaCS> 
	     (new LevelSigmaCS(m_amrGrids[lev], dx, sigmaCSGhost));
	   LevelSigmaCS& levelCoords_old = *vectCoords_old[lev];
	   const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
	   
	   
	   levelCoords_old.setIceDensity(levelCoords.iceDensity());
	   levelCoords_old.setWaterDensity(levelCoords.waterDensity());
	   levelCoords_old.setGravity(levelCoords.gravity());
	   // todo replace the copies below with a deepCopy of levelCoords
	   for (DataIterator dit( m_amrGrids[lev]); dit.ok(); ++dit)
	      {
		FArrayBox& oldH = levelCoords_old.getH()[dit];
		const FArrayBox& H = levelCoords.getH()[dit];
		oldH.copy(H);
	      }
	   levelCoords_old.setBaseHeight(levelCoords.getBaseHeight());
	   {
	     LevelSigmaCS* crseCoords = (lev > 0)?&(*vectCoords_old[lev-1]):NULL;
	     int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
	     levelCoords_old.recomputeGeometry( crseCoords, refRatio);
	   }
	   //levelCoords_old.setSurfaceHeight(levelCoords.getSurfaceHeight());
	   //levelCoords_old.setGradSurface(levelCoords.getGradSurface());
	   //levelCoords_old.setGradSurfaceFace(levelCoords.getGradSurfaceFace());
	   
#if BISICLES_Z == BISICLES_LAYERED
	   levelCoords_old.setFaceSigma(levelCoords.getFaceSigma());
#endif
	 }
    
      // compute div(F) and update solution
      

      for (int lev=0; lev <= finestTimestepLevel() ; lev++)
        {
          DisjointBoxLayout& levelGrids = m_amrGrids[lev];
          LevelData<FluxBox>& levelFlux = *vectFluxes[lev];
          LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
          LevelData<FArrayBox>& levelNewH = levelCoords.getH();
	  LevelData<FArrayBox>& levelOldH = (*vectCoords_old[lev]).getH();
	  LevelData<FArrayBox>& levelBalance = *m_balance[lev];
          const RealVect& dx = levelCoords.dx();              

          DataIterator dit = levelGrids.dataIterator();          
          
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box& gridBox = levelGrids[dit];
              FArrayBox& newH = levelNewH[dit];
              FArrayBox& oldH = levelOldH[dit];//(*m_old_thickness[lev])[dit];
              FluxBox& thisFlux = levelFlux[dit];
              newH.setVal(0.0);
              
              // loop over directions and increment with div(F)
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  // use the divergence from 
                  // Chombo/example/fourthOrderMappedGrids/util/DivergenceF.ChF
                  FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                                  CHF_FRA(newH),
                                  CHF_BOX(gridBox),
                                  CHF_CONST_REAL(dx[dir]),
                                  CHF_INT(dir));

		  
                }
              // add in thickness source
	      // if there are still diffusive fluxes to deal
	      // with, the source term will be included then
	      if (m_diffusionTreatment != IMPLICIT)
		{
		  newH.minus((*m_surfaceThicknessSource[lev])[dit], gridBox,0,0,1);
		  newH.minus((*m_basalThicknessSource[lev])[dit], gridBox,0,0,1);
		}
		  
	      levelBalance[dit].copy(newH);


              if (m_evolve_thickness)
                {



		  if (m_floating_ice_stable)
		    {
		      //keep floating ice stable if required
		      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
		      for (BoxIterator bit(gridBox); bit.ok(); ++bit)
			{
			  const IntVect& iv = bit();
			  if (mask(iv) == FLOATINGMASKVAL)
			    {
			      newH(iv) = 0.0;
			    }
			}
		    }
		  
		  if (m_grounded_ice_stable)
		    {
		      //keep grounded ice stable if required
		      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
		      for (BoxIterator bit(gridBox); bit.ok(); ++bit)
			{
			  const IntVect& iv = bit();
			  if (mask(iv) == GROUNDEDMASKVAL)
			    {
			      newH(iv) = 0.0;
			    }
			}
		    }


                  newH *= -1*a_dt;
		  
                }
              else 
                {
                  newH.setVal(0.0);
                }
	      
	      newH.plus(oldH, 0, 0, 1);

            } // end loop over grids
        } // end loop over levels
      

      //include any diffusive fluxes
      if (m_evolve_thickness && m_diffusionTreatment == IMPLICIT)
	implicitThicknessCorrection(a_dt,  m_surfaceThicknessSource, m_basalThicknessSource);

      // average down thickness to coarser levels and fill in ghost cells
      // before calling recomputeGeometry. 
      int Hghost = 2;
      Vector<LevelData<FArrayBox>* > vectH(m_finest_level+1, NULL);
      for (int lev=0; lev<vectH.size(); lev++)
        {
          IntVect HghostVect = Hghost*IntVect::Unit;
          LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
          vectH[lev] = &levelCoords.getH();
        }

      //average from the finest level down
      for (int lev =  finestTimestepLevel() ; lev > 0 ; lev--)
	{
	  CoarseAverage averager(m_amrGrids[lev],
                                 1, m_refinement_ratios[lev-1]);
          averager.averageToCoarse(*vectH[lev-1],
                                   *vectH[lev]);
          
	}

      // //interpolate H to any levels above finestTimestepLevel()
      // for (int lev=finestTimestepLevel() + 1; lev < vectH.size(); lev++)
      // 	{
      // 	  FineInterp interpolator(m_amrGrids[lev],1,
      // 				  m_refinement_ratios[lev-1],
      // 				  m_amrDomains[lev]);
	  
      // 	  interpolator.interpToFine(*vectH[lev], *vectH[lev-1]);

      // 	}
      // now pass back over and do PiecewiseLinearFillPatch
      for (int lev=1; lev<vectH.size(); lev++)
        {
      
          PiecewiseLinearFillPatch filler(m_amrGrids[lev],
                                          m_amrGrids[lev-1],
                                          1, 
                                          m_amrDomains[lev-1],
                                          m_refinement_ratios[lev-1],
                                          Hghost);

          Real interp_coef = 0;
          filler.fillInterp(*vectH[lev],
                            *vectH[lev-1],
                            *vectH[lev-1],
                            interp_coef,
                            0, 0, 1);
        }
    
      
      //re-fill ghost cells ouside the domain
      for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
	{
	  RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
	  m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
					    m_amrDomains[lev],levelDx, m_time, m_dt);
	}
  
      //allow calving model to modify geometry and velocity
      for (int lev=0; lev<= m_finest_level; lev++)
	{
	  m_calvingModelPtr->endTimeStepModifyState
	    (m_vect_coordSys[lev]->getH(), *this, lev);
	}

   


      //dont allow thickness to be negative
      for (int lev=0; lev<= m_finest_level; lev++)
	{
	  DisjointBoxLayout& levelGrids = m_amrGrids[lev];
          LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
          LevelData<FArrayBox>& levelH = levelCoords.getH();
          DataIterator dit = levelGrids.dataIterator();          
          
          for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      Real lim = 0.0;
	      FORT_MAXFAB(CHF_FRA(levelH[dit]), 
			  CHF_CONST_REAL(lim), 
			  CHF_BOX(levelH[dit].box()));
	    }
	}

      //average from the finest level down
      for (int lev =  finestTimestepLevel() ; lev > 0 ; lev--)
	{
	  CoarseAverage averager(m_amrGrids[lev],
                                 1, m_refinement_ratios[lev-1]);
          averager.averageToCoarse(*vectH[lev-1],
                                   *vectH[lev]);
          
	}

      for (int lev=1; lev<vectH.size(); lev++)
        {
      
          PiecewiseLinearFillPatch filler(m_amrGrids[lev],
                                          m_amrGrids[lev-1],
                                          1, 
                                          m_amrDomains[lev-1],
                                          m_refinement_ratios[lev-1],
                                          Hghost);

          Real interp_coef = 0;
          filler.fillInterp(*vectH[lev],
                            *vectH[lev-1],
                            *vectH[lev-1],
                            interp_coef,
                            0, 0, 1);
        }

      //interpolate levelSigmaCS to any levels above finestTimestepLevel()
      for (int lev = finestTimestepLevel()+1 ; lev<= m_finest_level; lev++)
	{
	  m_vect_coordSys[lev]->interpFromCoarse(*m_vect_coordSys[lev-1],
						 m_refinement_ratios[lev-1],
						 false , true);
	}



      // recompute thickness-derived data in SigmaCS
      for (int lev=0; lev<= m_finest_level; lev++)
        {
	  LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
	  //slc : I think we need an exchange here
	  levelCoords.getH().exchange();
	  LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
	  int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
	  levelCoords.recomputeGeometry(crseCoords, refRatio);            
        }
     
      


      // updateSurfaceGradient(m_vect_coordSys, m_amrGrids, 
      // 			    m_amrDomains, m_refinement_ratios, 
      // 			    m_amrDx, m_basalSlope, m_time, m_dt,
      // 			    (m_limitVelRHS)?m_gradLimitRadius:0,
      // 			    0, m_finest_level, m_thicknessIBCPtr);

      if (!m_isothermal)
	updateTemperature(layerTH_half, layerH_half, m_layerXYFaceXYVel,
			  m_layerSFaceXYVel,  a_dt, m_time,
			  m_vect_coordSys, vectCoords_old, 
			  m_surfaceThicknessSource, m_basalThicknessSource);
      

      // clean up temp storage
      for (int lev=0; lev<=m_finest_level; lev++)
        {
          if (grownVel[lev] != NULL)
            {
              delete grownVel[lev];
              grownVel[lev] = NULL;
            }
          
          
          if (H_half[lev] != NULL)
            {
              delete H_half[lev];
              H_half[lev] = NULL;
            }

	  if (layerTH_half[lev] != NULL)
            {
              delete layerTH_half[lev];
              layerTH_half[lev] = NULL;
            }
	   if (layerH_half[lev] != NULL)
            {
              delete layerH_half[lev];
              layerH_half[lev] = NULL;
            }
	  


          if (vectFluxes[lev] != NULL)
            {
              delete vectFluxes[lev];
              vectFluxes[lev] = NULL;
            }

        }

    } // end if 1st or 2nd-order temporal accuracy
  else if (m_temporalAccuracy == 3)
    {
      MayDay::Error("AmrIce::timestep -- RK3 not implemented");
    }
  else if (m_temporalAccuracy == 4)
    {
      // storage for fluxes 
      Vector<LevelData<FluxBox>* > vectFluxes(m_finest_level+1, NULL);
      
      // temp storage for div(flux)
      Vector<LevelData<FArrayBox>* > divFlux(m_finest_level+1, NULL);

      // temp storage for new thickness
      Vector<LevelData<FArrayBox>* > tempThickness(m_finest_level+1, NULL);
      Vector<LevelData<FArrayBox>* > newThickness(m_finest_level+1, NULL);

      for (int lev=0; lev<=m_finest_level; lev++)
        {
          // if we're doing AMR, we'll need to average these fluxes
          // down to coarser levels. As things stand now, 
          // CoarseAverageFace requires that the coarse LevelData<FluxBox>
          // have a ghost cell. 
          IntVect ghostVect = IntVect::Unit;
          vectFluxes[lev] = new LevelData<FluxBox>(m_amrGrids[lev],1, 
                                                   ghostVect);
          divFlux[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1,
                                                  IntVect::Zero);
          
          tempThickness[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1,
                                                        ghostVect);


          newThickness[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1,
                                                       ghostVect);

          m_old_thickness[lev]->copyTo(*tempThickness[lev]);
          m_old_thickness[lev]->copyTo(*newThickness[lev]);
        }
      

      Real time = m_time;
      
      // RK Stage 1: compute RHS.at old time
      computeDivThicknessFlux(divFlux,vectFluxes, 
                              m_old_thickness, time, a_dt/6);

      // first part of update
      incrementWithDivFlux(newThickness, divFlux, a_dt/6);
      
      // predict half-time value using divFlux
      incrementWithDivFlux(tempThickness, divFlux, a_dt/2);

#if 0      
      // RHSTmp gets RHS(tOld)
      // flux registers are also incremented if appropriate
      a_op.evalRHS(RHSTmp,a_oldSoln,
                   a_time,a_dt/6);
      
      // RK Stage 1: compute update {phi}_1, partial update to new vals.
      
      // first part of update...
      a_op.updateODE(a_newSoln,RHSTmp,a_dt/6);
      
      
      // predict half-time value using LofPhiTmp
      a_op.updateODE(solnTmp, RHSTmp,a_dt/2);
      
#endif
      
      
      // RK Stage 2: compute RHS based on first predicted value
      updateCoordSysWithNewThickness(tempThickness);
      solveVelocityField(m_velocity);
      
      time = m_time + a_dt/2;
      
      computeDivThicknessFlux(divFlux,vectFluxes, 
                              tempThickness, time, a_dt/3);
      
      
      // RK Stage 2: compute update {phi}_2, partial update to new vals.
      incrementWithDivFlux(newThickness, divFlux, a_dt/3);  
      
      // now predict half-time value using latest estimate of RHS
      for (int lev=0; lev<= m_finest_level; lev++)
        {
          m_old_thickness[lev]->copyTo(*tempThickness[lev]);
        }
      incrementWithDivFlux(tempThickness, divFlux, a_dt/2);
      
      // RK Stage 3: compute new RHS.
      updateCoordSysWithNewThickness(tempThickness);
      
      solveVelocityField(m_velocity);

      computeDivThicknessFlux(divFlux,vectFluxes, 
                              tempThickness, time, a_dt/3);
      
      
      // RK Stage 3: compute update {phi}_3, partial update to new vals.
      incrementWithDivFlux(newThickness, divFlux, a_dt/3);  
      
      // predict value at new time
      for (int lev=0; lev<= m_finest_level; lev++)
        {
          m_old_thickness[lev]->copyTo(*tempThickness[lev]);
        }
      incrementWithDivFlux(tempThickness, divFlux, a_dt);
      
      
      // RK Stage 4: compute RHS.
      
      time = m_time+a_dt;
      updateCoordSysWithNewThickness(tempThickness);
      
      solveVelocityField(m_velocity);

      computeDivThicknessFlux(divFlux,vectFluxes, 
                              tempThickness, time, a_dt/6);
      
      // RK Stage 4: compute final update of solution.
      
      incrementWithDivFlux(newThickness, divFlux, a_dt/6);
      
      
      // now update coordSys
      updateCoordSysWithNewThickness(newThickness);

      // clean up temp storage
      for (int lev=0; lev<= m_finest_level; lev++)
        {
          if (vectFluxes[lev] != NULL)
            {
              delete vectFluxes[lev];
              vectFluxes[lev] = NULL;
            }


          if (divFlux[lev] != NULL)
            {
              delete divFlux[lev];
              divFlux[lev]= NULL;
            }


          if (tempThickness[lev] != NULL)
            {
              delete tempThickness[lev];
              tempThickness[lev] = NULL;
            }


          if (newThickness[lev] != NULL)
            {
              delete newThickness[lev];
              newThickness[lev] = NULL;
            }
        }
    }
  else
    {
      MayDay::Error("AmrIce::timestep -- un-defined temporal accuracy");
    }

  // compute new ice velocity field
  solveVelocityField(m_velocity);



#if 0  
  // now average to faces
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      CellToEdge(*m_velocity[lev], *faceVel[lev]);
      faceVel[lev]->exchange();      
    }
#endif

  // finally, update to new time and increment current step
  m_dt = a_dt;
  m_time += a_dt;
  m_cur_step += 1;
  
  // write diagnostic info, like sum of ice
  Real sumIce = computeTotalIce();
  Real diffSum = sumIce - m_lastSumIce;
  Real totalDiffSum = sumIce - m_initialSumIce;
  
  Real sumGroundedIce, diffSumGrounded, totalDiffGrounded;
  if (m_report_grounded_ice)
    {
      sumGroundedIce = computeTotalGroundedIce();
      diffSumGrounded = sumGroundedIce - m_lastSumGroundedIce;
      totalDiffGrounded = sumGroundedIce - m_initialSumGroundedIce;      
      m_lastSumGroundedIce = sumGroundedIce;
    }

  if (s_verbosity > 0) 
    {
      pout() << "Step " << m_cur_step << ", time = " << m_time << " ( " << time() << " ) " 
             << ": sum(ice) = " << sumIce 
             << " (" << diffSum
             << " " << totalDiffSum
             << ")" << endl;
      
      if (m_report_grounded_ice)
        {
          pout() << "Step " << m_cur_step << ", time = " << m_time << " ( " << time() << " ) "
                 << ": sum(grounded ice) = " << sumGroundedIce 
                 << " (" << diffSumGrounded
                 << " " << totalDiffGrounded
                 << ")" << endl;
        }      
    }
  

  m_lastSumIce = sumIce;

  if (s_verbosity > 0) 
    {
      pout () << "AmrIce::timestep " << m_cur_step
              << " --     end time = " 
        //<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
              << m_time  << " ( " << time() << " )"
        //<< " (" << m_time/secondsperyear << " yr)"
              << ", dt = " 
        //<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
              << a_dt
        //<< " ( " << a_dt/secondsperyear << " yr )"
              << endl;

    }

  
  int totalCellsAdvanced = 0;
  for (int lev=0; lev<m_num_cells.size(); lev++) 
    {
      totalCellsAdvanced += m_num_cells[lev];
    }
     
  if (s_verbosity > 0) 
    {
      pout() << "Time = " << m_time  
             << " cells advanced = " 
             << totalCellsAdvanced << endl;

      for (int lev=0; lev<m_num_cells.size(); lev++) 
        {
          pout () << "Time = " << m_time 
                  << "  level " << lev << " cells advanced = " 
                  << m_num_cells[lev] << endl;
        }
    }
  
}





// do regridding
void
AmrIce::regrid()
{

  CH_TIME("AmrIce::regrid");

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::regrid" << endl;
    }

  // only do any of this if the max level > 0
  if (m_max_level > 0) 
    {

      m_n_regrids++;

      // in this code, lbase is always 0
      int lbase =0;
      
      // first generate tags
      Vector<IntVectSet> tagVect(m_max_level);
      tagCells(tagVect);
      
      {
	// now generate new boxes
	int top_level = min(m_finest_level, m_max_level-1);
	Vector<Vector<Box> > old_grids(m_finest_level+1);
	Vector<Vector<Box> > new_grids;
	
	// this is clunky, but i don't know of a better way to turn 
	// a DisjointBoxLayout into a Vector<Box>
	for (int lev=0; lev<= m_finest_level; lev++) 
	    {
	      const DisjointBoxLayout& levelDBL = m_amrGrids[lev];
	      old_grids[lev].resize(levelDBL.size());
	      LayoutIterator lit = levelDBL.layoutIterator();
	      int boxIndex = 0;
	      for (lit.begin(); lit.ok(); ++lit, ++boxIndex) 
		{
		  old_grids[lev][boxIndex] = levelDBL[lit()];
		}
	    }
	
	int new_finest_level;
	
	BRMeshRefine meshrefine(m_amrDomains[0], m_refinement_ratios,
				m_fill_ratio, m_block_factor, 
				m_nesting_radius, m_max_box_size);
	
	new_finest_level = meshrefine.regrid(new_grids, tagVect, 
					     lbase, top_level, 
					     old_grids);

	//test to see if grids have changed
	bool gridsSame = true;
	for (int lev=lbase+1; lev<= new_finest_level; ++lev)
	  {
	    int numGridsNew = new_grids[lev].size();
	    Vector<int> procIDs(numGridsNew);
	    LoadBalance(procIDs, new_grids[lev]);
	    const DisjointBoxLayout newDBL(new_grids[lev], procIDs,
					     m_amrDomains[lev]);
	    const DisjointBoxLayout oldDBL = m_amrGrids[lev];
	    gridsSame &= oldDBL.sameBoxes(newDBL);
	  }
	if (gridsSame)
	  {
	     if (s_verbosity > 3) 
	       { 
		 pout() << "AmrIce::regrid -- grids unchanged" << endl;
	       }
	     //return;
	  }

	// now loop through levels and redefine if necessary
	for (int lev=lbase+1; lev<= new_finest_level; ++lev)
	    {
	      int numGridsNew = new_grids[lev].size();
	      Vector<int> procIDs(numGridsNew);
	      LoadBalance(procIDs, new_grids[lev]);
	      
	      const DisjointBoxLayout newDBL(new_grids[lev], procIDs,
					     m_amrDomains[lev]);
	      
	      const DisjointBoxLayout oldDBL = m_amrGrids[lev];
	      
	      m_amrGrids[lev] = newDBL;
	      
	      // build new storage
	      LevelData<FArrayBox>* old_oldDataPtr = m_old_thickness[lev];
	      LevelData<FArrayBox>* old_velDataPtr = m_velocity[lev];
	      LevelData<FArrayBox>* old_tempDataPtr = m_temperature[lev];

	      RefCountedPtr<LevelData<FluxBox> > old_diffDataPtr = m_diffusivity[lev];
	     
	      LevelData<FArrayBox>* new_oldDataPtr = 
		new LevelData<FArrayBox>(newDBL, 1, m_old_thickness[0]->ghostVect());
	      
	      LevelData<FArrayBox>* new_velDataPtr = 
		new LevelData<FArrayBox>(newDBL, 2, m_velocity[0]->ghostVect());

	      LevelData<FArrayBox>* new_tempDataPtr = 
		new LevelData<FArrayBox>(newDBL, m_temperature[0]->nComp(), 
					 m_temperature[0]->ghostVect());
	      //since the temperature data has changed
	      m_A_valid = false;

	      RefCountedPtr<LevelData<FluxBox> > new_diffDataPtr 
		(new LevelData<FluxBox>(newDBL, m_diffusivity[0]->nComp(), 
					m_diffusivity[0]->ghostVect()));
	      
#if BISICLES_Z == BISICLES_LAYERED
	      LevelData<FArrayBox>* old_sTempDataPtr = m_sTemperature[lev];
	      LevelData<FArrayBox>* old_bTempDataPtr = m_bTemperature[lev];
	      LevelData<FArrayBox>* new_sTempDataPtr = 
		new LevelData<FArrayBox>(newDBL, m_sTemperature[0]->nComp(),
					 m_sTemperature[0]->ghostVect());
	      LevelData<FArrayBox>* new_bTempDataPtr = 
		new LevelData<FArrayBox>(newDBL, m_bTemperature[0]->nComp(),
					 m_bTemperature[0]->ghostVect());
	     
#endif	      
	      

	      // also need to handle LevelSigmaCS 
	      {
		IntVect sigmaCSGhost = IntVect::Unit;
		RealVect dx = m_amrDx[lev]*RealVect::Unit;
		RefCountedPtr<LevelSigmaCS > oldCoordSys = m_vect_coordSys[lev];
		
		RefCountedPtr<LevelSigmaCS > auxCoordSys = (lev > 0)?m_vect_coordSys[lev-1]:oldCoordSys;

		m_vect_coordSys[lev] = RefCountedPtr<LevelSigmaCS >
		  (new LevelSigmaCS(newDBL, dx, sigmaCSGhost));
		m_vect_coordSys[lev]->setIceDensity(auxCoordSys->iceDensity());
		m_vect_coordSys[lev]->setWaterDensity(auxCoordSys->waterDensity());
		m_vect_coordSys[lev]->setGravity(auxCoordSys->gravity());
#if BISICLES_Z == BISICLES_LAYERED
		m_vect_coordSys[lev]->setFaceSigma(auxCoordSys->getFaceSigma());
#endif		
		LevelSigmaCS* crsePtr = &(*m_vect_coordSys[lev-1]);
		int refRatio = m_refinement_ratios[lev-1];

		// initialize geometry from thicknessIBCPtr
		if (!m_interpolate_zb)
		  {
		    
		    m_thicknessIBCPtr->regridIceGeometry(*m_vect_coordSys[lev],
							 dx, 
							 m_domainSize, 
							 m_time, 
							 crsePtr,
							 refRatio);   
		  }
		
		//interpolate thickness & (maybe) topography
		m_vect_coordSys[lev]->interpFromCoarse(*m_vect_coordSys[lev-1],
						       m_refinement_ratios[lev-1],
						       m_interpolate_zb , true);

		//Defer to m_thicknessIBCPtr for boundary values - 
                //interpolation won't cut the mustard because it only fills
                //ghost cells overlying the valid regions.
		RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
		m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
						    m_amrDomains[lev],levelDx, m_time, m_dt);
						    
						    	 
		
	
		LevelData<FArrayBox>& thisLevelH = m_vect_coordSys[lev]->getH();
		LevelData<FArrayBox>& thisLevelB = m_vect_coordSys[lev]->getTopography();
		
		// overwrite interpolated fields in valid regiopns with such valid old data as there is
		if (oldDBL.isClosed()){	  
		  const LevelData<FArrayBox>& oldLevelH = oldCoordSys->getH();
                  oldLevelH.copyTo(thisLevelH);
		  const LevelData<FArrayBox>& oldLevelB = oldCoordSys->getTopography();
                  oldLevelB.copyTo(thisLevelB);
                }

		// exchange is necessary to fill periodic ghost cells
                // which aren't filled by the copyTo from oldLevelH
                thisLevelH.exchange();
		m_vect_coordSys[lev]->exchangeTopography();

		//m_vect_coordSys[lev]->recomputeGeometry(crsePtr, refRatio);
#if 0
                DataIterator dit = newDBL.dataIterator();
		for (dit.begin(); dit.ok(); ++dit)
		  {
                    // first need to set physical-domain ghost cells for H
                    // (we need to do this now because we'll be using those 
                    // ghost cells to average from cells->faces)
                    
                    // as a kluge, do extrapolation here
                    // this is probably better moved to a new function
                    // in IceThicknessIBC                    
                    {
                      FArrayBox& thisH = thisLevelH[dit];
                      Box testBox = thisH.box();
                      testBox &= m_amrDomains[lev];
                      if (testBox != thisH.box())
                        {
                          const Box& domainBox = m_amrDomains[lev].domainBox();
                          // we've got non-periodic ghost cells to fill...
                          for (int dir=0; dir<SpaceDim; dir++)
                            {
                              if (!m_amrDomains[lev].isPeriodic(dir))
                                {
                                  int rad = 1;
                                  
                                  // lo-side -- for now, only need one row of 
                                  // ghost cells
                                  int hiLo = 0;
                                  Box ghostBoxLo = adjCellLo(domainBox, 
                                                           dir, rad);
                                  // do this to try to catch corner cells
                                  ghostBoxLo.grow(1);
                                  ghostBoxLo.grow(dir,-1);
                                  ghostBoxLo &= thisH.box();
                                  if (!ghostBoxLo.isEmpty())
                                    {
                                      FORT_SIMPLEEXTRAPBC(CHF_FRA(thisH),
                                                          CHF_BOX(ghostBoxLo),
                                                          CHF_INT(dir), 
                                                          CHF_INT(hiLo));
                                    }

                                  // hi-side
                                  hiLo = 1;
                                  Box ghostBoxHi = adjCellHi(domainBox, 
                                                           dir, rad);
                                  // do this to try to catch corner cells
                                  ghostBoxHi.grow(1);
                                  ghostBoxHi.grow(dir,-1);
                                  ghostBoxHi &= thisH.box();
                                  if(!ghostBoxHi.isEmpty())
                                    {
                                      FORT_SIMPLEEXTRAPBC(CHF_FRA(thisH),
                                                          CHF_BOX(ghostBoxHi),
                                                          CHF_INT(dir), 
                                                          CHF_INT(hiLo));
                                    }
                                                      
                                } // end if not periodic in this direction
                            } // end loop over directions
                        } // end if there are non-periodic domain ghost cells
                    } // end scope for bc setting
                  } // end loop over grids
#endif
		{
		  LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
		  int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
		  m_vect_coordSys[lev]->recomputeGeometry(crseCoords,refRatio);
		}
              }
		
	      // first fill with interpolated data from coarser level
	      
	      {
		// may eventually want to do post-regrid smoothing on this!
		FineInterp interpolator(newDBL, 1, 
					m_refinement_ratios[lev-1],
					m_amrDomains[lev]);
	    
		interpolator.interpToFine(*new_oldDataPtr, *m_old_thickness[lev-1]);
	
		// now copy old-grid data into new holder
		if (old_oldDataPtr != NULL) 
		  {
		    if ( oldDBL.isClosed())
		      {
			old_oldDataPtr->copyTo(*new_oldDataPtr);
		      }
		    // can now delete old data 
		    delete old_oldDataPtr;
		  }
		
	      }
	      
	      {
		// may eventually want to do post-regrid smoothing on this!
		FineInterp interpolator(newDBL, 2, 
					m_refinement_ratios[lev-1],
					m_amrDomains[lev]);
		
		interpolator.interpToFine(*new_velDataPtr, *m_velocity[lev-1]);
		
		// now copy old-grid data into new holder
		if (old_velDataPtr != NULL)
		  {
		    if (oldDBL.isClosed()) 
		      {
			old_velDataPtr->copyTo(*new_velDataPtr);
		      }
		    // can now delete old data 
		    delete old_velDataPtr;
		  }
		
	      }

	      {
		// may eventually want to do post-regrid smoothing on this
		FineInterp interpolator(newDBL,m_temperature[0]->nComp(),
					m_refinement_ratios[lev-1],
					m_amrDomains[lev]);
		interpolator.interpToFine(*new_tempDataPtr, *m_temperature[lev-1]);

		
		PiecewiseLinearFillPatch ghostFiller
		  (m_amrGrids[lev],
		   m_amrGrids[lev-1],
		   m_temperature[lev-1]->nComp(),
		   m_amrDomains[lev-1],
		   m_refinement_ratios[lev-1],
		   m_temperature[lev-1]->ghostVect()[0]);

		ghostFiller.fillInterp(*new_tempDataPtr,*m_temperature[lev-1],
				      *m_temperature[lev-1],1.0,0,0,
				      m_temperature[lev-1]->nComp());

		if (old_tempDataPtr != NULL && oldDBL.isClosed())
		  {
		    old_tempDataPtr->copyTo(*new_tempDataPtr);
		  }
		delete old_tempDataPtr;
		new_tempDataPtr->exchange();
	      }
	      
#if BISICLES_Z == BISICLES_LAYERED
	      {
		// may eventually want to do post-regrid smoothing on this
		FineInterp interpolator(newDBL,m_sTemperature[0]->nComp(),
					m_refinement_ratios[lev-1],
					m_amrDomains[lev]);

		PiecewiseLinearFillPatch ghostFiller
		  (m_amrGrids[lev],
		   m_amrGrids[lev-1],
		   m_sTemperature[lev-1]->nComp(),
		   m_amrDomains[lev-1],
		   m_refinement_ratios[lev-1],
		   m_sTemperature[lev-1]->ghostVect()[0]);

	

		interpolator.interpToFine(*new_sTempDataPtr, *m_sTemperature[lev-1]);
		
		ghostFiller.fillInterp(*new_sTempDataPtr,*m_sTemperature[lev-1],
				       *m_sTemperature[lev-1],1.0,0,0,
				       m_sTemperature[lev-1]->nComp());


		if (old_sTempDataPtr != NULL && oldDBL.isClosed())
		  {
		    old_sTempDataPtr->copyTo(*new_sTempDataPtr);
		  }
		delete old_sTempDataPtr;
		new_sTempDataPtr->exchange();
		interpolator.interpToFine(*new_bTempDataPtr, *m_bTemperature[lev-1]);

		ghostFiller.fillInterp(*new_bTempDataPtr,*m_bTemperature[lev-1],
				       *m_bTemperature[lev-1],1.0,0,0,
				       m_bTemperature[lev-1]->nComp());
		
		if (old_bTempDataPtr != NULL && oldDBL.isClosed())
		  {
		    old_bTempDataPtr->copyTo(*new_bTempDataPtr);
		  }
		delete old_bTempDataPtr;
		new_bTempDataPtr->exchange();
	      }
#endif

	      // now copy old-grid data into new holder
	      if (old_diffDataPtr != NULL && oldDBL.isClosed())
		{
		  
		  (*old_diffDataPtr).copyTo(*new_diffDataPtr);
		  //old data will be destroyed when  old_diffDataPtr goes out of scope
		  
		}
	      
	      // now copy new holders into multilevel arrays
	      m_old_thickness[lev] = new_oldDataPtr;
	      m_velocity[lev] = new_velDataPtr;
	      m_diffusivity[lev] = new_diffDataPtr;
	      m_temperature[lev] = new_tempDataPtr;
#if BISICLES_Z == BISICLES_LAYERED
	      m_sTemperature[lev] = new_sTempDataPtr;
	      m_bTemperature[lev] = new_bTempDataPtr;
#endif


	      if (m_velBasalC[lev] != NULL)
		{
		  delete m_velBasalC[lev];
		}
	      m_velBasalC[lev] = new LevelData<FArrayBox>(newDBL, 1, IntVect::Unit);
	      
	      if (m_cellMuCoef[lev] != NULL)
		{
		  delete m_cellMuCoef[lev];
		}
	      m_cellMuCoef[lev] = new LevelData<FArrayBox>(newDBL, 1, IntVect::Unit);
  
	      if (m_faceMuCoef[lev] != NULL)
		{
		  delete m_faceMuCoef[lev];
		}
	      m_faceMuCoef[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);
	     
	      if (m_velRHS[lev] != NULL)
		{
		  delete m_velRHS[lev];
		}
	      m_velRHS[lev] = new LevelData<FArrayBox>(newDBL, 2, 
						       IntVect::Zero);

	      if (m_faceVel[lev] != NULL)
		{
		  delete m_faceVel[lev];
		}
	      m_faceVel[lev] = new LevelData<FluxBox>
		(newDBL, 1, IntVect::Unit);

	      if (m_surfaceThicknessSource[lev] != NULL)
		{
		  delete m_surfaceThicknessSource[lev];
		}
	      m_surfaceThicknessSource[lev] = 
		new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit) ;
	      
	      if (m_basalThicknessSource[lev] != NULL)
		{
		  delete m_basalThicknessSource[lev];
		}
	      m_basalThicknessSource[lev] = 
		new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit) ;
	      
	      if (m_balance[lev] != NULL)
		{
		  delete m_balance[lev];
		}
	      m_balance[lev] = 
		new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit) ;


#if BISICLES_Z == BISICLES_LAYERED
	      if (m_layerXYFaceXYVel[lev] != NULL)
		{
		  delete m_layerXYFaceXYVel[lev];
		}

	      m_layerXYFaceXYVel[lev] = new LevelData<FluxBox>
		(newDBL, m_nLayers, IntVect::Unit);
	      
	      if (m_layerSFaceXYVel[lev] != NULL)
		{
		  delete m_layerSFaceXYVel[lev];
		}
	      
	      m_layerSFaceXYVel[lev] = new LevelData<FArrayBox>
		(newDBL, SpaceDim*(m_nLayers + 1), IntVect::Unit);
#endif		
	      
	    } // end loop over currently defined levels
	  
	  // now ensure that any remaining levels are null pointers
	  // (in case of de-refinement)
	  for (int lev=new_finest_level+1; lev<m_old_thickness.size(); lev++)
	    {
	      if (m_old_thickness[lev] != NULL) 
		{
		  delete m_old_thickness[lev];
		  m_old_thickness[lev] = NULL;
		}


	      if (m_velocity[lev] != NULL) 
		{
		  delete m_velocity[lev];
		  m_velocity[lev] = NULL;
		}
	      
	      if (m_temperature[lev] != NULL) 
		{
		  delete m_temperature[lev];
		  m_temperature[lev] = NULL;
		}
#if BISICLES_Z == BISICLES_LAYERED
	      if (m_sTemperature[lev] != NULL) 
		{
		  delete m_sTemperature[lev];
		  m_sTemperature[lev] = NULL;
		}
	      if (m_bTemperature[lev] != NULL) 
		{
		  delete m_bTemperature[lev];
		  m_bTemperature[lev] = NULL;
		}
#endif	      	      
  
              if (m_velRHS[lev] != NULL)
                {
                  delete m_velRHS[lev];
                  m_velRHS[lev] = NULL;
                }
	      
	      if (m_velBasalC[lev] != NULL)
                {
                  delete m_velBasalC[lev];
                  m_velBasalC[lev] = NULL;
                }

	  
	      DisjointBoxLayout emptyDBL;
	      m_amrGrids[lev] = emptyDBL;
	    }
      
	  m_finest_level = new_finest_level;



	  // set up counter of number of cells
	  for (int lev=0; lev<=m_max_level; lev++)
	    {
	      m_num_cells[lev] = 0;
	      if (lev <= m_finest_level) 
		{
		  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
		  LayoutIterator lit = levelGrids.layoutIterator();
		  for (lit.begin(); lit.ok(); ++lit)
		    {
		      const Box& thisBox = levelGrids.get(lit());
		      m_num_cells[lev] += thisBox.numPts();
		    }
		} 
	    }
      
      
	  // finally, set up covered_level flags
	  m_covered_level.resize(m_max_level+1, 0);
	  // note that finest level can't be covered.
	  for (int lev=m_finest_level-1; lev>=0; lev--)
	    {
          
	      // if the next finer level is covered, then this one is too.
	      if (m_covered_level[lev+1] == 1)
		{
		  m_covered_level[lev] = 1;
		}
	      else
		{
		  // see if the grids finer than this level completely cover it
		  IntVectSet fineUncovered(m_amrDomains[lev+1].domainBox());
		  const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];
              
		  LayoutIterator lit = fineGrids.layoutIterator();
		  for (lit.begin(); lit.ok(); ++lit)
		    {
		      const Box& thisBox = fineGrids.get(lit());
		      fineUncovered.minus_box(thisBox);
		    }
              
		  if (fineUncovered.isEmpty()) 
		    {
		      m_covered_level[lev] = 1;
		    }
		}
	    } // end loop over levels to determine covered levels

	  //velocity solver needs to be re-defined
	  defineSolver();
	  //solve velocity field, but use the previous initial residual norm in place of this one
	  solveVelocityField(m_velocity, m_velocitySolveInitialResidualNorm);
	  
	} // end if tags changed
    } // end if max level > 0 in the first place
  
  m_groundingLineProximity_valid = false;
  m_viscousTensor_valid = false;
} 
      
                              
void 
AmrIce::tagCells(Vector<IntVectSet>& a_tags)
{
  
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::tagCells" << endl;
    }

  
  int top_level = a_tags.size();
  top_level = min(m_tag_cap,min(top_level-1, m_finest_level));
  // loop over levels
  for (int lev=0; lev<=top_level; lev++)
    {
      IntVectSet& levelTags = a_tags[lev];
      tagCellsLevel(levelTags, lev);
      IntVectSet& tagSubset = m_vectTagSubset[lev];
      if ( tagSubset.numPts() > 0)
	{
	  levelTags &= tagSubset;
	}
    }

  //throw away any coarse level tags outside m_tag_subset
  // if (s_verbosity > 3) 
  //   { 
  //     pout() << "AmrIce::tagCells, subset II" << endl;
  //   }
  // if (m_tag_subset.numPts() > 0)
  //   {
  //     IntVectSet tag_subset = m_tag_subset;
  //     a_tags[0] &= tag_subset;
  //     for (int lev = 1; lev <= top_level; lev++)
  // 	{
  // 	  tag_subset.refine(m_refinement_ratios[lev-1]);
  // 	  a_tags[lev] &= tag_subset;
  // 	}

  //   }

}

void
AmrIce::tagCellsLevel(IntVectSet& a_tags, int a_level)
{

  if (s_verbosity > 4) 
    { 
      pout() << "AmrIce::tagCellsLevel " << a_level << endl;
    }


  // base tags on undivided gradient of velocity
  // first stab -- don't do BC's; just do one-sided
  // stencils at box edges (hopefully good enough), 
  // since doing BC's properly is somewhat expensive.

  DataIterator dit = m_velocity[a_level]->dataIterator();
  
  LevelData<FArrayBox>& levelVel = *m_velocity[a_level];

  const DisjointBoxLayout& levelGrids = m_amrGrids[a_level];

  const LevelSigmaCS& levelCS = *m_vect_coordSys[a_level];

  // need to ensure that ghost cells are set properly
  levelVel.exchange(levelVel.interval());

  const LevelData<FluxBox>& levelFaceH = levelCS.getFaceH();

  LevelData<FArrayBox>& levelC = *m_velBasalC[a_level];

  IntVectSet local_tags;
  if (m_tagOnGradVel)
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          // note that we only need one component here
          // because the fortran subroutine stores the max(abs(grad)) 
          // over all components into the 0th position
          FArrayBox gradVel(levelGrids[dit()], 1);
          
          for (int dir=0; dir<SpaceDim; dir++)
            {
              const Box b = levelGrids[dit()];
              const Box bcenter = b & grow ( m_amrDomains[a_level], 
                                             -BASISV(dir) );
              const Box blo = b & adjCellLo( bcenter, dir );
              const Box bhi = b & adjCellHi( bcenter, dir );
              const int haslo = ! blo.isEmpty();
              const int hashi = ! bhi.isEmpty();
              FORT_UNDIVIDEDGRAD ( CHF_FRA1(gradVel,0),
                                   CHF_CONST_FRA(levelVel[dit()]),
                                   CHF_BOX(bcenter),
                                   CHF_BOX(blo),
                                   CHF_BOX(bhi),
                                   CHF_CONST_INT(dir),
                                   CHF_CONST_INT(haslo),
                                   CHF_CONST_INT(hashi));
              
              
              // now tag cells based on values
              BoxIterator bit(levelGrids[dit()]);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  if (abs(gradVel(iv,0)) > m_tagging_val) 
                    local_tags |= iv;
                } // end loop over cells
            } // end loop over directions
        } // end loop over grids
    } // end if tag on grad vel


  // tag on laplacian(velocity)     
  if (m_tagOnLapVel | m_tagOnGroundedLapVel)
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox lapVel(levelGrids[dit()], 2);
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
          lapVel.setVal(0.0);
          Real alpha = 0;
          Real beta = 1.0;
              
          // use undivided laplacian (set dx = 1)
          Real bogusDx = 1.0;
	  Box lapBox = levelVel[dit].box();
	  lapBox.grow(-2);
          lapBox &=  levelGrids[dit];
          // assumes that ghost cells boundary conditions are properly set
          FORT_OPERATORLAP(CHF_FRA(lapVel),
                           CHF_FRA(levelVel[dit]),
                           CHF_BOX(lapBox),
                           CHF_CONST_REAL(bogusDx),
                           CHF_CONST_REAL(alpha),
                           CHF_CONST_REAL(beta));
                            
          // now tag cells based on values
          BoxIterator bit(lapBox);
	  
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
	      for (int comp=0; comp<lapVel.nComp(); comp++)
		{
		  if ( (m_tagOnGroundedLapVel && mask(iv) == GROUNDEDMASKVAL) | m_tagOnLapVel )
		    {
		      if ( (abs(lapVel(iv,comp)) > m_laplacian_tagging_val) 
			   &&  (levelC[dit](iv) < m_laplacian_tagging_max_basal_friction_coef)) 
			local_tags |= iv;
		    }
		}
	      
            } // end loop over cells
        } // end loop over grids
    } // end if tag on laplacian(vel)
    

  // sometimes, it is easier to note where the grounding line is
  // and refine to the maximum level there
  if (m_tagGroundingLine)
    {
     
      for (dit.begin(); dit.ok(); ++dit)
	{
	  Box sbox = levelGrids[dit()];
	  //sbox.grow(-1);
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	  for (BoxIterator bit(sbox) ; bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		 int  tdir = (dir + 1)%SpaceDim;
		 const IntVect& ivm = iv - BASISV(dir);
		 const IntVect& ivp = iv + BASISV(dir);
		 
		 if (mask(iv) == GROUNDEDMASKVAL &&  levelC[dit](iv) <  m_groundingLineTaggingMaxBasalFrictionCoef)
		   {
		     if (mask(ivm) == FLOATINGMASKVAL || mask(ivm) == OPENSEAMASKVAL )
		       {
			 if (std::abs(levelVel[dit](iv,dir)) > m_groundingLineTaggingMinVel
			     || std::abs(levelVel[dit](ivm,dir)) > m_groundingLineTaggingMinVel
			     || std::abs(levelVel[dit](iv,tdir)) > m_groundingLineTaggingMinVel
			     || std::abs(levelVel[dit](ivm,tdir)) > m_groundingLineTaggingMinVel)
			   {
			     local_tags |= iv;  
			     local_tags |= ivm;
			   }   
		       }
			 
		     if ( mask(ivp) == FLOATINGMASKVAL || mask(ivp) == OPENSEAMASKVAL)
		       {
			 if (std::abs(levelVel[dit](iv,dir)) > m_groundingLineTaggingMinVel
			     || std::abs(levelVel[dit](ivp,dir)) > m_groundingLineTaggingMinVel
			     || std::abs(levelVel[dit](iv,tdir)) > m_groundingLineTaggingMinVel
			     || std::abs(levelVel[dit](ivp,tdir)) > m_groundingLineTaggingMinVel)
			   {
			     local_tags |= iv;  
			     local_tags |= ivp;
			   } 
		       }
		   }
		}
	    }
	}
    }
  
  // tag on div(H grad (vel)) 
  if (m_tagOndivHgradVel)
    {
      for (dit.begin(); dit.ok(); ++dit)
        {      
	  Box box = levelGrids[dit];
	  box.grow(-1);
	  const FluxBox& faceH = levelFaceH[dit];
	  const FArrayBox& vel = levelVel[dit];
	  BoxIterator bit(levelGrids[dit()]);
	  for (bit.begin(); bit.ok(); ++bit)
	  {
	    const IntVect& iv = bit();
	    for (int comp=0; comp < vel.nComp() ; comp++)
	      {
		Real t = 0.0;
		for (int dir=0; dir < SpaceDim ; dir++)
		  {		  
		    IntVect ivp = iv + BASISV(dir);
		    IntVect ivm = iv - BASISV(dir);
		  
		    t += faceH[dir](iv) * (vel(iv,comp)-vel(ivm,comp))
		      - faceH[dir](ivp) * (vel(ivp,comp)-vel(iv,comp));
		      
		  }
	
		if (abs(t) > m_divHGradVel_tagVal)
		  {
		   
		    local_tags |= iv;
		  }
	      }
	  }// end loop over cells
        } // end loop over grids
    } // end if tag on div(H grad (vel)) 
  
  if (m_tagOnEpsSqr)
    {
      IntVect tagsGhost = IntVect::Zero;
      LevelData<FArrayBox> epsSqr(levelGrids, 1, tagsGhost);
      LevelData<FArrayBox>* crseVelPtr = NULL;
      int nRefCrse = -1;
      if (a_level > 0)
        {
          crseVelPtr = m_velocity[a_level-1];
          nRefCrse = m_refinement_ratios[a_level-1];
        }

      m_constitutiveRelation->computeStrainRateInvariant(epsSqr,
                                                         levelVel,
                                                         crseVelPtr, 
                                                         nRefCrse,
                                                         levelCS,
                                                         tagsGhost);

      
      for (dit.begin(); dit.ok(); ++dit)
        {                    
          // now tag cells based on values
          // want undivided gradient
          epsSqr[dit] *= m_amrDx[a_level] * m_amrDx[a_level] ;
          Real levelTagVal = m_epsSqr_tagVal;
          BoxIterator bit(levelGrids[dit()]);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (abs(epsSqr[dit](iv,0)) > levelTagVal) 
                local_tags |= iv;
            } // end loop over cells
        } // end loop over grids
    } // end if tagging on strain rate invariant


  if (m_tagOnVelRHS)
    {
      for (dit.begin(); dit.ok(); ++dit)
        {                          
          const FArrayBox& thisVelRHS = (*m_velRHS[a_level])[dit];

          // now tag cells based on values
          // want RHS*dx (undivided gradient)
          Real levelTagVal = m_velRHS_tagVal/m_amrDx[a_level];
          BoxIterator bit(levelGrids[dit()]);
          for (int comp=0; comp<thisVelRHS.nComp(); comp++)
            {
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  if (abs(thisVelRHS(iv,comp)) > levelTagVal) 
                    local_tags |= iv;
                } // end loop over cells
            } // end loop over components
        } // end loop over grids
    } // end if tagging on velRHS
  
  // tag cells where thickness goes to zero
  if (m_tagMargin)
    {
      const LevelData<FArrayBox>& levelH = levelCS.getH();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box gridBox = levelGrids[dit];
          const FArrayBox& H = levelH[dit];

          for (BoxIterator bit(gridBox); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              // neglect diagonals for now...
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  IntVect ivm = iv - BASISV(dir);
                  IntVect ivp = iv + BASISV(dir);
                  if ( (H(iv,0) > 0) && (H(ivm,0) < TINY_THICKNESS) )
                    {
                      local_tags |= iv;
                      local_tags |= ivm;
                    } // end if low-side margin
                  if ( (H(iv,0) > 0) && (H(ivp,0) < TINY_THICKNESS) )
                    {
                      local_tags |= iv;
                      local_tags |= ivp;
                    } // end high-side margin
                } // end loop over directions
            } // end loop over cells
        } // end loop over boxes
    } // end if tagging on ice margins

  // tag anywhere there's ice
  if (m_tagAllIce)
    {
      const LevelData<FArrayBox>& levelH = levelCS.getH();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box gridBox = levelGrids[dit];
          const FArrayBox& H = levelH[dit];
          BoxIterator bit(gridBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (H(iv,0) > 0.0)
                {
                  local_tags |= iv;
                }
            } // end bit loop
        } // end loop over boxes
    } // end if tag all ice

  // now buffer tags
  local_tags.grow(m_tags_grow); 
  local_tags &= m_amrDomains[a_level];

  a_tags = local_tags;

}

void
AmrIce::tagCellsInit(Vector<IntVectSet>& a_tags)
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::tagCellsInit" << endl;
    }


  tagCells(a_tags);
  m_vectTags = a_tags;
  
}


void
AmrIce::initGrids(int a_finest_level)
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::initGrids" << endl;
    }


  m_finest_level = 0;
  // first create base level
  Vector<Box> baseBoxes;
  domainSplit(m_amrDomains[0], baseBoxes, m_max_base_grid_size, 
              m_block_factor);

  Vector<int> procAssign(baseBoxes.size());
  LoadBalance(procAssign,baseBoxes);
  
  DisjointBoxLayout baseGrids(baseBoxes, procAssign, m_amrDomains[0]);

  if (s_verbosity > 2) 
    {
      long long numCells0 = baseGrids.numCells();
      pout() << "Level 0: " << numCells0 << " cells: " << baseGrids << endl;
    }

  m_amrGrids.resize(m_max_level+1);
  m_amrGrids[0] = baseGrids;

  levelSetup(0,baseGrids);

  LevelData<FArrayBox>& baseLevelVel = *m_velocity[0];
  DataIterator baseDit = baseGrids.dataIterator();
  for (baseDit.begin(); baseDit.ok(); ++baseDit)
    {
      // initial guess at base-level velocity is zero
      baseLevelVel[baseDit].setVal(0.0);
    }

  // define solver before calling initData
  defineSolver();

  // initialize base level data
  initData(m_vect_coordSys,
           m_velocity);

  int numLevels = 1;
  bool moreLevels = (m_max_level > 0);

  int baseLevel = 0;
  //int topLevel = m_finest_level;
  
  
  BRMeshRefine meshrefine;
  if (moreLevels)
    {
      meshrefine.define(m_amrDomains[0], m_refinement_ratios,
                        m_fill_ratio, m_block_factor, 
                        m_nesting_radius, m_max_box_size);
    }
  
  Vector<IntVectSet> tagVect(m_max_level);
  
  Vector<Vector<Box> > oldBoxes(1);
  Vector<Vector<Box> > newBoxes;
  oldBoxes[0] = baseBoxes;
  newBoxes = oldBoxes;
  int new_finest_level = 0;

  while (moreLevels)
    {
      // default is moreLevels = false
      // (only repeat loop in the case where a new level is generated
      // which is still coarser than maxLevel)
      moreLevels = false;
      tagCellsInit(tagVect);
      
      // two possibilities -- need to generate grids
      // level-by-level, or we are refining all the
      // way up for the initial time.  check to 
      // see which it is by seeing if the finest-level
      // tags are empty
      if (tagVect[m_max_level-1].isEmpty())
        {
          int top_level = m_finest_level;
          int old_top_level = top_level;
          new_finest_level = meshrefine.regrid(newBoxes,
                                               tagVect, baseLevel,
                                               top_level,
                                               oldBoxes);

          if (new_finest_level > top_level) top_level++;
          oldBoxes = newBoxes;

          // now see if we need another pass through grid generation
          if ((top_level < m_max_level) && (top_level > old_top_level) && (new_finest_level <= m_tag_cap))
            {
              moreLevels = true;
            }
          
        }
      else 
        {
          
          // for now, define old_grids as just domains
          oldBoxes.resize(m_max_level+1);
          for (int lev=1; lev<=m_max_level; lev++) 
            {
              oldBoxes[lev].push_back(m_amrDomains[lev].domainBox());
            }
          
          int top_level = m_max_level -1;
          new_finest_level = meshrefine.regrid(newBoxes,
                                               tagVect, baseLevel,
                                               top_level,
                                               oldBoxes);
        }
      
      numLevels = Min(new_finest_level, m_max_level)+1;
      

      // now loop through levels and define
      for (int lev=baseLevel+1; lev<= new_finest_level; ++lev)
        {
          int numGridsNew = newBoxes[lev].size();
          Vector<int> procIDs(numGridsNew);
          LoadBalance(procIDs, newBoxes[lev]);
          const DisjointBoxLayout newDBL(newBoxes[lev], procIDs,
                                         m_amrDomains[lev]);
          m_amrGrids[lev] = newDBL;

          if (s_verbosity > 2)
            {
              long long levelNumCells = newDBL.numCells();          
              pout() << "   Level " << lev << ": " 
                     << levelNumCells << " cells: " 
                     << m_amrGrids[lev] << endl;
            }
              

          levelSetup(lev,m_amrGrids[lev]);
	  m_A_valid = false;
	  m_groundingLineProximity_valid = false;
	  m_viscousTensor_valid = false;

        } // end loop over levels

      m_finest_level = new_finest_level;
      
      // finally, initialize data on final hierarchy
      // only do this if we've created new levels
      if (m_finest_level > 0) 
        {
          defineSolver();

          initData(m_vect_coordSys,
                   m_velocity);
        }
    } // end while more levels to do

  


}


void
AmrIce::setupFixedGrids(const std::string& a_gridFile)
{
  Vector<Vector<Box> > gridvect;
  
  if (procID() == uniqueProc(SerialTask::compute))
    {
      gridvect.push_back(Vector<Box>(1,m_amrDomains[0].domainBox()));
    
      // read in predefined grids
      ifstream is(a_gridFile.c_str(), ios::in);
      
      if (is.fail())
        {
          MayDay::Error("Cannot open grids file");
        }

      // format of file:
      //   number of levels, then for each level (starting with level 1):
      //   number of grids on level, list of boxes
      int inNumLevels;
      is >> inNumLevels;

      CH_assert (inNumLevels <= m_max_level+1);

      if (s_verbosity > 3)
        {
          pout() << "numLevels = " << inNumLevels << endl;
        }

      while (is.get() != '\n');

      gridvect.resize(inNumLevels);

      // check to see if coarsest level needs to be broken up
      domainSplit(m_amrDomains[0],gridvect[0], m_max_base_grid_size, 
                  m_block_factor);

      if (s_verbosity >= 3)
        {
          pout() << "level 0: ";
          for (int n=0; n < gridvect[0].size(); n++)
            {
              pout() << gridvect[0][n] << endl;
            }
        }

      // now loop over levels, starting with level 1
      int numGrids = 0;
      for (int lev=1; lev<inNumLevels; lev++) 
        {
          is >> numGrids;

          if (s_verbosity >= 3)
            {
              pout() << "level " << lev << " numGrids = " 
                     << numGrids <<  endl;
              pout() << "Grids: ";
            }

          while (is.get() != '\n');

          gridvect[lev].resize(numGrids);

          for (int i=0; i<numGrids; i++)
            {
              Box bx;
              is >> bx;

              while (is.get() != '\n');

              // quick check on box size
              Box bxRef(bx);

              if (bxRef.longside() > m_max_box_size)
                {
                  pout() << "Grid " << bx << " too large" << endl;
                  MayDay::Error();
                }

              if (s_verbosity >= 3) 
                {
                  pout() << bx << endl;
                }

              gridvect[lev][i] = bx;
            } // end loop over boxes on this level
        } // end loop over levels
    } // end if serial proc

  // broadcast results
  broadcast(gridvect, uniqueProc(SerialTask::compute));

  // now create disjointBoxLayouts and allocate grids

  m_amrGrids.resize(m_max_level+1);
  IntVect sigmaCSGhost = m_num_thickness_ghost*IntVect::Unit;
  m_vect_coordSys.resize(m_max_level+1);
  
  // probably eventually want to do this differently
  RealVect dx = m_amrDx[0]*RealVect::Unit;

  for (int lev=0; lev<gridvect.size(); lev++)
    {
      int numGridsLev = gridvect[lev].size();
      Vector<int> procIDs(numGridsLev);
      LoadBalance(procIDs, gridvect[lev]);
      const DisjointBoxLayout newDBL(gridvect[lev],
                                     procIDs, 
                                     m_amrDomains[lev]);

      m_amrGrids[lev] = newDBL;

      // build storage for this level

      levelSetup(lev, m_amrGrids[lev]);
      if (lev < gridvect.size()-1)
        {
          dx /= m_refinement_ratios[lev];
        }
    }
  
  // finally set finest level and initialize data on hierarchy
  m_finest_level = gridvect.size() -1;

  // define solver before calling initData
  defineSolver();
  
  initData(m_vect_coordSys, m_velocity);

}
    

void
AmrIce::levelSetup(int a_level, const DisjointBoxLayout& a_grids)
{
  IntVect ghostVect = IntVect::Unit;
  // 4 ghost cells needed for advection. Could later go back and
  // make this a temporary if the additional storage becomes an issue...
  IntVect thicknessGhostVect = m_num_thickness_ghost*IntVect::Unit;
  m_old_thickness[a_level]->define(a_grids, 1,
                                   thicknessGhostVect);

  if (a_level == 0 || m_velocity[a_level] == NULL)
    {
      m_velocity[a_level] = new LevelData<FArrayBox>(a_grids, 2,
                                                     ghostVect);
    }
  else
    {
      // do velocity a bit differently in order to use previously 
      // computed velocity field as an initial guess
      {
        LevelData<FArrayBox>* newVelPtr = new LevelData<FArrayBox>(a_grids,
                                                                   2,
                                                                   ghostVect);
        
        // first do interp from coarser level
        FineInterp velInterp(a_grids, 2, 
                             m_refinement_ratios[a_level-1],
                             m_amrDomains[a_level]);
        
        velInterp.interpToFine(*newVelPtr, *m_velocity[a_level-1]);
        
        // can only copy from existing level if we're not on the
        // newly created level
        //if (a_level != new_finest_level)
        if (m_velocity[a_level]->isDefined())
          {
            m_velocity[a_level]->copyTo(*newVelPtr);
          }
        
        // finally, do an exchange (this may wind up being unnecessary)
        newVelPtr->exchange();
        
        delete (m_velocity[a_level]);
        m_velocity[a_level] = newVelPtr;
      }
    } // end interpolate/copy new velocity


  m_faceVel[a_level] = new LevelData<FluxBox>(a_grids,1,IntVect::Unit);
#if BISICLES_Z == BISICLES_LAYERED
  m_layerXYFaceXYVel[a_level] = new LevelData<FluxBox>(a_grids, m_nLayers, IntVect::Unit);
  m_layerSFaceXYVel[a_level] = new LevelData<FArrayBox>(a_grids, SpaceDim*(m_nLayers + 1), IntVect::Unit);
#endif

  m_velBasalC[a_level] = new LevelData<FArrayBox>(a_grids, 1, ghostVect);
  m_cellMuCoef[a_level] = new LevelData<FArrayBox>(a_grids, 1, ghostVect);
  m_faceMuCoef[a_level] = new LevelData<FluxBox>(a_grids, 1, ghostVect);

  m_diffusivity[a_level] = RefCountedPtr<LevelData<FluxBox> >
    (new LevelData<FluxBox>(a_grids, 1, IntVect::Zero));

  m_velRHS[a_level] = new LevelData<FArrayBox>(a_grids, 2, 
					 IntVect::Zero);
 
  m_surfaceThicknessSource[a_level] = 
    new LevelData<FArrayBox>(a_grids,   1, IntVect::Unit) ;
  m_basalThicknessSource[a_level] = 
    new LevelData<FArrayBox>(a_grids,   1, IntVect::Unit) ;
  m_balance[a_level] = 
    new LevelData<FArrayBox>(a_grids,   1, IntVect::Zero) ;

  // probably eventually want to do this differently
  RealVect dx = m_amrDx[a_level]*RealVect::Unit;

  //IntVect sigmaCSGhost = IntVect::Unit;
  IntVect sigmaCSGhost = thicknessGhostVect;
  m_vect_coordSys.resize(m_max_level+1);
  m_vect_coordSys[a_level] = RefCountedPtr<LevelSigmaCS >(new LevelSigmaCS(a_grids, 
                                                                     dx,
                                                                     sigmaCSGhost));
  
#if BISICLES_Z == BISICLES_LAYERED
  //in poor-man's multidim mode, use one FArrayBox component per layer
  //to hold the 3D temperature field
  m_temperature[a_level] = new LevelData<FArrayBox>(a_grids, m_nLayers, 
					      thicknessGhostVect);
  // plus base and surface temperatures
  m_sTemperature[a_level] = new LevelData<FArrayBox>(a_grids, 1, 
					       thicknessGhostVect);
  m_bTemperature[a_level] = new LevelData<FArrayBox>(a_grids, 1, 
					       thicknessGhostVect);
  
  m_vect_coordSys[a_level]->setFaceSigma(getFaceSigma());
  m_vect_coordSys[a_level]->setIceDensity(m_iceDensity);
  m_vect_coordSys[a_level]->setWaterDensity(m_seaWaterDensity);
  m_vect_coordSys[a_level]->setGravity(m_gravity);
#elif BISICLES_Z == BISICLES_FULLZ
  m_temperature[a_level] = new LevelData<FArrayBox>(a_grids, 1, thicknessGhostVect);	
#endif




}

void
AmrIce::initData(Vector<RefCountedPtr<LevelSigmaCS> >& a_vectCoordSys,
                 Vector<LevelData<FArrayBox>* >& a_velocity)
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::initData" << endl;
    }

  m_groundingLineProximity_valid = false;
  m_A_valid = false;

  for (int lev=0; lev<=m_finest_level; lev++)
    {

      if (lev > 0)
	{
	  // fill the ghost cells of a_vectCoordSys[lev]->getH();
	  LevelData<FArrayBox>& levelH = a_vectCoordSys[lev]->getH();
	  LevelData<FArrayBox>& coarseH = a_vectCoordSys[lev-1]->getH();
	  int nGhost = levelH.ghostVect()[0];
	  PiecewiseLinearFillPatch thicknessFiller
	    (m_amrGrids[lev],m_amrGrids[lev-1],1, m_amrDomains[lev-1],
	     m_refinement_ratios[lev-1], nGhost);
	  thicknessFiller.fillInterp(levelH,coarseH,coarseH,0.0, 0, 0, 1);

          // do the same with topography
          LevelData<FArrayBox>& levelTopog = a_vectCoordSys[lev]->getTopography();
          LevelData<FArrayBox>& crseTopog = a_vectCoordSys[lev-1]->getTopography();

          // should be able to use the same filler
          CH_assert(levelTopog.ghostVect() == levelH.ghostVect());
          thicknessFiller.fillInterp(levelTopog, crseTopog, crseTopog, 0.0, 0, 0, 1);
          
	}

      RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
      m_thicknessIBCPtr->define(m_amrDomains[lev],levelDx[0]);
      LevelSigmaCS* crsePtr = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
      int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:0;
      
      m_thicknessIBCPtr->initializeIceGeometry(*a_vectCoordSys[lev],
					       levelDx,
					       m_domainSize,
					       m_time,
					       crsePtr,
					       refRatio);
      

	


      a_vectCoordSys[lev]->recomputeGeometry(crsePtr, refRatio);

      //allow calving model to modify geometry and velocity
      m_calvingModelPtr->initialModifyState
      	(m_vect_coordSys[lev]->getH(), *this, lev);

      a_vectCoordSys[lev]->recomputeGeometry(crsePtr, refRatio);

      // initialize oldH to be the current value
      LevelData<FArrayBox>& currentH = a_vectCoordSys[lev]->getH();
      currentH.copyTo(*m_old_thickness[lev]);

#if BISICLES_Z == BISICLES_LAYERED
      m_temperatureIBCPtr->initializeIceTemperature
	(*m_temperature[lev], *m_sTemperature[lev], *m_bTemperature[lev],*a_vectCoordSys[lev] );

#elif BISICLES_Z == BISICLES_FULLZ
      m_temperatureIBCPtr->initializeIceTemperature(*m_temperature[lev],*a_vectCoordSys[lev]);
#endif
    }
  
  // now call velocity solver to initialize velocity field
  solveVelocityField(m_velocity);

  // may be necessary to average down here
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverage avgDown(m_amrGrids[lev],
                            2, m_refinement_ratios[lev-1]);
      avgDown.averageToCoarse(*m_velocity[lev-1], *m_velocity[lev]);
    }


#define writePlotsImmediately
#ifdef  writePlotsImmediately
  if (m_plot_interval >= 0)
    {
      writePlotFile();
    }
#endif


}

/// solve for velocity field
void
AmrIce::solveVelocityField(Vector<LevelData<FArrayBox>* >& a_velocity,
			   Real a_convergenceMetric)
{

  CH_TIME("AmrIce::solveVelocityField");

  //ensure A is up to date
#if BISICLES_Z == BISICLES_LAYERED
  if (!m_A_valid)
	  {
	    computeA(m_A,m_sA,m_bA,m_temperature,m_sTemperature,m_bTemperature,
		     m_vect_coordSys);
	    m_A_valid = true;
	  }
#else
  MayDay::Error("AmrIce::SolveVelocityField full z calculation of A not done"); 
#endif
  //certainly the viscous tensr field will need re-computing
  m_viscousTensor_valid = false;

  // define basal friction
  Vector<LevelData<FArrayBox>* > vectC(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectC0(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectRhs(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      vectRhs[lev] = m_velRHS[lev];
      vectC[lev] = m_velBasalC[lev];
      vectC0[lev] = new LevelData<FArrayBox>; vectC0[lev]->define(*vectC[lev]);
      LevelData<FArrayBox>& C0 = *vectC0[lev];
      for (DataIterator dit=C0.dataIterator(); dit.ok(); ++dit)
	{
	  C0[dit].setVal(0.0);
	}
    }

  //setBasalFriction(vectC);
  setMuCoefficient(m_cellMuCoef,m_faceMuCoef);
  //eliminateRemoteIce();
  setBasalFriction(vectC);
  //setMuCoefficient(m_cellMuCoef,m_faceMuCoef);

  // also sets beta=0 where ice is floating 
  defineVelRHS(vectRhs, vectC,vectC0);
  
  // write out sumRhs if appropriate
  if (s_verbosity > 3)
    {
      Real sumRhs = computeSum(vectRhs,
                               m_refinement_ratios,
                               m_amrDx[0],
                               Interval(0,0),
                               0);

      pout() << "Sum(rhs) for velocity solve = " << sumRhs << endl;

    }
  
  // put this in place to catch runs where plotfile writing is
  // going to hang _before_ I waste a few hours waiting for the 
  // velocity solve
  //#define writeTestPlots
#ifdef  writeTestPlots
  if (m_plot_interval >= 0 && m_cur_step == 0)
    {
      writePlotFile();
    }
#endif

  

  if (m_doInitialVelSolve) 
    {      
      if (m_finest_level == 0 && m_doInitialVelGuess)
        {
          // only really want or need to do this once
          m_doInitialVelGuess = false;
	  if (m_initialGuessType == SlidingLaw)
	    {
	      pout() << "computing an initial guess via a sliding law u = rhs/C "  << endl;
	      // compute initial guess as rhs/beta
	      LevelData<FArrayBox>& vel = *m_velocity[0];
	      LevelData<FArrayBox>& C = *m_velBasalC[0];
	      LevelData<FArrayBox>& rhs = *m_velRHS[0];
	      const DisjointBoxLayout& levelGrids = m_amrGrids[0];
	      
	      DataIterator dit = vel.dataIterator();
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  FORT_VELINITIALGUESS(CHF_FRA(vel[dit]),
				       CHF_FRA(rhs[dit]),
				       CHF_FRA1(C[dit],0),
				       CHF_BOX(levelGrids[dit]));
		}
	    }
	  else if (m_initialGuessType == ConstMu)
	    {
	      if (s_verbosity > 3) 
		{
		  pout() << "computing an initial guess by solving the velocity equations "
			 <<" with constant mu = " 
			 << m_initialGuessConstMu   
			 << "and constant initial velcocity = " << m_initialGuessConstVel
			 << endl;
		}
	      // compute initial guess by solving a linear problem with a
	      // modest constant viscosity
	      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;
	      constMuRelation* newPtr = new constMuRelation;
	      newPtr->setConstVal(m_initialGuessConstMu);
	      for (int lev=0; lev < m_finest_level + 1; lev++)
		{
		  for (DataIterator dit(m_amrGrids[lev]);dit.ok();++dit)
		    {
		      for (int dir = 0; dir < SpaceDim; dir++)
			{
			  (*m_velocity[lev])[dit].setVal(m_initialGuessConstVel[dir],dir);
			}
		    }

		}
              // do this by saving the exisiting velSolver and 
              // constitutiveRelation, re-calling defineSolver, then 
              // doing solve.
              IceVelocitySolver* velSolverSave = m_velSolver;
              ConstitutiveRelation* constRelSave = m_constitutiveRelation;
              int solverTypeSave = m_solverType;

              // new values prior to calling defineSolver
             
              m_constitutiveRelation = static_cast<ConstitutiveRelation*>(newPtr);
             
	      Real finalNorm = 0.0, initialNorm = 0.0, convergenceMetric = -1.0;
	      Vector<LevelData<FluxBox>* > muCoef(m_finest_level + 1,NULL);
	      int rc;
	      if (m_initialGuessSolverType == JFNK)
		{
		  //JFNK can be instructed to assume a linear solve
		  m_velSolver = NULL;
		  defineSolver();
		  JFNKSolver* jfnkSolver = dynamic_cast<JFNKSolver*>(m_velSolver);
		  CH_assert(jfnkSolver != NULL);
		  const bool linear = true;
		  rc = jfnkSolver->solve(m_velocity, initialNorm,finalNorm,convergenceMetric,
					  linear , m_velRHS, m_velBasalC, vectC0, m_A, muCoef,
					  m_vect_coordSys, m_time, 0, m_finest_level);
		}
	      else if (m_initialGuessSolverType == Picard)
		{
		  // since the constant-viscosity solve is a linear solve,
		  // Picard is the best option.
		  m_solverType = Picard;
		  m_velSolver = NULL;
		  defineSolver();
		  rc = m_velSolver->solve(m_velocity, initialNorm,finalNorm,convergenceMetric,
				     m_velRHS, m_velBasalC, vectC0, m_A, muCoef,
				     m_vect_coordSys, m_time, 0, m_finest_level);
		}
	      else
		{
		  MayDay::Error("unknown initial guess solver type");
		}


	      if (rc != 0)
		{
		  MayDay::Warning("constant mu solve failed");
		}

              // now put everything back the way it was...
	      delete m_constitutiveRelation;
              delete m_velSolver;
              m_velSolver = velSolverSave;
              m_constitutiveRelation = constRelSave;
              m_solverType = solverTypeSave;
              
#if 0	      
	      //put the solver back how it was
	      m_velSolver->define(m_amrDomains[0],
				  m_constitutiveRelation,
				  m_basalFrictionRelation,
				  m_amrGrids,
				  m_refinement_ratios,
				  dxCrse,
				  m_thicknessIBCPtr,
				  numLevels);
#endif
      
	    }
	  else if (m_initialGuessType == Function)
	    {
	      ParmParse pp("amr");
	      std::string functionType = "constant";
	      pp.query("initial_velocity_function_type", functionType );
	      if (functionType == "flowline")
		{
		  Real dx; 
		  std::string file, set;
		  pp.get("initial_velocity_function_flowline_dx", dx);
		  pp.get("initial_velocity_function_flowline_file", file);
		  pp.get("initial_velocity_function_flowline_set", set);
		  ExtrudedPieceWiseLinearFlowline f(file,set,dx);
		  for (int lev = 0; lev <  m_finest_level + 1; lev++)
		    {
		      LevelData<FArrayBox>& levelVel = *m_velocity[lev];
		      for (DataIterator dit(levelVel.disjointBoxLayout());
			   dit.ok(); ++dit)
			{
			  FArrayBox& vel = levelVel[dit];
			  const Box& box = vel.box(); 
			  for (BoxIterator bit(box); bit.ok(); ++bit)
			    {
			      const IntVect& iv = bit();
			      RealVect x = RealVect(iv) * m_amrDx[lev] 
				+ 0.5 * m_amrDx[lev];
			      vel(iv,0) = f(x);
			    }
			}
		    }
		}
	      
	    }
	  else
	    {
	      MayDay::Error("AmrIce::SolveVelocityField unknown initial guess type");
	    }
	}
     
      if (m_write_presolve_plotfiles)
        {
          string save_prefix = m_plot_prefix;
          m_plot_prefix.append("preSolve.");
	  bool t_write_fluxVel = m_write_fluxVel;
	  m_write_fluxVel = false; // turning this off in preSolve files for now
          writePlotFile();
          m_write_fluxVel = t_write_fluxVel;
	  m_plot_prefix = save_prefix;
        }
    

      int solverRetVal; 
      // if (m_isothermal)
      // 	{
      // 	  // no need to recompute A
      // 	  solverRetVal = m_velSolver->solve(m_velocity,
      // 					    m_velRHS, m_velBasalC,
      // 					    m_A,
      // 					    m_vect_coordSys,
      // 					    m_time,
      // 					    0, m_finest_level);
      // 	}
      // else
      {
	//Vector<LevelData<FluxBox>* > muCoef(m_finest_level + 1,NULL);
	{
	  
	  for (int lev=0; lev <= m_finest_level ; ++lev)
	    {
	      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
	      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
		{
		  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
		  FArrayBox& vel = (*m_velocity[lev])[dit];
		  for (BoxIterator bit(levelGrids[dit]); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      if (mask(iv) == OPENSEAMASKVAL || 
			  mask(iv) == OPENLANDMASKVAL )
			{
			  vel(iv,0) = 0.0; vel(iv,1) = 0.0;
			} 
		    }
		}
	    }
	}

	solverRetVal = m_velSolver->solve(m_velocity,
					  m_velocitySolveInitialResidualNorm, 
					  m_velocitySolveFinalResidualNorm,
					  a_convergenceMetric,
					  m_velRHS, m_velBasalC, vectC0,
					  m_A, m_faceMuCoef,
					  m_vect_coordSys,
					  m_time,
					  0, m_finest_level);
	  
      }
    }

  // set the thickness diffusivity: not really part
  // of the velocity solver but does depend on beta
  if (m_diffusionTreatment == IMPLICIT || m_diffusionTreatment == EXPLICIT)
    setThicknessDiffusivity(vectC);


  //calculate the face centred (flux) velocity 
  computeFaceVelocity(m_faceVel, m_layerXYFaceXYVel, m_layerSFaceXYVel);

  

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      if (vectC0[lev] != NULL)
	{
	  delete vectC0[lev]; vectC0[lev] = NULL;
	}
    }
#if 0  
  // debugging test -- redefine velocity as a constant field
  for (int lev=0; lev<m_velocity.size(); lev++)
    {
      DataIterator dit = m_velocity[lev]->dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*m_velocity[lev])[dit].setVal(1.0);
        }
    }
#endif

}

/// compute RHS for velocity field solve
void
AmrIce::defineVelRHS(Vector<LevelData<FArrayBox>* >& a_vectRhs,
                     Vector<LevelData<FArrayBox>* >& a_vectC,
		     Vector<LevelData<FArrayBox>* >& a_vectC0)
{
  // ice density from Pattyn(2003) in kg/m^3
  //Real rhoIce = 910.0;

  // gravitational acceleration in m/s^2
  //Real g = 9.81;
  
  
  // we will need ghost cells for H when we come to to
  // FORT_GLCORRECTION. 
  Vector<LevelData<FArrayBox>* >tempH(m_finest_level+1);
  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
     
      tempH[lev] = (&levelCoords.getH());
      LevelData<FArrayBox>& levelH = *tempH[lev];

      if (lev > 0)
	{
	  PiecewiseLinearFillPatch ghostFiller(levelGrids, 
					       m_amrGrids[lev-1],
					       1, 
					       m_amrDomains[lev-1],
					       m_refinement_ratios[lev-1],
					       1);

	  const LevelData<FArrayBox>& coarseH = *tempH[lev-1]; 
	  ghostFiller.fillInterp(levelH, coarseH, coarseH, 0.0, 0, 0, 1);
	}
      
      levelH.exchange();
    }


  // construct multilevel surface height
  Vector<LevelData<FArrayBox>* > zSurf(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const DisjointBoxLayout& grids = m_amrGrids[lev];
      zSurf[lev] = new LevelData<FArrayBox>(grids, 1, IntVect::Unit);
      levelCS.getSurfaceHeight(*zSurf[lev]);
    }
  
  // used for limiting RHS, if needed
  Real maxZs = 0.0;
  Real maxGrad = -1.0;
  if (m_limitVelRHS && (m_gradLimitRadius > 0) )
    {
      Interval comps(0,0);
      maxZs = computeMax(zSurf, m_refinement_ratios, comps, 0);
      maxGrad = maxZs/Real(m_gradLimitRadius);
    }
  

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const LevelData<FArrayBox>& levelH = levelCS.getH();
      
      const LevelData<FArrayBox>& levelGradS = levelCS.getGradSurface();
      const LevelData<BaseFab<int> >& levelMask = levelCS.getFloatingMask();
      const LayoutData<bool>& levelAnyFloating = levelCS.anyFloating();
      LevelData<FArrayBox>& levelC = (*a_vectC[lev]);
      LevelData<FArrayBox>& levelRhs = (*a_vectRhs[lev]);
      const DisjointBoxLayout& grids = m_amrGrids[lev];
      LevelData<FArrayBox>& levelZs = *zSurf[lev];

      const RealVect& dx = levelCS.dx();
      Real iceDensity = levelCS.iceDensity();
      Real gravity = levelCS.gravity();
      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& thisH = levelH[dit];
	  const BaseFab<int>& floatingMask = levelMask[dit];
	  bool  anyFloating = levelAnyFloating[dit];  
          // now take grad(z_s)
          FArrayBox& thisRhs = levelRhs[dit];
	  const FArrayBox& grad = levelGradS[dit];
	  const Box& gridBox = grids[dit];
	  
          

          // now add in background slope 
	  //and compute RHS.	  
          // break this into two steps and limit grad(z_s) if 
          // necessary
          thisRhs.copy(grad, gridBox);
          if (m_limitVelRHS && (m_gradLimitRadius > 0) )
            {
              RealVect maxGradDir;
              maxGradDir[0] = maxGrad / dx[0];
              maxGradDir[1] = maxGrad / dx[1];
              Box gridBoxPlus = gridBox;
              gridBoxPlus.grow(1);

              for (BoxIterator bit(gridBox); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  // don't limit in floating regions, since 
                  // grad(s) here represents the marine boundary condition
                  if  (floatingMask(iv) == GROUNDEDMASKVAL)
                    {
                      for (int dir=0; dir<2; dir++)
                        {
                          if (thisRhs(iv,dir) > maxGradDir[dir])
                            thisRhs(iv,dir) = maxGradDir[dir];
                          else if (thisRhs(iv,dir) < -maxGradDir[dir])
                            thisRhs(iv,dir) = -maxGradDir[dir];
                        }
                    }                          
                }
            }
          for (BoxIterator bit (gridBox); bit.ok(); ++bit)
            {
              IntVect iv = bit();

              for (int dir=0; dir<2; dir++)
                {     
                  thisRhs(iv,dir) *= iceDensity*gravity*thisH(iv,0);
                } // end loop over directions
            } // end loop over cells in this box



	  CH_assert(thisRhs.norm(0,0,SpaceDim) < HUGE_NORM);

	 

	  if(m_groundingLineCorrection && anyFloating == 1 
	     && 0.0 * m_basalLengthScale < dx[0])
	    {
	      Real rhog = iceDensity*gravity;
	      CH_assert(SpaceDim != 3);
	      for (int dir=0; dir<SpaceDim; dir++)
		{
		  const FArrayBox& thisH = (*tempH[lev])[dit];
		  const FArrayBox& thisZsurf = levelZs[dit];
		  FORT_GLCORRECTION(CHF_FRA1(thisRhs, dir),
				    CHF_CONST_FRA1(thisH,0),
				    CHF_CONST_FRA1(thisZsurf,0),
				    CHF_CONST_FIA1(floatingMask,0),
				    CHF_INT(dir),
				    CHF_CONST_REAL(dx[dir]),
				    CHF_CONST_REAL(rhog),
				    CHF_BOX(gridBox));
		}
	    }

	  
	  CH_assert(thisRhs.norm(0,0,SpaceDim) < HUGE_NORM);

	  
	  FArrayBox& thisC = levelC[dit];
	  FArrayBox& thisC0 = (*a_vectC0[lev])[dit];

	

	  // add drag due to ice in contact with ice-free rocky walls
	  
	  thisC0.setVal(0.0);
	  if (m_wallDrag)
	  {
	    IceVelocity::addWallDrag(thisC0, floatingMask, levelCS.getSurfaceHeight()[dit],
				     thisH, levelCS.getBaseHeight()[dit], thisC, m_wallDragExtra,
				     m_amrDx[lev], gridBox);
	  }
	  

	    // this is also a good place to set C=0 where
          // ice is floating, 
          if(anyFloating)
            {
	      
	      const FArrayBox& thisTopography = levelCS.getBaseHeight()[dit];
	      FArrayBox thisLowerSurface(gridBox,1);
	      if (m_basalFrictionDecay > 0.0)
		{
		  thisLowerSurface.copy(levelCS.getSurfaceHeight()[dit]);
		  thisLowerSurface -= thisH;
		}
              FORT_SETFLOATINGBETA(CHF_FRA1(thisC,0),
                                   CHF_CONST_FIA1(floatingMask,0),
				   CHF_CONST_FRA1(thisLowerSurface,0),
				   CHF_CONST_FRA1(thisTopography,0),
				   CHF_CONST_REAL(m_basalFrictionDecay),
                                   CHF_BOX(gridBox));

	      CH_assert(thisC.min(gridBox) >= 0.0); 
              // one-sided differencing near grounding line. 
              // Would slope limiting be sufficient?          

              // looping over directions not quite right in 3d
              CH_assert(SpaceDim != 3);
     
            } // end if anything is floating

	 

	   
        } // end loop over boxes

      // finally, modify RHS in problem-dependent ways,
      m_thicknessIBCPtr->modifyVelocityRHS(levelRhs, 
                                           *m_vect_coordSys[lev],
                                           m_amrDomains[lev],
                                           m_time, m_dt);

    } // end loop over levels


  
  if ( m_basalLengthScale > TINY_NORM ){
    if (s_verbosity > 1)
      {
	pout() << "Smoothing velocity basal C" << std::endl;
      }
    helmholtzSolve(m_velBasalC,Real(1.0),std::pow(m_basalLengthScale,2));
        if (s_verbosity > 1)
      {
	pout() << "Smoothing velocity rhs" << std::endl;
      }

    helmholtzSolve(m_velRHS,Real(1.0),std::pow(m_basalLengthScale,2));
  }

  // clean up
  for (int lev=0; lev<zSurf.size(); lev++)
    {
      if (zSurf[lev] != NULL)
        {
          delete zSurf[lev];
          zSurf[lev] = NULL;
        }
    }

}

/// set mu coefficient (phi) prior to velocity solve
void
AmrIce::setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef, 
			 Vector<LevelData<FluxBox>* >& a_faceMuCoef)
{
  CH_assert(m_muCoefficientPtr != NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      m_muCoefficientPtr->setMuCoefficient(*a_cellMuCoef[lev],
					   *a_faceMuCoef[lev],
					   *m_vect_coordSys[lev],
                                           m_time,
                                           m_dt);
    }

}


/// set basal friction coefficient (beta) prior to velocity solve
void
AmrIce::setBasalFriction(Vector<LevelData<FArrayBox>* >& a_vectBeta)
{
  CH_assert(m_basalFrictionPtr != NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      m_basalFrictionPtr->setBasalFriction(*a_vectBeta[lev],
                                           *m_vect_coordSys[lev],
                                           m_time,
                                           m_dt);

      // slc : not needed as velocity is to be measured in m/a
      // DataIterator dit = (a_vectBeta[lev]->getBoxes()).dataIterator();
      // // convert from Pa a/m -> Pa*s/m
      // for (dit.begin(); dit.ok(); ++dit)
      //   {
	  
      //     (*a_vectBeta[lev])[dit]  *= secondsperyear;
      //   }
      
      a_vectBeta[lev]->exchange();
    }
}

/// compute the thickness diffusivity D = H^2 / C
/// needs to change when we move on from linear Friction laws
void 
AmrIce::setThicknessDiffusivity(const Vector<LevelData<FArrayBox>* >& a_beta)
{
  
  Real rhog = m_vect_coordSys[0]->iceDensity() * m_vect_coordSys[0]->gravity();
  
  Vector<LevelData<FArrayBox>* > vectH(m_finest_level+1);

  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      LevelData<FluxBox>& levelD = *m_diffusivity[lev];
      const LevelData<FArrayBox>& levelBeta = *a_beta[lev];
      const LevelData<FArrayBox>& levelVel = *m_velocity[lev];
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
      vectH[lev] = const_cast<LevelData<FArrayBox>* >(&levelCoords.getH());

      //we need one ghost cell worth of H, which LevelSigmaCS
      // should hopefully give us.
      LevelData<FArrayBox>& levelH = *vectH[lev];
      CH_assert(levelH.ghostVect() > IntVect::Unit);
      
      if (lev > 0)
	{
	  PiecewiseLinearFillPatch ghostFiller(levelGrids, 
					       m_amrGrids[lev-1],
					       1, 
					       m_amrDomains[lev-1],
					       m_refinement_ratios[lev-1],
					       1);

	  const LevelData<FArrayBox>& coarseH = *vectH[lev-1]; 
	  ghostFiller.fillInterp(levelH, coarseH, coarseH, 0.0, 0, 0, 1);
	}
     
      levelH.exchange();

      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const FArrayBox& H = levelH[dit];
	  const FArrayBox& cellBeta = levelBeta[dit];
	  Box grownBox = levelGrids[dit];
	  grownBox.grow(1);
	  FArrayBox cellC(grownBox,1);;
	  
	  m_basalFrictionRelation->computeAlpha(cellC,levelVel[dit], 
						levelCoords.getThicknessOverFlotation()[dit],
						cellBeta, 
						levelCoords.getFloatingMask()[dit],
						grownBox);
	  //cellC *= cellBeta;

	  for (int dir = 0; dir < SpaceDim; ++dir)
	    {
	     
	      FArrayBox& faceD = levelD[dit][dir];
	      faceD.setVal(1.0e-12);
	      
	      Box faceBox = levelGrids[dit];
	      // Box inner = m_amrDomains[lev].domainBox();
	      // inner.grow(dir,-1);
	      // faceBox &= inner;
	      faceBox.surroundingNodes(dir);
	      FORT_SETTHICKDIFF(CHF_FRA1(faceD,0),
	      			CHF_CONST_FRA1(cellC,0),
	      			CHF_CONST_FRA1(H,0),
	      			CHF_CONST_FIA1(levelMask[dit],0),
	      			CHF_CONST_REAL(rhog),
	      			CHF_CONST_INT(dir),
	      			CHF_BOX(faceBox));
	      CH_assert(faceD.min(faceBox) >= 0.0);
	      
	      
	    }
	}
      
      levelD.exchange();
    }

}

// compute diffusive velocity (- D/H grad(H)) and source div(D grad H)
void 
AmrIce::computeDiffusionTerms(LevelData<FluxBox>& a_diffVel,
			      LevelData<FArrayBox>& a_diffSrc,
			      int a_lev) const
{
  const int& lev = a_lev;
  const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
  LevelData<FluxBox>& levelVel = a_diffVel;
  LevelData<FArrayBox>& levelSrc = a_diffSrc;
  const LevelData<FArrayBox>& levelH = levelCoords.getH();
  const LevelData<FluxBox>& levelFaceH = levelCoords.getFaceH();
  const LevelData<FluxBox>& levelD = *m_diffusivity[lev];  
  const RealVect& dx = levelCoords.dx();

  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
    {

      FArrayBox& cellSrc = levelSrc[dit];
      cellSrc.setVal(0.0);
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  FArrayBox& faceVel = levelVel[dit][dir];
	  faceVel.setVal(0.0);
	  const FArrayBox& faceD = levelD[dit][dir];
	  
	  const FArrayBox& cellH = levelH[dit];
	  const FArrayBox& faceH = levelFaceH[dit][dir];
	  Box faceBox = levelGrids[dit];
	  faceBox.surroundingNodes(dir);
	  FArrayBox gradHonH(faceBox,1);

	  Real oneOnDx = 1.0 / dx[dir];
	  for (BoxIterator bit(faceBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	  
	      gradHonH(iv) =  oneOnDx / std::max(faceH(iv),1.0) *
		(cellH(iv) - cellH(iv-BASISV(dir)));  

	      faceVel(iv) = -faceD(iv)*gradHonH(iv);
	    }
	 
	  const Box& cellBox = levelGrids[dit];
	  for (BoxIterator bit(cellBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      cellSrc(iv) += oneOnDx * (faceVel(iv+BASISV(dir))*faceH(iv+BASISV(dir))
					- faceVel(iv)*faceH(iv));
	    }
	}

      
    }

}



/// subtract - D grad(H) from a_flux
void 
AmrIce::subtractDiffusiveFlux(Vector<LevelData<FluxBox>*>& a_vectFluxes) const
{

  Vector<LevelData<FArrayBox>* > vectH(m_finest_level+1);
  for (int lev=0; lev <= finestTimestepLevel() ; ++lev)
    {
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      //we need one ghost cell worth of H, which SigmaCS doesn't give us yet,
      //though it could once we move to LevelSigmaCS
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      vectH[lev] = const_cast<LevelData<FArrayBox>*>(&levelCoords.getH());
      LevelData<FArrayBox>& levelH = *vectH[lev];
      
      if (lev > 0)
	{
	  PiecewiseLinearFillPatch ghostFiller(levelGrids, 
					       m_amrGrids[lev-1],
					       1, 
					       m_amrDomains[lev-1],
					       m_refinement_ratios[lev-1],
					       1);

	  const LevelData<FArrayBox>& coarseH = *vectH[lev-1]; 
	  ghostFiller.fillInterp(levelH, coarseH, coarseH, 0.0, 0, 0, 1);
	  
	  Real dx = m_amrDx[lev];
	  
	  QuadCFInterp qcfi(levelGrids, &m_amrGrids[lev-1], dx, 
			    m_refinement_ratios[lev-1],1, m_amrDomains[lev]);
	  qcfi.coarseFineInterp(levelH, *vectH[lev-1]);

	}
     
      levelH.exchange();

    }

  
  for (int lev=0; lev <= finestTimestepLevel() ; ++lev)
    {
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFluxes = *a_vectFluxes[lev];
      const LevelData<FArrayBox>& levelH = *vectH[lev];
      const LevelData<FluxBox>& levelD = *m_diffusivity[lev];
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      const RealVect& dx = levelCoords.dx();
      
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
       	  for (int dir = 0; dir < SpaceDim; ++dir)
	    {
	      FArrayBox& faceFlux = levelFluxes[dit][dir];
	      const FArrayBox& faceD = levelD[dit][dir];
	      const FArrayBox& cellH = levelH[dit];
	      Box faceBox = levelGrids[dit];
	      faceBox.surroundingNodes(dir);
	      FArrayBox gradH(faceBox,1);
	      Real oneOnDx = 1.0 / dx[dir];
	      for (BoxIterator bit(faceBox); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		  gradH(iv) =  oneOnDx * 
		    ( cellH(iv) - cellH(iv-BASISV(dir)));  
		}
	      gradH *= faceD;
	      faceFlux += gradH;
	      	      
	      //FORT_SUBTRACTDFLUX(CHF_FRA1(faceFlux,0),
	      //			 CHF_CONST_FRA1(cellH,0),
	      //			 CHF_CONST_FRA1(faceD,0),
	      //			 CHF_CONST_INT(dir),
	      //			 CHF_CONST_REAL(dx[dir]),
	      //			 CHF_BOX(faceBox));
	    }
	}
    }

}
/// given the current cell centred velocity field, compute a face centred velocity field
void 
AmrIce::computeFaceVelocity(Vector<LevelData<FluxBox>* >& a_faceVel, 
			    Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			    Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel) 
{
  for (int lev = 0; lev <= m_finest_level; lev++)
    {
      LevelData<FArrayBox>* crseVelPtr = (lev > 0)?m_velocity[lev-1]:NULL;
      int nRefCrse = (lev > 0)?m_refinement_ratios[lev-1]:1;
      IceVelocity::computeFaceVelocity
	(*a_faceVel[lev], *a_layerXYFaceXYVel[lev], *a_layerSFaceXYVel[lev],
	 *m_velocity[lev],  *m_vect_coordSys[lev],  
	 *m_A[lev], *m_sA[lev], *m_bA[lev], crseVelPtr, nRefCrse,
	 m_constitutiveRelation, m_additionalVelocity);
    }
}


/// compute div(vel*H) at a given time
void
AmrIce::computeDivThicknessFlux(Vector<LevelData<FArrayBox>* >& a_divFlux,
                                Vector<LevelData<FluxBox>* >& a_flux,
                                Vector<LevelData<FArrayBox>* >& a_thickness,
                                Real a_time, Real a_dt)
{

  //Vector<LevelData<FluxBox>* > faceVel(m_finest_level+1, NULL);
  
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      // construct face-centered velocity field
      LevelData<FluxBox> faceVel(levelGrids, 1, IntVect::Unit);
                                 
      if (lev > 0)
        {
          int nVelGhost = m_velocity[lev]->ghostVect()[0];
          
          PiecewiseLinearFillPatch velFiller(levelGrids, 
                                             m_amrGrids[lev-1],
                                             m_velocity[0]->nComp(), 
                                             m_amrDomains[lev-1],
                                             m_refinement_ratios[lev-1],
                                             nVelGhost);
          
          // since we're not subcycling, don't need to interpolate in time
          Real time_interp_coeff = 0.0;
          velFiller.fillInterp(*m_velocity[lev],
                               *m_velocity[lev-1],
                               *m_velocity[lev-1],
                               time_interp_coeff,
                               0, 0, m_velocity[0]->nComp());
          
        } // end if lev > 0
      m_velocity[lev]->exchange();

      // average velocities to faces
      CellToEdge(*m_velocity[lev],faceVel);
      faceVel.exchange();
      
      // flux = faceVel*faceH      
      LevelData<FluxBox>& levelFlux = *a_flux[lev];
      LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      LevelData<FluxBox>& faceH = levelCoords.getFaceH();

      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& thisflux = levelFlux[dit];
          FluxBox& thisVel = faceVel[dit];
          FluxBox& thisH = faceH[dit];

          for (int dir=0; dir<SpaceDim; dir++)
            {
              thisflux[dir].copy(thisVel[dir]);
              thisflux[dir].mult(thisH[dir]);
            }
        }

      // average fluxes to coarser levels, if needed
      if (lev>0)
        {
          CoarseAverageFace faceAverager(m_amrGrids[lev],
                                         1, m_refinement_ratios[lev-1]);
          faceAverager.averageToCoarse(*a_flux[lev-1], *a_flux[lev]);
          
        }
      
    } // end loop over levels

  // now compute div(flux)
  
  // compute div(F) and add source term
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFlux = *a_flux[lev];
      LevelData<FArrayBox>& levelDiv = *a_divFlux[lev];
      LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);

      LevelData<FArrayBox>& surfaceThicknessSource = *m_surfaceThicknessSource[lev];
      m_surfaceFluxPtr->surfaceThicknessFlux(surfaceThicknessSource, *this, lev, a_dt);
      
      LevelData<FArrayBox>& basalThicknessSource = *m_basalThicknessSource[lev];
      m_basalFluxPtr->surfaceThicknessFlux(basalThicknessSource, *this, lev, a_dt);
      m_calvingModelPtr->modifySurfaceThicknessFlux(basalThicknessSource, *this, lev, a_dt);

      const RealVect& dx = levelCoords.dx();          

      DataIterator dit = levelGrids.dataIterator();
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          FArrayBox& thisDiv = levelDiv[dit];
          
          FluxBox& thisFlux = levelFlux[dit];
          thisDiv.setVal(0.0);
          
          // loop over directions and increment with div(F)
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // use the divergence from 
              // Chombo/example/fourthOrderMappedGrids/util/DivergenceF.ChF
              FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                              CHF_FRA(thisDiv),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(dx[dir]),
                              CHF_INT(dir));
            }

          // add in thickness source here
          thisDiv.minus(surfaceThicknessSource[dit], gridBox,0,0,1);
	  thisDiv.minus(basalThicknessSource[dit], gridBox,0,0,1);
          //thisDiv *= -1*a_dt;
        } // end loop over grids
    } // end loop over levels
  
}

// increment phi := phi + dt*dphi
void
AmrIce::incrementWithDivFlux(Vector<LevelData<FArrayBox>* >& a_phi,
                             const Vector<LevelData<FArrayBox>* >& a_dphi,
                             Real a_dt)
{
  for (int lev=0; lev<a_phi.size(); lev++)
    {
      LevelData<FArrayBox>& levelPhi = *a_phi[lev];
      const LevelData<FArrayBox>& level_dPhi = *a_dphi[lev];

      DataIterator dit = levelPhi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          levelPhi[dit].plus(level_dPhi[dit], a_dt);
        }
    }
}
 

// increment coordSys with new thickness
void
AmrIce::updateCoordSysWithNewThickness(const Vector<LevelData<FArrayBox>* >& a_thickness)
{
  CH_assert(a_thickness.size() >= m_finest_level);
  
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      const LevelData<FArrayBox>& levelH = *a_thickness[lev];
      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      LevelData<FArrayBox>& levelCS_H = levelCS.getH();
      DataIterator dit = levelH.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisH = levelCS_H[dit];
          thisH.copy(levelH[dit]);
        }
      {
	LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
	int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
	levelCS.recomputeGeometry(crseCoords, refRatio);
      }
    } // end loop over levels      
}


// compute timestep
Real 
AmrIce::computeDt()
{
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::computeDt" << endl;
    }

  if (m_fixed_dt > TINY_NORM)
    return m_fixed_dt;

  Real dt = 1.0e50;
  for (int lev=0; lev<= finestTimestepLevel(); lev++)
    {
      // Real dtLev = dt;
      // const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      // //slc : this doesn't include any L1L2 contribution for now.
      // LevelData<FArrayBox>& levelVel = *m_velocity[lev];
      // DataIterator levelDit = levelVel.dataIterator();
      // for (levelDit.reset(); levelDit.ok(); ++levelDit)
      //   {
      //     int p = 0;
      //     Real maxVel = levelVel[levelDit].norm(levelGrids[levelDit],
      //                                           p, 0, 2);
      //     Real localDt = m_amrDx[lev]/maxVel;
      //     dtLev = min(dtLev, localDt);
      //   }
      Real dtLev = dt;
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      //since we now compute m_faceVel immediately after m_velocity, use that
      const LevelData<FluxBox>& levelVel = *m_faceVel[lev]; 
      DataIterator levelDit = levelVel.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
	for (int dir = 0; dir < SpaceDim; dir++)
	  {
	    int p = 0;
	    Box faceBox = levelGrids[levelDit];
	    faceBox.surroundingNodes(dir);
	    Real maxVel = 1.0 + levelVel[levelDit][dir].norm(faceBox,p, 0, 1);
	    Real localDt = m_amrDx[lev]/maxVel;
	    dtLev = min(dtLev, localDt);
	  }
      }
      
      if (m_diffusionTreatment == EXPLICIT){
	MayDay::Error("diffusion_treatment == explicit not supported now : use none");
      }
      dt = min(dt, dtLev);
    }

#ifdef CH_MPI
      Real tmp = 1.;
      int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL,
                                 MPI_MIN, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          MayDay::Error("communication error on norm");
        }
      dt = tmp;
#endif

  if (m_cur_step == 0)
    {
      dt *= m_initial_cfl;
    } 
  else 
    {
      dt *= m_cfl;
    }

  // also check to see if max grow rate applies
  // (m_dt > 0 test screens out initial time, when we set m_dt to a negative 
  // number by default)
  if (dt > m_max_dt_grow*m_dt && (m_dt > 0) )
    dt = m_max_dt_grow*m_dt;
  
  if (m_timeStepTicks){
    // reduce time step to integer power of two
    dt = std::pow(2.0, std::floor(std::log(dt)/std::log(two)));
    
  }
  
  

  return dt;// min(dt,2.0);

}

Real 
AmrIce::computeInitialDt()
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::computeInitialDt" << endl;
    }


  // for now, just call computeDt;
  Real dt = computeDt();
  return dt;
}



//determine the grouding line proximity 
/**

   Solves the elliptic problem 
      a * phi - b* grad^2 phi = 0;
   with natural boundary conditions.

   for grounded ice, a = 10^5 and b = 1
   for floating ice, s = 0 and b = 1
 */
void AmrIce::updateGroundingLineProximity() const
{

  CH_TIME("AmrIce::updateGroundingLineProximity");

  if (m_groundingLineProximity_valid)
    return;

  if (m_groundingLineProximity.size() < m_finest_level + 1)
    {
      m_groundingLineProximity.resize(m_finest_level + 1, NULL);
    }

  //Natural boundary conditions
  BCHolder bc(ConstDiriNeumBC(IntVect(0,0), RealVect(0.0,0.0),
  			      IntVect(0,0), RealVect(0.0,0.0)));

  //BCHolder bc(ConstDiriNeumBC(IntVect(0,0), RealVect(-1.0,-1.0),
  //			      IntVect(0,0), RealVect(1.0,1.0)));

  Vector<RefCountedPtr<LevelData<FArrayBox> > > a(m_finest_level + 1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > b(m_finest_level + 1);
  Vector<LevelData<FArrayBox>* > rhs(m_finest_level+ 1,NULL);
  Vector<DisjointBoxLayout> grids(finestTimestepLevel() + 1);
  Vector<ProblemDomain> domains(finestTimestepLevel() + 1);
  Vector<RealVect> dx(finestTimestepLevel() + 1);

  for (int lev=0; lev <= m_finest_level; ++lev)
    {
      dx[lev] = m_amrDx[lev]*RealVect::Unit;
      domains[lev] = m_amrDomains[lev];

      const LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const LevelData<BaseFab<int> >& levelMask = levelCS.getFloatingMask();
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      a[lev] = RefCountedPtr<LevelData<FArrayBox> >
 	(new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero));
      b[lev] = RefCountedPtr<LevelData<FluxBox> >
 	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Zero));
      rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
      
      grids[lev] = levelGrids;

     

      if (m_groundingLineProximity[lev] != NULL)
	{
	  delete m_groundingLineProximity[lev];
	  m_groundingLineProximity[lev] = NULL;
	}
      m_groundingLineProximity[lev] =  new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
      
      LevelData<FArrayBox>& levelPhi = *m_groundingLineProximity[lev];

      const Real& crseDx = m_amrDx[0];
      Real crseDxSq = crseDx*crseDx;

      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
 	{
	  FluxBox& B = (*b[lev])[dit];
	  
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      B[dir].setVal(crseDxSq);
	    }

 	  FArrayBox& r =  (*rhs[lev])[dit];
 	  r.setVal(0.0);
	  FArrayBox& A =  (*a[lev])[dit];
	  A.setVal(0.0);
	  FArrayBox& phi = levelPhi[dit];
	  phi.setVal(0.0);

 	  const BaseFab<int>& mask = levelMask[dit];
 	  const Box& gridBox = levelGrids[dit];
	  //	  const FArrayBox& u = (*m_velocity[lev])[dit];
 	  for (BoxIterator bit(gridBox);bit.ok();++bit)
 	    {
 	      const IntVect& iv = bit();
 	      if (mask(iv) == GROUNDEDMASKVAL )
 		{
 		  A(iv) = 1.0e+0;
		  r(iv) = 1.0e+0;
		  
 		} 
	      else
		{
		  A(iv) = crseDx / m_groundingLineProximityScale;
		}
	      
 	    }
 
 	}

      rhs[lev]->exchange();
      m_groundingLineProximity[lev]->exchange();
      a[lev]->exchange();
      b[lev]->exchange();
    }


  VCAMRPoissonOp2Factory* poissonOpFactory = new VCAMRPoissonOp2Factory;
  poissonOpFactory->define(domains[0], grids , m_refinement_ratios,
 			   m_amrDx[0], bc, 1.0, a,  1.0 , b);
  RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
    opFactoryPtr(poissonOpFactory);

  MultilevelLinearOp<FArrayBox> poissonOp;
  poissonOp.define(grids, m_refinement_ratios, domains, dx, opFactoryPtr, 0);
    
  RelaxSolver<Vector<LevelData<FArrayBox>* > >* relaxSolver
    = new RelaxSolver<Vector<LevelData<FArrayBox>* > >();

  relaxSolver->define(&poissonOp,false);
  relaxSolver->m_verbosity = s_verbosity;
  relaxSolver->m_normType = 0;
  relaxSolver->m_eps = 1.0e-8;
  relaxSolver->m_imax = 12;
  relaxSolver->m_hang = 0.05;
  relaxSolver->solve(m_groundingLineProximity,rhs);

  delete(relaxSolver);

#ifdef DUMP_PROXIMITY
  std::string file("proximity.2d.hdf5");
  Real dt = 0.0; 
  Real time = 0.0;
  Vector<std::string> names(1,"proximity");
  WriteAMRHierarchyHDF5(file ,grids, m_groundingLineProximity ,names, m_amrDomains[0].domainBox(),
  			m_amrDx[0], dt, m_time, m_refinement_ratios, m_groundingLineProximity.size());
#endif
  
  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      if (rhs[lev] != NULL)
 	{
 	  delete rhs[lev];
	  rhs[lev] = NULL;
 	}
    }

  m_groundingLineProximity_valid = true;
}

//access the viscous tensor (cell-centered)
const LevelData<FArrayBox>* AmrIce::viscousTensor(int a_level) const
{
  updateViscousTensor();
  if (!(m_viscousTensorCell.size() > a_level))
    {
      std::string msg("AmrIce::viscousTensor !(m_viscousTensorCell.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_viscousTensorCell.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_viscousTensorCell[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::viscousTensor m_viscousTensorCell[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;

}

//access the viscous tensor (cell-centered)
const LevelData<FArrayBox>* AmrIce::viscosityCoefficient(int a_level) const
{
  updateViscousTensor();
  if (!(m_viscosityCoefCell.size() > a_level))
    {
      std::string msg("AmrIce::viscosityCoef !(m_viscosityCoefCell.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_viscosityCoefCell.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_viscosityCoefCell[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::viscosityCoef m_viscosityCoefCell[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;

}

//access the drag coefficient (cell-centered)
const LevelData<FArrayBox>* AmrIce::dragCoefficient(int a_level) const
{
  updateViscousTensor();
  if (!(m_dragCoef.size() > a_level))
    {
      std::string msg("AmrIce::dragCoef !(m_dragCoef.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_dragCoef.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_dragCoef[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::dragCoef m_dragCoefCell[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;
}



//update the viscous tensor components
void AmrIce::updateViscousTensor() const
{
  CH_TIME("AmrIce::updateViscousTensor");

  if (m_viscousTensor_valid)
    return;
  
  if (m_viscousTensorCell.size() < m_finest_level + 1)
    {
      m_viscousTensorCell.resize(m_finest_level + 1, NULL);
    }
  if (m_viscosityCoefCell.size() < m_finest_level + 1)
    {
      m_viscosityCoefCell.resize(m_finest_level + 1, NULL);
    }
  if (m_dragCoef.size() < m_finest_level + 1)
    {
      m_dragCoef.resize(m_finest_level + 1, NULL);
    }

  if (m_viscousTensorFace.size() < m_finest_level + 1)
    {
      m_viscousTensorFace.resize(m_finest_level + 1, NULL);
    }

 
  Vector<LevelData<FluxBox>*> faceA(m_finest_level + 1,NULL);
  Vector<RefCountedPtr<LevelData<FluxBox> > > viscosityCoef;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > dragCoef;
  Vector<LevelData<FArrayBox>* > C0(m_finest_level + 1,  NULL);

  Vector<RealVect> vdx(m_finest_level + 1);
  for (int lev =0; lev <= m_finest_level; lev++)
    {
      faceA[lev] = new LevelData<FluxBox>(m_amrGrids[lev],m_A[lev]->nComp(),IntVect::Unit);
      CellToEdge(*m_A[lev],*faceA[lev]);

      if (m_viscousTensorFace[lev] != NULL)
	{
	  delete m_viscousTensorFace[lev];m_viscousTensorFace[lev]=NULL;
	}
      m_viscousTensorFace[lev] = new LevelData<FluxBox>(m_amrGrids[lev],SpaceDim,IntVect::Unit);

      if (m_viscousTensorCell[lev] != NULL)
	{
	  delete m_viscousTensorCell[lev];m_viscousTensorCell[lev]=NULL;
	}
      m_viscousTensorCell[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim*SpaceDim,IntVect::Zero);
      
      if (m_dragCoef[lev] != NULL)
	{
	  delete m_dragCoef[lev]; m_dragCoef[lev] = NULL;
	}
      m_dragCoef[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim,IntVect::Zero);
      
      if (m_viscosityCoefCell[lev] != NULL)
	{
	  delete m_viscosityCoefCell[lev]; m_viscosityCoefCell[lev] = NULL;
	}
      m_viscosityCoefCell[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim,IntVect::Zero);

      if (C0[lev] != NULL)
	{
	  delete C0[lev];C0[lev] = NULL;
	}
      C0[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1,m_velBasalC[0]->ghostVect());
      DataIterator dit = m_amrGrids[lev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*C0[lev])[dit].setVal(0.0);
        }
      vdx[lev] = RealVect::Unit*m_amrDx[lev];
    }

  //these parameters don't matter because we don't solve anything here. 
  Real vtopSafety = 1.0;
  int vtopRelaxMinIter = 4;
  Real vtopRelaxTol = 1.0;
  Real muMin = 0.0; 
  Real muMax = 1.23456789e+300;

  int numLevels = m_finest_level + 1;
  IceJFNKstate state(m_amrGrids, m_refinement_ratios, m_amrDomains, vdx, m_vect_coordSys, 
		     m_velocity, m_velBasalC, C0, numLevels-1, 
		     *m_constitutiveRelation,  *m_basalFrictionRelation, *m_thicknessIBCPtr,  
		     m_A, faceA, m_time, vtopSafety, vtopRelaxMinIter, vtopRelaxTol, 
		     muMin, muMax);
  state.setState(m_velocity);
  viscosityCoef = state.getViscosityCoef();
  dragCoef = state.getDragCoef();
  state.computeViscousTensorFace(m_viscousTensorFace);
  
  for (int lev =0; lev < numLevels; lev++)
    {
      
      //If a cell is adjacent to a calving front, we  set the (vertically integrated)
      //viscous tensor components at the intervening face to zero. That works well enough for the velocity solves,
      //but causes pain here because the cell average (in EdgeToCell) will end up half the value at the other face.
      for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
      	{

      	  const FArrayBox& thck = m_vect_coordSys[lev]->getH()[dit];
      	  const FArrayBox& dsdx = m_vect_coordSys[lev]->getGradSurface()[dit];
	  const FArrayBox& usrf = m_vect_coordSys[lev]->getSurfaceHeight()[dit];
      	  const BaseFab<int>& mask = m_vect_coordSys[lev]->getFloatingMask()[dit];
      	  const Real& rhoi = m_vect_coordSys[lev]->iceDensity();
      	  const Real& rhoo = m_vect_coordSys[lev]->waterDensity();
      	  const Real& gravity = m_vect_coordSys[lev]->gravity();
      	  const Real rgr = rhoi * gravity * (1.0-rhoi/rhoo);
      	  const RealVect& dx = m_vect_coordSys[lev]->dx();

      	  for (int dir = 0; dir < SpaceDim; dir++)
      	    {
      	      FArrayBox& facevt = (*m_viscousTensorFace[lev])[dit][dir];
      	      Real factor = dx[dir] * rgr;
      	      FORT_SETFRONTFACEVT(CHF_FRA1(facevt,dir),
      				  CHF_CONST_FRA1(thck,0),
      				  CHF_CONST_FRA1(usrf,0),
      				  CHF_CONST_FIA1(mask,0),
      				  CHF_CONST_INT(dir),
      				  CHF_CONST_REAL(factor),
      				  CHF_BOX(m_amrGrids[lev][dit]));
      	    }
      	}


      EdgeToCell(*m_viscousTensorFace[lev],*m_viscousTensorCell[lev]);
      EdgeToCell(*viscosityCoef[lev],*m_viscosityCoefCell[lev]);

      for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
      	{
	  const BaseFab<int>& mask = m_vect_coordSys[lev]->getFloatingMask()[dit];
	  FArrayBox& cellvt = (*m_viscousTensorCell[lev])[dit];
	  const Real z = 0.0;
	  for (int comp = 0; comp < SpaceDim * SpaceDim; comp++)
	    {
	      FORT_SETICEFREEVAL(CHF_FRA1(cellvt,comp), 
				 CHF_CONST_FIA1(mask,0),
				 CHF_CONST_REAL(z),
				 CHF_BOX(m_amrGrids[lev][dit]));
	    }
	}

      dragCoef[lev]->copyTo(Interval(0,0),*m_dragCoef[lev],Interval(0,0));

      if (faceA[lev] != NULL)
	{
	  delete faceA[lev]; faceA[lev] = NULL;
	}
      if (C0[lev] != NULL)
	{
	  delete C0[lev]; C0[lev] = NULL;
	}
    }

  m_viscousTensor_valid = true;

}

//access the grounding line proximity
const LevelData<FArrayBox>* AmrIce::groundingLineProximity(int a_level) const
{

  updateGroundingLineProximity();
  
  if (!(m_groundingLineProximity.size() > a_level))
    {
      std::string msg("AmrIce::groundingLineProximity !(m_groundingLineProximity.size() > a_level)");
      pout() << msg << endl;
      CH_assert((m_groundingLineProximity.size() > a_level));
      MayDay::Error(msg.c_str());
    }


  LevelData<FArrayBox>* ptr = m_groundingLineProximity[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::groundingLineProximity m_groundingLineProximity[a_level] == NULL)");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;
}



///Identify regions of floating ice that are remote
///from grounded ice and eliminate them.
/**
   Regions of floating ice unconnected to land 
   lead to an ill-posed problem, with a zero
   basal traction coefficient and Neumann boundary
   conditions. Here, we attempt to identify them
   by solving - div(A grad(phi)) + C phi = R with 
   natural (or periodic) BCs

   for ice, A = H * muCooef, elsewhere A = 0

   for grounded ice, C = 1 and R = 1
   for floating ice, C = 0 and R = 0
   elsewhere, C = 1 and R = -1

   The solution for the well-posed regions is
   phi = 1, and tensd to be -1 in the ill-posed
   regions.

 */ 
void AmrIce::eliminateRemoteIce()
{
  
  //Natural boundary conditions
  BCHolder bc(ConstDiriNeumBC(IntVect(0,0), RealVect(0.0,0.0),
 			      IntVect(0,0), RealVect(0.0,0.0)));
  
  Vector<RefCountedPtr<LevelData<FArrayBox> > > C(m_finest_level + 1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > D(m_finest_level + 1);
  Vector<LevelData<FArrayBox>* > rhs(m_finest_level+ 1,NULL);
  Vector<LevelData<FArrayBox>* > phi(m_finest_level + 1);
  Vector<DisjointBoxLayout> grids(finestTimestepLevel() + 1);
  Vector<ProblemDomain> domains(finestTimestepLevel() + 1);
  Vector<RealVect> dx(finestTimestepLevel() + 1);
  for (int lev=0; lev <= m_finest_level; ++lev)
    {
      dx[lev] = m_amrDx[lev]*RealVect::Unit;
      domains[lev] = m_amrDomains[lev];

      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const LevelData<BaseFab<int> >& levelMask = levelCS.getFloatingMask();
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      C[lev] = RefCountedPtr<LevelData<FArrayBox> >
 	(new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero));
      D[lev] = RefCountedPtr<LevelData<FluxBox> >
 	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Zero));
      rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
      phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
      grids[lev] = levelGrids;
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
 	{
	  FluxBox& d= (*D[lev])[dit];
	  
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      d[dir].copy(levelCS.getFaceH()[dit][dir]);
	      d[dir] *= (*m_faceMuCoef[lev])[dit][dir];
	      d[dir] += 1.0e+2;
	    }

 	  FArrayBox& r =  (*rhs[lev])[dit];
 	  r.setVal(0.0);

	  FArrayBox& p =  (*phi[lev])[dit];
 	  p.setVal(1.0);

	  FArrayBox& c =  (*C[lev])[dit];
	  c.setVal(0.0);

 	  const BaseFab<int>& mask = levelMask[dit];
 	  const Box& gridBox = levelGrids[dit];
 	  for (BoxIterator bit(gridBox);bit.ok();++bit)
 	    {
 	      const IntVect& iv = bit();
 	      if (mask(iv) == GROUNDEDMASKVAL)
 		{
 		  c(iv) = 1.0;
		  p(iv) = 1.0;
		  r(iv) = 1.0;
 		} 
	      else if (mask(iv) == OPENSEAMASKVAL || mask(iv) == OPENLANDMASKVAL)
 		{
 		  c(iv) = 1.0;
		  p(iv) = -1.0;
		  r(iv) = -1.0;
 		} 
	      else
		{
		  c(iv) = 0.0;
		  p(iv) = 0.0;
		  r(iv) = 0.0;
		}
	      
 	    }
 
 	}

      rhs[lev]->exchange();
      phi[lev]->exchange();
      D[lev]->exchange();
      C[lev]->exchange();
    }


  VCAMRPoissonOp2Factory* poissonOpFactory = new VCAMRPoissonOp2Factory;
  poissonOpFactory->define(domains[0], grids , m_refinement_ratios,
 			   m_amrDx[0], bc, 1.0, C,  1.0 , D);
  RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
    opFactoryPtr(poissonOpFactory);

  MultilevelLinearOp<FArrayBox> poissonOp;
  poissonOp.define(grids, m_refinement_ratios, domains, dx, opFactoryPtr, 0);
    
  RelaxSolver<Vector<LevelData<FArrayBox>* > >* relaxSolver
    = new RelaxSolver<Vector<LevelData<FArrayBox>* > >();

  relaxSolver->define(&poissonOp,false);
  relaxSolver->m_verbosity = s_verbosity;
  relaxSolver->m_normType = 0;
  relaxSolver->m_eps = 1.0e-10;
  relaxSolver->m_imax = 4;
  relaxSolver->m_hang = 0.05;
  relaxSolver->solve(phi,rhs);


  for (int lev = m_finest_level; lev > 0 ; --lev)
    {
      CoarseAverage averager(grids[lev],1,m_refinement_ratios[lev-1]);
      averager.averageToCoarse(*phi[lev-1],*phi[lev]);
    }

  //std::string file("connectivity.2d.hdf5");
  //Real dt = 0.0; //Real time = 0.0;
  //Vector<std::string> names(1,"conn");
  //WriteAMRHierarchyHDF5(file ,grids, phi ,names, m_amrDomains[0].domainBox(),
  //			m_amrDx[0], dt, m_time, m_refinement_ratios, phi.size());
  //MayDay::Error("stopping after connectivity.2d.hdf5");


  Real threshold = -0.99;

  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      const LevelData<FArrayBox>& levelPhi = *phi[lev];
      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
     
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const FArrayBox& thisPhi = levelPhi[dit];
	  FArrayBox& thisH = levelCS.getH()[dit];
	  FArrayBox& thisU = (*m_velocity[lev])[dit];
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	  for (BoxIterator bit(levelGrids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if (mask(iv) == FLOATINGMASKVAL && thisPhi(iv) < threshold)
	      	{
		  thisH(iv) = 0.0;
	      	}

	    }
	}

      LevelSigmaCS* crseCS = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
      int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
      levelCS.recomputeGeometry(crseCS, refRatio);     
    }



  delete(relaxSolver);

  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      if (rhs[lev] != NULL)
 	{
 	  delete rhs[lev];
	  rhs[lev] = NULL;
 	}
      if (phi[lev] != NULL)
 	{
 	  delete phi[lev];
 	  phi[lev] = NULL;
 	}
    }


}


void 
AmrIce::implicitThicknessCorrection(Real a_dt,
				    const Vector<LevelData<FArrayBox>* >& a_sts,
				    const Vector<LevelData<FArrayBox>* >& a_bts
				    )
{

  
 
  if  (m_temporalAccuracy == 1)
    {  
      //implicit Euler : solve (I - dt P) H = H_pred + dt * S
      
      //slc: at the moment, I'm setting eveything up every time-step,
      //pretending that diffusion is constant in time, and using the multi-grid
      //solver only. All these things are to be improved 
      
      //VCAMRPoissonOp2Factory opf;
      //RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > opfPtr 
      //	= RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > >(new VCAMRPoissonOp2Factory());
      //VCAMRPoissonOp2Factory& vcpFactory = static_cast<VCAMRPoissonOp2Factory> (*opfPtr);

      //BCHolder bc(doNothingBC); 
      // strictly, periodic problems only for now
      // although in many non-periodic cases there will no boundary
      // grad(H) and so no diffusive flux

      //Natural boundary conditions - OK for now, but ought to get 
      //moved into subclasses of IceThicknessIBC
      BCHolder bc(ConstDiriNeumBC(IntVect(0,0), RealVect(0.0,0.0),
      				  IntVect(0,0), RealVect(0.0,0.0)));

      Vector<RefCountedPtr<LevelData<FArrayBox> > > I(finestTimestepLevel() + 1);
      Vector<RefCountedPtr<LevelData<FluxBox> > > D(finestTimestepLevel() + 1);
      Vector<LevelData<FArrayBox>* > H(finestTimestepLevel() + 1);
      Vector<LevelData<FArrayBox>* > rhs(finestTimestepLevel()+ 1);
      Vector<DisjointBoxLayout> grids(finestTimestepLevel() + 1);

      for (int lev=0; lev <= finestTimestepLevel(); ++lev)
	{
	  LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
	  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

	  I[lev] = RefCountedPtr<LevelData<FArrayBox> >
	    (new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit));

	  //	  H[lev] = &levelCoords.getH();
	  H[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);

	  D[lev] = m_diffusivity[lev];

	  rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);

	  grids[lev] = levelGrids;

	  const LevelData<FArrayBox>& levelSTS = *a_sts[lev];
	  const LevelData<FArrayBox>& levelBTS = *a_bts[lev];
	  

	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      
	      (*I[lev])[dit].setVal(one);
	      (*H[lev])[dit].copy(levelCoords.getH()[dit] , 0 , 0, 1);
	      (*rhs[lev])[dit].copy( (*H[lev])[dit] , 0 , 0, 1);
	      (*rhs[lev])[dit].plus(levelSTS[dit],a_dt);
	      (*rhs[lev])[dit].plus(levelBTS[dit],a_dt); 
	      (*D[lev])[dit][0].plus(m_additionalDiffusivity);
	      (*D[lev])[dit][1].plus(m_additionalDiffusivity);
	      
	    }
	  rhs[lev]->exchange();
	  H[lev]->exchange();
	  m_diffusivity[lev]->exchange();
	  I[lev]->exchange();
	}

      VCAMRPoissonOp2Factory* poissonOpFactory = new VCAMRPoissonOp2Factory;
      poissonOpFactory->define(m_amrDomains[0], grids , m_refinement_ratios,
      		       m_amrDx[0], bc, 1.0, I,  a_dt, D);
      for (int lev=0; lev <= finestTimestepLevel(); ++lev)
	{
	  H[lev]->exchange();
	}
      
#if 0   
      //BiCGstab
      RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory
	= RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > >(poissonOpFactory);
      BiCGStabSolver<Vector<LevelData<FArrayBox>* > > krylovSolver;
      MultilevelLinearOp<FArrayBox> mlOp;
      Vector<RealVect> dx(m_amrDx.size());
      for (int i = 0; i < dx.size(); ++i)
	{
	  for (int dir =0; dir < SpaceDim; ++dir)
	    {
	      dx[i][dir] = m_amrDx[i];
	    }
	}
      mlOp.m_num_mg_smooth = 8 + 4 * m_finest_level ;
      mlOp.define(grids,m_refinement_ratios,m_amrDomains,dx,opFactory,0);
      krylovSolver.define(&mlOp , false);
      krylovSolver.m_verbosity = s_verbosity - 1;
      krylovSolver.m_reps = 1.0e-10;
      krylovSolver.m_imax = 50;
      krylovSolver.m_normType = 2;
      krylovSolver.solve(H, rhs);

#else
      //Plain MG
      BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
      AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
      mgSolver.define(m_amrDomains[0], *poissonOpFactory , &bottomSolver, finestTimestepLevel()+1);
      //parse these
      mgSolver.m_eps = 1.0e-10;
      mgSolver.m_normThresh = 1.0e-10;
    
      int numMGSmooth = 8 + 16 * m_finest_level ;
      mgSolver.m_pre = numMGSmooth;
      mgSolver.m_post = numMGSmooth;
      mgSolver.m_bottom = numMGSmooth;
      
      mgSolver.solve(H, rhs, finestTimestepLevel(), 0,  false);
      // for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
      // 	{
      // 	  RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
      // 	  m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
      // 					    m_amrDomains[lev],levelDx, m_time, m_dt);
      // 	}
      // mgSolver.m_eps = 1.0e-10;
#endif

      for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
	{
	  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	  LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
          LevelData<FArrayBox>& levelCoord_H = levelCoords.getH();
	  
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      Real maxNorm = (*H[lev])[dit].norm(0,0,1);
	      CH_assert(maxNorm < HUGE_THICKNESS);
	      //CH_assert((*H[lev])[dit].min() > 1.0);
	      levelCoord_H[dit].copy( (*H[lev])[dit], 0, 0, 1);

	      //put sensible values into the corners.
	      FArrayBox &thisH = levelCoord_H[dit];
	      Box sbox = thisH.box();
	      sbox.grow(-levelCoord_H.ghostVect()[0]);
	      FORT_EXTRAPCORNER2D(CHF_FRA(thisH),
      				  CHF_BOX(sbox));

	    }

	  delete rhs[lev];
	  delete H[lev];
	}
    }
  else 
    {    
      MayDay::Error("AmrIce::implicitThicknessCorrection, invalid temporal accuracy");
    }


  

}



#ifdef CH_USE_HDF5

  /// write hdf5 plotfile to the standard location
void 
AmrIce::writePlotFile() 
{
   
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::writePlotFile" << endl;
    }
  
  // plot comps: thickness + horizontal velocity + zb + zs
  // slc: + base velocity + 'flux' velocity
  int numPlotComps = 6;



  // may need a zvel for visit to do "3d" streamlines correctly
  bool writeZvel = true;
  if (writeZvel) numPlotComps+=1;

  if (m_write_fluxVel)
    {
      numPlotComps += 2;
      if (writeZvel) numPlotComps+=1;
    }
  
 
  if (m_write_baseVel)
    {
      numPlotComps += 2;
      if (writeZvel) numPlotComps+=1;
    }

  if (m_write_dHDt) numPlotComps += 1;
  if (m_write_solver_rhs) numPlotComps += (SpaceDim+2);

  if (m_write_temperature) 
    numPlotComps += m_temperature[0]->nComp();
#if BISICLES_Z == BISICLES_LAYERED
  if (m_write_temperature)
    numPlotComps += 2;// surface and basal temperatures

  //layer velocities
  if (m_write_layer_velocities)
    {
      numPlotComps += SpaceDim * (m_nLayers+1);
      if (writeZvel) numPlotComps += (m_nLayers+1);
    }

  if (m_write_viscousTensor)
    {
      numPlotComps += 1 + 1 + SpaceDim * SpaceDim;
    }
  
  if (m_write_thickness_sources)
    {
      numPlotComps += 3;  // surface and basal sources, plus the balance
    }


#endif
  // generate data names

  string thicknessName("thickness");
  string xVelName("xVel");
  string yVelName("yVel");
  string zVelName("zVel");
  string zsName("Z_surface");
  string zbName("Z_base");
  string zbottomName("Z_bottom");
  string dthicknessName("dThickness/dt");
  string betaName("basal_friction");
  string solverRhsxName("solverRHSx");
  string solverRhsyName("solverRHSy");
  string C0Name("C0");

  string xfVelName("xfVel");
  string yfVelName("yfVel");
  string zfVelName("zfVel");
  string xbVelName("xbVel");
  string ybVelName("ybVel");
  string zbVelName("zbVel");

  string temperatureName("temperature");
#if BISICLES_Z == BISICLES_LAYERED
  string xlayerVelName("xlayerVel");
  string ylayerVelName("ylayerVel");
  string zlayerVelName("zlayerVel");
#endif

  string xxVTname("xxViscousTensor");
  string xyVTname("xyViscousTensor");
  string xzVTname("xzViscousTensor");
  string yxVTname("yxViscousTensor");
  string yyVTname("yyViscousTensor");
  string yzVTname("yzViscousTensor");
  string zxVTname("zxViscousTensor");
  string zyVTname("zyViscousTensor");
  string zzVTname("zzViscousTensor");
  string viscosityCoefName("viscosityCoef");
  //string yViscosityCoefName("yViscosityCoef");
  //string zViscosityCoefName("zViscosityCoef");
  string dragCoefName("dragCoef");

  string basalThicknessSourceName("basalThicknessSource");
  string surfaceThicknessSourceName("surfaceThicknessSource");
  string surfaceThicknessBalanceName("surfaceThicknessBalance");

  Vector<string> vectName(numPlotComps);
  //int dThicknessComp;

  vectName[0] = thicknessName;
  vectName[1] = xVelName;
  vectName[2] = yVelName;
  int comp = 3;
  if (writeZvel) 
    {
      vectName[comp] = zVelName;
      comp++;
    }

  vectName[comp] = zsName;
  comp++;
  vectName[comp] = zbottomName;
  comp++;
  vectName[comp] = zbName;
  comp++;

  if (m_write_solver_rhs)
    {
      vectName[comp] = betaName;
      comp++;

      vectName[comp] = C0Name;
      comp++;

      if (SpaceDim == 2)
        {
          vectName[comp] = solverRhsxName;
          comp++;
          vectName[comp] = solverRhsyName;
          comp++;
        }
      else
        {
          MayDay::Error("writeSolverRHS undefined for this dimensionality");
        }
    }

  if (m_write_dHDt)
    {
      vectName[comp] = dthicknessName;      
      comp++;
    } 

  if (m_write_fluxVel)
    {
      vectName[comp] = xfVelName;
      comp++;
      vectName[comp] = yfVelName;
      comp++;
      
      if (writeZvel) 
        {
          vectName[comp] = zfVelName;
          comp++;
        }
    }

 

  if (m_write_baseVel)
    {
      vectName[comp] = xbVelName;
      comp++;
      vectName[comp] = ybVelName;
      comp++;
      
      if (writeZvel) 
        {
          vectName[comp] = zbVelName;
          comp++;
        }
    }

  if (m_write_temperature)
    {
#if BISICLES_Z == BISICLES_LAYERED
  vectName[comp] = temperatureName + string("Surface");
  comp++;
#endif    
      for (int l = 0; l < m_temperature[0]->nComp(); ++l)
	{
	  char idx[4]; sprintf(idx, "%04d", l);
	  vectName[comp] = temperatureName + string(idx);
	  comp++;
	}
    
#if BISICLES_Z == BISICLES_LAYERED
  vectName[comp] = temperatureName + string("Base");
  comp++;
#endif 
    }

#if BISICLES_Z == BISICLES_LAYERED
  if (m_write_layer_velocities){
    for (int l = 0; l < m_nLayers + 1; ++l)
      {
	char idx[4]; sprintf(idx, "%04d", l);
	vectName[comp] = xlayerVelName + string(idx);
	comp++;
	vectName[comp] = ylayerVelName + string(idx);
	comp++;
	if (writeZvel) 
	  {
	    vectName[comp] = zlayerVelName + string(idx);
	    comp++;
	  }
      }
  }
#endif

 
  if (m_write_viscousTensor)
    {
      vectName[comp] = dragCoefName; comp++;
      vectName[comp] = viscosityCoefName; comp++;
      
      vectName[comp] = xxVTname;comp++;
      if (SpaceDim > 1)
	{
	  vectName[comp] = yxVTname;comp++;
	  if (SpaceDim > 2)
	    {
	      vectName[comp] = zzVTname;comp++;
	    }
	  vectName[comp] = xyVTname;comp++;
	  vectName[comp] = yyVTname;comp++;
	  if (SpaceDim > 2)
	    {
	      vectName[comp] = zyVTname;comp++;
	      vectName[comp] = xzVTname;comp++;
	      vectName[comp] = yzVTname;comp++;
	      vectName[comp] = zzVTname;comp++;
	    }
	}
    }


  if (m_write_thickness_sources)
    {
      vectName[comp] = basalThicknessSourceName; comp++;
      vectName[comp] = surfaceThicknessSourceName; comp++;
      vectName[comp] = surfaceThicknessBalanceName; comp++;	
    }

  Box domain = m_amrDomains[0].domainBox();
  Real dt = 1.;
  int numLevels = m_finest_level +1;
  // compute plot data
  Vector<LevelData<FArrayBox>* > plotData(m_velocity.size(), NULL);

  // temp storage for C0
  Vector<LevelData<FArrayBox>* > vectC0(m_velocity.size(), NULL);

  // ghost vect makes things simpler
  IntVect ghostVect(IntVect::Unit);
  
  for (int lev=0; lev<numLevels; lev++)
    {
      // first allocate storage
      plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                               numPlotComps,
                                               ghostVect);

      vectC0[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                             1,
                                             m_velBasalC[0]->ghostVect());
      DataIterator dit = m_amrGrids[lev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*vectC0[lev])[dit].setVal(0.0);
        }
    }

  if (m_write_solver_rhs)
    {
      defineVelRHS(m_velRHS, m_velBasalC, vectC0);
    }




  // Vector<LevelData<FluxBox>*> viscousTensor(numLevels,NULL);
  // Vector<LevelData<FluxBox>*> faceA(numLevels,NULL);
  // Vector<RefCountedPtr<LevelData<FluxBox> > > viscosityCoef;
  // Vector<RefCountedPtr<LevelData<FArrayBox> > > dragCoef;
  // if (m_write_viscousTensor)
  //   {
  //     //need to compute the  Viscous Tensor, which might be expensive
  //     Vector<RealVect> vdx(numLevels);
  //     for (int lev =0; lev < numLevels; lev++)
  // 	{
  // 	  faceA[lev] = new LevelData<FluxBox>(m_amrGrids[lev],m_A[lev]->nComp(),IntVect::Unit);
  // 	  CellToEdge(*m_A[lev],*faceA[lev]);
  // 	  viscousTensor[lev] = new LevelData<FluxBox>(m_amrGrids[lev],SpaceDim,IntVect::Unit);
  // 	  vdx[lev] = RealVect::Unit*m_amrDx[lev];
  // 	}
      
  //     //these parameters don't matter because we don't solve anything here. 
  //     Real vtopSafety = 1.0;
  //     int vtopRelaxMinIter = 4;
  //     Real vtopRelaxTol = 1.0;
  //     Real muMin = 0.0; 
  //     Real muMax = 1.23456789e+300;

  //     IceJFNKstate state(m_amrGrids, m_refinement_ratios, m_amrDomains, vdx, m_vect_coordSys, 
  // 			 m_velocity, m_velBasalC, vectC0, numLevels-1, 
  // 			 *m_constitutiveRelation,  *m_basalFrictionRelation, *m_thicknessIBCPtr,  
  // 			 m_A, faceA, m_time, vtopSafety, vtopRelaxMinIter, vtopRelaxTol, 
  // 			 muMin, muMax);
  //     state.setState(m_velocity);
  //     viscosityCoef = state.getViscosityCoef();
  //     dragCoef = state.getDragCoef();
  //     state.computeViscousTensorFace(viscousTensor);
  //   }


  for (int lev=0; lev<numLevels; lev++)
    {
      // now copy new-time solution into plotData
      Interval thicknessComps(0,0);
      Interval velocityComps(1,2);

      LevelData<FArrayBox>& plotDataLev = *plotData[lev];

      const LevelSigmaCS& levelCS = (*m_vect_coordSys[lev]);
      const LevelData<FArrayBox>& levelH = levelCS.getH();
      const LevelData<FArrayBox>& levelZbase = levelCS.getBaseHeight();
      LevelData<FArrayBox> levelZsurf(m_amrGrids[lev], 1, ghostVect);
      levelCS.getSurfaceHeight(levelZsurf);

      if (m_write_thickness_sources)
	{
	  m_surfaceFluxPtr->surfaceThicknessFlux
	    (*m_surfaceThicknessSource[lev], *this, lev, m_dt);
	  m_basalFluxPtr->surfaceThicknessFlux
	    (*m_basalThicknessSource[lev], *this, lev, m_dt);
	  m_calvingModelPtr->modifySurfaceThicknessFlux
	    (*m_basalThicknessSource[lev],*this, lev, m_dt );
	}

      DataIterator dit = m_amrGrids[lev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = m_amrGrids[lev][dit];
          FArrayBox& thisPlotData = plotDataLev[dit];
          comp = 0;
          const FArrayBox& thisH = levelH[dit];
          
          thisPlotData.copy(thisH, 0, comp, 1);

          comp++;
          const FArrayBox& thisVel = (*m_velocity[lev])[dit];
          thisPlotData.copy(thisVel, 0, comp, 2);
          
          comp += 2;
	 
          if (writeZvel) 
            {
              // use zVel = zero for the moment
              Real zVel = 0.0;
              thisPlotData.setVal(zVel, comp);
              ++comp;
            }

          const FArrayBox& zBase = levelZbase[dit];
          
          // account for background slope of base 
          FArrayBox backgroundBase(thisPlotData.box(), 1);
          BoxIterator bit(thisPlotData.box());
	  const RealVect& basalSlope = levelCS.getBackgroundSlope();
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc(iv);
              loc += 0.5*RealVect::Unit;
              loc *= m_amrDx[lev];

              backgroundBase(iv,0) = D_TERM(loc[0]*basalSlope[0],
                                            +loc[1]*basalSlope[1],
                                            +loc[2]*basalSlope[2]);
            }
          // zsurface
          FArrayBox& zSurf = levelZsurf[dit];
          thisPlotData.copy(zSurf, 0, comp, 1);
          thisPlotData.plus(backgroundBase, 0, comp, 1);
          ++comp;

          // zbottom (bottom of ice
          thisPlotData.copy(zSurf, 0, comp, 1);
          thisPlotData.minus(thisH, 0, comp, 1);
          thisPlotData.plus(backgroundBase, 0, comp, 1);
          ++comp;

          // zbase 
          thisPlotData.copy(zBase, 0, comp, 1);
          thisPlotData.plus(backgroundBase, 0, comp, 1);
          ++comp;

          if (m_write_solver_rhs)
            {
              thisPlotData.copy((*m_velBasalC[lev])[dit],0,comp,1);
              comp++;
              thisPlotData.copy((*vectC0[lev])[dit],0,comp,1);
              comp++;
              thisPlotData.copy((*m_velRHS[lev])[dit],0,comp,SpaceDim);
              comp += SpaceDim;
            }

          // now copy for dthickness/dt 
          if (m_write_dHDt)
            {
              const FArrayBox& thisOldH = (*m_old_thickness[lev])[dit];
              thisPlotData.copy(thisH, 0, comp, 1);
              thisPlotData.minus(thisOldH, 0, comp, 1);
              if (m_dt > 0)
                {
                  thisPlotData.divide(m_dt, comp, 1);

                }              
              ++comp;

	    } // end if we are computing dHDt
      
	  // const FArrayBox& thisSurfaceVel = (*m_velocity[lev])[dit];
          // thisPlotData.copy(thisSurfaceVel, 0, comp, 2);
          
          // comp += 2;

	  

          if (m_write_fluxVel)
            {
              for (int dir = 0; dir < SpaceDim; ++dir)
                {
                  
                  const FArrayBox& thisVel = (*m_faceVel[lev])[dit][dir];
                  for (BoxIterator bit(gridBox); bit.ok(); ++bit)
                    {
                      const IntVect& iv = bit();
                      const IntVect ivp = iv + BASISV(dir);
                      thisPlotData(iv,comp) = half*(thisVel(iv) + thisVel(ivp));
                    }
                  comp++;
                }            

              if (writeZvel) 
                {
                  // use zVel = zero for the moment
                  Real zVel = 0.0;
                  thisPlotData.setVal(zVel, comp);
                  ++comp;
                }
            }



          if (m_write_baseVel)
            {
              const FArrayBox& thisBaseVel = (*m_velocity[lev])[dit];
              thisPlotData.copy(thisBaseVel, 0, comp, 2);
              
              comp += 2;
              
              if (writeZvel) 
                {
                  // use zVel = zero for the moment
                  Real zVel = 0.0;
                  thisPlotData.setVal(zVel, comp);
                  ++comp;
                }	  
            }
	  if (m_write_temperature)
	    {
#if BISICLES_Z == BISICLES_LAYERED
	      {
		const FArrayBox& thisTemp = (*m_sTemperature[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, thisTemp.nComp());
		comp++;
	      }
#endif
	      {
		const FArrayBox& thisTemp = (*m_temperature[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, thisTemp.nComp());
		comp += thisTemp.nComp();
	      }
#if BISICLES_Z == BISICLES_LAYERED
	      {
		const FArrayBox& thisTemp = (*m_bTemperature[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, thisTemp.nComp());
		comp++;
	      }
#endif
	    }
#if BISICLES_Z == BISICLES_LAYERED
	  if (m_write_layer_velocities)
	    {
	      const FArrayBox& thisVel = (*m_layerSFaceXYVel[lev])[dit];
	     
	      for (int j = 0; j < m_nLayers + 1; ++j)
		{
		  thisPlotData.copy(thisVel, j*SpaceDim, comp, SpaceDim);
		  
		  comp+= SpaceDim;
		  // end loop over components
		  if (writeZvel) 
		{
		  // use zVel = zero for the moment
		  Real zVel = 0.0;
		  thisPlotData.setVal(zVel, comp);
		  ++comp;
		} 
	      
		} // end loop over layers
	    }
#endif
	  if (m_write_viscousTensor)
	    {
	      thisPlotData.copy( (*dragCoefficient(lev))[dit],0,comp);
	      comp++;
	      thisPlotData.copy( (*viscosityCoefficient(lev))[dit],0,comp);
	      comp++;
	      thisPlotData.copy( (*viscousTensor(lev))[dit],0,comp, SpaceDim*SpaceDim);
	      comp += SpaceDim * SpaceDim;
	    }
	  

	  // if (m_write_viscousTensor)
	  //   {
	  //     thisPlotData.copy((*dragCoef[lev])[dit],0,comp);
	  //     comp++;

	  //     for (int dir = 0; dir < SpaceDim; ++dir)
	  // 	{
	  // 	  const FArrayBox& thisVisc = (*viscosityCoef[lev])[dit][dir];
		  
	  // 	  for (BoxIterator bit(gridBox); bit.ok(); ++bit)
	  // 	    {
	  // 	      const IntVect& iv = bit();
	  // 	      const IntVect ivp = iv + BASISV(dir);
	  // 	      thisPlotData(iv,comp) = half*(thisVisc(iv) + thisVisc(ivp));
	  // 	    }
	  // 	   comp++;
	  // 	}
	  //     for (int dir = 0; dir < SpaceDim; ++dir)
	  // 	{
	  // 	  const FArrayBox& thisVT = (*viscousTensor[lev])[dit][dir];
	  // 	  for (int vtcomp=0;vtcomp<SpaceDim;vtcomp++)
	  // 	    {
	  // 	      for (BoxIterator bit(gridBox); bit.ok(); ++bit)
	  // 		{
	  // 		  const IntVect& iv = bit();
	  // 		  const IntVect ivp = iv + BASISV(dir);
	  // 		  thisPlotData(iv,comp) = half*(thisVT(iv,vtcomp) + thisVT(ivp,vtcomp));
	  // 		}
	  // 	      comp++;
	  // 	    }
	  // 	}
	  //   }

	  if (m_write_thickness_sources)
	    {
	      thisPlotData.copy((*m_basalThicknessSource[lev])[dit], 0, comp, 1);
	      comp++;
	      thisPlotData.copy((*m_surfaceThicknessSource[lev])[dit], 0, comp, 1);
	      comp++;
	      thisPlotData.copy((*m_balance[lev])[dit], 0, comp, 1);
	      comp++;
	    }

	} // end loop over boxes on this level

      // this is just so that visit surface plots look right
      // fill coarse-fine ghost-cell values with interpolated data
      if (lev > 0)
        {
          PiecewiseLinearFillPatch interpolator(m_amrGrids[lev],
                                                m_amrGrids[lev-1],
                                                numPlotComps,
                                                m_amrDomains[lev-1],
                                                m_refinement_ratios[lev-1],
                                                ghostVect[0]);
          
          // no interpolation in time
          Real time_interp_coeff = 0.0;
          interpolator.fillInterp(*plotData[lev],
                                  *plotData[lev-1],
                                  *plotData[lev-1],
                                  time_interp_coeff,
                                  0, 0,  numPlotComps);
        }
      // just in case...
      //plotData[lev]->exchange();
    } // end loop over levels for computing plot data
  
  // generate plotfile name
  char iter_str[100];
  sprintf(iter_str, "%s%06d.", m_plot_prefix.c_str(), 
          m_cur_step );
  
  string filename(iter_str);
  
  // need to pull out SigmaCS pointers:
  Vector<const LevelSigmaCS* > vectCS(m_vect_coordSys.size(), NULL);
  for (int lev=0; lev<numLevels; lev++)
    {
      vectCS[lev] = dynamic_cast<const LevelSigmaCS* >(&(*m_vect_coordSys[lev]));
    }
  if (m_write_map_file)
    {
      WriteSigmaMappedAMRHierarchyHDF5(filename, m_amrGrids, plotData, vectName, 
                                       vectCS, domain, dt, m_time,
                                       m_refinement_ratios,
                                       numLevels);
    }
  else
    {
      if (SpaceDim == 2)
        {
          filename.append("2d.hdf5");
        } 
      else if (SpaceDim == 3)
        {
          filename.append("3d.hdf5");
        }

      WriteAMRHierarchyHDF5(filename, m_amrGrids, plotData, vectName, 
                            domain, m_amrDx[0], dt, time(), m_refinement_ratios, 
                            numLevels);
    }

  // need to delete plotData
  for (int lev=0; lev<numLevels; lev++)
    {
      if (plotData[lev] != NULL)
        {
          delete plotData[lev];
          plotData[lev] = NULL;
        }      

      if (vectC0[lev] != NULL)
        {
          delete vectC0[lev];
          vectC0[lev] = NULL;
        } 

      // if (faceA[lev] != NULL)
      // 	{
      // 	  delete faceA[lev];
      // 	  faceA[lev] = NULL;
      // 	}
    
      // if (viscousTensor[lev] != NULL)
      // 	{
      // 	  delete viscousTensor[lev];
      // 	  viscousTensor[lev] = NULL;
      // 	}
    }
}

  /// write checkpoint file out for later restarting
void 
AmrIce::writeCheckpointFile() const
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::writeCheckpointfile" << endl;
    }

#ifdef CH_USE_HDF5

  string thicknessName("thickness");
  Vector<string> vectName(1);
  for (int comp=0; comp<1; comp++)
    {
      char idx[4]; sprintf(idx, "%d", comp);
      vectName[comp] = thicknessName+string(idx);
    } 
  Box domain = m_amrDomains[0].domainBox();
  //int numLevels = m_finest_level +1;      

  // generate checkpointfile name
  char (iter_str[100]);
#if 0
  if (m_cur_step < 10)
    {
      sprintf(iter_str, "%s000%d.%dd.hdf5", m_check_prefix.c_str(), 
              m_cur_step, SpaceDim);
    } 
  else if (m_cur_step < 100)
    {
      sprintf(iter_str, "%s00%d.%dd.hdf5", m_check_prefix.c_str(), 
              m_cur_step, SpaceDim);
    }       
  else if (m_cur_step < 1000)
    {
      sprintf(iter_str, "%s0%d.%dd.hdf5", m_check_prefix.c_str(), 
              m_cur_step, SpaceDim);
    }       
  else
    {
      sprintf(iter_str, "%s%d.%dd.hdf5", m_check_prefix.c_str(), 
              m_cur_step, SpaceDim);
    }         
  #endif
  if (m_check_overwrite)
    {
      // overwrite the same checkpoint file, rather than re-writing them
      sprintf(iter_str, "%s.%dd.hdf5", m_check_prefix.c_str(), SpaceDim);
    }
  else 
    {
      // or hang on to them, if you are a bit sentimental. It's better than keeping
      // every core dump you generate.
      sprintf(iter_str, "%s%06d.%dd.hdf5", m_check_prefix.c_str(), m_cur_step, SpaceDim);
    }

  if (s_verbosity > 3) 
    {
      pout() << "checkpoint file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str, HDF5Handle::CREATE);

  // write amr data -- only dump out things which are essential
  // to restarting the computation (i.e. max_level, finest_level, 
  // time, refinement ratios, etc.).  Other paramters (regrid 
  // intervals, block-factor, etc can be changed by the inputs
  // file of the new run.
  // At the moment, the maximum level is not allowed to change,
  // although in principle, there is no real reason why it couldn't
  // 
  HDF5HeaderData header;
  header.m_int["max_level"] = m_max_level;
  header.m_int["finest_level"] = m_finest_level;
  header.m_int["current_step"] = m_cur_step;
  header.m_real["time"] = m_time;
  header.m_real["dt"] = m_dt;
  header.m_int["num_comps"] = 2 +  m_velocity[0]->nComp() 
    + m_temperature[0]->nComp();
#if BISICLES_Z == BISICLES_LAYERED
  header.m_int["num_comps"] +=2; // surface and base temperatures
#endif
  // at the moment, save cfl, but it can be changed by the inputs
  // file if desired.
  header.m_real["cfl"] = m_cfl;

  // periodicity info
  D_TERM(
         if (m_amrDomains[0].isPeriodic(0))
         header.m_int["is_periodic_0"] = 1;
         else
         header.m_int["is_periodic_0"] = 0; ,

         if (m_amrDomains[0].isPeriodic(1))
         header.m_int["is_periodic_1"] = 1;
         else
         header.m_int["is_periodic_1"] = 0; ,

         if (m_amrDomains[0].isPeriodic(2))
         header.m_int["is_periodic_2"] = 1;
         else
         header.m_int["is_periodic_2"] = 0; 
         );
         

  // set up component names
  char compStr[30];
  //string thicknessName("thickness");
  string compName;
  int nComp = 0;
  for (int comp=0; comp < 1; comp++)
    {
      // first generate component name
      char idx[5]; sprintf(idx, "%04d", comp);
      compName = thicknessName + string(idx);
      sprintf(compStr, "component_%04d", comp);
      header.m_string[compStr] = compName;
     
    }
  nComp++;

  string baseHeightName("bedHeight");
  for (int comp=0; comp < 1; comp++)
    {
      // first generate component name
      char idx[5]; sprintf(idx, "%04d", comp);
      compName = baseHeightName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
      
    }
  nComp++;

  string velocityName("velocity");
  for (int comp=0; comp < m_velocity[0]->nComp() ; comp++) 
    {
      // first generate component name
      char idx[5]; sprintf(idx, "%04d", comp);
      compName = velocityName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  nComp += m_velocity[0]->nComp() ;

  string temperatureName("temperature");
  for (int comp=0; comp < m_temperature[0]->nComp() ; comp++) 
    {
      char idx[5]; sprintf(idx, "%04d", comp);
      compName = temperatureName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  
  nComp += m_temperature[0]->nComp() ;

#if BISICLES_Z == BISICLES_LAYERED
  {
    sprintf(compStr, "component_%04d", nComp);
    compName = "sTemperature";
    header.m_string[compStr] = compName;
    nComp += 1;
    sprintf(compStr, "component_%04d", nComp);
    compName = "bTemperature";
    header.m_string[compStr] = compName;
    nComp += 1;
    //layer data
    const Vector<Real>& sigma = getFaceSigma();
    string s("sigma");
    for (int l =0; l < sigma.size(); ++l)
      {
	char idx[5]; sprintf(idx, "%04d", l);
	header.m_real[s + string(idx)] = sigma[l];  
      }
  }
#endif

  header.writeToFile(handle);

  // now loop over levels and write out each level's data
  // note that we loop over all allowed levels, even if they
  // are not defined at the moment.
  for (int lev=0; lev<= m_max_level; lev++)
    {
      // set up the level string
      char levelStr[20];
      sprintf(levelStr, "%d", lev);
      const std::string label = std::string("level_") + levelStr;
      
      handle.setGroup(label);
      
      // set up the header info
      HDF5HeaderData levelHeader;
      if (lev < m_max_level)
        {
          levelHeader.m_int["ref_ratio"] = m_refinement_ratios[lev];
        }
      levelHeader.m_real["dx"] = m_amrDx[lev];
      levelHeader.m_box["prob_domain"] = m_amrDomains[lev].domainBox();
      
      levelHeader.writeToFile(handle);
      
      // now write the data for this level
      // only try to write data if level is defined.
      if (lev <= m_finest_level)
        {
          write(handle, m_amrGrids[lev]);
	 
	 

	  // LevelData<FArrayBox> tmp;
	  const IntVect ghost = IntVect::Unit*2;
	  // tmp.define(m_amrGrids[lev],1,ghost);
	  
          const LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
          // const LevelData<FArrayBox>& levelH = levelCS.getH();
	  // for (DataIterator dit = tmp.dataIterator();
	  //      dit.ok(); ++dit){
	  //   tmp[dit].copy( levelH[dit]);
	  // }

	  write(handle, levelCS.getH() , "thicknessData", levelCS.getH().ghostVect());

          // const LevelData<FArrayBox>& levelZb = levelCS.getBaseHeight();
	  // for (DataIterator dit = tmp.dataIterator();
	  //      dit.ok(); ++dit){
	  //   tmp[dit].copy( levelZb[dit]);
	  // }
	  write(handle, levelCS.getBaseHeight() , "bedHeightData",
		levelCS.getBaseHeight().ghostVect()  );

	  write(handle, *m_velocity[lev], "velocityData", 
		m_velocity[lev]->ghostVect());

	  write(handle, *m_temperature[lev], "temperatureData", 
		m_temperature[lev]->ghostVect());

#if BISICLES_Z == BISICLES_LAYERED
	  write(handle, *m_sTemperature[lev], "sTemperatureData", 
		m_sTemperature[lev]->ghostVect());
	  write(handle, *m_bTemperature[lev], "bTemperatureData", 
		m_bTemperature[lev]->ghostVect());
#endif

        }
    }// end loop over levels
  
  handle.close();
#endif
}


/// read checkpoint file for restart 
void 
AmrIce::readCheckpointFile(HDF5Handle& a_handle)
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::readCheckpointFile" << endl;
    }

#ifndef CH_USE_HDF5
  MayDay::Error("code must be compiled with HDF5 to read checkpoint files");
#endif

#ifdef CH_USE_HDF5
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data: " << endl;
      pout() << header << endl;
    }

  // read max level
  if (header.m_int.find("max_level") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain max_level");
    }
  // we can change max level upon restart
  int max_level_check = header.m_int["max_level"];
  if (max_level_check != m_max_level)
    {
      if (s_verbosity > 0)
        {
          pout() << "Restart file has a different max level than inputs file"
                 << endl;
          pout() << "     max level from inputs file = " 
                 << m_max_level << endl;
          pout() << "     max level in checkpoint file = " 
                 << max_level_check << endl;                 
          pout() << "Using max level from inputs file" << endl;
        }
    }
  // read finest level
  if (header.m_int.find("finest_level") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain finest_level");
    }

  m_finest_level = header.m_int["finest_level"];
  if (m_finest_level > m_max_level)
    {
      MayDay::Error("finest level in restart file > max allowable level!");
    }

  // read current step
  if (header.m_int.find("current_step") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain current_step");
    }

  m_cur_step = header.m_int["current_step"];
  m_restart_step = m_cur_step;

  // read time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain time");
    }

  m_time = header.m_real["time"];

  // read timestep
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain dt");
    }

  m_dt = header.m_real["dt"];

  // read num comps
  if (header.m_int.find("num_comps") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain num_comps");
    }

  // read cfl
  if (header.m_real.find("cfl") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain cfl");
    }

  Real check_cfl = header.m_real["cfl"];
  ParmParse ppCheck("amr");

    if (ppCheck.contains("cfl"))
      { 
        // check for consistency and warn if different
        if (check_cfl != m_cfl)
          {
            if (s_verbosity > 0)
              {
                pout() << "CFL in checkpoint file different from inputs file" 
                       << endl;
                pout() << "     cfl in inputs file = " << m_cfl << endl;
                pout() << "     cfl in checkpoint file = " << check_cfl 
                       << endl;
                pout() << "Using cfl from inputs file" << endl;                
              }
          }  // end if cfl numbers differ
      } // end if cfl present in inputs file
    else
      {
        m_cfl = check_cfl;
      }          

  // read periodicity info
  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility 
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
           isPeriodic[0] = false; ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
           isPeriodic[1] = false; ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
           isPeriodic[2] = false;);

#if BISICLES_Z == BISICLES_LAYERED
  //retrieve sigma data
  Vector<Real> sigma;
  int l = 0;
  string s("sigma");
  bool found = false;
  do {
    char idx[6]; sprintf(idx, "%04d", l);
    string ss = s + string(idx);
    map<std::string, Real>::const_iterator it = header.m_real.find(ss);
    found = (it != header.m_real.end());
    if (found)
      sigma.push_back(it->second);
    ++l;
  } while (found);
  m_nLayers = sigma.size() - 1;
  CH_assert(m_nLayers > 0 && sigma[0] < TINY_NORM && abs(sigma[m_nLayers] - 1.0) < TINY_NORM);
#endif

  // now resize stuff 
  m_amrDomains.resize(m_max_level+1);
  m_amrGrids.resize(m_max_level+1);
  m_amrDx.resize(m_max_level+1);
  m_old_thickness.resize(m_max_level+1, NULL);
  m_velocity.resize(m_max_level+1, NULL);
  m_diffusivity.resize(m_max_level+1);
  m_vect_coordSys.resize(m_max_level+1);
  m_velRHS.resize(m_max_level+1);
  m_surfaceThicknessSource.resize(m_max_level+1,NULL);
  m_basalThicknessSource.resize(m_max_level+1,NULL);
  m_balance.resize(m_max_level+1,NULL);
  m_velBasalC.resize(m_max_level+1,NULL);
  m_cellMuCoef.resize(m_max_level+1,NULL);
  m_faceMuCoef.resize(m_max_level+1,NULL);
  m_faceVel.resize(m_max_level+1,NULL);
  m_temperature.resize(m_max_level+1,NULL);
#if BISCICLES_Z == BISICLES_LAYERED
  m_sTemperature.resize(m_max_level+1,NULL);
  m_bTemperature.resize(m_max_level+1,NULL);
  m_layerSFaceXYVel.resize(m_max_level+1,NULL);
  m_layerXYFaceXYVel.resize(m_max_level+1,NULL);
#endif
  IntVect sigmaCSGhost = m_num_thickness_ghost*IntVect::Unit;
	 

  // now read in level-by-level data
  for (int lev=0; lev<= m_max_level; lev++)
    {
      // set up the level string
      char levelStr[20];
      sprintf(levelStr, "%d", lev);
      const std::string label = std::string("level_") + levelStr;
      
      a_handle.setGroup(label);

      // read header info
      HDF5HeaderData header;
      header.readFromFile(a_handle);
      
      if (s_verbosity >= 3)
        {
          pout() << "level " << lev << " header data" << endl;
          pout() << header << endl;
        }

  // Get the refinement ratio
      if (lev < m_max_level)
        {
          int checkRefRatio;
          if (header.m_int.find("ref_ratio") == header.m_int.end())
            {
              MayDay::Error("checkpoint file does not contain ref_ratio");
            }
          checkRefRatio = header.m_int["ref_ratio"];

          // check for consistency
          if (checkRefRatio != m_refinement_ratios[lev])
            {
	      
	      MayDay::Error("inputs file and checkpoint file ref ratios inconsistent");
            }
        }
      
      // read dx
      if (header.m_real.find("dx") == header.m_real.end())
        {
          MayDay::Error("checkpoint file does not contain dx");
        }
      
      m_amrDx[lev] = header.m_real["dx"];
      
      // read problem domain box
      if (header.m_box.find("prob_domain") == header.m_box.end())
        {
          MayDay::Error("checkpoint file does not contain prob_domain");
        }
      Box domainBox = header.m_box["prob_domain"];

      m_amrDomains[lev] = ProblemDomain(domainBox, isPeriodic);


      // the rest is only applicable if this level is defined
      if (lev <= m_finest_level)
        {
          // read grids          
          Vector<Box> grids;
          const int grid_status = read(a_handle, grids);
          if (grid_status != 0) 
            {
              MayDay::Error("checkpoint file does not contain a Vector<Box>");
            }
          // do load balancing
          int numGrids = grids.size();
          Vector<int> procIDs(numGrids);
          LoadBalance(procIDs, grids);
          DisjointBoxLayout levelDBL(grids, procIDs, m_amrDomains[lev]);
          m_amrGrids[lev] = levelDBL;

          // allocate this level's storage
	  // 4 ghost cells needed for advection.
          m_old_thickness[lev] = new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
#if BISICLES_Z == BISICLES_LAYERED
	  m_temperature[lev] =  new LevelData<FArrayBox>
	    (levelDBL, m_nLayers, m_num_thickness_ghost*IntVect::Unit);
	  m_sTemperature[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
	  m_bTemperature[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
#elif BISICLES_Z == BISICLES_FULLZ
	  m_temperature[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
#endif
	  // other quantities need only one;
	  IntVect ghostVect(IntVect::Unit);
          m_velocity[lev] = new LevelData<FArrayBox>(levelDBL, 2, 
                                                     ghostVect);

	  m_faceVel[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Unit);
#if BISICLES_Z == BISICLES_LAYERED
	  m_layerXYFaceXYVel[lev] = new LevelData<FluxBox>
	    (m_amrGrids[lev], m_nLayers, IntVect::Unit);
	  m_layerSFaceXYVel[lev] = new LevelData<FArrayBox>
	    (m_amrGrids[lev], SpaceDim*(m_nLayers + 1), IntVect::Unit);
#endif

	  m_velBasalC[lev] = new LevelData<FArrayBox>(levelDBL, 1, ghostVect);
	  m_cellMuCoef[lev] = new LevelData<FArrayBox>(levelDBL, 1, ghostVect);
	  m_faceMuCoef[lev] = new LevelData<FluxBox>(levelDBL, 1, ghostVect);
	  m_velRHS[lev] = new LevelData<FArrayBox>(levelDBL, 2, IntVect::Zero);


	  m_surfaceThicknessSource[lev] = 
	    new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_basalThicknessSource[lev] = 
	    new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_balance[lev] = 
	    new LevelData<FArrayBox>(levelDBL,   1, IntVect::Zero) ;
	  m_diffusivity[lev] =  RefCountedPtr<LevelData<FluxBox> >
	    (new LevelData<FluxBox>(levelDBL, 1,  ghostVect));

          // read this level's data
	  
          LevelData<FArrayBox>& old_thickness = *m_old_thickness[lev];  
	  //IntVect ghost = old_thickness.ghostVect();
          int dataStatus = read<FArrayBox>(a_handle,
                                           old_thickness,
                                           "thicknessData",
                                           levelDBL);
	  //CH_assert( old_thickness.ghostVect() == ghost);

          if (dataStatus != 0)
            {
              MayDay::Error("checkpoint file does not contain thickness data");
            }


	  LevelData<FArrayBox> bedHeight;
	  bedHeight.define(old_thickness);
	  dataStatus = read<FArrayBox>(a_handle,
				       bedHeight,
				       "bedHeightData",
				       levelDBL);

	  if (dataStatus != 0)
            {
              MayDay::Error("checkpoint file does not contain bed height data");
            }

	  //having read thickness and base data, we can define
          //the co-ordinate system 

	  RealVect dx = m_amrDx[lev]*RealVect::Unit;
          m_vect_coordSys[lev] = RefCountedPtr<LevelSigmaCS >
            (new LevelSigmaCS(m_amrGrids[lev], dx, sigmaCSGhost));
	  m_vect_coordSys[lev]->setIceDensity(m_iceDensity);
	  m_vect_coordSys[lev]->setWaterDensity(m_seaWaterDensity);
	  m_vect_coordSys[lev]->setGravity(m_gravity);
#if BISICLES_Z == BISICLES_LAYERED
	  m_vect_coordSys[lev]->setFaceSigma(sigma);
#endif
          LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
          LevelData<FArrayBox>& levelH = levelCS.getH();

          DataIterator dit = levelH.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              levelH[dit].copy((*m_old_thickness[lev])[dit]);
            }
          levelCS.setBaseHeight(bedHeight);

	  {
	    LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
	    int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
	    levelCS.recomputeGeometry(crseCoords, refRatio);
	  }
          LevelData<FArrayBox>& velData = *m_velocity[lev];
           dataStatus = read<FArrayBox>(a_handle,
                                       velData,
                                       "velocityData",
				       levelDBL);
	   m_velocitySolveInitialResidualNorm = 1.0e+6; //\todo fix this
	   // dit.reset();
	   // for (dit.begin(); dit.ok(); ++dit)
	   //   {
	   //     (*m_velocity[lev])[dit].setVal(0.0);
	   //   }
	  
	   //this check doesn't work, because HDF5::SetGroup attempts
	   //to create a group if it doesn't exist. And since the file has been
	   // opened readonly, the previous call generates an exception
 			       
          if (dataStatus != 0)
            {
	      MayDay::Error("checkpoint file does not contain velocity data");
	      
            }

	  LevelData<FArrayBox>& temperatureData = *m_temperature[lev];
	  dataStatus = read<FArrayBox>(a_handle,
				       temperatureData,
                                       "temperatureData",
				       levelDBL);

	  if (dataStatus != 0)
	    {
	      MayDay::Error("checkpoint file does not contain temperature data"); 
            }

#if BISICLES_Z == BISICLES_LAYERED
	  
	  CH_assert(m_temperature[lev]->nComp() == sigma.size() -1);
	  
	  LevelData<FArrayBox>& sTemperatureData = *m_sTemperature[lev];
	  dataStatus = read<FArrayBox>(a_handle,
				       sTemperatureData,
                                       "sTemperatureData",
				       levelDBL);

	  if (dataStatus != 0)
	    {
	      MayDay::Error("checkpoint file does not contain sTemperature data"); 
            }

	  
	  LevelData<FArrayBox>& bTemperatureData = *m_bTemperature[lev];
	  dataStatus = read<FArrayBox>(a_handle,
				       bTemperatureData,
                                       "bTemperatureData",
				       levelDBL);

	  if (dataStatus != 0)
	    {
	      MayDay::Error("checkpoint file does not contain bTemperature data"); 
	    }


  
#endif

        } // end if this level is defined
    } // end loop over levels                                    
          
  
  // do we need to close the handle?
 
  //this is just to make sure the diffusivity is computed
  //(so I should improve that)
  defineSolver();
  m_doInitialVelSolve = false; // since we have just read the velocity field
  m_doInitialVelGuess = false; // ditto
  solveVelocityField(m_velocity);
  m_doInitialVelSolve = true;

#endif
  
}

/// set up for restart
void 
AmrIce::restart(string& a_restart_file)
{
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::restart" << endl;
    }

  HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
  // first read in data from checkpoint file
  readCheckpointFile(handle);
  handle.close();
  // don't think I need to do anything else, do I?


}
#endif

Real AmrIce::computeTotalGroundedIce() const
{
  
  Real totalGroundedIce = 0;

  Vector<LevelData<FArrayBox>* > vectGroundedThickness(m_finest_level+1, NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const LevelData<FArrayBox>& levelThickness = m_vect_coordSys[lev]->getH();
      // temporary with only ungrounded ice
      vectGroundedThickness[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1,
                                                            IntVect::Zero);

      LevelData<FArrayBox>& levelGroundedThickness = *vectGroundedThickness[lev];
      // now copy thickness to       
      levelThickness.copyTo(levelGroundedThickness);

      const LevelData<BaseFab<int> >& levelMask = m_vect_coordSys[lev]->getFloatingMask();
      // now loop through and set to zero where we don't have grounded ice.
      // do this the slow way, for now
      DataIterator dit=levelGroundedThickness.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const BaseFab<int>& thisMask = levelMask[dit];
          FArrayBox& thisThick = levelGroundedThickness[dit];
          BoxIterator bit(thisThick.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              if (thisMask(iv,0) != GROUNDEDMASKVAL)
                {
                  thisThick(iv,0) = 0.0;
                }
            }
        }
    

    }

  // now compute sum
  Interval thicknessInt(0,0);
  totalGroundedIce = computeSum(vectGroundedThickness, m_refinement_ratios,
                                m_amrDx[0], thicknessInt, 0);

  
  // clean up temp storage
  for (int lev=0; lev<vectGroundedThickness.size(); lev++)
    {
      if (vectGroundedThickness[lev] != NULL)
        {
          delete vectGroundedThickness[lev];
          vectGroundedThickness[lev] = NULL;
        }
    }

  return totalGroundedIce;

}

void AmrIce::helmholtzSolve
(Vector<LevelData<FArrayBox>* >& a_phi,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 Real a_alpha, Real a_beta) const
{

  // AMRPoissonOp supports only one component of phi
  // if m_finest_level > 0 (its LevelFluxRegisters are
  // defined with only one component, so we will do
  // one component at a time. should try to avoid some
  // of the rhs copies...
  for (int icomp = 0; icomp < a_phi[0]->nComp(); ++icomp)
    {
      
      //make a copy of a_phi with one ghost cell
      Vector<LevelData<FArrayBox>* > phi(m_finest_level + 1, NULL);
      Vector<LevelData<FArrayBox>* > rhs(m_finest_level + 1, NULL);
      Vector<DisjointBoxLayout> grids(m_finest_level + 1);
      for (int lev=0; lev < m_finest_level + 1; ++lev)
	{
	  grids[lev] = m_amrGrids[lev];

	  const LevelData<FArrayBox>& levelPhi = *a_phi[lev]; 
	  phi[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 
					      1, IntVect::Unit);
	  levelPhi.copyTo(Interval(icomp,icomp),*phi[lev], Interval(0,0));
	  phi[lev]->exchange();
	  

	  const LevelData<FArrayBox>& levelRhs = *a_rhs[lev];
	  rhs[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 
					      1, IntVect::Zero);
	  levelRhs.copyTo(Interval(icomp,icomp),*rhs[lev], Interval(0,0));
	  rhs[lev]->exchange();
      
    }


      //Natural boundary conditions
      BCHolder bc(ConstDiriNeumBC(IntVect(0,0), RealVect(0.0,0.0),
				  IntVect(0,0), RealVect(0.0,0.0)));
      
      
      AMRPoissonOpFactory opf;
      opf.define(m_amrDomains[0],  grids , m_refinement_ratios,
		 m_amrDx[0], bc, a_alpha, -a_beta );
      
      AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
      BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
      mgSolver.define(m_amrDomains[0], opf, &bottomSolver, m_finest_level+1);
      mgSolver.m_eps = TINY_NORM;
      mgSolver.m_normThresh = TINY_NORM;
      mgSolver.m_iterMax = 8;
      int numMGSmooth = 4;
      mgSolver.m_pre = numMGSmooth;
      mgSolver.m_post = numMGSmooth;
      mgSolver.m_bottom = numMGSmooth;
      mgSolver.m_verbosity = s_verbosity - 1;
      
      mgSolver.solve(phi, rhs, m_finest_level, 0,  false);
      
      for (int lev=0; lev < m_finest_level + 1; ++lev)
	{
	  LevelData<FArrayBox>& levelPhi = *a_phi[lev];
	  phi[lev]->copyTo(Interval(0,0), levelPhi, Interval(icomp, icomp));
	  
	  if (phi[lev] != NULL){
	    delete phi[lev];
	    phi[lev] = NULL;
	  }
	}
    }
}


void AmrIce::helmholtzSolve
(Vector<LevelData<FArrayBox>* >& a_phi, Real a_alpha, Real a_beta) const
{
  
  Vector<LevelData<FArrayBox>* > rhs(m_finest_level + 1, NULL);
 
  for (int lev=0; lev < m_finest_level + 1; ++lev)
    {
      const LevelData<FArrayBox>& levelPhi = *a_phi[lev]; 
      rhs[lev] = new LevelData<FArrayBox>
	(m_amrGrids[lev], levelPhi.nComp(), IntVect::Zero);
      levelPhi.copyTo(*rhs[lev]);
    }
 
  helmholtzSolve(a_phi, rhs, a_alpha, a_beta);

   for (int lev=0; lev < m_finest_level + 1; ++lev)
     {
       if (rhs[lev] != NULL){
	 delete rhs[lev];
	 rhs[lev] = NULL;
       }
     }

}



// void AmrIce::updateSurfaceGradient(Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys,
// 				   const Vector<DisjointBoxLayout>& a_grids, 
// 				   const Vector<ProblemDomain>& a_domains,
// 				   const Vector<int>& a_refinementRatios,
// 				   const Vector<Real>& a_dx, const RealVect a_basalSlope,
// 				   Real a_time, Real a_dt,
// 				   int a_clipRadius,
// 				   const int a_lmin, const int a_lmax,
// 				   IceThicknessIBC* a_thicknessIBCPtr)
// {
//   CH_assert(a_lmin >= 0 && a_lmin <  a_coordSys.size());
//   CH_assert(a_lmax >= 0 && a_lmax <  a_coordSys.size());
//   CH_assert(2 == SpaceDim);
//   MayDay::Error("AmrIce::updateSurfaceGradient is now doomed");

//   Vector<LevelData<FArrayBox>* > tempZs(a_lmax + 1, NULL);
  
//   //first, construct surface elevation : 
//   // \todo this probably ought to be part of levelSigmaCS, or at least
//   // live in its own function but for now
//   for (int lev= a_lmin ; lev <= a_lmax ; lev++)
//     {
//       const DisjointBoxLayout& levelGrids = a_grids[lev];
//       //two ghost cells are needed by the L1L2 model, which needs
//       //grad(s) in the first ghost cell.
//       tempZs[lev] = new LevelData<FArrayBox>(levelGrids, 1,  2 * IntVect::Unit);
//       LevelData<FArrayBox>& levelZs = (*tempZs[lev]);
//       LevelSigmaCS& levelCS = *a_coordSys[lev];
//       levelCS.getSurfaceHeight(levelZs);
      
//       if (lev > 0)
// 	{
// 	  //fill ghost regions of levelZs
// 	  int nGhost = levelZs.ghostVect()[0];
// 	  PiecewiseLinearFillPatch surfaceFiller
// 	    (levelGrids,a_grids[lev-1],1, a_domains[lev-1],
// 	     a_refinementRatios[lev-1],nGhost);
	  
// 	  Real time_interp_coeff = 0.0;
// 	  surfaceFiller.fillInterp
// 	    (levelZs,  *tempZs[lev-1],*tempZs[lev-1],
// 	     time_interp_coeff, 0, 0, 1);
	  
// 	}
      
//       // set any non-periodic BC's on surface height
//       RealVect levelDx = a_dx[lev]*RealVect::Unit;
//       a_thicknessIBCPtr->setSurfaceHeightBCs
// 	(*tempZs[lev],*(a_coordSys[lev]),a_domains[lev],
// 	 levelDx, a_time, a_dt);

//       //tempZs[lev]->exchange();
      
//       // put some sensible values into the corners
//       if (SpaceDim == 2)
//       	{
//       	  for (DataIterator dit(levelZs.dataIterator()); dit.ok(); ++dit)
//       	    {
//       	      FArrayBox &thisZs = levelZs[dit];
//       	      Box sbox = thisZs.box();
//       	      sbox.grow(-levelZs.ghostVect()[0]);
//       	      FORT_EXTRAPCORNER2D(CHF_FRA(thisZs),
//       				  CHF_BOX(sbox));
//       	    }
	  
//       	}

//       levelCS.setSurfaceHeight(levelZs);
//     } // end construct surface elevation
  
//   // next, compute grad(s) at cell centers and faces, and store 
//   for (int lev=a_lmin; lev<= a_lmax; lev++)
//     {
//       LevelSigmaCS& levelCS = *a_coordSys[lev];
      
//       LevelData<FArrayBox>& levelZs = (*tempZs[lev]);
//       const DisjointBoxLayout& levelGrids = a_grids[lev];
//       //temporaries \todo eliminate these
//       LevelData<FArrayBox> levelGrad(levelGrids, SpaceDim, IntVect::Unit);
//       LevelData<FluxBox> levelFaceGrad(levelGrids, SpaceDim, IntVect::Zero);
//       const RealVect& dx = levelCS.dx();
      
//       // compute maxGrad , to later limit grad(s) to max(s)/(gradRadius*dx)
//       Real maxZs, maxGrad;
//       if (a_clipRadius > 0)
// 	{
// 	  Interval comps(0,0);
// 	  maxZs = computeMax(tempZs, a_refinementRatios[lev], comps, 0);
// 	  maxGrad = maxZs/Real(a_clipRadius);
// 	}
	 
//       for (DataIterator dit(levelGrids); dit.ok(); ++dit)
// 	{
// 	  Real ratio = levelCS.iceDensity() / levelCS.waterDensity();
// 	  const Real& seaLevel = levelCS.seaLevel();
// 	  const BaseFab<int>& thisMask = levelCS.getFloatingMask()[dit];
// 	  const FArrayBox& thisTopg = levelCS.getTopography()[dit];
// 	  const FArrayBox& thisH = levelCS.getH()[dit];
// 	  FArrayBox& thisSurf = levelZs[dit];
// 	  const Box gridBox = levelGrids[dit];
// 	  //we need grad(s) in the first ghost, so
// 	  Box gridBoxPlus = gridBox;
// 	  gridBoxPlus.grow(1);
// 	  FArrayBox& grad = levelGrad[dit];

// 	  // set surface elevation for any open sea or land regions
// 	  FORT_SETOPENSURFACE(CHF_FRA1(thisSurf,0),
// 			      CHF_CONST_FIA1(thisMask,0),
// 			      CHF_CONST_FRA1(thisTopg,0),
// 			      CHF_CONST_REAL(seaLevel),
// 			      CHF_BOX(thisSurf.box()));

// 	  for (int dir =0; dir < SpaceDim; ++dir)
// 	    {
// 	      Box faceBox = gridBoxPlus;
// 	      faceBox.growHi(dir,1);
// 	      // FORT_GLGRADS needs a workspace at cell faces
// 	      FArrayBox faceS(faceBox,1);
	      
// 	      FORT_GLGRADS(CHF_FRA1(grad,dir),
// 			   CHF_FRA1(faceS,0),
// 			   CHF_CONST_FRA1(thisH,0),
// 			   CHF_CONST_FRA1(thisSurf,0),
// 			   CHF_CONST_FRA1(thisTopg,0),
// 			   CHF_CONST_FIA1(thisMask,0),
// 			   CHF_CONST_REAL(ratio),
// 			   CHF_CONST_REAL(seaLevel),
// 			   CHF_CONST_REAL(dx[dir]),
// 			   CHF_CONST_INT(dir),
// 			   CHF_BOX(gridBoxPlus),
// 			   CHF_BOX(faceBox));

// 	      //limit grad if required
// 	      if (a_clipRadius > 0)
// 		{
// 		  Real maxGradDir = maxGrad / dx[dir];
// 		  for (BoxIterator bit(gridBoxPlus); bit.ok(); ++bit)
// 		    {
// 		      IntVect iv = bit();
		      
// 		      if (grad(iv,dir) > maxGradDir)
// 			grad(iv,dir) = maxGradDir;
// 		      else if (grad(iv,dir) < -maxGradDir)
// 			grad(iv,dir) = -maxGradDir;
// 		    }
// 		}
	      
// 	      CH_assert(grad.norm(0) < HUGE_NORM);

// 	    }  

// 	  //grad(s) at cell faces
// 	  for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
// 	    {
// 	      Real maxGradDir = maxGrad / dx[faceDir];
// 	      FArrayBox& thisFaceGrad = levelFaceGrad[dit][faceDir];
// 	      const Box& faceBox = thisFaceGrad.box();
// 	      Real oneOnDx = 1.0 / dx[faceDir];
// 	      for (BoxIterator fbit(faceBox); fbit.ok(); ++fbit)
// 		{
// 		  IntVect iv = fbit(); 
// 		  thisFaceGrad(iv, faceDir) = oneOnDx * (thisSurf(iv,0) 
// 			       - thisSurf(iv - BASISV(faceDir),0));

// 		  //limiting
// 		  if (a_clipRadius > 0)
// 		    {
			 
// 		      if (thisFaceGrad(iv,faceDir) > maxGradDir)
// 			thisFaceGrad(iv,faceDir) = maxGradDir;
// 		      else if (thisFaceGrad(iv,faceDir) < -maxGradDir)
// 			thisFaceGrad(iv,faceDir) = -maxGradDir;
// 		    }
		     

// 		}
// 	      if (SpaceDim == 2)
// 		{
// 		  int transDir = (faceDir+1)%SpaceDim;
// 		  Real oneOnFourDx = 0.25 / dx[faceDir];
// 		  Real maxGradDir = maxGrad / dx[transDir];
// 		  for (BoxIterator fbit(faceBox) ; fbit.ok(); ++fbit)
// 		    {
// 		      IntVect iv = fbit(); 
// 		      thisFaceGrad(iv, transDir) =  oneOnFourDx
// 			* (thisSurf(iv + BASISV(transDir),0) 
// 			   - thisSurf(iv - BASISV(transDir),0)
// 			   + thisSurf(iv + BASISV(transDir) 
// 				       - BASISV(faceDir),0) 
// 			   - thisSurf(iv - BASISV(transDir) 
// 				       - BASISV(faceDir),0));

// 		      //limiting
// 		      if (a_clipRadius > 0)
// 			{ 
			 
// 			  if (thisFaceGrad(iv,transDir) > maxGradDir)
// 			    thisFaceGrad(iv,transDir) = maxGradDir;
// 			  else if (thisFaceGrad(iv,transDir) < -maxGradDir)
// 			    thisFaceGrad(iv,transDir) = -maxGradDir;
// 			}

// 		    }
// 		} 
// 	    }
	     

// 	} // end loop over boxes
	 
//       //store the results
//       levelGrad.exchange();
//       levelCS.setGradSurface(levelGrad); 
//       //      levelFaceGrad.exchange();
//       levelCS.setGradSurfaceFace(levelFaceGrad);

//     } // end loop over levels

//   // clean up temp storage
//   for (int lev=0; lev<tempZs.size(); lev++)
//     {
//       if (tempZs[lev] != NULL)
// 	{
// 	  delete tempZs[lev];
// 	  tempZs[lev] = NULL;
// 	}
//     }

// } // end setSurfaceGradient

/// during AmrIce::regrid, we interpolate thickness a_oldH onto the
/// newly created higher levels. This has an unfortunate side effect : 
/// because gradients of H are steep glose to the grounding line, we tend
/// to turn floating ice into grounded ice there. 
///
///  | G | F | F   : grouded(G) and floating(F) ice prior to interpolation
///  |G|G|G|F|F|F  : grouded(G) and floating(F) ice post interpolation
///
/// This has a pretty severe effect, moving the turning point in the velocity
/// a cell seaward
/// 
/// Here, we attempt correct this without resorting to pure piecewise constant
/// interpolation
// void AmrIce::postInterpolationReFloat(LevelData<FArrayBox>& a_H,
// 				      const LevelData<FArrayBox>& a_coarseH,
// 				      const LevelData<FArrayBox>& a_coarseBed,
// 				      const DisjointBoxLayout a_newDBL,
// 				      const ProblemDomain& a_domain,
// 				      int a_refRatio,
// 				      Real a_seaLevel, 
// 				      Real a_waterDensity, 
// 				      Real a_iceDensity)
// {
  
//   // this is a little unpleasant, but for now. The idea is simple: 
//   // wherever piecwise constant intepolation of gives floating (or grounded)
//   // ice, adjust a_H to ensure the same. Probably ought to do 
//   // a_refRatio * a_refRatio blocks of cells to conserve H.
//   // We ought to get passed both a_coarseBed and bed, rather than
//   // doing piecewise linear interpolation here.
  
  
//   FineInterp interpolator(a_newDBL, 1, a_refRatio, a_domain);
  
//   LevelData<FArrayBox> bed(a_newDBL, 1, a_coarseBed.ghostVect());
//   interpolator.interpToFine(bed, a_coarseBed);
  
//   LevelData<FArrayBox> pwcBed(a_newDBL, 1, a_coarseBed.ghostVect());
//   interpolator.pwcinterpToFine(pwcBed, a_coarseBed);
  
//   LevelData<FArrayBox> pwcH(a_newDBL, 1, a_coarseH.ghostVect());
//   interpolator.pwcinterpToFine(pwcH , a_coarseH);
  
//   Real densityRatio = a_waterDensity / a_iceDensity;

//   Real smallH = 1.0; 
//   // floating ice will be no thicker than fH - smallH

//   for (DataIterator dit(a_newDBL); dit.ok(); ++dit)
//     {
//       for (BoxIterator bit(a_newDBL[dit]); bit.ok(); ++bit)
// 	{
// 	  const IntVect iv = bit();
// 	  Real pwcFH = (a_seaLevel - pwcBed[dit](iv))*densityRatio;
// 	  Real fH = (a_seaLevel - bed[dit](iv))*densityRatio;
// 	  if (pwcH[dit](iv) <  pwcFH)
// 	    {
// 	      //floating ice on the coarse grid
// 	      if ( !(a_H[dit](iv) < fH))
// 		//grounded ice on the fine grid
// 		{ 
// 		  a_H[dit](iv) = fH - smallH ;
// 		}
// 	    }
// 	  else if (pwcH[dit](iv) >  pwcFH)
// 	    {
// 	      //floating ice on the coarse grid
// 	      if ( !(a_H[dit](iv) > fH))
// 		//grounded ice on the fine grid
// 		{		  
// 		  a_H[dit](iv) = fH + smallH ; 
// 		}
// 	    }

	  
// 	}
//     }
  
// }



#if BISICLES_Z == BISICLES_LAYERED

/// update the flow law coefficient A
void AmrIce::computeA(Vector<LevelData<FArrayBox>* >& a_A, 
		      Vector<LevelData<FArrayBox>* >& a_sA,
		      Vector<LevelData<FArrayBox>* >& a_bA,
		      const Vector<LevelData<FArrayBox>* >& a_temperature, 
		      const Vector<LevelData<FArrayBox>* >& a_sTemperature,
		      const Vector<LevelData<FArrayBox>* >& a_bTemperature,
		      const Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys) const
		      
{
  if (s_verbosity > 0)
    {
      pout() <<  "AmrIce::computeA" <<  std::endl;
    }

  //for now, throw a_A etc away and recompute
  for (int lev = 0; lev < a_A.size(); ++lev)
    {
      if (a_A[lev] != NULL)
	{
	  delete a_A[lev]; a_A[lev] = NULL;
	}
      
      if (a_sA[lev] != NULL)
	{
	  delete a_sA[lev]; a_sA[lev] = NULL;
	}
      if (a_bA[lev] != NULL)
	{
	  delete a_bA[lev]; a_bA[lev] = NULL;
	}
      
    }
  a_A.resize(m_finest_level+1,NULL);
  a_sA.resize(m_finest_level+1,NULL);
  a_bA.resize(m_finest_level+1,NULL);
	
  for (int lev = 0; lev <= m_finest_level; ++lev)
    {
      
      a_A[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], m_nLayers, IntVect::Unit);
      a_sA[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1, IntVect::Unit);
      a_bA[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1, IntVect::Unit);
      
      const LevelSigmaCS& levelCoords = *a_coordSys[lev];
      const Vector<Real>& sigma = levelCoords.getSigma();
      for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
	{
	  // compute A(T)
	  // need a temperature field corrected to the pressure melting point,
	  // \theta ^* = \min(\theta,\theta _r) + a * p)
	  // a is a constant, p is pressure, \thera _r is the melting point in standard pressure
	  // using p = \rho * g * \sigma * H 
          // (used by Glimmer, even with higher order stresses)
	  // should be p = T_xx + T_yy + \rho * g * \sigma * H
	  Real Tmax = triplepoint - TINY_NORM;
	  Real fbase = levelCoords.iceDensity() * levelCoords.gravity() * icepmeltfactor;
	  const FArrayBox& theta = (*a_temperature[lev])[dit];
	  // if (s_verbosity > 0)
    // 	    {
    // 	      Real Amin = theta.min();
    // 	      Real Amax = theta.max();
    // 	      pout() << Amin << " <= theta(x,y,sigma) <= " << Amax << std::endl;
	      
    // }
	  for (int layer = 0; layer < m_nLayers; ++layer)
	    {
	      Real f = fbase * sigma[layer];
	      FArrayBox layerA;
	      layerA.define(Interval(layer,layer),(*a_A[lev])[dit]);
	      const Box& box = layerA.box();
	      FArrayBox thetaStar(box,1);
	      
	      

	      FORT_FABMINPLUS(CHF_FRA1(thetaStar,0),
			      CHF_FRA1(theta,layer),
			      CHF_FRA1(levelCoords.getH()[dit],0),
			      CHF_CONST_REAL(f),
			      CHF_CONST_REAL(Tmax),
			      CHF_BOX(box));
	

	      CH_assert(0.0 < thetaStar.min(box));
	      CH_assert(thetaStar.max(box) < triplepoint); 

	      m_rateFactor->computeA(layerA,
				     thetaStar, 
				     layerA.box());




	    } // end loop over layers
	  
	  //surface
	  {
	    m_rateFactor->computeA((*a_sA[lev])[dit],
				   (*a_sTemperature[lev])[dit] , 
				   (*a_sA[lev])[dit].box());
	  }
	  //base
	  {
	    const Box& box = (*m_bA[lev])[dit].box();
	    FArrayBox thetaStar((*a_bA[lev])[dit].box() ,1);
	    FORT_FABMINPLUS(CHF_FRA1(thetaStar,0),
			    CHF_FRA1((*a_bTemperature[lev])[dit],0),
			    CHF_FRA1(levelCoords.getH()[dit],0),
			    CHF_CONST_REAL(fbase),
			    CHF_CONST_REAL(Tmax),
			    CHF_BOX(box));
	    m_rateFactor->computeA((*a_bA[lev])[dit],
				   thetaStar , box);
	  }
	} // end loop over boxes
      
    }//end loop over AMR levels

  if (s_verbosity > 0)
    {
      Real Amin = computeMin(a_A,  m_refinement_ratios, Interval(0,a_A[0]->nComp()-1));
      Real Amax = computeMax(a_A,  m_refinement_ratios, Interval(0,a_A[0]->nComp()-1));
      pout() << Amin << " <= A(x,y,sigma) <= " << Amax << std::endl;

    }
}


//compute the temperature field and the bulk dissipation at the half-time step
void AmrIce::computeTHalf(Vector<LevelData<FluxBox>* >& a_layerTH_half,
			  Vector<LevelData<FluxBox>* >& a_layerH_half,
			  const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel, 
			  const Real a_dt, const Real a_time)
{

  
  //first fill the ghost regions of m_temperature
  for (int lev = 1 ; lev <= m_finest_level; lev++)
    {
      PiecewiseLinearFillPatch tempFiller(m_amrGrids[lev],
					  m_amrGrids[lev-1],
					  m_temperature[lev]->nComp(),
					  m_amrDomains[lev-1],
					  m_refinement_ratios[lev-1],
					  m_temperature[lev]->ghostVect()[0]);
      tempFiller.fillInterp(*m_temperature[lev],*m_temperature[lev-1],
			    *m_temperature[lev-1],1.0,0,0,m_temperature[lev]->nComp());
			    
    }


  //in the 2D case (with poor man's multidim) this
  //is a little pained using AdvectPhysics, but for the time being
  //we need to construct a single component thisLayerT_Half for each layer, 
  //given a temperature and horizontal velocity and then copy it into 
  //the multicomponent T_half[lev]
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {

      if (a_layerTH_half[lev] != NULL)
	delete(a_layerTH_half[lev]);

      a_layerTH_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 
					     m_temperature[lev]->nComp(), 
					     IntVect::Unit);

      if (a_layerH_half[lev] != NULL)
	delete(a_layerH_half[lev]);
      a_layerH_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 
					     m_temperature[lev]->nComp(), 
					     IntVect::Unit);

      
      PatchGodunov* patchGoduTPtr = m_temperaturePatchGodVect[lev];
      AdvectPhysics* advectPhysTPtr = dynamic_cast<AdvectPhysics*>(patchGoduTPtr->getGodunovPhysicsPtr());
      if (advectPhysTPtr == NULL)
	{
	  MayDay::Error("AmrIce::timestep -- unable to upcast GodunovPhysics to AdvectPhysics");
	}
      patchGoduTPtr->setCurrentTime(m_time);
      
      const LevelData<FArrayBox>& levelTemperature = *m_temperature[lev]; 
      const LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev]; 
      const LevelData<FluxBox>& levelLayerXYFaceXYVel = *a_layerXYFaceXYVel[lev]; 
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      for (int layer = 0; layer < m_nLayers; ++layer)
	{

	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      const Box& box = levelTemperature[dit].box(); // grid box plus ghost cells
	      
	      FluxBox layerXYFaceXYVel(box,1);
	      for (int dir = 0; dir < SpaceDim; ++dir){
		layerXYFaceXYVel[dir].copy(levelLayerXYFaceXYVel[dit][dir],layer,0,1);
	      }

	      FArrayBox layerCellXYVel(box,SpaceDim);
	      EdgeToCell(layerXYFaceXYVel,layerCellXYVel);

	      //\todo compute bulk heat sources
	      FArrayBox heatSource(levelGrids[dit], 1);
	      heatSource.setVal(0.0);
	      patchGoduTPtr->setCurrentBox(levelGrids[dit]);
	      

	      advectPhysTPtr->setVelocities(&layerCellXYVel,&layerXYFaceXYVel);
	      FArrayBox WGdnv(box,1);

	      //HT at half time and cell faces
	      WGdnv.copy(levelTemperature[dit],layer,0,1);
	      WGdnv *= levelOldThickness[dit];
	      //CH_assert(0.0 < WGdnv.min(levelGrids[dit]));
	      //CH_assert(WGdnv.max(levelGrids[dit]) < levelOldThickness[dit].max()*triplepoint);
	      Box grownBox = levelGrids[dit];
	      grownBox.grow(1);
	      FluxBox HThalf(grownBox,1);
	      patchGoduTPtr->computeWHalf(HThalf,
					  WGdnv,
					  heatSource,
					  a_dt,
					  levelGrids[dit]);
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		  Box faceBox(levelGrids[dit]);
		  faceBox.surroundingNodes(dir);
		  CH_assert(HThalf[dir].norm(faceBox,0) < HUGE_NORM);
		  (*a_layerTH_half[lev])[dit][dir].copy(HThalf[dir],0,layer,1);
		}
	      
	      //H at half time and cell faces
	      WGdnv.copy(levelOldThickness[dit]);
	      FluxBox Hhalf(grownBox,1);
	      patchGoduTPtr->computeWHalf(Hhalf,
					  WGdnv,
					  heatSource,
					  a_dt,
					  levelGrids[dit]);
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		  Box faceBox(levelGrids[dit]);
		  faceBox.surroundingNodes(dir);
		  CH_assert(Hhalf[dir].norm(faceBox,0) < HUGE_NORM);
		  (*a_layerH_half[lev])[dit][dir].copy(Hhalf[dir],0,layer,1);
		}
	      
	    }
	  
	}
    }
      
  // coarse average new T-Half to covered regions
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],
				     a_layerTH_half[lev]->nComp(), m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*a_layerTH_half[lev-1], *a_layerTH_half[lev]);
      faceAverager.averageToCoarse(*a_layerH_half[lev-1], *a_layerH_half[lev]);
    }

}

void AmrIce::updateTemperature(Vector<LevelData<FluxBox>* >& a_layerTH_half, 
			       Vector<LevelData<FluxBox>* >& a_layerH_half,
			       const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			       const Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel,
			       const Real a_dt, const Real a_time,
			       Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysNew,
			       Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysOld,
			       const Vector<LevelData<FArrayBox>*>& a_surfaceThicknessSource,
			       const Vector<LevelData<FArrayBox>*>& a_basalThicknessSource)
{

  //update the temperature fields, 2D case
  Vector<LevelData<FluxBox>* > vectLayerFluxes(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      LevelData<FluxBox>& levelXYFaceXYVel = *a_layerXYFaceXYVel[lev];
      LevelData<FluxBox>& levelFaceTH = *a_layerTH_half[lev];
      
      IntVect ghostVect = IntVect::Unit;//CoarseAverageFace requires a ghost cell
      vectLayerFluxes[lev] = new LevelData<FluxBox>
	(m_amrGrids[lev], levelXYFaceXYVel.nComp() , ghostVect);
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  FluxBox& faceVel = levelXYFaceXYVel[dit];
	  FluxBox& faceTH = levelFaceTH[dit];
	  FluxBox& flux = (*vectLayerFluxes[lev])[dit];
	  const Box& gridBox = levelGrids[dit];
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      Box faceBox(gridBox);
	      faceBox.surroundingNodes(dir);
	      flux[dir].copy(faceTH[dir], faceBox);
	      flux[dir].mult(faceVel[dir], faceBox, 0, 0, faceVel[dir].nComp());
	      CH_assert(flux[dir].norm(faceBox,0) < HUGE_NORM);
	    }
	}
    }
  // average fine fluxes down to coarse levels
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],
				     vectLayerFluxes[lev]->nComp(), m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*vectLayerFluxes[lev-1], *vectLayerFluxes[lev]);
    }

      

  // compute rhs =a_dt *(H*dissipation - div(u H T)) and update solution
  for (int lev=0; lev<= m_finest_level; lev++)
    {
	  
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFlux = *vectLayerFluxes[lev];
      LevelSigmaCS& levelCoordsNew = *(a_coordSysNew[lev]);
      LevelSigmaCS& levelCoordsOld = *(a_coordSysOld[lev]);
      LevelData<FArrayBox>& levelOldT = *m_temperature[lev];
      LevelData<FArrayBox> dissipation(levelGrids,m_nLayers,IntVect::Zero);
      //dissipation.define(levelOldT);
      //we really ought to compute this with T_half, but for now..
      LevelData<FArrayBox>* crseVelPtr = NULL;
      int nRefCrse = -1;
      if (lev > 0)
        {
          crseVelPtr = m_velocity[lev-1];
          nRefCrse = m_refinement_ratios[lev-1];
        }
      m_constitutiveRelation->computeDissipation
	(dissipation,*m_velocity[lev],  crseVelPtr,
	 nRefCrse, *m_A[lev],
	 levelCoordsOld , m_amrDomains[lev], IntVect::Zero);
						 

      // grad(H) will be needed to evaulate metric terms Hg^{i\sigma},
      // perhaps this ought to live in LevelSigmaCS
	  
      LevelData<FArrayBox> levelGradHNew(m_amrGrids[lev], SpaceDim, IntVect::Zero);
      computeCCDerivatives(levelGradHNew, levelCoordsNew.getH(), levelCoordsNew,
			   Interval(0,0),Interval(0,SpaceDim-1));
      LevelData<FArrayBox> levelGradHOld(m_amrGrids[lev], SpaceDim, IntVect::Zero);
      computeCCDerivatives(levelGradHOld, levelCoordsOld.getH(), levelCoordsOld,
			   Interval(0,0),Interval(0,SpaceDim-1));

      const RealVect& dx = levelCoordsNew.dx(); 
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{

	  // first, do the ordinary fluxes : if we just had
	  // horizontal advection and grad(H) = grad(S) = 0., 
	  // this would be the lot
	       
	  const Box& gridBox = levelGrids[dit];
	  FArrayBox& oldT = levelOldT[dit];
	  FArrayBox& oldH = levelCoordsOld.getH()[dit];
	  FArrayBox& newH = levelCoordsNew.getH()[dit];
	  FArrayBox Hhalf(gridBox, 1);
	  Hhalf.copy(oldH);
	  Hhalf.plus(newH);
	  Hhalf *= 0.5;

	  FArrayBox rhs(gridBox, oldT.nComp());
	  rhs.setVal(0.0);

	  // horizontal fluxes across faces with constant x or y (usual terms)
	  FluxBox& thisFlux = levelFlux[dit];
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
			      CHF_FRA(rhs),
			      CHF_BOX(gridBox),
			      CHF_CONST_REAL(dx[dir]),
			      CHF_INT(dir));
	     
		
	    }
	  CH_assert(rhs.norm(0) < HUGE_NORM);
	  for (int layer = 0; layer < dissipation.nComp(); ++layer)
	    {
	      dissipation[dit].mult(newH,0,layer,1);
	      dissipation[dit].mult(1.0/(iceheatcapacity * levelCoordsNew.iceDensity()));
	      //CH_assert(dissipation[dit].norm(0) < 10.0);
	      //pout() << dissipation[dit].max(layer) 
	      //	     << " > H Phi > " 
	      //	     << dissipation[dit].min(layer) << std::endl;
	    } 

	  rhs -= dissipation[dit]; 
	  rhs *= -a_dt;
	  CH_assert(rhs.norm(0) < HUGE_NORM);
#ifdef TARMAC
	  {
	    CH_assert(SpaceDim == 2); // we'll probably never care about 1D

	    FArrayBox gradH(rhs.box(), SpaceDim);
	    gradH.copy(levelGradHNew[dit]);
	    gradH.plus(levelGradHOld[dit]);
	    gradH*=0.5;
	    
	    FArrayBox gradS(rhs.box(), SpaceDim);
	    gradS.copy(levelCoordsOld.getGradSurface()[dit]);
	    gradS.plus(levelCoordsNew.getGradSurface()[dit]);
	    gradS*=0.5;
	    
	    const FArrayBox& sFaceXYVel = (*a_layerSFaceXYVel[lev])[dit];
		 
	    // this copy perhaps indicates layer should run faster than
	    // dir in sFaceXYVel, but for now ...
	    FArrayBox uX(rhs.box(), m_nLayers+1);
	    FArrayBox uY(rhs.box(), m_nLayers+1);
		 
	    for (int l = 0; l < m_nLayers+1; l++)
	      {
		uX.copy(sFaceXYVel, l*SpaceDim, l);
		uY.copy(sFaceXYVel, l*SpaceDim + 1, l);
	      }

	    const Vector<Real>& dSigma = levelCoordsNew.getDSigma();
	    FArrayBox& sT = (*m_sTemperature[lev])[dit];
	    FArrayBox& bT = (*m_bTemperature[lev])[dit];
	    Real time = a_time; Real dt = a_dt;
	    int n = m_nLayers;
		 
	   

	    //horizontal contribution to div(Hu) at cell centres, 
            // viz d(Hu_x)/dx' + d(Hu_y)/dy'
	    FArrayBox divUHxy(rhs.box(), m_nLayers);
	    divUHxy.setVal(0.0);
	   
	    
	    const FluxBox& thisU = (*a_layerXYFaceXYVel[lev])[dit];
	    const FluxBox& faceHNew =levelCoordsNew.getFaceH()[dit];
	    const FluxBox& faceHOld =levelCoordsOld.getFaceH()[dit];
	      
	    for (int dir =0; dir < SpaceDim; dir++)
	      {
		FArrayBox uH(thisU[dir].box(),m_nLayers);
		for (int l = 0; l < m_nLayers; ++l)
		  {
		    uH.copy( (*a_layerH_half[lev])[dit][dir],l,l);
		    
		    //uH.copy(faceHOld[dir],0,l);
		    //uH.plus(faceHNew[dir],0,l);
		  }
		//uH *= 0.5;
		uH *= thisU[dir];
		FORT_DIVERGENCE(CHF_CONST_FRA(uH),
				CHF_FRA(divUHxy),
				CHF_BOX(rhs.box()),
				CHF_CONST_REAL(dx[dir]),
				CHF_INT(dir));
		
	
	      }
	  
	    

	    //basal thickness source
	    FArrayBox& bts = (*a_basalThicknessSource[lev])[dit];
	    //surface thickness source
	    const FArrayBox& sts = (*a_surfaceThicknessSource[lev])[dit];

	    FArrayBox sumdivUH(rhs.box(),1);
	   
	    sumdivUH.setVal(0.0);
	   
	    for (int l = 0; l < m_nLayers; ++l){
	      sumdivUH.plus(divUHxy,dSigma[l],l,0,1);
	      
	    }
	    

	     
	    FArrayBox dHdta(rhs.box(),1);  
	    dHdta.copy(newH);
	    dHdta.plus(oldH,-1.0,0,0,1);
	    dHdta *= 1.0/a_dt;

	    //sumdivUH += dHdt;
	    
	    
	    FArrayBox dHdtb(rhs.box(),1); 
	    dHdtb.copy(sumdivUH);
	    dHdtb*=-1;
	    dHdtb+=sts;
	    dHdtb+=bts;
	      
	    FArrayBox tmp(rhs.box(),1);
	    tmp.copy(dHdta);
	    tmp -= dHdtb;

	    pout() << "err dHdt = " <<  tmp.norm(0) << std::endl;

	    //calculation of dSdt assumes surface elevation is up to date
            //in LevelSigmaCS
	    FArrayBox dSdt(rhs.box(),1); 
	    dSdt.copy(levelCoordsNew.getSurfaceHeight()[dit]);
	    dSdt -= levelCoordsOld.getSurfaceHeight()[dit];
	    dSdt *= 1.0/a_dt;

	    // z-component of velocity at layer faces
	    FArrayBox uZ(rhs.box(),m_nLayers + 1); 
	    // sigma-component of velocity at layer faces
	    FArrayBox uSigma(rhs.box(),m_nLayers + 1); 
	    // z-component of velocity at surface (by bc)
	    FArrayBox uZs(rhs.box(), 1);

	    FORT_COMPUTEZVEL(CHF_FRA(uZ),
			     CHF_FRA1(uZs,0),
			     CHF_FRA(uSigma),
			     CHF_CONST_FRA(uX),
			     CHF_CONST_FRA(uY),
			     CHF_CONST_FRA(divUHxy),
			     CHF_CONST_VR(levelCoordsNew.getFaceSigma()),
			     CHF_CONST_VR(levelCoordsNew.getSigma()),
			     CHF_CONST_VR(dSigma),
			     CHF_CONST_FRA1(Hhalf,0),
			     CHF_CONST_FRA1(gradS,0), 
			     CHF_CONST_FRA1(gradH,0),
			     CHF_CONST_FRA1(gradS,1), 
			     CHF_CONST_FRA1(gradH,1),
			     CHF_CONST_FRA1(dSdt,0), 
			     CHF_CONST_FRA1(dHdta,0),
			     CHF_CONST_FRA1(sts,0),
			     CHF_CONST_FRA1(bts,0),
			     CHF_INT(n),
			     CHF_BOX(rhs.box()));

	    
	     if (s_verbosity > 3)
	       {
		 FArrayBox differ(rhs.box(), 1);
		 differ.copy(uZs);
		 differ.minus(uZ,0,0,1);
		 pout() 
		   << " max | u_z( integrated)| = " 
		   << uZ.norm(0)
		   << " max | u_z(surface, integrated) - u_z(surface, bc)| = " 
		   << differ.norm(0) 
		   << " max | u^sigma( integrated)| = " 
		   << uSigma.norm(0)
		   << std::endl;
	       }

	     //compute heat flux across base due to basal dissipation & geothermal heat 
	     FArrayBox basalHeatFlux(rhs.box(),1);
	     m_basalFrictionRelation->computeDissipation
	       (basalHeatFlux , (*m_velocity[lev])[dit] , (*m_velBasalC[lev])[dit], rhs.box());
	     CH_assert(0.0 <=  basalHeatFlux.min() &&  basalHeatFlux.max() < HUGE_NORM);
	     basalHeatFlux /= (iceheatcapacity * levelCoordsNew.iceDensity());
	     pout() << basalHeatFlux.max() << " > bflux > " << basalHeatFlux.min() << std::endl;
	     // update the layer, surface, and base temperature 
	     FORT_UPDATETEMPERATURE(CHF_FRA(oldT), 
				    CHF_FRA1(sT,0), 
				    CHF_FRA1(bT,0),
				    CHF_FRA1(basalHeatFlux,0),
				    CHF_CONST_FIA1(levelCoordsNew.getFloatingMask()[dit],0),
				    CHF_CONST_FRA(rhs),
				    CHF_CONST_FRA1(oldH,0),
				    CHF_CONST_FRA1(newH,0),
				    CHF_CONST_FRA(uSigma),
				    CHF_CONST_VR(levelCoordsNew.getFaceSigma()),
				    CHF_CONST_VR(dSigma),
				    CHF_REAL(time), 
				    CHF_REAL(dt),
				    CHF_INT(n),
				    CHF_BOX(rhs.box()));

	     CH_assert(0.0 < oldT.min(rhs.box()) &&  oldT.max(rhs.box()) < triplepoint); 


	  }
#endif

	}
    }
 
   if ( m_basalLengthScale > TINY_NORM ){
  //smoothing...
     helmholtzSolve(m_temperature,Real(1.0),std::pow(m_basalLengthScale,2));
   }
     
   //coarse average from finer levels & exhange
   for (int lev = m_finest_level; lev >= 0 ; --lev)
   {
     if (lev > 0)
       {
	 CoarseAverage an(m_amrGrids[lev],
			  m_temperature[lev]->nComp(),
			  m_refinement_ratios[lev-1]);
	 
	 an.averageToCoarse(*m_temperature[lev-1], *m_temperature[lev]);
	 
	 CoarseAverage aone(m_amrGrids[lev],
			    1,m_refinement_ratios[lev-1]);

	 aone.averageToCoarse(*m_sTemperature[lev-1], *m_sTemperature[lev]);
	 aone.averageToCoarse(*m_bTemperature[lev-1], *m_bTemperature[lev]);

       }


     m_temperature[lev]->exchange();
     m_sTemperature[lev]->exchange();
     m_bTemperature[lev]->exchange();
   }

  for (int lev = 0; lev < vectLayerFluxes.size(); ++lev)
    {
      

      if (vectLayerFluxes[lev] != NULL)
	{
	  delete vectLayerFluxes[lev];
	  vectLayerFluxes[lev] = NULL;
	}
    }

  //finally, A is no longer valid 
  m_A_valid = false;
  

}
#endif

#include "NamespaceFooter.H"
