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
using std::string;

#include "AmrIce.H"
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
#include "EBAmrIce.H"
#include "computeNorm.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "PiecewiseLinearFillPatch.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "DerivativesF_F.H"
#include "DivergenceF_F.H"
#include "computeSum.H"
#include "CONSTANTS.H"
#include "IceConstants.H"
#include "ExtrapBCF_F.H"
#include "amrIceF_F.H"
#include "EBAmrIceF_F.H"
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

#include "GroundingLineIF.H"
#include "HyperPlaneIF.H"

#ifdef CH_USE_EB
#include "NewGeometryShop.H"
#include "EBIndexSpace.H"
#else
#include "Notation.H"
#include "RootFinder.H"
#include "BaseIF.H"
#include "DataFileIF.H"
#endif

#include "areaWeightedFriction.H"

#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

#include "NamespaceHeader.H"
void  EBAmrIce::define(ConstitutiveRelation * a_constRelPtr,
                       RateFactor           * a_rateFactor,
                       SurfaceFlux          * a_surfFluxPtr,
                       SurfaceFlux          * a_basalFluxPtr,
                       MuCoefficient        * a_muCoefficientPtr,
                       BasalFrictionRelation* a_basalFrictionRelationPtr,
                       BasalFriction        * a_basalFrictionPtr,
                       IceThicknessIBC      * a_thicknessIBC,
                       IceTemperatureIBC    * a_temperatureIBC,
                       const RealVect       & a_domainSize)
{
  setConstitutiveRelation (a_constRelPtr);
  setRateFactor           (a_rateFactor);
  setSurfaceFlux          (a_surfFluxPtr);
  setBasalFlux            (a_basalFluxPtr); 
  setMuCoefficient        (a_muCoefficientPtr);
  setBasalFrictionRelation(a_basalFrictionRelationPtr);
  setBasalFriction        (a_basalFrictionPtr);
  setThicknessBC          (a_thicknessIBC);
  setTemperatureBC        (a_temperatureIBC);
  setDomainSize           (a_domainSize);
  
  initialize();
}

/// define nonlinear ellipic solver for computing velocity field
void EBAmrIce::defineSolver()
{
  AmrIce::defineSolver();
}

/// solve for velocity field (or just set up some auxilliary quantities)
void EBAmrIce::velInitialGuess(Vector<LevelData<FArrayBox>* >& a_vectC0,
                               const Real                    & a_convergenceMetric)
{
 if (m_finest_level == 0 && m_doInitialVelGuess)
   {
     // only need to do this once
     m_doInitialVelGuess = false;
              
     if (m_initialGuessType == SlidingLaw)
       {
         pout() << "computing an initial guess via a sliding law u = rhs/C "  << endl;
                  
         // compute initial guess as rhs/beta
         LevelData<FArrayBox>& vel = *m_velocity[0];
         LevelData<FArrayBox>& C   = *m_velBasalC[0];
         LevelData<FArrayBox>& rhs = *m_velRHS[0];
                  
         // dbl
         const DisjointBoxLayout& levelGrids = m_amrGrids[0];
                  
         DataIterator dit = vel.dataIterator();
         for (dit.begin(); dit.ok(); ++dit)
           {
             FORT_VELINITIALGUESS(CHF_FRA (vel       [dit]),
                                  CHF_FRA (rhs       [dit]),
                                  CHF_FRA1(C         [dit],0),
                                  CHF_BOX (levelGrids[dit]));
           }
       }
     else if (m_initialGuessType == ConstMu)
       {
         // compute initial guess by solving a linear problem with modest constant viscosity
         if (s_verbosity > 3) 
           {
             pout() << "computing an initial guess by solving the velocity equations "
                    <<" with constant mu = " 
                    << m_initialGuessConstMu   
                    << "and constant initial velcocity = " << m_initialGuessConstVel
                    << endl;
           }
              
         // dx on the coarsest level
         RealVect dxCrse = m_amrDx[0]*RealVect::Unit;
                  
         // viscosity
         constMuRelation* newPtr = new constMuRelation;
         newPtr->setConstVal(m_initialGuessConstMu);
                  
         // set velocity to m_initialGuessConstVel
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
                  
         // cache velSolver and constitutiveRelation
         IceVelocitySolver   * velSolverSave = m_velSolver;
         ConstitutiveRelation* constRelSave  = m_constitutiveRelation;
         int solverTypeSave                  = m_solverType;
                  
         // temp value
         m_constitutiveRelation = static_cast<ConstitutiveRelation*>(newPtr);
                  
         //solver parameters
         Real finalNorm         = 0.0; 
         Real initialNorm       = 0.0;
         Real convergenceMetric = -1.0;
                  
         // mucoef
         Vector<LevelData<FluxBox>* > muCoef(m_finest_level + 1,NULL);
                  
         // return condition
         int rc;
                  
         if (m_initialGuessSolverType == JFNK)
           {
             //JFNK can be instructed to assume a linear solve
             m_velSolver = NULL;
                      
             // define solver
             defineSolver();
             JFNKSolver* jfnkSolver = dynamic_cast<JFNKSolver*>(m_velSolver);
             CH_assert(jfnkSolver != NULL);
                      
             const bool linear = true;
             rc = jfnkSolver->solve(m_velocity, 
                                    initialNorm,
                                    finalNorm,
                                    convergenceMetric,
                                    linear, 
                                    m_velRHS, 
                                    m_velBasalC, 
                                    a_vectC0, 
                                    m_A, 
                                    muCoef,
                                    m_vect_coordSys, 
                                    m_time, 
                                    0, 
                                    m_finest_level);
           }
         else if (m_initialGuessSolverType == Picard)
           {
             // since the constant-viscosity solve is a linear solve, Picard is the best option.
             m_solverType = Picard;
             m_velSolver  = NULL;
                      
             // define solver
             defineSolver();

#if CH_SPACEDIM== 1
#if CUH_SPACEDIM == 1     
             m_velSolver->setGroundingLineData(m_groundingLineIv,
                                               m_physCoordGroundingPt,
                                               m_lengthFraction);
#endif
#endif
             // solve
             rc = m_velSolver->solve(m_velocity, 
                                     initialNorm,
                                     finalNorm,
                                     convergenceMetric,
                                     m_velRHS, 
                                     m_velBasalC, 
                                     a_vectC0, 
                                     m_A, 
                                     muCoef,
                                     m_vect_coordSys, 
                                     m_time, 
                                     0, 
                                     m_finest_level);
           }
         else
           {
             MayDay::Error("unknown initial guess solver type");
           }
                  
                  
         if (rc != 0)
           {
             MayDay::Warning("constant mu solve failed");
           }
                  
         // Put everything back the way it was for solves after the initial solve
         delete m_constitutiveRelation;
         delete m_velSolver;
                  
         // member data
         m_velSolver            = velSolverSave;
         m_constitutiveRelation = constRelSave;
         m_solverType           = solverTypeSave;
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
             pp.get("initial_velocity_function_flowline_dx"  , dx);
             pp.get("initial_velocity_function_flowline_file", file);
             pp.get("initial_velocity_function_flowline_set" , set);
                      
             // function for initial guess
             ExtrudedPieceWiseLinearFlowline f(file,set,dx);
                      
             // iterate over every grid location, every level setting velocity = f(x)
             for (int lev = 0; lev <  m_finest_level + 1; lev++)
               {
                 LevelData<FArrayBox>& levelVel = *m_velocity[lev];
                 for (DataIterator dit(levelVel.disjointBoxLayout());dit.ok(); ++dit)
                   {
                     // velocity
                     FArrayBox& vel = levelVel[dit];
                              
                     //iterate over box
                     const Box& box = vel.box(); 
                     for (BoxIterator bit(box); bit.ok(); ++bit)
                       {
                         // IntVect
                         const IntVect& iv = bit();
                                  
                         // physical coordinates with origin hard-wired to 0
                         RealVect x = RealVect(iv) * m_amrDx[lev] + 0.5 * m_amrDx[lev];
                                  
                         // velocity = f(x)
                         vel(iv,0) = f(x);
                       }
                   }
               }
           }
                  
       }
     else
       {
         MayDay::Error("EBAmrIce::SolveVelocityField unknown initial guess type");
       }
   }
          
 // I/O
 if (m_write_presolve_plotfiles)
   {
     string save_prefix = m_plot_prefix;
     m_plot_prefix.append("preSolve.");
              
     bool t_write_fluxVel = m_write_fluxVel;
              
     // turning this off in preSolve files 
     m_write_fluxVel = false; 
              
     // write
     writePlotFile();
              
     // assign member data
     m_write_fluxVel = t_write_fluxVel;
     m_plot_prefix   = save_prefix;
   }
          
 // set velocity to zero over open sea or open land
 for (int lev=0; lev <= m_finest_level ; ++lev)
   {
     const DisjointBoxLayout& levelGrids = m_amrGrids      [lev];
     LevelSigmaCS           & levelCS    = *m_vect_coordSys[lev];
              
     for (DataIterator dit(levelGrids); dit.ok(); ++dit)
       {
         const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
         FArrayBox         & vel  = (*m_velocity[lev])       [dit];
                  
         for (BoxIterator bit(levelGrids[dit]); bit.ok(); ++bit)
           {
             const IntVect& iv = bit();
             if (mask(iv) == OPENSEAMASKVAL || mask(iv) == OPENLANDMASKVAL )
               {
                 vel(iv,0) = 0.0; 
                 vel(iv,1) = 0.0;
               } 
           }
       }
   }
}

/// solve for velocity field (or just set up some auxilliary quantities)
void EBAmrIce::solveVelocityField(Real a_convergenceMetric)
{
  CH_TIME("EBAmrIce::solveVelocityField");
  
  createEBIndexSpace();

  //ensure A is up to date
#if BISICLES_Z == BISICLES_LAYERED
  if (!m_A_valid)
    {
      computeA(m_A,
               m_sA,
               m_bA,
               m_temperature,
               m_sTemperature,
               m_bTemperature,
               m_vect_coordSys);
      m_A_valid = true;
    }
#else
  MayDay::Error("EBAmrIce::SolveVelocityField full z calculation of A not done"); 
#endif
  
  // viscous tensor field will need updating at the conclusion of solveVelocityField
  m_viscousTensor_valid = false;

  // basal friction
  Vector<LevelData<FArrayBox>* > vectC  (m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectC0 (m_finest_level+1, NULL);

  // right hand side
  Vector<LevelData<FArrayBox>* > vectRhs(m_finest_level+1, NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      vectRhs[lev] = m_velRHS[lev];
      vectC  [lev] = m_velBasalC[lev];

      vectC0 [lev] = new LevelData<FArrayBox>; 
      vectC0 [lev]->define(*vectC[lev]);

      // set vectC0 to zero
      LevelData<FArrayBox>& C0 = *vectC0[lev];
      for (DataIterator dit    = C0.dataIterator(); dit.ok(); ++dit)
	{
	  C0[dit].setVal(0.0);
	}
    }

  setMuCoefficient(m_cellMuCoef,
                   m_faceMuCoef);

  setBasalFriction(vectC);

  // also sets beta = 0 where ice is floating 
  defineVelRHS(vectRhs, 
               vectC,
               vectC0);
  
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
  
  // catch runs where plotfile writing is going to hang
#define writeTestPlot
#ifdef writeTestPlot
  if (m_plot_interval >= 0 && m_cur_step == 0)
    {
      writePlotFile();
    }
#endif

  if (m_doInitialVelSolve) 
    {
      if(m_doInitialVelGuess)
        {
          velInitialGuess(vectC0,
                          a_convergenceMetric);
        }

#ifndef CH_USE_EB
#if CH_SPACEDIM== 1
      m_velSolver->setGroundingLineData(m_groundingLineIv,
                                        m_physCoordGroundingPt,
                                        m_lengthFraction);
#endif
#endif
      // multi-fluid viscous tensor solve
      int solverRetVal = m_velSolver->solve(m_velocity,
                                            m_velocitySolveInitialResidualNorm, 
                                            m_velocitySolveFinalResidualNorm,
                                            a_convergenceMetric,
                                            m_velRHS, 
                                            m_velBasalC, 
                                            vectC0,
                                            m_A, 
                                            m_faceMuCoef,
                                            m_vect_coordSys,
                                            m_time,
                                            0, 
                                            m_finest_level);
          if (solverRetVal != 0)
            {
              MayDay::Abort("velocity solve failed");
            }
    }

  //calculate the face centered velocity and diffusion coefficients
  computeFaceVelocity(m_faceVelAdvection,
                      m_faceVelTotal,
                      m_diffusivity,
                      m_layerXYFaceXYVel, 
                      m_layerSFaceXYVel);
      
  for (int lev=0; lev <= m_finest_level; lev++)
    {
      if (vectC0[lev] != NULL)
        {
          delete vectC0[lev]; 
          vectC0[lev] = NULL;
        }
    }
}

/// compute RHS for velocity field solve
/** also sets beta to zero where ice is floating
 */
void EBAmrIce::defineVelRHS(Vector<LevelData<FArrayBox>* >& a_vectRhs,
                            Vector<LevelData<FArrayBox>* >& a_vectC,
                            Vector<LevelData<FArrayBox>* >& a_vectC0)
{
  // fill ghost cells for H  
  Vector<LevelData<FArrayBox>* >tempH(m_finest_level+1);
  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      const DisjointBoxLayout& levelGrids  = m_amrGrids      [lev];
      LevelSigmaCS           & levelCoords = *m_vect_coordSys[lev];
     
      // use levelH to call fillInterp
      tempH[lev]                   = (&levelCoords.getH());
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
	  ghostFiller.fillInterp(levelH, 
                                 coarseH, 
                                 coarseH, 
                                 0.0, 
                                 0, 
                                 0, 
                                 1);
	}
      
      // exchange
      levelH.exchange();
    }

  // multilevel surface height
  Vector<LevelData<FArrayBox>* > zSurf(m_finest_level+1, NULL);
  
  // use levelCS to read surface height into zSurf
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const DisjointBoxLayout& grids   = m_amrGrids[lev];
      LevelSigmaCS           & levelCS = *m_vect_coordSys[lev];
      
      zSurf[lev] = new LevelData<FArrayBox>(grids, 1, IntVect::Unit);
      levelCS.getSurfaceHeight(*zSurf[lev]);
    }
  
  // used for limiting RHS, if needed
  Real maxZs   = 0.0;
  Real maxGrad = -1.0;
  if (m_limitVelRHS && (m_gradLimitRadius > 0) )
    {
      Interval comps(0,0);
      maxZs = computeMax(zSurf, 
                         m_refinement_ratios, 
                         comps, 
                         0);

      maxGrad = maxZs/Real(m_gradLimitRadius);
    }
  
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      // surface
      LevelData<FArrayBox>& levelZs = *zSurf[lev];

      // levelCS parameters
      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const RealVect& dx    = levelCS.dx();
      Real iceDensity       = levelCS.iceDensity();
      Real gravity          = levelCS.gravity();

      // grad surface and depth from levelCS
      const LevelData<FArrayBox>& levelH     = levelCS.getH();
      const LevelData<FArrayBox>& levelGradS = levelCS.getGradSurface();
 
      // mask for open water, open ground, floating, grounded
      const LevelData<BaseFab<int> >& levelMask        = levelCS.getFloatingMask();
      const LayoutData<bool>        & levelAnyFloating = levelCS.anyFloating();

      //friction
      LevelData<FArrayBox>& levelC   = (*a_vectC[lev]);

      // rhs
      LevelData<FArrayBox>& levelRhs = (*a_vectRhs[lev]);
          
      // data iterate
      const DisjointBoxLayout& grids = m_amrGrids[lev];
      DataIterator dit               = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& thisH           = levelH          [dit];
	  const BaseFab<int>& floatingMask = levelMask       [dit];
	  bool  anyFloating                = levelAnyFloating[dit];  
          
          // now take grad(z_s)
          FArrayBox& thisRhs    = levelRhs  [dit];
	  const FArrayBox& grad = levelGradS[dit];
	  const Box& gridBox    = grids     [dit];

          // add in background slope and compute RHS.	  
          // break this into two steps and limit grad(z_s) if necessary
          thisRhs.copy(grad, gridBox);
          
          if (m_limitVelRHS && (m_gradLimitRadius > 0) )
            {
              RealVect maxGradDir;
              maxGradDir[0] = maxGrad / dx[0];
              maxGradDir[1] = maxGrad / dx[1];
              
              // grown dbl box
              Box gridBoxPlus = gridBox;
              gridBoxPlus.grow(1);

              for (BoxIterator bit(gridBox); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  
                  // don't limit grad(s) in floating regions where grad(s) represents the marine boundary condition
                  if  (floatingMask(iv) == GROUNDEDMASKVAL)
                    {
                      for (int dir=0; dir<2; dir++)
                        {
                          if (thisRhs(iv,dir) > maxGradDir[dir])
                            {
                              // thisRhs = grad(s) 
                              thisRhs(iv,dir) = maxGradDir[dir];
                            }
                          else if (thisRhs(iv,dir) < -maxGradDir[dir])
                            {
                              // thisRhs = grad(s) 
                              thisRhs(iv,dir) = -maxGradDir[dir];
                            }
                        }                          
                    }
                }
            }

          for (BoxIterator bit (gridBox); bit.ok(); ++bit)
            {
              IntVect iv = bit();

              for (int idir = 0; idir < SpaceDim; ++idir)
                {     
                  thisRhs(iv,idir) *= iceDensity*gravity*thisH(iv,0);
                } 
            }
          
	  if(m_groundingLineCorrection && anyFloating)
	    {
	      Real rhog = iceDensity*gravity;
	      CH_assert(SpaceDim != 3);
	      for (int dir=0; dir<SpaceDim; dir++)
		{
		  const FArrayBox& thisH = (*tempH[lev])[dit];
		  const FArrayBox& thisZsurf = levelZs[dit];
		  FORT_GLCORRECTION(CHF_FRA1      (thisRhs, dir),
				    CHF_CONST_FRA1(thisH,0),
				    CHF_CONST_FRA1(thisZsurf,0),
				    CHF_CONST_FIA1(floatingMask,0),
				    CHF_INT       (dir),
				    CHF_CONST_REAL(dx[dir]),
				    CHF_CONST_REAL(rhog),
				    CHF_BOX       (gridBox));
		}
	    }

	  // check that rhs not too large
	  CH_assert(thisRhs.norm(0,0,SpaceDim) < HUGE_NORM);
	  
	  FArrayBox& thisC  = levelC          [dit];
	  FArrayBox& thisC0 = (*a_vectC0[lev])[dit];

	  // add drag due to ice in contact with ice-free rocky walls
	  thisC0.setVal(0.0);
	  if (m_wallDrag)
	  {
	    IceVelocity::addWallDrag(thisC0, 
                                     floatingMask, 
                                     levelCS.getSurfaceHeight()[dit],
				     thisH, 
                                     levelCS.getBaseHeight()[dit], 
                                     thisC, 
                                     m_wallDragExtra,
				     RealVect::Unit*m_amrDx[lev], 
                                     gridBox);
	  }
	  
          // set beta to 100 on open land and open sea; set betqa = 0 where ice is floating;
          if(anyFloating)
            {
              FORT_SETOPENLANDOPENSEABETA(CHF_FRA1      (thisC,0),
                                          CHF_CONST_FIA1(floatingMask,0),
                                          CHF_BOX       (gridBox));

	      // friction non-negative
              CH_assert(thisC.min(gridBox) >= 0.0); 
            } 
        } 

      // Modify RHS in problem-dependent ways,
      m_thicknessIBCPtr->modifyVelocityRHS(levelRhs, 
                                           *m_vect_coordSys[lev],
                                           m_amrDomains[lev],
                                           m_time, 
                                           m_dt);
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

/// set basal friction coefficient (beta) prior to velocity solve
void EBAmrIce::setBasalFriction(Vector<LevelData<FArrayBox>* >& a_vectBeta)
{
  AmrIce::setBasalFriction(a_vectBeta);

  // for debugging
  ParmParse pp2("main");
  bool useAreaWeightedFriction;
  pp2.get("useAreaWeightedFriction",useAreaWeightedFriction);
  if (!useAreaWeightedFriction)
    {
      for (int lev=0; lev<=m_finest_level; lev++)
        {
          // geometry information
          LevelSigmaCS& levelCS = *m_vect_coordSys[lev];

          // masks
          const LevelData<BaseFab<int> >& levelMask        = levelCS.getFloatingMask();
          const LayoutData<bool        >& levelAnyFloating = levelCS.anyFloating    ();
          
          // level data
          LevelData<FArrayBox>& levelBeta = (*a_vectBeta[lev]);
          
          // dbl
          const DisjointBoxLayout& dbl = m_amrGrids[lev];
          
          // iterate through dbl
          DataIterator dit = dbl.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              Box dblBox = dbl[dit()];

              // masks for floating ice
              bool                anyFloating  = levelAnyFloating[dit()];  
              const BaseFab<int>& floatingMask = levelMask       [dit()];
  
              // friction
              FArrayBox& beta = levelBeta[dit];
                           
              // set friction on open land and open sea to 100. Set friction on floating ice to 0
              if(anyFloating && m_reset_floating_friction_to_zero)
                {
                  FORT_SETFLOATINGBETA(CHF_FRA1       (beta       ,0),
                                       CHF_CONST_FIA1(floatingMask,0),
                                       CHF_BOX        (dblBox));
                  
                  // friction must be non-negative
                  CH_assert(beta.min(dblBox) >= 0.0); 
                } 
            }
        }
    }
}

void  EBAmrIce::setBasalFriction(const BasalFriction* a_basalFrictionPtr)
{
  AmrIce::setBasalFriction(a_basalFrictionPtr);
}

/// set mu coefficient (phi) prior to velocity solve
void EBAmrIce::setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef, 
                                Vector<LevelData<FluxBox>* >& a_faceMuCoef)
{
  AmrIce::setMuCoefficient(a_cellMuCoef, 
                           a_faceMuCoef);
}

void  EBAmrIce::setMuCoefficient(const MuCoefficient* a_muCoefficientPtr)
{
  AmrIce::setMuCoefficient(a_muCoefficientPtr);
}

void EBAmrIce::createEBIndexSpace()
{
  // index into coarsest amr level
  int baseLevel = 0;

  // level sigma-coordinate-system
  RefCountedPtr<LevelSigmaCS> levelCoordPtr = m_vect_coordSys[baseLevel];
  
  // physical parameters repeated on each level
  Real iceDensity   = levelCoordPtr->iceDensity  ();
  Real waterDensity = levelCoordPtr->waterDensity();
  Real seaLevel     = levelCoordPtr->seaLevel    ();
  
  // level data of ice thickness
  LevelData<FArrayBox>& iceThicknessLD = levelCoordPtr->getH();
  iceThicknessLD.exchange();
  
  // todo: not valid for AMR
  Box iceThicknessBox      = m_amrDomains[baseLevel].domainBox();
  Box grownIceThicknessBox = iceThicknessBox;
  grownIceThicknessBox.grow(2);
  
  // GroundingLineIF and DataIF use ref counted pointers
  RefCountedPtr<FArrayBox> iceThicknessFab(new FArrayBox(grownIceThicknessBox,1));
  
  // copy iceThicknessLD onto iceThicknessFab
  for (DataIterator dit = iceThicknessLD.dataIterator(); dit.ok(); ++dit)
    {
      iceThicknessFab->copy(iceThicknessLD[dit()]);
    }
    
  // bad data value
  Real noDataValue = LARGEREALVAL;
  
  // number of data point in each direction
  IntVect num      = IntVect::Unit + m_amrDomains[baseLevel].domainBox().bigEnd();
  
  // dx for the interface data, which in this case equals dx, defined below
  RealVect spacing = m_amrDx[baseLevel]*RealVect::Unit;
  
  // origin
  setOrigin(IndexTM<Real,SpaceDim>::Zero);
  
  // value of the IF that corresponds to the interface
  Real value       = 0.0;
  
  // which side of the interface is in the calculation
  bool thicknessInside  = true;
  
  // otherwise bi-linear
  bool useCubicInterp = true;
  
  // not using binary data
  RefCountedPtr<BaseFab<unsigned char> > binaryData(NULL);
  
  // data File IF used even though we actually have a state variable IF
  RefCountedPtr<DataFileIF> iceThicknessDataIF(new DataFileIF(iceThicknessFab,
                                                              binaryData,
                                                              noDataValue,
                                                              num,
                                                              spacing,
                                                              m_origin,
                                                              value,
                                                              thicknessInside,
                                                              useCubicInterp));
  
  // data for synthetic "topography" generated by a plane
  ParmParse pp("marineIceSheet");
  
  // slope of plane
  Vector<Real >vecBasalSlope;
  pp.getarr("basal_slope",vecBasalSlope,0,SpaceDim);
  
  // convert vector to RealVect
  IndexTM<Real,SpaceDim> basalSlope;
  for(int idir = 0; idir < SpaceDim; ++idir)
    {
      basalSlope[idir] = vecBasalSlope[idir];
    }
  
  // point on the plane
  Real originElevation;
  pp.get("originElevation",originElevation);
  IndexTM<Real,SpaceDim> planeOrigin;
  CH_assert(basalSlope[0] != 0.0);
  planeOrigin[0] = -originElevation/basalSlope[0];

#if (CH_SPACEDIM == 2)
  {
    planeOrigin[1] = 0.0;
  }
#endif

  // todo: why is inside correct here?
  bool inside = false;
  
  // topography 
  RefCountedPtr<HyperPlaneIF> topographyIF(new HyperPlaneIF(basalSlope, planeOrigin, inside));
  
  // Archimedes principle determines the interface
  GroundingLineIF implicitFunction(iceDensity,
                                   waterDensity,
                                   seaLevel,
                                   topographyIF);
  
  implicitFunction.setDepth(iceThicknessDataIF);
  
  // this won't work with amr
  RealVect dx = RealVect::Unit*m_amrDx[baseLevel];
  
#ifdef CH_USE_EB
  {
    // todo: what constraints make sense
    int order = 1;
    int degree = 2;
    bool useConstraints = false;
    
    ProblemDomain ebProblemDomain(m_amrDomains[baseLevel].domainBox());
    
    // class that uses implicit function to create geometry data for EBIS
    NewGeometryShop workshopPtr(implicitFunction, 
                                m_origin,
                                dx,
                                ebProblemDomain,
                                order,
                                degree,
                                useConstraints);
    
    // todo: do something more clear than hard-wiring these values
    int ebMaxSize      = 32;
    int maxCoarsenings = -1;  
    
    // build EBIS
    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->define(ebProblemDomain, 
                    m_origin, 
                    m_amrDx[baseLevel], 
                    workshopPtr, 
                    ebMaxSize, 
                    maxCoarsenings);
    
  }
#elif (CH_SPACEDIM ==1)
  {
    ParmParse pp2("main");
    bool useAreaWeightedFriction;
    pp2.get("useAreaWeightedFriction",useAreaWeightedFriction);

    if (useAreaWeightedFriction)    
      {
        bool foundGroundingLine = false;
        
        //iterate through every iv
        for (BoxIterator bit(iceThicknessBox); bit.ok() && !foundGroundingLine; ++bit)
          {
            // IntVect
            const IntVect& iv = bit();
            
            // physical coordinates 
            IndexTM<Real,SpaceDim> physCoordIv;
            for (int idir = 0; idir <SpaceDim; ++ idir)
              {
                physCoordIv[idir] = m_origin[idir] + (iv[idir] + 0.5)* m_amrDx[baseLevel];
              }
            
            IndexTM<Real,SpaceDim> loEnd = physCoordIv;
            loEnd[0] -= 0.5* m_amrDx[baseLevel];
            
            IndexTM<Real,SpaceDim> hiEnd = physCoordIv;
            hiEnd[0] += 0.5*m_amrDx[baseLevel];
            
            //check endpoints
            Real loValue = implicitFunction.value(IndexTM<int,SpaceDim>::Zero,loEnd);
            Real hiValue = implicitFunction.value(IndexTM<int,SpaceDim>::Zero,hiEnd);
            
            if (loValue * hiValue <= 0)
              {
                foundGroundingLine = true;
                RootFinder rootFinder(&implicitFunction);
                // Brent returns a value between [-1,1]
                Real lengthFraction = rootFinder.Brent(loEnd,hiEnd);
                
                // scale the intersection to be in (0.0,1.0)
                lengthFraction = 0.5*(1.0 + lengthFraction);
                
                IndexTM<Real,SpaceDim>  physCoordGroundingPt = loEnd*lengthFraction + (1.0 - lengthFraction)*hiEnd;  
                setGroundingLineData(iv,physCoordGroundingPt,lengthFraction);
                
                pout()<< "physCoordGroundingPt = "<<physCoordGroundingPt<<endl;
                pout()<< "iv = "<<iv<<endl;
                
                // modify BasalFriction
                areaWeightedFriction* areaWeightedFrictionPtr = dynamic_cast<areaWeightedFriction*>(m_basalFrictionPtr);
                CH_assert(areaWeightedFrictionPtr != NULL);
                
                areaWeightedFrictionPtr->setGroundingLineData(iv,physCoordGroundingPt,lengthFraction);
              }
          }
      }
  }
#else
  {
    MayDay::Abort("Confusion about CH_USE_EB and CH_SPACEDIM == 1");
  } 
#endif
}

#ifndef CH_USE_EB
#if CH_SPACEDIM==1
void EBAmrIce::setGroundingLineData(const IntVect                & a_groundingLineIv,
                                    const IndexTM<Real,SpaceDim> & a_physCoordGroundingPt,
                                    const Real                   & a_lengthFraction)

{
  m_groundingLineIv      = a_groundingLineIv;
  m_physCoordGroundingPt = a_physCoordGroundingPt;
  
  m_lengthFraction       = a_lengthFraction;
}                
#endif
#endif

void EBAmrIce::setOrigin(const IndexTM<Real,SpaceDim>& a_origin)
{
  m_origin = a_origin;
}

#include "NamespaceFooter.H"
