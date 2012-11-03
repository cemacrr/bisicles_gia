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
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

#include "NamespaceHeader.H"

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

  //use levelsigmaCS to get H and then ebisPtr
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
#ifdef  writeTestPlots
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

              for (int dir=0; dir<2; dir++)
                {     
                  thisRhs(iv,dir) *= iceDensity*gravity*thisH(iv,0);
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
              FORT_SETFLOATINGBETA(CHF_FRA1      (thisC,0),
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

}

/// set mu coefficient (phi) prior to velocity solve
void EBAmrIce::setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef, 
                                Vector<LevelData<FluxBox>* >& a_faceMuCoef)
{
  AmrIce::setMuCoefficient(a_cellMuCoef, 
                           a_faceMuCoef);
}

#include "NamespaceFooter.H"
