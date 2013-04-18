#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "PicardSolver.H"
#include "RelaxSolver.H"
#include "NoOpSolver.H"
#include "ViscousTensorOp.H"
#include "ParmParse.H"
#include "QuadCFInterp.H"
#include "CornerCopier.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"

#ifdef CH_USE_EB
#include "LSquares.H"
#include "Notation.H"
#endif

#if    CH_SPACEDIM == 1 
#ifdef CH_USE_ICE_MF

#include "LSquares.H"
#include "Notation.H"

#endif
#endif

#include "AMRIO.H"
#include "NamespaceHeader.H"

PicardSolver::PicardSolver()
{
  // default constructor leaves things in an undefined state
  setDefaultValues();
  //m_opFactoryPtr = NULL;
  m_bc = NULL;
  m_mlOpPtr = NULL;
  m_biCGStabSolverPtr = NULL;
  m_GMRESSolverPtr = NULL;
  m_MGSolverPtr = NULL;
#ifdef CH_USE_PETSC
  m_petscSolver = NULL;
#endif
  m_bottomSolverPtr = NULL;
  m_isOpDefined = false;
  m_isSolverDefined = false;
  m_isDefined = false;
  m_vtopSafety = VTOP_DEFAULT_SAFETY;
}

PicardSolver::~PicardSolver() 
{
  if (m_biCGStabSolverPtr != NULL)
    {
      delete m_biCGStabSolverPtr;
      m_biCGStabSolverPtr = NULL;
    }

  if (m_GMRESSolverPtr != NULL)
    {
      delete m_GMRESSolverPtr;
      m_GMRESSolverPtr = NULL;
    }
#ifdef CH_USE_PETSC
  if (m_petscSolver != NULL)
    {
      delete m_petscSolver;
      m_petscSolver = NULL;
    }
#endif

  if(m_mlOpPtr != NULL)
    {
      delete m_mlOpPtr;
      m_mlOpPtr = NULL;
    }

  if (m_MGSolverPtr != NULL)
    {
      delete m_MGSolverPtr;
      m_MGSolverPtr = NULL;
    }
  
  if (m_bottomSolverPtr != NULL)
    {
      delete m_bottomSolverPtr;
      m_bottomSolverPtr = NULL;
    }

  if (m_bc != NULL)
    {
      delete m_bc;
      m_bc = NULL;
    }
}

void 
PicardSolver::define(const ProblemDomain& a_coarseDomain,
                     ConstitutiveRelation* a_constRelPtr,
		     BasalFrictionRelation* a_basalFrictionRelPtr,
                     const Vector<DisjointBoxLayout>& a_vectGrids,
                     const Vector<int>& a_vectRefRatio,
                     const RealVect& a_dxCrse,
                     IceThicknessIBC* a_bc,
                     int a_numLevels)
{
  m_coarseDomain = a_coarseDomain;
  m_domains.resize(a_numLevels);
  m_domains[0] = a_coarseDomain;
  for (int lev = 1; lev < a_numLevels; ++lev)
    {
      //m_dxs[lev] = m_dxs[lev-1] / m_refRatios[lev-1];
      m_domains[lev] = m_domains[lev-1];
      m_domains[lev].refine(a_vectRefRatio[lev-1]);
      //m_grids[lev] = a_grids[lev];
      
    }


  m_constRelPtr = a_constRelPtr;
  m_basalFrictionRelPtr = a_basalFrictionRelPtr;
  m_vectGrids = a_vectGrids;
  m_vectRefRatio = a_vectRefRatio;
  m_dxCrse = a_dxCrse;
  m_bc = a_bc->new_thicknessIBC();
  m_numLevels = a_numLevels;  
  
  m_vectMu.resize(m_numLevels);
  m_vectLambda.resize(m_numLevels);
  m_vectBeta.resize(m_numLevels);
  m_vectC.resize(m_numLevels);

  for (int lev=0; lev<m_numLevels; lev++)
    {
      m_vectMu[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_vectGrids[lev], 
                                                                                1, 
                                                                                IntVect::Unit) );


      m_vectLambda[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_vectGrids[lev], 
                                                                                    1, 
                                                                                    IntVect::Zero) );

      // beta only has one component...
      m_vectBeta[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_vectGrids[lev], 
                                                                                      1, 
                                                                                      IntVect::Zero));

      // C only has one component...
      m_vectC[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_vectGrids[lev],
									       1, 
									       IntVect::Zero));

    
  
    }



  // get this from ParmParse for now)
  ParmParse picardPP("picardSolver");
  picardPP.query("writeResidualPlots", m_writeResidToPlotFile);

  picardPP.query("absoluteTolerance", m_absolute_tolerance);

  picardPP.query("max_picard_iterations",m_max_iter);


  //defineOpFactory();

  //defineLinearSolver();
}
 
/// solve for isothermal ice
/** beta scales sliding coefficient C -- acoef in terms of the ViscousTensorOp
 */
int PicardSolver::solve(Vector<LevelData<FArrayBox>* >       & a_horizontalVel,
                        Real                                 & a_initialResidualNorm, 
                        Real                                 & a_finalResidualNorm,
                        const Real                             a_convergenceMetric,
                        const Vector<LevelData<FArrayBox>*  >& a_rhs,
                        const Vector<LevelData<FArrayBox>*  >& a_beta,
                        const Vector<LevelData<FArrayBox>*  >& a_beta0, // not used
                        const Vector<LevelData<FArrayBox>*  >& a_A,
                        const Vector<LevelData<FluxBox>  *  >& a_muCoef,
                        Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                        Real                                   a_time,
                        int                                    a_lbase, 
                        int                                    a_maxLevel)
{

  int returnCode = 0;


  bool initVelToZero = false;
  
  // initial implementation -- redefine solver for every Picard
  // iteration. not ideal, but a good place to start. Eventually grab
  // coefficients from operators and reset in existing solver
  //cell face A, needed to compute mu
  Vector<LevelData<FluxBox>* > faceA(a_maxLevel+1, NULL);
  for (int lev = 0; lev < a_maxLevel + 1; ++lev)
    {
      faceA[lev] = new LevelData<FluxBox>(m_vectGrids[lev], 
                                          a_A[lev]->nComp(), 
                                          IntVect::Zero);
      CellToEdge(*a_A[lev] , *faceA[lev]);
    }
 
  // copy beta into local storage (also not terribly efficient; we'll fix that later as well)
  for (int lev=0; lev<= a_maxLevel; lev++)
    {
      LevelData<FArrayBox>& localBeta = *m_vectBeta[lev];
      LevelData<FArrayBox>& argBeta = *a_beta[lev];
      CH_assert(localBeta.nComp() == argBeta.nComp());

      DataIterator dit = localBeta.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          localBeta[dit].copy(argBeta[dit]);
        }
    }


  computeMu(a_horizontalVel, 
            faceA, 
            a_coordSys, 
            a_time);
  
  Real& residNorm    = a_finalResidualNorm;
  Real& initialResid = a_initialResidualNorm;
  int numIter        = 0;
  Real convergenceMetric = 1.0;

#ifdef CH_USE_ICE_MF
  {
    residNorm = LARGEREALVAL;  
  }
#endif

  // section of the code for AMRViscousTensorOp
#ifndef CH_USE_EB       
#ifndef CH_USE_ICE_MF 
     
  if (m_isSolverDefined)
    {
      // simplest thing to do here is to clear the solvers; come back later and implement solver re-use
      clearLinearSolvers();
    }
  
  defineLinearSolver();
  
  m_MGSolverPtr->init(a_horizontalVel, a_rhs, a_maxLevel, a_lbase);
  
  if (initVelToZero)
    {
      for (int lev=0; lev<=a_maxLevel; lev++)
        {
          DataIterator dit = a_horizontalVel[lev]->dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*a_horizontalVel[lev])[dit].setVal(0.0);
              }
        }
    }
        
  // create local storage for AMR residual to make it easy to write the AMR residual to a plotfile if needed
  Vector<LevelData<FArrayBox>* > amrResid(a_maxLevel+1, NULL);
  Vector<LevelData<FArrayBox>* > plotData(a_maxLevel+1, NULL);
  
  Vector<string> vectName(4);
  vectName[0] = "x-residual";
  vectName[1] = "y-residual";
  vectName[2] = "x-RHS";
  vectName[3] = "y-RHS";
  
  Interval residComps(0,1);
  Interval rhsComps(2,3);
  Real bogusDt = 0;
  
  int numLevels = a_maxLevel +1;
  
  if (m_writeResidToPlotFile)
    {
      for (int lev=0; lev<=a_maxLevel; lev++)
        {
          IntVect residGhost = IntVect::Zero;
          amrResid[lev] = new LevelData<FArrayBox>(a_rhs[lev]->getBoxes(),
                                                   a_rhs[lev]->nComp(),
                                                   residGhost);
          
          plotData[lev] = new LevelData<FArrayBox>(a_rhs[lev]->getBoxes(),
                                                   2*a_rhs[lev]->nComp(),
                                                   residGhost);
        }
      
      // Compute initial residual
      // Note that whichever solver we're using, we've got the AMRMultigrid solver, so can use that to compute the residual
      initialResid = m_MGSolverPtr->computeAMRResidual(amrResid,
                                                       a_horizontalVel, 
                                                       a_rhs, 
                                                       a_maxLevel, 
                                                       a_lbase);
      for (int lev=0; lev < plotData.size(); lev++)
        {
          amrResid[lev]->copyTo(amrResid[lev]->interval(),
                                *(plotData[lev]),
                                residComps);
          
          a_rhs[lev]->copyTo(a_rhs[lev]->interval(),
                             *(plotData[lev]),
                             rhsComps);
        }
      
      string filename = "picardResidual";
      char suffix[100];
      sprintf(suffix,"%06d.%dd.hdf5", numIter, SpaceDim);
      filename.append(suffix);
      
      WriteAMRHierarchyHDF5(filename, m_vectGrids, plotData, vectName, 
                            m_coarseDomain.domainBox(), m_dxCrse[0], 
                            bogusDt, a_time, 
                            m_vectRefRatio, 
                            numLevels);
    } 
    else
      {
        
        initialResid = m_MGSolverPtr->computeAMRResidual(a_horizontalVel, 
                                                         a_rhs, 
                                                         a_maxLevel, 
                                                         a_lbase);
      }
    
  residNorm = initialResid;
  convergenceMetric = (a_convergenceMetric > 0.0)?a_convergenceMetric:initialResid;
  CH_assert(residNorm < HUGE_NORM);

  // end of section relevant only for AMRVicousTensorOP
#endif
#endif  

  if (m_verbosity >= 3)
    {
      pout() << "Picard iteration 0 max(resid) = "
             << initialResid << endl;
    }
  
  bool done = false;  
  if (residNorm < m_solver_tolerance*convergenceMetric) 
    {
      // this case will only be true if the convergence metric
      // is passed in from outside -- prevents the solver
      // from doing any work if the initial residual is small enough
      done = true;
    }
  else if (residNorm < m_absolute_tolerance) 
    {
      done = true;
    }
  else 
    {
      while (!done)
        {
          // do linear solve
          if (m_solver_type == multigrid)
            {
              // may eventually want to switch to solveNoInit
              bool zeroPhi = false;
                
              //if (numIter == 0) zeroPhi = true;
              m_MGSolverPtr->solveNoInit(a_horizontalVel, 
                                         a_rhs,
                                         a_maxLevel, 
                                         a_lbase, 
                                         zeroPhi);
            }
          else 
            {
              // in case horizontalVel has more levels than we're 
              // currently solving for (common for the case where
              // we're not refined all the way to the maximum allowable
              // level), create local Vectors with only the levels we're
              // solving for so that the MultilevelLinearOp and BiCGStab 
              // don't get confused. Kind of a kluge, but a fairly harmless
              // one, I think (DFM, 7/12/10)
              Vector<LevelData<FArrayBox>*> localVel(a_maxLevel+1, NULL);
              Vector<LevelData<FArrayBox>*> localRHS(a_maxLevel+1, NULL);
              for (int lev=0; lev<= a_maxLevel; lev++)
                {
                  localVel[lev] = a_horizontalVel[lev];
                  localRHS[lev] = a_rhs[lev];
                }
                
              if (m_solver_type == BiCGStab)
                {
                  m_biCGStabSolverPtr->solve(localVel, localRHS);
                }
              else if (m_solver_type == GMRES)
                {
                  m_GMRESSolverPtr->solve(localVel, localRHS);
                }
              else if (m_solver_type == petsc)
                {
#ifdef CH_USE_PETSC
                  // at the moment, this only works for single-level
                  CH_assert(localVel.size() == 1);
                  m_petscSolver->setInitialGuessNonzero(true);
                  m_petscSolver->solve(*localVel[0], *localRHS[0]);	
                  // this signals for next solve that its new (nonlinear)
		  m_petscSolver->resetOperator(); 
#endif
                }

#ifdef CH_USE_ICE_MF
#if    CH_SPACEDIM == 1
              else if (m_solver_type == multiFluidSolver)
                {
                  // allocate matrix
                  Real** EBMatrix;
                  
                  // domain info
                  Box domainBox  = m_domains[a_maxLevel].domainBox();
                  int xDir = 0;

                  int rows = domainBox.size()[xDir] + 1;
                  int cols = rows;
                  
                  // allocArray initializes EBmatrix to 0.0
                  LSquares matrixUtility;
                  matrixUtility.allocArray(rows,cols,EBMatrix);
                  Vector<Real> vectorVel(rows);
                  Vector<Real> vectorRhs(rows);

                  bool okay = MultiFluidOp1D(EBMatrix,
                                            a_beta,
                                            a_beta0,
                                            a_A,
                                            a_muCoef,
                                            m_vectMu,
                                            a_time,
                                            a_lbase, 
                                            a_maxLevel,
                                            m_vectGrids[a_maxLevel],
                                            m_coarseDomain.domainBox(),
                                            m_dxCrse[0],
                                            m_groundingLineIv,
                                            m_lengthFraction);

                  if (!okay)
                    {
                      MayDay::Abort("Matrix building failed.");
                    }

                  residNorm         = LARGEREALVAL;
                  Real oldResidNorm = LARGEREALVAL;
                  Real tolerance = 1.0e-8;
                  int iterationNumber = -1;
                  while (residNorm > tolerance && iterationNumber<100)
                    {
                      iterationNumber += 1;
                      okay = InvertMultiFluidOp1D(a_horizontalVel,
                                                  vectorVel,
                                                  vectorRhs,
                                                  EBMatrix,
                                                  rows,
                                                  a_rhs,
                                                  a_maxLevel,
                                                  m_groundingLineIv,
                                                  m_lengthFraction,
                                                  m_dxCrse[0]);
                      if (!okay)
                        {
                          MayDay::Abort("Matrix inversion failed");
                        }

                      // recompute mu
                      computeMu(a_horizontalVel, 
                                faceA, 
                                a_coordSys,  
                                a_time );
                      
                      okay = MultiFluidOp1D(EBMatrix,
                                            a_beta,
                                            a_beta0,
                                            a_A,
                                            a_muCoef,
                                            m_vectMu,
                                            a_time,
                                            a_lbase, 
                                            a_maxLevel,
                                            m_vectGrids[a_maxLevel],
                                            m_coarseDomain.domainBox(),
                                            m_dxCrse[0],
                                            m_groundingLineIv,
                                            m_lengthFraction);
                      if (!okay)
                        {
                          MayDay::Abort("Matrix building failed.");
                        }

                      // used for stopping criterion
                      residNorm = computeResidualNorm(EBMatrix,
                                                      vectorVel,
                                                      vectorRhs,
                                                      iterationNumber);

                      Real rate = (oldResidNorm/residNorm);
                      pout()<<endl<<"old/new = "<<rate<<endl; 
                      
                      oldResidNorm = residNorm;
                    }
                  
                  pout()<<endl<<"Number of Picard iterations = "<<iterationNumber<<endl; 
                  
                  // free matrix memory
                  matrixUtility.freeArray(rows,cols,EBMatrix);
                
                }
#endif 
#endif
              else
                {
                  MayDay::Error("Invalid solver type");
                }
            }

          Real oldResidNorm = residNorm;

#ifdef CH_USE_ICE_MF
#if    CH_SPACEDIM == 1 
          
          {
            residNorm = 0.5*m_solver_tolerance*convergenceMetric;
          }
#endif
#else
          {
            // compute new residual -- take advantage of the fact that 
            // mg solver is using the RefCountedPtr's of the coefficients, so 
            // we can just call computeAMRResidual w/o recomputing anything
            m_MGSolverPtr->m_convergenceMetric = initialResid;
              
              
            if (m_writeResidToPlotFile)
              {
                  
                residNorm = m_MGSolverPtr->computeAMRResidual(amrResid,
                                                              a_horizontalVel, 
                                                              a_rhs, 
                                                              a_maxLevel, 
                                                              a_lbase);
                CH_assert(residNorm < 1e+10);
                  
                string filename = "picardResidual";
                char suffix[100];
                sprintf(suffix,"%06d.%dd.hdf5", numIter+1, SpaceDim);
                filename.append(suffix);
                  
                for (int lev=0; lev < plotData.size(); lev++)
                  {
                    amrResid[lev]->copyTo(amrResid[lev]->interval(),
                                          *(plotData[lev]),
                                          residComps);
                      
                    a_rhs[lev]->copyTo(a_rhs[lev]->interval(),
                                       *(plotData[lev]),
                                       rhsComps);
                  }
                  
                  
                WriteAMRHierarchyHDF5(filename, m_vectGrids, plotData, vectName, 
                                      m_coarseDomain.domainBox(), m_dxCrse[0], 
                                      bogusDt, a_time, 
                                      m_vectRefRatio, 
                                      numLevels);
              }
            else
              {          
                residNorm = m_MGSolverPtr->computeAMRResidual(a_horizontalVel, 
                                                              a_rhs, 
                                                              a_maxLevel, 
                                                              a_lbase);
              }
            
          }
#endif 
          Real rate = oldResidNorm/residNorm;
            
          ++numIter;
          if (m_verbosity >= 3)
            {
              pout() << "Picard iteration " << numIter 
                     << " max(resid) = " << residNorm
                     << " ------- Rate = " << rate <<  endl;
            }
            
          if (residNorm < m_solver_tolerance*convergenceMetric) 
            {
              done = true;
              returnCode = 0;
            }
          else if (residNorm < m_absolute_tolerance) 
            {
              done = true;
              returnCode = 0;
            }
          else if (numIter > m_max_iter) 
            {
              done = true;
              if (m_verbosity >= 3)
                {
                  pout() << "Picard Solver reached max number of iterations"
                         << endl;
                }
              returnCode = 1;
            }
          
          if (!done)
            {
              // re-initialize solver with new coefficients
              // for the moment, delete existing solver and re-initialize
              // -- we'll eventually do this more intelligently
              clearLinearSolvers();
              
              defineLinearSolver();
              m_MGSolverPtr->init(a_horizontalVel, a_rhs, a_maxLevel, a_lbase);
            }
        }
      
      if (m_verbosity >= 1)
        {
          if (residNorm < m_solver_tolerance*initialResid || residNorm < m_absolute_tolerance)
            {
              pout() << "PicardSolver converged -- final norm(resid) = "
                     << residNorm << " after " << numIter << " iterations"
                     << endl;
            }
          else 
            {
              pout() << "PicardSolver NOT CONVERGED -- final norm(resid) = "
                     << residNorm << " after " << numIter << " iterations"
                     << endl;          
            }
        } // end if verbosity >= 1

      
    } // end if initial residual not zero

#ifndef CH_USE_ICE_MF
  // clean up temporary storage
  for (int lev=0; lev<amrResid.size(); lev++)
    {

      if (amrResid[lev] != NULL)
        {
          delete amrResid[lev];
          amrResid[lev] = NULL;
        }
    }
#endif
  for (int lev=0; lev<faceA.size(); lev++)
    {
      if (faceA[lev] != NULL)
	{
	  delete faceA[lev];
	  faceA[lev] = NULL;
	}
    }

  return returnCode;
}



void PicardSolver::setDefaultValues()
{
  m_solver_tolerance = 1.0e-10;
  m_absolute_tolerance = 5e-10;
  m_max_iter = 200;

#ifdef CH_USE_ICE_MF
  m_solver_type = multiFluidSolver;
#else
  m_solver_type = multigrid;
#endif
  m_bottom_solver_type = BiCGStab;
  //m_solver_type = BiCGStab;
  //m_solver_type = GMRES;
  // use isothermal ice temp from Pattyn(2003)
  m_constThetaVal = 238.15;
  m_linearsolver_tolerance = 1e-5;

  m_writeResidToPlotFile = false;
}
 
void
PicardSolver::defineLinearSolver()
{
  // for now, make this definable from ParmParse
  ParmParse picardPP("picardSolver");
   
  if (!m_isOpDefined) {
    picardPP.query("vtopSafety", m_vtopSafety);
    defineOpFactory(); 
  }
    
  if (picardPP.contains("linearSolver"))
    {
      std::string solverTypeString;
      picardPP.get("linearSolver", solverTypeString);
      if (solverTypeString == "multigrid")
        {
          m_solver_type = multigrid;
        }
#ifdef CH_USE_ICE_MF
      else if (solverTypeString == "multiFluidSolver")
        {
          m_solver_type = multiFluidSolver;
        }
#endif
      else if (solverTypeString == "BiCGStab")
        {
          m_solver_type = BiCGStab;
        }
      else if (solverTypeString == "GMRES")
        {
          m_solver_type = GMRES;
        }
      else if (solverTypeString == "petsc")
        {
#ifdef CH_USE_PETSC
          m_solver_type = petsc;
#else
          // go to MG if petsc not compiled in (this is to simplfy
          // petsc vs. MG comparison by allowing us to use the same 
          // inputs file)
          m_solver_type = multigrid;
#endif
        }
      
      else 
        {
          MayDay::Error("PicardSolver -- unknown linear solver type");
        }
    }
  
  if (picardPP.contains("mgBottomSolver"))
    {
      std::string solverTypeString;
      picardPP.get("mgBottomSolver", solverTypeString);
      if (solverTypeString == "BiCGStab")
        {
          m_bottom_solver_type = BiCGStab;
        }
      else if (solverTypeString == "relaxSolver")
        {
          m_bottom_solver_type = relaxSolver;
        }
      else if (solverTypeString == "noOpSolver")
        {
          m_bottom_solver_type = noOpSolver;
        }
      else 
        {
          MayDay::Error("PicardSolver -- unknown bottom solver type");
        }
    }

  // default has already been set
  picardPP.query("linearsolver_tolerance", m_linearsolver_tolerance);

  // default is to coarsen as much as possible
  int maxMGdepth = -1;
  picardPP.query("maxMGdepth", maxMGdepth);

  // need this for multigrid (where it's the primary linear solver)
  // or for BiCGStab (where it's the preconditioner)
  CH_assert(m_MGSolverPtr == NULL);
  m_MGSolverPtr = new AMRMultiGrid<LevelData<FArrayBox> >;
  m_MGSolverPtr->m_maxDepth = maxMGdepth;
  m_MGSolverPtr->m_verbosity = m_verbosity - 2;

  CH_assert(m_bottomSolverPtr == NULL);

  if (m_bottom_solver_type == BiCGStab)
    {
      BiCGStabSolver<LevelData<FArrayBox> >* biCGStabPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
      biCGStabPtr->m_verbosity = (m_verbosity - 3);
      m_bottomSolverPtr = biCGStabPtr;
    }
  else if (m_bottom_solver_type == relaxSolver)
    {
      RelaxSolver<LevelData<FArrayBox> >* relaxSolverPtr = new RelaxSolver<LevelData<FArrayBox> >;
      relaxSolverPtr->m_verbosity = (m_verbosity - 3);
      int numMGSmooth = 16;
      picardPP.query("numMGSmooth", numMGSmooth);
      relaxSolverPtr->m_imax = numMGSmooth;
      m_bottomSolverPtr = relaxSolverPtr;
    }
  else if (m_bottom_solver_type == noOpSolver)
    {
      NoOpSolver<LevelData<FArrayBox> >* noOpSolverPtr = new NoOpSolver<LevelData<FArrayBox> >;
      m_bottomSolverPtr = noOpSolverPtr;
    }
      
  defineMGSolver(*m_MGSolverPtr, m_bottomSolverPtr);

  CH_assert(m_mlOpPtr == NULL);
  if (m_solver_type == multigrid)
    {
      // don't need to do anything else
    }
  else if (m_solver_type == BiCGStab)
    {

      CH_assert(m_biCGStabSolverPtr == NULL);
      m_biCGStabSolverPtr = new BiCGStabSolver<Vector<LevelData<FArrayBox>* > >;

      m_mlOpPtr = defineBiCGStabSolver(*m_biCGStabSolverPtr, *m_MGSolverPtr);
      
    }
  else if (m_solver_type == GMRES)
    {

      CH_assert(m_GMRESSolverPtr == NULL);
      m_GMRESSolverPtr = new GMRESSolver<Vector<LevelData<FArrayBox>* > >;

      m_mlOpPtr = defineGMRESSolver(*m_GMRESSolverPtr, *m_MGSolverPtr);
      
    }
#ifdef CH_USE_PETSC
  else if (m_solver_type == petsc)
    {
      Real opAlpha, opBeta;
      getOperatorScaleFactors(opAlpha, opBeta);
      
      // single-level only
      if( !m_petscSolver )
	{
	  // we could delete this everytime and redo coarse grids -- should make an interval!!!
	  m_petscSolver = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
	  LinearOp<LevelData<FArrayBox> >* op = m_opFactoryPtr->AMRnewOp(m_domains[0]);      
	  m_petscSolver->define( op, false ); // just sets dx & crdx
	}      
      m_petscSolver->setVTParams( opAlpha, opBeta, &(*m_vectC[0]), &(*m_vectMu[0]), &(*m_vectLambda[0]) );
    }
#endif
#ifdef CH_USE_ICE_MF
  else if (m_solver_type == multiFluidSolver)
    {
      // do nothing
    }
#endif
  else
    {
      MayDay::Error("invalid solver type");
    }  

  m_isSolverDefined = true;
}



void PicardSolver::clearLinearSolvers()
{
  if (m_MGSolverPtr != NULL)
    {
      delete m_MGSolverPtr;
      m_MGSolverPtr = NULL;
    }
  
  if (m_bottomSolverPtr != NULL)
    {
      delete m_bottomSolverPtr;
      m_bottomSolverPtr = NULL;
    }
  
  if (m_mlOpPtr != NULL)
    {
      delete m_mlOpPtr;
      m_mlOpPtr = NULL;
    }
  
  if (m_biCGStabSolverPtr != NULL)
    {
      delete m_biCGStabSolverPtr;
      m_biCGStabSolverPtr = NULL;
    }

  if (m_GMRESSolverPtr != NULL)
    {
      delete m_GMRESSolverPtr;
      m_GMRESSolverPtr = NULL;
    }
  
  m_isSolverDefined = false;
}



void
PicardSolver::defineMGSolver(AMRMultiGrid<LevelData<FArrayBox> >& a_mgSolver,
                             LinearSolver<LevelData<FArrayBox> >* a_bottomSolverPtr)
{
  int numMGSmooth = 16;
  ParmParse pp2("picardSolver");
  pp2.query("numMGSmooth", numMGSmooth);
  // switched to "numMGSmooth" for compatibility with JFNK, etc, 
  // keeping num_smooth around for backward compatibility
  pp2.query("num_smooth", numMGSmooth);  

  // set default to be arithmetic averaging of viscosity
  int mgAverageType = CoarseAverageFace::arithmetic;
  pp2.query("mgCoefficientAverageType", mgAverageType);
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  pp2.query("mgProlongType", mgProlongType);
  ViscousTensorOp::s_prolongType = mgProlongType;

  // for symmetry with the BiCGStab case
  int max_iter = 15;
  pp2.query("max_iterations", max_iter);

  a_mgSolver.m_iterMax = max_iter;
  a_mgSolver.m_pre = numMGSmooth;
  a_mgSolver.m_post = numMGSmooth;
  a_mgSolver.m_bottom = numMGSmooth;

  a_mgSolver.define(m_coarseDomain,
                    *m_opFactoryPtr,
                    a_bottomSolverPtr,
                    m_numLevels);
  int numMG = 1;
  pp2.query("numMGCycle", numMG);
  a_mgSolver.setMGCycle(numMG);

  pp2.query("linearsolver_tolerance", m_linearsolver_tolerance);
  a_mgSolver.m_eps = m_linearsolver_tolerance;
}
 
/// define BiCGStab solver
/** returns pointer to the MultilevelLinearOp used by the solver. 
    As things stand now, we'll need to manage that one ourselves.
*/
MultilevelIceVelOp*  
PicardSolver::defineBiCGStabSolver(BiCGStabSolver<Vector<LevelData<FArrayBox>* > >& a_solver,
                                   AMRMultiGrid<LevelData<FArrayBox> >& a_mgSolver)
{
  int lBase = 0;
  CH_assert(m_opFactoryPtr != NULL);

  MultilevelIceVelOp* mlOpPtr = new MultilevelIceVelOp();
  MultilevelIceVelOp& mlOp = *mlOpPtr;

  ParmParse pp2("picardSolver");

  int numMGIter = 1;
  // num_mg is deprecated, numMGIter is better for agreement with JFNKSolver
  pp2.query("num_mg", numMGIter);
  pp2.query("numMGIter", numMGIter);

  mlOp.m_num_mg_iterations = numMGIter;
  int numMGSmooth = 4;
  pp2.query("num_smooth", numMGSmooth);
  mlOp.m_num_mg_smooth = numMGSmooth;
  int preCondSolverDepth = -1;
  pp2.query("preCondSolverDepth", preCondSolverDepth);
  mlOp.m_preCondSolverDepth = preCondSolverDepth;
  
  Real tolerance = 1.0e-7;
  pp2.query("tolerance", tolerance);
  
  int max_iter = 10;
  pp2.query("max_iterations", max_iter);
  
  Vector<RealVect> vectDx(m_numLevels);
  Vector<ProblemDomain> vectDomain(m_numLevels);
  vectDx[0] = m_dxCrse;
  vectDomain[0] = m_coarseDomain;

  // do this in case m_vectGrids has more (undefined) levels
  // than m_numLevels, which confuses the MultilevelLinearOp
  Vector<DisjointBoxLayout> localVectGrids(m_numLevels);
  localVectGrids[0] = m_vectGrids[0];

  for (int lev=1; lev<m_numLevels; lev++)
    {
      int nRef = m_vectRefRatio[lev-1];
      vectDx[lev] = vectDx[lev-1]/nRef;
      vectDomain[lev] = vectDomain[lev-1];
      vectDomain[lev].refine(nRef);
      localVectGrids[lev] = m_vectGrids[lev];
    }


  
  mlOp.define(localVectGrids, m_vectRefRatio, vectDomain,
              vectDx, m_opFactoryPtr, lBase);
  bool homogeneousBC = false;
  a_solver.define(mlOpPtr, homogeneousBC);
  a_solver.m_verbosity = m_verbosity-1;
  a_solver.m_normType = 0;
  a_solver.m_eps = tolerance;
  a_solver.m_imax = max_iter;
  
  return mlOpPtr;
}


/// define GMRES solver
/** returns pointer to the MultilevelLinearOp used by the solver. 
    As things stand now, we'll need to manage that one ourselves.
*/
MultilevelIceVelOp*  
PicardSolver::defineGMRESSolver(GMRESSolver<Vector<LevelData<FArrayBox>* > >& a_solver,
                                AMRMultiGrid<LevelData<FArrayBox> >& a_mgSolver)
{
  int lBase = 0;
  CH_assert(m_opFactoryPtr != NULL);

  MultilevelIceVelOp* mlOpPtr = new MultilevelIceVelOp();
  MultilevelIceVelOp& mlOp = *mlOpPtr;

  ParmParse pp2("picardSolver");

  int numMGIter = 1;
  pp2.query("num_mg", numMGIter);

  mlOp.m_num_mg_iterations = numMGIter;
  int numMGSmooth = 4;
  pp2.query("num_smooth", numMGSmooth);
  mlOp.m_num_mg_smooth = numMGSmooth;
  int preCondSolverDepth = -1;
  pp2.query("preCondSolverDepth", preCondSolverDepth);
  mlOp.m_preCondSolverDepth = preCondSolverDepth;
  
  Real tolerance = 1.0e-7;
  pp2.query("tolerance", tolerance);
  
  int max_iter = 10;
  pp2.query("max_iterations", max_iter);
  
  Vector<RealVect> vectDx(m_numLevels);
  Vector<ProblemDomain> vectDomain(m_numLevels);
  vectDx[0] = m_dxCrse;
  vectDomain[0] = m_coarseDomain;

  // do this in case m_vectGrids has more (undefined) levels
  // than m_numLevels, which confuses the MultilevelLinearOp
  Vector<DisjointBoxLayout> localVectGrids(m_numLevels);
  localVectGrids[0] = m_vectGrids[0];

  for (int lev=1; lev<m_numLevels; lev++)
    {
      int nRef = m_vectRefRatio[lev-1];
      vectDx[lev] = vectDx[lev-1]/nRef;
      vectDomain[lev] = vectDomain[lev-1];
      vectDomain[lev].refine(nRef);
      localVectGrids[lev] = m_vectGrids[lev];
    }


  
  mlOp.define(localVectGrids, m_vectRefRatio, vectDomain,
              vectDx, m_opFactoryPtr, lBase);
  bool homogeneousBC = false;
  a_solver.define(mlOpPtr, homogeneousBC);
  a_solver.m_verbosity = m_verbosity-1;
  //  a_solver.m_normType = 0;
  a_solver.m_eps = tolerance;
  a_solver.m_imax = max_iter;
  // maxnorm for consistency
  a_solver.m_normType = 0;
  
  return mlOpPtr;
}



void PicardSolver::defineOpFactory()
{
  if (SpaceDim <= 2)
    {
      CH_assert(m_opFactoryPtr == NULL);
     
      Real alpha, beta;
      getOperatorScaleFactors(alpha, beta);

      // for the moment, at least, this only works for dx = dy:
      if (SpaceDim > 1)
        CH_assert(m_dxCrse[0] == m_dxCrse[1]);
      
     

      BCHolder velSolveBC = m_bc->velocitySolveBC();
      m_opFactoryPtr = 
        RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
	(new ViscousTensorOpFactory(m_vectGrids, m_vectMu, m_vectLambda, m_vectC, alpha, 
				    beta, m_vectRefRatio, m_coarseDomain, m_dxCrse[0], velSolveBC, m_vtopSafety) );
      
    }
  else 
    {
      MayDay::Error("PicardSolver::defineOpFactory not implemented for dim = SpaceDim");
    }
  m_isOpDefined = true;
}

void
PicardSolver::getOperatorScaleFactors(Real& a_alpha, Real& a_beta) const
{
  // in 2d, we can use the existing Chombo ViscousTensorOp
  a_alpha = -1.0;
  //alpha = 0.0;
  // beta = 0.5 for similarity with algorithm doc...
  //Real beta = 0.5;
  // beta = 1.0 for smiliarity with Pattyn etc
  a_beta = 1.0; 
}
 
// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
void 
PicardSolver::computeMu(Vector<LevelData<FArrayBox>*        >& a_horizontalVel,
			Vector<LevelData<FluxBox>*          >& a_A, 
                        Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                        Real                                   a_time)
{
  Vector<LevelData<FArrayBox>* >& avgDownVel = a_horizontalVel;
    
  ProblemDomain levelDomain = m_coarseDomain;
  for (int lev=0; lev<m_numLevels; lev++)
    {
      const DisjointBoxLayout& levelGrids   = m_vectGrids  [lev];
      const LevelSigmaCS     & levelCS      = *a_coordSys  [lev];

      LevelData      <FArrayBox>& levelVel  = *avgDownVel  [lev];
      LevelData      <FArrayBox>& levelC    = *m_vectC     [lev];
      const LevelData<FArrayBox>& levelBeta = *m_vectBeta  [lev];

      LevelData<FluxBox>& levelMu           = *m_vectMu    [lev];
      LevelData<FluxBox>& levelA            = *a_A         [lev];
      LevelData<FluxBox>& levelLambda       = *m_vectLambda[lev];
      
      // if there is a finer level, average-down the current velocity field
      if (lev < (m_numLevels -1) )
        {
          LevelData<FArrayBox>& finerVel = *avgDownVel[lev+1];

          CoarseAverage averager(finerVel.getBoxes(),
                                 levelGrids,
                                 finerVel.nComp(),
                                 m_vectRefRatio[lev]);

          averager.averageToCoarse(levelVel, finerVel);
        }

      // just in case, add an exchange here
       levelVel.exchange();

      // first set BC's on vel
      m_bc->velocityGhostBC(levelVel,
                            levelCS,
                            levelDomain, a_time);

      // this is a bit of a strange way to do this, but we're avoiding an 
      // assertion in the case where we don't have any grids on this level
      // (in which case we don't care what dx is anyway)
      Real dxLevel = levelCS.dx()[0];

      if (lev > 0) 
        {
          QuadCFInterp qcfi(levelGrids, 
                            &m_vectGrids[lev-1],
                            dxLevel, 
                            m_vectRefRatio[lev-1], 
                            2, 
                            levelDomain);

          qcfi.coarseFineInterp(levelVel, *avgDownVel[lev-1]);
        }

      //slc : qcfi.coarseFineInterp fills the edges of lev > 0 cells
      //but not the corners. We need them filled to compute the
      //rate-of-strain invariant, so here is a bodge for now
      DataIterator dit = levelGrids.dataIterator();
      if (SpaceDim == 2)
      	{
      	  for (dit.begin(); dit.ok(); ++dit)
      	    {
      	      Box sbox = levelVel[dit].box();
      	      sbox.grow(-1);
      	      FORT_EXTRAPCORNER2D(CHF_FRA(levelVel[dit]),
      				  CHF_BOX(sbox));
      	    }
	  
      	}

      // actually need to use a cornerCopier, too...
      CornerCopier cornerCopier(levelGrids, levelGrids, 
                                levelDomain,levelVel.ghostVect(),
                                true);
      levelVel.exchange(cornerCopier);      

      // constant temperature...
      LevelData<FluxBox> theta(levelGrids, 1, IntVect::Unit);   
      LevelData<FArrayBox> cellTheta(levelGrids, 1, IntVect::Unit);         
      for (dit.begin(); dit.ok(); ++dit)
        {
          theta[dit].setVal(m_constThetaVal);
          cellTheta[dit].setVal(m_constThetaVal);          
        }

      LevelData<FArrayBox>* crseVelPtr = NULL;
      int nRefCrse = -1;

      if (lev > 0)
        {
          crseVelPtr = a_horizontalVel[lev-1];
          nRefCrse = m_vectRefRatio[lev-1];
        }

      IntVect muGhost = IntVect::Zero;
      (*m_constRelPtr).computeFaceMu(levelMu,
                                     levelVel,
                                     crseVelPtr,
                                     nRefCrse,
                                     levelA,
                                     levelCS,
				     m_domains[lev],
                                     muGhost);
      levelMu.exchange();
      // now multiply by ice thickness H
      const LevelData<FluxBox>& faceH = levelCS.getFaceH();
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
              levelMu[dit][dir].mult(faceH[dit][dir],
                                     levelMu[dit][dir].box(),
                                     0,
                                     0,
                                     1);
	    }
	            
          // lambda = 2*mu
          FluxBox& lambda = levelLambda[dit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              lambda[dir].copy(levelMu[dit][dir]);
              lambda[dir] *= 2.0;
            }
	  
	  // also update alpha (or C)
          const Box& gridBox = levelGrids[dit];
	  m_basalFrictionRelPtr->computeAlpha(levelC                             [dit], 
                                              levelVel                           [dit], 
                                              levelCS.getThicknessOverFlotation()[dit], 
                                              levelBeta                          [dit],
                                              levelCS.getFloatingMask()          [dit],
                                              gridBox);

        }

      if (lev != m_numLevels-1)
        {
          levelDomain.refine(m_vectRefRatio[lev]);
        }
    }
}
#ifdef CH_USE_ICE_MF 
#if    CH_SPACEDIM == 1 
Real PicardSolver::computeResidualNorm(Real**                                 a_EBMatrix,
                                       const Vector<Real>                   & a_vectorVel,
                                       const Vector<Real>                   & a_vectorRhs,
                                       const int                            & a_iterationNumber)
{ 
  int xDir = 0;
  int rows = a_vectorVel.size();
  
  // Apply the operator 
  LSquares matrixUtility;
  Vector<Real>  Ax(rows);
  matrixUtility.AtimesX(a_EBMatrix,
                        a_vectorVel,
                        rows,
                        Ax);
  
  // calculate residual
  Vector<Real> residual(rows);
  setResidual1D(residual,Ax,a_vectorRhs);
  
  Real maxNorm;
  int maxAchieved;
  Real L1Norm;
  Real L2Norm;
  
  // calculate norm of residual
  residualNorm1D(maxNorm,maxAchieved,L1Norm,L2Norm,residual,m_dxCrse[0]);
   
  
  
  pout()<<"maxNorm = "<<maxNorm<<endl;
  pout()<<"L1Norm  = "<<L1Norm<<endl;
  pout()<<"L2Norm  = "<<L2Norm<<endl;
  pout()<<endl;
  
  pout()<<"Iteration number "<<a_iterationNumber<<endl;
  pout()<<"Maximum value achieved at cell "<<maxAchieved<<endl;
  pout()<<endl;
  // used for stopping criterion
  return(maxNorm);
}

int PicardSolver::InvertMultiFluidOp1D(Vector<LevelData<FArrayBox>* >       & a_horizontalVel,
                                       Vector<Real>                         & a_vectorVel,
                                       Vector<Real>                         & a_vectorRhs,
                                       Real **                                a_EBMatrix,
                                       const int                            & a_EBMAtrixRows,
                                       const Vector<LevelData<FArrayBox>* > & a_rhs,
                                       const int                            & a_maxLevel,
                                       const IntVect                        & a_groundingLineIv,
                                       const Real                           & a_lengthFraction,
                                       const Real                           & a_dx)
{
  // allocate matrix
  Real** EBMatrixCopy;
  
  // domain info
  DisjointBoxLayout dbl = m_vectGrids[a_maxLevel];
  Box domainBox  = m_domains[a_maxLevel].domainBox();
  int xDir = 0;
  
  int rows = domainBox.size()[xDir] + 1;
  int cols = rows;
  
  // allocArray initializes EBmatrix to 0.0
  LSquares matrixUtility;
  matrixUtility.allocArray(rows,cols,EBMatrixCopy);
  for (int irow = 0; irow < rows; ++irow)
    {
      for (int icol = 0; icol < cols; ++icol)
        {
          EBMatrixCopy[irow][icol] = a_EBMatrix[irow][icol];
        }
    }
   
  // cell-centered rhs
  const LevelData<FArrayBox>* finestLevelRhs  = a_rhs   [a_maxLevel];
  const LevelData<FArrayBox>& rhsLD           = (*finestLevelRhs);
  setVector(a_vectorRhs,
            rhsLD,
            a_groundingLineIv,
            a_lengthFraction,
            a_dx);

  //debug statement
#if 1
  outputVector(a_rhs,
               a_groundingLineIv,
               a_lengthFraction,
               a_dx);
#endif

  // multiply rhs by -1.0, because of diference with main code over sign convention
  Vector<Real> rhsCopy(a_vectorRhs.size());
  for (int irow = 0; irow < a_vectorRhs.size(); ++ irow)
    {
      a_vectorRhs[irow] *= -1.0;
      rhsCopy[irow]      = a_vectorRhs[irow];
    }
  

  // gaussian elimination with partial pivoting
  matrixUtility.gaussElim(EBMatrixCopy,rhsCopy);
     
  // solve for unknown velocity and copy unknown onto vel
  matrixUtility.backSolve(EBMatrixCopy,rhsCopy,a_EBMAtrixRows,a_vectorVel);
 
 // free matrix memory
  matrixUtility.freeArray(rows,cols,EBMatrixCopy);

 // unknown velocity for the matrix equation
  LevelData<FArrayBox>* finestLevelVel = a_horizontalVel[a_maxLevel];

  //assign unknown to velocity in a dataiterator, boxiterator
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      // cell-centered data
      FArrayBox& finestLevelVelFab  = (*finestLevelVel)[dit()];
      
      // box iterate
      Box dblBox = dbl[dit()];
      for (BoxIterator bit(dblBox); bit.ok(); ++bit)
        {
          // IntVect
          const IntVect& iv = bit();
          
          // assign finestLevelRhsFab
          if (iv[0] < a_groundingLineIv[0])
            {
              finestLevelVelFab(iv,0) = a_vectorVel[iv[0]];
            }
          else if (iv[0] == a_groundingLineIv[0])
            {
              finestLevelVelFab(iv) = a_lengthFraction*a_vectorVel[iv[0]] + (1.0 - a_lengthFraction)*a_vectorVel[iv[0] + 1];
            }
          else if (iv[0] > a_groundingLineIv[0])
            {
              finestLevelVelFab(iv,0) = a_vectorVel[iv[0] + 1];
            }
        }
    }
  
#if 1
  pout()<<"velocity vector"<<endl;
  outputVector(a_horizontalVel,
               a_groundingLineIv,
               a_lengthFraction,
               a_dx);
#endif
  return true;
}

int PicardSolver::checkWork(Vector<LevelData<FArrayBox>* >       & a_horizontalVel,
                            Real **                                a_EBMatrix,
                            const int                            & a_EBMAtrixRows,
                            const Vector<LevelData<FArrayBox>* > & a_rhs,
                            const int                            & a_maxLevel,
                            const IntVect                        & a_groundingLineIv,
                            const Real                           & a_lengthFraction)
{
  
  LSquares matrixUtility;

  LevelData      <FArrayBox>* finestLevelVel  = a_horizontalVel[a_maxLevel];
  const LevelData<FArrayBox>* finestLevelRhs  = a_rhs          [a_maxLevel];

  // domain info
  DisjointBoxLayout dbl = m_vectGrids[a_maxLevel];
  Box domainBox         = m_domains  [a_maxLevel].domainBox();

  int xDir = 0;
  int rows = domainBox.size()[xDir] + 1;

  // rhs for the matrix equation
  Vector<Real> rhs;
  rhs.resize(rows);

 // unknown velocity for the matrix equation
  Vector<Real> unknown;
  unknown.resize(rows);
  
  // create a vector rhs
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      // cell-centered data
      const FArrayBox& rhsFab = (*finestLevelRhs)[dit()];
                         
      // box iterate
      Box dblBox = dbl[dit()];
      
      for (BoxIterator bit(dblBox); bit.ok(); ++bit)
        {
          // IntVect
          const IntVect& iv = bit();
          const int irow    = iv[0];

          // assign rhs
          if (iv[0] < a_groundingLineIv[0])
            {
              rhs[irow]     = rhsFab(iv,0);
            }
          else if (iv[0] == a_groundingLineIv[0])
            {
              Real loCentroidFrac = 0.5*a_lengthFraction;
              Real hiCentroidFrac = a_lengthFraction + 0.5*(1.0 - a_lengthFraction); 
              rhs[irow]     = (1.0 - loCentroidFrac)*rhsFab(iv,0) + loCentroidFrac*rhsFab(iv + BASISV(xDir),0);
              rhs[irow + 1] = (1.0 - hiCentroidFrac)*rhsFab(iv,0) + hiCentroidFrac*rhsFab(iv + BASISV(xDir),0);
            }
          else
            {
              rhs[irow + 1] = rhsFab(iv);
            }
        }
    }
  
  // save rhs
  Vector<Real> originalRhs(rhs.size());
  for (int irow = 0; irow < rhs.size(); ++ irow)
    {
      originalRhs[irow] = rhs[irow];
    }
  
  // gaussian elimination with partial pivoting
   matrixUtility.gaussElim(a_EBMatrix,rhs);
  
  // solve for unknown velocity and copy unknown onto vel
   matrixUtility.backSolve(a_EBMatrix,rhs,a_EBMAtrixRows,unknown);
   for (int irow = 0; irow < rhs.size(); ++ irow)
     {
       if (irow < 10)
         {
           pout()<<"velocity["<<irow<<"] = "<<unknown[irow]<<endl;
         }
       

       if (irow > 45 && irow < 55)
         {
           pout()<<"velocity["<<irow<<"] = "<<unknown[irow]<<endl;
         }
       if (irow > 118)
         {
           pout()<<"velocity["<<irow<<"] = "<<unknown[irow]<<endl;
         }
     }
   // recover rhs
   Vector<Real>  Ax(rows);
   matrixUtility.AtimesX(a_EBMatrix,
                         unknown,
                         rows,
                         Ax);

   //compare A*velocity to original rhs
  Vector<Real> residual(rhs.size());
  for (int irow = 0; irow < rhs.size(); ++ irow)
    {
      residual[irow] = originalRhs[irow] - Ax[irow];
    }

#if 1
  for (int irow = 0; irow < rows; ++irow)
    {
      pout()<<"Ax      ["<<irow<<"]="<<Ax      [irow]<<endl;
      pout()<<"rhs     ["<<irow<<"]="<<rhs     [irow]<<endl;
      pout()<<"residual["<<irow<<"]="<<residual[irow]<<endl;
      pout()<<endl;
    }
#endif     
  
#if 0  
  matrixUtility.output(rhs.size(),
                       rhs.size(),
                       a_EBMatrix,
                       "matrix");
#endif
  return true;
}

int PicardSolver::MultiFluidOp1D(Real**                                             & a_EBMatrix,   
                                 const Vector<LevelData              <FArrayBox>* > & a_beta,
                                 const Vector<LevelData              <FArrayBox>* > & a_beta0,
                                 const Vector<LevelData              <FArrayBox>* > & a_A,
                                 const Vector<LevelData              <FluxBox>* >   & a_muCoef,
                                 const Vector<RefCountedPtr<LevelData<FluxBox> > >  & a_mu,
                                 const Real                                         & a_time,
                                 const int                                          & a_lbase, 
                                 const int                                          & a_maxLevel,
                                 const DisjointBoxLayout                            & a_dbl,
                                 const Box                                          & a_domainBox,
                                 const Real                                         & a_dx,
                                 const IntVect                                      & a_groundingLineIv,
                                 const Real                                         & a_lengthFraction)
{
  const LevelData<FArrayBox>* finestLevelBeta = a_beta[a_maxLevel];
  
  // mu and muCoef are face-centered
  const LevelData<FluxBox>* finestLevelMuCoef = a_muCoef[a_maxLevel]; 
  const LevelData<FluxBox>* finestLevelMu     = a_mu    [a_maxLevel];
    
  // dx
  Real dx2 = a_dx*a_dx;
  
  // indexing macros
  int xDir       = 0;
  int zerothComp = 0;
    
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      // cell-centered data
      const FArrayBox& finestLevelBetaFab = (*finestLevelBeta )[dit()];
      
      // face-centered data
      const FluxBox& finestLevelMuCoefFab = (*finestLevelMuCoef)[dit()];
      const FluxBox& finestLevelMuFab     = (*finestLevelMu)    [dit()];
      
      // box iterate
      Box dblBox = a_dbl[dit()];
      for (BoxIterator bit(dblBox); bit.ok(); ++bit)
        {
          // IntVect
          const IntVect& iv = bit();
          
          // matrix elements
          int irow = iv[0];

          // typical interior cells
          Real diagonalElement = LARGEREALVAL;
          Real hiSideElement   = LARGEREALVAL;
          Real loSideElement   = LARGEREALVAL;
          
          // groundingLineIv
          Vector<Real> groundingPtDiagonal(2);
          Vector<Real> singleOffsetHi     (2);
          Vector<Real> singleOffsetLo     (2);
          Vector<Real> doubleOffsetHi     (2);
          Vector<Real> doubleOffsetLo     (2);
          
          // friction
          Real beta   = finestLevelBetaFab(iv,zerothComp);
          Real betaHi = beta;
          Real betaLo = beta;
          
          if (iv[0] != a_domainBox.bigEnd()[0])
            {
              betaHi = finestLevelBetaFab (iv + BASISV(xDir),zerothComp);
            }
          if (iv[0] != a_domainBox.smallEnd()[0])
            {
              betaLo = finestLevelBetaFab (iv - BASISV(xDir),zerothComp);
            }
          
          // coefficient of viscosity
          Real muCoef   = finestLevelMuCoefFab[xDir](iv,zerothComp);
          Real muCoefLo = muCoef;
          Real muCoefHi = muCoef;

          // debugging concern about box boundaries. Remember that I added a ghostCell to mu
          muCoefHi = finestLevelMuCoefFab[xDir](iv + BASISV(xDir),zerothComp);
          
          // viscosity
          Real mu   = finestLevelMuFab[xDir](iv,zerothComp);
          Real muLo = mu;
          Real muHi = mu;
          
          // debugging concern about box boundaries. Remember that I added a ghostCell to mu
          muHi = finestLevelMuFab[xDir](iv + BASISV(xDir),zerothComp);
                   
          // effective viscosity
          Real totalMuLo = muCoefLo*muLo;
          Real totalMuHi = muCoefHi*muHi;
          
          // common factor in matrix elements                          
          Real preFactor = 4.0/dx2;
          
          if ((iv[0] != a_groundingLineIv[0]) && (iv[0] != a_domainBox.smallEnd()[0])&& (iv[0] != a_domainBox.bigEnd()[0]))
            {
              // neither the hi nor the lo are in the groundingLineIv
              if( ((iv[0] + 1) != a_groundingLineIv[0]) &&  ( (iv[0] - 1) != a_groundingLineIv[0]))
                {
                  // diagonal
                  diagonalElement = beta + preFactor*(totalMuHi+totalMuLo);
                  
                  // hi side
                  hiSideElement   = -preFactor*totalMuHi;

                  // lo side
                  loSideElement   = -preFactor*totalMuLo;
                }
              else if ((iv[0] + 1) == a_groundingLineIv[0])
                {
                  Real totalMuExtrapHi = LARGEREALVAL;

                  Real extrapP    = a_dx*(iv[0] + 1);
                  
                  Real pLo        = LARGEREALVAL;
                  Real valPLo     = LARGEREALVAL;
                  Real valPHi     = LARGEREALVAL;
                  Real distPLoPHi  = LARGEREALVAL;
              
                  pLo        = a_dx*(iv[0] - 1);
                  valPLo     = finestLevelMuFab[xDir](iv - BASISV(xDir),zerothComp);
                  valPHi     = finestLevelMuFab[xDir](iv               ,zerothComp);
                  distPLoPHi = a_dx;

                  // extrapolate mu to a (lo) face of the grounding Line iv
                  totalMuExtrapHi = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                                    
                  Real distFracHiSide = 0.5*(1.0 + a_lengthFraction);
                  Real distFrac2D     = 0.5*distFracHiSide + 0.5;
                    
                  // diagonal
                  diagonalElement = beta + preFactor*(totalMuExtrapHi/(distFracHiSide*distFrac2D) + totalMuLo/distFrac2D);
                  
                  // hi side in groundingLineIv
                  hiSideElement   = -preFactor*totalMuExtrapHi/(distFracHiSide*distFrac2D);
                                    
                  // lo side
                  loSideElement   = -preFactor*totalMuLo/distFrac2D;
                }
              else
                {
                  Real totalMuExtrapLo = LARGEREALVAL;

                  Real extrapP    = a_dx*(iv[0]);
                  
                  Real pLo        = LARGEREALVAL;
                  Real valPLo     = LARGEREALVAL;
                  Real valPHi     = LARGEREALVAL;
                  Real distPLoPHi  = LARGEREALVAL;
              
                  pLo        = a_dx*(iv[0] + 1);
                  valPLo     = finestLevelMuFab[xDir](iv + BASISV(xDir)               ,zerothComp);
                  valPHi     = finestLevelMuFab[xDir](iv + BASISV(xDir) + BASISV(xDir),zerothComp);
                  distPLoPHi = a_dx;

                  // extrapolate mu to a  (hi) face of the grounding Line iv
                  totalMuExtrapLo = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                      
                  Real distFracLoSide = 0.5*(2.0 - a_lengthFraction);
                  Real distFrac2D     = 0.5*distFracLoSide + 0.5;
          
                  // diagonal
                  diagonalElement = beta + preFactor*(totalMuHi/distFrac2D + totalMuExtrapLo/(distFrac2D*distFracLoSide));
                  
                  // lo side in groundingLineIv
                  loSideElement   = -preFactor*totalMuExtrapLo/(distFrac2D*distFracLoSide);
                  
                  // lo side
                  hiSideElement   = -preFactor*totalMuHi/distFrac2D;
                }
            }
          else if(iv[0] == a_domainBox.smallEnd()[0])
            {
              // u = 0
              CH_assert(iv[0] != a_groundingLineIv[0]);
              
              // diagonal
              diagonalElement = beta + preFactor*(totalMuHi + 2.0*totalMuLo);
              
              // hi side
              hiSideElement      = -preFactor*totalMuHi;
                          
              // lo side
              loSideElement      = 0.0;
            }
          else if (iv[0] == a_domainBox.bigEnd()[0])
            {
              // du/dx = 0
              CH_assert(iv[0] != a_groundingLineIv[0]);
              
              Real muBar  = 0.5*(totalMuHi + totalMuLo);
              Real muDiff =     (totalMuHi - totalMuLo);
              Real dMuDx  = muDiff/a_dx;
              
              // diagonal
              diagonalElement = beta - (4.0/a_dx)*(muBar - 0.5*dMuDx);
              
              // lo side
              loSideElement      =   + (4.0/a_dx)*(muBar - 0.5*dMuDx); 
                          
              // lo side
              hiSideElement      = 0.0;
            }
          else
            {
              // viscosity at the groundling line
              Real totalMuBdLo = LARGEREALVAL;
              Real totalMuBdHi = LARGEREALVAL;

              Real extrapP    = a_dx*(iv[0] + a_lengthFraction);
              
              Real pLo        = LARGEREALVAL;
              Real valPLo     = LARGEREALVAL;
              Real valPHi     = LARGEREALVAL;
              Real distPLoPHi  = LARGEREALVAL;
              
              // face centered data
              pLo        = a_dx*(iv[0] - 2);
              valPLo     = finestLevelMuFab[xDir](iv  - BASISV(xDir) - BASISV(xDir),zerothComp);
              valPHi     = finestLevelMuFab[xDir](iv  - BASISV(xDir)                  ,zerothComp);
              distPLoPHi = a_dx;
              totalMuBdLo = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
              
              pLo        = a_dx*(iv[0] + 2);
              valPLo     = finestLevelMuFab[xDir](iv + BASISV(xDir) + BASISV(xDir)               ,zerothComp);
              valPHi     = finestLevelMuFab[xDir](iv + BASISV(xDir) + BASISV(xDir) + BASISV(xDir),zerothComp);
              distPLoPHi = a_dx;
              totalMuBdHi = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                          
              // matrix elements associated with the the grounding line.
              Real loDist = a_dx*       a_lengthFraction;
              Real hiDist = a_dx*(1.0 - a_lengthFraction);

              Real hiD1 = 0.5*a_dx + hiDist;
              Real hiD2 = a_dx + hiD1;
              
              Real loD1 = 0.5*a_dx + loDist;
              Real loD2 = a_dx + loD1;
              
              // R.K. Crockett 2010 JCP:/partial /partial u|_Bd = (1/(d2-d1))*(d2/d1 - d1/d2)*u_Bd + (1/(d2-d1))*((d1/d2)*u_2 - (d2/d1)*u_1)
              // this normal points toward the fluid (on both sides)
              Real hi1 = (hiD2/hiD1)/(hiD2 - hiD1); 
              Real hi2 = (hiD1/hiD2)/(hiD2 - hiD1);
              
              Real lo1 = (loD2/loD1)/(loD2 - loD1); 
              Real lo2 = (loD1/loD2)/(loD2 - loD1); 
              
              Real hiDiff = hi1 - hi2;
              Real loDiff = lo1 - lo2;
             
              Real denominator = totalMuBdHi*hiDiff + totalMuBdLo*loDiff;
              Real outerTermLo = -4.0/loDist;
              Real outerTermHi = -4.0/hiDist;
              
              Real u7 = (iv[0]+2.5)*a_dx                      ;
              Real u6 = (iv[0]+1.5)*a_dx                      ;

              Real u5 = iv[0]*a_dx + loDist + 0.5*hiDist;
              Real u4 = iv[0]*a_dx +          0.5*loDist;

              Real u3 = (iv[0] - 0.5)*a_dx                      ;
              Real u2 = (iv[0] - 1.5)*a_dx                      ;

              Real uBdExact      = sin(iv[0]*a_dx + loDist);
              Real uBdDerivExact = cos(iv[0]*a_dx + loDist);

              Real uDerivLo = loDiff*uBdExact - lo1*sin(u3) + lo2*sin(u2);
              Real uDerivHi = hiDiff*uBdExact - hi1*sin(u6) + hi2*sin(u7);

              Real uBest = (totalMuBdLo*(lo1*sin(u3) - lo2*sin(u2)) + totalMuBdHi*(hi1*sin(u6) - hi2*sin(u7)))/denominator;
              
              Real lou3Coef = -lo1 + totalMuBdLo*lo1*loDiff/denominator;
              Real lou2Coef =  lo2 - totalMuBdLo*lo2*loDiff/denominator;
              Real lou6Coef =  totalMuBdHi*hi1*loDiff      /denominator;
              Real lou7Coef = -totalMuBdHi*hi2*loDiff      /denominator;

              Real hiu3Coef = -lo1 + totalMuBdLo*lo1*loDiff/denominator;
              Real hiu2Coef =  lo2 - totalMuBdLo*lo2*loDiff/denominator;
              Real hiu6Coef =  totalMuBdHi*hi1*loDiff      /denominator;
              Real hiu7Coef = -totalMuBdHi*hi2*loDiff      /denominator;
              
              Real uDerivLoest = lou2Coef*sin(u2) + lou3Coef*sin(u3) + lou6Coef*sin(u6) + lou7Coef*sin(u7);
              

              if (denominator != 0.0)
                {
                  groundingPtDiagonal[0] = betaLo - outerTermLo*(totalMuLo*2.0/(loDist  + a_dx)                                                          ); 
                             
                  singleOffsetHi     [0] =          outerTermLo*(                                  totalMuBdLo*(      loDiff*totalMuBdHi*hi1/denominator));
                  singleOffsetLo     [0] =          outerTermLo*((totalMuLo*2.0/(loDist + a_dx)) - totalMuBdLo*(lo1 - loDiff*totalMuBdLo*lo1/denominator));
                  
                  doubleOffsetHi     [0] =          outerTermLo*(                                - totalMuBdLo*(      loDiff*totalMuBdHi*hi2/denominator));
                  doubleOffsetLo     [0] =          outerTermLo*(                                  totalMuBdLo*(lo2 - loDiff*totalMuBdLo*lo2/denominator));
                  
                  groundingPtDiagonal[1] = betaHi - outerTermHi*( totalMuHi*2.0/(hiDist + a_dx)                                                            ); 
                             
                  singleOffsetHi     [1] =          outerTermHi*( totalMuHi*2.0/(hiDist + a_dx)    -totalMuBdHi*(hi1  - hiDiff*totalMuBdHi*hi1/denominator));
                  singleOffsetLo     [1] =          outerTermHi*(                                   totalMuBdHi*(       hiDiff*totalMuBdLo*lo1/denominator));
                  
                  doubleOffsetHi     [1] =          outerTermHi*(                                   totalMuBdHi*(hi2  - hiDiff*totalMuBdHi*hi2/denominator));
                  doubleOffsetLo     [1] =          outerTermHi*(                                 - totalMuBdHi*(       hiDiff*totalMuBdLo*lo2/denominator));
                }
              else
                {
                  MayDay::Abort("Need to think more about this case");
                  // the derivatives on each side of the grounding line are equal.
                  groundingPtDiagonal[0] = beta + preFactor*(totalMuHi+totalMuLo);
                  groundingPtDiagonal[1] = beta + preFactor*(totalMuHi+totalMuLo);
                  
                  singleOffsetHi     [0] = -preFactor*totalMuHi;
                  singleOffsetLo     [0] = -preFactor*totalMuLo;
                  
                  doubleOffsetHi     [0] = 0.0;
                  doubleOffsetLo     [0] = 0.0;
                  
                  singleOffsetHi     [1] = -preFactor*totalMuHi;
                  singleOffsetLo     [1] = -preFactor*totalMuLo;;
                  
                  doubleOffsetHi     [1] = 0.0;
                  doubleOffsetLo     [1] = 0.0;
                }
            }
          
          // assign matrix elements
          if (iv[0] < a_groundingLineIv[0])
            {
              a_EBMatrix[irow][irow    ] = diagonalElement;
              a_EBMatrix[irow][irow + 1] = hiSideElement;
          
              // check for low end of domain
              if (irow > a_domainBox.smallEnd()[0])
                {
                  a_EBMatrix[irow][irow - 1] = loSideElement;
                }
            }
          else if (iv[0] ==  a_groundingLineIv[0])
            {
              a_EBMatrix[irow][irow    ] = groundingPtDiagonal[0];
              a_EBMatrix[irow][irow + 2] = singleOffsetHi     [0];
              a_EBMatrix[irow][irow - 1] = singleOffsetLo     [0];
              a_EBMatrix[irow][irow + 3] = doubleOffsetHi     [0];
              a_EBMatrix[irow][irow - 2] = doubleOffsetLo     [0];
                  
              a_EBMatrix[irow + 1][irow     + 1] = groundingPtDiagonal[1];
              a_EBMatrix[irow + 1][irow + 1 + 1] = singleOffsetHi     [1];
              a_EBMatrix[irow + 1][irow - 2 + 1] = singleOffsetLo     [1];
              a_EBMatrix[irow + 1][irow + 2 + 1] = doubleOffsetHi     [1];
              a_EBMatrix[irow + 1][irow - 3 + 1] = doubleOffsetLo     [1];
            }
          else
            {
              a_EBMatrix[irow + 1][irow + 1    ] = diagonalElement;
              
              // check for hi end of domain
              if (irow < a_domainBox.bigEnd()[0])
                {
                  a_EBMatrix[irow + 1][irow + 1 + 1] = hiSideElement;
                }
              
              a_EBMatrix[irow + 1][irow - 1 + 1] = loSideElement;
            }
        }
    }
  
  return true;
}
void PicardSolver::setResidual1D(Vector<Real>       & a_residual,
                                 const Vector<Real> & a_Ax,
                                 const Vector<Real> & a_exactX)
{
  for (int iVal = 0; iVal < a_residual.size();++iVal)
    {
      a_residual[iVal] = a_Ax[iVal] - a_exactX[iVal];
    }
}

void PicardSolver::residualNorm1D(Real               & a_maxNorm,
                                  int                & a_maxNormAchieved,
                                  Real               & a_L1Norm,
                                  Real               & a_L2Norm,
                                  const Vector<Real> & a_residual,
                                  const Real         & a_deltaX)
{
  a_maxNorm = 0.0;
  a_L1Norm  = 0.0;
  a_L2Norm  = 0.0;

  for(int iRes = 0; iRes < a_residual.size();++ iRes)
    {
      if (Abs(a_residual[iRes]) > a_maxNorm)
        {
          a_maxNorm = Abs(a_residual[iRes]);
          a_maxNormAchieved = iRes;
        }

      a_L1Norm += Abs(a_residual[iRes])*a_deltaX;
      a_L2Norm += Abs(a_residual[iRes])*Abs(a_residual[iRes])*a_deltaX;
    }

  a_L1Norm  /= (a_residual.size()*a_deltaX);
  a_L2Norm  /= (a_residual.size()*a_deltaX);
  a_L2Norm = sqrt(a_L2Norm);
}
Real PicardSolver::linearExtrapolation(const Real & a_extrapP,
                                       const Real & a_pLo,
                                       const Real & a_valPLo,
                                       const Real & a_valPHi,
                                       const Real & a_distPLoPHi)
{
  Real slope = (a_valPHi - a_valPLo)/a_distPLoPHi;

  // equation of line through a_valPLo and a_valPHi evaluated at an increment of extrapP away from PLo
  Real retval = (a_extrapP-a_pLo)*slope + a_valPLo;

  return retval;
}

void PicardSolver::setVector(Vector<Real>               & a_vector,
                             const LevelData<FArrayBox> & a_levelData,
                             const IntVect              & a_groundingLineIv,
                             const Real                 & a_lengthFraction,
                             const Real                 & a_dx)
{
  int maxLevel = 0;
  int xDir = 0;

  DisjointBoxLayout dbl = m_vectGrids[maxLevel];
  Box domainBox         = m_domains  [maxLevel].domainBox();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& fab = a_levelData[dit()];

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit)
        {
          // IntVect
          const IntVect& iv = bit();
          const int irow = iv[0];
           
          // assign rhs
          if (iv[0] < a_groundingLineIv[0])
            {
              a_vector[irow] = fab(iv,0);
            }
          else if (iv[0] == a_groundingLineIv[0])
            {
              if (a_lengthFraction < 0.5)
                {
                  // loCentroid extrapolation info
                  IntVect loVofLoCentroid = iv - BASISV(xDir) - BASISV(xDir);
                  IntVect hiVofLoCentroid = iv - BASISV(xDir)               ;
                  Real extrapP            = a_dx*(iv[0] + 0.5*a_lengthFraction); 
                  Real pLo                = a_dx*(iv[0] - 2 + 0.5);
                  Real valPLo             = fab(loVofLoCentroid,0);
                  Real valPHi             = fab(hiVofLoCentroid,0);
                  Real distPLoPHi         = a_dx;
                           
                  // hiCentroid interp info
                  Real loSideCoefHiCentroid =       0.5*a_lengthFraction;
                  Real hiSideCoefHiCentroid = 1.0 - 0.5*a_lengthFraction;
                  
                  IntVect loVofHiCentroid = iv               ;
                  IntVect hiVofHiCentroid = iv + BASISV(xDir); 
                                   
                  a_vector[irow    ] = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                  a_vector[irow + 1] = loSideCoefHiCentroid*fab(loVofHiCentroid,0) + hiSideCoefHiCentroid*fab(hiVofHiCentroid,0);
                }
              else
                {
                   // // loCentroid interp info
                  Real loSideCoefLoCentroid = 0.5*(1.0 + a_lengthFraction); 
                  Real hiSideCoefLoCentroid = 0.5*(1.0 - a_lengthFraction);
                                    
                  IntVect loVofLoCentroid = iv - BASISV(xDir);
                  IntVect hiVofLoCentroid = iv                ;
                                    
                  // hiCentroid extrapolation info
                  IntVect loVofHiCentroid = iv + BASISV(xDir)               ;
                  IntVect hiVofHiCentroid = iv + BASISV(xDir) + BASISV(xDir);             ;
                  Real extrapP            = a_dx*(iv[0] + 0.5*(1.0 + a_lengthFraction)); 
                  Real pLo                = a_dx*(iv[0] + 1 + 0.5);
                  Real valPLo             = fab(loVofHiCentroid,0);
                  Real valPHi             = fab(hiVofHiCentroid,0);
                  Real distPLoPHi         = a_dx;
                               
                  a_vector[irow    ] = loSideCoefLoCentroid*fab(loVofLoCentroid,0) + hiSideCoefLoCentroid*fab(hiVofLoCentroid,0);
                  a_vector[irow + 1] = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                }
            }
          else
            {
              a_vector[irow + 1] = fab(iv);
            }
         }
    }
}
void PicardSolver::setVector(Vector<Real>               & a_vector,
                             const LevelData<FluxBox> & a_levelData,
                             const IntVect              & a_groundingLineIv,
                             const Real                 & a_lengthFraction,
                             const Real                 & a_dx)
{
  int maxLevel = 0;
  int xDir = 0;

  DisjointBoxLayout dbl = m_vectGrids[maxLevel];
  Box domainBox         = m_domains  [maxLevel].domainBox();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FluxBox& fab = a_levelData[dit()];

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit)
        {
          // IntVect
          const IntVect& iv = bit();
          const int irow = iv[0];
           
          // assign rhs
          if (iv[0] < a_groundingLineIv[0])
            {
              a_vector[irow] = fab[xDir](iv,0);
            }
          else if (iv[0] == a_groundingLineIv[0])
            {
              if (a_lengthFraction < 0.5)
                {
                  // loCentroid extrapolation info
                  IntVect loVofLoCentroid = iv - BASISV(xDir) - BASISV(xDir);
                  IntVect hiVofLoCentroid = iv - BASISV(xDir)               ;
                  Real extrapP            = a_dx*(iv[0] + 0.5*a_lengthFraction); 
                  Real pLo                = a_dx*(iv[0] - 2 + 0.5);
                  Real valPLo             = fab[xDir](loVofLoCentroid,0);
                  Real valPHi             = fab[xDir](hiVofLoCentroid,0);
                  Real distPLoPHi         = a_dx;
                           
                  // hiCentroid interp info
                  Real loSideCoefHiCentroid =       0.5*a_lengthFraction;
                  Real hiSideCoefHiCentroid = 1.0 - 0.5*a_lengthFraction;
                  
                  IntVect loVofHiCentroid = iv               ;
                  IntVect hiVofHiCentroid = iv + BASISV(xDir); 
                                   
                  a_vector[irow    ] = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                  a_vector[irow + 1] = loSideCoefHiCentroid*fab[xDir](loVofHiCentroid,0) + hiSideCoefHiCentroid*fab[xDir](hiVofHiCentroid,0);
                }
              else
                {
                   // // loCentroid interp info
                  Real loSideCoefLoCentroid = 0.5*(1.0 + a_lengthFraction); 
                  Real hiSideCoefLoCentroid = 0.5*(1.0 - a_lengthFraction);
                                    
                  IntVect loVofLoCentroid = iv - BASISV(xDir);
                  IntVect hiVofLoCentroid = iv                ;
                                    
                  // hiCentroid extrapolation info
                  IntVect loVofHiCentroid = iv + BASISV(xDir)               ;
                  IntVect hiVofHiCentroid = iv + BASISV(xDir) + BASISV(xDir);             ;
                  Real extrapP            = a_dx*(iv[0] + 0.5*(1.0 + a_lengthFraction)); 
                  Real pLo                = a_dx*(iv[0] + 1 + 0.5);
                  Real valPLo             = fab[xDir](loVofHiCentroid,0);
                  Real valPHi             = fab[xDir](hiVofHiCentroid,0);
                  Real distPLoPHi         = a_dx;
                               
                  a_vector[irow    ] = loSideCoefLoCentroid*fab[xDir](loVofLoCentroid,0) + hiSideCoefLoCentroid*fab[xDir](hiVofLoCentroid,0);
                  a_vector[irow + 1] = linearExtrapolation(extrapP,pLo,valPLo,valPHi,distPLoPHi);
                }
            }
          else
            {
              a_vector[irow + 1] = fab[xDir](iv);
            }
         }
    }
}

void PicardSolver::outputVector(const Vector<LevelData<FArrayBox>* > & a_vecLDPtr,
                                const IntVect                        & a_groundingLineIv,
                                const Real                           & a_lengthFraction,
                                const Real                           & a_dx)
{
  int maxLevel = 0;
  const LevelData<FArrayBox>* LDPtr = a_vecLDPtr[maxLevel];
  
  // domain info
  DisjointBoxLayout dbl = m_vectGrids[maxLevel];
  Box domainBox         = m_domains  [maxLevel].domainBox();
  
  int xDir = 0;
  int rows = domainBox.size()[xDir] + 1;
  
  // rhs for the matrix equation
  Vector<Real> vec(rows);
                   
  // cell-centered data
  const LevelData<FArrayBox>& LDRef = *LDPtr;
  
  // create vectors from Fabs
  setVector(vec,LDRef,a_groundingLineIv,a_lengthFraction,a_dx);
  for(int irow = 0; irow <rows; ++irow)
    {
      pout()<<"vec["<<irow<<"] = "<<vec[irow]<<endl;
    }

  pout()<<endl;
}

void PicardSolver::outputVector(const Vector<LevelData<FluxBox>* > & a_vecLDPtr,
                                const IntVect                      & a_groundingLineIv,
                                const Real                         & a_lengthFraction,
                                const Real                         & a_dx)
{
  int maxLevel = 0;
  const LevelData<FluxBox  >* LDPtr = a_vecLDPtr[maxLevel];
  
  // domain info
  DisjointBoxLayout dbl = m_vectGrids[maxLevel];
  Box domainBox         = m_domains  [maxLevel].domainBox();
  
  int xDir = 0;
  int rows = domainBox.size()[xDir] + 1;
  
  // rhs for the matrix equation
  Vector<Real> vec(rows);
                   
  // cell-centered data
  const LevelData<FluxBox>& LDRef = *LDPtr;
  
  // create vectors from Fabs
  setVector(vec,LDRef,a_groundingLineIv,a_lengthFraction,a_dx);
  for(int irow = 0; irow <rows; ++irow)
    {
      pout()<<"vec["<<irow<<"] = "<<vec[irow]<<endl;
    }
}
#endif
#endif

  
#include "NamespaceFooter.H"

