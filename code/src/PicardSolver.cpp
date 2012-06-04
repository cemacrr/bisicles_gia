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
                                                                                IntVect::Zero) );


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
int
PicardSolver::solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		    Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		    const Real a_convergenceMetric,
                    const Vector<LevelData<FArrayBox>* >& a_rhs,
                    const Vector<LevelData<FArrayBox>* >& a_beta,
		    const Vector<LevelData<FArrayBox>* >& a_beta0, // not used
		    const Vector<LevelData<FArrayBox>* >& a_A,
		    const Vector<LevelData<FluxBox>* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel)
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
      faceA[lev] = new LevelData<FluxBox>
	(m_vectGrids[lev], a_A[lev]->nComp(), IntVect::Zero);
      CellToEdge(*a_A[lev] , *faceA[lev]);
    }
 
  // copy beta into local storage (also not terribly efficient; we'll
  // fix that later as well)
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
        
  //if (!m_isOpDefined)
  //  {
  //    defineOpFactory();
  //  }

  
  computeMu(a_horizontalVel, faceA, a_coordSys, a_time);
  
  Real& residNorm = a_finalResidualNorm;
  Real& initialResid = a_initialResidualNorm;
  int numIter = 0;


  if (m_isSolverDefined)
    {
      // simplest thing to do here is to clear the solvers
      // come back later and implement solver re-use
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
  // compute initial residual
  // note that whichever solver we're using, we've got the
  // AMRMultigrid solver, so can use that to compute the residual
  //bool homoBCs = false;
  
  // create local storage for AMR residual to make it easy to write 
  // the AMR residual to a plotfile if needed

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
  
  Real convergenceMetric = (a_convergenceMetric > 0.0)?a_convergenceMetric:initialResid;

  CH_assert(residNorm < HUGE_NORM);
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
              m_MGSolverPtr->solveNoInit(a_horizontalVel, a_rhs,
                                         a_maxLevel, a_lbase, zeroPhi);
            }
          else 
            {
              // in case horizontalVel has more levels than we're 
              // currently solving for (common for the case where
              // we're not refined all the way to the maximum allowable
              // level, create local Vectors with only the levels we're
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
              else
                {
                  MayDay::Error("Invalid solver type");
                }
            }

          // recompute mu
          computeMu( a_horizontalVel, faceA, a_coordSys,  a_time );
          
          // compute new residual -- take advantage of the fact that 
          // mg solver is using the RefCountedPtr's of the coefficients, so 
          // we can just call computeAMRResidual w/o recomputing anything
          m_MGSolverPtr->m_convergenceMetric = initialResid;
          Real oldResidNorm = residNorm;

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
  
  // clean up temporary storage
  for (int lev=0; lev<amrResid.size(); lev++)
    {
      if (amrResid[lev] != NULL)
        {
          delete amrResid[lev];
          amrResid[lev] = NULL;
        }
      if (faceA[lev] != NULL)
	{
	  delete faceA[lev];
	  faceA[lev] = NULL;
	}
    }

  return returnCode;
}
  
  

void
PicardSolver::setDefaultValues()
{
  m_solver_tolerance = 1.0e-10;
  m_absolute_tolerance = 5e-10;
  m_max_iter = 100;
  m_solver_type = multigrid;
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
  else
    {
      MayDay::Error("invalid solver type");
    }  

  m_isSolverDefined = true;
}


void
PicardSolver::clearLinearSolvers()
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
  if (SpaceDim == 2)
    {
      CH_assert(m_opFactoryPtr == NULL);
     
      Real alpha, beta;
      getOperatorScaleFactors(alpha, beta);

      // for the moment, at least, this only works for dx = dy:
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
PicardSolver::computeMu(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
			Vector<LevelData<FluxBox>* >& a_A, 
                        Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                        Real a_time)
{
  Vector<LevelData<FArrayBox>* >& avgDownVel = a_horizontalVel;
    
  ProblemDomain levelDomain = m_coarseDomain;
  for (int lev=0; lev<m_numLevels; lev++)

    {
      const DisjointBoxLayout& levelGrids = m_vectGrids[lev];
      const LevelSigmaCS& levelCS = *a_coordSys[lev];
      LevelData<FArrayBox>& levelVel = *avgDownVel[lev];
      LevelData<FluxBox>& levelMu = *m_vectMu[lev];
      LevelData<FluxBox>& levelA = *a_A[lev];
      LevelData<FluxBox>& levelLambda = *m_vectLambda[lev];
      const LevelData<FArrayBox>& levelBeta = *m_vectBeta[lev];
      LevelData<FArrayBox>& levelC = *m_vectC[lev];



      // // first thing, if there is a finer level, average-down
      // // the current velocity field
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
      //Real dxLevel = 1;
      Real dxLevel = levelCS.dx()[0];

      if (lev > 0) 
        {
          QuadCFInterp qcfi(levelGrids, &m_vectGrids[lev-1],
                            dxLevel, m_vectRefRatio[lev-1], 
                            2, levelDomain);
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

      // now multiply by ice thickness H
      const LevelData<FluxBox>& faceH = levelCS.getFaceH();
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          
          for (int dir=0; dir<SpaceDim; dir++)
            {
              levelMu[dit][dir].mult(faceH[dit][dir],
                                     levelMu[dit][dir].box(),0,0,1);
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
	  m_basalFrictionRelPtr->computeAlpha
	    (levelC[dit], levelVel[dit], levelCS.getThicknessOverFlotation()[dit], levelBeta[dit] ,
	     levelCS.getFloatingMask()[dit],gridBox);

        }
      if (lev != m_numLevels-1) levelDomain.refine(m_vectRefRatio[lev]);
    }
}


// // compute face-centered coefficients for tensor solver (really
// // winds up being H*mu) -- non-isothermal version...
// void
// PicardSolver::computeMu(Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
//                         Vector<LevelData<FArrayBox>* >& a_A, 
//                         Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys, 
//                         Real a_time)
// {
  
// }

#include "NamespaceFooter.H"

