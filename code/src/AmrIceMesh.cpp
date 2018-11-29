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
#include "BISICLES_VERSION.H"
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
#include "IceThermodynamics.H"
#include "PicardSolver.H"
#include "JFNKSolver.H"
#include "InverseVerticallyIntegratedVelocitySolver.H"
#include "PetscIceSolver.H"
#include "RelaxSolver.H"
#ifdef CH_USE_FAS
#include "FASIceSolverI.H"
#endif
#include "KnownVelocitySolver.H"
#include "VCAMRPoissonOp2.H"
#include "AMRPoissonOpF_F.H"
#include "CH_HDF5.H"
#include "IceUtility.H"
#include "LevelMappedDerivatives.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

#include "NamespaceHeader.H"


// do regridding
void
AmrIce::regrid()
{

  CH_TIME("AmrIce::regrid");

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::regrid" << endl;
    }
  
  //first part of a conservation of volume check
  Real volumeBefore = computeTotalIce();

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

#ifdef REGRID_EH
	// We need to regrid the internal energy, but the conserved quantity is H*E
	// Set E <- E*H now (and E <- E/H later)
	for (int lev=0; lev <= m_finest_level ; ++lev)
	  {
	    for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
	      {
		FArrayBox& E = (*m_internalEnergy[lev])[dit];
		FArrayBox H(E.box(),1);
		H.copy((*m_old_thickness[lev])[dit]);
		H += 1.0e-10;
		for (int comp  = 0; comp < m_internalEnergy[0]->nComp(); comp++)
		  {
		    E.mult( H,0,comp,1);
		  }
	      }
	  }
#endif

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
	    LevelData<FArrayBox>* old_old_thicknessDataPtr = m_old_thickness[lev];
	    LevelData<FArrayBox>* old_velDataPtr = m_velocity[lev];
	    LevelData<FArrayBox>* old_tempDataPtr = m_internalEnergy[lev];
	    LevelData<FArrayBox>* old_calvDataPtr = m_calvedIceThickness[lev];
	    LevelData<FArrayBox>* old_removedDataPtr = m_removedIceThickness[lev];
	    LevelData<FArrayBox>* old_addedDataPtr = m_addedIceThickness[lev];

	    LevelData<FArrayBox>* old_deltaTopographyDataPtr = m_deltaTopography[lev];
            LevelData<FArrayBox>* old_iceFracDataPtr = m_iceFrac[lev];
	     
	    LevelData<FArrayBox>* new_old_thicknessDataPtr = 
	      new LevelData<FArrayBox>(newDBL, 1, m_old_thickness[0]->ghostVect());
	      
	    LevelData<FArrayBox>* new_velDataPtr = 
	      new LevelData<FArrayBox>(newDBL, SpaceDim, m_velocity[0]->ghostVect());

	    LevelData<FArrayBox>* new_tempDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_internalEnergy[0]->nComp(), 
				       m_internalEnergy[0]->ghostVect());
	    //since the internalEnergy data has changed
	    m_A_valid = false;

            LevelData<FArrayBox>* new_iceFracDataPtr = 
              new LevelData<FArrayBox>(newDBL, 1, m_iceFrac[0]->ghostVect());
            
	    LevelData<FArrayBox>* new_calvDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_calvedIceThickness[0]->nComp(), 
				       m_calvedIceThickness[0]->ghostVect());
	    LevelData<FArrayBox>* new_removedDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_removedIceThickness[0]->nComp(), 
				       m_removedIceThickness[0]->ghostVect());
	    LevelData<FArrayBox>* new_addedDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_addedIceThickness[0]->nComp(), 
				       m_addedIceThickness[0]->ghostVect());
	      
	    LevelData<FArrayBox>* new_deltaTopographyDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_deltaTopography[0]->nComp(), 
				       m_deltaTopography[0]->ghostVect());

	      
#if BISICLES_Z == BISICLES_LAYERED
	    LevelData<FArrayBox>* old_sTempDataPtr = m_sInternalEnergy[lev];
	    LevelData<FArrayBox>* old_bTempDataPtr = m_bInternalEnergy[lev];
	    LevelData<FArrayBox>* new_sTempDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_sInternalEnergy[0]->nComp(),
				       m_sInternalEnergy[0]->ghostVect());
	    LevelData<FArrayBox>* new_bTempDataPtr = 
	      new LevelData<FArrayBox>(newDBL, m_bInternalEnergy[0]->nComp(),
				       m_bInternalEnergy[0]->ghostVect());
	     
#endif	      
	      
	    {
	      // first we need to regrid m_deltaTopography, it will be needed to 
	      // regrid the bedrock topography & hence LevelSigmaCS
	      // may eventually want to do post-regrid smoothing on this
	      FineInterp interpolator(newDBL,m_deltaTopography[0]->nComp(),
				      m_refinement_ratios[lev-1],
				      m_amrDomains[lev]);
	      interpolator.interpToFine(*new_deltaTopographyDataPtr, *m_deltaTopography[lev-1]);

		
	      PiecewiseLinearFillPatch ghostFiller
		(m_amrGrids[lev],
		 m_amrGrids[lev-1],
		 m_deltaTopography[lev-1]->nComp(),
		 m_amrDomains[lev-1],
		 m_refinement_ratios[lev-1],
		 m_deltaTopography[lev-1]->ghostVect()[0]);

	      ghostFiller.fillInterp(*new_deltaTopographyDataPtr,*m_deltaTopography[lev-1],
				     *m_deltaTopography[lev-1],1.0,0,0,
				     m_deltaTopography[lev-1]->nComp());

	      if (old_deltaTopographyDataPtr != NULL && oldDBL.isClosed())
		{
		  old_deltaTopographyDataPtr->copyTo(*new_deltaTopographyDataPtr);
		}
	      delete old_deltaTopographyDataPtr;
		
	    }



	    // also need to handle LevelSigmaCS 

	    // assume level 0 has correct ghosting
	    IntVect sigmaCSGhost = m_vect_coordSys[0]->ghostVect();
	    {
	      RealVect dx = m_amrDx[lev]*RealVect::Unit;
	      RefCountedPtr<LevelSigmaCS > oldCoordSys = m_vect_coordSys[lev];
		
	      RefCountedPtr<LevelSigmaCS > auxCoordSys = (lev > 0)?m_vect_coordSys[lev-1]:oldCoordSys;

	      m_vect_coordSys[lev] = RefCountedPtr<LevelSigmaCS >
		(new LevelSigmaCS(newDBL, dx, sigmaCSGhost));
	      m_vect_coordSys[lev]->setIceDensity(auxCoordSys->iceDensity());
	      m_vect_coordSys[lev]->setWaterDensity(auxCoordSys->waterDensity());
	      m_vect_coordSys[lev]->setGravity(auxCoordSys->gravity());
	      m_vect_coordSys[lev]->setBackgroundSlope(auxCoordSys->getBackgroundSlope());
#if BISICLES_Z == BISICLES_LAYERED
	      m_vect_coordSys[lev]->setFaceSigma(auxCoordSys->getFaceSigma());
#endif		
	      LevelSigmaCS* crsePtr = &(*m_vect_coordSys[lev-1]);
	      int refRatio = m_refinement_ratios[lev-1];

	      bool interpolate_zb = (m_interpolate_zb ||
				     !m_thicknessIBCPtr->regridIceGeometry
				     (*m_vect_coordSys[lev],dx,  m_domainSize, 
				      m_time,  crsePtr,refRatio ) );
		
	      if (!interpolate_zb)
		{
		  // need to re-apply accumulated bedrock (GIA). Could be optional?
		  for (DataIterator dit(newDBL); dit.ok(); ++dit)
		    {
		      m_vect_coordSys[lev]->getTopography()[dit] += (*new_deltaTopographyDataPtr)[dit];
		    }
		}

		{
		  //interpolate thickness & (maybe) topography
		  bool interpolateThickness(true);
		  bool preserveMask(true);
		  bool interpolateTopographyGhost(true); 
		  bool interpolateThicknessGhost(true); 
		  bool preserveMaskGhost(true);
		  m_vect_coordSys[lev]->interpFromCoarse(*m_vect_coordSys[lev-1],
							 m_refinement_ratios[lev-1],
							 interpolate_zb,
							 interpolateThickness, 
							 preserveMask,
							 interpolateTopographyGhost, 
							 interpolateThicknessGhost, 
							 preserveMaskGhost, 
							 m_regrid_thickness_interpolation_method);
		}


	      LevelData<FArrayBox>& thisLevelH = m_vect_coordSys[lev]->getH();
	      LevelData<FArrayBox>& thisLevelB = m_vect_coordSys[lev]->getTopography();
		
	      // overwrite interpolated fields in valid regiopns with such valid old data as there is
	      if (oldDBL.isClosed()){	  
		const LevelData<FArrayBox>& oldLevelH = oldCoordSys->getH();
		oldLevelH.copyTo(thisLevelH);
		const LevelData<FArrayBox>& oldLevelB = oldCoordSys->getTopography();
		oldLevelB.copyTo(thisLevelB);
	      }

	      //Defer to m_thicknessIBCPtr for boundary values - 
	      //interpolation won't cut the mustard because it only fills
	      //ghost cells overlying the valid regions.
	      RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
	      m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
						m_amrDomains[lev],levelDx, m_time, m_dt);



	      // exchange is necessary to fill periodic ghost cells
	      // which aren't filled by the copyTo from oldLevelH
	      thisLevelH.exchange();
	      m_vect_coordSys[lev]->exchangeTopography();

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
	    
	      interpolator.interpToFine(*new_old_thicknessDataPtr, *m_old_thickness[lev-1]);
	
	      // now copy old-grid data into new holder
	      if (old_old_thicknessDataPtr != NULL) 
		{
		  if ( oldDBL.isClosed())
		    {
		      old_old_thicknessDataPtr->copyTo(*new_old_thicknessDataPtr);
		    }
		  
		}

		interpolator.interpToFine(*new_iceFracDataPtr, *m_iceFrac[lev-1]);
	
		// now copy old-grid data into new holder
		if (old_iceFracDataPtr != NULL) 
		  {
		    if ( oldDBL.isClosed())
		      {
			old_iceFracDataPtr->copyTo(*new_iceFracDataPtr);
		      }
		    // can now delete old data 
		    delete old_iceFracDataPtr;
		  }

		
	    }
	      
	    {
	      // may eventually want to do post-regrid smoothing on this!
	      FineInterp interpolator(newDBL, SpaceDim, 
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

	      //handle ghost cells on the coarse-fine interface
	      QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
				m_amrDx[lev], m_refinement_ratios[lev-1],
				2, m_amrDomains[lev]);
	      qcfi.coarseFineInterp(*new_velDataPtr, *m_velocity[lev-1]);
		
	      //boundary ghost cells
	      m_thicknessIBCPtr->velocityGhostBC
		(*new_velDataPtr, *m_vect_coordSys[lev],m_amrDomains[lev],m_time);
		

	    }

	    {
	      // may eventually want to do post-regrid smoothing on this
	    
	      FineInterp interpolator(newDBL,m_internalEnergy[0]->nComp(),
				      m_refinement_ratios[lev-1],
				      m_amrDomains[lev]);
	      interpolator.interpToFine(*new_tempDataPtr, *m_internalEnergy[lev-1]);

		
	      PiecewiseLinearFillPatch ghostFiller
		(m_amrGrids[lev],
		 m_amrGrids[lev-1],
		 m_internalEnergy[lev-1]->nComp(),
		 m_amrDomains[lev-1],
		 m_refinement_ratios[lev-1],
		 m_internalEnergy[lev-1]->ghostVect()[0]);

	      ghostFiller.fillInterp(*new_tempDataPtr,*m_internalEnergy[lev-1],
				     *m_internalEnergy[lev-1],1.0,0,0,
				     m_internalEnergy[lev-1]->nComp());

	      if (old_tempDataPtr != NULL && oldDBL.isClosed())
		{
		  old_tempDataPtr->copyTo(*new_tempDataPtr);
		}
	      delete old_tempDataPtr;	
	    }
	      
	   
	    {
	      FineInterp interpolator(newDBL,m_calvedIceThickness[0]->nComp(),
				      m_refinement_ratios[lev-1],
				      m_amrDomains[lev]);
	      interpolator.interpToFine(*new_calvDataPtr, *m_calvedIceThickness[lev-1]);
		
	      PiecewiseLinearFillPatch ghostFiller
		(m_amrGrids[lev],
		 m_amrGrids[lev-1],
		 m_calvedIceThickness[lev-1]->nComp(),
		 m_amrDomains[lev-1],
		 m_refinement_ratios[lev-1],
		 m_calvedIceThickness[lev-1]->ghostVect()[0]);

	      ghostFiller.fillInterp(*new_calvDataPtr,*m_calvedIceThickness[lev-1],
				     *m_calvedIceThickness[lev-1],1.0,0,0,
				     m_calvedIceThickness[lev-1]->nComp());

	      if (old_calvDataPtr != NULL && oldDBL.isClosed())
		{
		  old_calvDataPtr->copyTo(*new_calvDataPtr);
		}
	      delete old_calvDataPtr;
		
	    }
	    {
	      FineInterp interpolator(newDBL,m_removedIceThickness[0]->nComp(),
				      m_refinement_ratios[lev-1],
				      m_amrDomains[lev]);
	      interpolator.interpToFine(*new_removedDataPtr, *m_removedIceThickness[lev-1]);
		
	      PiecewiseLinearFillPatch ghostFiller
		(m_amrGrids[lev],
		 m_amrGrids[lev-1],
		 m_removedIceThickness[lev-1]->nComp(),
		 m_amrDomains[lev-1],
		 m_refinement_ratios[lev-1],
		 m_removedIceThickness[lev-1]->ghostVect()[0]);

	      ghostFiller.fillInterp(*new_removedDataPtr,*m_removedIceThickness[lev-1],
				     *m_removedIceThickness[lev-1],1.0,0,0,
				     m_removedIceThickness[lev-1]->nComp());

	      if (old_removedDataPtr != NULL && oldDBL.isClosed())
		{
		  old_removedDataPtr->copyTo(*new_removedDataPtr);
		}
	      delete old_removedDataPtr;
		
	    }
	    {
	      FineInterp interpolator(newDBL,m_addedIceThickness[0]->nComp(),
				      m_refinement_ratios[lev-1],
				      m_amrDomains[lev]);
	      interpolator.interpToFine(*new_addedDataPtr, *m_addedIceThickness[lev-1]);
		
	      PiecewiseLinearFillPatch ghostFiller
		(m_amrGrids[lev],
		 m_amrGrids[lev-1],
		 m_addedIceThickness[lev-1]->nComp(),
		 m_amrDomains[lev-1],
		 m_refinement_ratios[lev-1],
		 m_addedIceThickness[lev-1]->ghostVect()[0]);

	      ghostFiller.fillInterp(*new_addedDataPtr,*m_addedIceThickness[lev-1],
				     *m_addedIceThickness[lev-1],1.0,0,0,
				     m_addedIceThickness[lev-1]->nComp());

	      if (old_addedDataPtr != NULL && oldDBL.isClosed())
		{
		  old_addedDataPtr->copyTo(*new_addedDataPtr);
		}
	      delete old_addedDataPtr;
		
	    }

 

#if BISICLES_Z == BISICLES_LAYERED
	    {
	      // may eventually want to do post-regrid smoothing on this
	      FineInterp interpolator(newDBL,m_sInternalEnergy[0]->nComp(),
				      m_refinement_ratios[lev-1],
				      m_amrDomains[lev]);

	      PiecewiseLinearFillPatch ghostFiller
		(m_amrGrids[lev],
		 m_amrGrids[lev-1],
		 m_sInternalEnergy[lev-1]->nComp(),
		 m_amrDomains[lev-1],
		 m_refinement_ratios[lev-1],
		 m_sInternalEnergy[lev-1]->ghostVect()[0]);

	

	      interpolator.interpToFine(*new_sTempDataPtr, *m_sInternalEnergy[lev-1]);
		
	      ghostFiller.fillInterp(*new_sTempDataPtr,*m_sInternalEnergy[lev-1],
				     *m_sInternalEnergy[lev-1],1.0,0,0,
				     m_sInternalEnergy[lev-1]->nComp());


	      if (old_sTempDataPtr != NULL && oldDBL.isClosed())
		{
		  old_sTempDataPtr->copyTo(*new_sTempDataPtr);
		}
	      delete old_sTempDataPtr;
		
	      interpolator.interpToFine(*new_bTempDataPtr, *m_bInternalEnergy[lev-1]);

	      ghostFiller.fillInterp(*new_bTempDataPtr,*m_bInternalEnergy[lev-1],
				     *m_bInternalEnergy[lev-1],1.0,0,0,
				     m_bInternalEnergy[lev-1]->nComp());
		
	      if (old_bTempDataPtr != NULL && oldDBL.isClosed())
		{
		  old_bTempDataPtr->copyTo(*new_bTempDataPtr);
		}
	      delete old_bTempDataPtr;

	      new_tempDataPtr->exchange();
	      new_sTempDataPtr->exchange();
	      new_bTempDataPtr->exchange();
	      //set boundary for non-periodic cases values
	      m_internalEnergyIBCPtr->setIceInternalEnergyBC
		(*new_tempDataPtr,*new_sTempDataPtr,*new_bTempDataPtr,
		 *m_vect_coordSys[lev] );

	    }
#endif
	    // can now delete old data 
	    delete old_old_thicknessDataPtr;

	    // now copy new holders into multilevel arrays
	    m_old_thickness[lev] = new_old_thicknessDataPtr;
	    m_velocity[lev] = new_velDataPtr;
	    m_internalEnergy[lev] = new_tempDataPtr;
	    m_calvedIceThickness[lev] = new_calvDataPtr;
	    m_removedIceThickness[lev] = new_removedDataPtr;
	    m_addedIceThickness[lev] = new_addedDataPtr;
	    m_deltaTopography[lev] = new_deltaTopographyDataPtr;
            m_iceFrac[lev] = new_iceFracDataPtr;      
#if BISICLES_Z == BISICLES_LAYERED
	    m_sInternalEnergy[lev] = new_sTempDataPtr;
	    m_bInternalEnergy[lev] = new_bTempDataPtr;
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
	     
	    if (m_velRHS[lev] != NULL)
	      {
		delete m_velRHS[lev];
	      }
	    m_velRHS[lev] = new LevelData<FArrayBox>(newDBL, SpaceDim, 
						     IntVect::Zero);

	    if (m_faceVelAdvection[lev] != NULL)
	      {
		delete m_faceVelAdvection[lev];
	      }
	    m_faceVelAdvection[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);

	    if (m_faceVelTotal[lev] != NULL)
	      {
		delete m_faceVelTotal[lev];
	      }
	    m_faceVelTotal[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);


	    if (m_diffusivity[lev] != NULL)
	      {
		delete m_diffusivity[lev];
	      }
	    m_diffusivity[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);


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


	    if (m_divThicknessFlux[lev] != NULL)
	      {
		delete m_divThicknessFlux[lev];
	      }
	    m_divThicknessFlux[lev] = 
	      new LevelData<FArrayBox>(newDBL,   1, IntVect::Zero) ;


	    if (m_bHeatFlux[lev] != NULL)
	      {
		delete m_bHeatFlux[lev];
	      }
	    m_bHeatFlux[lev] = 
	      new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit);

	    if (m_sHeatFlux[lev] != NULL)
	      {
		delete m_sHeatFlux[lev];
	      }
	    m_sHeatFlux[lev] = 
	      new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit);



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
	      
	    if (m_internalEnergy[lev] != NULL) 
	      {
		delete m_internalEnergy[lev];
		m_internalEnergy[lev] = NULL;
	      }

	      if (m_iceFrac[lev] != NULL) 
		{
		  delete m_iceFrac[lev];
		  m_iceFrac[lev] = NULL;
		}

#if BISICLES_Z == BISICLES_LAYERED
	    if (m_sInternalEnergy[lev] != NULL) 
	      {
		delete m_sInternalEnergy[lev];
		m_sInternalEnergy[lev] = NULL;
	      }
	    if (m_bInternalEnergy[lev] != NULL) 
	      {
		delete m_bInternalEnergy[lev];
		m_bInternalEnergy[lev] = NULL;
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

	// this is a good time to check for remote ice
	if ((m_eliminate_remote_ice_after_regrid) 
	    && !(m_eliminate_remote_ice))
	  eliminateRemoteIce();
#ifdef REGRID_EH 	
	// Since we set E <- E*H earlier, set E <- E/H later now
	for (int lev=0; lev<= new_finest_level; ++lev)
	  {
	    for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
	      {
		FArrayBox& E = (*m_internalEnergy[lev])[dit];
		FArrayBox H(E.box(),1);
		H.copy((*m_old_thickness[lev])[dit]);
		H += 1.0e-10;
		for (int comp  = 0; comp < m_internalEnergy[0]->nComp(); comp++)
		  {
		    E.divide( H,0,comp,1);
		  }
	      }
	  }
#endif

	
	//applyCalvingCriterion(CalvingModel::PostRegrid);

	if (m_evolve_velocity)
	  {
	    //velocity solver needs to be re-defined
	    defineSolver();
	    //solve velocity field, but use the previous initial residual norm in place of this one
	    //and force a solve even if other conditions (e.g the timestep interval condition) are not met
	    solveVelocityField(true, m_velocitySolveInitialResidualNorm);
	  }
	else
	  {
	    CH_assert(m_evolve_velocity);
	    MayDay::Error("AmrIce::regrid() not implemented for !m_evolve_velocity");
	  }
         
	  
      } // end if tags changed
    } // end if max level > 0 in the first place
  
  Real volumeAfter = computeTotalIce();
  Real volumeDifference = volumeAfter - volumeBefore;
  if (s_verbosity > 3) 
    { 
      
      pout() << "AmrIce::regrid: volume on input,output,difference =  " 
	     << volumeBefore << "," << volumeAfter << "," << volumeDifference << " m^3" << endl;
    }


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
          FArrayBox lapVel(levelGrids[dit()], SpaceDim);
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
  
  // tag on |vel| * dx > m_velDx_tagVal. This style of tagging is used in the AMRControl.cpp
  // (but with the observed field), so can be used to construct similar meshes
  if (m_tagVelDx)
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
	  const Box& box = levelGrids[dit()];
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	  const FArrayBox& vel = levelVel[dit];
	  for (BoxIterator bit(box) ; bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if ((mask(iv) == GROUNDEDMASKVAL && a_level < m_velDx_tagVal_finestLevelGrounded) ||
		  (mask(iv) == FLOATINGMASKVAL && a_level < m_velDx_tagVal_finestLevelFloating) )
		{
		  Real v = 0.0;
		  for (int dir = 0; dir < SpaceDim; ++dir)
		    {
		      v += std::pow(vel(iv,dir),2);
		    }
		  if (sqrt(v)*dx(a_level)[0] > m_velDx_tagVal)
		    {
		      local_tags |= iv;
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
  
  // tag cells  with thin cavities
  if (m_tag_thin_cavity)
    {
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	  Box gridBox = levelGrids[dit];
	  const FArrayBox& H  = levelCS.getH()[dit];
	  for (BoxIterator bit(gridBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if (mask(iv) == FLOATINGMASKVAL && 
		  H(iv) <  m_tag_thin_cavity_thickness)
		{
		  local_tags |= iv;
		}
	    }
	}
    }

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
	      if ( a_level < m_margin_tagVal_finestLevel )
		{
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
		}
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

  


  // tag anywhere and everywhere
  if (m_tagEntireDomain)
    {
      // this is super-simple...
      Box domainBox = m_amrDomains[a_level].domainBox();
      local_tags |= domainBox;
          
    } // end if tag entire domain



#ifdef HAVE_PYTHON
  if (m_tagPython)
    {
      //tag via a python function f(x,y,dx,H,R) (more args to come)
      // 
      Vector<Real> args(SpaceDim + 3);
      Vector<Real> rval;
      for (dit.begin(); dit.ok(); ++dit)
        {
	  for (BoxIterator bit(levelGrids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      int i = 0;
	      for (int dir=0; dir < SpaceDim; dir++)
		{
		  args[i++] = (Real(iv[dir]) + 0.5)*m_amrDx[a_level];
		}
	      args[i++]  = m_amrDx[a_level];
	      args[i++] = levelCS.getH()[dit](iv,0);
	      args[i++] = levelCS.getTopography()[dit](iv,0);				    
	      PythonInterface::PythonEval(m_tagPythonFunction, rval,  args);
	      if (rval[0] > 0.0)
		local_tags |= iv;
	    } // end bit loop
	} // end loop over boxes
    } // end if tag via python
#endif

  // now buffer tags
  
  local_tags.grow(m_tags_grow);
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      if (m_tags_grow_dir[dir] > m_tags_grow)
	local_tags.grow(dir, std::max(0,m_tags_grow_dir[dir]-m_tags_grow));
    }
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

  CH_TIME("AmrIce::initGrids");
  
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

  if (s_verbosity > 3) 
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

  bool moreLevels = (m_max_level > 0);
  int baseLevel = 0;
  
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
    


#include "NamespaceFooter.H"
