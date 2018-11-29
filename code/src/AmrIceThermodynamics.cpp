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


#if BISICLES_Z == BISICLES_LAYERED
void AmrIce::updateInternalEnergy(Vector<LevelData<FluxBox>* >& a_layerEH_half, 
				  Vector<LevelData<FluxBox>* >& a_layerH_half,
				  const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
				  const Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel,
				  const Real a_dt, const Real a_time,
				  Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysNew,
				  Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysOld,
				  const Vector<LevelData<FArrayBox>*>& a_surfaceThicknessSource,
				  const Vector<LevelData<FArrayBox>*>& a_basalThicknessSource)
{

  CH_TIME("AmrIce::updateInternalEnergy");
  //update the internalEnergy fields, 2D case
  Vector<LevelData<FluxBox>* > vectLayerFluxes(m_finest_level+1, NULL);
  Vector<LevelData<FluxBox>* > vectLayerThicknessFluxes(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectUSigma(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectDivUHxy(m_finest_level+1, NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      LevelData<FluxBox>& levelXYFaceXYVel = *a_layerXYFaceXYVel[lev];
      LevelData<FluxBox>& levelFaceEH = *a_layerEH_half[lev];
      LevelData<FluxBox>& levelFaceH = *a_layerH_half[lev];
      IntVect ghostVect = IntVect::Unit;//CoarseAverageFace requires a ghost cell

      vectUSigma[lev] = new LevelData<FArrayBox>
	(m_amrGrids[lev], m_nLayers + 1 , IntVect::Zero);
       
      vectDivUHxy[lev] = new LevelData<FArrayBox>
	(m_amrGrids[lev], m_nLayers + 1 , ghostVect);
     
      vectLayerFluxes[lev] = new LevelData<FluxBox>
	(m_amrGrids[lev], levelXYFaceXYVel.nComp() , ghostVect);

      vectLayerThicknessFluxes[lev] = new LevelData<FluxBox>
	(m_amrGrids[lev], levelXYFaceXYVel.nComp() , ghostVect);

      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  FluxBox& faceVel = levelXYFaceXYVel[dit];
	  FluxBox& faceEH = levelFaceEH[dit];
	  FluxBox& faceH = levelFaceH[dit];
	  FluxBox& flux = (*vectLayerFluxes[lev])[dit];
	  FluxBox& thicknessFlux = (*vectLayerThicknessFluxes[lev])[dit];

	  const Box& gridBox = levelGrids[dit];
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      Box faceBox(gridBox);
	      faceBox.surroundingNodes(dir);
	      flux[dir].copy(faceEH[dir], faceBox);
	      flux[dir].mult(faceVel[dir], faceBox, 0, 0, faceVel[dir].nComp());

	      

	      thicknessFlux[dir].copy(faceH[dir],faceBox);
	      thicknessFlux[dir].mult(faceVel[dir], faceBox, 0, 0, faceVel[dir].nComp());
	        
	      // CH_assert(flux[dir].norm(faceBox,0,0,flux[dir].nComp()) < HUGE_NORM);
	      // CH_assert(thicknessFlux[dir].norm(faceBox,0,0,thicknessFlux[dir].nComp()) < HUGE_NORM);
		
	    }
	}
    }
  // average fine fluxes down to coarse levels
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],
				     vectLayerFluxes[lev]->nComp(), m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*vectLayerFluxes[lev-1], *vectLayerFluxes[lev]);
      faceAverager.averageToCoarse(*vectLayerThicknessFluxes[lev-1], *vectLayerThicknessFluxes[lev]);
    }
 

     
  //vertical and cross-layer velocity components (u[z] and u^sigma)
  for (int lev=0; lev <= m_finest_level; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
     
      LevelSigmaCS& levelCoordsNew = *(a_coordSysNew[lev]);
      LevelSigmaCS& levelCoordsOld = *(a_coordSysOld[lev]);
      const Vector<Real>& dSigma = levelCoordsNew.getDSigma();

      LevelData<FArrayBox> levelGradHNew(m_amrGrids[lev], SpaceDim, IntVect::Zero);
      computeCCDerivatives(levelGradHNew, levelCoordsNew.getH(), levelCoordsNew,
			   Interval(0,0),Interval(0,SpaceDim-1));
      LevelData<FArrayBox> levelGradHOld(m_amrGrids[lev], SpaceDim, IntVect::Zero);
      computeCCDerivatives(levelGradHOld, levelCoordsOld.getH(), levelCoordsOld,
			   Interval(0,0),Interval(0,SpaceDim-1));


      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const Box& box = levelGrids[dit];
	   
	  // this copy perhaps indicates layer should run faster than
	  // dir in sFaceXYVel, but for now ...
	  const FArrayBox& sFaceXYVel = (*a_layerSFaceXYVel[lev])[dit];
	  FArrayBox uX(box, m_nLayers+1);
	  FArrayBox uY(box, m_nLayers+1);
	  
	  for (int l = 0; l < m_nLayers+1; l++)
	    {
	      uX.copy(sFaceXYVel, l*SpaceDim, l);
	      uY.copy(sFaceXYVel, l*SpaceDim + 1, l);
	    }

	  FArrayBox& oldH = levelCoordsOld.getH()[dit];
	  FArrayBox& newH = levelCoordsNew.getH()[dit];
	  //cell centered thickness at t + dt/2
	  FArrayBox Hhalf(box, 1);
	  Hhalf.copy(oldH);
	  Hhalf.plus(newH);
	  Hhalf *= 0.5;

	  //cell centered grad(thickness) at t + dt/2
	  FArrayBox gradH(box, SpaceDim);
	  gradH.copy(levelGradHNew[dit]);
	  gradH.plus(levelGradHOld[dit]);
	  gradH*=0.5;
	  //cell centered grad(surface) at t + dt/2
	  FArrayBox gradS(box, SpaceDim);
	  gradS.copy(levelCoordsOld.getGradSurface()[dit]);
	  gradS.plus(levelCoordsNew.getGradSurface()[dit]);
	  gradS*=0.5;
	  //horizontal contribution to div(Hu) at cell centres, 
	  // viz d(Hu_x)/dx' + d(Hu_y)/dy'
	  FArrayBox divUHxy(box, m_nLayers);
	  {
	    divUHxy.setVal(0.0);
	    
	    const RealVect& dx = levelCoordsNew.dx(); 
	    for (int dir =0; dir < SpaceDim; dir++)
	      {
		const FArrayBox& uH = (*vectLayerThicknessFluxes[lev])[dit][dir];
		FORT_DIVERGENCE(CHF_CONST_FRA(uH),
				CHF_FRA(divUHxy),
				CHF_BOX(box),
				CHF_CONST_REAL(dx[dir]),
				CHF_INT(dir));
	      }
	    
	  }

	  //dH / dt
	  FArrayBox dHdt(box,1);  
	  dHdt.copy(newH);
	  dHdt.plus(oldH,-1.0,0,0,1);
	  dHdt *= 1.0/a_dt;
	    
	  //calculation of dS/dt assumes surface elevation is up to date
	  //in LevelSigmaCS
	  FArrayBox dSdt(box,1); 
	  dSdt.copy(levelCoordsNew.getSurfaceHeight()[dit]);
	  dSdt -= levelCoordsOld.getSurfaceHeight()[dit];
	  dSdt *= 1.0/a_dt;

	  //surface and basal thickness source
	  const FArrayBox& bts = (*a_basalThicknessSource[lev])[dit];
	  const FArrayBox& sts = (*a_surfaceThicknessSource[lev])[dit];
	  // z-component of velocity at layer faces
	  FArrayBox uZ(box,m_nLayers + 1); 
	  // sigma-componnet of velocity at layer faces
	  FArrayBox& uSigma = (*vectUSigma[lev])[dit]; 
	  // z-component of velocity at surface 
	  FArrayBox uZs(box, 1);
	  //divUHxy.setVal(0.0);
	  int nLayers = m_nLayers;
	  uSigma.setVal(0.0);

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
			   CHF_CONST_FRA1(dHdt,0),
			   CHF_CONST_FRA1(sts,0),
			   CHF_CONST_FRA1(bts,0),
			   CHF_CONST_INT(nLayers),
			   CHF_BOX(box));

	  CH_assert(uSigma.norm(0) < HUGE_NORM);
	} //end compute vertical velocity loop over boxes

      vectUSigma[lev]->exchange();

    }//end compute vertical velocity loop over levels


  // compute rhs =a_dt *(H*dissipation - div(u H T)) and update solution
  for (int lev=0; lev <= m_finest_level; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFlux = *vectLayerFluxes[lev];
      LevelSigmaCS& levelCoordsOld = *(a_coordSysOld[lev]);
      LevelSigmaCS& levelCoordsNew = *(a_coordSysNew[lev]);
      const Vector<Real>& dSigma = levelCoordsNew.getDSigma();
      //caculate dissipation due to internal stresses
      LevelData<FArrayBox> dissipation(levelGrids,m_nLayers,IntVect::Zero);
      {
	LevelData<FArrayBox>* crseVelPtr = NULL;
	int nRefCrse = -1;
	if (lev > 0)
	  {
	    crseVelPtr = m_velocity[lev-1];
	    nRefCrse = m_refinement_ratios[lev-1];
	  }

	m_velocity[lev]->exchange();

	m_constitutiveRelation->computeDissipation
	  (dissipation,*m_velocity[lev],  crseVelPtr,
	   nRefCrse, *m_A[lev],
	   levelCoordsOld , m_amrDomains[lev], IntVect::Zero);
      }
      
      LevelData<FArrayBox>& surfaceHeatFlux = *m_sHeatFlux[lev];
      if (surfaceHeatBoundaryDirichlett())
	{
	  surfaceHeatBoundaryData().evaluate(*m_sInternalEnergy[lev], *this, lev, a_dt);
	  if (surfaceHeatBoundaryTemperature())
	    {
	      //convert surface temperature to internal energy
	      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
		{
		  FArrayBox& E = (*m_sInternalEnergy[lev])[dit];
		  FArrayBox T(E.box(),1);
		  T.copy(E);
		  FArrayBox W(E.box(),1);
		  W.setVal(0.0);
		  IceThermodynamics::composeInternalEnergy(E, T, W, E.box());
		}
	    }
	}
      else
	{
	  m_surfaceHeatBoundaryDataPtr->evaluate(surfaceHeatFlux, *this, lev, a_dt);
	}
      
      LevelData<FArrayBox>& basalHeatFlux = *m_bHeatFlux[lev];
      basalHeatBoundaryData().evaluate(basalHeatFlux, *this, lev, a_dt);

      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const Box& box = levelGrids[dit];
	  const FArrayBox& oldH = levelCoordsOld.getH()[dit];
	  const FArrayBox& newH = levelCoordsNew.getH()[dit];
	  FArrayBox& E = (*m_internalEnergy[lev])[dit];
	  FArrayBox& sT = (*m_sInternalEnergy[lev])[dit];	
	  FArrayBox& bT = (*m_bInternalEnergy[lev])[dit];
	  

	  // first, do the ordinary fluxes : if we just had
	  // horizontal advection and grad(H) = grad(S) = 0., 
	  // this would be the lot
	       
	  FArrayBox rhs(box, m_nLayers);
	  rhs.setVal(0.0);
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      Real dx = levelCoordsOld.dx()[dir];              
	   
	      FORT_DIVERGENCE(CHF_CONST_FRA(levelFlux[dit][dir]),
			      CHF_FRA(rhs),
			      CHF_BOX(box),
			      CHF_CONST_REAL(dx),
			      CHF_INT(dir));
	     

	    }
	  for (int layer = 0; layer < dissipation.nComp(); ++layer)
	    {
	      dissipation[dit].mult(newH,0,layer,1);
	      
	    } 
	  dissipation[dit] /= levelCoordsOld.iceDensity();
	  rhs -= dissipation[dit]; 
	  rhs *= -a_dt;

	  //compute heat flux across base due to basal dissipation
	  FArrayBox basalDissipation(rhs.box(),1);
	  m_basalFrictionRelation->computeDissipation
	    (basalDissipation , (*m_velocity[lev])[dit] , (*m_velBasalC[lev])[dit],
	     levelCoordsOld , dit ,rhs.box());
	  

	  //add to user set (e.g geothermal) heat flux
	  basalHeatFlux[dit] += basalDissipation;

	  //zero heat flux outside grounded ice
	  for (BoxIterator bit(rhs.box());bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (levelCoordsOld.getFloatingMask()[dit](iv) != GROUNDEDMASKVAL)
		{
		  basalHeatFlux[dit](iv) = 0.0;
		}
	    }
	  
	  //basalHeatFlux[dit] /= levelCoordsNew.iceDensity()); // scale conversion
	  FArrayBox scaledBasalHeatFlux(basalHeatFlux[dit].box(),basalHeatFlux[dit].nComp());
	  scaledBasalHeatFlux.copy(basalHeatFlux[dit]);
	  scaledBasalHeatFlux /=  levelCoordsNew.iceDensity();

	  //surfaceHeatFlux[dit] /= (levelCoordsNew.iceDensity()); // scale conversion
	  FArrayBox scaledSurfaceHeatFlux(surfaceHeatFlux[dit].box(),surfaceHeatFlux[dit].nComp());
	  scaledSurfaceHeatFlux.copy(surfaceHeatFlux[dit]);
	  scaledSurfaceHeatFlux /= levelCoordsNew.iceDensity();

	  //solve H(t+dt)E(t+dt) + vertical transport terms = H(t)E(t) - rhs(t+/dt)
	  //with either a Dirichlett or flux boundary condition at the upper surface and a flux condition at base
	 
	  Real halftime = time() + 0.5*a_dt;
	  int nLayers = m_nLayers;
	  const Real& rhoi = levelCoordsNew.iceDensity();
	  const Real& rhoo = levelCoordsNew.waterDensity();
	  const Real& gravity = levelCoordsNew.gravity();
	 
	  int surfaceTempDirichlett = surfaceHeatBoundaryDirichlett()?1:0;
	  
	  FORT_UPDATEINTERNALENERGY
	    (CHF_FRA(E), 
	     CHF_FRA1(sT,0), 
	     CHF_FRA1(bT,0),
	     CHF_CONST_FRA1(scaledSurfaceHeatFlux,0),
	     CHF_CONST_FRA1(scaledBasalHeatFlux,0),
	     CHF_CONST_FIA1(levelCoordsOld.getFloatingMask()[dit],0),
	     CHF_CONST_FIA1(levelCoordsNew.getFloatingMask()[dit],0),
	     CHF_CONST_FRA(rhs),
	     CHF_CONST_FRA1(oldH,0),
	     CHF_CONST_FRA1(newH,0),
	     CHF_CONST_FRA((*vectUSigma[lev])[dit]),
	     CHF_CONST_VR(levelCoordsOld.getFaceSigma()),
	     CHF_CONST_VR(dSigma),
	     CHF_CONST_REAL(halftime), 
	     CHF_CONST_REAL(a_dt),
	     CHF_CONST_REAL(rhoi),
	     CHF_CONST_REAL(rhoo),
	     CHF_CONST_REAL(gravity),
	     CHF_CONST_INT(nLayers),
	     CHF_CONST_INT(surfaceTempDirichlett),
	     CHF_BOX(box));

	  scaledBasalHeatFlux *= (levelCoordsNew.iceDensity());
	  basalHeatFlux[dit].copy(scaledBasalHeatFlux);
	  scaledSurfaceHeatFlux *= (levelCoordsNew.iceDensity());
	  surfaceHeatFlux[dit].copy(scaledSurfaceHeatFlux);	    
	} // end update internal energy loop over grids
    } // end update internal energy loop over levels

  //coarse average from finer levels & exchange
  for (int lev = m_finest_level; lev >= 0 ; --lev)
    {
      if (lev > 0)
	{
	  CoarseAverage avN(m_amrGrids[lev],
			    m_amrGrids[lev-1],
			    m_internalEnergy[lev]->nComp(),
			    m_refinement_ratios[lev-1], 
			    IntVect::Zero);
	  
	  
	  
	  avN.averageToCoarse(*m_internalEnergy[lev-1], *m_internalEnergy[lev]);
	
	  
	  CoarseAverage avOne(m_amrGrids[lev],m_amrGrids[lev-1],
			      1,m_refinement_ratios[lev-1], IntVect::Zero);
	  
	  avOne.averageToCoarse(*m_sInternalEnergy[lev-1], *m_sInternalEnergy[lev]);
	  avOne.averageToCoarse(*m_bInternalEnergy[lev-1], *m_bInternalEnergy[lev]);
	  avOne.averageToCoarse(*m_sHeatFlux[lev-1], *m_sHeatFlux[lev]);
	  avOne.averageToCoarse(*m_bHeatFlux[lev-1], *m_bHeatFlux[lev]);
	}
      
      m_internalEnergy[lev]->exchange();
      m_sInternalEnergy[lev]->exchange();
      m_bInternalEnergy[lev]->exchange();
      m_sHeatFlux[lev]->exchange();
      m_bHeatFlux[lev]->exchange();
    }
  
  for (int lev = 0; lev < vectLayerFluxes.size(); ++lev)
    {
      if (vectUSigma[lev] != NULL)
	{
	  delete vectUSigma[lev]; vectUSigma[lev] = NULL;
	}

      if (vectDivUHxy[lev] != NULL)
	{
	  delete  vectDivUHxy[lev]; vectDivUHxy[lev] = NULL;
	}

      if (vectLayerFluxes[lev] != NULL)
	{
	  delete vectLayerFluxes[lev];vectLayerFluxes[lev] = NULL;
	}

      if (vectLayerThicknessFluxes[lev] != NULL)
	{
	  delete vectLayerThicknessFluxes[lev];vectLayerThicknessFluxes[lev] = NULL;
	}

    }

  //finally, A is no longer valid 
  m_A_valid = false;
  //#endif

}
#endif


#if BISICLES_Z == BISICLES_LAYERED
//compute the face- and layer- centered internal energy (a_layerEH_half)
//and thickness (a_layerH_half) at time a_time + 1/2 * a_dt
void AmrIce::computeInternalEnergyHalf(Vector<LevelData<FluxBox>* >& a_layerEH_half,
				       Vector<LevelData<FluxBox>* >& a_layerH_half,
				       const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel, 
				       const Real a_dt, const Real a_time)
{

  CH_TIME("AmrIce::computeInternalEnergyHalf");
  
  //delete and re-create storage for a_layerEH_half and a_layerH_half.
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      
      if (a_layerEH_half[lev] != NULL)
	delete(a_layerEH_half[lev]);

      a_layerEH_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 
						   m_internalEnergy[lev]->nComp(), 
						   IntVect::Unit);
      if (a_layerH_half[lev] != NULL)
	delete(a_layerH_half[lev]);
      
      a_layerH_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 
						  m_internalEnergy[lev]->nComp(), 
						  IntVect::Unit);
    }


  //assume the ghost regions of m_internalEnergy are not correct
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      if (lev > 0)
	{
	  PiecewiseLinearFillPatch pwl(m_amrGrids[lev],
				       m_amrGrids[lev-1],
				       m_internalEnergy[lev]->nComp(),
				       m_amrDomains[lev-1],
				       m_refinement_ratios[lev-1],
				       m_internalEnergy[lev]->ghostVect()[0]);
	  pwl.fillInterp(*m_internalEnergy[lev],*m_internalEnergy[lev-1],
			 *m_internalEnergy[lev-1],1.0,0,0,m_internalEnergy[lev]->nComp());
	}
      m_internalEnergy[lev]->exchange();
    }
   
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      //in the 2D case (with poor man's multidim) this
      //is a little pained using AdvectPhysics, but for the time being
      //we need to construct a single component thisLayerEH_Half for each layer, 
      //given a internalEnergy and horizontal velocity and then copy it into 
      //the multicomponent EH_half[lev]
     
      // PatchGodunov object for layer thickness/energy advection
      PatchGodunov patchGodunov;
      {
	int normalPredOrder = 1;
	bool useFourthOrderSlopes = false;
	bool usePrimLimiting = false;
	bool useCharLimiting = false;
	bool useFlattening = false;
	bool useArtificialViscosity = false;
	Real artificialViscosity = 0.0;
	AdvectPhysics advectPhys;
	advectPhys.setPhysIBC(m_internalEnergyIBCPtr);

	patchGodunov.define(m_amrDomains[lev], m_amrDx[lev],
			    &advectPhys, normalPredOrder,
			    useFourthOrderSlopes,usePrimLimiting,
			    useCharLimiting,useFlattening,
			    useArtificialViscosity,artificialViscosity);
	patchGodunov.setCurrentTime(m_time);
      }
      
      AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>(patchGodunov.getGodunovPhysicsPtr());
      if (advectPhysPtr == NULL)
	{
	  MayDay::Error("AmrIce::computeInternalEnergyHalf -- unable to upcast GodunovPhysics to AdvectPhysics");
	}

      const LevelData<FArrayBox>& levelInternalEnergy = *m_internalEnergy[lev]; 
      const LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev]; 
      const LevelData<FluxBox>& levelLayerXYFaceXYVel = *a_layerXYFaceXYVel[lev]; 
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      for (int layer = 0; layer < m_nLayers; ++layer)
	{
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      const Box& box = levelInternalEnergy[dit].box(); // grid box plus ghost cells
	      
	      FluxBox layerXYFaceXYVel(box,1);
	      layerXYFaceXYVel.setVal(0.0);
	      for (int dir = 0; dir < SpaceDim; ++dir){
		layerXYFaceXYVel[dir].copy(levelLayerXYFaceXYVel[dit][dir],layer,0,1);
		Box faceBox = levelGrids[dit].surroundingNodes(dir);
		CH_assert(layerXYFaceXYVel[dir].norm(faceBox,0) < HUGE_NORM);
	      }

	      FArrayBox layerCellXYVel(box,SpaceDim);
	      EdgeToCell(layerXYFaceXYVel,layerCellXYVel);

	      //\todo compute bulk heat sources
	      FArrayBox heatSource(levelGrids[dit], 1);
	      heatSource.setVal(0.0);

	      patchGodunov.setCurrentBox(levelGrids[dit]);
	      advectPhysPtr->setVelocities(&layerCellXYVel,&layerXYFaceXYVel);

	      FArrayBox WGdnv(box,1);

	      //HE at half time and cell faces
	      WGdnv.copy(levelInternalEnergy[dit],layer,0,1);
	      WGdnv *= levelOldThickness[dit];
	      Box grownBox = levelGrids[dit];
	      grownBox.grow(1);
	      FluxBox HEhalf(grownBox,1);
	      patchGodunov.computeWHalf(HEhalf,
					WGdnv,
					heatSource,
					a_dt,
					levelGrids[dit]);
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		  Box faceBox(levelGrids[dit]);
		  faceBox.surroundingNodes(dir);
		  CH_assert(HEhalf[dir].norm(faceBox,0) < HUGE_NORM);
		  (*a_layerEH_half[lev])[dit][dir].copy(HEhalf[dir],0,layer,1);
		}
	      
	      //H at half time and cell faces
	      WGdnv.copy(levelOldThickness[dit]);
	      FluxBox Hhalf(grownBox,1);
	      //\todo compute layer thickness sources
	      FArrayBox HSource(levelGrids[dit], 1);
	      HSource.setVal(0.0);
	      patchGodunov.computeWHalf(Hhalf,
					WGdnv,
					HSource,
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
      
  // coarse average new EH-Half to covered regions
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],a_layerEH_half[lev]->nComp(), m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*a_layerEH_half[lev-1], *a_layerEH_half[lev]);
      faceAverager.averageToCoarse(*a_layerH_half[lev-1], *a_layerH_half[lev]);
    }

}

#endif

#include "NamespaceFooter.H"
