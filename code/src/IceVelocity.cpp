
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IceVelocity.H"
#include "IceConstants.H"
#include "amrIceF_F.H"
#include "AMRPoissonOpF_F.H"
#include "QuadCFInterp.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "PiecewiseLinearFillPatch.H"
#include "L1L2ConstitutiveRelation.H"
#include "NamespaceHeader.H"

/// compute RHS for velocity field solve
void IceVelocity::defineRHS(Vector<LevelData<FArrayBox>* >& a_rhs,
			    const Vector<RefCountedPtr<LevelSigmaCS > >& a_CS,
			    const Vector<DisjointBoxLayout>& a_grids,
			    const Vector<RealVect>& a_dx)
			    
{
  CH_assert(SpaceDim == 2);
  int xDir = SpaceDim-2;
  int yDir = SpaceDim-1;
  CH_assert(a_rhs.size() <= a_CS.size());
  for (int lev = 0; lev < a_rhs.size(); ++lev)
    {
      CH_assert(a_rhs[lev] != NULL && a_CS[lev] != NULL);
      LevelData<FArrayBox>& levelRhs = *a_rhs[lev];
      const LevelSigmaCS& levelCS = *a_CS[lev];
      const DisjointBoxLayout& levelGrids = a_grids[lev];
      Real rhog = levelCS.iceDensity()*levelCS.gravity();
      const RealVect& levelDx = a_dx[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  FArrayBox& rhs = levelRhs[dit];
	  const FArrayBox& gradS = levelCS.getGradSurface()[dit];
	  const FArrayBox& thck = levelCS.getH()[dit];
	  
	  rhs.copy(gradS);
	  rhs.mult(thck,0,xDir);
	  rhs.mult(thck,0,yDir);
	  rhs *= rhog;

	  const bool& anyFloating = levelCS.anyFloating()[dit];
	  
	  if (anyFloating)
	    {
	      const Box& box = levelGrids[dit];
	      const FArrayBox& usrf = levelCS.getSurfaceHeight()[dit];
	      const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	      
	      for (int dir=0; dir<SpaceDim; dir++)
		{
		  FORT_GLCORRECTION(CHF_FRA1(rhs, dir),
				    CHF_CONST_FRA1(thck,0),
				    CHF_CONST_FRA1(usrf,0),
				    CHF_CONST_FIA1(mask,0),
				    CHF_INT(dir),
				    CHF_CONST_REAL(levelDx[dir]),
				    CHF_CONST_REAL(rhog),
				    CHF_BOX(box));
		}
	    }
	}
    }
}

// compute RHS for velocity field solve
void IceVelocity::setFloatingC(Vector<LevelData<FArrayBox>* >& a_C,
			       const Vector<RefCountedPtr<LevelSigmaCS > >& a_CS,
			       const Vector<DisjointBoxLayout>& a_grids,
			       const Real& a_basalFrictionDecay)
{
   CH_assert(SpaceDim == 2);
   CH_assert(a_C.size() <= a_CS.size());
  for (int lev = 0; lev < a_C.size(); ++lev)
    {
      CH_assert(a_C[lev] != NULL && a_CS[lev] != NULL);
      LevelData<FArrayBox>& levelC = *a_C[lev];
      const LevelSigmaCS& levelCS = *a_CS[lev];
      const DisjointBoxLayout& levelGrids = a_grids[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	   const bool& anyFloating = levelCS.anyFloating()[dit];
	   if (anyFloating)
	     {
	       const Box& box = levelGrids[dit];
	       FArrayBox& C = levelC[dit];
	       const FArrayBox& usrf = levelCS.getSurfaceHeight()[dit];
	       const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	       const FArrayBox& topg = levelCS.getBaseHeight()[dit];

	       FArrayBox lsrf(box,1);
	       lsrf.copy(usrf);
	       lsrf -= levelCS.getH()[dit];
	       
	       FORT_SETFLOATINGBETA(CHF_FRA1(C,0),
				    CHF_CONST_FIA1(mask,0),
				    CHF_CONST_FRA1(lsrf,0),
				    CHF_CONST_FRA1(topg,0),
				    CHF_CONST_REAL(a_basalFrictionDecay),
				    CHF_BOX(box));
	     }
	}
    }
}

//apply cell-centred helmholtz operator a phi + b grad^2(phi) 
//to cell-centred phi. Assumes that boundary values have been
//set.
void IceVelocity::applyHelmOp
(LevelData<FArrayBox>& a_lapPhi,
 const LevelData<FArrayBox>& a_phi, 
 const Real& a_a, const Real& a_b,
 const DisjointBoxLayout& a_grids,
 const RealVect& a_dx )
{
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      FORT_OPERATORLAP(CHF_FRA(a_lapPhi[dit]),
		       CHF_FRA(a_phi[dit]),
		       CHF_BOX(a_grids[dit]),
		       CHF_CONST_REAL(a_dx[0]),
		       CHF_CONST_REAL(a_a),
		       CHF_CONST_REAL(a_b));
    }
}



//apply cell-centred  operator div(u)
//to face-centred u. Assumes that ghost cells
//have been set
void IceVelocity::applyDiv
(LevelData<FArrayBox>& a_divU,
 const LevelData<FluxBox>& a_u, 
 const DisjointBoxLayout& a_grids,
 const RealVect& a_dx
)
{
 
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      a_divU[dit].setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  FORT_DIVERGENCE(CHF_CONST_FRA(a_u[dit][dir]),
			  CHF_FRA(a_divU[dit]),
			  CHF_BOX(a_grids[dit]),
			  CHF_CONST_REAL(a_dx[dir]),
			  CHF_INT(dir));
	}
    }
}




//compute face-centered us given cell-centered
//vector u and scalar s 
//Assumes that boundary values have been set.
//USES CENTERED DIFFERENCES
void IceVelocity::computeFaceFlux
(LevelData<FluxBox>& a_us,
 const LevelData<FArrayBox>& a_u, 
 const LevelData<FArrayBox>& a_s, 
 const DisjointBoxLayout& a_grids)
{
  CH_assert(a_s.nComp() == 1);
  CH_assert(a_us.nComp() == 1);
  CH_assert(a_u.nComp() == SpaceDim);
  LevelData<FArrayBox> ccus(a_grids,SpaceDim,IntVect::Unit);
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      ccus[dit].copy(a_u[dit]);
      for (int dir = 0; dir < SpaceDim; ++dir) 
	ccus[dit].mult(a_s[dit],0,dir,1);
    }
  CellToEdge(ccus,a_us);
}


//apply cell-centred gradient operator grad (phi) 
//to cell-centred phi.Assumes that ghost cells
//have been set
void IceVelocity::applyGrad
(LevelData<FArrayBox>& a_gradPhi,
 const LevelData<FArrayBox>& a_phi, 
 const DisjointBoxLayout& a_grids,
 const RealVect& a_dx)
{
	
  for (DataIterator dit(a_grids); dit.ok(); ++dit)
    {
      for (int icomp = 0; icomp < a_phi.nComp(); icomp++) 
	{
	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      Real oneOnTwoDx = 1.0 / (2.0 * a_dx[dir]);
	      for (BoxIterator bit(a_grids[dit]);bit.ok();++bit)
		{
		  
		  const IntVect& iv = bit();
		  a_gradPhi[dit](iv,icomp) = 
		    oneOnTwoDx * 
		    (a_phi[dit](iv + BASISV(dir),icomp)
		     -a_phi[dit](iv - BASISV(dir),icomp));
		}
	    }
	}
    }
}

//compute A(x,y,sigma) given temperature, geometry etc
void IceVelocity::computeA
(LevelData<FArrayBox>& a_A,
 const Vector<Real>& a_sigma,
 const LevelSigmaCS& a_coordSys,
 const RateFactor* a_rateFactor,
 const LevelData<FArrayBox>& a_temperature)
{

  const DisjointBoxLayout grids =  a_coordSys.grids();

  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      // compute A(T)
      // need a temperature field corrected to the pressure melting point,
      // \theta ^* = \min(\theta,\theta _r) + a * p)
      // a is a constant, p is pressure, \thera _r is the melting point in standard pressure
      // using p = \rho * g * \sigma * H 
      // (used by Glimmer, even with higher order stresses)
      // should be p = T_xx + T_yy + \rho * g * \sigma * H
      Real Tmax = triplepoint - TINY_NORM;
      Real fbase = a_coordSys.iceDensity() * a_coordSys.gravity() * icepmeltfactor;
      for (int layer = 0; layer < a_sigma.size(); ++layer)
	{
	  Real f = fbase * a_sigma[layer];
	  FArrayBox layerA;
	  layerA.define(Interval(layer,layer),a_A[dit]);
	  const Box& box = layerA.box();
	  FArrayBox thetaStar(box,1);
	  
	  FORT_FABMINPLUS(CHF_FRA1(thetaStar,0),
			  CHF_FRA1(a_temperature[dit],layer),
			  CHF_FRA1(a_coordSys.getH()[dit],0),
			  CHF_CONST_REAL(f),
			  CHF_CONST_REAL(Tmax),
			  CHF_BOX(box));
	
	  CH_assert(0.0 < thetaStar.min(box));
	  CH_assert(thetaStar.max(box) < triplepoint); 
	  a_rateFactor->computeA(layerA,thetaStar,box);
	}
    }

}

void IceVelocity::addWallDrag(FArrayBox& a_drag, 
			      const BaseFab<int>& a_mask,
			      const FArrayBox& a_usrf,
			      const FArrayBox& a_thk,
			      const FArrayBox& a_topg,
			      const FArrayBox& a_beta,
			      const Real& a_extra, 
			      const RealVect& a_dx,
			      const Box& a_box)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      for (int sign = -1; sign <= 1; sign+=2)
	{

	  for (BoxIterator bit(a_box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if (a_mask(iv) == GROUNDEDMASKVAL || a_mask(iv) == FLOATINGMASKVAL)
		{
		  const IntVect ivp = iv + sign*BASISV(dir);
		  if (a_mask(ivp) != GROUNDEDMASKVAL && a_mask(ivp) != FLOATINGMASKVAL)
		    {
		      Real contact = 
			std::min(a_thk(iv) * 0.5,  a_topg(ivp) - (a_usrf(iv)-a_thk(iv)));
		      if (contact > 0.0)
			{
			  a_drag(iv) += (a_beta(iv) + a_extra) * contact / a_dx[dir];
			}
		    }
		}
	    }
	}
    }
}

void IceVelocity::computeFaceVelocity
(LevelData<FluxBox>& a_faceVelAdvection,
 LevelData<FluxBox>& a_faceVelTotal,
 LevelData<FluxBox>& a_faceDiffusivity,
 LevelData<FArrayBox>& a_cellDiffusivity,
#if BISCICLES_Z == BISICLES_LAYERED
 LevelData<FluxBox>& a_layerXYFaceXYVel,
 LevelData<FArrayBox>& a_layerSFaceXYVel,
#endif
 const LevelData<FArrayBox>& a_velocity,
 const LevelSigmaCS& a_coordSys,
 const LevelData<FArrayBox>& a_A,
#if BISCICLES_Z == BISICLES_LAYERED
 const LevelData<FArrayBox>& a_sA,
 const LevelData<FArrayBox>& a_bA,
#endif			 
 const LevelData<FArrayBox>* a_crseVelocity,
 const LevelData<FArrayBox>* a_crseDiffusivity,
 const int a_nRefCrse,
 const ConstitutiveRelation* a_constitutiveRelation,
 const bool a_additionalVelocity) 
{			   
  
  //We need to copy the cell-centered velocity onto
  //LevelData<FArrayBox> with 2 ghost cells, because we will 
  //need one ghost cell's worth of face-centred data.
  IntVect grownVelGhost(2*IntVect::Unit);
  const DisjointBoxLayout& grids = a_velocity.disjointBoxLayout();
  LevelData<FArrayBox> grownVel(grids,a_velocity.nComp(), grownVelGhost);
  
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      grownVel[dit].setVal(0.0);
      grownVel[dit].copy(a_velocity[dit], a_velocity[dit].box());
      a_cellDiffusivity[dit].setVal(0.0);
      for (int dir = 0; dir < SpaceDim; dir++)
	a_faceDiffusivity[dit][dir].setVal(0.0);
    }
  
  if (a_crseVelocity != NULL)
    {
      const DisjointBoxLayout& crseGrids = a_crseVelocity->disjointBoxLayout();
      PiecewiseLinearFillPatch velFiller(grids , crseGrids, a_velocity.nComp(), 
					 crseGrids.physDomain(), a_nRefCrse,
					 grownVelGhost[0]);
	
      Real time_interp_coeff = 0.0;
      velFiller.fillInterp(grownVel, *a_crseVelocity, *a_crseVelocity,
			   time_interp_coeff,0, 0, a_velocity.nComp());
	
      //slc. In L1L2 we will compute second derivatives of v, hence
      //we need quadratic coarse-fine interpolation here.
      Real dx = 1.0;
      QuadCFInterp qcfi(grids , &crseGrids, dx, a_nRefCrse, 
			a_velocity.nComp(), grids.physDomain());
      qcfi.coarseFineInterp(grownVel, *a_crseVelocity);

    }

  grownVel.exchange();
      
  //default calculation : average to faces 
  CellToEdge(grownVel, a_faceVelAdvection);
#if BISCICLES_Z == BISICLES_LAYERED
  //for layered models (SSA,L1L2) assume du/dz = 0
  for (int j = 0; j < a_layerXYFaceXYVel.nComp(); ++j)
    a_faceVelAdvection.copyTo(Interval(0,0), a_layerXYFaceXYVel, Interval(j,j));
  
  for (int j = 0; j < a_layerSFaceXYVel.nComp(); j+=SpaceDim)
    grownVel.copyTo(Interval(0,SpaceDim-1), a_layerSFaceXYVel,Interval(j,j+SpaceDim-1)); 
#endif

  //modification to fluxes at the margins, that is where mask changes to open sea or land.
  //
  // -----|-----|-----|-----|-----|
  //      |     |     |     |     |
  //   x  o  x  o  x  F     |     |
  //      |     |     |     |     |
  // -----|-----|-----|-----|-----|
  //
  //  n-2   n-2    n    n+1
  //
  // 2D case (L1L2 / SSA)
  // there are valid values of the basal velocity at points x and
  // and the z-varying velocity at points o. The points at o 
  // have been interpolated from the x, and should vary little vertically
  // We need the flux at F, but since there is no x_n+1, the interpolated
  // value makes no sense. Using the margin boundary condition to get
  // du/dx at the face is perhaps the ideal approach , but for now we just
  // take x_{n-1} and o_{n-1} and extrapolate. 
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  Box faceBox = grids[dit];
	  faceBox.surroundingNodes(dir);
	    
	  FArrayBox& faceVel = a_faceVelAdvection[dit][dir];
	  CH_assert(faceVel.box().contains(faceBox));
	  const FArrayBox& cellVel = grownVel[dit];
	  const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
	  //CH_assert(faceVel.norm(faceBox,0) < HUGE_NORM);
	  FORT_EXTRAPTOMARGIN(CHF_FRA1(faceVel,0),
			      CHF_CONST_FRA1(cellVel,dir),
			      CHF_CONST_FIA1(mask,0),
			      CHF_CONST_INT(dir),
			      CHF_BOX(faceBox));
	  //pout() << "FORT_EXTRAPTOMARGIN" << std::endl;
	  CH_assert(faceVel.norm(faceBox,0) < HUGE_NORM);

	  FArrayBox& faceVelTotal = a_faceVelTotal[dit][dir];
	  faceVelTotal.copy(faceVel);
	}
    }
    
  
  


#if BISCICLES_Z == BISICLES_LAYERED
  //copy the basal velocity into the vertically varying velocities.
  {
    
    for (DataIterator dit(grids); dit.ok(); ++dit)
      {

	

	FArrayBox& SFaceXYVel = a_layerSFaceXYVel[dit];
	const FArrayBox& cellVel = grownVel[dit];
	for (int ic = 0; ic < SFaceXYVel.nComp()-1; ic+=SpaceDim)
	  {
	    SFaceXYVel.copy( cellVel , 0, ic, SpaceDim);
	  }

	for (int dir = 0; dir < SpaceDim; ++dir)
	  {
	    const FArrayBox& faceVel = a_faceVelAdvection[dit][dir];
	    FArrayBox& XYFaceXYVel = a_layerXYFaceXYVel[dit][dir];
	    
	    for (int ic = 0; ic < XYFaceXYVel.nComp(); ic++)
	      {
		XYFaceXYVel.copy( faceVel , 0, ic, 1);
	      }
	  }
      }
  }
#endif 


  if (a_additionalVelocity)
    {
      const L1L2ConstitutiveRelation* L1L2Ptr = dynamic_cast<const L1L2ConstitutiveRelation*>(a_constitutiveRelation);
      if (L1L2Ptr != NULL)
	{
	  //L1L2Ptr->computeFaceFluxVelocity(grownVel, a_crseVelocity, a_nRefCrse, a_coordSys, 
	  //				   grids,  grids.physDomain(), a_A, a_sA, a_bA,
	  //				   a_faceVelTotal ,a_layerXYFaceXYVel, a_layerSFaceXYVel);
	 

	  
	  L1L2Ptr->modifyTransportCoefficients(grownVel, a_crseVelocity, a_crseDiffusivity, a_nRefCrse, a_coordSys, 
					       grids,  grids.physDomain(), a_A, a_sA, a_bA,
					       a_faceVelAdvection, a_faceVelTotal, a_faceDiffusivity,
					       a_cellDiffusivity, a_layerXYFaceXYVel, a_layerSFaceXYVel);

	}
    }
  
  a_faceVelAdvection.exchange();
  a_faceVelTotal.exchange();
#if BISCICLES_Z == BISICLES_LAYERED
  a_layerXYFaceXYVel.exchange();
  a_layerSFaceXYVel.exchange();
#endif
}



#include "NamespaceFooter.H"
