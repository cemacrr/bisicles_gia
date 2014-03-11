#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _ICEUTILITY_H_
#define _ICEUTILITY_H_
#include "Vector.H"
#include "LevelSigmaCS.H"
#include "ConstitutiveRelation.H"
#include "IceThicknessIBC.H"
#include "NamespaceHeader.H"

/// IceUtility namespace : general purpose ice sheet functions common to the 'big' classes (AmrIce and AMRIceControl for now)
namespace IceUtility
{

  void addWallDrag(FArrayBox& a_drag, 
		   const BaseFab<int>& a_mask,
		   const FArrayBox& a_usrf,
		   const FArrayBox& a_thk,
		   const FArrayBox& a_topg,
		   const FArrayBox& a_beta,
		   const Real& a_extra, 
		   const RealVect& a_dx,
		   const Box& a_box);

  void defineRHS(Vector<LevelData<FArrayBox>* >& a_rhs,
		 const Vector<RefCountedPtr<LevelSigmaCS > >& a_CS,
		 const Vector<DisjointBoxLayout>& a_grids,
		 const Vector<RealVect>& a_dx);


  void  setFloatingC(Vector<LevelData<FArrayBox>* >& a_C,
		    const Vector<RefCountedPtr<LevelSigmaCS > >& a_CS,
		    const Vector<DisjointBoxLayout>& a_grids,
		    const Real& a_basalFrictionDecay);

  /// compute the operator L phi = a * phi + b * grad^2 (phi) 
  void applyHelmOp
  (LevelData<FArrayBox>& a_lapPhi,
   const LevelData<FArrayBox>& a_phi, 
   const Real& a_a, const Real& a_b,
   const DisjointBoxLayout& a_grids,
   const RealVect& a_dx);

  ///apply cell-centred  operator div(phi) to face-centred u. 
  ///Assumes that ghost cells have been set
  void applyDiv
  (LevelData<FArrayBox>& a_divU,
   const LevelData<FluxBox>& a_u, 
   const DisjointBoxLayout& a_grids,
   const RealVect& a_dx);

  ///apply cell-centred  operator [grad(phi)]^2 to face-centred u. 
  ///Assumes that ghost cells have been set
  void applyGradSq
  (LevelData<FArrayBox>& a_gradPhiSq,
   const LevelData<FArrayBox>& a_phi, 
   const DisjointBoxLayout& a_grids, 
   const RealVect& a_dx);
 
  ///compute face-centered us given cell-centered
  ///vector u and scalar s 
  ///Assumes that boundary values have been set.
  ///USES CENTERED SCHEME
  void computeFaceFluxCentered
  (LevelData<FluxBox>& a_us,
   const LevelData<FArrayBox>& a_u, 
   const LevelData<FArrayBox>& a_s, 
   const DisjointBoxLayout& a_grids);

  ///compute face-centered us given cell-centered
  ///vector u and scalar s 
  ///Assumes that boundary values have been set.
  ///USES  FIRST ORDER UPWIND SCHEME
  void computeFaceFluxUpwind
  (LevelData<FluxBox>& a_us,
   const LevelData<FluxBox>& a_u, 
   const LevelData<FArrayBox>& a_s, 
   const DisjointBoxLayout& a_grids);

  ///compute cell centered rate factor A from the temperature 
  void computeA
  (LevelData<FArrayBox>& a_A,
   const Vector<Real>& a_sigma,
   const LevelSigmaCS& a_coordSys,
   const RateFactor* a_rateFactor,
   const LevelData<FArrayBox>& a_temperature);

  ///compute face-centered velocity and thickness diffusion coefficients from cell-centered velocity
  void computeFaceVelocity
  (LevelData<FluxBox>& a_faceVelAdvection,
   LevelData<FluxBox>& a_faceVelTotal,
   LevelData<FluxBox>& a_faceDiffusivity,
   LevelData<FArrayBox>& a_cellDiffusivity,
#if BISICLES_Z == BISICLES_LAYERED
   LevelData<FluxBox>& a_layerXYFaceXYVel,
   LevelData<FArrayBox>& a_layerSFaceXYVel,
#endif
   const LevelData<FArrayBox>& a_velocity,
   const LevelSigmaCS& a_coordSys,
   const IceThicknessIBC* a_bc,
   const LevelData<FArrayBox>& a_A,
#if BISICLES_Z == BISICLES_LAYERED
   const LevelData<FArrayBox>& a_sA,
   const LevelData<FArrayBox>& a_bA,
#endif			 
   const LevelData<FArrayBox>* a_crseVelocity,
   const LevelData<FArrayBox>* a_crseDiffusivity,
   const int a_nRefCrse,
   const ConstitutiveRelation* a_constitutiveRelation,
   const bool a_additionalVelocity); 

  ///Identify regions of floating ice that are remote
  ///from grounded ice and eliminate them.
  void eliminateRemoteIce
  (Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const Vector<DisjointBoxLayout>& a_grids,
   const Vector<ProblemDomain>& a_domain,
   const Vector<int>& a_refRatio, Real a_crseDx,
   int a_finestLevel, int a_maxIter);
  
}

#include "NamespaceFooter.H"
#endif