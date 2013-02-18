#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _ICEVELOCITY_H_
#define _ICEVELOCITY_H_
#include "Vector.H"
#include "LevelSigmaCS.H"
#include "ConstitutiveRelation.H"
#include "NamespaceHeader.H"

namespace IceVelocity
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

  void applyHelmOp
  (LevelData<FArrayBox>& a_lapPhi,
   const LevelData<FArrayBox>& a_phi, 
   const Real& a_a, const Real& a_b,
   const DisjointBoxLayout& a_grids,
   const RealVect& a_dx);

  //apply cell-centred  operator div(u)
  //to face-centred u. Assumes that ghost cells
  //have been set
  void applyDiv
  (LevelData<FArrayBox>& a_divU,
   const LevelData<FluxBox>& a_u, 
   const DisjointBoxLayout& a_grids,
   const RealVect& a_dx);

  void applyGradSq
  (LevelData<FArrayBox>& a_gradPhiSq,
   const LevelData<FArrayBox>& a_phi, 
   const DisjointBoxLayout& a_grids, 
   const RealVect& a_dx);
 
//compute face-centered us given cell-centered
//vector u and scalar s 
//Assumes that boundary values have been set.
//USES CENTERED DIFFERENCES
  void computeFaceFlux
  (LevelData<FluxBox>& a_us,
   const LevelData<FArrayBox>& a_u, 
   const LevelData<FArrayBox>& a_s, 
   const DisjointBoxLayout& a_grids);


  void computeA
  (LevelData<FArrayBox>& a_A,
   const Vector<Real>& a_sigma,
   const LevelSigmaCS& a_coordSys,
   const RateFactor* a_rateFactor,
   const LevelData<FArrayBox>& a_temperature);


  void computeFaceVelocity
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
   const bool a_additionalVelocity); 

}

#include "NamespaceFooter.H"
#endif
