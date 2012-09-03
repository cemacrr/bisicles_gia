#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _L1L2CONSTITUTIVERELATION_H_
#define _L1L2CONSTITUTIVERELATION_H_

#include "ConstitutiveRelation.H"
#include "IceThicknessIBC.H"
#include "NamespaceHeader.H"


#if BISICLES_Z == BISICLES_LAYERED
/// class which encapsulates the L1L2 constitutive relation for ice
/** 
    Described by Schoof, C and Hindmarsh, R C A, ...
    Like GlensFlowRelation, this class computes \f$\mu_{AS}\f$ such that (stress =
    \f$\mu_{AS}\f$*strainrate ) assuming that \f$\mu_{AS} = \mu(u,\theta)\f$.
    But mu is now found by solving a nonlinear problem 
    \f$1 - 2 \mu A (4 \mu ^2 \epsilon^2 + \tau^2)^((n-1)/2) = 0 \f$
    numerically for mu at each point and integrating vertically.
    Although an analytic solution does exist when n = 3, it is expensive/unreliable to compute
*/
class L1L2ConstitutiveRelation : public ConstitutiveRelation
{

  GlensFlowRelation glensFlowRelation;

public:
  ///
  L1L2ConstitutiveRelation() :m_solverTol(1.0e-4) {;}
  ///
  virtual ~L1L2ConstitutiveRelation() {;}

  GlensFlowRelation* getGlensFlowRelationPtr()
  { return &glensFlowRelation;}
  GlensFlowRelation& getGlensFlowRelation()
  { return glensFlowRelation;}

  virtual Real power() const 
  { return glensFlowRelation.power(); }

  /// creates a new copy of this ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

  /// computes cell-centered \f$\mu_{AS}\f$ based on the cell-centered velocity
  /**
     \param a_mu -- mu_{AS} based on the local velocity field.
     \param a_vel -- Cell-centered velocity field.
     \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     \param a_nRefCrse -- refinement ratio to next coarser level     
     \param a_A: Cell-centered flow law coefficient (Glenn's A) field
     \param a_coordSys: LevelSigmaCS object containing the geometry.
     \param a_ghostVect: how the boxes on which we want to compute mu compare
                         with the DisjointBoxLayout. (can be negative)
  */
  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel, 
                         const LevelData<FArrayBox>* a_crseVelPtr,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const;
  
 
  /// Compute a cell centred bulk dissipation \f$\Phi = 4\mu\epsilon^2\f$ 
  /// (heat source) at the cell centres. This ought to have the same dependencies 
  /// as cell centred mu
  /**
     a_dissipation -- \f$\Phi\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
     (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (Glenn's A) field
     a_coordSys:  SigmaCS object containing the geometry of this patch.
     a_box: cell-centered box over which to do this computation
  **/
  virtual void computeDissipation(LevelData<FArrayBox>& a_dissipation,
				  const LevelData<FArrayBox>& a_vel, 
				  const LevelData<FArrayBox>* a_crseVel,
				  int a_nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const;
  
  /// computes face-centered \f$\mu_{AS}\f$ based on cell-centered velocity
  /** \param a_mu: face-centered mu_{AS} based on the local velocity field.
      \param a_vel: Cell-centered velocity field.
      \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                          (NULL if no coarser level)
      \param a_nRefCrse -- refinement ratio to next coarser level     
      \param a_A: Cell-centered flow law coefficient (Glenn's A) field
      \param a_coordSys: LevelSigmaCS object containing the geometry.
      \param a_ghostVect: how the boxes on which we want to compute mu compare
                          with the DisjointBoxLayout. (can be negative)
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, 
                             const LevelData<FArrayBox>* a_crseVelPtr,
                             int nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;

  
  ///override method from ConstitutiveRelation
  virtual void computeFaceFluxVelocity(const LevelData<FArrayBox>& a_cellVel,
                                       const LevelData<FArrayBox>* a_crseVelPtr,
                                       int a_nRefCrse,
				       const LevelSigmaCS& a_coordSys,
				       const DisjointBoxLayout& a_grids,
				       const ProblemDomain& a_domain,
				       const LevelData<FArrayBox>& a_A,
				       const LevelData<FArrayBox>& a_surfaceA,
				       const LevelData<FArrayBox>& a_baseA,
				       LevelData<FluxBox>& a_faceVel,
				       LevelData<FluxBox>& a_layerXYFaceXYVel,
				       LevelData<FArrayBox>& a_layerSFaceXYVel
				       ) const;


  virtual void modifyTransportCoefficients
  (const LevelData<FArrayBox>& a_cellVel,
   const LevelData<FArrayBox>* a_crseVelPtr,
   int a_nRefCrse,
   const LevelSigmaCS& a_coordSys,
   const DisjointBoxLayout& a_grids,
   const ProblemDomain& a_domain,
   const LevelData<FArrayBox>& a_A,
   const LevelData<FArrayBox>& a_sA,
   const LevelData<FArrayBox>& a_bA,
   LevelData<FluxBox>& a_faceVelAdvection,
   LevelData<FluxBox>& a_faceVelTotal,
   LevelData<FluxBox>& a_faceDiffusivity,
   LevelData<FluxBox>& a_layerXYFaceXYVel,
   LevelData<FArrayBox>& a_layerSFaceXYVel) const;
    
  virtual void modifyThicknessDiffusivity(const LevelData<FArrayBox>& a_cellVel,
                                          const LevelData<FArrayBox>* a_crseVelPtr,
                                          int a_nRefCrse,
					  const LevelSigmaCS& a_coordSys,
					  const DisjointBoxLayout& a_grids,
					  const LevelData<FluxBox>& a_faceA,
					  LevelData<FluxBox>& a_diffusivity) const;

  /// sets tolerance for secant solve
  void solverTolerance(Real a_solverTol) 
  { m_solverTol = a_solverTol; }


private:
  ///compute mu for given sigma at a_box cell centres,
  /// plus velocity gradient & rate-of-strain invariant
  /// at base
  ///a_mu has one component for each layer
  virtual void computeMuZ(LevelData<FArrayBox>& a_mu,
			  LevelData<FArrayBox>& a_gradVel,
			  LevelData<FArrayBox>& a_epsSqr,
			  const Vector<Real>& a_sigma,
			  const LevelData<FArrayBox>& a_vel, 
                          const LevelData<FArrayBox>* a_crseVelPtr,
                          int a_nRefCrse,
			  const LevelData<FArrayBox>& a_A,
			  const LevelSigmaCS& a_coordSys,
			  const ProblemDomain& a_domain,
			  const IntVect& a_ghostVect = IntVect::Unit) const;
  
  

  //compute mu for given sigma at a_faceDir cell faces
  /// plus velocity gradient & rate-of-strain invariant
  /// at base
  ///a_mu has one component for each layer
  virtual void computeFaceMuZ(LevelData<FluxBox>& a_mu,
			      LevelData<FluxBox>& a_gradVel,
			      LevelData<FluxBox>& a_epsSqr,
			      const Vector<Real>& a_sigma,
			      LevelData<FArrayBox>& a_vel, 
                              const LevelData<FArrayBox>* a_crseVelPtr,
                              int a_nRefCrse,
			      const LevelData<FluxBox>& a_A,
			      const LevelSigmaCS& a_coordSys,
			      const ProblemDomain& a_domain,
			      const IntVect& a_ghostVect) const;
  

  /// calculate mu for given sigma at a cell faces or centers
  /// plus velocity gradient & rate-of-strain invariant
  /// at base
  ///a_mu has one component for each layer
  virtual void computeEitherMuZ(FArrayBox& a_mu,
				const Vector<Real>& a_sigma,
				const FArrayBox& a_grads,
				const FArrayBox& a_epsSqr,
				const FArrayBox& a_A,
				const FArrayBox& a_H,
				const Real& a_rhog,
				const Box& a_box) const;

  /// computes the vertically averaged, face-centred, face direction
  /// component of the velocity field
  /** 
      \param a_vel: the velocity field at the ice base at cell centers
      \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                          (NULL if no coarser level)
      \param a_nRefCrse -- refinement ratio to next coarser level     \
      \param a_faceVel: the vertically averaged velocity field at cell faces
      \param a_faceDir : face direction
      \param a_cellBox : box of cell centres: computes face flux for each of these
      The L1L2 constutive relation implies a z-varying
      velocity field, which is constructed from a 2D velocity
      field computed at the ice base. 
  */  
  

  virtual void computeFaceFluxVelocity(LevelData<FArrayBox>& a_vel,
                                       const LevelData<FArrayBox>* a_crseVelPtr,
                                       int a_nRefCrse,
                                       const LevelData<FArrayBox>& a_A,
                                       const LevelData<FArrayBox>& a_thickness,
                                       const RealVect& a_dx,
                                       LevelData<FluxBox>& a_fluxVel,
				       LevelData<FluxBox>& a_layerXYFaceXYVel,
				       LevelData<FArrayBox>& a_layerSFaceXYVel,
                                       const LevelSigmaCS& a_coordSys,
				       const ProblemDomain& a_domain,
                                       const IntVect& a_cellGhost = IntVect::Zero) const;

  

  /// thickness diffusivity at one face 
  virtual void modifyThicknessDiffusivity(const FArrayBox& a_vel,
                                          const FArrayBox& a_A,
                                          const FArrayBox& a_thickness,
                                          const FArrayBox& a_bedrock,
                                          const FArrayBox& a_faceThickness,
                                          const FArrayBox& a_faceGradSurface,
                                          const BaseFab<int>& a_floatingMask,
                                          const bool& a_anyFloating,
                                          const Real a_seaLevel,
					  const Real a_iceDensity,
					  const Real a_seaWaterDensity,
					  const Real a_gravity,
					  const Vector<Real>& a_sigma,
					  const Vector<Real>& a_dSigma,
                                          const Box& a_cellBox,
                                          FArrayBox& a_diffusivity,
                                          int a_faceDir) const;

  Real m_solverTol;



};
#endif 
#include "NamespaceFooter.H"
#endif
