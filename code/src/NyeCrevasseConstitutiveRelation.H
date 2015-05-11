#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _NYECREVASSECONSTITUTIVERELATION_H_
#define _NYECREVASSECONSTITUTIVERELATION_H_

#include "ConstitutiveRelation.H"
#include "NamespaceHeader.H"
/// Computes the effective viscosity \f$\mu\f$ including a factor due to crevassing
/**
   The effective viscosity \f$\mu\f$ is computed by multiplying
   an effective viscosity \f$\mu_0\f$ as calculated by another subclass
   of ConstitutiveRelation by a factor \f$(1-d)/H\f$, where \f$d\f$
   is the total depth of crevasses (surface plus base). \f$d\f$ is
   set to \f$d = a d_{\rm NYE} + b d_{\rm OLD}\f$. \f$d_{\rm NYE}\f$
   is the Nye crevasse depth, and \f$d_{\rm OLD}\f$ is an 'old' crevasse
   depth, either the Nye depth ffrom a previous time step, os some other
   addition due to e.g water filling.
 */
class NyeCrevasseConstitutiveRelation : public ConstitutiveRelation
{
  ConstitutiveRelation* m_uncrevassedConstitutiveRelation;
public:

  NyeCrevasseConstitutiveRelation(ConstitutiveRelation* a_ptr,
  				  const Real& a_NyeWaterDepth, 
  				  const Real& a_NyeA)
  : m_NyeWaterDepth(a_NyeWaterDepth),
    m_NyeA(a_NyeA)
  {
    m_uncrevassedConstitutiveRelation = a_ptr->getNewConstitutiveRelation();
  }


  virtual ~NyeCrevasseConstitutiveRelation();
  
  /// computes cell-centered \f$\mu_{AS}\f$ based on the cell-centered velocity
  /**
     a_mu -- \f$\mu_{AS}\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
     (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (glenn's A) field
     a_coordSys:  SigmaCS object containing the geometry of this patch.
     a_ghostVect -- how the boxes over which we want to compute mu relate 
                    to those in the DisjointBoxLayout (can be negative)
  */
  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel, 
                         const LevelData<FArrayBox>* a_crseVel,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const;
  
  /// Compute a cell centred bulk dissipation 
  /// \f$\Phi/(\rho _i c _i) = \sigma_ij \epsilon _ji /(\rho _i c _i) \f$ 
  /// (heat source) at the cell centres. This ought to have the same dependencies 
  /// as cell centred mu
  /**
     a_dissipation -- \f$\Phi\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (glenn's A) field
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

  // computes face-centered \f$\mu_{AS}\f$ based on cell-centered velocity
  /** a_mu: face-centered \f$\mu_{AS}\f$ based on the local velocity field.
      a_vel: Cell-centered velocity field.
      a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
      a_nRefCrse -- refinement ratio to next coarser level
      a_A: Cell-centered flow law coefficient (glenn's A) field
      a_coordSys: SigmaCS object containing the geometry of this patch.
      a_ghostVect -- how the boxes over which we want to compute mu relate 
                     to those in the DisjointBoxLayout (can be negative) 
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, 
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;

  /// creates a new copy of this ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

private:

  Real m_NyeWaterDepth;
  Real m_NyeA;
  NyeCrevasseConstitutiveRelation() : m_uncrevassedConstitutiveRelation(NULL)
				      
  {;}

};

#include "NamespaceFooter.H"
#endif
