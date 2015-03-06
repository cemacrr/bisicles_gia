#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CONSTITUTIVERELATION_H_
#define _CONSTITUTIVERELATION_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelSigmaCS.H"
#include "CellToEdge.H"
#include "NamespaceHeader.H"



/// class which encapsulates the constitutive relation for ice
/** This class is computes \f$\mu_{AS}\f$ such that (stress =
    \f$\mu_{AS}\f$*strainrate ) assuming that \f$\mu_{AS}\f$ = mu(vel,theta)
*/
class 
ConstitutiveRelation
{

public:

  static ConstitutiveRelation* parse(const char* a_prefix);

  /// default constructor
  ConstitutiveRelation() {;}

  /// 
  virtual ~ConstitutiveRelation() {;}
  
  /// computes cell-centered \f$\mu_{AS}\f$ based on the cell-centered velocity
  /** 
     \param a_mu -- \f$\mu_{AS}\f$ based on the local velocity field. 
     \param a_vel -- Cell-centered velocity field.
     \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     \param a_nRefCrse -- refinement ratio to next coarser level
     \param a_A: Cell-centered flow law coefficient (glenn's A) field
     \param a_coordSys: LevelSigmaCS object containing the ice-sheet geometry
     \param a_ghostVect: how the boxes on which we want to compute mu compare
                         with the DisjointBoxLayout. (can be negative)
  */
  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel, 
                         const LevelData<FArrayBox>* a_crseVel,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const = 0;

  
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
				  const LevelData<FArrayBox>& a_cellVel,
				  const LevelData<FArrayBox>* a_crseVel,
				  int nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const = 0; 

  // computes face-centered \f$\mu_{AS}\f$ based on cell-centered velocity
  /** 
      \param a_mu: face-centered \f$\mu_{AS}\f$ based on the local velocity field.
      \param a_vel: Cell-centered velocity field.
      \param a_A: Cell-centered flow law coefficient (glenn's A) field
      \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                          (NULL if no coarser level)
      \param a_nRefCrse -- refinement ratio to next coarser level
      \param a_coordSys: LevelSigmaCS object containing the geometry of this patch
      \param a_ghostVect: how the boxes on which we want to compute mu compare
                         with the DisjointBoxLayout. (can be negative)
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, 
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const = 0;


  /// creates a new copy of this  ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const = 0;

  /// compute cell-centered epsilon^2
  void computeStrainRateInvariant(LevelData<FArrayBox>& a_epsilonSquared,
                                  const LevelData<FArrayBox>& a_velocity,
                                  const LevelData<FArrayBox>* a_crseVel,
                                  int nRefCrse,
                                  const LevelSigmaCS& a_coordSys,
                                  const IntVect& a_ghostVect = IntVect::Zero) const;

  /// compute face-centered epsilon^2 based on cell-centered velocity
  void computeStrainRateInvariantFace(LevelData<FluxBox>& a_epsilonSquared,
                                      LevelData<FArrayBox>& a_velocity,
                                      const LevelData<FArrayBox>* a_crseVel,
                                      int a_nRefCrse,
                                      const LevelSigmaCS& a_coordSys,
                                      const IntVect& a_ghostVect = IntVect::Zero) const;  

  /// compute cell-centered epsilon^2 and velocity gradients.
  /// no more costly than just computing epsilon^2
  void computeStrainRateInvariant(LevelData<FArrayBox>& a_epsilonSquared,
				  LevelData<FArrayBox>& a_gradVelocity,
                                  const LevelData<FArrayBox>& a_velocity,
                                  const LevelData<FArrayBox>* a_crseVel,
                                  int a_nRefCrse,
                                  const LevelSigmaCS& a_coordSys,
                                  const IntVect& a_ghostVect = IntVect::Zero) const;

  /// compute face-centered epsilon^2 and velocity gradients
  /// based on cell-centered velocity. no more costly than just computing epsilon^2
  void computeStrainRateInvariantFace(LevelData<FluxBox>& a_epsilonSquared,
				      LevelData<FluxBox>& a_gradVelocity,
                                      LevelData<FArrayBox>& a_velocity,
                                      const LevelData<FArrayBox>* a_crseVel,
                                      int a_nRefCrse,
                                      const LevelSigmaCS& a_coordSys,
                                      const IntVect& a_ghostVect = IntVect::Zero) const;  

  //virtual Real power() const = 0;

protected:

};
 



/// rate factor A(T) in (e.g) Glen's law  
class RateFactor 
{

public:
  
  virtual ~RateFactor(){;}

  virtual void computeA(FArrayBox& a_A, 
			const FArrayBox& a_thetaStar,
			const FArrayBox& a_pressure,
			const Box& a_box) const = 0;

  virtual RateFactor* getNewRateFactor() const = 0;


} ;

///Constant Rate Factor
/**
   Allows A to be specified as a constant, useful for isothermal
   ice sheets.
*/
class ConstantRateFactor : public RateFactor  
{
  Real m_A; // the constant rate factor
public:
  
  ConstantRateFactor() : m_A(9.2e-18) // A measured in Pa^{n} a^{-1}
  {
  }
 
  ConstantRateFactor(Real a_A) : m_A(a_A)
  {
  }


  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_theta, 
		const FArrayBox& a_pressure,
		const Box& a_box) const
  {
    a_A.setVal(m_A, a_box, 0);
  }

  RateFactor* getNewRateFactor() const 
  {
    ConstantRateFactor* newPtr = new ConstantRateFactor(m_A);
    return static_cast<RateFactor*>(newPtr);
  }


};

///Arrhenius Rate Factor
/**
   The Arrhenius rate factor is used widely in ice sheet modelling, e.g 
   (Pattyn, 2003)
 */
class ArrheniusRateFactor : public RateFactor  
{
  /// power law exponent
  Real m_n;
  /// Pattyn's `"enhancement factor"
  Real m_enhance;
  /// flow rate factor
  Real m_B0;
  /// limit temperature in flow-rate factor
  Real m_theta_r;
  /// flow rate exponent
  Real m_K;
  /// flow rate factor
  Real m_C;
  /// universal gas constant
  Real m_R;
  /// activation energy for creep 
  // (insert funny comment about lazy creeps here...)
  Real m_Q;

public:
  ArrheniusRateFactor();
 
  void setDefaultParameters();
 
  
  void setParameters(Real a_n,
		     Real a_enhance,
		     Real a_B0,
		     Real a_theta_r,
		     Real a_K,
		     Real a_C,
		     Real a_R,
		     Real a_Q);

  //compute the (temperature-dependent) rate factor A
  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_thetaStar,
		const FArrayBox& a_pressure,
		const Box& a_box) const;

  RateFactor* getNewRateFactor() const; 
  

};

/// Paterson Rate Factor
/**
   Another  Arrhenius law, this one from
   Cuffey and Paterson,  2010, 'The physics of glaciers', p72-74
   (a revision of earlier books by Paterson) 

   A = E*A0*exp(-Qc / R * [1/T* - 1/T0*]

   Th is the usual pressure-corrected temperature
   Th = T + 7e−8 * P
   T0* = T0  + 7e−8 * P (T0 = 263, P = pressure)
   Qc = Qm = if Th < T∗; Qc = Qp if Th > T∗.
  
   Defaults: 

   A0 = 3.5e-25 Pa−3 s−1 * secondsperyear
   Qp = 115 kJ mol−1
   Qm = 6e4 kJ mol−1
   R=8.314 J

 */
class PatersonRateFactor : public RateFactor  
{
  
  /// enhancement factor
  Real m_E;
  /// flow rate factor
  Real m_A0;
  /// limit temperature in flow-rate factor
  Real m_T0;
  /// universal gas constant
  Real m_R;
  /// activation energies for creep 
  // (insert funny comment about two lazy creeps here...)
  Real m_Qm, m_Qp;

public:

  PatersonRateFactor();
 
  void setDefaultParameters();
 
  void setA0(Real a_A0)
  {
    m_A0 = a_A0;
  }

  void setParameters(Real a_E, Real a_A0, Real a_T0, Real a_R, 
		     Real a_Qm, Real a_Qp);

  //compute the (pressure-corrected temperature- 
  //and pressure- dependent) rate factor A
  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_thetaPC, 
		const FArrayBox& a_pressure,
		const Box& a_box) const;

  RateFactor* getNewRateFactor() const; 
  

};


/// Zwinger Rate Factor
/**
   Another  Arrhenius law, this one from
   Zwinger (2007) Annals of Glaciology vol 45 p 29
   and used in (older versions of?) Elmer/Ice. 

   A = E*A0(T)*exp(-Qc / (R * T) ]
   where Qc = Qm = if theta < theta0; Qc = Qp otherwise.
   and   A0 = Am = if theta < theta0; A0 = Ap otherwise.
   Am is computed so that A is continuous with T 

   Th is the usual pressure-corrected temperature
   Th = T + 7e−8 * P
   T0* = T0  + 7e−8 * P (T0 = 263, P = pressure)
   Qc = Qm = if Th < T∗; Qc = Qp if Th > T∗.
  
   Defaults: 

   A0 = 1.916e3 Pa−3 s−1 * secondsperyear
   Qp = 139 kJ mol−1
   Qm = 6e4 kJ mol−1
   R=8.314 J

 */
class ZwingerRateFactor : public RateFactor  
{
  
  /// enhancement factor
  Real m_E;
  /// flow rate factor
  Real m_A0;
  /// limit temperature in flow-rate factor
  Real m_T0;
  /// universal gas constant
  Real m_R;
  /// activation energies for creep 
  // (insert funny comment about two lazy creeps here...)
  Real m_Qm, m_Qp;

public:

  ZwingerRateFactor();
 
  void setDefaultParameters();
 

  void setParameters(Real a_E, Real a_A0, Real a_T0, Real a_R, 
		     Real a_Qm, Real a_Qp);

  //compute the (pressure-corrected temperature- 
  //and pressure- dependent) rate factor A
  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_thetaPC, 
		const FArrayBox& a_pressure,
		const Box& a_box) const;

  RateFactor* getNewRateFactor() const; 
  

};


/// implementation of Glen's flow law constitutive relation
/** The  GlensFlowRelation class is publicly derived from the 
    ConstitutiveRelation Class and uses Glen's flow law, as described
    in Dukowicz, Et Al (2009)
*/
class GlensFlowRelation: public ConstitutiveRelation
{
public:
  ///
  GlensFlowRelation();

  virtual ~GlensFlowRelation();


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

#if 0 // getting rid of computeFaceFluxVelocity
#if BISICLES_Z == BISICLES_LAYERED
  // in the 2D case, we need two velocities : a 2D, 2 component field which appears
  // in the thickness evolution equation and a 3D field for temperature advection.
  // The 3D field is split into a LevelData<FluxBox> a_layerXYFaceXYVel  of 
  // horizontal velocities at horizontal faces,
  // each FluxBox having n components, and a LevelData<FArrayBox> a_layerSFaceXYVel of 
  // horizontal velocities at sigma faces, each FArrayBox having (n+1)*SpaceDim 
  // components, and 0 < direction index < SpaceDim running faster than layer index 
  
  // In the SSA case, the velocity doesn't vary in z, so we just can copy it
  virtual void computeFaceFluxVelocity(const LevelData<FArrayBox>& a_cellVel,
                                       const LevelData<FArrayBox>* a_crseVel,
                                       int nRefCrse,
				       const LevelSigmaCS& a_coordSys,
				       const DisjointBoxLayout& a_grids,
				       const ProblemDomain& a_domain,
				       const LevelData<FArrayBox>& a_A,
				       const LevelData<FArrayBox>& a_sTheta,
				       const LevelData<FArrayBox>& a_bTheta,
				       LevelData<FluxBox>& a_faceVel,
				       LevelData<FluxBox>& a_layerXYFaceXYVel,
				       LevelData<FArrayBox>& a_layerSFaceXYVel) const
  {
    CellToEdge(a_cellVel, a_faceVel);
    
    for (int j = 0; j < a_layerXYFaceXYVel.nComp(); ++j)
      {
	a_faceVel.copyTo(Interval(0,0), a_layerXYFaceXYVel, Interval(j,j));
      }

    for (int j = 0; j < a_layerSFaceXYVel.nComp(); j+=SpaceDim)
      {
	a_cellVel.copyTo(Interval(0,SpaceDim-1), a_layerSFaceXYVel, 
			 Interval(j,j+SpaceDim-1));
      }


  }  
#elif BISICLES_Z == BISICLES_FULLZ
  virtual void computeFaceFluxVelocity(const LevelData<FArrayBox>& a_cellVel,
                                       const LevelData<FArrayBox>* a_crseVel,
                                       int nRefCrse,
				       const LevelSigmaCS& a_coordSys,
				       const DisjointBoxLayout& a_grids,
				       const ProblemDomain& a_domain,
				       const LevelData<FArrayBox>& a_A,
				       LevelData<FluxBox>& a_faceVel) const
  {
    MayDay::Error(" 3D face velocities not done");
  }
  
#endif
#endif  // getting rid of computeFaceFluxVelocity


#if 0
  /// computes face-centered \f$\mu_{AS}\f$ for a single face direction 
  /**
     a_mu: face-centered \f$\mu_{AS}\f$ based on the local velocity field.
     a_vel: Cell-centered velocity field.
     a_coordSys: SigmaCS object containing the geometry of this patch.
     a_box: face-centered box over which to do this computation
     a_dir: face direction
  */
  virtual void computeFaceMu(FArrayBox& a_mu,
                             const FArrayBox& a_vel, 
                             const FArrayBox& a_A,
                             const LevelSigmaCS& a_coordSys,
                             const Box& a_box,
                             int a_dir) const;
#endif

  virtual void setDefaultParameters();

  virtual void setParameters(Real a_n,
                             Real a_epsSqr0, 
			     Real a_delta);

  /// creates a new copy of this ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

  /* these are all parameters used to compute the temprature-dependent
     part of the effective viscosity mu_0, based on (Pattyn, 2003)
  */
  
  /// Power law exponent
  Real m_n;

  virtual Real power() const {return m_n;}

  /// Rate factor object
  //RateFactor* m_rateFactor;

  /// small numbers to ensure that viscosity remains finite as
  /// strain rates approach 0 or inifinity
  Real m_epsSqr0, m_delta;

protected:
  /// utility function to compute flow law coefficient dependent part of Glens's
  /// law (ie A^{1/n})
  void computeMu0(FArrayBox& a_mu0,
                  const FArrayBox& a_A,
                  const Box& a_box) const;
  

};


/// implementation of constant-mu constitutive relation
/** The constMuRelation class is publicly derived from the 
    ConstitutiveRelation Class and sets mu to be a constant value
*/
class constMuRelation: public ConstitutiveRelation
{
public:
  ///
  constMuRelation();

  virtual ~constMuRelation();


  /// computes cell-centered \f$\mu_{AS}\f$ based on the cell-centered velocity
  /**
     a_mu -- \f$\mu_{AS}\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (glenn's A) field
     a_coordSys:  SigmaCS object containing the geometry of this patch.
     a_box: cell-centered box over which to do this computation
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
     \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     \param a_nRefCrse -- refinement ratio to next coarser level
      a_A: Cell-centered flow law coefficient (glenn's A) field
      a_coordSys: SigmaCS object containing the geometry of this patch.
      a_box: cell-centered box over which to do this computation
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, 
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse, 
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;


#if 0
  /// computes face-centered \f$\mu_{AS}\f$ for a single face direction 
  /**
     a_mu: face-centered \f$\mu_{AS}\f$ based on the local velocity field.
     a_vel: Cell-centered velocity field.
     \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     \param a_nRefCrse -- refinement ratio to next coarser level
     a_coordSys: SigmaCS object containing the geometry of this patch.
     a_box: face-centered box over which to do this computation
     a_dir: face direction
  */
  virtual void computeFaceMu(FArrayBox& a_mu,
                             const FArrayBox& a_vel, 
                             const FArrayBox& a_A,
                             const LevelSigmaCS& a_coordSys,
                             const Box& a_box,
                             int a_dir) const;
#endif
#if 0 // getting rid of computeFaceFluxVelocity
  #if BISICLES_Z == BISICLES_LAYERED
  // in the 2D case, we need two velocities : a 2D, 2 component field which appears
  // in the thickness evolution equation and a 3D field for temperature advection.
  // The 3D field is split into a LevelData<FluxBox> a_layerXYFaceXYVel  of 
  // horizontal velocities at horizontal faces,
  // each FluxBox having n components, and a LevelData<FArrayBox> a_layerSFaceXYVel of 
  // horizontal velocities at sigma faces, each FArrayBox having (n+1)*SpaceDim 
  // components, and 0 < direction index < SpaceDim running faster than layer index 
  
  // In the const mu case, 
  // the velocity doesn't vary in z, so we just can copy it
  virtual void computeFaceFluxVelocity(const LevelData<FArrayBox>& a_cellVel,
                                       const LevelData<FArrayBox>* a_crseVel,
                                       int nRefCrse,
				       const LevelSigmaCS& a_coordSys,
				       const DisjointBoxLayout& a_grids,
				       const ProblemDomain& a_domain,
				       const LevelData<FArrayBox>& a_A,
				       const LevelData<FArrayBox>& a_sTheta,
				       const LevelData<FArrayBox>& a_bTheta,
				       LevelData<FluxBox>& a_faceVel,
				       LevelData<FluxBox>& a_layerXYFaceXYVel,
				       LevelData<FArrayBox>& a_layerSFaceXYVel) const
  {
    CellToEdge(a_cellVel, a_faceVel);
    
    for (int j = 0; j < a_layerXYFaceXYVel.nComp(); ++j)
      {
	a_faceVel.copyTo(Interval(0,0), a_layerXYFaceXYVel, Interval(j,j));
      }

    for (int j = 0; j < a_layerSFaceXYVel.nComp(); j+=SpaceDim)
      {
	a_cellVel.copyTo(Interval(0,SpaceDim-1), a_layerSFaceXYVel, 
			 Interval(j,j+SpaceDim-1));
      }


  }  
#elif BISICLES_Z == BISICLES_FULLZ
  virtual void computeFaceFluxVelocity(const LevelData<FArrayBox>& a_cellVel,
                                       const LevelData<FArrayBox>* a_crseVel,
                                       int nRefCrse,
				       const LevelSigmaCS& a_coordSys,
				       const DisjointBoxLayout& a_grids,
				       const ProblemDomain& a_domain,
				       const LevelData<FArrayBox>& a_A,
				       LevelData<FluxBox>& a_faceVel) const
  {
    MayDay::Error(" 3D face velocities not done");
  }
  
#endif
#endif  // getting rid of computeFaceFluxVelocity
  virtual void setDefaultParameters();

  virtual void setConstVal(Real a_mu) {m_mu = a_mu;}


  /// creates a new copy of this  ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

  /* these are all parameters used to compute the temprature-dependent
     part of the effective viscosity mu_0, based on (Pattyn, 2003)
  */
  
  // constant value of mu
  Real m_mu;

  virtual Real power() const{return 1.0;}

protected:

};



#include "NamespaceFooter.H"

#endif
