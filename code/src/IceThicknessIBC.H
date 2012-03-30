#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  IceThicknessIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _ICETHICKNESSIBC_H_
#define _ICETHICKNESSIBC_H_

#include "PhysIBC.H"
#include "LevelSigmaCS.H"
#include "BCFunc.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions for ice-sheet problems
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for an ice sheet. Derived from PhysIBC, adds
   functions for initializing SigmaCS.
 */
class IceThicknessIBC : public PhysIBC
{
public:
  /// Constructor
  /**
   */
  IceThicknessIBC() {}

  /// Destructor
  /**
   */
  virtual ~IceThicknessIBC() {}

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object. Default is to just call
     PhysIBC::define
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx)
  {PhysIBC::define(a_domain, a_dx);}

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysIBC* new_physIBC() 
  {return static_cast<PhysIBC*>(new_thicknessIBC());}

  /// same as new_physIBC, except pointer cast to an IceThicknessIBC
  virtual IceThicknessIBC* new_thicknessIBC() = 0;
 

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

  /// set up initial ice state
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio) = 0;
  

  /// set up topography, etc at regrid time
  /*  just like initializeIceGeometry except that we don't set 
      ice thickness, since it's assumed to be set using current state.
      
      default is to do nothing, in which case we assume that the current
      state will be interpolated in AmrIce::regrid (to enable isostasy, 
      for example)
  */
  virtual void regridIceGeometry(LevelSigmaCS& a_coords,
                                 const RealVect& a_dx,
                                 const RealVect& a_domainSize,
                                 const Real& a_time,
				 const LevelSigmaCS* a_crseCoords,
				 const int a_refRatio) {;}
  /// Set boundary fluxes
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time) = 0;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  virtual
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time) = 0;

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time) = 0;

  /// return boundary condition for Ice velocity solve
  /** 
   */
  virtual BCHolder velocitySolveBC() = 0;


  /// fill ghost cells on velocity
  /** this is a separate function from the BCHolder returned by 
      the velocitySolveBC function to take into account the case 
      where the boundary condition used in the solve is different 
      from the one used when working with the velocity. One example 
      of this is the case where we incorporate the normal-stress boundary 
      condition into the RHS for the velocity solve; in that case, the 
      BCHolder will be a homogeneous BC. However we want to use the 
      non-homogeneous form of the BC's when we do things like compute 
      derivatives of the velocity at physical domain boundaries 
      (say, for example, when computing the strain invariant in order 
      to compute the viscosity).  This function fills ghost cells on the 
      boundary in the generic case. Default implementation simply calls 
      the BCHolder implementation... */             
  virtual void velocityGhostBC(LevelData<FArrayBox>& a_velocity,
                               const LevelSigmaCS& a_coords,
                               const ProblemDomain& a_domain,
                               Real a_time)
  {
    BCHolder velBC = velocitySolveBC();
    const RealVect& dx = a_coords.dx();
    for (DataIterator dit=a_velocity.dataIterator(); dit.ok(); ++dit)
      velBC(a_velocity[dit],a_velocity.getBoxes()[dit],
            a_domain, dx[0], false, dit());
  }

      
  
  /// if appropriate, modify velocity solve RHS in a problem-dependent way. 
  /** default is to do nothing
   */
  virtual void modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
                                 LevelSigmaCS& a_coords,
                                 const ProblemDomain& a_domain,
                                 Real a_time, Real a_dt) {;}

  /// set non-periodic ghost cells for surface height z_s. 
  /** Default is to do nothing
      note that we pass in ProblemDomain and dx
   */
  virtual void setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                   LevelSigmaCS& a_coords,
                                   const ProblemDomain& a_domain,
                                   const RealVect& a_dx,
                                   Real a_time, Real a_dt) {;}

  
  /// set non-periodic ghost cells for geometry
  /** Default is to do nothing
      note that we pass in ProblemDomain and dx
  */
  virtual void setGeometryBCs(LevelSigmaCS& a_coords,
			      const ProblemDomain& a_domain,
			      const RealVect& a_dx,
			      Real a_time, Real a_dt) {;}
  

  

  /// set AMR grid hierarchy (for the BC's which need this)
  virtual void setGridHierarchy(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                                Vector<ProblemDomain>& a_vectDomain) {;}


  // this is just here to so that FortranInterfaceIBC can re-implement it
  // and I can pester it on a frequent basis to work out who is 
  // corrupting its memory
  virtual void checkOK() const {;}

protected:

  

private:
  // Disallowed for all the usual reasons
  void operator=(const IceThicknessIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  IceThicknessIBC(const IceThicknessIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


/// BCFunction-derived class to handle velocity-solve BC's 
/** in particular, this BCFunction contains the 
    Vector<LevelSigmaCS* > which we need to define the marine boundary 
    condition */
class  IceVelBCFunction : public BCFunction
{
public:
  IceVelBCFunction() {;}  

  IceVelBCFunction(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                   Vector<ProblemDomain>& a_vectDomain) {;}

  
  ~IceVelBCFunction() {;}

  virtual void define(Vector<RefCountedPtr<LevelSigmaCS> >& a_vectCS,
                      Vector<ProblemDomain>& a_vectDomain) =0;

  
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous) =0;


protected:
  Vector<RefCountedPtr<LevelSigmaCS > > m_vectCS;

  Vector<ProblemDomain>                       m_vectDomain;


};


#include "NamespaceFooter.H"
#endif


