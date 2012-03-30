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
//  LevelDataIBC.H
// ============
//
// PhysIBC-derived class which stores initial topography and thickness data
// and imposes either periodic or reflection boundary conditions

#ifndef _LEVELDATAIBC_H_
#define _LEVELDATAIBC_H_

#include "IceThicknessIBC.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   
 */
class LevelDataIBC : public IceThicknessIBC
{
public:
  /// Constructor
  /**
   */
  LevelDataIBC(RefCountedPtr<LevelData<FArrayBox> > a_thck, 
	       RefCountedPtr<LevelData<FArrayBox> > a_topg,
	       const RealVect& a_dx);

  LevelDataIBC();

  /// Destructor
  /**
   */
  virtual ~LevelDataIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new IceThicknessIBC object with 
     m_isDefined = false (i.e., its define() must be called before it
     is used). 
   */
  virtual IceThicknessIBC* new_thicknessIBC();


  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U);

  void regridIceGeometry(LevelSigmaCS& a_coords,
			 const RealVect& a_dx,
			 const RealVect& a_domainSize,
			 const Real& a_time,
			 const LevelSigmaCS* a_crseCoords,
			 const int a_refRatio);

  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
   */
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time,
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio);


  /// Set boundary fluxes
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time);

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
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  /// return boundary condition for Ice velocity solve
  /** 
   */
  virtual BCHolder velocitySolveBC();

  /// if appropriate, modify velocity solve RHS in a problem-dependent way. 
  virtual void modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
                                 LevelSigmaCS& a_coords,
                                 const ProblemDomain& a_domain,
                                 Real a_time, Real a_dt)
  {;}


  /// set non-periodic ghost cells for surface height z_s. 
  virtual void setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                   LevelSigmaCS& a_coords,
                                   const ProblemDomain& a_domain,
                                   const RealVect& a_dx, 
                                   Real a_time, Real a_dt);
  /// set non-periodic geometry (thickness, topography) ghost cells
  virtual void setGeometryBCs(LevelSigmaCS& a_coords,
			      const ProblemDomain& a_domain,
			      const RealVect& a_dx,
			      Real a_time, Real a_dt);

private:

  RefCountedPtr<LevelData<FArrayBox> > m_thck;
  RefCountedPtr<LevelData<FArrayBox> > m_topg;
  RealVect m_dx;
  BCFunc m_velBCs;
  bool m_isBCsetUp;
  bool m_verbose;
  void setupBCs();

private:
  // Disallowed for all the usual reasons
  void operator=(const LevelDataIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelDataIBC(const LevelDataIBC& a_input)
    :m_thck(a_input.m_thck),m_topg(a_input.m_topg),m_dx(a_input.m_dx)
  {
    MayDay::Error("invalid operator");
  }

  void groundCoarse(LevelSigmaCS& a_coords, int a_refRatio);

};

#include "NamespaceFooter.H"
#endif
