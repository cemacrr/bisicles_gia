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
//  FortranInterfaceIBC.H
// ============
//
// PhysIBC-derived class for interfacing with fortran -- IBC fields like 
// initial ice thickness and topography are passed in as fortran arrays, which 
// are then used to set values in the BISICLES code
//

#ifndef _FORTRANTHICKNESSIBC_H_
#define _FORTRANTHICKNESSIBC_H_

#include "IceThicknessIBC.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"

class AmrIce;

/// Physical/domain initial and boundary conditions
/**
   
 */
class FortranInterfaceIBC : public IceThicknessIBC
{

 
public:
  /// Constructor
  /**
   */
  FortranInterfaceIBC();

  /// Destructor
  /**
   */
  virtual ~FortranInterfaceIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  virtual void  setFAB(Real* a_data_ptr,
		       const int* a_dimInfo,
                       const int* boxlo, const int* boxhi, 
		       const Real* a_dew, const Real* a_dns,
		       const IntVect& a_nGhost,
		       FArrayBox& a_fab,
		       const bool a_nodal = false);

  /// nGhost is the number of ghost cells in in the input data (which are ignored)
  virtual void setThickness(Real* a_data_ptr,
			    const int* a_dimInfo,
                            const int* boxlo, const int* boxhi, 
			    const Real* a_dew, const Real* a_dns,
                            const IntVect& a_nGhost = IntVect::Zero, 
			    const bool a_nodal = false);

  virtual void setTopography(Real* a_data_ptr,
			     const int* a_dimInfo,
                             const int* boxlo, const int* boxhi, 
			     const Real* a_dew, const Real* a_dns,
                             const IntVect& a_nGhost = IntVect::Zero, 
			     const bool a_nodal = false);

  /// regions where we artificially set thickness to zero
  /** this is done in setThickness, for lack of a better place, so 
      this needs to be set before setThickness is called. 
      a_clearRegions defines logically-rectangular regions where 
      the thickness is artificially set to zero. These regions are 
      defined relative to the original input data (i.e. before any shifting 
      due to ghost cells)
  */
  virtual void setThicknessClearRegions(const Vector<Box>& a_clearRegions);

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

  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
   */
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio) ;

  /// update the topgraphy in LevelSigmaCS but *not* the thickness
  virtual void regridIceGeometry(LevelSigmaCS& a_coords,
				 const RealVect& a_dx,
				 const RealVect& a_domainSize,
				 const Real& a_time, 
				 const LevelSigmaCS* a_crseCoords,
				 const int a_refRatio) ;

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

  void setVerbose(bool a_verbose) {m_verbose = a_verbose;}

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
  
  /// set grids using Boxes passed in from Glimmer-CISM
  /** creates a DisjointBoxLayout using the grid boxes and 
      processor distribution used by CISM. */
  virtual void setGrids(const Box& a_gridBox);

  /// has setGrids been called?
  bool gridsSet() const {return m_gridsSet;}


  /// utility function to fill in holes in topography
  /** looks for isolated values of holeVal and replaces then
      with average of neighbors */
  virtual void fillTopographyHoles(Real holeVal);

  virtual void checkOK() const;

protected:

  // ice sheet thickness at domain boundaries
  Real          m_boundaryThickness;

  RealVect m_domainSize;

  bool m_verbose;

  // have bc's been set up?
  bool m_isBCsetUp;

  /// set up bc's based on parmparse inputs
  void setupBCs();

  bool m_extrapBoundary; //extrapolate the boundary geoemtry? or reflect it?

  BCFunc m_velBCs;

  DisjointBoxLayout m_grids;
  bool m_gridsSet;

  // ice sheet thickness from glimmer
  /// alias to data array passed in from glimmer-CISM
  FArrayBox m_inputThickness;
  RealVect m_inputThicknessDx;
  IntVect m_thicknessGhost; 
  /// distributed version of thickness
  RefCountedPtr<LevelData<FArrayBox> > m_inputThicknessLDF;

  /// regions where we will artificially set ice thickness to zero
  Vector<Box> m_thicknessClearRegions;

  /// alias to data array passed in from glimmer-CISM
  FArrayBox m_inputTopography;
  RealVect m_inputTopographyDx;
  IntVect m_topographyGhost;
  /// distributed version of topography
  RefCountedPtr<LevelData<FArrayBox> > m_inputTopographyLDF;

private:
  // Disallowed for all the usual reasons
  void operator=(const FortranInterfaceIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FortranInterfaceIBC(const FortranInterfaceIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
