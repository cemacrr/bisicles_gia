#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef EBAMR_ICE_H
#define EBAMR_ICE_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "IceVelocitySolver.H"
#include "RealVect.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "AmrIce.H"
#include "SurfaceFlux.H"
#include "BasalFriction.H"
#include "IceThicknessIBC.H"
#include "IceTemperatureIBC.H"
#include "CalvingModel.H"
#include "MuCoefficient.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "NamespaceHeader.H"

/// class to manage non-subcycled EBAMR ice-sheet model
/**
 */
class EBAmrIce:public AmrIce
{

public:
  /// Default constructor
  /** At the moment, there is only one constructor, which defines itself
      based on an inputs file.  (This may change, of course)
  */
  EBAmrIce();

  /// destructor
  /** destructor
   */
  ~EBAmrIce();

protected:
  /// define nonlinear ellipic solver for computing velocity field
  void defineSolver();

  void velInitialGuess( Vector<LevelData<FArrayBox>* >& a_vectC0,
                       const Real                    & a_convergenceMetric);

  /// solve for velocity field (or just set up some auxilliary quantities)
  void solveVelocityField(Real a_convergenceMetric = -1.0);

  /// compute RHS for velocity field solve
  /** also sets beta to zero where ice is floating
   */
  void defineVelRHS(Vector<LevelData<FArrayBox>* >& a_vectRhs,
                    Vector<LevelData<FArrayBox>* >& a_vectC,
		    Vector<LevelData<FArrayBox>* >& a_vectC0);

  /// set basal friction coefficient (beta) prior to velocity solve
  void setBasalFriction(Vector<LevelData<FArrayBox>* >& a_vectBeta);

  /// set mu coefficient (phi) prior to velocity solve
  void setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef, 
			Vector<LevelData<FluxBox>* >& a_faceMuCoef);


};  
  

#include "NamespaceFooter.H"

#endif


