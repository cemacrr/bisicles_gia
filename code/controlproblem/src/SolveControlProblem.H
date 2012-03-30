#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// SolveControlProblem.H
// contruct and solve a control problem
// by a method related to Vieli and Payne, Ann Glaciol, 2003. vol 26
// assumes a 2D SSA model or L1L2 model
//===========================================================================

#ifndef SOLVE_CONTROL_PROBLEM_H
#define SOLVE_CONTROL_PROBLEM_H

#include "ConstitutiveRelation.H"
#include "BasalFrictionRelation.H"
#include "LevelDataIBC.H"
#include "Vector.H"
#include "AMRIO.H"
#include "LevelDataOps.H"

//set up an AMR hierarchy, from a single level of data
//and a vector of refinement ratios.
//The level data is read from an hdf5 file and
//resides on level a_level of the hierarchy. 
void initAMR( Vector<RefCountedPtr<LevelSigmaCS > >& a_coords,
	      Vector<LevelData<FArrayBox>*>& a_C,
	      Vector<LevelData<FArrayBox>*>& a_A,
	      Vector<LevelData<FArrayBox>*>& a_velObs,
	      Vector<LevelData<FArrayBox>*>& a_vel,
	      Vector<LevelData<FArrayBox>*>& a_DivUHObs,
	      Vector<DisjointBoxLayout>& a_grids,
	      Vector<ProblemDomain>& a_domain,
	      Vector<RealVect>& a_dx,
	      LevelDataIBC* a_bc,
	      const std::string a_filename,
	      const Vector<int>& a_refRatio,
	      const int& a_level,
	      const ProblemDomain& a_crseDomain);


//set up and solve the forward  problem
void solveForward(const std::string a_infile,
	   const std::string a_outfile,
	   const ProblemDomain& a_crseDomain,
	   ConstitutiveRelation* a_constRelPtr, 
	   BasalFrictionRelation* a_basalFrictionRelPtr);



/// define and solve a control problem .
/// 
/// miminise L = (u - u_obs)^2 + q (fwd(u)) with respect to C
///
/// the foward problem fwd(u)=0 is solved for a velcity field u 
/// 
/// div(H*mu(u)*(grad(u)+...)) - C nu(u)*u = rho * g * H * grad(s) 
///
/// the adjoint problem adj(q) = 0 is solved for a largange multiplier q
/// 
/// div(H*mu(u)*(grad(q)+...)) - C nu(u)*q = rhs_adj
///
/// then dL/dC = q * u 

/// a_constRelPtr : constitutuve relation for mu calculation
/// a_basalFrictionRealPtr : constitutuve relation for nu calculation
/// a_levelThck : LevelData<FArrayBox>& of ice thickness, spans the domain
/// a_levelTopg : topography data
/// a_levelVelObs : velocity observations
/// a_levelC : on entry, initial guess for C. on succesful exit, 
///           a C which reduces (u - u_obs)^2.  
void solveControl(LevelData<FArrayBox>& a_levelC,
		  const LevelData<FArrayBox>& a_levelVelObs,
		  const LevelData<FArrayBox>& a_levelVelCoef,
		  const LevelData<FArrayBox>& a_levelDivUHObs,
		  const ProblemDomain& a_problemDomain,
		  const DisjointBoxLayout& a_grids,
#if BISICLES_Z == BISICLES_LAYERED
		  const Vector<Real>& a_faceSigma,
#endif
		  IceThicknessIBC* a_bc,
		  ConstitutiveRelation* a_constRelPtr, 
		  BasalFrictionRelation* a_basalFrictionRelPtr,
		  const RealVect& a_baseDx, const RealVect& a_dataDx);

//construct and solve, for u, the forward problem
// div(H*mu(u)*(grad(u)+...)) - C nu(u)*u = rho * g * H * grad(s)
void solveForward(LevelData<FArrayBox>& a_levelVel,
		  LevelData<FArrayBox>& levelDivUH,
		  const LevelData<FArrayBox>& a_levelC,
		  const ProblemDomain& a_problemDomain,
		  const DisjointBoxLayout& a_grids,
#if BISICLES_Z == BISICLES_LAYERED
		  const Vector<Real>& a_faceSigma,
#endif
		  IceThicknessIBC* a_bc,
		  ConstitutiveRelation* a_constRelPtr, 
		  BasalFrictionRelation* a_basalFrictionRelPtr,
		  const RealVect& a_baseDx, const RealVect& a_dataDx);

//construct and solve, for u, the forward problem
// div(H*mu(u)*(grad(u)+...)) - C nu(u)*u = rho * g * H * grad(s)
void solveForward(Vector<LevelData<FArrayBox>* >& a_u,
		  bool a_linear,
		  const Vector<LevelData<FArrayBox>* >& a_C,
		  const Vector<LevelData<FArrayBox>* >& a_C0,
		  const Vector<LevelData<FArrayBox>* >& a_A,
		  const Vector<LevelData<FArrayBox>* >& a_rhs,
		  Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		  const ProblemDomain& a_crseDomain,
		  const Vector<DisjointBoxLayout>& a_grids,
		  const Vector<int> a_refRatios,
		  IceThicknessIBC* a_bc,
		  ConstitutiveRelation* a_constRelPtr, 
		  BasalFrictionRelation* a_basalFrictionRelPtr,
		  const RealVect& a_crseDx);



#endif
