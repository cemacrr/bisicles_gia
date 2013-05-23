#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCAMRSOLVERSOLVER_H_
#define _PETSCAMRSOLVERSOLVER_H_

#include "PetscCompGridVTO.H"
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"

#include "NamespaceHeader.H"

class JFNKOp;
// solver class for VTO operators
class PetscAMRSolver
{
public:
  PetscAMRSolver();
  ///
  /**
     Solve a_op*a_phi = a_rhs using the PETSC matrix free functions
     The preconditioner (for which a matrix is formed) need not be 
     the same as the actual operator.
  */
  virtual PetscErrorCode solve_mfree( Vector<LevelData<FArrayBox>*>& a_phi, 
				      const Vector<LevelData<FArrayBox>*>& a_rhs, 
				      JFNKOp *a_op );

  virtual PetscErrorCode solve( Vector<LevelData<FArrayBox>*>& a_phi, 
				const Vector<LevelData<FArrayBox>*>& a_rhs );
  
  static PetscErrorCode apply_mfree(Mat A, Vec x, Vec f);  
  
  JFNKOp *m_op_mfree;
  Vector<LevelData<FArrayBox>*> m_phi_mfree;
  Vector<LevelData<FArrayBox>*> m_Lphi_mfree;
  bool m_mfree_homogeneous;
  PetscCompGridVTO m_petscCompMat; // this is specific for VTO for now, should make PetscCompGrid*
};

#include "NamespaceFooter.H"
#endif
