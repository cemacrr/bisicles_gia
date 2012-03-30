#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCAMRSOLVER_H_
#define _PETSCAMRSOLVER_H_

#include "PetscLinearSolver.H"
#include "ViscousTensorOp.H"
#include "LinearSolver.H"

#include "parstream.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

/// very simple implementation of Petsc AMR solver
/** Basic implementation of a PETSC AMR solver which takes a Vector 
    single-level PETSC-based solvers and a Chombo-based AMR solver
    and does a composite solve. For now, it is up to the calling context 
    to ensure that the solvers are consistently defined.
*/

template <class T> 
class PetscAMRSolver<T>
{

};

#include "NamespaceFooter.H"
