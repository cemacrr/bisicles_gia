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
//  FillFromReference.H
// ============
//
// function to 
//

#ifndef _FILLFROMREFERENCE_H_
#define _FILLFROMREFERENCE_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "RealVect.H"


#include "NamespaceHeader.H"

/// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference FArrayBox
/**
   
 */

void FillFromReference(LevelData<FArrayBox>& a_destData,
		       const FArrayBox& a_srcData,
		       const RealVect&  a_destDx,
		       const RealVect&  a_srcDx,
		       const IntVect& a_srcGhost,
		       const bool a_verbose);



/// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference  LevelData<FArrayBox>
/**
   
 */
void FillFromReference(LevelData<FArrayBox>& a_destData,
		       const LevelData<FArrayBox>&  a_srcData,
		       const RealVect&  a_destDx,
		       const RealVect&  a_srcDx,
		       const bool a_verbose);


#include "NamespaceFooter.H"

#endif
