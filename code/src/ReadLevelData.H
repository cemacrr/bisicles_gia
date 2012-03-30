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
// ReadLevelData.H
//===========================================================================

#ifndef _READLEVELDATA_H_
#define _READLEVELDATA_H_
#include "AMRIO.H"
#include "NamespaceHeader.H"

//fill a LevelData<FArrayBox> from data stored in an AMR file 
void readLevelData(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_data,
		   Real& a_dx,
		   const std::string a_file,
		   const Vector<std::string>& a_names,
		   int a_nComp);

#include "NamespaceFooter.H"

#endif
