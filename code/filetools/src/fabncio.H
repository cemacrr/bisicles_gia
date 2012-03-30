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
// ncio.H
// read/write FABs from netcdf files
//===========================================================================

#include "FArrayBox.H"
#include "NamespaceHeader.H"

#ifndef _FABNCIO_H_
#define _FABNCIO_H_

void writeNetCDF(const std::string& a_file,
		 const Vector<std::string>& a_names,
		 const FArrayBox& a_fab, 
		 const Real& a_dx);

void readNetCDF(const std::string& a_file,
		const Vector<std::string>& a_var,
		FArrayBox& a_fab,
		Real& a_dx);

#include "NamespaceFooter.H"
#endif
