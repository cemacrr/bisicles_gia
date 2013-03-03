#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ROOTFINDER_H_
#define _ROOTFINDER_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "REAL.H"
#include "IndexTM.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"
class RootFinder
{
public:
  
  RootFinder(BaseIF* a_function);

Real Brent(const IndexTM<Real,SpaceDim>& a_loPt,
           const IndexTM<Real,SpaceDim>& a_hiPt);

BaseIF* m_function;

};

#include "NamespaceFooter.H"
#endif
