
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )

#include <unistd.h>
#define _GLIBCPP_USE_C99 
#endif

#include <cmath>

#include "BaseIF.H"
#include "IndexTM.H"
#include "CH_Attach.H"

#include "RootFinder.H"

#include "NamespaceHeader.H"

RootFinder::RootFinder(BaseIF* a_function)
  :m_function(a_function)
{
}
// returns a number between -1.0 and 1.0
Real RootFinder::Brent(const IndexTM<Real,SpaceDim>& a_loPt,
                       const IndexTM<Real,SpaceDim>& a_hiPt) 

{
  // Max allowed iterations and floating point precision
  const unsigned int MAXITER = 100;
  const Real         EPS   = 3.0e-15;

  unsigned int i;
  Real aPt;
  Real bPt;
  Real c, fa, fb, fc;
  Real d, e;
  Real tol1, xm;
  Real p, q, r, s;

  aPt = -1.0; // a_smallestRoot;
  bPt =  1.0; // a_biggestRoot;

  IndexTM<Real,SpaceDim> physCoordAPt = a_loPt;
  IndexTM<Real,SpaceDim> physCoordBPt = a_hiPt;

  fa = -m_function->value(IndexTM<int,SpaceDim>::Zero,physCoordAPt);
  fb = -m_function->value(IndexTM<int,SpaceDim>::Zero,physCoordBPt);

  //  Init these to be safe
  c = d = e = 0.0;

  if (fb*fa > 0)
  {
    pout() << "fa " << fa << " fb " << fb <<endl;
    MayDay::Abort("IFData::BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
  }

  fc = fb;

  for (i = 0; i < MAXITER; i++)
  {
    if (fb*fc > 0)
    {
      //  Rename a, b, c and adjust bounding interval d
      c = aPt;
      fc  = fa;
      d = bPt - aPt;
      e = d;
    }

    if (Abs(fc) < Abs(fb))
    {
      aPt = bPt;
      bPt = c;
      c = aPt;
      fa  = fb;
      fb  = fc;
      fc  = fa;
    }

    //  Convergence check
    tol1  = 2.0 * EPS * Abs(bPt) + 0.5 * TOLERANCE;
    xm    = 0.5 * (c - bPt);

    if (Abs(xm) <= tol1 || fb == 0.0)
    {
      break;
    }

    if (Abs(e) >= tol1 && Abs(fa) > Abs(fb))
    {
      //  Attempt inverse quadratic interpolation
      s = fb / fa;
      if (aPt == c)
      {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q-r) - (bPt-aPt) * (r-1.0));
        q = (q-1.0) * (r-1.0) * (s-1.0);
      }

      //  Check whether in bounds
      if (p > 0) q = -q;

      p = Abs(p);

      if (2.0 * p < Min(((Real)3.0)*xm*q-Abs(tol1*q), Abs(e*q)))
      {
        //  Accept interpolation
        e = d;
        d = p / q;
      }
      else
      {
        //  Interpolation failed, use bisection
        d = xm;
        e = d;
      }
    }
    else
    {
      //  Bounds decreasing too slowly, use bisection
      d = xm;
      e = d;
    }

    //  Move last best guess to a
    aPt = bPt;
    fa  = fb;

    //  Evaluate new trial root
    if (Abs(d) > tol1)
    {
      bPt = bPt + d;
    }
    else
    {
      if (xm < 0) bPt = bPt - tol1;
      else        bPt = bPt + tol1;
    }

    physCoordBPt[0] = ((1.0 - bPt)/2.0) * a_loPt[0] + ((1.0 + bPt)/2.0) * a_hiPt[0];
    fb = -m_function->value(IndexTM<int,SpaceDim>::Zero,physCoordBPt);
  }

  if (i >= MAXITER)
  {
    std::cerr << "BrentRootFinder: exceeding maximum iterations: "
              << MAXITER
              << "\n";

  }

  return bPt;
}
