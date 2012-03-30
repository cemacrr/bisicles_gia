#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// SLC, Nov 9, 2011

#ifndef _CGOPTIMIZE_H_
#define _CGOPTIMIZE_H_


#include "parstream.H"
#include "REAL.H"
#include <cmath>
#include "NamespaceHeader.H"
#include "CH_assert.H"

/**
   Nonlinear Conjugate-gradient method, used to mimimise F(x) w.r.t x 
*/
inline Real sqr(const Real& a)
{
  return a*a;
}

template <class O, class X>
void CGOptimize(O& a_F, X& a_x, 
		const int  a_maxIter, 
		const Real a_tol,
		const Real a_secantParameter,
		const int  a_maxSecantIter,
		const Real a_secantTol)
{

  //Preconditioned Nonlinear Conjugate Gradients 
  //with Secant line search and Polak-Ribiere updates
  //from J. R. Shewchuk, 1994 
  //"An Introduction to the Conjugate Gradient Method
  // Without the Agonizing Pain"

  X& x = a_x;
  X r,s,d,y,p;
  a_F.create(r,x);
  a_F.create(s,x);
  a_F.create(d,x);
  a_F.create(y,x);
  a_F.create(p,x);

  Real f;
  a_F.computeObjectiveAndGradient(f,r,x,false);
  a_F.scale(r,-1.0);
  a_F.preCond(s,r);
  a_F.assign(d,s);
  Real deltaNew = a_F.dotProduct(r,s);
  Real deltaZero = deltaNew;
  int iter = 0; int k = 0;
  while (iter < a_maxIter && deltaNew > sqr(a_tol) * deltaZero)
    {
      
      pout() << "CGOptimize iteration " << iter << " ||f'(x)||^2 = " 
	     << deltaNew << " f(x) = " << f << std::endl;


      //secant line search
      Real deltaD = a_F.dotProduct(r,d);
      Real eta = -deltaD;
      Real alpha = -a_secantParameter;
      a_F.assign(y,x);
      a_F.incr(y,d,a_secantParameter);
      Real fs;
      a_F.computeObjectiveAndGradient(fs,p,y,true);
      Real etaPrev = a_F.dotProduct(p,d);
      pout() << " ... initial secant ||f'(x+sd).f'(x)||^2 = " 
	     << -etaPrev << " s = " << a_secantParameter << " f(x+sd) = " << fs << std::endl;
      
      int j = 0;
      Real q = 0.0;
      
      //lowest known positive value of f(x+qd).f(x)
      bool lposKnown = false;
      Real lposQ = 1e+100,lposEta = -1.0e+100;
      //highest known negative value of f(x+qd).f(x)
      bool unegKnown = false;
      Real unegQ = 1e+100,unegEta = 1.0e+100;
	
      do {
	if (j > 0)
	  {
	    a_F.computeObjectiveAndGradient(f,r,x,true);
	    eta = a_F.dotProduct(r,d);
	  }
	if (eta > 0)
	  {
	    if (!lposKnown || (lposKnown && eta < lposEta))
	      {
		lposEta = eta;
		lposQ = q;
		lposKnown = true;
	      }
	  }
	if (eta < 0)
	  {
	    if (!unegKnown || (unegKnown && eta > unegEta))
	      {
		unegEta = eta;
		unegQ = q;
		unegKnown = true;
	      }
	  }


	pout() <<  " ... secant iteration j = " << j << " f'(x[0]+qd).f'(x[0]) = " << -eta
	       << " a = " << alpha << " q =  " << q <<  " f(x[j]+ad) = " << f << std::endl;
	alpha = alpha * eta / (etaPrev-eta);
	q += alpha;
	pout() << " ... testing a = " << alpha << "  q = " << q  << std::endl;
	
	if (lposKnown && unegKnown)
	  {
	    if (q > lposQ && q > unegQ)
	      {
		
		Real newq = std::max(lposQ,unegQ);
		alpha -= (q - newq);
		q = newq;
		pout() << " (q > lposQ = " << lposQ << " &&  q  >  unegQ = " << unegQ 
		       << " ) so setting q = " << q << std::endl;
	      }
	    else if (q < lposQ && q < unegQ)
	      {
		Real newq = std::min(lposQ,unegQ);
		alpha += (newq - q);
		q = newq;
		pout() << " (q < lposQ = " << lposQ << " &&  q  <  unegQ = " << unegQ 
		       << " ) so setting q = " << q << std::endl;
	      }
	  }

	a_F.incr(x,d,alpha);
	etaPrev = eta;
	j++;
      } while (j < a_maxSecantIter 
	       && sqr(alpha) * deltaD > sqr(a_secantTol));
    
      
      a_F.computeObjectiveAndGradient(f,r,x,false);
      a_F.scale(r,-1.0);
      Real deltaOld = deltaNew;
      Real deltaMid = a_F.dotProduct(r,s);
      a_F.preCond(s,r);
      deltaNew = a_F.dotProduct(r,s);
      Real beta = (deltaNew - deltaMid)/deltaOld;
      k++;
      if (k == a_F.nDoF(x) || beta <= 0.0)
	{
	  pout() << "CGOptimize restart k = " << k 
		 << " beta = " << beta << std::endl;
	  a_F.assign(d,s);
	  k = 0;
	}
      else
	{
	  a_F.scale(d,beta);
	  a_F.incr(d,s,1.0);
	}
      iter++;
    }
  pout() << "CGOptimize iteration " << iter 
	 << " ||f'(x)||^2 = " 
	 << deltaNew << " f(x) = " << f << std::endl;
  a_F.free(r);
  a_F.free(s);
  a_F.free(d);
  a_F.free(y);
  a_F.free(p);


}


#include "NamespaceFooter.H"
#endif /*_CGOPTIM_H_*/
