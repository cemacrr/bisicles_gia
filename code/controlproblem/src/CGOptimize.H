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
		int  a_maxIter, 
		Real a_tol,
		Real a_hang,
		Real a_secantParameter,
	        Real a_secantStepMaxGrow,
	        int  a_maxSecantIter,
		Real a_secantTol)
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

  Real fm,fp, fmOld, fpOld; // these are just for retreiving the the misfit norm and penalty norm, 
              // and play no part in the method
  a_F.computeObjectiveAndGradient(fm,fp,r,x,false);
  a_F.scale(r,-1.0);
  a_F.preCond(s,r);
  a_F.assign(d,s);
  Real deltaNew = a_F.dotProduct(r,s);
  Real deltaZero = deltaNew;
  int iter = 0; int k = 0;
  while (   (iter < a_maxIter) 
	    && (deltaNew > sqr(a_tol) * deltaZero) 
	    && (fm > 1.0e-16) 
	    && ( (iter < 1) || (fm/fmOld < a_hang) ) )
    {
      
      pout() << "CGOptimize iteration " << iter 
	     << " ||f'(x)||^2 = " << deltaNew 
	     << " ||fm(x)||^2 = " << fm 
	     << " ||fp(x)||^2 = " << fp; 
      if (iter > 0)
	{
	  pout() << " ||fm(x)||^2/||fm_old(x)||^2 = " << fm/fmOld;
	}
      pout() << std::endl;

      fmOld = fm;
      fpOld = fp;

      //secant line search
      Real deltaD = a_F.dotProduct(r,d);
      Real eta = -deltaD;
      Real alpha = -a_secantParameter;
      a_F.assign(y,x);
      a_F.incr(y,d,a_secantParameter);
      Real fms,fps;
      a_F.computeObjectiveAndGradient(fms,fps,p,y,true);
      Real etaPrev = a_F.dotProduct(p,d);
      pout() << " ... initial secant ||f'(x+sd).f'(x)||^2 = " 
	     << -etaPrev << " s = " << a_secantParameter << " f(x+sd) = " << fms << std::endl;
      

      {
	int n = 1;
	while ( (fms + fps < fm + fp) && (-etaPrev > -eta))
	  {
	    //we have a problem, because although we have a descent direction in the  objective
	    //function f , we have an ascent direction for f'(x+sd).f'(x). For now, just walk along d till
	    // we do have descent in both, then start the secant method offset by alpha + a_secantParameter
	    pout() << " d is a descent direction for f(x+sd) but not for f'(x+sd).f'(x)  " 
		   << " so walking the line (remind SLC to look up a better method) " << std::endl;
	    alpha -= a_secantParameter;
	    a_F.incr(y,d,a_secantParameter);
	    a_F.computeObjectiveAndGradient(fms,fps,p,y,true);
	    eta = etaPrev;
	    etaPrev = a_F.dotProduct(p,d);
	    pout() << " ... initial secant f'(x+ sd).f'(x) = " 
		   << -etaPrev << " s = " << alpha << " f(x+ sd) = " << fms + fps << std::endl;
	  }
      }
      
      int j = 0;
      Real q = -(alpha + a_secantParameter); // normally alpha = - a_secantParameter, but sometimes we needed to walk the line
      if (q > 0.0)
	a_F.incr(x,d,q);
      
      alpha = -a_secantParameter;
       
      do {
	
	if (j > 0)
	  {
	    a_F.computeObjectiveAndGradient(fm,fp,r,x,true);
	    eta = a_F.dotProduct(r,d);
	  }
	
	pout() <<  " ... secant iteration j = " << j << " f'(x[0]+qd).f'(x[0]) = " << -eta
	       << " a = " << alpha << " q =  " << q << std::endl;
	
	Real absAlphaPrev = std::abs(alpha);
	alpha = alpha * eta / (etaPrev-eta);
	
	//limit the rate at which alpha grows
	if ( (a_secantStepMaxGrow > 0.0) && (std::abs(alpha) > a_secantStepMaxGrow * absAlphaPrev) )
	  {
	    alpha = ((alpha > 0.0)?(1.0):(-1.0)) * a_secantStepMaxGrow * absAlphaPrev;
	  }
	
	q += alpha;
	pout() << " ... testing a = " << alpha << "  q = " << q  << std::endl;
	
	a_F.incr(x,d,alpha);
	etaPrev = eta;
	j++;
      } while (j < a_maxSecantIter 
	       && sqr(alpha) * deltaD > sqr(a_secantTol));
      
     
      a_F.computeObjectiveAndGradient(fm,fp,r,x,false);
      a_F.scale(r,-1.0);
      Real deltaOld = deltaNew;
      Real deltaMid = a_F.dotProduct(r,s);
      a_F.preCond(s,r);
      deltaNew = a_F.dotProduct(r,s);
      Real beta = (deltaNew - deltaMid)/deltaOld;
      k++;
      if (k == a_F.nDoF(x) || beta <= 0.0 )
	{
	  pout() << "CGOptimize restart k = " << k 
		 << " beta = " << beta << std::endl;
	  
	  a_F.assign(d,s);
	  a_F.restart();
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
	 << " ||f'(x)||^2 = " << deltaNew 
	 << " ||fm(x)||^2 = " << fm 
	 << " ||fp(x)||^2 = " << fp << std::endl;
  a_F.free(r);
  a_F.free(s);
  a_F.free(d);
  a_F.free(y);
  a_F.free(p);


}


#include "NamespaceFooter.H"
#endif /*_CGOPTIM_H_*/
