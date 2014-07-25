#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
// VelocityBC.H
// ============
// Collection of BCFunction classes

#ifndef _VELOCITYBC_H_
#define _VELOCITYBC_H_

#include "BCFunc.H"
#include "ReflectGhostCells.H"
#include "ExtrapGhostCells.H"
#include "NamespaceHeader.H"

/// Apply no-slip conditions along one boundary
/**
   Normal velocity component = 0
   Tangent velocity component = 0
*/
class NoSlipBCFunction : public BCFunction
{

public:

  NoSlipBCFunction(int a_dir, Side::LoHiSide a_side)
  {
    m_dir = a_dir;
    m_side = a_side;
  }

  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    const IntVect ghostVect = IntVect::Unit;

    if (!(a_domain.isPeriodic(m_dir)))
      {
	ReflectGhostCells(a_vel, a_domain,ghostVect, m_dir, m_side);
	Box ghostBox = adjCellBox(a_valid, m_dir, m_side, 1);
	if(!a_domain.domainBox().contains(ghostBox))
	  {
	    ghostBox &= a_vel.box();
	    a_vel.mult(-1.0,ghostBox,0,SpaceDim);
	  }
      }
  }

private:
  int m_dir;
  Side::LoHiSide m_side;
  NoSlipBCFunction(){;}
  NoSlipBCFunction(const NoSlipBCFunction&)
  {
    MayDay::Error(" NoSlipBCFunction(const NoSlipBCFunction&) not implemented");
  }

  NoSlipBCFunction& operator=(const NoSlipBCFunction&)
  {
    MayDay::Error(" NoSlipBCFunction& operator= (const NoSlipBCFunction&) not implemented");
    return *this;
  }
};

/// Apply ice divide conditions along one boundary
/**
   Normal velocity component = 0
   Tangent velocity component normal gradient = 0
*/
class IceDivideBCFunction : public BCFunction
{

public:

  IceDivideBCFunction(int a_dir, Side::LoHiSide a_side)
  {
    m_dir = a_dir;
    m_side = a_side;
  }

  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    const IntVect ghostVect = IntVect::Unit;

    if (!(a_domain.isPeriodic(m_dir)))
      {
	ReflectGhostCells(a_vel, a_domain,ghostVect, m_dir, m_side);
	Box ghostBox = adjCellBox(a_valid, m_dir, m_side, 1);
	if(!a_domain.domainBox().contains(ghostBox))
	  {
	    ghostBox &= a_vel.box();
	    a_vel.mult(-1.0,ghostBox,m_dir,1);
	  }
      }
  }

private:
  int m_dir;
  Side::LoHiSide m_side;
  IceDivideBCFunction(){;}
  IceDivideBCFunction(const IceDivideBCFunction&)
  {
    MayDay::Error(" IceDivideBCFunction(const IceDivideBCFunction&) not implemented");
  }

  IceDivideBCFunction& operator=(const IceDivideBCFunction&)
  {
    MayDay::Error(" IceDivideBCFunction& operator= (const IceDivideBCFunction&) not implemented");
    return *this;
  }
};



/// Apply a different BCFunction along each boundary 
class PerBoundaryBCFunction : public BCFunction
{
  Vector<RefCountedPtr<BCFunction> > m_loBC;
  Vector<RefCountedPtr<BCFunction> > m_hiBC;

public:

  PerBoundaryBCFunction(Vector<RefCountedPtr<BCFunction> > a_loBC,
			Vector<RefCountedPtr<BCFunction> > a_hiBC)
    :m_loBC(a_loBC), m_hiBC(a_hiBC)
  {
    CH_assert(m_loBC.size() == SpaceDim);
    CH_assert(m_hiBC.size() == SpaceDim);
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	CH_assert(m_loBC[dir] != NULL);
	CH_assert(m_hiBC[dir] != NULL);
      }
  }

  
  PerBoundaryBCFunction()
  {
  }

  PerBoundaryBCFunction(const PerBoundaryBCFunction&)
  {
    MayDay::Error(" PerBoundaryBCFunction(const PerBoundaryBCFunction&) not implemented");
  }

  PerBoundaryBCFunction& operator=(const PerBoundaryBCFunction&)
  {
    MayDay::Error(" PerBoundaryBCFunction& operator= (const PerBoundaryBCFunction&) not implemented");
    return *this;
  }


  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	(*m_loBC[dir])(a_vel, a_valid, a_domain, a_dx, a_homogeneous);
	(*m_hiBC[dir])(a_vel, a_valid, a_domain, a_dx, a_homogeneous);
      }
  }

private:
  int m_dir;
  Side::LoHiSide m_side;

  
};

#include "NamespaceFooter.H"

#endif
