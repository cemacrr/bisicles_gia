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
// MultiLevelDataIBC.cpp
// ============
//
// PhysIBC-derived class which stores initial topography and thickness data
// and imposes either periodic or reflection boundary conditions


#include "MultiLevelDataIBC.H"
#include "ReflectGhostCells.H"
#include "FillFromReference.H"
#include "NamespaceHeader.H"

MultiLevelDataIBC::MultiLevelDataIBC
(const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_thck,
 const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_topg, 
 const RealVect& a_dxCrse, const Vector<int> & a_refRatio)
  :m_isBCsetUp(false)
{
  MayDay::Error("MultiLevelDataIBC::MultiLevelDataIBC not implemented yet");
}

MultiLevelDataIBC::~MultiLevelDataIBC()
{
  
}

void MultiLevelDataIBC::define(const ProblemDomain& a_domain,
			       const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}


void  reflectionBC_MLDIBC(FArrayBox& a_vel,
			  const Box& a_valid,
			  const ProblemDomain& a_domain,
			  Real a_dx,
			  bool a_homogeneous)
{
  
  const IntVect ghostVect = IntVect::Unit;
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_vel, a_domain,ghostVect, dir, Side::Lo);
	  ReflectGhostCells(a_vel, a_domain,ghostVect, dir, Side::Hi);
	  
	  Box ghostBoxLo = adjCellBox(a_valid, dir, Side::Lo, 1);
              
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      ghostBoxLo &= a_vel.box();
	      a_vel.mult(-1.0,ghostBoxLo,dir,1);
		  
	    }
	  Box ghostBoxHi = adjCellBox(a_valid, dir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      ghostBoxHi &= a_vel.box();
	      a_vel.mult(-1.0,ghostBoxHi,dir,1);
	    }

	}
    }
}


BCHolder MultiLevelDataIBC::velocitySolveBC()
{
  if (!m_isBCsetUp)
    {
      m_velBCs = reflectionBC_MLDIBC;
      m_isBCsetUp = true;
    }
  return m_velBCs;
}

void MultiLevelDataIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
					 const RealVect& a_dx,
					 const RealVect& a_domainSize,
					 const Real& a_time, 
					 const LevelSigmaCS* a_crseCoords,
					 const int a_refRatio)
{
  MayDay::Error("MultiLevelDataIBC::initializeIceGeometry not implemented yet");
}

void MultiLevelDataIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				     const RealVect& a_dx,
				     const RealVect& a_domainSize,
				     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio)
{
  MayDay::Error("MultiLevelDataIBC::regridIceGeometry not implemented yet");
}

IceThicknessIBC* MultiLevelDataIBC::new_thicknessIBC()
{
  MayDay::Error("MultiLevelDataIBC::new_thicknessIBC not implemented yet");
  MultiLevelDataIBC* retval = NULL;
  return static_cast<IceThicknessIBC*>(retval);
}


/// set non-periodic ghost cells for surface height z_s. 
void MultiLevelDataIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
				       LevelSigmaCS& a_coords,
				       const ProblemDomain& a_domain,
				       const RealVect& a_dx,
				       Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Hi);
	}
    }
}

/// set non-periodic ghost cells for thickness & topography
void MultiLevelDataIBC::setGeometryBCs(LevelSigmaCS& a_coords,
				  const ProblemDomain& a_domain,
				  const RealVect& a_dx,
				  Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir))){
	ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Lo);
	ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Hi);
	ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Lo);
	ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Hi);
      }
    }
}

void MultiLevelDataIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("MultiLevelDataIBC::initialize not implemented");
}

// set boundary fluxes
void MultiLevelDataIBC::primBC(FArrayBox&            a_WGdnv,
			   const FArrayBox&      a_Wextrap,
			   const FArrayBox&      a_W,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const Real&           a_time)
{
 if (!m_domain.isPeriodic(a_dir))
   {
     int lohisign;
     Box tmp = a_WGdnv.box();
     // Determine which side and thus shifting directions
     lohisign = sign(a_side);
     tmp.shiftHalf(a_dir,lohisign);
     // (DFM - 5/28/10) this little dance with the ghostBox is a bit 
     // of a kluge to handle the case where a_WGdnv has more than one layer   
     // of ghosting, in which case just testing agains tmp isn't 
     // sufficient to determine whether you're up against the domain edge
     Box ghostBox = (a_side == Side::Lo)?
       adjCellLo(m_domain.domainBox(),a_dir, 1):
       adjCellHi(m_domain.domainBox(),a_dir, 1);
     ghostBox &= tmp;
     // Is there a domain boundary next to this grid
     if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
       {
	 tmp &= m_domain;
	 Box boundaryBox = (a_side == Side::Lo)?
	   bdryLo(tmp,a_dir):bdryHi(tmp,a_dir);
	 BoxIterator bit(boundaryBox);
	 for (bit.begin(); bit.ok(); ++bit){
	   const IntVect& i = bit();
	   a_WGdnv(i,0) = std::max(0.0,a_Wextrap(i,0));
	 }
       }
   }
}
#include "NamespaceFooter.H"
