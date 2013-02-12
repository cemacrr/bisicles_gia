#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "FortranInterfaceIBC.H"
#include "ParmParse.H"
#include "FillFromReference.H"
#include "ExtrapGhostCells.H"
#include "ExtrapBCF_F.H"
#include "ReflectGhostCells.H"
#include "FIBCF_F.H"
#include "NamespaceHeader.H"


void zeroBCValue_FIBC(Real* pos,
		     int* dir,
		     Side::LoHiSide* side,
		     Real* a_values)
{
  a_values[0]=0.0;
  a_values[1]=0.0;
}



void iceNeumannBC_FIBC(FArrayBox& a_state,
		       const Box& a_valid,
		       const ProblemDomain& a_domain,
		       Real a_dx,
		       bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
		  Interval NeumInterval(0,SpaceDim-1);
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue_FIBC,
                         dir,
                         Side::Lo,
			 NeumInterval);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
		  Interval NeumInterval(0,SpaceDim-1);
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue_FIBC,
                         dir,
                         Side::Hi, 
			 NeumInterval);
                }

            } // end if is not periodic in ith direction
        }
    }
}


void iceDirichletBC_FIBC(FArrayBox& a_state,
			 const Box& a_valid,
			 const ProblemDomain& a_domain,
			 Real a_dx,
			 bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue_FIBC,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue_FIBC,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}

//set all components of u to zero in ghost regions
void iceDivideBC_FIBC(FArrayBox& a_state,
		      const Box& a_valid,
		      const ProblemDomain& a_domain,
		      Real a_dx,
		      bool a_homogeneous)
{
if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
		  ghostBoxLo &= a_state.box();
		  a_state.setVal(0.0,  ghostBoxLo, 0, a_state.nComp());
		  
                }
	      Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
		  ghostBoxHi &= a_state.box();
		  a_state.setVal(0.0,  ghostBoxHi, 0, a_state.nComp());
                }

            } // end if is not periodic in ith direction
        }
    }

}


//ice stream parallel to x : homogeneous Neumann conditions
//along the y faces, and homogeneous Dirichlet conditions
//along the x faces. 
void iceStreamXBC_FIBC(FArrayBox& a_state,
		      const Box& a_valid,
		      const ProblemDomain& a_domain,
		      Real a_dx,
		      bool a_homogeneous)
{

  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      
      int xDir = 0;
      if (!a_domain.isPeriodic(xDir))
	{
	  Box ghostBoxLo = adjCellBox(valid, xDir, Side::Lo, 1);
	  
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      ghostBoxLo &= a_state.box();
	      a_state.setVal(0.0,  ghostBoxLo, 0, a_state.nComp());
	      
	    }
	  Box ghostBoxHi = adjCellBox(valid, xDir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      ghostBoxHi &= a_state.box();
	      a_state.setVal(0.0,  ghostBoxHi, 0, a_state.nComp());
	    }
	}

      if (SpaceDim > 1)
	{
	  int yDir = 1;
	  Box ghostBoxLo = adjCellBox(valid, yDir, Side::Lo, 1);
	  Box ghostBoxHi = adjCellBox(valid, yDir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      //Real bcVal = 0.0;
	      Interval NeumInterval(0,SpaceDim-1);
	      NeumBC(a_state,
		     valid,
		     a_dx,
		     a_homogeneous,
		     zeroBCValue_FIBC,
		     yDir,
		     Side::Lo,
		     NeumInterval);
	    }
	  
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      Interval NeumInterval(0,SpaceDim-1);
	      NeumBC(a_state,
		     valid,
		     a_dx,
		     a_homogeneous,
		     zeroBCValue_FIBC,
		     yDir,
		     Side::Hi, 
		     NeumInterval);
	    }
	  
	}
    }
}




// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
FortranInterfaceIBC::FortranInterfaceIBC() : m_boundaryThickness(0.0)
{
  m_isBCsetUp = false;
  m_isDefined = false;
  m_gridsSet = false;
  // set default to be true
  m_verbose = true;
  m_thicknessGhost = IntVect::Zero;
  m_topographyGhost = IntVect::Zero;
  // default is an empty vector here 
  // (don't set thickness to zero anywhere)
  m_thicknessClearRegions = Vector<Box>();
  
}

FortranInterfaceIBC::~FortranInterfaceIBC()
{
}


/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
FortranInterfaceIBC::define(const ProblemDomain& a_domain,
			    const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
  m_isDefined = true;
}

// define cell centred a_fab, directly with cell centered data
// a_data_ptr if !a_nodal, or by averaging nodal data a_data_ptr.
// if the data is nodal, also compute the change in value across the box
void FortranInterfaceIBC::setFAB(Real* a_data_ptr,
				 const int* a_dimInfo,
                                 const int* a_boxlo, const int* a_boxhi, 
				 const Real* a_dew, const Real* a_dns,
				 const IntVect& a_offset,
				 const IntVect& a_nGhost,
				 FArrayBox & a_fab,
				 const bool a_nodal)
{
  IntVect loVect(D_DECL(a_boxlo[0],a_boxlo[1], a_boxlo[2]-1));
  IntVect hiVect(D_DECL(a_boxhi[0],a_boxhi[1], a_boxhi[2]-1));
  Box fabBox(loVect, hiVect);
  fabBox.shift(-a_nGhost);
  fabBox.shift(a_offset);
  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      Box gridBox(fabBox);
      gridBox.grow(-a_nGhost);
      setGrids(gridBox);
    }
  if (a_nodal)
    {
      a_fab.define(fabBox, 1);
      Box nodeBox(IntVect::Zero, hiVect + IntVect::Unit);
      nodeBox.shift(-a_nGhost);
      FArrayBox nodeFAB;
      nodeFAB.define(nodeBox,1,a_data_ptr);
      FORT_NODETOCELL(CHF_CONST_FRA1(nodeFAB,0),
		      CHF_FRA1(a_fab,0),
		      CHF_BOX(fabBox));
    }
  else 
    {
      a_fab.define(fabBox, 1, a_data_ptr);
    }
}


void
FortranInterfaceIBC::setThickness(Real* a_data_ptr,
				  const int* a_dimInfo,
                                  const int* a_boxlo, const int* a_boxhi, 
                                  const Real* a_dew, const Real* a_dns,
				  const IntVect& a_offset,
                                  const IntVect& a_nGhost,
				  const bool a_nodal)
{

  m_thicknessGhost = a_nGhost;

  // dimInfo is (SPACEDIM, nz, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering
  //cout << "a_dimonfo" << a_dimInfo[0] << a_dimInfo[1] << endl;  

  setFAB(a_data_ptr, a_dimInfo,a_boxlo, a_boxhi,
         a_dew,a_dns,a_offset,a_nGhost,
	 m_inputThickness,a_nodal);
  m_inputThicknessDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // now define LevelData and copy from FAB->LevelData 
  // (at some point will likely change this to be an  aliased 
  // constructor for the LevelData, but this should be fine for now....
  RefCountedPtr<LevelData<FArrayBox> > localLDFPtr(new LevelData<FArrayBox>(m_grids, 1, m_thicknessGhost) );
  m_inputThicknessLDF = localLDFPtr;
  // fundamental assumption that there is only one box per processor here
  DataIterator dit = m_grids.dataIterator();
  dit.begin();
  Box copyBox = (*m_inputThicknessLDF)[dit].box();
  copyBox &= m_inputThickness.box();
  (*m_inputThicknessLDF)[dit].copy(m_inputThickness, copyBox);


  // if necessary, set thickness to zero where needed
  if (m_thicknessClearRegions.size() > 0) 
    {
      for (int i=0; i<m_thicknessClearRegions.size(); i++)
        {
          Box region = m_thicknessClearRegions[i];
          // adjust for ghosting
          region.shift(-a_nGhost);
          region &= m_inputThickness.box();
          if (!region.isEmpty())
            {
              // do this in distributed copy, rather than in the original
              (*m_inputThicknessLDF)[dit].setVal(0.0, region, 0, 1);
            }
        } // end loop over thickness clear regions
    }

}

void
FortranInterfaceIBC::setTopography(Real* a_data_ptr,
				   const int* a_dimInfo,
                                   const int* a_boxlo, const int* a_boxhi, 
				   const Real* a_dew, const Real* a_dns,
				   const IntVect& a_offset,
                                   const IntVect& a_nGhost, 
				   const bool a_nodal)
{

  // dimInfo is (SPACEDIM, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering

  m_topographyGhost = a_nGhost;
  setFAB(a_data_ptr, a_dimInfo,a_boxlo, a_boxhi, a_dew,a_dns,
         a_offset, a_nGhost,m_inputTopography,a_nodal);
  m_inputTopographyDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // now define LevelData and copy from FAB->LevelData 
  // (at some point will likely change this to be an aliased 
  // constructor for the LevelData, but this should be fine for now....
  RefCountedPtr<LevelData<FArrayBox> > localLDFPtr(new LevelData<FArrayBox>(m_grids, 1, m_topographyGhost));
  m_inputTopographyLDF = localLDFPtr;
  // fundamental assumption that there is only one box per processor here
  DataIterator dit = m_grids.dataIterator();
  dit.begin();
  Box copyBox = (*m_inputTopographyLDF)[dit].box();
  copyBox &= m_inputTopography.box();
  (*m_inputTopographyLDF)[dit].copy(m_inputTopography, copyBox);



}


  /// regions where we artificially set thickness to zero
  /** this is done in setThickness, for lack of a better place, so 
      this needs to be set before setThickness is called. 
      a_clearRegions defines logically-rectangular regions where 
      the thickness is artificially set to zero. These regions are 
      defined relative to the original input data (i.e. before any shifting 
      due to ghost cells)
  */
void
FortranInterfaceIBC::setThicknessClearRegions(const Vector<Box>& a_clearRegions)
{
  m_thicknessClearRegions = a_clearRegions;
}


/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
IceThicknessIBC* 
FortranInterfaceIBC::new_thicknessIBC()
{
  FortranInterfaceIBC* retval = new FortranInterfaceIBC();

  retval->m_grids = m_grids;
  retval->m_gridsSet = m_gridsSet;

  // only do this if these guys are actually defined
  if (!m_inputThickness.box().isEmpty())
    {
      //retval->m_inputThickness.define(m_inputThickness.interval(), m_inputThickness);
      retval->m_inputThickness.define(m_inputThickness.box(), m_inputThickness.nComp());
      retval->m_inputThickness.copy(m_inputThickness);

      retval->m_inputThicknessLDF = m_inputThicknessLDF;
    }
  retval->m_thicknessGhost = m_thicknessGhost;
  retval->m_inputThicknessDx = m_inputThicknessDx;

  if (!m_inputTopography.box().isEmpty())
    {
      //retval->m_inputTopography.define(m_inputTopography.interval(), m_inputTopography);
      retval->m_inputTopography.define(m_inputTopography.box(), m_inputTopography.nComp());
      retval->m_inputTopography.copy(m_inputTopography);

      retval->m_inputTopographyLDF = m_inputTopographyLDF;
    }
  retval->m_inputTopographyDx = m_inputTopographyDx;
  retval->m_topographyGhost = m_topographyGhost;

  retval->m_verbose = m_verbose;

  return static_cast<IceThicknessIBC*>(retval);
}

/// Set up initial conditions
/**
 */
void
FortranInterfaceIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("FortranInterfaceIBC::initialize not implemented");
}

/// Set boundary fluxes
/**
 */
void 
FortranInterfaceIBC::primBC(FArrayBox&            a_WGdnv,
                          const FArrayBox&      a_Wextrap,
                          const FArrayBox&      a_W,
                          const int&            a_dir,
                          const Side::LoHiSide& a_side,
                          const Real&           a_time)
{
  // do nothing in periodic case
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
	    //a_WGdnv(i,0) = 0.0;
	  }
	}
    }
  
}

/// Set boundary slopes
/**
   The boundary slopes in a_dW are already set to one sided difference
   approximations.  If this function doesn't change them they will be
   used for the slopes at the boundaries.
*/

void 
FortranInterfaceIBC::setBdrySlopes(FArrayBox&       a_dW,
                                 const FArrayBox& a_W,
                                 const int&       a_dir,
                                 const Real&      a_time)
{
  // one-sided differences sounds fine with me, so do nothing...
}

/// Adjust boundary fluxes to account for artificial viscosity
/**
 */
void 
FortranInterfaceIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here -- if we wind up here, need to
  // give it some thought
  MayDay::Error("FortranInterfaceIBC::artViscBC not implemented");
}


/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
BCHolder
FortranInterfaceIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs;
}


/// set non-periodic ghost cells for surface height z_s. 
void
FortranInterfaceIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                         LevelSigmaCS& a_coords,
                                         const ProblemDomain& a_domain,
                                         const RealVect& a_dx,
                                         Real a_time, Real a_dt)
{
  if (m_extrapBoundary)
    {
      ExtrapGhostCells(a_zSurface, a_domain);
    } 
  else
    {
      //slc : i think reflection is a better default : works for isolated
      //      islands (Antarctica, Greenland) and divides, 
      //      where grad(s) = grad(H) = grad(topo) = 0 (Pine Island).
      //      (and in that case, this ought to be unnecessary)
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Hi);
	}
    }
}

/// set non-periodic ghost cells for thickness & topography
void
FortranInterfaceIBC::setGeometryBCs(LevelSigmaCS& a_coords,
				    const ProblemDomain& a_domain,
				    const RealVect& a_dx,
				    Real a_time, Real a_dt)
{
  if (m_extrapBoundary)
    {
      ExtrapGhostCells(a_coords.getH(), a_domain);
      ExtrapGhostCells(a_coords.getTopography(), a_domain);
    } 
  else
    {
      //slc : i think reflection is a better default : works for isolated
      //      islands (Antarctica, Greenland) and divides, 
      //      where grad(s) = grad(H) = grad(topo) = 0 (Pine Island).
      //      (and in that case, this ought to be unnecessary)
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Hi);
	  ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Hi);
	}
    }

}


/// set grids using Boxes passed in from Glimmer-CISM
/** creates a DisjointBoxLayout using the grid boxes and 
    processor distribution used by CISM. */
void 
FortranInterfaceIBC::setGrids(const Box& a_gridBox)
{

  if (m_verbose)
    {
      pout() << "in setGrids, box = " << a_gridBox << endl;
    }

  CH_assert(m_isDefined);

  const int numBox = numProc();
  Vector<Box> boxes(numBox);

  if (numBox == 1) boxes[0] = a_gridBox;

  Box* nonConstBox = const_cast<Box*>(&a_gridBox);
#ifdef CH_MPI
  int boxSize = sizeof(Box);

  pout() << "entering allGather -- sending " << *nonConstBox << endl;
  pout () << "numBox = " << numBox << ", boxSize = " << boxSize << endl;
  
  MPI_Allgather(nonConstBox,  boxSize,  MPI_BYTE, &(boxes[0]), 
                boxSize , MPI_BYTE , Chombo_MPI::comm);

  pout () << "after allGather" << endl;
  pout () << "nonConstbox = " << *nonConstBox << endl;


#endif

  pout () << "numBoxes = " << boxes.size() << endl;
  for (int i=0; i< boxes.size(); i++)
    {
      pout () << "box " << i << ": " << boxes[i] << endl;
    }

  Vector<int> procAssign(boxes.size());
  for (int i=0; i< boxes.size(); i++)
    {
      procAssign[i] = i;
    }
  
  //if (m_verbose)
    {
      pout() << "processor " << procID() << ": grids: " << endl;
      for (int i=0; i<boxes.size(); i++)
        {
          pout () << procAssign[i] << ": " << boxes[i] << endl;
        }
    }

    pout () << " before DBL define" << endl;

  // define DisjointBoxLayout
  m_grids.define(boxes, procAssign, m_domain);

  pout () << "after DBL define" << endl;

  m_gridsSet = true;
}

/// utility function to fill in holes in topography
/** looks for isolated values of holeVal and replaces then
    with average of neighbors */
void
FortranInterfaceIBC::fillTopographyHoles(Real a_holeVal)
{

  if (m_verbose) 
    {
      pout() << "FortranInterfaceIBC::fillTopographyHoles" << endl;
    }

  CH_assert(!m_inputTopography.box().isEmpty());
  
  Real compareThreshold = 1.0e-10;
  // neighborCutoff is number of neighbors which are not the 
  // holeValue required before we do averaging. (max number of 
  // neigbors in 2d is 8
  int neighborCutoff = 5;

  Box testBox = m_inputTopography.box();
  testBox.grow(-1);

  Box neighborBox(-IntVect::Unit, IntVect::Unit);
  IntVectSet neighborIVS(neighborBox);
  neighborIVS -= IntVect::Zero;

  IVSIterator neighborIt(neighborIVS);

  BoxIterator bit(testBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (abs(m_inputTopography(iv, 0) - a_holeVal) < compareThreshold)
        {
          Real sum = 0;
          int num = 0;
          for (neighborIt.begin(); neighborIt.ok(); ++neighborIt)
            {
              IntVect newIV = iv + neighborIt();
              if (abs(m_inputTopography(newIV, 0) - a_holeVal) > compareThreshold)
                {
                  sum += m_inputTopography(newIV, 0);
                  num++;
                }
            }

          if (num > neighborCutoff)
            {
              m_inputTopography(iv,0) = sum / num;
            }
        }
    }
          

#if 0
  // try this in c++ first
  FORT_FILLHOLES(CHF_FRA1(m_inputTopography,0),
                 CHF_REAL(a_holeVal),
                 CHF_REAL(compareThreshold),
                 CHF_INT(neighborCutoff),
                 CHF_BOX(testBox));
#endif
}


  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
   */
void
FortranInterfaceIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                                           const RealVect& a_dx,
                                           const RealVect& a_domainSize,
                                           const Real& a_time, 
					   const LevelSigmaCS* a_crseCoords,
					   const int a_refRatio)
{
  if (m_verbose) 
    {
      pout() << "FortranInterfaceIBC::initializeIceGeometry" << endl;
    }

  
  ParmParse geomPP("geometry");

  //CH_assert(m_isDefined);

  // assume that thickness and topography from glimmer are defined on the same 
  // mesh
  CH_assert(m_inputThickness.box() == m_inputTopography.box());
  CH_assert(m_inputThicknessDx == m_inputTopographyDx);

  // get scaling factor for ice thickness -- not sure if I really need this
  Real thicknessScale =1;
  geomPP.query("thickness_scale", thicknessScale);

  if (m_verbose)
    {
      pout() << " ...setting up thickness LevelData..." << endl;
    }

  
  // grab periodicity info from a_coords
  const DisjointBoxLayout& destGrids = a_coords.grids();
  const ProblemDomain& destDomain = destGrids.physDomain();
  if (geomPP.contains("basalSlope") 
      && (destDomain.isPeriodic(0) || destDomain.isPeriodic(1)) )
    {
      //slc : periodic domains need not have periodic topography b,
      //instead b = f(x,y) + Ax + By where A,B are constants and f is periodic.
      //for now, read A,B from the input file - periodic problems
      //are rare, so it is possibly not worth the effort of devising
      //something robust, e.g working out A and B from glimmer's model%geometry%topg, 
      //which is node-centred data. It would be easy to do that in serial, 
      //but in parallel the details would depend on the
      //interaction between glimmer's domain decomposition and ours.
      RealVect basalSlope(RealVect::Zero);
      Vector<Real> t(SpaceDim, 0.0);
      geomPP.getarr("basalSlope", t, 0, SpaceDim);
      D_TERM(basalSlope[0] = t[0];,
	     basalSlope[1] = t[1];,
	     basalSlope[2] = t[2];);
      if (m_verbose)
	{      
	  pout() << " ... adding constant background slope " 
		 << basalSlope ;
	}
      if (basalSlope.dotProduct(basalSlope) > 1.0e-10)
	{
	  MayDay::Error("periodic problems with non-zero basal slope not yet implemented");
	}
      a_coords.setBackgroundSlope(basalSlope);

      RealVect unitShift = basalSlope * a_domainSize;
      if (m_verbose)
	{      
	  pout() << " a_domainSize[0] = "  
		 << a_domainSize[0]
		 << " ... " << std::endl;
	}

      
    }
  // assume here that topography and thickness defined on same box
  // passed-in ice topography on a LevelData

#if 0
  if (m_verbose)
    {
      pout() << " ...setting up topography LevelData..." << endl;
    }
  CH_assert(m_thicknessGhost == m_topographyGhost);

  DisjointBoxLayout& topographyDBL = thicknessDBL;
  LevelData<FArrayBox> topographyLD(topographyDBL, 1);
  
  {
    DataIterator dit = topographyDBL.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
	topographyLD[dit].copy(m_inputTopography);
      }
  }
#endif



  DisjointBoxLayout ChomboGrids = a_coords.grids();
  LevelData<FArrayBox>& ChomboThickness = a_coords.getH();

  IntVect thicknessGhostVect =  ChomboThickness.ghostVect();
  IntVect sigmaGhostVect = thicknessGhostVect - IntVect::Unit;

  //LevelData<FArrayBox> ChomboTopography(ChomboGrids, 1, thicknessGhostVect);
  LevelData<FArrayBox>& ChomboTopography = a_coords.getTopography();
  Real tolerance = 1.0e-6;
  if (a_crseCoords != NULL &&
      a_refRatio * a_dx[0] <= m_inputThicknessDx[0] * (1.0 + tolerance))
    {
      // in this (common) case, interpolation from a_crseCoords is as good as it gets
      if (m_verbose)
	{
	  pout() << " ...interpolating data from coarse LevelSigmaCS with refinement ratio = " 
		 << a_refRatio << endl;
	}

      a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
    }
  else 
    {

      FillFromReference(ChomboThickness,
			(*m_inputThicknessLDF),
			a_dx,
			m_inputThicknessDx,
			m_verbose);

      FillFromReference(ChomboTopography,
			(*m_inputTopographyLDF),
			a_dx,
			m_inputTopographyDx,
			m_verbose);
      if (a_crseCoords!= NULL)
	{
	  // fill ghost cells
	  a_coords.interpFromCoarse(*a_crseCoords, a_refRatio, 
				    false, false, false);
	}

    }
  
  ChomboThickness.exchange();
  //ChomboTopography.exchange();
  //LevelSigmaCS.exchangeTopography() takes care of constant topography slopes
  a_coords.exchangeTopography();
  
  // if we have ghost cells, then use extrapolation BC's (for now, at least)
  // to fill in reasonable values for topography and thickness
  m_extrapBoundary = true;
  geomPP.query("extrap_boundary", m_extrapBoundary );
  if (m_extrapBoundary)
    {
      if (m_verbose)
	{      
	  pout() << "extrapolating boundary thickness and topography" << endl;
	}
      ExtrapGhostCells(ChomboThickness, destDomain);
      ExtrapGhostCells(ChomboTopography, destDomain);
    }
  else 
    {
      //slc : I think reflection would make a better default : works for isolated
      //      islands (Antarctica, Greenland) and divides, 
      //      where grad(s) = grad(H) = grad(topo) = 0 (Pine Island).
      if (m_verbose)
	{      
	  pout() << "reflecting boundary thickness and topography" << endl;
	}
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  ReflectGhostCells(ChomboThickness, destDomain, dir, Side::Lo);
	  ReflectGhostCells(ChomboThickness, destDomain, dir, Side::Hi);
	  ReflectGhostCells(ChomboTopography, destDomain, dir,  Side::Lo);
	  ReflectGhostCells(ChomboTopography, destDomain, dir, Side::Hi);
	}
    }
  //this can wait
  //a_coords.recomputeGeometry();
  if (m_verbose)
    {      
      pout() << "leaving FortranInterfaceIBC::initializeIceGeometry" << endl;
    }
}

void
FortranInterfaceIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				       const RealVect& a_dx,
				       const RealVect& a_domainSize,
				       const Real& a_time, 
				       const LevelSigmaCS* a_crseCoords,
				       const int a_refRatio)
{

   if (m_verbose)
    {      
      pout() << "entering FortranInterfaceIBC::regridIceGeometry" << endl;
    }

   Real tolerance = 1.0e-6;
   if (a_crseCoords != NULL &&
       a_refRatio * a_dx[0] <= m_inputThicknessDx[0] * (1.0 + tolerance))
     {
       // in this (common) case, interpolation from a_crseCoords is as good as it gets
       if (m_verbose)
	 {
	   pout() << " ...interpolating data from coarse LevelSigmaCS with refinement ratio = " 
		  << a_refRatio << endl;
	 }
       
       a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
     }
   else
     {
       LevelData<FArrayBox>& ChomboTopography = a_coords.getTopography();
       FillFromReference(ChomboTopography,
			 m_inputTopography,
			 a_dx,
			 m_inputTopographyDx,
			 m_topographyGhost,
			 m_verbose);
       if (a_crseCoords!= NULL)
	 {
	   // fill ghost cells
	   a_coords.interpFromCoarse(*a_crseCoords, a_refRatio, 
				     false, false, false);
	 }
     }
   if (m_verbose)
     {      
       pout() << "leaving FortranInterfaceIBC::regridIceGeometry" << endl;
     }
   
}


void 
FortranInterfaceIBC::setupBCs()
{
  ParmParse ppBC("bc");

  // get boundary conditions 
  Vector<int> loBCvect(SpaceDim), hiBCvect(SpaceDim);
  ppBC.getarr("lo_bc", loBCvect, 0, SpaceDim);
  ppBC.getarr("hi_bc", hiBCvect, 0, SpaceDim);

  // this is a placeholder until I can get a BCHolder to work...
  // require all directions to have the same BC for now
  CH_assert(loBCvect[0] == loBCvect[1]);
  CH_assert(hiBCvect[0] == hiBCvect[1]);
  CH_assert(loBCvect[0] == hiBCvect[0]);

  if (loBCvect[0] == 0)
    {
      m_velBCs = iceDirichletBC_FIBC;
    }
  else if (loBCvect[0] == 1)
    {
      m_velBCs = iceNeumannBC_FIBC;
    }
  else if (loBCvect[0] == 2)
    {
      m_velBCs = iceDivideBC_FIBC;
    }
  else if (loBCvect[0] == 3)
    {
      m_velBCs = iceStreamXBC_FIBC;
    }

  else
    {
      MayDay::Error("bad BC type");
    }
  m_isBCsetUp = true;
}
      
void FortranInterfaceIBC::checkOK() const
{
  pout() << "FortranInterfaceIBC::checkOK() &m_inputTopography = " 
	 << &m_inputTopography << std::endl;
  CH_assert(m_inputTopography.norm(0) < 1.0e+5);
}
#include "NamespaceFooter.H"
