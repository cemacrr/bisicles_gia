#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifdef HAVE_PYTHON

// PythonInterface.cpp
// ============
// BISICLES embedded Python interface. Includes
//
// PythonInterface, a namespace for common parts of the interface
//
// PythonInterface::PythonIBC, a PhysIBC-derived class which allows the bedrock and topography to be set by
// a user-defined python function. Reflection boundary conditions are imposed on
// each edge. 
//
// PythonInterface::PythonSurfaceFlux , a SurfaceFlux derived class which allows a thickness source to be set
// by a python function f(x,y,thickness,topography)


#include "PythonInterface.H"
#include "ReflectGhostCells.H"
#include "FillFromReference.H"
#include "CoarseAverage.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "signal.h"
#include "NamespaceHeader.H"

void PythonInterface::PythonEval(PyObject* a_pFunc, 
				 Vector<Real>& a_value, 
				 Vector<Real>& a_args)
{
  PyObject *pArgs, *pValue;
  pArgs = PyTuple_New(a_args.size());
  for (int j = 0; j < a_args.size(); ++j)
    {
      pValue = PyFloat_FromDouble(a_args[j]);
      if (!pValue)
	{
	  pout() << "PythonInterface::PythonEval failed to construct Python args" << endl;
	  CH_assert(pValue != NULL);
	  MayDay::Error("PythonInterface::PythonEval failed to construct Python args");
	}
      PyTuple_SetItem(pArgs, j, pValue);
    }

  pValue = PyObject_CallObject(a_pFunc, pArgs);
  Py_DECREF(pArgs);
  if (!pValue)
    {
      pout() << "PythonInterface::PythonEval  python call failed" << endl;
      PyErr_Print();
      CH_assert(pValue != NULL);
      MayDay::Error("PythonInterface::PythonEval python call failed");
    }
  
  if (PyFloat_CheckExact(pValue))
    {
      a_value.resize(1,0.0);
      a_value[0] =  PyFloat_AS_DOUBLE(pValue);
    }
  else
    {
      Py_DECREF(pValue);
      CH_assert(PyFloat_CheckExact(pValue));
      MayDay::Error("PythonInterface::PythonEval python return value is not a float");
    }
}

void PythonInterface::InitializePythonModule(PyObject **a_pModule,
					     const std::string& a_pyModuleName)
{
  //initialize Python without overriding SIGINT  
  PyOS_sighandler_t sigint = PyOS_getsig(SIGINT);
  Py_Initialize();
  sigint =  PyOS_setsig(SIGINT, sigint);
  
  PyObject* pName = PyString_FromString(a_pyModuleName.c_str());
  *a_pModule =  PyImport_Import(pName);

  if (*a_pModule == NULL)
    {
      PyErr_Print();
      pout() << "failed to import Python module " << a_pyModuleName << endl;
      CH_assert(a_pModule != NULL);
      MayDay::Error("failed to import Python module ");
    }

}

void PythonInterface::InitializePythonFunction(  PyObject **a_pFunc,
						 PyObject*  a_pModule,
						 const std::string& a_pyFuncName)
{

 
  pout() << "PythonInterface::InitializePythonFunction initializing " 
	 << a_pyFuncName << endl;
  
  if ( a_pModule == NULL)
    {
      pout() << "null python module " << a_pyFuncName << endl;
      CH_assert(a_pModule != NULL);
      MayDay::Error("null python module");
    }

  char* t = const_cast<char*>(a_pyFuncName.c_str());
  *a_pFunc = PyObject_GetAttrString(a_pModule,t);
  
  if (a_pFunc == NULL)
    {
      pout() << "failed to initialize Python function " << a_pyFuncName << endl;
      CH_assert(*a_pFunc != NULL);
      MayDay::Error("failed to initialize Python function");
    }

  if (!PyCallable_Check(*a_pFunc))
    {
      pout() << "Python function " << a_pyFuncName << " not callable " <<  endl;
      CH_assert(PyCallable_Check(*a_pFunc));
      MayDay::Error("Python function not callable");
    }


}


PythonInterface::PythonIBC::PythonIBC(const std::string& a_pyModuleName,
				      const std::string& a_pyFuncThckName,
				      const std::string& a_pyFuncTopgName)
  :m_isBCsetUp(false),m_verbose(true)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFuncThck, m_pModule,  a_pyFuncThckName);
  InitializePythonFunction(&m_pFuncTopg, m_pModule,  a_pyFuncTopgName);
}


PythonInterface::PythonIBC::PythonIBC()
  :m_isBCsetUp(false),m_verbose(true)
{

}

PythonInterface::PythonIBC::~PythonIBC()
{

}

void PythonInterface::PythonIBC::define(const ProblemDomain& a_domain,
		       const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
  //m_dx = a_dx * RealVect::Unit;
}

void PythonInterface::PythonIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("PythonIBC::initialize not implemented");
}

// set boundary fluxes
void PythonInterface::PythonIBC::primBC(FArrayBox&            a_WGdnv,
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

void  PythonInterface::PythonIBC::setBdrySlopes(FArrayBox&       a_dW,
				  const FArrayBox& a_W,
				  const int&       a_dir,
				  const Real&      a_time)
{
  // one-sided differences sounds fine with me, so do nothing...
}


void PythonInterface::PythonIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here
  MayDay::Error("PythonIBC::artViscBC not implemented");
}

/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
BCHolder PythonInterface::PythonIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs; 
}

void  iceDivideBC_PyBC(FArrayBox& a_vel,
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


void  PythonInterface::PythonIBC::setupBCs()
{
  m_velBCs = iceDivideBC_PyBC;
  m_isBCsetUp = true;
}

/// set non-periodic ghost cells for surface height z_s. 
void  PythonInterface::PythonIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
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
  a_zSurface.exchange();
}

/// set non-periodic ghost cells for thickness & topography
void  PythonInterface::PythonIBC::setGeometryBCs(LevelSigmaCS& a_coords,
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
  a_coords.getTopography().exchange();
  a_coords.getH().exchange();
}

void  PythonInterface::PythonIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
					 const RealVect& a_dx,
					 const RealVect& a_domainSize,
					 const Real& a_time, 
					 const LevelSigmaCS* a_crseCoords,
					 const int a_refRatio)
{
  CH_TIME("PythonInterface::PythonIBC::initializeIceGeometry");
  if (m_verbose)
    {
      pout() << " PythonIBC::initializeIceGeometry" << endl;
    }

  
  Vector<Real> args(SpaceDim);
  Vector<Real> rval;
  const DisjointBoxLayout& grids = a_coords.grids();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      FArrayBox& topg = a_coords.getTopography()[dit];
      FArrayBox& thck = a_coords.getH()[dit];


      for (BoxIterator bit(grids[dit]);bit.ok();++bit)
	{
	   IntVect iv = bit();
	   int i = 0;
	   for (int dir=0; dir < SpaceDim; dir++)
	     {
	       args[i] = (Real(iv[dir]) + 0.5)*a_dx[dir];
	       i++;
	       PythonEval(m_pFuncThck, rval,  args);
	       thck(iv) = rval[0];
	       PythonEval(m_pFuncTopg, rval,  args);
	       topg(iv) = rval[0]; 	 
	     }
	}
    }
  a_coords.getTopography().exchange();
  a_coords.getH().exchange();
  Real dt = 0.0;
  setGeometryBCs(a_coords, grids.physDomain(), a_dx, a_time, dt);  
}

void  PythonInterface::PythonIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				     const RealVect& a_dx,
				     const RealVect& a_domainSize,
				     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio)
{
  CH_TIME("PythonInterface::PythonIBC::regridIceGeometry");
  Vector<Real> args(SpaceDim);
  Vector<Real> rval;
  const DisjointBoxLayout& grids = a_coords.grids();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      FArrayBox& topg = a_coords.getTopography()[dit];
      for (BoxIterator bit(grids[dit]);bit.ok();++bit)
	{
	   IntVect iv = bit();
	   int i = 0;
	   for (int dir=0; dir < SpaceDim; dir++)
	     {
	       args[i] = (Real(iv[dir]) + 0.5)*a_dx[dir];
	       i++;
	       
	       PythonEval(m_pFuncTopg, rval,  args);
	       topg(iv) = rval[0]; 	 
	     }
	}
    }
  Real dt = 0.0;
  //setGeometryBCs(a_coords, grids.physDomain(), a_dx, a_time, dt);
}


IceThicknessIBC*  PythonInterface::PythonIBC::new_thicknessIBC()
{
  PythonIBC* retval = new PythonIBC( m_pModule, m_pFuncThck, m_pFuncTopg);
  return static_cast<IceThicknessIBC*>(retval);
}


PythonInterface::PythonSurfaceFlux::PythonSurfaceFlux(const std::string& a_pyModuleName,
						      const std::string& a_pyFuncName)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFunc, m_pModule,  a_pyFuncName);
}

PythonInterface::PythonSurfaceFlux::~PythonSurfaceFlux()
{
  Py_DECREF(m_pModule);
  Py_DECREF(m_pFunc);
}

SurfaceFlux* PythonInterface::PythonSurfaceFlux::new_surfaceFlux()
{


  PythonSurfaceFlux* ptr = new PythonSurfaceFlux(m_pModule, m_pFunc);
  return static_cast<SurfaceFlux*>( ptr);
}


void PythonInterface::PythonSurfaceFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
							      const AmrIce& a_amrIce, 
							      int a_level, Real a_dt)
{
   CH_TIME("PythonInterface::PythonSurfaceFlux::surfaceThicknessFlux");
 
  Vector<Real> args(SpaceDim + 3);
  Vector<Real> rval;
  const RefCountedPtr<LevelSigmaCS> levelCS = a_amrIce.geometry(a_level);

  const DisjointBoxLayout& grids = a_flux.disjointBoxLayout();
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      a_flux[dit].setVal(0.0);
      const Box& b = a_flux[dit].box();//grids[dit];
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  int i = 0;
	  for (int dir=0; dir < SpaceDim; dir++)
	    {
	      args[i] = (Real(iv[dir]) + 0.5)*levelCS->dx()[dir];
	      i++;
	    }
	  args[i] = a_amrIce.time();i++;
	  args[i] = levelCS->getH()[dit](iv);i++;
	  args[i] = levelCS->getTopography()[dit](iv);i++;
	  PythonEval(m_pFunc, rval,  args);
	  a_flux[dit](iv) =  rval[0];

	}
    }
}

PythonInterface::PythonBasalFriction::PythonBasalFriction
(const std::string& a_pyModuleName,
 const std::string& a_pyFuncName)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFunc, m_pModule,  a_pyFuncName);
}

PythonInterface::PythonBasalFriction::~PythonBasalFriction()
{
  Py_DECREF(m_pModule);
  Py_DECREF(m_pFunc);
}

BasalFriction* PythonInterface::PythonBasalFriction::new_basalFriction() const
{
  PythonBasalFriction* ptr = new PythonBasalFriction(m_pModule, m_pFunc);
  return static_cast<BasalFriction*>(ptr);
}


void PythonInterface::PythonBasalFriction::setBasalFriction
(LevelData<FArrayBox>& a_C,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{
  CH_TIME("PythonInterface::PythonBasalFriction::setBasalFriction");
  Vector<Real> args(SpaceDim + 3);
  Vector<Real> rval;
  
  const DisjointBoxLayout& grids = a_C.disjointBoxLayout();
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      a_C[dit].setVal(0.0);
      const Box& b = a_C[dit].box();//grids[dit];
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  int i = 0;
	  for (int dir=0; dir < SpaceDim; dir++)
	    {
	      args[i] = (Real(iv[dir]) + 0.5)*a_coordSys.dx()[dir];
	      i++;
	    }
	  args[i] = a_time;i++;
	  args[i] = a_coordSys.getH()[dit](iv);i++;
	  args[i] = a_coordSys.getTopography()[dit](iv);i++;
	  PythonEval(m_pFunc, rval,  args);
	  a_C[dit](iv) =  rval[0];
	  
	}
    }
}



#include "NamespaceFooter.H"

#endif
