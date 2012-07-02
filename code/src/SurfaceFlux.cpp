#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SurfaceFlux.H"
#include "LevelDataSurfaceFlux.H"
#include "GroundingLineLocalizedFlux.H"
#include "IceConstants.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "BisiclesF_F.H"
#include "ParmParse.H"
#include "AmrIce.H"
#include "NamespaceHeader.H"

  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */

SurfaceFlux* 
zeroFlux::new_surfaceFlux()
{
  zeroFlux* newPtr = new zeroFlux;
  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
zeroFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
			       const AmrIce& a_amrIce, 
			       int a_level, Real a_dt)
{
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }
}


constantFlux::constantFlux() : m_isValSet(false)
{
}

SurfaceFlux* 
constantFlux::new_surfaceFlux()
{
  constantFlux* newPtr = new constantFlux;
  newPtr->m_fluxVal = m_fluxVal;
  newPtr->m_isValSet = m_isValSet;
  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
constantFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				   const AmrIce& a_amrIce, 
				   int a_level, Real a_dt)
{
  CH_assert(m_isValSet);
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(m_fluxVal);
    }
}


///
void
constantFlux::setFluxVal(const Real& a_fluxVal) 
{
  m_fluxVal = a_fluxVal; 
  // input value is in meters/year divide by secondsperyear 
  // to get flux in meters/second
  // slc: switch to flux in m/a
  //m_fluxVal/= secondsperyear;
  
  m_isValSet = true;
}


// --------------------------------------------------------------
// fortran interface surface flux
// --------------------------------------------------------------

/// class which takes an input fortran array 
/** averages or interpolates as necessary to fill the flux
 */

/// constructor
fortranInterfaceFlux::fortranInterfaceFlux()
  : m_isValSet(false)
{
}

/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* 
fortranInterfaceFlux::new_surfaceFlux()
{
  fortranInterfaceFlux* newPtr = new fortranInterfaceFlux;
  newPtr->m_fluxVal.define(m_fluxVal.interval(), m_fluxVal);
  newPtr->m_isValSet = m_isValSet;
  return static_cast<SurfaceFlux*>(newPtr);  
}

/// define source term for thickness evolution and place it in flux
/** dt is included in case one needs integrals or averages over a
    timestep. flux should be defined in meters/second in the current 
    implementation. 
*/
void 
fortranInterfaceFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					   const AmrIce& a_amrIce, 
					   int a_level, Real a_dt)
{
  CH_assert(m_isValSet);

  // this looks a lot like the code in FortranInterfaceIBC
  Vector<Box> importBoxes(1, m_fluxVal.box());
  Vector<int> procAssign(1,0);
  DisjointBoxLayout importDBL(importBoxes, procAssign);
  LevelData<FArrayBox> importFlux(importDBL, 1, IntVect::Zero);

  DataIterator importDit = importDBL.dataIterator();
  for (importDit.begin(); importDit.ok(); ++importDit)
    {
      importFlux[importDit].copy(m_fluxVal);
    }

  DisjointBoxLayout levelGrids = a_flux.getBoxes();
  RealVect dx = a_amrIce.dx(a_level);

  // refinement ratio for flux
  Real refRatio = m_fluxDx[0]/dx[0];
 
  Real tolerance = 1.0e-6;

  if (refRatio > 1 + tolerance)
    {
      // importFlux coarser than what we want, have to interpolate
      //int nRef = (int)(refRatio + tolerance);

    }
  else if (refRatio < 1 + tolerance)
    {
      // importFlux finer than what we want, have to average
      //int nRef = (int)(1.0/refRatio + tolerance);
    }
  else
    {
      // same size, just copy  
      importFlux.copyTo(a_flux);
    }
  


}

/// set fortran array-valued surface flux
void
fortranInterfaceFlux::setFluxVal(Real* a_data_ptr,
                                 const int* a_dimInfo,
                                 const Real* a_dew, const Real* a_dns)
{
  // dimInfo is (SPACEDIM, nz, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering
  IntVect hiVect(D_DECL(a_dimInfo[2]-1,a_dimInfo[3]-1, a_dimInfo[1]-1));
  Box fabBox(IntVect::Zero, hiVect);

  //cout << "hiVect" << hiVect << endl;
  
  m_fluxVal.define(fabBox, 1, a_data_ptr);
  m_fluxDx = RealVect(D_DECL(*a_dew, *a_dns, 1));  
}

 /// constructor
MaskedFlux::MaskedFlux(SurfaceFlux* a_groundedIceFlux, SurfaceFlux* a_floatingIceFlux,
		       SurfaceFlux* a_openSeaFlux, SurfaceFlux* a_openLandFlux)
  :m_groundedIceFlux(a_groundedIceFlux),m_floatingIceFlux(a_floatingIceFlux),
   m_openSeaFlux(a_openSeaFlux),m_openLandFlux(a_openLandFlux)
{
  CH_assert(a_groundedIceFlux);
  CH_assert(a_floatingIceFlux);
  CH_assert(a_openSeaFlux);
  CH_assert(a_openLandFlux);
}
/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* MaskedFlux::new_surfaceFlux()
{
  SurfaceFlux* f = m_floatingIceFlux->new_surfaceFlux();
  SurfaceFlux* g = m_groundedIceFlux->new_surfaceFlux();
  SurfaceFlux* s = m_openSeaFlux->new_surfaceFlux();
  SurfaceFlux* l = m_openLandFlux->new_surfaceFlux();
  return static_cast<SurfaceFlux*>(new MaskedFlux(g,f,s,l));
}

void MaskedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				      const AmrIce& a_amrIce, 
				      int a_level, Real a_dt)
{

  //somewhat ineffcient, because we compute all fluxes everywhere. 
  //At some point, come back and only compute (say) grounded ice flux
  //in boxes where all the ice is grounded.

  //first, grounded ice values
  if (m_groundedIceFlux)
    m_groundedIceFlux->surfaceThicknessFlux(a_flux,a_amrIce,a_level,a_dt);

  //floating ice values
  LevelData<FArrayBox> floatingFlux;
  floatingFlux.define(a_flux);

  if (m_floatingIceFlux)
    m_floatingIceFlux->surfaceThicknessFlux(floatingFlux, a_amrIce,a_level,a_dt);

  for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask =  a_amrIce.geometry(a_level)->getFloatingMask()[dit];
      
      Box box = mask.box();
      box &= a_flux[dit].box();

      int m = FLOATINGMASKVAL;
      FORT_MASKEDREPLACE(CHF_FRA1(a_flux[dit],0),
			 CHF_CONST_FRA1(floatingFlux[dit],0),
			 CHF_CONST_FIA1(mask,0),
			 CHF_CONST_INT(m),
			 CHF_BOX(box));

      
      m = OPENSEAMASKVAL;
      FORT_MASKEDREPLACE(CHF_FRA1(a_flux[dit],0),
			 CHF_CONST_FRA1(floatingFlux[dit],0),
			 CHF_CONST_FIA1(mask,0),
			 CHF_CONST_INT(m),
			 CHF_BOX(box));

    }

}


/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* CompositeFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new CompositeFlux(m_fluxes));
}

CompositeFlux::CompositeFlux(const Vector<SurfaceFlux*>& a_fluxes)
{
  m_fluxes.resize(a_fluxes.size());
  for (int i =0; i < a_fluxes.size(); i++)
    {
      CH_assert(a_fluxes[i] != NULL);
      m_fluxes[i] =  a_fluxes[i]->new_surfaceFlux();
    }
}

CompositeFlux::~CompositeFlux()
{
  for (int i =0; i < m_fluxes.size(); i++)
    {
      if (m_fluxes[i] != NULL)
	{
	  delete m_fluxes[i];
	  m_fluxes[i] = NULL;
	}
    }
}

void CompositeFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIce& a_amrIce, 
					  int a_level, Real a_dt)
{
  m_fluxes[0]->surfaceThicknessFlux(a_flux, a_amrIce, a_level,a_dt );

  // this is hardly effcient... but it is convenient
  LevelData<FArrayBox> tmpFlux(a_flux.disjointBoxLayout(),1,a_flux.ghostVect());
  for (int i = 1; i <  m_fluxes.size(); i++)
    {
      m_fluxes[i]->surfaceThicknessFlux(tmpFlux, a_amrIce, a_level,a_dt );
      for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  a_flux[dit] += tmpFlux[dit];
	}
    }

}



SurfaceFlux* BoxBoundedFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>( new BoxBoundedFlux(m_lo, m_hi,m_startTime,m_endTime,m_fluxPtr));
}


void BoxBoundedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIce& a_amrIce, 
					  int a_level, Real a_dt)
{

  Real time = a_amrIce.time();
  

  for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }

  if (time >= m_startTime && time < m_endTime)
    {
      // this is hardly efficient... but it is convenient
      LevelData<FArrayBox> tmpFlux(a_flux.disjointBoxLayout(),1,a_flux.ghostVect());
      m_fluxPtr->surfaceThicknessFlux(tmpFlux, a_amrIce, a_level,a_dt);
      const RealVect& dx = a_amrIce.dx(a_level);
      
      IntVect ilo,ihi;
      for (int dir =0; dir < SpaceDim; dir++)
	{
	  ilo[dir] = int(m_lo[dir]/dx[dir] - 0.5);
	  ihi[dir] = int(m_hi[dir]/dx[dir] - 0.5);
	}
      
    
      for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  const Box& b = a_flux[dit].box();
	  if (b.intersects(Box(ilo,ihi)))
	    { 
	      Box sub(max(b.smallEnd(),ilo), min(b.bigEnd(),ihi));
	      a_flux[dit].plus(tmpFlux[dit],sub,0,0,1);
	    }
	}
    }

}

PiecewiseLinearFlux::PiecewiseLinearFlux(const Vector<Real>& a_abscissae, 
					 const Vector<Real>& a_ordinates)
  :m_abscissae(a_abscissae),m_ordinates(a_ordinates)
{
  CH_assert(m_abscissae.size() == m_ordinates.size());
}


SurfaceFlux* PiecewiseLinearFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new PiecewiseLinearFlux(m_abscissae,m_ordinates));
}

void PiecewiseLinearFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					       const AmrIce& a_amrIce, 
					       int a_level, Real a_dt)
{
  Vector<Real> dx(m_abscissae.size());
  Vector<Real> db(m_abscissae.size());

  const LevelData<FArrayBox>& levelH = a_amrIce.geometry(a_level)->getH();

  for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
    {
      FORT_PWLFILL(CHF_FRA1(a_flux[dit],0),
		   CHF_CONST_FRA1(levelH[dit],0),
		   CHF_CONST_VR(m_abscissae),
		   CHF_CONST_VR(m_ordinates),
		   CHF_VR(dx),CHF_VR(db),
		   CHF_BOX(a_flux[dit].box()));
    }

}


SurfaceFlux* SurfaceFlux::parseSurfaceFlux(const char* a_prefix)
{
  
  SurfaceFlux* ptr = NULL;
  std::string type = "";
  
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "zeroFlux")
    {
      ptr = new zeroFlux;
    }
  else if (type == "constantFlux")
    {
      constantFlux* constFluxPtr = new constantFlux;
      Real fluxVal;
      pp.get("flux_value", fluxVal);
      constFluxPtr->setFluxVal(fluxVal);
      ptr = static_cast<SurfaceFlux*>(constFluxPtr);
    }
  else if (type == "LevelData")
    {
      std::string fileFormat;
      pp.get("fileFormat",fileFormat);
      int n;
      pp.get("n",n);
      Real startTime = 0.0, timeStep = 1.0;
      pp.query("startTime", startTime);
      pp.query("timeStep", timeStep);
      std::string name = "flux";
      pp.query("name", name);
      	
      RefCountedPtr<std::map<Real,std::string> > tf
	(new std::map<Real,std::string>);
      
      for (int i =0; i < n; i++)
	{
	  char* file = new char[fileFormat.length()+32];
	  sprintf(file, fileFormat.c_str(),i);
	  tf->insert(make_pair(startTime + Real(i)*timeStep, file));
	}
      
      LevelDataSurfaceFlux* ldptr = new LevelDataSurfaceFlux(tf,name);
      ptr = static_cast<SurfaceFlux*>(ldptr);
    }
  else if (type == "piecewiseLinearFlux")
    {
      int n = 1;  
      pp.query("n",n);
      Vector<Real> vabs(n,0.0);
      Vector<Real> vord(n,0.0);
      pp.getarr("abscissae",vabs,0,n);
      pp.getarr("ordinates",vord,0,n);
      PiecewiseLinearFlux* pptr = new PiecewiseLinearFlux(vabs,vord);
      ptr = static_cast<SurfaceFlux*>(pptr);
    }
  else if (type == "maskedFlux")
    {
      std::string groundedPrefix(a_prefix);
      groundedPrefix += ".grounded";
      SurfaceFlux* groundedPtr = parseSurfaceFlux(groundedPrefix.c_str());
      if (groundedPtr == NULL)
	{
	  groundedPtr = new zeroFlux;
	}
 
      std::string floatingPrefix(a_prefix);
      floatingPrefix += ".floating";
      SurfaceFlux* floatingPtr = parseSurfaceFlux(floatingPrefix.c_str());
      if (floatingPtr == NULL)
	{
	  floatingPtr = new zeroFlux;
	}

      std::string openLandPrefix(a_prefix);
      openLandPrefix += ".openLand";
      SurfaceFlux* openLandPtr = parseSurfaceFlux(openLandPrefix.c_str());
      if (openLandPtr == NULL)
	{
	  openLandPtr = new zeroFlux;
	}

      
      std::string openSeaPrefix(a_prefix);
      openSeaPrefix += ".openSea";
      SurfaceFlux* openSeaPtr = parseSurfaceFlux(openSeaPrefix.c_str());
      if (openSeaPtr == NULL)
	{
	  openSeaPtr = new zeroFlux;
	}

      ptr = static_cast<SurfaceFlux*>
	(new MaskedFlux(groundedPtr->new_surfaceFlux(),
			floatingPtr->new_surfaceFlux(),
			openSeaPtr->new_surfaceFlux(),
			openLandPtr->new_surfaceFlux()));
      
      delete groundedPtr;
      delete floatingPtr;
      delete openSeaPtr;
      delete openLandPtr;
    }
  else if (type == "boxBoundedFlux")
    {
      Vector<Real> tmp(SpaceDim); 
      pp.getarr("lo",tmp,0,SpaceDim);
      RealVect lo (D_DECL(tmp[0], tmp[1],tmp[2]));
      pp.getarr("hi",tmp,0,SpaceDim);
      RealVect hi (D_DECL(tmp[0], tmp[1],tmp[2]));

      Vector<Real> time(2);
      time[0] = -1.2345678e+300;
      time[1] = 1.2345678e+300;
      pp.queryarr("time",time,0,2);
           
      std::string prefix(a_prefix);
      prefix += ".flux";
      SurfaceFlux* fluxPtr = parseSurfaceFlux(prefix.c_str());

      BoxBoundedFlux bbf(lo,hi,time[0],time[1],fluxPtr);
      ptr = static_cast<SurfaceFlux*>(bbf.new_surfaceFlux());

    }
  else if (type == "compositeFlux")
   {
     
     
     int nElements;
     pp.get("nElements",nElements);
     
     std::string elementPrefix(a_prefix);
     elementPrefix += ".element";

     Vector<SurfaceFlux*> elements(nElements);
     for (int i = 0; i < nElements; i++)
       {
	 std::string prefix(elementPrefix);
	 char s[32];
	 sprintf(s,"%i",i);
	 prefix += s;
	 ParmParse pe(prefix.c_str());
	 elements[i] = parseSurfaceFlux(prefix.c_str());
	 CH_assert(elements[i] != NULL);
	 
       }
     CompositeFlux compositeFlux(elements);
     ptr = static_cast<SurfaceFlux*>(compositeFlux.new_surfaceFlux());
   
   }
  else if (type == "groundingLineLocalizedFlux")
    {
      std::string glPrefix(a_prefix);
      glPrefix += ".groundingLine";
      SurfaceFlux* glPtr = parseSurfaceFlux(glPrefix.c_str());
      if (glPtr == NULL)
	{
	  glPtr = new zeroFlux;
	}
       
      std::string ambientPrefix(a_prefix);
      ambientPrefix += ".ambient";
      SurfaceFlux* ambientPtr = parseSurfaceFlux(ambientPrefix.c_str());
      if (ambientPtr == NULL)
	{
	  ambientPtr = new zeroFlux;
	}
      
      ptr = static_cast<SurfaceFlux*>
	(new GroundingLineLocalizedFlux(glPtr->new_surfaceFlux(),
					ambientPtr->new_surfaceFlux()));
	 
      delete glPtr;
      delete ambientPtr;

    }
#ifdef HAVE_PYTHON
  else if (type == "pythonFlux") {
    
    std::string module;
    pp.get("module",module);
    std::string function;
    pp.get("function",function);
    PythonSurfaceFlux pythonFlux(module, function);
    ptr = static_cast<SurfaceFlux*>(pythonFlux.new_surfaceFlux());

  }
#endif
  return ptr;
  
}


#ifdef HAVE_PYTHON
#include "signal.h"

PythonSurfaceFlux::PythonSurfaceFlux(const std::string& a_pyModule,
				     const std::string& a_pyFunc)
{
  
  //initialize Python without overriding SIGINT  
  PyOS_sighandler_t sigint = PyOS_getsig(SIGINT);
  Py_Initialize();
  sigint =  PyOS_setsig(SIGINT, sigint);
  
  PyObject* pName = PyString_FromString(a_pyModule.c_str());
  m_pModule =  PyImport_Import(pName);

  if (m_pModule == NULL)
    {
      pout() << "failed to import Python module " << a_pyModule << endl;
      CH_assert(m_pModule != NULL);
      MayDay::Error("failed to import Python module ");
    }
  char* t = const_cast<char*>(a_pyFunc.c_str());
  m_pFunc = PyObject_GetAttrString(m_pModule,t);
  
  if (m_pFunc == NULL)
    {
      pout() << "failed to initialize Python function " << a_pyFunc << endl;
      CH_assert(m_pFunc != NULL);
      MayDay::Error("failed to initialize Python function");
    }

  if (!PyCallable_Check(m_pFunc))
    {
      pout() << "Python function " << a_pyFunc << " not callable " <<  endl;
      CH_assert(PyCallable_Check(m_pFunc));
      MayDay::Error("Python function not callable");
    }


}

PythonSurfaceFlux::~PythonSurfaceFlux()
{
  Py_DECREF(m_pModule);
  Py_DECREF(m_pFunc);
}

SurfaceFlux* PythonSurfaceFlux::new_surfaceFlux()
{


  PythonSurfaceFlux* ptr = new PythonSurfaceFlux(m_pModule, m_pFunc);
  return static_cast<SurfaceFlux*>( ptr);
}


void PythonSurfaceFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					     const AmrIce& a_amrIce, 
					     int a_level, Real a_dt)
{
  
  PyObject *pArgs, *pValue;
  Vector<Real> args(SpaceDim + 3);

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
	  
	  //construct a python tuple of args
	  pArgs = PyTuple_New(args.size());
	  for (int j = 0; j < args.size(); ++j)
	    {
	      pValue = PyFloat_FromDouble(args[j]);
	      if (!pValue)
		{
		  pout() << "PythonSurfaceFlux::surfaceThicknessFlux failed to construct Python args" << endl;
		  CH_assert(pValue != NULL);
		  MayDay::Error("PythonSurfaceFlux::surfaceThicknessFlux failed to construct Python args");
		}
	      PyTuple_SetItem(pArgs, j, pValue);
	    }
	  pValue = PyObject_CallObject(m_pFunc, pArgs);
	  Py_DECREF(pArgs);
	  if (!pValue)
	    {
	      pout() << "PythonSurfaceFlux::surfaceThicknessFlux python call failed" << endl;
	      PyErr_Print();
	      CH_assert(pValue != NULL);
	      MayDay::Error("PythonSurfaceFlux::surfaceThicknessFlux python call failed");
	    }
	  
	  if (PyFloat_CheckExact(pValue))
	    {
	      a_flux[dit](iv) =  PyFloat_AS_DOUBLE(pValue);
	    }
	  else
	    {
	      Py_DECREF(pValue);
	      CH_assert(PyFloat_CheckExact(pValue));
	      MayDay::Error("PythonSurfaceFlux::surfaceThicknessFlux python return value is not a float");
	    }


	}
    }
}
#endif
#include "NamespaceFooter.H"
