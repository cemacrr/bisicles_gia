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
#include "PiecewiseLinearFlux.H"
#include "IceConstants.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "BisiclesF_F.H"
#include "ParmParse.H"
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
                               LevelSigmaCS& a_coordSys,
                               Real a_time,
                               Real a_dt)
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
                                   LevelSigmaCS& a_coordSys,
                                   Real a_time,
                                   Real a_dt)
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
                                           LevelSigmaCS& a_coordSys,
                                           Real a_time,
                                           Real a_dt)
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
  RealVect dx = a_coordSys.dx();

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
MaskedFlux::MaskedFlux(SurfaceFlux* a_groundedIceFlux, SurfaceFlux* a_floatingIceFlux)
  :m_groundedIceFlux(a_groundedIceFlux),m_floatingIceFlux(a_floatingIceFlux)
{

}
/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* MaskedFlux::new_surfaceFlux()
{
  SurfaceFlux* f = m_floatingIceFlux->new_surfaceFlux();
  SurfaceFlux* g = m_groundedIceFlux->new_surfaceFlux();
  return static_cast<SurfaceFlux*>(new MaskedFlux(g,f));
}

void MaskedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				      LevelSigmaCS& a_coordSys,
				      Real a_time,
				      Real a_dt)
{

  //somewhat ineffcient, because we compute all fluxes everywhere. 
  //At some point, come back and only compute (say) grounded ice flux
  //in boxes where all the ice is grounded.

  //first, grounded ice values
  m_groundedIceFlux->surfaceThicknessFlux(a_flux, a_coordSys, a_time, a_dt);

  //floating ice values
  LevelData<FArrayBox> floatingFlux;
  floatingFlux.define(a_flux);
  m_floatingIceFlux->surfaceThicknessFlux(floatingFlux, a_coordSys, a_time, a_dt);

  for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask =  a_coordSys.getFloatingMask()[dit];
      
      Box box = mask.box();
      box &= a_flux[dit].box();
      int m = FLOATINGMASKVAL;

      FORT_MASKEDREPLACE(CHF_FRA1(a_flux[dit],0),
			 CHF_CONST_FRA1(floatingFlux[dit],0),
			 CHF_CONST_FIA1(mask,0),
			 CHF_CONST_INT(m),
			 CHF_BOX(box));
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

      ptr = static_cast<SurfaceFlux*>
	(new MaskedFlux(groundedPtr->new_surfaceFlux(),floatingPtr->new_surfaceFlux()));
      
      delete groundedPtr;
      delete floatingPtr;
    }

  return ptr;
  
}



#include "NamespaceFooter.H"
