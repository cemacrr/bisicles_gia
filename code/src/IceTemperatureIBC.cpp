 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "IceTemperatureIBC.H"
#include "IceConstants.H"
#include "BoxIterator.H"
#include "ExtrapGhostCells.H"
#include "NamespaceHeader.H"

ConstantIceTemperatureIBC* 
ConstantIceTemperatureIBC::new_temperatureIBC()
{
  return new ConstantIceTemperatureIBC(m_T);
}

/// set a basal heat flux to zero. units are Joules / Year
void ConstantIceTemperatureIBC::basalHeatFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIce& a_amrIce, 
 int a_level, Real a_dt)
{
  for (DataIterator dit = a_flux.dataIterator();dit.ok();++dit)
    {
      a_flux[dit].setVal(0.0);
    }
}

#if BISICLES_Z == BISICLES_LAYERED
void 
ConstantIceTemperatureIBC::initializeIceTemperature
(LevelData<FArrayBox>& a_T,
 LevelData<FArrayBox>& a_surfaceT, 
 LevelData<FArrayBox>& a_basalT,
 const LevelSigmaCS& a_coordSys)
{
  
  for (DataIterator dit(a_T.disjointBoxLayout()); dit.ok() ; ++dit)
    {
      for (int l = 0; l < a_T[dit].nComp(); ++l)
	{
	  a_T[dit].setVal(m_T, l);
	}	    
      a_surfaceT[dit].setVal(m_T);
      a_basalT[dit].setVal(m_T);
    }
}

void 
ConstantIceTemperatureIBC::setIceTemperatureBC
(LevelData<FArrayBox>& a_T,
 LevelData<FArrayBox>& a_surfaceT, 
 LevelData<FArrayBox>& a_basalT,
 const LevelSigmaCS& a_coordSys)
{
  const ProblemDomain& domain = a_coordSys.grids().physDomain();
  ExtrapGhostCells( a_T, domain);
  ExtrapGhostCells( a_basalT, domain);
}


#elif BISICLES_Z == BISICLES_FULLZ
void 
ConstantIceTemperatureIBC::initializeIceTemperature
(LevelData<FArrayBox>& a_T,
 const LevelSigmaCS& a_coordSys)
{
   
  for (DataIterator dit(a_T.disjointBoxLayout()); dit.ok() ; ++dit)
    {  
       a_T[dit].setVal(m_T, l);
    }
}

void 
ConstantIceTemperatureIBC::setIceTemperatureBC
(LevelData<FArrayBox>& a_T,
 const LevelSigmaCS& a_coordSys)
{
  const ProblemDomain& domain = a_coordSys.grids().physDomain();
  ExtrapGhostCells( a_T, domain);
}

#endif

void 
IceTemperatureIBC::initialize(LevelData<FArrayBox>& a_U)
{
  //we shouldn't get here
  MayDay::Error("IceTemperatureIBC::initialize not implemented");
}

void  IceTemperatureIBC::primBC
(FArrayBox&            a_WGdnv,
 const FArrayBox&      a_Wextrap,
 const FArrayBox&      a_W,
 const int&            a_dir,
 const Side::LoHiSide& a_side,
 const Real&           a_time)
{

  // do nothing in periodic case
  if (!m_domain.isPeriodic(a_dir))
    {
      
#if CH_SPACEDIM == 2
      // 2D case : We are at either an ice divide (where u = 0)
      // or an outflow (where extrapolation is fine)
      // \todo : support an inflow with a known temperature?
   int lohisign;
   
   //find the strip of cells just inside the domain. DFM notes that 
   //a_WGdnv might have more than one layer of ghost cells, so this
   //is slightly more complicated than just testing tmp is at the domain
   //edge
   Box tmp = a_WGdnv.box();
   lohisign = sign(a_side);
   tmp.shiftHalf(a_dir,lohisign);
   const Box& dbox = m_domain.domainBox();
   Box ghostBox = (a_side == Side::Lo)?adjCellLo(dbox,a_dir,1):adjCellHi(dbox,a_dir,1);
   ghostBox &= tmp;
 
   if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
     {
       tmp &= m_domain;
       Box boundaryBox = (a_side == Side::Lo)?bdryLo(tmp,a_dir):bdryHi(tmp,a_dir);
       BoxIterator bit(boundaryBox);
       for (bit.begin(); bit.ok(); ++bit){
	 const IntVect& i = bit();
	 a_WGdnv(i,0) = a_Wextrap(i,0);
	 
       }
     }

#elif CH_SPACEDIM == 3
  MayDay::Error("IceThicknessIBC::primBC DIM = 3 not yet impplemented");
#else
  MayDay::Error("IceThicknessIBC::primBC only supports DIM = 2 and DIM = 3");
#endif

    }
}



#include "NamespaceFooter.H"
