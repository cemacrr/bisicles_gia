#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CalvingModel.H"
#include "BennCalvingModel.H"
#include "FlotationCalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

void 
DeglaciationCalvingModelA::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	}
    }
}

void DomainEdgeCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  //const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
  const IntVect ghost = a_thickness.ghostVect();
  //const LevelData<FArrayBox>& vt  = *a_amrIce.viscousTensor(a_level);
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      //const Box& gridBox = grids[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  if (!domain.isPeriodic(dir))
	    {

	      if (m_frontLo[dir] > 0)
		{
		  Box loBox = adjCellLo(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  loBox.grow(transverseVect);
		  loBox &= a_thickness[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv + BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
			a_thickness[dit](iv) = 0.0;
		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  hiBox.grow(transverseVect);
		  hiBox &= a_thickness[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv - BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
			a_thickness[dit](iv) = 0.0;
		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  thck(iv) = std::max(thck(iv),0.0);
	}

    } // end loop over boxes

}

void ProximityCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{

  Real time = a_amrIce.time();
  bool calvingActive = (time >= m_startTime && time < m_endTime);
  calvingActive = false;
  pout() << " time = " << time 
	 << " m_startTime = " <<  m_startTime
	 << " m_endTime = " <<  m_endTime
	 << "calvingActive = " << calvingActive
	 << std::endl;
  if (true || calvingActive)
    {
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& proximity = *a_amrIce.groundingLineProximity(a_level);
      const LevelData<FArrayBox>& velocity = *a_amrIce.velocity(a_level);
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& thck = a_thickness[dit];
	  const FArrayBox& prox = proximity[dit];
	  const FArrayBox& vel = velocity[dit];
	  Box b = thck.box();b &= prox.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real vmod = std::sqrt(vel(iv,0)*vel(iv,0) + vel(iv,1)*vel(iv,1));
	      if (prox(iv) < m_proximity && calvingActive && vmod > m_velocity)
		{
		  //thck(iv) *= 0.5; thck(iv) = max(thck(iv),10.0);
		  thck(iv) = 0.0;
		}
	      if (mask(iv) == OPENSEAMASKVAL)
		{
		   thck(iv) = 0.0;
		}
	      if (mask(iv) == FLOATINGMASKVAL)
		{
		  thck(iv) = max(thck(iv),1.0);
		}
	    }
	}
    }
}



void ProximityCalvingModel::modifySurfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 LevelData<FArrayBox>& a_calvingMelt,
 const AmrIce& a_amrIce,
 int a_level,
 Real a_dt)
{

  Real time = a_amrIce.time();
  bool calvingActive = (time >= m_startTime && time < m_endTime);
  pout() << " time = " << time 
	 << " m_startTime = " <<  m_startTime
	 << " m_endTime = " <<  m_endTime
	 << " calvingActive = " << calvingActive
	 << std::endl;
  if (calvingActive)
    {
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& proximity = *a_amrIce.groundingLineProximity(a_level);
      const LevelData<FArrayBox>& velocity = *a_amrIce.velocity(a_level);
      
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  //const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& flux = a_flux[dit];
	  FArrayBox& melt = a_calvingMelt[dit];
	  const FArrayBox& prox = proximity[dit];
	  const FArrayBox& vel = velocity[dit];
	  const FArrayBox& thck = levelCoords.getH()[dit];
	  Box b = flux.box();b &= prox.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real extraMelt = 0.0;
	      Real vmod = std::sqrt(vel(iv,0)*vel(iv,0) + vel(iv,1)*vel(iv,1));
	      if (prox(iv) < m_proximity && calvingActive && vmod > m_velocity)
	  	{
	  	  //flux(iv) -= 0.95 * thck(iv)/a_dt;
		  extraMelt = 0.95 * thck(iv)/a_dt;
		  flux(iv) -= extraMelt;
	  	}
	      melt(iv) = extraMelt;
	    }
	}
    }

  
}



CalvingModel* CalvingModel::parseCalvingModel(const char* a_prefix)
{

  CalvingModel* ptr = NULL;
  std::string type = "";
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "NoCalvingModel")
    {
      ptr = new NoCalvingModel;
    }
  else if (type == "DomainEdgeCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new DomainEdgeCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "FixedFrontCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      ptr = new DeglaciationCalvingModelA
	(0.0,  1.0e+10, minThickness, -1.2345678e+300, 1.2345678e+300);
    }
  else if (type == "DeglaciationCalvingModelA")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.get("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelA
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "DeglaciationCalvingModelB")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelB
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "ProximityCalvingModel")
    {
      Real proximity = 0.0;
      pp.get("proximity", proximity );
      Real velocity = 0.0;
      pp.query("velocity", velocity );
      Real startTime = -1.2345678e+300;
      pp.get("startTime",  startTime);
      Real endTime = 1.2345678e+300;
      pp.get("endTime",  endTime);
      ptr = new ProximityCalvingModel(proximity,velocity, startTime, endTime);
    }
  else if (type == "FlotationCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new FlotationCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "BennCalvingModel")
    {
      Real waterDepth = 0.0;
      pp.get("waterDepth",waterDepth);
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      bool inclBasalCrev = false;
      pp.query("inclBasalCrev",inclBasalCrev);
      bool oneDimCrev = false;
      pp.query("oneDimCrev",oneDimCrev);
      bool RedFreqCalv = false;
      pp.query("RedFreqCalv",RedFreqCalv);
      bool setZeroThck = false;
      pp.query("setZeroThck",setZeroThck);
      bool oldMeltRate = false;
      pp.query("oldMeltRate",oldMeltRate);
      Real NoiseScale = 0.0;
      pp.query("NoiseScale",NoiseScale);
      Real calvingFreq = 1.0;
      pp.query("calvingFreq",calvingFreq);
      Real decay = 3.0e-1;
      pp.query("decay",decay);
      Real timescale = 0.25;
      pp.query("timescale",timescale);
      ptr = new BennCalvingModel(waterDepth, frontLo, frontHi, preserveSea,
				 preserveLand, inclBasalCrev, oneDimCrev, 
				 RedFreqCalv, setZeroThck,
				 oldMeltRate, NoiseScale, calvingFreq,
				 decay, timescale);
    }
  else if (type == "ThicknessCalvingModel")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new ThicknessCalvingModel
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "MaxiumumExtentCalvingModel")  
    {

      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);

      Vector<Real> vect(SpaceDim,0.0);

      pp.getarr("lowLoc",vect,0,SpaceDim);
      RealVect lowLoc(D_DECL(vect[0], vect[1],vect[2]));      

      pp.getarr("highLoc",vect,0,SpaceDim);
      RealVect highLoc(D_DECL(vect[0], vect[1],vect[2]));      

      MaximumExtentCalvingModel* Ptr = new MaximumExtentCalvingModel(highLoc,
                                                                     lowLoc,
                                                                     startTime,
                                                                     endTime);
      ptr = static_cast<CalvingModel*>(Ptr);

    }
  else if (type == "CompositeCalvingModel")
    {
      int nElements;
      pp.get("nElements",nElements);
     
      std::string elementPrefix(a_prefix);
      elementPrefix += ".element";

      Vector<CalvingModel*> elements(nElements);
      for (int i = 0; i < nElements; i++)
        {
          std::string prefix(elementPrefix);
          char s[32];
          sprintf(s,"%i",i);
          prefix += s;
          ParmParse pe(prefix.c_str());
          elements[i] = parseCalvingModel(prefix.c_str());
          CH_assert(elements[i] != NULL);
        }
      CompositeCalvingModel* compositePtr = new CompositeCalvingModel(elements);
      ptr = static_cast<CalvingModel*>(compositePtr);
    }   

//   else
//     {
//       pout() << "Unknown calving model" << type << std::endl;
//       MayDay::Error("Unknown calving model");
//     }
  
  return ptr;
}


void 
DeglaciationCalvingModelB::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          else if ((mask(iv) == FLOATINGMASKVAL) && (thck(iv) < m_calvingThickness))
            {
	      thck(iv) = m_minThickness;              
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	}
    }
}


void 
ThicknessCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{
  // actually need a non-const AmrIce here because we also want to 
  // change the iceMask. 
  AmrIce* nonConstAmrIce = const_cast<AmrIce*>(&a_amrIce);
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  LevelData<FArrayBox>& levelIceMask = *nonConstAmrIce->iceMask(a_level);

  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& iceMask = levelIceMask[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox effectiveThickness(thck.box(), 1);
      effectiveThickness.copy(thck);
      Box b = thck.box();
      b &= iceMask.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();          
          // if iceMask > 0, then rescale effectiveThickness
          // by dividing by iceMask value, which gives "actual" thickness
          // in the partial cell. Probably eventually want to move this to 
          // fortran
          if (iceMask(iv,0) > 0.0)
            {
              effectiveThickness(iv,0) /= iceMask(iv,0);
            }
            
          if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          // allow ice to spread into open sea regions too, if appropriate
          else if (((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
                   && (effectiveThickness(iv) < m_calvingThickness))
            {
              // note that we're setting thck here, not effectiveThickness, 
              // which is a temporary
              // also set the iceMask to zero in these cells
	      thck(iv) = m_minThickness; 
              iceMask(iv,0) = 0.0;
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	}
    }
}



  
//alter the thickness field at the end of a time step
void
MaximumExtentCalvingModel::endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
                                                  const AmrIce& a_amrIce,
                                                  int a_level)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const Real dx = a_amrIce.amrDx()[a_level];
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();

          // compute location of cell center
          RealVect loc(iv);          
          loc += 0.5*RealVect::Unit;
          loc *= dx;
          
          // check high and low extents
          if ((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
            {
              if (loc[0] <= m_lowLoc[0])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] <= m_lowLoc[1])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[0] > m_highLoc[0]) 
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] > m_highLoc[1])
                {
                  thck(iv) = 0.0;
                }
              
            } // end if floating or opensea
        } // end loop over cells
  
    }

}


//alter the thickness field prior to the first time step. 
void
CompositeCalvingModel::initialModifyState(LevelData<FArrayBox>& a_thickness, 
                                          const AmrIce& a_amrIce,
                                          int a_level)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->initialModifyState(a_thickness, a_amrIce, a_level);
    }
}


//alter the thickness field at the end of a time step
void 
CompositeCalvingModel::endTimeStepModifyState(LevelData<FArrayBox>& a_thickness, 
                                              const AmrIce& a_amrIce,
                                              int a_level)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->endTimeStepModifyState(a_thickness, a_amrIce, 
                                              a_level);
    }
}

  
//alter the thickness field after a regrid
void 
CompositeCalvingModel::regridModifyState(LevelData<FArrayBox>& a_thickness, 
                                         const AmrIce& a_amrIce,
                                         int a_level)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->regridModifyState(a_thickness, a_amrIce, a_level);
    }

}

//modify a thickness flux
void 
CompositeCalvingModel::modifySurfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                                  LevelData<FArrayBox>& a_calvingMelt,
                                                  const AmrIce& a_amrIce,
                                                  int a_level,
                                                  Real a_dt)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->modifySurfaceThicknessFlux(a_flux, a_calvingMelt, 
                                                  a_amrIce, a_level, a_dt);
                                                  
    }  
}


CompositeCalvingModel::~CompositeCalvingModel()
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      delete m_vectModels[n];
      m_vectModels[n] = NULL;
    }
}


#include "NamespaceFooter.H"
