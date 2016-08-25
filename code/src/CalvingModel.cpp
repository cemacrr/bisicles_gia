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
#include "CrevasseCalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

void 
DeglaciationCalvingModelA::applyCriterion
(LevelData<FArrayBox>& a_thickness, 
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce, 
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();

      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
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
	      
	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}

void DomainEdgeCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
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
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;

		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
					  a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

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
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;

		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
				      a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelMask[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  thck(iv) = std::max(thck(iv),0.0);

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}

    } // end loop over boxes

}

void ProximityCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
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
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& prox = proximity[dit];
	  const FArrayBox& vel = velocity[dit];
	  Box b = thck.box();b &= prox.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real prevThck = thck(iv);
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

	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}

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
      ptr = new BennCalvingModel(pp);
    }
  else if (type == "VanDerVeenCalvingModel")
    {
      ptr = new VdVCalvingModel(pp);
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
DeglaciationCalvingModelB::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
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

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}


void 
ThicknessCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& iceFrac = a_iceFrac[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox effectiveThickness(thck.box(), 1);
      effectiveThickness.copy(thck);
      Box b = thck.box();
      b &= iceFrac.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();          
          // if iceFrac > 0, then rescale effectiveThickness
          // by dividing by iceFrac value, which gives "actual" thickness
          // in the partial cell. Probably eventually want to move this to 
          // fortran
	  Real prevThck = thck(iv);
          if (iceFrac(iv,0) > 0.0)
            {
              effectiveThickness(iv,0) /= iceFrac(iv,0);
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
              // also set the iceFrac to zero in these cells
	      thck(iv) = m_minThickness; 
              iceFrac(iv,0) = 0.0;
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



  
//alter the thickness field at the end of a time step
void
MaximumExtentCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const Real dx = a_amrIce.amrDx()[a_level];
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
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

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

        } // end loop over cells
  
    }

}




//alter the thickness field at the end of a time step
void 
CompositeCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce, 
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
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

void FlotationCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  m_domainEdgeCalvingModel.applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  Real extraMelt = 0.0;
	  if (mask(iv) == FLOATINGMASKVAL)
	    {
	      thck(iv) = 0.0; 
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}

void
CalvingModel::updateCalvedIce(const Real& a_thck, const Real a_prevThck, const int a_mask, Real& a_added, Real& a_calved, Real& a_removed)
{

  if (a_thck > a_prevThck)
    {
      a_added += (a_prevThck-a_thck);
    }
  else 
    {
      if ((a_mask == OPENSEAMASKVAL) || (a_mask == FLOATINGMASKVAL))
	{
	  a_calved += (a_prevThck-a_thck);
	}
      else
	{
	  a_removed += (a_prevThck-a_thck);
	}
    } 

}



#include "NamespaceFooter.H"
