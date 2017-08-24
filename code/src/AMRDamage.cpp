#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRDamage.H"
#include "DamageConstitutiveRelation.H"
#include "AmrIce.H"
#include "FineInterp.H"
#include "CoarseAverageFace.H"
#include "AdvectPhysics.H"
#include "PatchGodunov.H"
#include "PiecewiseLinearFillPatch.H"
#include "DivergenceF_F.H"
#include "ParmParse.H"
#include "LevelMappedDerivatives.H"
#include "NyeCrevasseF_F.H"
#include "SigmaCSF_F.H"
#include "CH_HDF5.H"
#include "NamespaceHeader.H"

DamageIceObserver::DamageIceObserver()
{
  m_damagePtr = new AMRDamage();
}

DamageIceObserver::~DamageIceObserver()
{
  if (m_damagePtr != NULL)
    {
      delete m_damagePtr;
      m_damagePtr = NULL;
    }
}

AMRDamage& DamageIceObserver::damage() const
{
  return *m_damagePtr;
}

void DamageIceObserver::notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce)
{

  pout() <<  "DamageIceObserver::notify" << std::endl;

  if (a_n == AmrIce::Observer::PreVelocitySolve)
    {
      // the velocity solve will need the damage field
      m_damagePtr->define(a_amrIce.grids(), a_amrIce.refRatios(), 
			  a_amrIce.finestLevel(), a_amrIce.dx(0));
    }
  else if (a_n == AmrIce::Observer::PostGeometryUpdate)
    {
      //m_damagePtr->define(a_amrIce.grids(), a_amrIce.refRatios(), a_amrIce.finestLevel(), 
      //			  a_amrIce.dx(0));
      const Vector<LevelData<FluxBox>* >& faceVel = a_amrIce.faceVelocities();
      const Vector<LevelData<FArrayBox>* >& cellVel = a_amrIce.amrVelocity();


      Vector<LevelData<FArrayBox>const * > visTensor(a_amrIce.finestLevel() + 1, NULL);
      for (int lev =0; lev <= a_amrIce.finestLevel(); lev++)
	{
	  visTensor[lev] = a_amrIce.viscousTensor(lev);
	}

      Vector<LevelData<FArrayBox>const * > surfThickSource(a_amrIce.finestLevel() + 1, NULL);
      for (int lev =0; lev <= a_amrIce.finestLevel(); lev++)
	{
	  surfThickSource[lev] = a_amrIce.surfaceThicknessSource(lev);
	}

      Vector<LevelData<FArrayBox>const * > baseThickSource(a_amrIce.finestLevel() + 1, NULL);
      for (int lev =0; lev <= a_amrIce.finestLevel(); lev++)
	{
	  baseThickSource[lev] = a_amrIce.basalThicknessSource(lev);
	}

      const Vector<RefCountedPtr<LevelSigmaCS> >& geometry = a_amrIce.amrGeometry();

 
      m_damagePtr->timestep(a_amrIce.dt(), visTensor, geometry, surfThickSource, baseThickSource, cellVel, faceVel);
     
      

    }


}
#ifdef CH_USE_HDF5
void DamageIceObserver::addPlotVars(Vector<std::string>& a_vars)
{
  m_damagePtr->addPlotVars(a_vars);
}

void DamageIceObserver::writePlotData(LevelData<FArrayBox>& a_data, int a_level)
{
  m_damagePtr->writePlotData(a_data, a_level);
}

/// fill a_var with the names of variables to add to the checkpoint file 
void DamageIceObserver::addCheckVars(Vector<std::string>& a_vars)
{
  m_damagePtr->addCheckVars(a_vars);
}
  
/// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
void DamageIceObserver::writeCheckData(HDF5Handle& a_handle, int a_level)
{
  m_damagePtr->writeCheckData(a_handle, a_level);
}
  
/// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
void DamageIceObserver::readCheckData(HDF5Handle& a_handle, HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids)
{
  m_damagePtr->readCheckData(a_handle, a_header, a_level, a_grids);
}
#endif


AMRDamage::~AMRDamage()
{
  
  for (int lev = 0; lev < m_damage.size(); lev++)
    {
      if (m_damage[lev] != NULL)
	{
	  delete m_damage[lev]; m_damage[lev] = NULL;
	}
    }

  for (int lev =0; lev < m_faceMinDamage.size(); lev++)
    {
      if (m_faceMinDamage[lev] != NULL)
	{
	  delete m_faceMinDamage[lev];m_faceMinDamage[lev]=NULL; 
	}
    }
 
}
AMRDamage::AMRDamage()
{
  m_time_step = 0;
  m_time = 0.0;
}

const LevelData<FArrayBox>* AMRDamage::damage(int a_level) const
{
  if (!(m_damage.size() > a_level))
    {
      std::string msg("AMRDamage::damage !(m_damage.size() > a_level)");
      pout() << msg <<endl;
      CH_assert((m_damage.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_damage[a_level];
  if (ptr == NULL)
    {
      std::string msg("AMRDamage::damage m_damage[a_level] == NULL");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }
  return ptr;
}

void AMRDamage::define
(const Vector<DisjointBoxLayout>& a_grids, 
 const Vector<int>& a_ratio,
 int a_finestLevel,
 const RealVect& a_crseDx)
{

  if (m_dx.size() > 0)
    {
      if ( !(m_dx[0] == a_crseDx) )
	{
	  std::string msg("AMRDamage::define, incompatible mesh");
	  pout() << msg << std::endl;
	  CH_assert(m_dx[0] == a_crseDx);
	  MayDay::Error(msg.c_str());
	}
    }

  //update the mesh hierarchy
  m_finestLevel = a_finestLevel;
  m_grids.resize(m_finestLevel  + 1);
  m_ratio.resize(m_finestLevel + 1, 1);
  for (int lev = 0; lev <= m_finestLevel ; lev++)
    {
      m_grids[lev] = a_grids[lev];
      if (lev < m_finestLevel )
	{
	  m_ratio[lev] = a_ratio[lev];
	}
    }
  
  m_dx.resize(m_finestLevel + 1);
  m_dx[0] = a_crseDx;
  for (int lev = 1; lev <= m_finestLevel;  lev++)
    {
      m_dx[lev] = m_dx[lev-1] / Real(m_ratio[lev-1]);
    }


  //copy any previous damage data pointers
  Vector<LevelData<FArrayBox>* > prevDamage;
  prevDamage.resize(m_damage.size());
  for (int lev = 0; lev < m_damage.size(); lev++)
    {
      prevDamage[lev] = m_damage[lev];
    }

  //initialize the damage data
  m_damage.resize(m_finestLevel + 1);
  for (int lev = 0; lev <= m_finestLevel; lev++)
    {
      m_damage[lev] = new LevelData<FArrayBox>(m_grids[lev],  DAMAGE_N_COMP , DAMAGE_N_GHOST * IntVect::Unit);
    }
  if (prevDamage.size() == 0)
    {
     
      // for now, set to zero, but need to set initial conditions from input data of some sort, 
      // or read checkpoints
      for (int lev = 0; lev <= m_finestLevel; lev++)
       	{
       	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       	    {
       	      FArrayBox& damage = (*m_damage[lev])[dit];
	      damage.setVal(0.0);
	      
	    }
	}
    }

  if (prevDamage.size() > 0)
    {
      // if previous damage data exists, interpolate onto the new mesh
      for (int lev = 0; lev <= m_finestLevel; lev++)
	{

	  //FIXME : only doing this to avoid stupid value outside the domain, should be a BC thing....
	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       	    {
       	      FArrayBox& damage = (*m_damage[lev])[dit];
	      damage.setVal(0.0);  
	    }


	  if (lev > 0)
	    {
	      FineInterp fi(m_grids[lev], m_damage[lev]->nComp(), 
			    m_ratio[lev-1], m_grids[lev].physDomain());
	      fi.interpToFine(*m_damage[lev], *m_damage[lev-1]);

	      //need to fill the ghost cells on CF-interfaces now, because *m_damage[lev] has to be
	      //in good shape at all times if the interface with DamageConstitutiveRelation is to work
	      
	      PiecewiseLinearFillPatch ghostFiller(m_grids[lev],m_grids[lev-1],m_damage[lev]->nComp(),
						   m_grids[lev-1].physDomain(), m_ratio[lev-1], m_damage[lev]->ghostVect()[0]);
	      ghostFiller.fillInterp(*m_damage[lev], *m_damage[lev-1],*m_damage[lev-1],1.0,0,0,m_damage[lev]->nComp());

	    }
	  Interval ival(0,m_damage[lev]->nComp()-1);

	  if (prevDamage.size() > lev && prevDamage[lev])
	    prevDamage[lev]->copyTo(ival,  *m_damage[lev], ival);
	  m_damage[lev]->exchange();

	}
      //free old data
      for (int lev =0; lev < prevDamage.size(); lev++)
	{
	  if (prevDamage[lev] != NULL)
	    {
	      delete prevDamage[lev]; prevDamage[lev] = NULL;
	    }	
	}

    }

  //update storage for m_faceMinDamage : I think we don't care what it contains....
   for (int lev =0; lev < m_faceMinDamage.size(); lev++)
     {
       if (m_faceMinDamage[lev] != NULL)
	 {
	   delete m_faceMinDamage[lev];m_faceMinDamage[lev]=NULL; 
	 }
     }
   m_faceMinDamage.resize(m_finestLevel + 1);
    for (int lev =0; lev < m_faceMinDamage.size(); lev++)
     {
       m_faceMinDamage[lev] = new LevelData<FluxBox>(m_grids[lev],1, IntVect::Zero);
     }
  


}

void AMRDamage::timestep(Real a_dt,
			 const Vector<LevelData<FArrayBox>const * >& a_visTensor,
			 const Vector<RefCountedPtr<LevelSigmaCS> >& geometry,
			 const Vector<LevelData<FArrayBox>const * >& a_surfThickSource,
			 const Vector<LevelData<FArrayBox>const * >& a_baseThickSource,
			 const Vector<LevelData<FArrayBox>* >& a_cellVel, 
			 const Vector<LevelData<FluxBox>* >& a_faceVel)
{

  // fill ghost data
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (lev > 0)
	{
	  int nGhost = m_damage[lev]->ghostVect()[0];
	  PiecewiseLinearFillPatch pwl(m_grids[lev],  m_grids[lev-1], 1, 
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost);    
	  // since we're not subcycling, don't need to interpolate in time
	  Real time_interp_coeff = 0.0;
	  pwl.fillInterp(*m_damage[lev], *m_damage[lev-1], *m_damage[lev-1],
			 time_interp_coeff,0, 0, 1);
	}  
      m_damage[lev]->exchange();
    }
      
  Vector<LevelData<FArrayBox>* >  source;
  source.resize(m_finestLevel + 1);
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      source[lev] = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
    }
  computeSource(source, m_damage,geometry, a_visTensor, a_surfThickSource, a_baseThickSource, a_dt);
  
  Vector<LevelData<FluxBox>* >  faceFlux;
  faceFlux.resize(m_finestLevel + 1);
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      faceFlux[lev] = new LevelData<FluxBox>(m_grids[lev],1,IntVect::Unit);
    }
  
  computeFlux(faceFlux, m_damage, source, a_cellVel, a_faceVel, a_dt );

  updateDamage(m_damage, faceFlux, m_faceMinDamage , a_faceVel, source, geometry, a_visTensor,  a_dt);


  //free temporary data.
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (source[lev] != NULL)
	{
	  delete source[lev]; source[lev] = NULL;
	}
    
      if (faceFlux[lev] != NULL)
	{
	  delete faceFlux[lev]; faceFlux[lev] = NULL;
	}

     
    }

  m_time_step++;
  m_time += a_dt;

}

///Compute the source part of equation dD/dt + div(uD) = source.
/**
   There are two obvious source/sink terms : reduction in dh damage removal due to the 
   loss of new surface or basal ice, and damage removal due to time dependent healing.
   For now we have only the first. Damage reduction due to surface/basal accretion
   does not need a source: D will remian constant while h grows, so D/h drops.
 */
void AMRDamage::computeSource(Vector<LevelData<FArrayBox>* >& a_source, 
			      const Vector<LevelData<FArrayBox>* >& a_damage,
			      const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
			      const Vector<LevelData<FArrayBox> const * >& a_visTensor,
			      const Vector<LevelData<FArrayBox>const * >& a_surfThickSource,
			      const Vector<LevelData<FArrayBox>const * >& a_baseThickSource,
			      Real a_dt)
{
 
  for (int lev = 0 ; lev <= m_finestLevel; lev++)
    {
      const LevelSigmaCS& levelCoords = *a_geometry[lev];
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  FArrayBox& source = (*a_source[lev])[dit];
	  FArrayBox& damage = (*a_damage[lev])[dit];
	  const FArrayBox& smb = (*a_surfThickSource[lev])[dit];
	  const FArrayBox& bmb = (*a_baseThickSource[lev])[dit];
	  const FArrayBox& thick = levelCoords.getH()[dit];
	  for (BoxIterator bit(source.box());bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      Real surfDamageSource;
	      surfDamageSource = 0.0;
	      if (smb(iv)<0.0){surfDamageSource = smb(iv);}
	      Real baseDamageSource;
	      baseDamageSource = 0.0;
	      if (bmb(iv)<0.0){baseDamageSource = bmb(iv);}
	      //damage source removed by surface and basal melting
	      source(iv) = (surfDamageSource + baseDamageSource)*damage(iv)/thick(iv);
	    }
	}
    }
}
//compute the half time damage face flux
void AMRDamage::computeFlux
(Vector<LevelData<FluxBox>* >& a_faceFlux,
 const Vector<LevelData<FArrayBox>* >& a_damage,
 const Vector<LevelData<FArrayBox>* >& a_source,
 const Vector<LevelData<FArrayBox>* >& a_cellVel, 
 const Vector<LevelData<FluxBox>* >& a_faceVel,
 Real a_dt)
{

  for (int lev = 0 ; lev <= m_finestLevel; lev++)
    { 
      Real dx = m_dx[lev][0];
      
      AdvectPhysics advectPhys;
      advectPhys.define( m_grids[lev].physDomain(),dx );
      DamagePhysIBC* ptr = new DamagePhysIBC;
      advectPhys.setPhysIBC(ptr);
      delete ptr;

      PatchGodunov patchG;
      int normalPredOrder = 2;
      bool useFourthOrderSlopes = true;
      bool usePrimLimiting = true;
      bool useCharLimiting = false;
      bool useFlattening = false;
      bool useArtificialViscosity = false;
      Real artificialViscosity = 0.0;
      patchG.define(m_grids[lev].physDomain(), 
		    dx, &advectPhys, 
		    normalPredOrder,
		    useFourthOrderSlopes,
		    usePrimLimiting,
		    useCharLimiting,
		    useFlattening,
		    useArtificialViscosity,
		    artificialViscosity);

      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  const FArrayBox& cellVel = (*a_cellVel[lev])[dit];
	  const FArrayBox& damage = (*a_damage[lev])[dit]; 
	  const FluxBox& faceVel = (*a_faceVel[lev])[dit];
	  FluxBox& flux = (*a_faceFlux[lev])[dit];

	  patchG.setCurrentBox(m_grids[lev][dit]);
	  patchG.setCurrentTime(0.0);
	  
	  FArrayBox* cv = const_cast<FArrayBox*>(&cellVel);
	  FluxBox* fv = const_cast<FluxBox*>(&faceVel);

	  AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>
	    (patchG.getGodunovPhysicsPtr());
	  advectPhysPtr->setVelocities( cv, fv  );
	  
	  //todo : compute source terms
	  //FArrayBox& source = (*a_source[lev])[dit]; 
	  FArrayBox source(m_grids[lev][dit],1); source.setVal(0.0);

	  FluxBox faceDamage(flux.box(),1);
	  patchG.computeWHalf(faceDamage, damage, source, 
			      a_dt, m_grids[lev][dit]);

	  //compute flux from face-centered values
	  
	  for (int dir=0; dir<SpaceDim; dir++)
            {
	      Box faceBox(m_grids[lev][dit]);
	      faceBox.surroundingNodes(dir);
	      flux[dir].copy(faceDamage[dir], faceBox);
              flux[dir].mult(faceVel[dir], faceBox, 0, 0, 1);
	    }
	    
	}
    }

 


  for (int lev = m_finestLevel ; lev > 0; lev--)
    {
      CoarseAverageFace faceAverager(m_grids[lev],1, m_ratio[lev-1]);
      faceAverager.averageToCoarse(*a_faceFlux[lev-1], *a_faceFlux[lev]);
    }

}



void AMRDamage::updateDamage
(Vector<LevelData<FArrayBox>* >& a_damage,
 const Vector<LevelData<FluxBox>* >& a_faceFlux,
 const Vector<LevelData<FluxBox>* >& a_faceMinDamage,
 const Vector<LevelData<FluxBox>* >& a_faceVel,
 const Vector<LevelData<FArrayBox>* >& a_source,
 const Vector<RefCountedPtr<LevelSigmaCS> >& geometry,
 const Vector<LevelData<FArrayBox>const * >& a_vt,
 Real a_dt)
{
 
  //compute damage = damage - div(faceFlux)*a_dt (+ source*dt);
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      //usual advection / source transport
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       {
	 const Box& box = m_grids[lev][dit];
	 FArrayBox& damage = (*a_damage[lev])[dit]; 
	 const FluxBox& flux =  (*a_faceFlux[lev])[dit];

	 FArrayBox& source = (*a_source[lev])[dit];
       
	 FArrayBox div(box,1);
	 div.setVal(0.0);
	 for (int dir=0; dir<SpaceDim; dir++)
	   {
	     FORT_DIVERGENCE(CHF_CONST_FRA(flux[dir]),
			     CHF_FRA(div),
			     CHF_BOX(m_grids[lev][dit]),
			     CHF_CONST_REAL(m_dx[lev][dir]),
			     CHF_INT(dir));
	   }
	  
	  div *= -a_dt;
	  source *= a_dt;
	  damage += div;
	  damage += source;
	  
	  FArrayBox localDamage(box,1);
	  DamageConstitutiveRelation::computeLocalDamageVT
	    (localDamage, (*a_vt[lev])[dit], geometry[lev]->getH()[dit], 
	     geometry[lev]->getTopography()[dit],  geometry[lev]->iceDensity(),
	     geometry[lev]->waterDensity(), geometry[lev]->gravity(), 
	     geometry[lev]->seaLevel(), box);
	
	  // set D = max(D_tranport,D_local) 
	  FORT_FABMAX(CHF_FRA1(damage, 0), CHF_FRA1(localDamage, 0), CHF_BOX(box));

	  //limit to thickness
	  FArrayBox max_damage(damage.box(), 1);
	  max_damage.copy(geometry[lev]->getH()[dit]);
	  max_damage *= DAMAGE_MAX_VAL;
	  FORT_FABMIN(CHF_FRA1(damage, 0), CHF_FRA1(max_damage, 0), CHF_BOX(damage.box()));

       } // end loop over boxes
    } // end loop over levels
       
}
#ifdef CH_USE_HDF5
void AMRDamage::addPlotVars(Vector<std::string>& a_vars)
{
  a_vars.push_back("VIDamage");
}

void AMRDamage::writePlotData(LevelData<FArrayBox>& a_data, int a_level)
{
  for (DataIterator dit(m_grids[a_level]);dit.ok();++dit)
    {
      a_data[dit].copy( (*m_damage[a_level])[dit],0,0,1);
    }
}



/// fill a_var with the names of variables to add to the checkpoint file 
void AMRDamage::addCheckVars(Vector<std::string>& a_vars)
{
  a_vars.push_back("VIDamage");
}
  
/// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
void AMRDamage::writeCheckData(HDF5Handle& a_handle, int a_level)
{
  write(a_handle, *m_damage[a_level], "DamageData", m_damage[a_level]->ghostVect());
}
  
/// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
void AMRDamage::readCheckData(HDF5Handle& a_handle, HDF5HeaderData&  a_header, int a_level, const DisjointBoxLayout& a_grids)
{
  bool containsDamageData(false);
  map<std::string, std::string>::const_iterator i;
  for (i = a_header.m_string.begin(); i!= a_header.m_string.end(); ++i)
    {
      containsDamageData |= (i->second == "VIDamage");
    }

  if (containsDamageData)
    {

      if (m_damage.size() <= a_level)
	{
	  m_damage.resize(a_level + 1, NULL);
	  if (m_damage[a_level] != NULL)
	    {
	      delete m_damage[a_level];
	    }
	  
	}
      m_damage[a_level] = new LevelData<FArrayBox>(a_grids, DAMAGE_N_COMP, DAMAGE_N_GHOST * IntVect::Unit); 
      int dataStatus =  read<FArrayBox>(a_handle, *m_damage[a_level], "DamageData", a_grids);
      if (dataStatus != 0)
	{
	  MayDay::Error("failed to read damage data from checkpoint, but the header indicated its presence");
	}
    }
}
#endif


void DamagePhysIBC::define(const ProblemDomain& a_domain,
			  const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

PhysIBC* DamagePhysIBC::new_physIBC()
{
  DamagePhysIBC *ptr = new DamagePhysIBC();
  ptr->define(PhysIBC::m_domain, PhysIBC::m_dx);
  return static_cast<PhysIBC*>(ptr);
}

void DamagePhysIBC::initialize(LevelData<FArrayBox>& a_U)
{
  MayDay::Error("DamagePhysIBC::initialize not implemented");
}

void DamagePhysIBC::primBC(FArrayBox&            a_WGdnv,
			   const FArrayBox&      a_Wextrap,
			   const FArrayBox&      a_W,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const Real&           a_time)
{

  //first off, set WGdnv on the domain boundaries to a_Wextrap
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
	  for (BoxIterator bit(boundaryBox); bit.ok(); ++bit)
	    {
	      const IntVect& i = bit();
	      //All boundaries are outflows
	      a_WGdnv(i,0) = std::max(0.0,a_Wextrap(i,0));
	    }
	}
    }
 
}

void DamagePhysIBC::setBdrySlopes(FArrayBox&       a_dW,
				 const FArrayBox& a_W,
				 const int&       a_dir,
				 const Real&      a_time)
{
  //one-sided differences are fine, so do nothing
}

void DamagePhysIBC:: artViscBC(FArrayBox&       a_F,
		 const FArrayBox& a_U,
		 const FArrayBox& a_divVel,
		 const int&       a_dir,
		 const Real&      a_time)
{
  MayDay::Error("DamagePhysIBC:artViscBC not implemented");
}
