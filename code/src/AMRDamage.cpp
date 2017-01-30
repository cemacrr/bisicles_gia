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
#include "FineInterp.H"
#include "CoarseAverageFace.H"
#include "AdvectPhysics.H"
#include "PatchGodunov.H"
#include "PiecewiseLinearFillPatch.H"
#include "DivergenceF_F.H"
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

const AMRDamage& DamageIceObserver::damage() const
{
  return *m_damagePtr;
}

void DamageIceObserver::notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce)
{

  pout() <<  "DamageIceObserver::notify" << std::endl;

  if (a_n == AmrIce::Observer::PostGeometryUpdate)
    {
      m_damagePtr->define(a_amrIce.grids(), a_amrIce.refRatios(), a_amrIce.finestLevel(), 
			  a_amrIce.dx(0));
      const Vector<LevelData<FluxBox>* >& faceVel = a_amrIce.faceVelocities();
      const Vector<LevelData<FArrayBox>* >& cellVel = a_amrIce.amrVelocity();
      m_damagePtr->timestep(a_amrIce.dt(), cellVel, faceVel);

    }


}

AMRDamage::AMRDamage()
{
  m_time_step = 0;
  m_time = 0.0;
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
      m_damage[lev] = new LevelData<FArrayBox>(m_grids[lev],1,4 * IntVect::Unit);
    }
  if (prevDamage.size() == 0)
    {
      // for now, make up some stripes. needs to be something sensible.
      for (int lev = 0; lev <= m_finestLevel; lev++)
	{
	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	    {
	      FArrayBox& damage = (*m_damage[lev])[dit];
	      for (BoxIterator bit(damage.box());bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  Real x = (Real(iv[0]) + 0.5)*m_dx[lev][0];
		  damage(iv) = sin(2.0*M_PI*x/(40.0e+3));
		}
	    }
	}
    }

  if (prevDamage.size() > 0)
    {
      // if previous damage data exists, interpolate onto the new mesh
      for (int lev = 0; lev <= m_finestLevel; lev++)
	{
	  if (lev > 0)
	    {
	      FineInterp fi(m_grids[lev], m_damage[lev]->nComp(), 
			    m_ratio[lev-1], m_grids[lev].physDomain());
	      fi.interpToFine(*m_damage[lev], *m_damage[lev-1]);
	    }
	  Interval ival(0,m_damage[lev]->nComp()-1);
	  prevDamage[lev]->copyTo(ival,  *m_damage[lev], ival);
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

}

void AMRDamage::timestep(Real a_dt, 
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
      

  computeSource(a_dt);
  
  Vector<LevelData<FluxBox>* >  faceFlux;
  faceFlux.resize(m_finestLevel + 1);
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      faceFlux[lev] = new LevelData<FluxBox>(m_grids[lev],1,IntVect::Unit);
    }
  
  computeFlux(faceFlux, m_damage, a_cellVel, a_faceVel, a_dt );
 
  updateDamage(m_damage, faceFlux, a_dt);

  //free flux data
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (faceFlux[lev] != NULL)
	{
	  delete faceFlux[lev]; faceFlux[lev] = NULL;
	}
    }


  //just for now, some output
  {
    Vector<std::string> name(1,"damage");
    std::stringstream ss;
    ss << "damage.";
    ss.width(6);
    ss.fill('0');
    ss << m_time_step;
    std::string filename(ss.str());
    filename.append(".2d.hdf5");
    WriteAMRHierarchyHDF5(filename, m_grids, m_damage, name, 
			  m_grids[0].physDomain().domainBox(), 
			  m_dx[0][0], a_dt, m_time, m_ratio, 
			  m_finestLevel + 1);
  }
  m_time_step++;
  m_time += a_dt;

}

void AMRDamage::computeSource(Real a_dt)
{

}

void AMRDamage::computeFlux
(Vector<LevelData<FluxBox>* >& a_faceFlux,
 const Vector<LevelData<FArrayBox>* >& a_damage,
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
 Real a_dt)
{

  //compute damage = damage - div(faceFlux)*a_dt (+ source*dt);
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       {
	  FArrayBox& damage = (*a_damage[lev])[dit]; 
	  const FluxBox& flux =  (*a_faceFlux[lev])[dit];
	  
	  FArrayBox div(m_grids[lev][dit],1);
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
	  damage += div;
	  
       }
    }
 
}


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
