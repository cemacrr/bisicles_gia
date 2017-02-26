#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#include "NyeCrevasseConstitutiveRelation.H"
#include "NyeCrevasseF_F.H"
#include "SigmaCSF_F.H"
#include "LevelMappedDerivatives.H"
#include "NamespaceHeader.H"

//NyeCrevasseConstitutiveRelation::NyeCrevasseConstitutiveRelation(ConstitutiveRe//lation* a_ptr, const Real& a_NyeWaterDepth, const Real& a_NyeA)
//{
//  m_uncrevassedConstitutiveRelation = a_ptr->getNewConstitutiveRelation();
//}

NyeCrevasseConstitutiveRelation::~NyeCrevasseConstitutiveRelation()
{
  if ( m_uncrevassedConstitutiveRelation != NULL)
    {
      delete m_uncrevassedConstitutiveRelation;
      m_uncrevassedConstitutiveRelation = NULL;
    }
}

ConstitutiveRelation* 
NyeCrevasseConstitutiveRelation::getNewConstitutiveRelation() const
{
  NyeCrevasseConstitutiveRelation* newPtr = new NyeCrevasseConstitutiveRelation(m_uncrevassedConstitutiveRelation,m_NyeWaterDepth,m_NyeA);

  return static_cast<ConstitutiveRelation*>(newPtr);
}

void
NyeCrevasseConstitutiveRelation::computeMu(LevelData<FArrayBox>& a_mu,
                                    const LevelData<FArrayBox>& a_vel, 
                                    const LevelData<FArrayBox>* a_crseVelPtr,
                                    int a_nRefCrse,
                                    const LevelData<FArrayBox>& a_A,
                                    const LevelSigmaCS& a_coordSys,
				    const ProblemDomain& a_domain,
                                    const IntVect& a_ghostVect) const
{

  CH_TIME("NyeCrevasseConstitutiveRelation::computeMu");
  MayDay::Error("NyeCrevasseConstitutiveRelation::computeMu");
} 

void
NyeCrevasseConstitutiveRelation::computeDissipation(LevelData<FArrayBox>& a_dissipation,
					     const LevelData<FArrayBox>& a_vel, 
					     const LevelData<FArrayBox>* a_crseVelPtr,
					     int a_nRefCrse,
					     const LevelData<FArrayBox>& a_A,
					     const LevelSigmaCS& a_coordSys,
					     const ProblemDomain& a_domain,
					     const IntVect& a_ghostVect) const
{
  CH_TIME("NyeCrevasseConstitutiveRelation::computeDissipation");
  MayDay::Error("NyeCrevasseConstitutiveRelation::computeDissipation");
}

void  
NyeCrevasseConstitutiveRelation::computeFaceMu(LevelData<FluxBox>& a_mu,
                                        LevelData<FArrayBox>& a_velocity, 
                                        const LevelData<FArrayBox>* a_crseVelPtr,
                                        int a_nRefCrse,
                                        const LevelData<FluxBox>& a_A,
                                        const LevelSigmaCS& a_coordSys,
					const ProblemDomain& a_domain,
                                        const IntVect& a_ghostVect) const
{
  CH_TIME("NyeCrevasseConstitutiveRelation::computeFaceMu");
  
  //All components of the velocity gradient;
  Interval DerivDir(0,SpaceDim-1);
  Interval DerivComps(0,SpaceDim-1);
  const DisjointBoxLayout& grids = a_mu.disjointBoxLayout();
  LevelData<FluxBox> epsSqr(grids, 1, a_ghostVect);
  LevelData<FluxBox> gradU(grids, SpaceDim*SpaceDim, a_ghostVect);
  LevelData<FluxBox> faceTopography(grids, 1, a_ghostVect);

  //we need face-centered topography to compute face-centered thickness above flotation
  CellToEdge(a_coordSys.getTopography(),  faceTopography);

  computeStrainRateInvariantFace
    (epsSqr,gradU,a_velocity, a_crseVelPtr, a_nRefCrse, a_coordSys,
     a_ghostVect);
  
  m_uncrevassedConstitutiveRelation->computeFaceMu
    ( a_mu, a_velocity,  a_crseVelPtr, a_nRefCrse, a_A, a_coordSys, 
      a_domain, a_ghostVect);
    

  for (DataIterator dit(grids);dit.ok();++dit)
    {
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
	{
	  Box box = grids[dit];
	  box.surroundingNodes(faceDir);
	  int dudxComp = derivComponent(0,0);
	  int dudyComp = derivComponent(1,0);
	  int dvdxComp = derivComponent(0,1);
	  int dvdyComp = derivComponent(1,1);
	  
	  //strain rate
	  FArrayBox e(box, SpaceDim*SpaceDim);
	  FORT_UPLUSUT(CHF_FRA(e), 
		       CHF_CONST_FRA(gradU[dit][faceDir]),
		       CHF_INT(dudxComp),
		       CHF_INT(dudyComp),
		       CHF_INT(dvdxComp),
		       CHF_INT(dvdyComp),
		       CHF_BOX(box));
	  
	  //principal strain rate components
	  FArrayBox lambda(box, SpaceDim);
	  FORT_SYMTEIGEN(CHF_FRA(lambda), 
			 CHF_CONST_FRA(e),
			 CHF_INT(dudxComp),
			 CHF_INT(dudyComp),
			 CHF_INT(dvdxComp),
			 CHF_INT(dvdyComp),
			 CHF_BOX(box));

	 

	  //principal stress components
	  lambda.mult(a_mu[dit][faceDir],0,0,1); // only need the first principal stress
	 
	  const FArrayBox& thck = a_coordSys.getFaceH()[dit][faceDir];
	  lambda.mult(thck,0,0,1);

	  //crevasse depth
	  FArrayBox depth(box, 1);
	  Real rhoi =  a_coordSys.iceDensity();
	  Real rhoo =  a_coordSys.waterDensity();
	  Real sea =  a_coordSys.seaLevel();
	  Real gravity =  a_coordSys.gravity();

	  //todo : these all need to be class args of some sort]
	  Real waterDepth = m_NyeWaterDepth; // water depth
	  FArrayBox depth0(box, 1); depth0.setVal(waterDepth);
	  Real a = m_NyeA; // weight for old crevasses / water filling of crevasses
	  Real b = 1.0e0; // weight for nye crevasses
	  Real eps = 1.0e-10; // minimum denominator
	 
	  FArrayBox thckab(box,1);
	  FORT_THICKNESSOVERFLOTATION(CHF_FRA1(thckab,0),
				      CHF_CONST_FRA1(thck,0),
				      CHF_CONST_FRA1(faceTopography[dit][faceDir],0),
				      CHF_CONST_REAL(rhoi),
				      CHF_CONST_REAL(rhoo),
				      CHF_CONST_REAL(sea),
				      CHF_BOX(box));


	  FORT_CREVASSEDEPTH(CHF_FRA1(depth,0), 
			     CHF_CONST_FRA1(depth0,0),
			     CHF_CONST_FRA1(thck,0),
			     CHF_CONST_FRA1(lambda,0),
			     CHF_CONST_FRA1(thckab,0),
			     CHF_CONST_REAL(rhoi),
			     CHF_CONST_REAL(rhoo),
			     CHF_CONST_REAL(gravity),
			     CHF_CONST_REAL(a),
			     CHF_CONST_REAL(b),
			     CHF_CONST_REAL(eps),
			     CHF_BOX(box));

	  //h * mu *= (1-d/h)
	  FORT_CREVASSEMU(CHF_FRA1(a_mu[dit][faceDir],0),
			  CHF_CONST_FRA1(thck,0),
			  CHF_CONST_FRA1(depth,0),
			  CHF_CONST_REAL(eps),
			  CHF_BOX(box));

	} //end loop of face directions
    } //end loop over boxes

}




// void
// NyeCrevasseConstitutiveRelation::modifyTransportCoefficients
// (const LevelData<FArrayBox>& a_cellVel,
//  const LevelData<FArrayBox>* a_crseVelPtr,
//  const LevelData<FArrayBox>* a_crseDiffusivityPtr,
//  int a_nRefCrse,
//  const LevelSigmaCS& a_coordSys,
//  const DisjointBoxLayout& a_grids,
//  const ProblemDomain& a_domain,
//  const LevelData<FArrayBox>& a_A,
//  const LevelData<FArrayBox>& a_sA,
//  const LevelData<FArrayBox>& a_bA,
//  LevelData<FluxBox>& a_faceVelAdvection,
//  LevelData<FluxBox>& a_faceVelTotal,
//  LevelData<FluxBox>& a_faceDiffusivity,
//  LevelData<FArrayBox>& a_cellDiffusivity,
//  LevelData<FluxBox>& a_layerXYFaceXYVel,
//  LevelData<FArrayBox>& a_layerSFaceXYVel) const
// {
//   CH_TIME("NyeCrevasseConstitutiveRelation::modifyTransportCoefficients");
//   MayDay::Error("NyeCrevasseConstitutiveRelation::modifyTransportCoefficients not implemented");
// }


#include "NamespaceFooter.H"
