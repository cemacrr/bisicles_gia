#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NONLINEARVISCOUSTENSOR_H_
#define _NONLINEARVISCOUSTENSOR_H_


#include "IceVelocitySolver.H"
#include "AMRMultiGrid.H"
#include "MultilevelIceVelOp.H"
#include "ViscousTensorOp.H"
#ifdef CH_USE_PETSC
#include "PetscAMRSolver.H"

#endif

#include "NamespaceHeader.H"

/** virtual base for classes which encapsulate a nonlinear viscosity tensor. 
    States are defined by the velocity, and the important dependents are the 
    coefficients in the linear ViscousTensorOp: alpha, mu, lambda. 
*/ 
class NonlinearViscousTensor
{
public:
  virtual ~NonlinearViscousTensor(){};
  //factory method
  virtual NonlinearViscousTensor* newNonlinearViscousTensor()=0;
  // set the state 
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_state)=0;
  // fill ghost cells of a_state as per the boundary conditions
  virtual void applyBC(Vector<LevelData<FArrayBox>*>& a_phi)=0;
  // retrieve the state 
  virtual const Vector<LevelData<FArrayBox>*>& getState() = 0;
  //retrieve the coefficient alpha
  virtual Vector<RefCountedPtr< LevelData<FArrayBox> > >& alpha() =0;
  //retrieve the coefficient mu 
  virtual Vector<RefCountedPtr< LevelData<FluxBox> > >& mu()  =0;
  //retrieve the coefficient lambda
  virtual Vector<RefCountedPtr< LevelData<FluxBox> > >& lambda() =0;
  // provide a pointer to an AMRLevelOpFactory<LevelData<FArrayBox> >
  // defined at the current state
  virtual RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > opFactoryPtr()=0;
};

//manage coeffients mu(u), lambda(u) = 2*mu(u) and alpha(u) 
//in L[u] = div( mu(theta) * (grad(u) + grad(u)^T) + lambda * div(u)*I)  - alpha*u
//and provide the appropriate ViscousTensorOpFactory
//alpha(u) and mu(u),lambda(u) are determined by BasalFrictionRelation and 
//ConstitutiveRelation members.
class IceNonlinearViscousTensor : public NonlinearViscousTensor
{
 
  Vector<LevelData<FArrayBox>*> m_u;
  Vector<LevelData<FArrayBox>*> m_C;
  Vector<LevelData<FArrayBox>*> m_C0;
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > m_opFactoryPtr;
  
  //mesh parameters
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_refRatio;
  Vector<ProblemDomain> m_domains;
  Vector<RealVect> m_dxs;
  int m_finestLevel;
  Vector<RefCountedPtr<LevelSigmaCS > > m_coordSys;
  
  //constitutive relations
  const ConstitutiveRelation* m_constRelPtr;
  const BasalFrictionRelation* m_basalFrictionRelPtr;
  
  //boundary conditions
  IceThicknessIBC* m_bcPtr;

  //pointers to temperature fields
  Vector<LevelData<FArrayBox>*> m_A;
  Vector<LevelData<FluxBox>*> m_faceA;

  //additional face-centred coefficient of H*mu, usually NULL
  Vector<LevelData<FluxBox>*> m_muCoef;

  // time : this class doesn't care, but m_bcPtr does
  Real m_time;

  // safety factor for VicousTensorOp
  Real m_vtopSafety;
  //  minimum number of smnooths in viscous tensor op mg relax
  int m_vtopRelaxMinIter;
   // relax tolerance for  viscous tensor op
  Real m_vtopRelaxTol;

  //lower and upper bounds on the viscosity
  Real m_muMin, m_muMax;
  
  // BCHolder
  RefCountedPtr<BCFunction> m_velSolveBC;

  //coefficient storage. 
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_mu;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_lambda;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_alpha;

  void setupCoeffs();

public:

  ~IceNonlinearViscousTensor();

  virtual NonlinearViscousTensor* newNonlinearViscousTensor();

  IceNonlinearViscousTensor(const IceNonlinearViscousTensor&);

  IceNonlinearViscousTensor(const Vector<DisjointBoxLayout>& a_grids,
	       const Vector<int>& a_refRatio,
	       const Vector<ProblemDomain>& a_domains,
	       const Vector<RealVect>& a_dxs,
	       const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
	       const Vector<LevelData<FArrayBox>*>& a_u,
	       const Vector<LevelData<FArrayBox>*>& a_C,
	       const Vector<LevelData<FArrayBox>*>& a_C0,
	       const int a_finestLevel,
	       const ConstitutiveRelation& a_constRel,
	       const BasalFrictionRelation& a_basalFrictionRel,
	       IceThicknessIBC& a_bc,
	       const Vector<LevelData<FArrayBox>*>& a_A,
	       const Vector<LevelData<FluxBox>*>& a_faceFlowLawCoefficient,
	       const Real a_time, 
	       const Real a_vtopSafety,
	       const int a_vtopRelaxMinIter,
	       const Real a_vtopRelaxTol,
	       const Real a_muMin = 0.0,
	       const Real a_muMax = 1.234567e+300);

 
  

  //implementation of NonlinearViscousTensor interface 
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_u);

  // fill ghost cells of a_state as per the boundary conditions
  virtual void applyBC(Vector<LevelData<FArrayBox>*>& a_phi)
  {
    for (int lev =0; lev < a_phi.size(); lev++)
      {
	for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	  {
	    (*m_velSolveBC)((*a_phi[lev])[dit], m_grids[lev][dit], m_domains[lev], 
                            m_dxs[lev][0], dit(), false);
	  }
      }
  }


  virtual Vector<LevelData<FArrayBox>*>& getState() 
  {
    return m_u;
  }
  virtual Vector<LevelData<FArrayBox>*>& getC() 
  {
    return m_C;
  }
  virtual Vector<RefCountedPtr<LevelData<FluxBox> > >& mu()
  {
    return m_mu;
  }
  virtual Vector<RefCountedPtr<LevelData<FluxBox> > >& lambda()
  {
    return m_lambda;
  }
  virtual Vector<RefCountedPtr<LevelData<FArrayBox> > >& alpha() 
  {
    return m_alpha;
  }

  virtual RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
  opFactoryPtr() {
    return m_opFactoryPtr;
  } 

  //an overloaded setState, because beta changes with time
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_u, 
			const Vector<LevelData<FArrayBox>*>& a_basalFrictionCoefficent)
  {
    m_C = a_basalFrictionCoefficent;
    setState(a_u);
  }

  //allows an additional coefficient f of the face-centered viscosity to be specified,
  // so that H*mu*(grad(u) ... ) is replaced with  f*H*mu*(grad(u) ... ).
  virtual void setFaceViscCoef(const Vector<LevelData<FluxBox>*>& a_muCoef)
  {
    int nlev = min(a_muCoef.size(),m_muCoef.size());

    for (int lev=0; lev < nlev; lev++)
      {
	m_muCoef[lev] = a_muCoef[lev];
      }
  }
  virtual void unsetFaceViscCoef()
  {
    for (int lev=0; lev < m_muCoef.size(); lev++)
      {
	m_muCoef[lev] = NULL;
      }
  }

  //provide the components of the Visocus Tensor G[ij] at cell faces, such that
  //at j-faces we have G[ij] for all i and so on
  virtual void computeViscousTensorFace(const Vector<LevelData<FluxBox>*>& a_viscousTensor);
    
};

#include "NamespaceFooter.H"
#endif
