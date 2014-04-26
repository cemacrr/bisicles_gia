#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef AMR_ICE_H
#define AMR_ICE_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "IceVelocitySolver.H"
#include "RealVect.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "SurfaceFlux.H"
#include "BasalFriction.H"
#include "IceThicknessIBC.H"
#include "IceTemperatureIBC.H"
#include "CalvingModel.H"
#include "MuCoefficient.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "NamespaceHeader.H"

/// class to manage non-subcycled AMR ice-sheet model
/**
 */
class AmrIce
{

public:
  /// Default constructor
  /** At the moment, there is only one constructor, which defines itself
      based on an inputs file.  (This may change, of course)
  */
  AmrIce();

  /// destructor
  /** destructor
   */
  ~AmrIce();

  /// set default values before initialization
  void setDefaults();

  /// initializes object based on inputs data passed in through parmParse
  /** initializes new object based on data in ParmParse database.  This 
      may change to a more explicit initialization if it's thought to be
      a good idea.
  */
  void initialize();

  /// set domain size (in meters)
  void setDomainSize(const RealVect& a_domainSize)
  { m_domainSize = a_domainSize;}
    

  /// set constitutive relation
  /** needs to be done before calling initialize
   */
  void setConstitutiveRelation(ConstitutiveRelation* a_constRelPtr) 
  {m_constitutiveRelation = a_constRelPtr->getNewConstitutiveRelation(); }

   /// set rate factor
  /** needs to be done before calling initialize
   */
  void setRateFactor(RateFactor* a_rateFactorPtr)
  {m_rateFactor = a_rateFactorPtr->getNewRateFactor(); }


   /// set rate factor at base
  /** needs to be done before calling initialize
   */
  void setBasalRateFactor(RateFactor* a_rateFactorPtr)
  {m_basalRateFactor = a_rateFactorPtr->getNewRateFactor(); }
  

  /// set basal Friction relation
  /** needs to be done before calling initialize
   */
  void setBasalFrictionRelation(BasalFrictionRelation* a_basalFrictionRelPtr) 
  {m_basalFrictionRelation = a_basalFrictionRelPtr->getNewBasalFrictionRelation(); }


  /// set BC for thickness advection
  void setThicknessBC( IceThicknessIBC* a_thicknessIBC);
  
  /// set BC for temperature advection
  void setTemperatureBC( IceTemperatureIBC* a_temperatureIBC);

  /// set surfact flux for ice sheet
  void setSurfaceFlux(SurfaceFlux* a_surfaceFluxPtr)
  {
    if (m_surfaceFluxPtr != NULL)
      {
	delete(m_surfaceFluxPtr);
	m_surfaceFluxPtr = NULL;
      }
    m_surfaceFluxPtr = a_surfaceFluxPtr->new_surfaceFlux(); 
  }
  
  /// set basal flux for ice sheet
  void setBasalFlux(SurfaceFlux* a_basalFluxPtr)
  {
    if (m_basalFluxPtr != NULL)
      {
	delete(m_basalFluxPtr);
	m_basalFluxPtr = NULL;
      }
    m_basalFluxPtr = a_basalFluxPtr->new_surfaceFlux(); 
  }

  /// set surface Heat flux for ice sheet
  void setSurfaceHeatBoundaryData(SurfaceFlux* a_surfaceHeatBoundaryDataPtr, bool a_dirichlett)
  {
    if (m_surfaceHeatBoundaryDataPtr != NULL)
      {
	delete(m_surfaceHeatBoundaryDataPtr);
	m_surfaceHeatBoundaryDataPtr = NULL;
      }
    
    //special case for backward compatibility : allow a NULL m_surfaceHeatBoundaryDataPtr
    //with the Dirichlett boundary condition. Data will then be determined by the TemperatureIBC
    CH_assert(a_dirichlett | (a_surfaceHeatBoundaryDataPtr != NULL));
    m_surfaceTemperatureDirichlett = a_dirichlett;
    if (a_surfaceHeatBoundaryDataPtr != NULL)
      {
	m_surfaceHeatBoundaryDataPtr = a_surfaceHeatBoundaryDataPtr->new_surfaceFlux(); 
      }
  }
  
  /// set basal Heat flux for ice sheet
  void setBasalHeatBoundaryData(SurfaceFlux* a_basalHeatBoundaryDataPtr)
  {
    if (m_basalHeatBoundaryDataPtr != NULL)
      {
	delete(m_basalHeatBoundaryDataPtr);
	m_basalHeatBoundaryDataPtr = NULL;
      }
    //special case for backward compatibility : allow a NULL m_basalHeatBoundaryDataPtr and determine the flux from m_temperatureIBCPtr
    if (a_basalHeatBoundaryDataPtr != NULL)
      {
	m_basalHeatBoundaryDataPtr = a_basalHeatBoundaryDataPtr->new_surfaceFlux(); }
  }


  /// set basal friction
  void setBasalFriction(const BasalFriction* a_basalFrictionPtr)
  {m_basalFrictionPtr = a_basalFrictionPtr->new_basalFriction();}
   
  void setMuCoefficient(const MuCoefficient* a_muCoefficientPtr)
  {m_muCoefficientPtr = a_muCoefficientPtr->new_muCoefficient();}

  /// advance solution until either max_time or max_step are reached
  /**
   */
  void run(Real a_max_time, int a_max_step);

  /// compute one timestep
  void timeStep(Real a_dt);

  // compute face-centered half-time ice thickness using ppm
  void computeH_half(Vector<LevelData<FluxBox>* >& a_H_half, Real a_dt);

  // compute thickness fluxes
  void computeThicknessFluxes(Vector<LevelData<FluxBox>* >& a_vectFluxes,
                              const Vector<LevelData<FluxBox>* >& a_H_half,
                              const Vector<LevelData<FluxBox>* >& a_faceVelAdvection);

// update ice thickness
  void updateThickness(Vector<RefCountedPtr<LevelSigmaCS> >& a_vect_coordSys_new, 
                       Vector<RefCountedPtr<LevelSigmaCS> >& a_vectCoords_old, 
                       const Vector<LevelData<FluxBox>* >& a_vectFluxes, 
                       Real a_dt);
    
  /// do regridding
  void regrid();

  /// diagnostic function -- integrates thickness over domain
  Real computeTotalIce() const;


  /// fill flattened Fortran array of data with ice thickness
  void getIceThickness(Real* a_data_ptr, int* a_dim_info, 
		       Real* a_dew, Real* a_dns) const;

  /// is this object defined and initialized?
  /**
   */
  bool isDefined() const;

#ifdef CH_USE_HDF5

  /// write hdf5 plotfile to the standard location
  /** writes plotfile to the standard location; filename is
      <plot_root>.<step>.<DIM>d.hdf5
   */
  void writePlotFile();


  /// write checkpoint file out for later restarting
  /** writes checkpoint data to the standard location; filename is
      <check_root>.<step>.<DIM>d.hdf5
   */
  void writeCheckpointFile() const;

  /// write checkpoint file out for later restarting
  /** writes checkpoint data to the location specified by a_file
   */
  void writeCheckpointFile(const string& a_file) const;

  /// read checkpoint file for restart 
  /** read checkpoint data from file pointed to by a_handle
   */
  void readCheckpointFile(HDF5Handle& a_handle);

  /// set up for restart
  void restart(const string& a_restart_file);

#endif

  /// compute total ice and total grounded ice (volume, at this point)
  Real computeTotalGroundedIce() const;

  /// compute volume over flotation
  Real computeVolumeAboveFlotation() const;

  /// compute total flux over ice
  Real computeFluxOverIce(const Vector<LevelData<FArrayBox>*> a_flux);

  /// compute crevasses depths
  void computeCrevasseDepths(LevelData<FArrayBox>& a_surfaceCrevasse,
			     LevelData<FArrayBox>& a_basalCrevasse,
			     int a_level) const;

  enum velSolverTypes { Picard = 0,
			JFNK = 1,
			KnownVelocity = 2,
			PetscNLSolver = 3,
			FASMGAMR = 4,
                        NUM_SOLVER_TYPES};

  /// read only access members
  Real time() const
  {
    return m_time + m_offsetTime;
  };

  int finestLevel() const
  {
    return m_finest_level;
  }

  const RealVect dx(int a_level) const
  {
    return m_vect_coordSys[a_level]->dx();
  };
  
  const RefCountedPtr<LevelSigmaCS> geometry(int a_level) const 
  {
    return m_vect_coordSys[a_level];
  };
  const DisjointBoxLayout& grids(int a_level) const
  {
    return m_amrGrids[a_level];
  }
  
  const Vector<RefCountedPtr<LevelSigmaCS> >& amrGeometry() const
  {
    return m_vect_coordSys;
  }

  const Vector<int>& refRatios() const { return m_refinement_ratios;}

  const Vector<Real>& amrDx() const {return m_amrDx;}

  const LevelData<FArrayBox>* velocity(int a_level) const
  {
    return m_velocity[a_level];
  };

  const Vector<LevelData<FArrayBox>* >& amrVelocity() const
  {
    return m_velocity;
  }


  const LevelData<FArrayBox>* groundingLineProximity(int a_level) const;
  
  enum viscousTensorComponent{xx,xy,yx,yy,max};

  const LevelData<FArrayBox>* viscousTensor(int a_level) const;
  const LevelData<FArrayBox>* dragCoefficient(int a_level) const;
  const LevelData<FArrayBox>* viscosityCoefficient(int a_level) const;

  /// access function for face-centered velocities
  const Vector<LevelData<FluxBox>* >& faceVelocities() const {return m_faceVelTotal;}

  /// access function for layered face-centered velocities
  const Vector<LevelData<FluxBox>* >& faceVelocitiesLayered() const {return m_layerXYFaceXYVel;} 

  /// access function for calved ice thickness
  const  Vector<LevelData<FArrayBox>* >& calvedIceThickness() const {return m_calvedIceThickness; }
  
  /// access function for surface temperature
  const  Vector<LevelData<FArrayBox>* >& surfaceTemperature() const {return m_sTemperature; }

  /// access function for surface heat flux
  const  Vector<LevelData<FArrayBox>* >& surfaceHeatFlux() const {return m_sHeatFlux; }

protected:
  /// compute tags for regridding
  void tagCells(Vector<IntVectSet>& a_tags);

  /// compute tags for the level a_level
  void tagCellsLevel(IntVectSet& a_tags, int a_level);
 
  /// compute tags at initial time
  void tagCellsInit(Vector<IntVectSet>& a_tags);

  /// initialize grids at initial time
  void initGrids(int a_finest_level);
  
  /// set up grids from grids file
  void setupFixedGrids(const std::string& a_gridFile);

  /// set up storage, etc on an AMR level
  void levelSetup(int a_level, const DisjointBoxLayout& a_grids);
                       
  /// initialize data on hierarchy
  void initData(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCoordSys,
                Vector<LevelData<FArrayBox>* >& a_initialVelocity);


  /// define nonlinear ellipic solver for computing velocity field
  virtual void defineSolver();

  /// solve for velocity field (or just set up some auxilliary quantities)
  virtual void solveVelocityField(Real a_convergenceMetric = -1.0);

  /// compute RHS for velocity field solve
  /** also sets beta to zero where ice is floating
   */
  virtual void defineVelRHS(Vector<LevelData<FArrayBox>* >& a_vectRhs,
                            Vector<LevelData<FArrayBox>* >& a_vectC,
                            Vector<LevelData<FArrayBox>* >& a_vectC0);

  /// set basal friction coefficient (beta) prior to velocity solve
  virtual void setBasalFriction(Vector<LevelData<FArrayBox>* >& a_vectBeta);

  /// set mu coefficient (phi) prior to velocity solve
  virtual void setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef, 
			Vector<LevelData<FluxBox>* >& a_faceMuCoef);

  /// given the current cell centred velocity field, compute a face centred velocity field
  void computeFaceVelocity(Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
			   Vector<LevelData<FluxBox>* >& a_faceVelTotal,
			   Vector<LevelData<FluxBox>* >& a_diffusivity,
			   Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			   Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel);


  /// compute div(vel*H) at a given time
  void computeDivThicknessFlux(Vector<LevelData<FArrayBox>* >& a_divFlux,
                               Vector<LevelData<FluxBox>* >& a_flux,
                               Vector<LevelData<FArrayBox>* >& a_thickness,
                               Real a_time, Real a_dt);

  // increment phi := phi + dt*dphi
  void incrementWithDivFlux(Vector<LevelData<FArrayBox>* >& a_phi,
                            const Vector<LevelData<FArrayBox>* >& a_dphi,
                            Real a_dt);

  // incremente coordSys with new thickness
  void updateCoordSysWithNewThickness(const Vector<LevelData<FArrayBox>* >& a_thickness);

  /// compute timestep
  Real computeDt();

  /// compute timestep at initial time
  Real computeInitialDt();

  //set the calved ice thickness to zero
  void zeroCalvedIceThickness();

  /// implicit solve for diffusive fluxes
  void implicitThicknessCorrection(Real a_dt,
				   const Vector<LevelData<FArrayBox>* >& a_sts,
				   const Vector<LevelData<FArrayBox>* >& a_bts
				   );
  
  void eliminateRemoteIce();

  /// solve a_alpha * phi - a_beta * laplacian(phi) = a_rhs
  /// and grad(phi).n = 0 on the domain boundaries
  void helmholtzSolve(Vector<LevelData<FArrayBox>* >& a_phi,
		      const Vector<LevelData<FArrayBox>* >& a_rhs, 
		      Real a_alpha, Real a_beta) const;

  /// solve a_alpha* phi - a_beta * laplacian(phi)  = phi'
  /// and grad(phi).n = 0 on the domain boundaries
  void helmholtzSolve(Vector<LevelData<FArrayBox>* >& a_phi,
		      Real a_alpha, Real a_beta) const;
 

  // 1D SSA flux across the grounding line from  Schoof (2007) J Geophys. Res. vol 112, p F03S28
  Real glFluxSchoof(Real a_H, Real a_A, Real a_C, Real a_n, Real a_m, Real a_rho, Real a_rhow, Real a_g)
  {
    return std::pow(a_H, (a_n+a_m+3)/(a_m+1.0)) 
      * std::pow( a_A * std::pow(a_rho*a_g,a_n+1.0) 
		  * std::pow(1.0-a_rho/a_rhow,a_n) 
		  / (std::pow(4.0,a_n)*a_C), 1.0/(a_m+1.0));
  }






#if BISICLES_Z == BISICLES_LAYERED
 

  /// compute the flow rate coefficient A given temperature & geometry
  void computeA(Vector<LevelData<FArrayBox>* >& a_A, 
			Vector<LevelData<FArrayBox>* >& a_sA,
			Vector<LevelData<FArrayBox>* >& a_bA,
			const Vector<LevelData<FArrayBox>* >& a_temperature, 
			const Vector<LevelData<FArrayBox>* >& a_sTemperature,
			const Vector<LevelData<FArrayBox>* >& a_bTemperature,
			 const Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys) const;


  //compute the face temperature at t + dt/2
  void computeTHalf(Vector<LevelData<FluxBox>* >& a_layerTH_half, 
		    Vector<LevelData<FluxBox>* >& a_layerH_half,
		    const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel, 
		    const Real a_dt, const Real a_time);

  //update the temperature fields
  void updateTemperature(Vector<LevelData<FluxBox>* >& a_layerTH_half, 
			 Vector<LevelData<FluxBox>* >& a_layerH_half,
			 const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			 const Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel,
			 const Real a_dt, const Real a_time, 
			 Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysNew,
			 Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysOld,
			 const Vector<LevelData<FArrayBox>*>& a_surfaceThicknessSource,
			 const Vector<LevelData<FArrayBox>*>& a_basalThicknessSource);

#endif
  // /// update the surface gradient in level lmin to lmax of a_vectCoordsys
  // void updateSurfaceGradient(Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys,
  // 				    const Vector<DisjointBoxLayout>& a_grids,
  // 				    const Vector<ProblemDomain>&  a_domains,
  // 				    const Vector<int>& a_refinementRatios,
  // 				    const Vector<Real>& a_dx, 
  // 				    const RealVect a_basalSlope,
  // 				    Real a_time, Real a_dt,
  // 				    int a_clipRadius,
  // 				    const int a_lmin, const int a_lmax,
  // 				    IceThicknessIBC* a_thicknessIBCPtr);

  static void postInterpolationReFloat(LevelData<FArrayBox>& a_H,
				       const LevelData<FArrayBox>& a_coarseH,
				       const LevelData<FArrayBox>& a_coarseBed,
				       const DisjointBoxLayout a_newDBL,
				       const ProblemDomain& a_domain,
				       int a_refRatio,
				       Real a_seaLevel, 
				       Real a_waterDensity, 
				       Real a_iceDensity);


#if 0
  /// set up BC's based on parmparse inputs
  void setupBCs();
#endif

  /// temporal accuracy 
  /** 1 = 1st-order (use old-time velocity for flux computation)
      2 = 2nd-order (recompute velocity at half-time for flux computation)
      3 = 3rd-order Runge-Kutta update (if eventually implemented)
      4 = 4th-order Runge-Kutta update
  */
  int m_temporalAccuracy;
  
  /// number of ghost cells we need for m_old_thickness
  /** this is dependent on the advection scheme we're using
      AMRGodunov-based -- 4
      RK4 -- 1
   */
  int m_num_thickness_ghost;

  /// what type of nonlinear solver to use for velocity field
  int m_solverType;

  /// max number of iterations for velocity solver
  int m_maxSolverIterations;

  /// solver tolerance for velSolver.
  Real  m_velocity_solver_tolerance;

  /// solve the velocity problem every m_velocitySolveInterval timesteps
  int m_velocity_solve_interval;

  /// Residual norms at the start and end of the last velocity solve
  Real  m_velocitySolveInitialResidualNorm;
  Real  m_velocitySolveFinalResidualNorm;
  /// nonlinear elliptic solver for velocity field
  IceVelocitySolver* m_velSolver;
  

  /// constitutive relation
  ConstitutiveRelation* m_constitutiveRelation;
  /// rate factor
  RateFactor* m_rateFactor;
  /// basal Friction relation
  BasalFrictionRelation* m_basalFrictionRelation;
  /// rate factor at base
  RateFactor* m_basalRateFactor;

  /// if false, don't reset basal friction for floating ice to zero (default is true)
  /** This option exists primarly so that we can run the CISM circular 
      shelf test problem
  */
  bool  m_reset_floating_friction_to_zero;

  bool m_wallDrag;
  Real m_wallDragExtra;

  // Attempt to limit velocity
  // gradients close to the grounding line by adding a drag
  // term -rHu to the LHS and a source term -rQ to the RHS.
  // r = m_groundingLineRegularization
  // Q is the analytic expression for the GL flux  of 
  // Schoof (2007) J Geophys. Res. vol 112, p F03S28, which
  // is only valid for 1D SSA models.
  Real m_groundingLineRegularization;

  //use one sided diffs at the grounding line?
  bool m_groundingLineCorrection;
#if 0
  /// BCFunc for velSolver. Eventually would like this to be a BCHolder
  BCFunc m_velBCs;
#endif

  /// pointer to a GodunovPhysics for thickness advection
  GodunovPhysics* m_thicknessPhysPtr;

  /// PatchGodunov used for thickness advection -- one for each level
  Vector<PatchGodunov*> m_thicknessPatchGodVect;

  /// IBC for thickness advection
  IceThicknessIBC* m_thicknessIBCPtr;

  /// pointer to a GodunovPhysics for temperature advection
  GodunovPhysics* m_temperaturePhysPtr;

  /// PatchGodunov used for horizontal temperature advection -- one for each level
  Vector<PatchGodunov*> m_temperaturePatchGodVect;
  
  /// IBC for temperature advection
  IceTemperatureIBC* m_temperatureIBCPtr;

  /// Calving model
  CalvingModel* m_calvingModelPtr;

  /// sets surface flux for ice sheet
  SurfaceFlux* m_surfaceFluxPtr;
  /// sets basal flux for ice sheet
  SurfaceFlux* m_basalFluxPtr;
  
   /// sets surface heat flux or temperature for ice sheet
  SurfaceFlux* m_surfaceHeatBoundaryDataPtr;
  /// sets basal heat flux or temperature for ice sheet
  SurfaceFlux* m_basalHeatBoundaryDataPtr;
  

  /// sets basal friction coefficient (beta^2)
  BasalFriction* m_basalFrictionPtr;

  /// sets mu coefficient (phi)
  MuCoefficient* m_muCoefficientPtr;
  
  /// max number ofSflux levels
  int m_max_level;
  
  /// current finest level
  int m_finest_level;
  
  /// if 0 <  m_finest_timestep_level < m_finest_level and >= 0, 
  /// timestep only levels 0 to m_finest_timestep_level, and interpolate
  /// any levels above that. The idea is to avoid being hammered too hard by 
  /// diminishing stability limits on dt
  int m_finest_timestep_level;
  int finestTimestepLevel() const
  {
    return (0 <= m_finest_timestep_level)?
      std::min(m_finest_timestep_level,m_finest_level):m_finest_level;
  }


  /// maximum level at which refinement tags are computed
  int m_tag_cap;

  /// blocking factor
  int m_block_factor;

  /// only allow tagging within a subset of the coarse domain
  IntVectSet m_tag_subset;

  Vector<IntVectSet> m_vectTagSubset;

  /// grid efficiency
  Real m_fill_ratio;

  /// proper nesting radius
  int m_nesting_radius;

  /// max box size
  int m_max_box_size;

  /// max box size for level 0 (defaults to max_box_size)
  /** this allows a different max box size for the level 0 domain,
      which can be useful since the number of boxes on level 0 is often 
      driven by load balancing concerns, while max_box_size on other levels 
      is generally driven by issues of cache size, etc) */
  int m_max_base_grid_size;

  /// regrid interval
  int m_regrid_interval;
  int m_n_regrids;

  /// if true, limit velocity RHS to try to deal with noisy initial conditions
  /** for now, use max(abs(RHS)) = max(Z_s/gradLimitRadius/dx)
   */
  bool m_limitVelRHS;
  int m_gradLimitRadius;

  /// if m_limitFluxSpeed > 0, impose a maximum on the magnitude of  *flux* velocity,
  /// that is, the vertically averaged face-centered velocity which appears in the 
  /// thickness transport equation . 
  Real m_limitFluxSpeed;

  /// if false, don't do initial velocity solve during problem setup. 
  /// (default is true)
  bool m_doInitialVelSolve;

  /// if true, use nonzero initial guess for velocity field (default is false)
  /** 
      if true, make semi-informed initial guess for first pass through 
      ice velocity solver
  */
  bool m_doInitialVelGuess;
  // types of initial guess
  enum initialGuessTypes {SlidingLaw=0, // sliding law with the form Cu = rhs 
			  ConstMu=1, // solve the elliptic equations with a supplied mu
			  Function=2, //use a RealFunction<RealVect> 
			  NUM_INITIAL_GUESS_TYPES};
  int m_initialGuessType;

  //paramamters which specify the linear solver in the case m_initialGuessType = ConstMu
  Real m_initialGuessConstMu; // domain wide viscosity
  int m_initialGuessSolverType; // solver type, defaults to the Picard solver
  RealVect m_initialGuessConstVel; // initial velocity, useful when the basal friction relation is nonlinear

  /// if true, interpolate things like base topography at regrid time 
  /** if false, use IceThicknesIBC::regridIceGeometry function to set values
   */
  bool m_interpolate_zb;


  // if > 0, then the coefficient of friction at the base
  // of the ice and the RHS will be smoothed.
  Real m_basalLengthScale;
 
  /// stores tags from previous regrids
  Vector<IntVectSet> m_vectTags; 

  /// tag on grad(velocity)?
  bool m_tagOnGradVel;
  
  /// tagging value (undivided gradient(vel) threshold for regridding)
  Real m_tagging_val;

  /// tag on laplacian(velocity)
  bool m_tagOnLapVel;

  /// tag on laplacian(velocity) in grounded regions only
  bool m_tagOnGroundedLapVel;

  /// tagging threshold value for undivided lap(vel)
  Real m_laplacian_tagging_val;
  /// when tagging on undivided lap(vel), ignore cells where the basal friction coeffcient exceeds 
  Real m_laplacian_tagging_max_basal_friction_coef;
  /// tag on strain rate invariant?
  bool m_tagOnEpsSqr;

  /// tagging value for strain rate invariant
  /** tag if epsSqr*dx > than this value
   */
  Real m_epsSqr_tagVal;

  /// tag on velocity solve RHS
  bool m_tagOnVelRHS;

  /// tagging value for velocity RHS
  Real m_velRHS_tagVal;

  /// tag on div(H grad (vel)) 
  bool m_tagOndivHgradVel;

  /// tagging threshold value for undivided div(H grad(vel))
  Real m_divHGradVel_tagVal;

  /// tag the grounding line
  bool m_tagGroundingLine;
  /// threshold velocity at the grounding line : tag if |u| > value
  Real m_groundingLineTaggingMinVel;
  /// threshold Basal friction coefficient at the grounding line : tag if C < value
  Real m_groundingLineTaggingMaxBasalFrictionCoef;

  /// tag at the ice margin
  bool m_tagMargin;

  /// tag any cell with ice in it
  bool m_tagAllIce;

  //attempt to eliminate floating ice that is unconnected to grounded ice 
  bool m_eliminate_remote_ice;
  int m_eliminate_remote_ice_max_iter;
  Real m_eliminate_remote_ice_tol;

  //attempt to eliminate unconnected floating ice after regridding
  bool m_eliminate_remote_ice_after_regrid;
  
#ifdef HAVE_PYTHON
  /// tag cells according to a python function
  bool m_tagPython;
  PyObject* m_tagPythonModule;
  PyObject* m_tagPythonFunction;
#endif
  /// amount to buffer tags used in regridding
  // grow tags in all directions by m_tags_grow
  int m_tags_grow; 
  //grow tags in all direction dir by max(0,m_tags_grow_dir[dir] -  m_tags_grow)
  IntVect m_tags_grow_dir; 

  /// refinement ratios
  Vector<int>  m_refinement_ratios;

  /// cell spacing at each level
  Vector<Real> m_amrDx;

  /// problem domains at each level
  Vector<ProblemDomain> m_amrDomains;
  
  // problem domain size
  RealVect m_domainSize;

  /// current grids
  Vector<DisjointBoxLayout> m_amrGrids;

  /// keeps track of which levels are completely covered
  Vector<int> m_covered_level;

  /// book-keeping; keeps track of number of cells per level
  Vector<int> m_num_cells;

  /// current time
  Real m_time; Real m_offsetTime;

  /// most recent timestep
  Real m_dt;

  /// most recent stable timestep (> m_dt if we took a step just to hit a plot interval)
  Real m_stable_dt;

  /// timestep scaling
  Real m_cfl;

  /// cfl number for initial timestep (useful if initial data needs small cfl)
  Real m_initial_cfl;

  /// maximum amount cfl number may grow in one timestep
  Real m_max_dt_grow;

  // current step
  int m_cur_step;

  // isothermal model ?
  bool m_isothermal;
  
  // isothermal temperature
  Real m_isothermalTemperature;

  //depth of water in surface crevasses
  Real m_waterDepth;

  /// should a Dirichlett condition be imposed at the upper surface?
  bool m_surfaceTemperatureDirichlett;

  //ice density
  Real m_iceDensity;

  //sea water density
  Real m_seaWaterDensity;

 //sea water density
  Real m_gravity;

  /// current old-time data
  Vector<LevelData<FArrayBox>*> m_old_thickness;

#if 0  // I think I can drop this one, to be honest...
  /// new-time data
  Vector<LevelData<FArrayBox>*> m_new_thickness;
#endif

  /// horizontal velocity field (found by an IceVelocitySolver)
  Vector<LevelData<FArrayBox>*> m_velocity;

  /// velocity solve RHS (in case we want to save it)
  Vector<LevelData<FArrayBox>* > m_velRHS;

  // coeffient of |u|^m-1 u in the basal traction
  Vector<LevelData<FArrayBox>* > m_velBasalC; 

  //coefficient of mu in the stress tensor
  Vector<LevelData<FArrayBox>* > m_cellMuCoef;
  Vector<LevelData<FluxBox>* > m_faceMuCoef;

  // thickness diffusion coefficient
  Vector<LevelData<FluxBox>* > m_diffusivity;
  Real m_additionalDiffusivity;
  
  //face velocities
  Vector<LevelData<FluxBox>* > m_faceVelAdvection;
  Vector<LevelData<FluxBox>* > m_faceVelTotal; // m_faceVelTotal = m_faceVelAdvection - m_diffusivity * grad(H)
#if BISICLES_Z == BISICLES_LAYERED
  Vector<LevelData<FArrayBox>* > m_layerSFaceXYVel;
  Vector<LevelData<FluxBox>* > m_layerXYFaceXYVel;
#endif

 

  // SigmaCS mapping
  Vector<RefCountedPtr<LevelSigmaCS > > m_vect_coordSys;
  

  // Cell-centered temperature field
  Vector<LevelData<FArrayBox>* > m_temperature;
  // Cell-centered A field
  Vector<LevelData<FArrayBox>* > m_A;
  // used to decide whether A needs updating (after a regrid, or a temperature update)
  mutable bool m_A_valid;
#if BISICLES_Z == BISICLES_LAYERED
  // Temperature and A at the top and bottom surfaces of 
  // the ice, horizontally cell centred (but at vertical faces)
  Vector<LevelData<FArrayBox>* > m_sTemperature;
  Vector<LevelData<FArrayBox>* > m_bTemperature;
  Vector<LevelData<FArrayBox>* > m_sHeatFlux;
  Vector<LevelData<FArrayBox>* > m_bHeatFlux;
  Vector<LevelData<FArrayBox>* > m_sA;
  Vector<LevelData<FArrayBox>* > m_bA;
#elif BISICLES_Z == BISICLES_FULL
  // possibly won't even need surface and basal temperature
  // field, but if we do they will need to be 2D
#endif
  //the thickness source will have two components, from the surface and base.
  //Although they are combined for thickness evolution, they must
  //be kept seperate for the temperature calculation. 
  Vector<LevelData<FArrayBox>* > m_surfaceThicknessSource;
  Vector<LevelData<FArrayBox>* > m_basalThicknessSource;
  
  Vector<LevelData<FArrayBox>* > m_balance; // set to -dH/dt
  //(accumulated) thickness of calved ice 
  Vector<LevelData<FArrayBox>* > m_calvedIceThickness;

  //cell centered storage for the stress-balance equation coefficients
  mutable Vector<LevelData<FArrayBox>* > m_dragCoef;
  mutable Vector<LevelData<FArrayBox>* > m_viscosityCoefCell;

  //cell centered storage for the viscous tensor components
  mutable Vector<LevelData<FArrayBox>* > m_viscousTensorCell;
  //face centered storage for the viscous tensor components
  mutable Vector<LevelData<FluxBox>* > m_viscousTensorFace;
  // cache validity flag
  mutable bool m_viscousTensor_valid;
  void updateViscousTensor() const;

 
  

  //Scalar field that varies in the ice shelf from 0 at the grounding line
  mutable Vector<LevelData<FArrayBox>* > m_groundingLineProximity;
  // cache validity flag;
  mutable bool m_groundingLineProximity_valid;
  // scale length for grounding line proximity
  Real m_groundingLineProximityScale;
  int m_groundingLineProximityCalcType;
  void updateGroundingLineProximity() const;
  
  /// diagnostic quantity -- sum(ice) at initial time
  Real m_initialSumIce;

  /// diagnostic quantity -- sum(ice) at last timestepa
  Real m_lastSumIce;

  /// if true, report sum(grounded ice) as well as sum(total ice)
  bool m_report_grounded_ice;

  /// diagnostic quantity -- sum(grounded ice) at initial time
  Real m_initialSumGroundedIce;

  /// 
  /// diagnostic quantity -- sum(grounded ice) at last timestepa
  Real m_lastSumGroundedIce;

  /// diagnostic quantity -- sum(ice over flotation) at initial time
  Real m_initialVolumeAboveFlotation;

  /// 
  /// diagnostic quantity -- sum(ice over flotation) at last timestepa
  Real m_lastVolumeAboveFlotation;
 
  /// if true, report sum(basal flux) and sum(surface flux)
  bool m_report_total_flux;

  /// type of basal friction distribution
  int m_beta_type;

  /// background value for beta
  Real m_betaVal;

  /// small parameter for use in computing beta, if needed
  Real m_betaEps;

  /// background slope of bottom topography
  /** for various reasons (mostly having to do with periodicity)
      it's simpler to just keep this as a parameter here in AmrIce, 
      where it's added to the grad(z_surf) when it comes time to 
      compute the RHS for the velocity solve, rather than include it 
      in the sigmaCS definition.

      The way this is implemented, the basal height at a realVect
      location X is defined as 
      
      z_b = sigmaCS.getTopography() +basalSlope[0]*X[0] +basalSlope[1]*X[1]
  */
  RealVect m_basalSlope;


  /// how verbose should we be?
  static int s_verbosity;

  /// is this object initialized?
  bool m_is_defined;
  
  /// 
  bool m_do_restart;

  /// if starting from a restart, timestep of restart
  int m_restart_step;

  /// 
  string m_plot_prefix;
  
  ///
  string m_check_prefix;

  bool m_check_overwrite; // overwrite check point files
  bool m_check_exit; // exit immediately after writing checkpoint?


  /// 
  int m_plot_interval;  
  /// write a plot file every m_plot_time_interval years (rather than steps)
  Real m_plot_time_interval;


  /// if true, include mask in plotfiles
  bool m_write_mask;

  /// if true, include dH/dt in plotfiles
  bool m_write_dHDt;

  // if true, include (cell-averaged) flux velocity in plotfiles
  bool m_write_fluxVel;

  // if true, include (cell-averaged) viscous tensor components in plotfiles
  bool m_write_viscousTensor;

  // if true, include base velocity in plotfiles
  bool m_write_baseVel;

  /// if true, include velocity solver RHS in plotfiles
  bool m_write_solver_rhs;

  //if true, include temperature in plot files
  bool m_write_temperature;

  //if true, include thickness source terms in plot files
  bool m_write_thickness_sources;

  /// if true, write out .map.hdf5 file
  bool m_write_map_file;
  
  /// if true, write out plotfile before doing velocity solves
  /** this can be convenient when working with the solverBenchmark
   */
  bool m_write_presolve_plotfiles;

  /// 
  int m_check_interval;

  // if true, allow thickness to change in time
  bool m_evolve_thickness; 

  //if true, allow velocity field to change in time
  bool m_evolve_velocity;

  //if true, keep floating or grounded ice stable;
  bool m_grounded_ice_stable ;
  bool m_floating_ice_stable ;

  enum DiffusionTreatment {NONE,IMPLICIT,EXPLICIT};
  DiffusionTreatment m_diffusionTreatment;
 
  // if true, time step are always an integer power of two
  bool m_timeStepTicks;

  //ised to store a fixed time step > 0 , if desired
  Real m_fixed_dt;

  // if true, report sum of grounded ice (volume) -- default is false
  bool m_reportGroundedIce;

  



#if BISICLES_Z == BISICLES_LAYERED
  //additional members used in poor-mans multidim mode
public:
  const Vector<Real>& getFaceSigma() const
  {
    return m_faceSigma;
  }
protected:
  void setLayers(const Vector<Real>& a_sigma)
  {
    if (m_sigmaSet)
      MayDay::Error("AmrIce::SetLayers should only be called once");
    m_faceSigma = a_sigma;
    m_nLayers = m_faceSigma.size()-1;
    m_sigmaSet = true;
  }

  bool m_sigmaSet;
  unsigned int m_nLayers;
  Vector<Real> m_faceSigma;
  //if true, include layer velocities in plot files
  bool m_write_layer_velocities;
  bool m_additionalVelocity;

#endif
private: 
  template <class T>
  void levelAllocate(LevelData<T>** a_ptrPtr, const DisjointBoxLayout& a_grids, 
		     int a_nComp, const IntVect& a_ghost)
  {
    CH_assert(a_ptrPtr != NULL);
    if (*a_ptrPtr != NULL)
      delete *a_ptrPtr;
    *a_ptrPtr = new LevelData<T>(a_grids, a_nComp, a_ghost);
  }

};  
  

#include "NamespaceFooter.H"

#endif


