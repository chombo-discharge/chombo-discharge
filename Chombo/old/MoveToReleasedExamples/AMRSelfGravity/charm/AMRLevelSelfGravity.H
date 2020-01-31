#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELSELFGRAVITY_H_
#define _AMRLEVELSELFGRAVITY_H_

#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"

#include "FArrayBox.H"
#include "AMRLevel.H"
#include "LevelData.H"
#include "LevelGodunov.H"

#include "CoarseAverage.H"
#include "FineInterp.H"
#include "QuadCFInterp.H"
#include "LevelFluxRegister.H"
#include "RefCellTagger.H"

#include "Gradient.H"
#include "SelfGravityPhysics.H"

#include "CH_HDF5.H"
#include "CONSTANTS_C.H"

/// AMRLevel-derived class for AMR self gravity code
class AMRLevelSelfGravity : public AMRLevel
{
public:
  /// Constructor
  AMRLevelSelfGravity();

  /// Destructor
  virtual ~AMRLevelSelfGravity();

  /// This instance should never get called - historical
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMR level
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Define the parameters the object needs
  void defineParams(const Real&                 a_cfl,
                    const Real&                 a_domainLength,
                    const int&                  a_verbosity,
                    const int&                  a_tagBufferSize,
                    const int&                  a_maxInitRefLevel,
                    const Real&                 a_initialDtMultiplier,
                    const GodunovPhysics* const a_godunovPhysics,
                    const int&                  a_normalPredOrder,
                    const bool&                 a_useFourthOrderSlopes,
                    const bool&                 a_usePrimLimiting,
                    const bool&                 a_useCharLimiting,
                    const bool&                 a_useFlattening,
                    const bool&                 a_useArtificialViscosity,
                    const Real&                 a_artificialViscosity,
                    const RefCellTagger* const  a_refCellTagger,
                    const bool&                 a_useDeltaPhiCorr,
                    const StencilType&          a_stencil);

  /// Advance by one timestep
  virtual Real advance();

  /// Things to do after a timestep
  virtual void postTimeStep();

  /// Create tags for regridding
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// Set up data on this level after regridding
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// postRegridding ops
  virtual void postRegrid(int a_base_level);

  /// Initialize grids
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Initialize data
  virtual void initialData();

  /// Things to do after initialization
  virtual void postInitialize();

  /// Define Cosmology class for cosmological simulations
  //virtual void setCosmology(const Cosmology& a_cosmology);

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

#endif

  /// Returns the dt computed earlier for this level
  virtual Real computeDt();

  /// Compute dt using initial data
  virtual Real computeInitialDt();

  ///
  const LevelData<FArrayBox>& getStateNew() const;

  ///
  const LevelData<FArrayBox>& getStateOld() const;

  /// get the RHS to Poisson eq., interpolated to a_time
  virtual LevelData<FArrayBox>* getPoissonRhs(const Real& a_time);

  /// get the Phi, the gravitational potential
  virtual LevelData<FArrayBox>* getPhi(const Real& a_time);

  /// Solves Poisson eq. on the hierarchy of this and higher levels.
  virtual void gravity(const int  a_baseLevel, const bool a_srceCorr = true);

  // solves elliptic (Poisson's) equation
  virtual void ellipticSolver(const int  a_baseLevel,
                              const bool a_isLevelSolve);

  /// Overloaded version (see below) to be called by gravity()
  virtual void computeForce();

  /// compute 2nd order corrections for gravity source term
  virtual void secondOrderCorrection();

  /// set boolean stating whether or not this is the finest lev.
  virtual void isThisFinestLev(const bool a_isThisFinestLev);

  /// return bool stating whether or not this is the finest level.
  virtual bool isThisFinestLev() const;

  ///
  bool allDefined() const;

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelSelfGravity* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelSelfGravity* getFinerLevel() const;

  //
  void makePoissonRhs(LevelData<FArrayBox>&       a_rhs,
                      const LevelData<FArrayBox>& a_U,
                      const Real&                 a_time);

  // Compute force from a potential.
  void computeForce(LevelData<FArrayBox>&  a_force,
                    LevelData<FArrayBox>& a_phi,
                    const Real&           a_time);

  // compute source term
  void setupSourceTerm(LevelData<FArrayBox>&       a_source,
                       const LevelData<FArrayBox>& a_U,
                       const LevelData<FArrayBox>&  a_force,
                       const Real&                 a_time,
                       const Real&                 a_dt);

protected:
  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;

  // Gravitational potential Phi for projection at old and new time
  LevelData<FArrayBox> m_phiOld,m_phiNew,m_phiInt,m_deltaPhi;

  //
  LevelData<FArrayBox> m_forceOld, m_forceNew;

  // grad(Phi)=Rhs = 1.5*rho_gas
  LevelData<FArrayBox> m_rhs;

  //
  Real m_rhsOffset;

  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  // New and Old time step
  Real m_dtNew, m_dtOld;

  // Number of converved states
  int m_numStates;

  // Names of conserved states
  Vector<string> m_stateNames;

  // Names of plot variables
  Vector<string> m_plotNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Number of ghost cells (in each direction) for RHS to Poisson eq.
  int m_numRhsGhost;

  // Number of ghost for the potential; always 1 except for 4 pts stencil
  int m_numPhiGhost;

  // Number of ghost cells where it's necessary to compute the force
  int m_numForceGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // whether or not the deltaPhi correction should be applied
  bool m_useDeltaPhiCorr;

  // Interpolation from coarse to fine level
  FineInterp m_fineInterp;

  // Averaging from fine to coarse level
  CoarseAverage m_coarseAverage;

  // Interpolation from coarse to fine level
  FineInterp m_fineInterpPhi;

  // Averaging from fine to coarse level
  CoarseAverage m_coarseAveragePhi;

  // Interpolator for filling in ghost cells from the next coarser level
  PiecewiseLinearFillPatch m_patcher;

  // 2nd order interpolant for filling in ghost cells from next coarser level
  QuadCFInterp m_quadCFInterp;

  // Interpolator for filling in force ghost cells from the next coarser level
  PiecewiseLinearFillPatch     m_forcePatcher;

  // Level integrator for hydro
  LevelGodunov m_levelGodunov;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // gradient operator
  Gradient* m_gradient;

  // Patch integrator
  SelfGravityPhysics* m_gdnvPhysics;

  // Class that tags cells for refinement
  RefCellTagger* m_refCellTagger;

  // max level of refinement at start
  int m_maxInitRefLevel;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  // boolean stating whether or not this is the FINEst lev.
  bool m_isThisFinestLev;

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  int m_normalPredOrder;

  // Use 4th order slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  //
  DisjointBoxLayout m_grids;

  // True if all the parameters for this object are defined
  bool m_paramsDefined;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelSelfGravity& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelSelfGravity(const AMRLevelSelfGravity& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
