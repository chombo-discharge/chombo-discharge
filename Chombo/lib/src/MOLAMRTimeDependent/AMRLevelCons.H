#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELCONS_H_
#define _AMRLEVELCONS_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FourthOrderFineInterp.H"
#include "LevelFluxRegister.H"

#include "LevelConsOperator.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"

#include "NamespaceHeader.H"

/// AMR Level for fourth-order solution to conservation law
/**
 */
class AMRLevelCons : public AMRLevel
{
  friend class AMRLevelConsFactory ;
public:
  /// Constructor
  /**
   */
  AMRLevelCons();

  /// Destructor
  /**
   */
  virtual ~AMRLevelCons();

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep by RK4, and return new timestep.
  /**
   */
  virtual Real advance();

  /// add dissipation to m_Unew:  Kreiss-Oliger artificial viscosity
  /**
   */
  virtual void addDissipation();

  /// Things to do after a timestep -- reflux
  /**
     If (!m_hasFiner) then this does nothing except write out some
     diagnostics depending on the verbosity level.
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// Set up data on this level after regridding
  /**
   */
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  /**
   */
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Define m_Unew and m_Uold.
  /**
   */
  virtual void defineDataHolders();

  /// Initialize data:  Set m_Unew at time 0.
  /**
   */
  virtual void initialData();

  /// Set a_U at time a_time.
  /**
   */
  virtual void setData(LevelData<FArrayBox>& a_U,
                       Real a_time) const;

  /// Things to do after initialization
  /**
   */
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  /**
   */
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  /**
   */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  /**
   */
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();

  /// Return maximum speed from conserved variables
  /**
   */
  virtual Real getMaxWaveSpeed(const LevelData<FArrayBox>& a_U);

  /// Set the CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Set the spatial order of accuracy
  /**
     can be 2 or 4 (default)
  */
  virtual void spaceOrder(int a_spaceOrder);

  /// sets whether to limit face values in operator
  virtual void limitFaceValues(bool a_limitFaceValues);

  // sets whether to use high-order limiter
  virtual void highOrderLimiter(bool a_highOrderLimiter);

  /// sets whether initial data is average; default false
  virtual void initialAverage(bool a_initialAverage);

  /// sets whether to flatten extrapolations to faces
  virtual void useFlattening(bool a_useFlattening);

  /// sets whether to use artificial viscosity
  virtual void useArtVisc(bool a_useArtVisc);

  /// sets whether to turn off PPM
  virtual void noPPM(bool a_noPPM);

  /// sets m_doDeconvolution
  virtual void doDeconvolution(bool a_doDeconvolution);

  /// sets m_doFaceDeconvolution
  virtual void doFaceDeconvolution(bool a_doFaceDeconvolution);

  /// sets whether to use artificial viscosity (from divergence)
  virtual void useArtificialViscosity(bool a_useArtificialViscosity);

  /// sets coefficient of artificial viscosity (from divergence)
  virtual void artificialViscosity(Real a_artificialViscosity);

  /// sets coefficient of artificial viscosity
  virtual void ratioArtVisc(Real a_ratioArtVisc);

  /// sets whether forward Euler is used (instead of RK4)
  virtual void forwardEuler(bool a_forwardEuler);

  /// sets whether to enforce a min value
  virtual void enforceMinVal(bool a_enforceMinVal, Real a_minVal);

  /// Set the physical dimension of the longest side of the domain

  /**
   */
  virtual void domainLength(Real a_domainLength);

  /// Set the refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Set whether refinement threshold is scaled with dx
  /**
   */
  virtual void refinementIsScaled(bool a_refinementIsScaled);

  /// Set whether to tag on pressure instead of on density
  /**
   */
  virtual void tagPressure(bool a_tagPressure);

  /// Set whether to tag on vorticity instead of on density
  /**
   */
  virtual void tagVorticity(bool a_tagVorticity);

  /// Set the tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  /// Set the MOLPhysics class
  /**
   */
  void molPhysics(const MOLPhysics* const a_molPhysics);

  /// Accessors
  /**
   */
  LevelData<FArrayBox> * U()
  {
    return &m_Unew;
  }

  // return high-order estimate of vorticity of m_Unew
  void computeVorticity(LevelData<FArrayBox>& a_vorticity,
                        const LevelData<FArrayBox>& a_U) const;

protected:

  virtual void setDefaultValues();

  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  virtual void levelSetup();

  // Transfer settings from this to m_levelConsOperator:  called by levelSetup()
  virtual void transferSettingsToLevelOp();

  // Transfer settings from a_amrConsPtr to this:  called by define()
  virtual void transferSettingsFromAMRLevel(AMRLevelCons* a_amrConsPtr);

  // Get the next coarser level
  AMRLevelCons* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelCons* getFinerLevel() const;

  /// Advance by one timestep by RK4.
  virtual void advanceU();

  /// Compute new timestep from m_Unew, and store it in m_dtNew.
  virtual Real computeNewDt();

  /// Compute max vorticity of m_Unew, and store their magnitudes in a_vecMag.
  virtual Real maxVorticity(LevelData<FArrayBox>&  a_vecMag);

  /// Index within primitive variables for tagging cells
  virtual int indexForTagging();

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelCons& a_input);

  // Disallowed for all the usual reasons
  AMRLevelCons(const AMRLevelCons& a_input);

protected:

  // member variables

  // conserved variables U at old and new time
  LevelData<FArrayBox> m_Uold, m_Unew;

  // lower corner of problem domain in physical space
  Real m_r0;

  // CFL number
  Real m_cfl;

  // spatial order of accuracy (default is 4)
  int m_spaceOrder;

  // if true, limit face values in advection operator
  bool m_limitFaceValues;

  /// if true, use high-order limiter
  bool m_highOrderLimiter;

  // if true, initial values are average
  bool m_initialAverage;

  // if true, use flattening
  bool m_useFlattening;

  // if true, use artificial viscosity
  bool m_useArtVisc;

  // if m_useArtVisc, artificial viscosity coefficient
  Real m_ratioArtVisc;

  // use forward Euler if true, or RK4 if false
  bool m_forwardEuler;

  // if true, enforce min value on advected quantity
  bool m_enforceMinVal;

  // Avoid PPM?  Default false.
  bool m_noPPM;

  // Do deconvolution?  Default yes.
  bool m_doDeconvolution;

  // Do face deconvolution?  Default yes.
  bool m_doFaceDeconvolution;

  // use artificial viscosity (from divergence)?
  bool m_useArtificialViscosity;

  // coefficient of artificial viscosity (from divergence)
  Real m_artificialViscosity;

  // if enforcing minval, what value to enforce
  Real m_minVal;

  // Grid spacing
  Real m_dx;

  // Averaging from coarse to fine level
  CoarseAverage* m_coarseAveragePtr;

  // New time step
  Real m_dtNew;

  // Number of field variables (U) = number of conserved states
  int m_numStates;

  // Names of conserved states
  Vector<std::string> m_stateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Ghost cell vector
  IntVect m_ghostVect;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Interpolation from fine to coarse level
  FourthOrderFineInterp m_fineInterp;

  // Level integrator
  LevelConsOperator* m_levelConsOperatorPtr;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Pointer to the class defining the physics of the problem
  MOLPhysics* m_molPhysics;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Whether refinement threshold is scaled with dx
  bool m_refinementIsScaled;

  // Tag on pressure instead of on density?
  bool m_tagPressure;

  // Tag on vorticity instead of on density?
  bool m_tagVorticity;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  DisjointBoxLayout m_grids;
};

#include "NamespaceFooter.H"

#endif
