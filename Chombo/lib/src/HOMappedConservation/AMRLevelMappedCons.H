#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMAPPEDCONS_H_
#define _AMRLEVELMAPPEDCONS_H_

#include "Vector.H"
#include "IntVectSet.H"
#include "Box.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "LevelFluxRegister.H"
#include "CoarseAverage.H"
#include "AMRLevelCons.H"
#include "MappedLevelData.H"
#include "MultiBlockCoordSys.H"
#include "MultiBlockUtil.H"
#include "AMRLevelMappedStabilityStrategy.H"
#include "AMRLevelMappedTaggingStrategy.H"

#include "LevelMappedConsOperator.H"

#include "NamespaceHeader.H"

//--Forward declarations

class LevelGridMetrics;


/*******************************************************************************
 */
/// AMR Level for fourth-order solution to conservation law
/**
 *//*+*************************************************************************/

class AMRLevelMappedCons : public AMRLevelCons
{
  friend class AMRLevelMappedConsFactory;
  friend class LevelGridMetrics;

public:

  /// Factory constructor
  AMRLevelMappedCons(MultiBlockCoordSysFactory *const a_coordSysFact,
                     const AMRLevelMappedStabilityStrategy* const a_stabilityStrategy,
                     const AMRLevelMappedTaggingStrategy* const a_taggingStategy);

  /// Destructor
  virtual ~AMRLevelMappedCons();


/*==============================================================================
 * Member functions
 *============================================================================*/

  /// Define new AMR level
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep by RK4, and return new timestep.
  virtual Real advance();

  /// add dissipation to m_Unew:  Kreiss-Oliger artificial viscosity
//   virtual void addDissipation();

  /// Things to do after a timestep -- reflux
  virtual void postTimeStep();

  /// Create tags for regridding
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// Pre-regrid operations necessary for updating metric terms
  virtual void preRegrid(int a_base_level,
                         const Vector<Vector<Box> >& a_new_grids);

  /// Set up data on this level after regridding
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Compute the grid metrics for mapped grids
  virtual void postInitialGrid(const bool a_restart);

  /// Initialize data
  virtual void initialData();

  /// Things to do after initialization
  virtual void postInitialize();

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

  /// Write mapped-grid info
  virtual void writeMappedPlotFile() const;


#endif

  /// Conclude by testing conservation if necessary
  virtual void conclude(int a_step) const;

  /// Compute dt using initial data
  virtual Real computeInitialDt();

  void coordinateSystem(MultiBlockCoordSysFactory* a_coordSysFact);

  /// Return maximum speed from conserved variables
  virtual Real getMaxWaveSpeed(const LevelData<FArrayBox>& a_U);

  /// Set whether to use a source term
  void useSourceTerm(bool a_useSourceTerm);

  /// Set the source term
  void sourceTerm(const LevelSourceTerm* const a_sourceTerm);

  /// Set whether to add artificial dissipation
  void useArtificialDissipation(bool a_useArtificialDissipation);

  /// Set coefficient of artificial dissipation
  void artificialDissipation(Real a_artificialDissipation);

  /// Set the MOLPhysics class
/*
  void molPhysics(const MOLPhysics* const a_molPhysics);
*/

#if 0
  /// return high-order estimate of vorticity of m_Unew
  void computeVorticity(LevelData<FArrayBox>& a_vorticity,
                        const LevelData<FArrayBox>& a_U) const;
#endif

protected:

  virtual void setDefaultValues();

/*
  /// Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);
*/

  /// Setup menagerie of data structures
  virtual void levelSetup();

  /// Transfer settings from this to m_levelMappedConsOperator: called by levelSetup()
  virtual void transferSettingsToLevelOp();

  /// Transfer settings from a_amrConsPtr to this:  called by define()
  virtual void transferSettingsFromAMRLevelMapped(AMRLevelMappedCons* a_amrConsPtr);

  /// Get the next coarser level
  AMRLevelMappedCons* getCoarserLevel() const;

  /// Get the next finer level
  AMRLevelMappedCons* getFinerLevel() const;

  /// Advance by one timestep by RK4.
  virtual void advanceU();

  /// Compute new timestep from m_Unew, and store it in m_dtNew.
  virtual Real computeNewDt();

  /// Index within primitive variables for tagging cells
  virtual int indexForTagging();

  /// Compute a norm for some reason
  void myComputeNorm(const int a_type, Vector<Real>& a_norm) const;

  /// Compute and write out norm of conserved quantities
  void reportNorm() const;

  /// Report norms of errors (requires exact solution)
  void reportError() const;

  /// Initialize <U> and <JU> at time m_time.
  virtual void setData(LevelData<FArrayBox>& a_U,
                       LevelData<FArrayBox>& a_JU) const;

  void writeJ(bool a_writeJ)
  {
    s_writeJ = a_writeJ;
  }

  void writeExact(bool a_writeExact)
  {
    s_writeExact = a_writeExact;
  }

  void writeError(bool a_writeError)
  {
    s_writeError = a_writeError;
  }

  void writeMap(bool a_writeMap)
  {
    s_writeMap = a_writeMap;
  }

private:

  // Copy and assignment not allowed
  AMRLevelMappedCons(const AMRLevelMappedCons& a_input);
  void operator=(const AMRLevelMappedCons& a_input);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

//--Parameters initialized by factory constructor

  // Factory for the coordinate system in each block
  MultiBlockCoordSysFactory* m_coordSysFact;
  // A pointer to the instance of the coordinate system for this level.
  MultiBlockCoordSys* m_coordSys;
  // A pointer to the instance of the coordinate system for finer level.
  MultiBlockCoordSys* m_coordSysFine;

  MultiBlockUtil m_mbUtil;

//--Parameters initialized by define

  // Grid metrics for the level
  LevelGridMetrics *m_levelGridMetricsPtr;

  // Level data (conserved variables <JU> and <U>)
  MappedLevelData m_data;

  // Level integrator
  LevelMappedConsOperator* m_levelMappedConsOperatorPtr;

  // Multiblock-capable flux register
  LevelFluxRegister* m_fluxRegisterPtr;

  int m_levelStep;

  bool m_useSourceTerm;

  bool m_useArtificialDissipation;

  Real m_artificialDissipation;

  // Strategy for computing a stable time step.
  AMRLevelMappedStabilityStrategy* m_stabilityStrategy;

  // Strategy for determining which cells to tag for refinement.
  AMRLevelMappedTaggingStrategy* m_taggingStrategy;

  // Prefix for the plot file name
  std::string m_plotPrefix;

  LevelSourceTerm* m_sourceTermPtr;

  static int s_finestLevel;           ///< The finest level from the initial
                                      ///< grid or latest regrid.

//--Testing/validation

  static Vector<Real> s_L1NormJUConsvRef;
                                      ///< L_1 norm of JU at start of solution
                                      ///< used for checking conservation

  /// whether to write out J
  static bool s_writeJ;

  /// whether to write out exact value from time
  static bool s_writeExact;

  /// whether to write out error from time
  static bool s_writeError;

  /// whether to write out mapping file
  static bool s_writeMap;

  mutable RefCountedPtr< LevelData<FArrayBox> > m_UexactPtr;

  mutable RefCountedPtr< LevelData<FArrayBox> > m_JUexactPtr;
};

#include "NamespaceFooter.H"

#endif
