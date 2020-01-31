#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELVOLTEST_H_
#define _AMRLEVELVOLTEST_H_

#include "VolTestParameters.H"
#include "FArrayBox.H"
#include "BoxLayoutData.H"
#include "CoarseAverage.H"
#include "FourthOrderFineInterp.H"
#include "NewCoordSys.H"
#include "LevelGridMetrics.H"

#include "NamespaceHeader.H"

/*******************************************************************************
 */
///  Implementation of AMRLevel for testing volumes of mapped grids
/**
 *
 ******************************************************************************/

class AMRLevelVolTest : public AMRLevel
{
public:

  /// Default constructor
  AMRLevelVolTest();

  /// Factory constructor
  AMRLevelVolTest(MultiBlockCoordSysFactory *const a_coordSysFact,
                  const VT::ICType                 a_icType,
                  const VT::RefType                a_refType,
                  const int                        a_tagBufferSize,
                  const std::string&               a_plotPrefix);

  /// Destructor
  virtual
  ~AMRLevelVolTest();

  /// Defines this AMRLevel.
  virtual
  void define(AMRLevel*            a_coarser_level_ptr,
              const ProblemDomain& a_problem_domain,
              int                  a_level,
              int                  a_ref_ratio);

/*--------------------------------------------------------------------*
 * Member functions
 *--------------------------------------------------------------------*/

  /// Advances this level by one time step.
  virtual
  Real advance();

  ///
  virtual
  void postTimeStep();

  ///
  virtual
  void tagCells(IntVectSet& a_tags);

  ///
  virtual
  void tagCellsInit(IntVectSet& a_tags);

  ///
  virtual
  void preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids);

  ///
  virtual
  void regrid(const Vector<Box>& a_new_grids);

  /// ?
  virtual
  void postRegrid(int a_base_level);

  ///
  virtual
  void initialGrid(const Vector<Box>& a_new_grids);

  ///
  virtual
  void postInitialGrid(const bool a_restart);

  ///
  virtual
  void initialData();

  ///
  virtual
  void postInitialize();

  ///
  virtual
  Real computeDt();

  ///
  virtual
  Real computeInitialDt();

//--HDF5 Routines

#ifdef CH_USE_HDF5

  ///
  virtual
  void writeCheckpointHeader(HDF5Handle& a_handle) const;

  ///
  virtual
  void writeCheckpointLevel(HDF5Handle& a_handle) const;

  ///
  virtual
  void readCheckpointHeader(HDF5Handle& a_handle);

  ///
  virtual
  void readCheckpointLevel(HDF5Handle& a_handle);

  ///
  virtual
  void writePlotHeader(HDF5Handle& a_handle) const;

  ///
  virtual
  void writePlotLevel(HDF5Handle& a_handle) const;

  ///
  void writeMappedPlotFile() const;
#endif

  ///
  void conclude(int a_step) const;

protected:

  /// Create a load-balanced DisjointBoxLayout from a collection of
  /// boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids) const;

  /// Set up data structures for the level
  void levelSetup();

  /// Fill ghost cells and compute \<U\>
  void fillGhostsAndComputeU();

  /// Get the next coarser level
  AMRLevelVolTest* getCoarserLevel() const;

  /// Get the next finer level
  AMRLevelVolTest* getFinerLevel() const;

  /// Compute a norm (conservation or grid-convergence testing)
  void computeNorm(const int a_type, Real a_norm[3]) const;

  /// Check the result in each cell (free-stream preservation testing)
  int allCellsTest(Real a_testVal) const;

private:

  // No copy or assignment
  AMRLevelVolTest(const AMRLevelVolTest&);
  AMRLevelVolTest &operator=(const AMRLevelVolTest&);

/*--------------------------------------------------------------------*
 * Data members
 *--------------------------------------------------------------------*/

protected:

//--Parameters

  // Constructor initialized
  const int m_numStates;              ///< Number of conserved states
  int m_levelStep;                    ///< Iteration count for the level

  // Factory initialized
  MultiBlockCoordSysFactory* m_coordSysFact;
                                      ///< Factory for the coordinate system in
                                      ///< each block
  const VT::ICType m_icType;          ///< Type of IC
  const VT::RefType m_refType;        ///< Type of refinement
  int m_tagBufferSize;                ///< Tag buffer size
  std::string m_plotPrefix;           ///< Prefix for the plot file name

  // Initialized by define()
  Vector<string> m_stateNames;        ///< Names of conserved states
  int m_numGhost;                     ///< Number of ghost cells (for all
                                      ///< directions)
  int m_numGhostU;                    //** FIXME comment
  IntVect m_ghostVect;                ///< Number of ghost cells (in each
                                      ///< direction)
  IntVect m_ghostVectU;               //** FIXME comment

//--Utilities

  CoarseAverage m_coarseAverage;      ///< Averaging from fine to coarse level
  Copier m_exchangeCopier1;           ///< Exchanging 1 layer of valid ghosts
  Copier m_exchangeCopier;            ///< Exchanging for valid ghosts of \<U\>

//--Solution data

  LevelGridMetrics m_levelGridMetrics;
                                      ///< Provides grid metrics for this level
  LevelData<FArrayBox> m_U;           ///< Conserved state, U, in physical space
  LevelData<FArrayBox> m_JU;          ///< Conserved state, JU

  bool m_hasCoarser;                  ///< T - Coarser level exists
  bool m_hasFiner;                    ///< T - Finer level exists

public:

  friend class LevelGridMetrics;

};

#include "NamespaceFooter.H"

#endif
