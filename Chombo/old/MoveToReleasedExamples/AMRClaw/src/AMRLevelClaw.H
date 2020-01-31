#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELCLAW_H_
#define _AMRLEVELCLAW_H_

#include <string>
using std::string;

// Base class
#include "AMRLevel.H"

// And members
#include "CoarseAverageClaw.H"
#include "FineInterpClaw.H"
#include "PiecewiseLinearFillPatch.H"
#include "LevelFluxRegister.H"
#include "ClawPatch.H"

#include "NamespaceHeader.H"

/// AMR Level for second-order wave equation.
/**
 */
class AMRLevelClaw : public AMRLevel
{
    friend class AMRLevelClawFactory ;
public:
    /// Constructor
    /**
     */
    AMRLevelClaw();

    /// Destructor
    /**
     */
    virtual ~AMRLevelClaw();

    /// Define new AMR level
    /**
     */
    virtual void define(AMRLevel*            a_coarserLevelPtr,
                        const ProblemDomain& a_problemDomain,
                        int                  a_level,
                        int                  a_refRatio);

    /**
     */

    virtual Real advance();

    /// Things to do after a timestep -- reflux
    /**
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

    /// Initialize data
    /**
     */
    virtual void initialData();

    /// Things to do after initialization
    /**
     */
    virtual void postInitialize();

    virtual void postRegrid(int a_base_level);

    virtual void clawPatch(const ClawPatch& a_clawPatch);

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

  /// Set the CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Set the physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  /// Set the refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Set the tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelClaw* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelClaw* getFinerLevel() const;


private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelClaw& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelClaw(const AMRLevelClaw& a_input)
  {
    MayDay::Error("invalid operator");
  }

protected:

    // member variables

    // lower corner of problem domain in physical space

    // Solution phi at old and new time
    LevelData<FArrayBox> m_phiOld,m_phiNew;
    LevelData<FArrayBox> m_aux;
    LevelData<FArrayBox> m_kappa;
    LevelData<FArrayBox> m_coarsened_kappa;
    int m_maux;  // Number of actual aux array variables.
    int m_mcapa;
    bool m_hasKappa;

    // CFL number
    Real m_cfl;

  // Grid spacing
  Real m_dx;

  // Interpolation from fine to coarse level
  FineInterpClaw m_fineInterp;

  // Averaging from coarse to fine level
  CoarseAverageClaw m_coarseAverage;

  // Boundary filler.
  PiecewiseLinearFillPatch m_patcher;

  // New time step
  Real m_dtNew;

  Real m_coarse_cfl;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  DisjointBoxLayout m_grids;

  // Initial / boundary condition interface class.

  // Patch integrator (factory and object)
    ClawPatch m_clawPatch;

    // Names of conserved states
    std::vector<string> m_stateNames;

    // Physical dimension of the longest side of the domain
    Real m_domainLength;

    int m_numGhost;  // set to m_mbc,  from clawpack input, above. See setClawParams.
    int m_numStates; // set to m_meqn, from clawpack input, above. See setClawParams.
};

#include "NamespaceFooter.H"

#endif
