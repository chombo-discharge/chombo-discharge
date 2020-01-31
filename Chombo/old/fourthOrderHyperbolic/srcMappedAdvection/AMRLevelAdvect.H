#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELADVECT_H_
#define _AMRLEVELADVECT_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"

#include "LevelAdvectOperator.H"
#include "PatchGodunov.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "BasicIBC.H"
#include "CoordSys.H"
#include "UsingNamespace.H"

/// AMR Level for fourth-order advection
/**
 */
class AMRLevelAdvect : public AMRLevel
{
  friend class AMRLevelAdvectFactory ;
public:
  /// Constructor
  /**
   */
  AMRLevelAdvect();

  /// Destructor
  /**
   */
  virtual ~AMRLevelAdvect();

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep -- RK4
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

  /// write mapped-grid info
  void writeMappedPlotFile() const;
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

  /// Set the spatial order of accuracy
  /**
     can be 2 or 4 (default)
  */
  virtual void spaceOrder(int a_spaceOrder);

  /// sets whether to limit face values in advection operator
  virtual void limitFaceValues(bool a_limitFaceValues);

  /// sets whether to enforce a min value
  virtual void enforceMinVal(bool a_enforceMinVal, Real a_minVal);

  /// sets the low-order flux for min val enforcement [0:DCU,1:CTU]
  virtual void lowOrderFluxScheme(std::string a_lowOrderFluxScheme);

  /// sets whether to use redistribution to preserve positivity
  virtual void redistributeNegativeVal(bool a_redistributeNegativeVal,
                                       int a_maxRedistributionPasses);

  /// sets whether to use hyperviscous stabilization (default is false)
  virtual void useHyperviscosity(bool a_useHyperviscosity);

  /// sets hyperviscosity coefficient
  virtual void hyperviscosity(Real a_hyperviscosity);

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

  void IBC(BasicIBC* a_ibc)
  {
    m_basicIBCPtr = a_ibc;
  }

  void coordinateSystem(CoordSysFactory<FArrayBox,FluxBox>* a_coordSysFact);

  /// Accessors
  /**
   */
  LevelData<FArrayBox> * phi()
  {
    return &m_phiNew;
  }

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelAdvect* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelAdvect* getFinerLevel() const;

  /// compute max- and min-preserving low-order fluxes
  /** note that this takes the face-centered N^T*vel (which is the
      area-weighted normal velocity on faces)
  */
  void computeLowOrderFlux(LevelData<FluxBox>& a_CTUflux,
                           const LevelData<FArrayBox>& a_phi,
                           const LevelData<FluxBox>& a_NTvel,
                           Real a_time,
                           Real a_dt);
#if 0
  /// apply Zalesak limiter to enforce minVal
  /** a_flux comes in with high-order flux and is modified in place using
      the low-order flux
  */
  void applyZalesakLimiter(LevelData<FluxBox>& a_flux,
                           const LevelData<FluxBox>& a_lowOrderFlux,
                           const LevelData<FArrayBox>& a_oldPhi,
                           Real minVal);
#endif

  /// write out error norms
  void writeErrorNorms() const;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelAdvect& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelAdvect(const AMRLevelAdvect& a_input)
  {
    MayDay::Error("invalid operator");
  }

protected:

    // member variables

  //
  int m_levelStep;

    // lower corner of problem domain in physical space
    Real m_r0 ;


    // phi at old and new time
    LevelData<FArrayBox> m_phiOld,m_phiNew;

    // cell-centered advection velocity
    LevelData<FArrayBox> m_advVel;

  // CFL number
  Real m_cfl;

  // spatial order of accuracy (default is 4)
  int m_spaceOrder;

  // if true, limit face values in advection operator
  bool m_limitFaceValues;

  // if true, enforce min value on advected quantity
  bool m_enforceMinVal;

  // if enforcing minval, what value to enforce
  Real m_minVal;

  // if enforcing minval, low order scheme to use
  int m_lowOrderFluxScheme;

  // if true, enforce postivity with redistribution
  bool m_redistributeNegativeVal;

  // if enforcing positivity with redistribution, number of passes over which to redistribute
  int m_maxRedistributionPasses;

  /// if true, proceed as if to do Zalesak limiting, but use high-order fluxes
  /**
      useful for debugging -- essentially the same as minVal = -infinity
  */
  bool m_highOrderFluxesOnly;

  /// if true, use low-order fluxes only
  /**
      useful for debugging low-order scheme
  */
  bool m_lowOrderFluxesOnly;

  // if true, add hyperviscous fluxes to each flux for stabilization
  bool m_useHyperviscosity;
  Real m_hyperviscosity;

  // Grid spacing
  Real m_dx;

  // Interpolation from fine to coarse level
  FineInterp m_fineInterp;

  // Averaging from coarse to fine level
  CoarseAverage m_coarseAverage;

  // New time step
  Real m_dtNew;

  // Number of field variables (phi) = number of conserved states
  int m_numStates;

  // Names of conserved states
  Vector<std::string> m_stateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Level integrator
  LevelAdvectOperator m_levelAdvectOperator;

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

  BasicIBC* m_basicIBCPtr;

  // PatchGodunov  object for doing CTU (needed for postivity preservation)
  PatchGodunov m_CTUpatchGod;

  // Coordinate system
  CoordSys<FArrayBox,FluxBox>* m_coordSysPtr;
  CoordSysFactory<FArrayBox,FluxBox>* m_coordSysFactPtr;

  /// need this in AMRLevel in order to generate mapped-grid output
  std::string m_plotfile_prefix;

};

#endif
