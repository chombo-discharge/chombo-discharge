#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELMAPPEDCONSOPERATOR_H_
#define _LEVELMAPPEDCONSOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "Tuple.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "LevelFluxRegister.H"
#include "MOLPhysics.H"
#include "PatchConsOperator.H"
#include "TimeInterpolatorRK4.H"
#include "LevelSourceTerm.H"

#include "NamespaceHeader.H"

//--Forward declarations

class LevelGridMetrics;
class NewFourthOrderCoordSys;
class BlockRegister;

//! \class LevelMappedConsOperator
//! This class defines a flux divergence operator for hyperbolic conservation
//! laws on a single refinement level of a set of structured grids.
class LevelMappedConsOperator
{
public:

  /// Default constructor
  LevelMappedConsOperator();

  /// Constructs a Level operator with the given PatchConsOperator instance.
  explicit LevelMappedConsOperator(PatchConsOperator* a_operator);

  /// Destructor
  virtual ~LevelMappedConsOperator();


/*==============================================================================
 * Member functions
 *============================================================================*/

  /// Define the object so that time stepping can begin (actual constructor)
  virtual void define(LevelGridMetrics *const     a_levelGridMetricsPtr,
                      LevelData<FArrayBox> *const a_UPtr,
                      const Copier *const         a_UExchangeCopierPtr,
                      const Copier *const         a_JUExchangeCopierPtr,
                      const DisjointBoxLayout*    a_coarserGrids,
                      const ProblemDomain&        a_domain,
                      const int                   a_numGhost,
                      const int                   a_numInterpolatorCrFnGhost,
                      const int                   a_refineCoarse,
                      const Real                  a_dx,
                      const MOLPhysics* const     a_molPhysics,
                      const int                   a_numStates,
                      const bool                  a_hasCoarser);

  /// Define indices in to the metrics matrix
  void defineMetricsIndices(const NewFourthOrderCoordSys *const a_coordSysPtr);

  /// Find unit normals for applying the Riemann problem on mapped grids
  void defineUnitNormals(LevelData<FluxBox>& a_NLev);

  /// Evaluate the operator \f$(-div(F))\f$ at a given time.
  virtual void evalRHS(LevelData<FArrayBox>&       a_LofU,
                       LevelData<FArrayBox>&       a_JU,
                       LevelFluxRegister&          a_finerFluxRegister,
                       LevelFluxRegister&          a_coarserFluxRegister,
                       const LevelData<FArrayBox>& a_JUcoarseOld,
                       const Real&                 a_timeCoarseOld,
                       const LevelData<FArrayBox>& a_JUcoarseNew,
                       const Real&                 a_timeCoarseNew,
                       Real                        a_time,
                       Real                        a_fluxweight);


  /// update solution -- soln += dt*rhs
  virtual void updateODE(LevelData<FArrayBox>&       a_soln,
                         const LevelData<FArrayBox>& a_rhs,
                         Real                        a_dt);

  /// define newSoln to match existingSoln, including ghost cells
  virtual void defineSolnData(LevelData<FArrayBox>&       a_newSoln,
                              const LevelData<FArrayBox>& a_existingSoln);


  /// define RHS data based on existingSoln (in this case, w/o ghost cells)
  virtual void defineRHSData(LevelData<FArrayBox>&       a_newRHS,
                             const LevelData<FArrayBox>& a_existingSoln);

  /// copy data from src->dest
  void copySolnData(LevelData<FArrayBox>&       a_dest,
                    const LevelData<FArrayBox>& a_src);

  /// set spatial order of accuracy
  void spaceOrder(int a_spaceOrder);

  /// if true, limit face values
  void limitFaceValues(bool a_limitFaceValues);

  /// sets whether to use high-order limiter
  void highOrderLimiter(bool a_highOrderLimiter);

  /// sets whether to flatten extrapolations to faces
  void useFlattening(bool a_useFlattening);

  /// sets m_noPPM
  void noPPM(bool a_noPPM);

  /// sets m_doDeconvolution
  void doDeconvolution(bool a_doDeconvolution);

  /// sets m_doFaceDeconvolution
  void doFaceDeconvolution(bool a_doFaceDeconvolution);

  /// sets whether to use artificial viscosity (from divergence)
  void useArtificialViscosity(bool a_useArtificialViscosity);

  /// sets coefficient of artificial viscosity (from divergence)
  void artificialViscosity(Real a_artificialViscosity);

  /// sets whether to use source term
  void useSourceTerm(bool a_useSourceTerm);

  /// sets m_numGhost, the number of ghost cells
  void numGhost(int a_numGhost);

  /// set number of calls that will be made to evalRHS
  void evalCountMax(int a_evalCountMax);

  /// reset m_evalCount to 0
  void resetEvalCount();

  /// add artificial viscosity to a_Unew
  virtual void addArtificialViscosity(
    LevelData<FArrayBox>&       a_Unew,
    const LevelData<FArrayBox>& a_Uold,
    LevelFluxRegister&          a_finerFluxRegister,
    LevelFluxRegister&          a_coarserFluxRegister,
    Real                        a_oldTime,
    Real                        a_weight);

  /// update a_finerFluxRegister and a_coarserFluxRegister
  virtual void updateFluxTotalsAndRegisters(
    FluxBox&           a_Fface,
    LevelFluxRegister& a_finerFluxRegister,
    LevelFluxRegister& a_coarserFluxRegister,
    const DataIndex&   a_dataIndex,
    Real               a_weight);

  /// Return time interpolator object within patch-fill object.
  TimeInterpolatorRK4& getTimeInterpolator();

  /// fill in ghost cells of a_U to RK4 intermediates
  virtual void fillGhostsRK4AndComputeU(
    LevelData<FArrayBox>&   a_JU,
    const Real&             a_time,
    int                     a_stage,
    const Real&             a_timeCoarseOld,
    const Real&             a_timeCoarseNew);

  /// set source term object
  void setSourceTerm(LevelSourceTerm* a_sourceTermPtr);

protected:

  /// sets m_flattening
  void defineFlattening();

  /// sets m_patchConsOperator
  virtual void definePatch(const MOLPhysics* const   a_molPhysics);

  virtual void setPatchIndex(const DataIndex&  a_ind) const;

//--Multiblock only (defined in LevelMappedConsOperator_MB.cpp)

  /// Define members required for multiblock grids
  void defineMultiblockMbrs();

  /// Sets single-valued fluxes at block boundaries
  // Using the contents of a BlockRegister, set the single-valued
  // fluxes at block boundaries.
  void setCommonFlux(LevelData<FluxBox>&  a_flux,
                     const BlockRegister& a_blockRegister) const;

//--End multiblock only

private:

  // Disallowed for all the usual reasons
  void operator=(const LevelMappedConsOperator&);
  LevelMappedConsOperator(const LevelMappedConsOperator&);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  // Metrics for the level
  LevelGridMetrics* m_levelGridMetricsPtr;

  // <U>
  LevelData<FArrayBox>* m_UPtr;

  // Unit normals for applying Riemann solvers on mapped grids
  LayoutData<FluxBox> m_unitNormalLay;

  // Exchange copier for <U>
  const Copier* m_UExchangeCopierPtr;

  // Exchange copier for <JU>
  const Copier* m_JUExchangeCopierPtr;

  // Patch integrator
  PatchConsOperator* m_patchConsOperatorPtr;

  // Number of ghost cells need locally for this level
  int m_numGhost;

  // The time interpolator
  TimeInterpolatorRK4 m_timeInterpolator;

  // Indexing of metrics components [face direction][component]
  Tuple<IntVect, SpaceDim> m_metricTermComponents;

  // Grid spacing
  Real m_dx;

  // Problem domain - index space for this level
  // NOTE: This is mutable so that it can be changed within setPatchIndex(), which
  // NOTE: is a const method.
  mutable ProblemDomain m_domain;

  // Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  // Number of conserved variables
  int m_numCons;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  // A velocity interval, if it exists in the conserved variables
  Interval m_velocityIntv;

  // Flags for the existence of coarser and finer levels, respectively
  bool m_hasCoarser;

  // m_levelGridMetricsPtr->getBoxes() grown by 1 then intersected with domain
  BoxLayout m_grow1inDomainLayout;

  // flattening coefficients, if m_useFlattening
  BoxLayoutData<FArrayBox> m_flattening;

  // spatial order of accuracy
  int m_spaceOrder;

  // if true, limit face values
  bool m_limitFaceValues;

  /// if true, use high-order limiter
  bool m_highOrderLimiter;

  // if true, use flattening
  bool m_useFlattening;

  int m_numFields;

  // Avoid PPM?  Default false.
  bool m_noPPM;

  // Do deconvolution?  Default yes.
  bool m_doDeconvolution;

  // Do face deconvolution?  Default yes.
  bool m_doFaceDeconvolution;

  // number of calls that will be made to evalRHS
  int m_evalCountMax;

  // number of calls that have been made to evalRHS so far
  int m_evalCount;

  // use artificial viscosity (from divergence)?
  bool m_useArtificialViscosity;

  // coefficient of artificial viscosity (from divergence)
  Real m_artificialViscosity;

  // use source term?  Default false.  If true, in m_sourceTermPtr.
  bool m_useSourceTerm;

  // define() has been called
  bool m_defined;

  // Source term
  LevelSourceTerm* m_sourceTermPtr;

//--Multiblock only

  Tuple< LayoutData<IntVectSet>, 2*SpaceDim> m_facesToFill;
};

#include "NamespaceFooter.H"

#endif
