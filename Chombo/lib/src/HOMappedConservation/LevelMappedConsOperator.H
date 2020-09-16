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
#include "PatchMappedConsOperator.H"
#include "TimeInterpolatorRK4.H"
#include "LevelSourceTerm.H"
#include "LevelConsOperator.H"

#include "NamespaceHeader.H"

//--Forward declarations

class LevelGridMetrics;
class NewFourthOrderCoordSys;
class BlockRegister;

//! \class LevelMappedConsOperator
//! This class defines a flux divergence operator for hyperbolic conservation
//! laws on a single refinement level of a set of structured grids.
class LevelMappedConsOperator : public LevelConsOperator
{
public:

  /// Default constructor
  LevelMappedConsOperator();

  /// Constructs a Level operator with the given PatchMappedConsOperator instance.
  explicit LevelMappedConsOperator(PatchMappedConsOperator* a_operator);

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



  /// sets whether to use source term
  void useSourceTerm(bool a_useSourceTerm);

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

  /// sets whether to use artificial dissipation
  void useArtificialDissipation(bool a_useArtificialDissipation);
  
  /// sets coefficient of artificial dissipation
  void artificialDissipation(Real a_artificialDissipation);
  
protected:

  /// sets m_flattening
  virtual void defineFlattening();

  /// sets m_patchMappedConsOperator
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

  /// return maximum speed from m_UPtr
  virtual Real getMaxSpeed();

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
  PatchMappedConsOperator* m_patchMappedConsOperatorPtr;

  // The time interpolator
  TimeInterpolatorRK4 m_timeInterpolator;

  // Indexing of metrics components [face direction][component]
  Tuple<IntVect, SpaceDim> m_metricTermComponents;


  // A velocity interval, if it exists in the conserved variables
  Interval m_velocityIntv;

  // Interval within fluxes of vector components
  Interval m_vectorFluxIntv;

  // m_levelGridMetricsPtr->getBoxes() grown by 1 then intersected with domain
  BoxLayout m_grow1inDomainLayout;

  // use source term?  Default false.  If true, in m_sourceTermPtr.
  bool m_useSourceTerm;

  // Source term
  LevelSourceTerm* m_sourceTermPtr;

  // use artificial dissipation?  Default false.
  bool m_useArtificialDissipation;

  // coefficient of artificial dissipation
  Real m_artificialDissipation;

//--Multiblock only

  Tuple< LayoutData<IntVectSet>, 2*SpaceDim> m_facesToFill;
};

#include "NamespaceFooter.H"

#endif
