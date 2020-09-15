#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHMAPPEDCONSOPERATOR_H_
#define _PATCHMAPPEDCONSOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "FArrayBox.H"
#include "FluxBox.H"
#include "ProblemDomain.H"
#include "MultiBlockCoordSys.H"
#include "MOLMappedPhysics.H"
#include "MOLUtilities.H"
#include "LevelGridMetrics.H"
#include "PatchConsOperator.H"

#include "NamespaceHeader.H"

//--Forward declarations

template<class T> class LayoutData;

/// Patch Mapped Cons Operator - a class for evaluating hyperbolic conservation fluxes
/// on mapped multiblock grids

/**
 */
class PatchMappedConsOperator : public PatchConsOperator
{
public:
  /// Constructor
  /**
     Object requires define() to be called before all other functions.
  */
  PatchMappedConsOperator(const int a_numFluxesPerField = 1);

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  virtual ~PatchMappedConsOperator();

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).  a_refine is the refinement ratio between
     this level and the next coarser level.  For the coarsest level, an
     empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.
   */
  virtual void define(const ProblemDomain&      a_domain,
                      const Real&               a_dx,
                      LevelGridMetrics*         a_levelGridMetrics,
                      const MOLPhysics* const   a_molPhysics, // MOLMappedPhysics
                      const int&                a_numStates);

  /// Set the current time before calling evalRHS or addArtificialViscosity.
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Set the current data index and box before calling
  /** Required before calling evalRHS or addArtificialViscosity.
   */
  virtual void setCurrentBox(const DataIndex& a_currentDataIdx,
                             const Box&       a_currentBox);

  /// Evaluate the flux at a given time.
  virtual void evalFlux(const FArrayBox&    a_UavgFab,
                        FluxBox&            a_FfaceAvg,
                        FluxBox&            a_FfromWfaceAvg,
                        Real                a_weight,
                        bool                a_setFlattening,
                        FArrayBox&          a_flatteningFab);

  /// Evaluate the operator (div(flux)) at a given time.
  virtual void evalRHS(FArrayBox&          a_LofU,
                       const FluxBox&      a_FfaceAvg);

  /// finds face-centered W and limits it using Colella&Sekora limiter; also solves Riemann problem
  virtual void getFaceAvg(// we'll apply limiter to a_faceW, face-averaged primitive variables
                          FluxBox& a_faceW,
                          /// cell-averaged primitive variables
                          const FArrayBox& a_cellW,
                          const FArrayBox& a_WofUavg,
                          FArrayBox& a_flatteningFab,
                          bool a_setFlattening);

  /// sets m_unitNormalLayPtr
  void unitNormals(const LayoutData<FluxBox> *const a_unitNormalLayPtr);

  /// sets m_maxSpeed
  void setMaxSpeed(Real a_maxSpeed);

  /// sets whether to use artificial dissipation
  void useArtificialDissipation(bool a_useArtificialDissipation);

  /// sets coefficient of artificial dissipation
  void artificialDissipation(Real a_artificialDissipation);

  /// add artificial viscosity to a_Unew
  virtual void addArtificialViscosity(FArrayBox&         a_UnewFab,
                                      const FArrayBox&   a_UoldFab,
                                      FluxBox&           a_flux,
                                      Real  a_weight);

  /// Add artificial viscosity to a_JUnew for mapped grids
  void addMappedArtificialViscosity(FArrayBox&       a_JUnewFab,
                                    const FArrayBox& a_UoldFab,
                                    FluxBox&         a_flux,
                                    const FluxBox&   a_N,
                                    const FArrayBox& a_J,
                                    const Box&       a_box,
                                    Real             a_weight);

  virtual Real getMaxWaveSpeed(const FArrayBox&   a_U,
                               const Box&         a_bx);

protected:

  /// may be used in derived classes just before calling Riemann solver in direction a_dir
  virtual void preRiemann(FArrayBox&  a_WLeft,
                          FArrayBox&  a_WRight,
                          int         a_dir,
                          const Box&  a_box);

  /// may be used in derived classes just after calling Riemann solver in direction a_dir
  virtual void postRiemann(FArrayBox&  a_Wface,
                           int         a_dir,
                           const Box&  a_box);

  virtual void forwardBasisTransform(FArrayBox& a_W,
                                     const int a_dir,
                                     const Box& a_box);

  virtual void reverseBasisTransform(FArrayBox& a_W,
                                     const int a_dir,
                                     const Box& a_box);

  virtual void getAllFluxes(FluxBox&        a_FfaceAvg,
                            FluxBox&        a_FfaceCen,
                            const FluxBox&  a_WfaceAvg,
                            const FluxBox&  a_WfaceCen);

  virtual void fluxCorrection(FluxBox&         a_fluxAvg,
                              const FArrayBox& a_UavgFab);

  // Domain on a block -- maintains ghosts cells for interpolation between blocks,
  // but none at the domain boundary.
  ProblemDomain m_blockDomain;

  // Multi-block coordinate system.
  // NOTE: We don't "own" this, so we don't deallocate it upon destruction.
  const MultiBlockCoordSys* m_coordSys;

  MOLMappedPhysics* m_molMappedPhysics;
  // MOLPhysics* m_molPhysics;

  // m_levelGridMetrics is not const, because blockDomain modifies it.
  LevelGridMetrics* m_levelGridMetrics;

  // Unit normals (for the level) if using mapped grids
  const LayoutData<FluxBox>* m_unitNormalLayPtr;

  // Whether to add artificial dissipation to fluxes
  bool m_useArtificialDissipation;

  // coefficient of artificial dissipation
  Real m_artificialDissipation;

  // maximum speed, used with artificial dissipation
  Real m_maxSpeed;

  // whether m_maxSpeed has been set
  bool m_isMaxSpeedSet;

  // Current data index
  DataIndex m_currentDataIdx;

private:
  // Disallowed for all the usual reasons
  PatchMappedConsOperator& operator=(const PatchMappedConsOperator& a_input);
  PatchMappedConsOperator(const PatchMappedConsOperator& a_input);
};

#include "NamespaceFooter.H"

#endif
