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
#include "GodunovPhysics.H"
#include "GodunovUtilities.H"
#include "PhysIBC.H"
#include "Tuple.H"
#include "UsingNamespace.H"

/// Patch Cons Operator - a class for evaluating the Laplacian.

/**
 */
class PatchMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  PatchMappedConsOperator();

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
  void define(const ProblemDomain&      a_domain,
              const Real&               a_dx,
              const GodunovPhysics* const  a_gdnvPhysics,
              const int&                a_numStates);

  /// Set the current time before calling evalRHS or addArtificialViscosity.
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Set the current box before calling evalRHS or addArtificialViscosity.
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Evaluate the operator (div(flux) ) at a given time.
  /**
     For the coarsest level, UcoarseOld and UcoarseNew are empty
     LevelData<FArrayBox> objects.  Also, either UcoarseOld or
     UcoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  U must be defined on the same
     DisjointBoxLayouts as were used in define().  Coarse flux register
     is flux register with the next coarser level.  Fine flux register is
     the flux register with the next finer level.  To the finer level FR,
     this level is the coarse level.  To the coarser level FR, this level
     is the fine level.
   */
  void evalRHS(FArrayBox&          a_LofU,
               const FArrayBox&    a_JUavgFab,
               const FArrayBox&    a_JFab,
               const FluxBox&      a_faceMetricTerms,
               const Tuple<IntVect, SpaceDim>& a_metricTermComponents,
               FluxBox&            a_FfaceAvg,
               Real                a_weight,
               const Box&          a_gridBox,
               bool                a_setFlattening,
               FArrayBox&          a_flatteningFab);

  /// update solution -- soln += dt*rhs (required by LevelRK4)
  void updateODE(FArrayBox&        a_solnFab,
                 const FArrayBox&  a_rhsFab,
                 Real a_dt);

  /// limits the face-centered W using Colella&Sekora limiter; also solves Riemann problem
  void getFaceAvg(// we'll apply limiter to a_faceW, face-averaged primitive variables
                  FluxBox& a_faceW,
                  /// cell-averaged primitive variables
                  const FArrayBox& a_cellW,
                  const FArrayBox& a_WofUavg,
                  const FluxBox& a_unitNormals,
                  // CELL-centered base box
                  const Box& a_gridBox,
                  FArrayBox& a_flatteningFab,
                  bool a_setFlattening);

  /// set spatial order of accuracy
  /*
    Can be 2 or 4 (default)
   */
  void spaceOrder(int a_spaceOrder);


  /// if true, limit face values
  void limitFaceValues(bool a_limitFaceValues);

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

  /// add artificial viscosity to a_Unew
  void addArtificialViscosity(FArrayBox&         a_UnewFab,
                              const FArrayBox&   a_UoldFab,
                              FluxBox&           a_flux,
                              const Box&         a_gridBox,
                              Real  a_weight);

  void computeCompFaceFluxes( FluxBox& a_uTimesV,
                              const FluxBox& a_u,
                              const FluxBox& a_v) const;

  void cellUJToCellU(
                     FArrayBox& a_UavgFab,
                     const FArrayBox& a_UJavgFab,
                     const FArrayBox& a_JFab) const;

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
   */
  virtual bool isDefined() const;

protected:

  void forwardBasisTransform(FArrayBox& a_W,
                             const FArrayBox& a_unitNormalFab) const;

  void reverseBasisTransform(FArrayBox& a_W,
                             const FArrayBox& a_unitNormalFab) const;

  // define() has been called
  bool m_isDefined;

  // Grid spacing
  Real m_dx;

  // Problem domain - index space for this level
  ProblemDomain m_domain;

  GodunovPhysics* m_gdnvPhysics;

  GodunovUtilities m_util;

  // spatial order of accuracy
  int m_spaceOrder;

  // if true, limit face values
  bool m_limitFaceValues;

  // if true, use flattening
  bool m_useFlattening;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  int m_numGhost;
  int m_numFields;

  // Flags for the existence of coarser and finer levels, respectively
  bool m_hasCoarser;
  bool m_hasFiner;

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

  // Current time and has it been set
  Real m_currentTime;
  bool m_isCurrentTimeSet;

  // Current box and has it been set
  Box  m_currentBox;
  bool m_isCurrentBoxSet;

private:
  // Disallowed for all the usual reasons
  //  void operator=(const PatchMappedConsOperator& a_input)
  //  {
  //    MayDay::Error("invalid operator");
  //  }

  // Disallowed for all the usual reasons
  //  PatchMappedConsOperator(const PatchMappedConsOperator& a_input)
  //  {
  //    MayDay::Error("invalid operator");
  //  }
};

#endif
