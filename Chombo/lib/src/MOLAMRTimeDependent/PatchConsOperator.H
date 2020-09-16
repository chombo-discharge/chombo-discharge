#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHCONSOPERATOR_H_
#define _PATCHCONSOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "FArrayBox.H"
#include "FluxBox.H"
#include "ProblemDomain.H"
#include "MOLPhysics.H"
#include "MOLUtilities.H"

#include "NamespaceHeader.H"

/// Patch Cons Operator - a class for evaluating the Laplacian.

/**
 */
class PatchConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  PatchConsOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  virtual ~PatchConsOperator();

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
                      const MOLPhysics* const  a_molPhysics,
                      const int&                a_numStates);

  /// Set the current time before calling evalRHS or addArtificialViscosity.
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Set the current box before calling evalRHS or addArtificialViscosity.
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Evaluate the normal flux.
  virtual void getNormalFlux(FluxBox&            a_FfaceAvg,
                             const FArrayBox&    a_UavgFab,
                             Real                a_weight,
                             bool                a_setFlattening,
                             FArrayBox&          a_flatteningFab);

  /// update solution -- soln += dt*rhs (required by LevelRK4)
  virtual void updateODE(FArrayBox&        a_solnFab,
                         const FArrayBox&  a_rhsFab,
                         Real a_dt);

  /// finds face-centered W and limits it using Colella&Sekora limiter; also solves Riemann problem
  virtual void getFaceAvg(// we'll apply limiter to a_faceW, face-averaged primitive variables
                          FluxBox& a_faceW,
                          /// cell-averaged primitive variables
                          const FArrayBox& a_cellW,
                          const FArrayBox& a_WofUavg,
                          FArrayBox& a_flatteningFab,
                          bool a_setFlattening);

  /// set spatial order of accuracy
  /*
    Can be 2 or 4 (default)
   */
  void spaceOrder(int a_spaceOrder);


  /// if true, limit face values
  void limitFaceValues(bool a_limitFaceValues);

  // sets whether to use high-order limiter
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

  /// sets m_numGhost, the number of ghost cells
  void numGhost(int a_numGhost);

  /// add artificial viscosity to a_Unew
  virtual void addArtificialViscosity(FArrayBox&         a_UnewFab,
                                      const FArrayBox&   a_UoldFab,
                                      FluxBox&           a_flux,
                                      Real  a_weight);

  virtual void getFluxDivergence(FArrayBox&   a_LofU,
                                 const FluxBox&     a_FfaceAvg);

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
   */
  virtual bool isDefined() const;

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

  virtual void getAllFluxes(FluxBox&        a_FfaceAvg,
                            FluxBox&        a_FfaceCen,
                            const FluxBox&  a_WfaceAvg,
                            const FluxBox&  a_WfaceCen);

  // may be used by derived classes
  virtual void reduceFlux(FluxBox&     a_FfaceAvg,
                          FluxBox&     a_FfaceForGrad);

  // define() has been called
  bool m_isDefined;

  // Grid spacing
  Real m_dx;

  // Problem domain - index space for this level
  ProblemDomain m_domain;

  MOLPhysics* m_molPhysics;

  MOLUtilities m_util;

  // spatial order of accuracy
  int m_spaceOrder;

  // if true, limit face values
  bool m_limitFaceValues;

  // if true, use high-order limiter
  bool m_highOrderLimiter;

  // if true, use flattening
  bool m_useFlattening;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  int m_numGhost;
  int m_numFields;

  int m_numFluxesPerField;

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
  void operator=(const PatchConsOperator& a_input);

  // Disallowed for all the usual reasons
  PatchConsOperator(const PatchConsOperator& a_input);
};

#include "NamespaceFooter.H"

#endif
