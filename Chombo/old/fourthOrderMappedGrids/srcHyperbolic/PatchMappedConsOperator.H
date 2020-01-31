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

#include "PatchConsOperator.H"
#include "MOLPhysics.H"
#include "Tuple.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Patch Cons Operator - a class for evaluating the Laplacian.

/**
 */
class PatchMappedConsOperator : public PatchConsOperator
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
  virtual void define(const ProblemDomain&      a_domain,
                      const Real&               a_dx,
                      const MOLPhysics* const  a_molPhysics,
                      const int&                a_numStates);

  /// Set the current box before calling evalRHS or addArtificialViscosity.
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Set coordinate system.
  /**
   */
  virtual void setCoordSys(MultiBlockCoordSys* a_coordSysPtr);


  /// Evaluate the normal flux.
  virtual void getNormalFlux(FluxBox&            a_FfaceAvg,
                             const FArrayBox&    a_UavgFab,
                             Real                a_weight,
                             bool                a_setFlattening,
                             FArrayBox&          a_flatteningFab);

  virtual void computeCompFaceFluxes( FluxBox& a_uTimesV,
                                      const FluxBox& a_u,
                                      const FluxBox& a_v) const;

  virtual void cellUJToCellU(FArrayBox& a_UavgFab,
                             const FArrayBox& a_UJavgFab) const;

  virtual void cellUToCellUJ(FArrayBox& a_UJavgFab,
                             const FArrayBox& a_UavgFab) const;

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

  virtual void setBasisVectors(const Box& a_box,
                               int a_dir);

  virtual void unsetBasisVectors();

  virtual void forwardBasisTransform(FArrayBox& a_W);

  virtual void reverseBasisTransform(FArrayBox& a_W);

  virtual void getAllFluxes(FluxBox&        a_FfaceAvg,
                            FluxBox&        a_FfaceCen,
                            const FluxBox&  a_WfaceAvg,
                            const FluxBox&  a_WfaceCen);

  virtual void getDirFluxes(FArrayBox&  a_Fface,
                            const FArrayBox&  a_Wface,
                            const Box&  a_box);

  virtual void reduceFlux(FluxBox&     a_FfaceAvg,
                          FluxBox&     a_FfaceForGrad);

  //  MOLPhysics* m_molPhysics;

  // [face direction][component]
  Tuple<IntVect, SpaceDim> m_metricTermComponents;

  MultiBlockCoordSys* m_coordSysPtr;

  // NewFourthOrderCoordSys* m_coordSysPtr;
  NewCoordSys* m_coordSysCurrentPtr;

  Vector<Box> m_mappingBlocks;

  ProblemDomain m_blockDomain;

  FluxBox* m_faceMetricTerms;

  FArrayBox* m_unitNormalFabPtr;
};

#include "NamespaceFooter.H"

#endif
