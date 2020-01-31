#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHSHALLOWWATERMAPPEDOPERATOR_H_
#define _PATCHSHALLOWWATERMAPPEDOPERATOR_H_

#include "PatchMappedConsOperator.H"

#include "NamespaceHeader.H"

/// Patch Cons Operator - a class for evaluating the Laplacian.

/**
 */
class PatchShallowWaterMappedOperator : public PatchMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  PatchShallowWaterMappedOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  virtual ~PatchShallowWaterMappedOperator();

  virtual void preRiemann(FArrayBox&  a_WLeft,
                          FArrayBox&  a_WRight,
                          int         a_dir,
                          const Box&  a_box);

  virtual void postRiemann(FArrayBox&  a_Wface,
                           int         a_dir,
                           const Box&  a_box);

  virtual void setCurrentBox(const Box& a_currentBox);

  virtual void getNormalFlux(
                             FluxBox&            a_FfaceAvg,
                             FluxBox&            a_WfaceOrthoCen,
                             const FArrayBox&    a_UavgFab,
                             Real                a_weight,
                             bool                a_setFlattening,
                             FArrayBox&          a_flatteningFab);

  /// finds face-centered W and limits it using Colella&Sekora limiter; also solves Riemann problem
  virtual void getFaceAvg(// we'll apply limiter to a_faceW, face-averaged primitive variables
                          FluxBox& a_faceW,
                          /// cell-averaged primitive variables
                          const FArrayBox& a_cellW,
                          const FArrayBox& a_WofUavg,
                          FArrayBox& a_flatteningFab,
                          bool a_setFlattening);

  void setContravariantMetricFace(const FluxBox* a_contravariantMetricFacePtr);

  void setOrthoMatrix(const FluxBox* a_orthoMatrixFacePtr);

protected:

  virtual void getAllFluxes(FluxBox&        a_FfaceAvg,
                            FluxBox&        a_FfaceCen,
                            const FluxBox&  a_WfaceAvg,
                            const FluxBox&  a_WfaceCen);

  FluxBox* m_contravariantMetricFacePtr;

  FluxBox* m_orthoMatrixFacePtr;
};

#include "NamespaceFooter.H"

#endif
