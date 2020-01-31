#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHADVECTMAPPEDOPERATOR_H_
#define _PATCHADVECTMAPPEDOPERATOR_H_

#include "PatchMappedConsOperator.H"

#include "NamespaceHeader.H"

/// Patch Cons Operator - a class for evaluating the Laplacian.

/**
 */
class PatchAdvectMappedOperator : public PatchMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  PatchAdvectMappedOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  virtual ~PatchAdvectMappedOperator();

  virtual void preRiemann(FArrayBox&  a_WLeft,
                          FArrayBox&  a_WRight,
                          int         a_dir,
                          const Box&  a_box);

  virtual void postRiemann(FArrayBox&  a_Wface,
                           int         a_dir,
                           const Box&  a_box);

  void setAdvVelFace(const FluxBox* a_advVelFacePtr);

protected:

  virtual void getAllFluxes(FluxBox&        a_FfaceAvg,
                            FluxBox&        a_FfaceCen,
                            const FluxBox&  a_WfaceAvg,
                            const FluxBox&  a_WfaceCen);

  FluxBox* m_advVelFacePtr;
};

#include "NamespaceFooter.H"

#endif
