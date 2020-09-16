#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ANALYTICLONLATADVECTIONIBC_H_
#define _ANALYTICLONLATADVECTIONIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "PhysMappedIBC.H"
#include "RefCountedPtr.H"
#include "ScalarFunction.H"

#include "UsingNamespace.H"

//! \class AnalyticLonLatAdvectionIBC
//! This class uses analytic functions to set up initial conditions for
//! mapped advection.
class AnalyticLonLatAdvectionIBC: public PhysMappedIBC
{
  public:

  //! Constructs initial/boundary conditions for a mapped advection algorithm.
  explicit AnalyticLonLatAdvectionIBC(RefCountedPtr<ScalarFunction> a_initialState);

  //! Destructor
  ~AnalyticLonLatAdvectionIBC();

  // Overridden methods.

  PhysMappedIBC *new_physIBC();

  void initialize(LevelData<FArrayBox>& a_U);

  void primBC(FArrayBox&             a_WGdnv,
              const FArrayBox&       a_Wextrap,
              const FArrayBox&       a_W,
              const FArrayBox *const a_unitNormalBasisPtr,
              const Interval&        a_velIntv,
              const int&             a_dir,
              const Side::LoHiSide&  a_side,
              const Real&            a_time);

  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  void artViscBC(FArrayBox&                                   a_NtFdir,
                 const CHArray<Real, SpaceDim+1, ArRangeCol>& a_Nctg,
                 const FArrayBox&                             a_U,
                 const FArrayBox&                             a_unitNormalBasis,
                 const FArrayBox&                             a_divVel,
                 const FArrayBox&                             a_csq,
                 const FArrayBox&                             a_dxFace,
                 const Interval&                              a_momIntv,
                 const Real                                   a_alpha,
                 const Real                                   a_beta,
                 const Box&                                   a_loFaceBox,
                 const int                                    a_hasLo,
                 const Box&                                   a_hiFaceBox,
                 const int                                    a_hasHi,
                 const int                                    a_dir)
  {
    // No-op
  }

  private:

  RefCountedPtr<ScalarFunction> m_initialState;
};

#endif
