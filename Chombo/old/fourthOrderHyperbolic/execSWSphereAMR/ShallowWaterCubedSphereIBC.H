#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SHALLOWWATERCUBEDSPHEREIBC_H_
#define _SHALLOWWATERCUBEDSPHEREIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "PhysMappedIBC.H"
#include "RefCountedPtr.H"
#include "ScalarFunction.H"
#include "VectorFunction.H"
// #include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

class ShallowWaterCubedSphereIBC: public PhysMappedIBC
{
public:


  //! Constructor with initial solution, ShallowWater velocity.
  /**
   */
  explicit ShallowWaterCubedSphereIBC(RefCountedPtr<ScalarFunction> a_height,
                                      RefCountedPtr<VectorFunction> a_velocity,
                                      Real a_gravity,
                                      Real a_omega,
                                      Real a_alpha);

  //! Destructor
  ~ShallowWaterCubedSphereIBC();

  // Overridden methods.

  PhysMappedIBC *new_physIBC(); // Factory method.

  void initialize(LevelData<FArrayBox>& a_U);

  // void initializeWithJ(LevelData<FArrayBox>& a_U);

  void primBC(FArrayBox&             a_WGdnv,
              const FArrayBox&       a_Wextrap,
              const FArrayBox&       a_W,
              const FArrayBox *const a_unitNormalBasisPtr,
              const Interval&        a_velIntv,
              const int&             a_dir,
              const Side::LoHiSide&  a_side,
              const Real&            a_time);

  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time)
  {
    // no-op, because there are no boundaries
  };

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
                 const int                                    a_dir);

protected:

  //  void initializeUnified(LevelData<FArrayBox>& a_U,
  //                         bool a_includeJ);

  // Functions.
  RefCountedPtr<ScalarFunction> m_height;
  RefCountedPtr<VectorFunction> m_velocity;

  Real m_gravity;
  Real m_omega;
  Real m_alpha;

private:

  // Forbidden.
  //  ShallowWaterCubedSphereIBC();
  //  ShallowWaterCubedSphereIBC(const ShallowWaterCubedSphereIBC&);
  //  ShallowWaterCubedSphereIBC& operator=(const ShallowWaterCubedSphereIBC&);
};

#include "NamespaceFooter.H"

#endif
