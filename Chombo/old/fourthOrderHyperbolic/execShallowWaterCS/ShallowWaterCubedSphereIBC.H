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
#include "Vector.H"
#include "RealVect.H"
#include "ScalarFunction.H"
#include "VectorFunction.H"
#include "RefCountedPtr.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

class ShallowWaterCubedSphereIBC: public PhysMappedIBC
{
  public:


  //! Constructor with initial solution, ShallowWater velocity.
  /**
   */
  ShallowWaterCubedSphereIBC(RefCountedPtr<ScalarFunction> a_height,
                             RefCountedPtr<VectorFunction> a_velocity);

  //! Destructor
  ~ShallowWaterCubedSphereIBC();

  // Overridden methods.

  PhysMappedIBC *new_physIBC(); // Factory method.

  void initialize(LevelData<FArrayBox>& a_U);

  void initializeWithJ(LevelData<FArrayBox>& a_U);

  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  protected:

  void initializeUnified(LevelData<FArrayBox>& a_U,
                         bool a_includeJ);

  // Functions.
  RefCountedPtr<ScalarFunction> m_height;
  RefCountedPtr<VectorFunction> m_velocity;

  private:

  // Forbidden.
  ShallowWaterCubedSphereIBC();
  ShallowWaterCubedSphereIBC(const ShallowWaterCubedSphereIBC&);
  ShallowWaterCubedSphereIBC& operator=(const ShallowWaterCubedSphereIBC&);
};

#include "NamespaceFooter.H"

#endif
