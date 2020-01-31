#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ROSSBYHAURWITZIBC_H_
#define _ROSSBYHAURWITZIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"
#include "PhysShallowWaterMappedIBC.H"

#include "NamespaceHeader.H"

/// Rossby-Haurwitz wave.
/**
   This test was taken from Case 6 in Williamson et. al.,
   "A Standard Test Set for Numerical Approximations to the Shallow Water
   Equations in Spherical Geometry", JCP 102:211-224 (1992).
 */
class RossbyHaurwitzIBC: public PhysShallowWaterMappedIBC
{
public:
  /// Null constructor
  /**
   */
  RossbyHaurwitzIBC();

  /// Constructor which defines parameters
  /**
   */
  RossbyHaurwitzIBC(
                    Real a_w,
                    Real a_K,
                    Real a_h0,
                    int a_R);

  /// Destructor
  /**
   */
  virtual ~RossbyHaurwitzIBC();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysMappedIBC *new_physIBC();

  /// Set up initial conditions
  /**
   */
  void initializeUnified(LevelData<FArrayBox>& a_U,
                         bool a_multiplyJ);

  void initialize(LevelData<FArrayBox>& a_U);

  void initializeWithJ(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values.
  /**
   */
  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

protected:

  Real m_w;

  Real m_K;

  Real m_h0;

  int m_R;
};

#include "NamespaceFooter.H"

#endif
