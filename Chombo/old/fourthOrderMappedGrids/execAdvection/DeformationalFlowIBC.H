#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DEFORMATIONALFLOWIBC_H_
#define _DEFORMATIONALFLOWIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

/// Example of a cosine bell distribution
/**
   Parameters:
     a_ambientHeight  - Ambient height
     a_deltaHeight    - Height of the cosine bell
     a_center         - Center of the cosine bell (vector)
     a_size           - Standard deviation of the gaussian
     a_advectVelocity - Advection velocity
     a_advectAngle    - Advection angle
         a_evalTime       - Time to evaluate initial profile
 */
class DeformationalFlowIBC: public PhysAdvectMappedIBC
{
public:
  /// Null constructor
  /**
   */
  DeformationalFlowIBC();

  /// Constructor which defines parameters
  /**
   */
  DeformationalFlowIBC(
                    const Real&     a_ambientDensity,
                    const Real&     a_deltaDensity,
                    const RealVect& a_center,
                    const Real&     a_size,
                    const Real&     a_advectVelocity,
                    const Real&     a_advectAngle,
                    const Real&     a_evalTime);

  /// Destructor
  /**
   */
  virtual ~DeformationalFlowIBC();

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

  ///
  /* Return advection velocity.
   */
  virtual void advVel(LevelData<FArrayBox>& a_advVel,
                      Real a_time = 0.);

protected:

  void initializeUnified(LevelData<FArrayBox>& a_U,
                         bool a_includeJ);

  //  bool m_params_are_set;

  /// Ambient density of the advected field
  Real m_ambientDensity;

  /// Height of the cosine bell
  Real m_deltaDensity;

  /// Center of the cosine bell as (latitude, longitude)
  RealVect m_center;

  /// Initial time to evaluate profile
  Real m_evalTime;

  /// Size of the cosine bell
  Real m_size;

  /// Advection velocity
  Real m_advectVelocity;

  /// Advection angle
  Real m_advectAngle;
};

#include "NamespaceFooter.H"

#endif
