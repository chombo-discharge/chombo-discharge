#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTCONSTANTIBC_H_
#define _ADVECTCONSTANTIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

/// Example of a constant distribution
/**
   Parameters:
     a_fieldVal       - Constant initial value of the advected quantity
     a_advectVelocity - Advection velocity
     a_advectAngle    - Advection angle
 */
class AdvectConstantIBC: public PhysAdvectMappedIBC
{
public:
  /// Null constructor
  /**
   */
  AdvectConstantIBC();

  /// Constructor which defines parameters
  /**
   */
  AdvectConstantIBC(
                    Real a_fieldVal,
                    Real a_advectVelocity,
                    Real a_advectAngle);

  /// Destructor
  /**
   */
  virtual ~AdvectConstantIBC();

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

  /// Constant initial field
  Real m_fieldVal;

  /// Advection velocity
  Real m_advectVelocity;

  /// Advection angle
  Real m_advectAngle;

  LevelData<FArrayBox> m_advVel;

  /// Is m_advVel filled in?
  bool m_gotAdvVel;
};

#include "NamespaceFooter.H"

#endif
