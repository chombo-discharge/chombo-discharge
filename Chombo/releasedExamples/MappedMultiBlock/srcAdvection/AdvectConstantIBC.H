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
  virtual ~AdvectConstantIBC();

  // Overridden methods.

  PhysMappedIBC *new_physIBC();

  void initialize(LevelData<FArrayBox>& a_U);

  ///
  /* Return advection velocity.
   */
  virtual void advVel(FArrayBox& a_advVelFab,
                      Real a_time = 0.);

protected:

  //  bool m_params_are_set;

  /// Constant initial field
  Real m_fieldVal;

  /// Advection velocity
  Real m_advectVelocity;

  /// Advection angle
  Real m_advectAngle;

  LevelData<FArrayBox> m_advVel;

  /// Is m_advVel filled in?
  bool m_haveAdvVel;
};

#include "NamespaceFooter.H"

#endif
