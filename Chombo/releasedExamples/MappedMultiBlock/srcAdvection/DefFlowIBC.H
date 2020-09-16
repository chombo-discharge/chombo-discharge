#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DEFFLOWIBC_H_
#define _DEFFLOWIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

/// Example of a deformational flow
class DefFlowIBC: public PhysAdvectMappedIBC
{
public:
  /// Null constructor
  /**
   */
  DefFlowIBC();

  /// Constructor which defines parameters
  /**
   */
  DefFlowIBC(
             const Real&  a_period,
             const Real&  a_kappa,
             const Real&  a_evalTime);

  /// Destructor
  /**
   */
  virtual ~DefFlowIBC();

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

  ///
  /* Return advection velocity.
   */
  virtual void advVel(FArrayBox& a_advVel,
                      Real a_time = 0.);

protected:

  /// Period, in days (set to 5)
  Real m_period;

  /// Scalar multiple for velocity (set to 2)
  Real m_kappa;

  /// Initial time to evaluate profile (default 0)
  Real m_evalTime;
};

#include "NamespaceFooter.H"

#endif
