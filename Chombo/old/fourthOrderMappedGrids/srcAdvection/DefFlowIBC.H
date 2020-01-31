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
  virtual void initialize(LevelData<FArrayBox>& a_U);

  virtual void initializeWithJ(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values.
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
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
  virtual void setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual void artViscBC(FArrayBox&       a_F,
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

  virtual void initializeUnified(LevelData<FArrayBox>& a_U,
                                 bool a_includeJ) = 0;

  /// Period, in days (set to 5)
  Real m_period;

  /// Scalar multiple for velocity (set to 2)
  Real m_kappa;

  /// Initial time to evaluate profile (default 0)
  Real m_evalTime;
};

#include "NamespaceFooter.H"

#endif
