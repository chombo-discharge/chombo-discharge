#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BAROTROPICINSTABILITYIBC_H_
#define _BAROTROPICINSTABILITYIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"
#include "PhysShallowWaterMappedIBC.H"

#include "NamespaceHeader.H"

/// barotropic instability test case
/**
   This test was taken from Galewsky, Scott, and Polvani,
   "An initial-value problem for testing numerical models of the global
   shallow-water equations", Tellus 56A:429-440 (2004).
 */
class BarotropicInstabilityIBC: public PhysShallowWaterMappedIBC
{
public:
  /// Null constructor
  /**
   */
  BarotropicInstabilityIBC();

  /// Constructor which defines parameters
  /**
   */
  BarotropicInstabilityIBC(
                           Real a_phi0,
                           Real a_phi1,
                           Real a_phi2,
                           Real a_umax,
                           Real a_hmean,
                           Real a_hhat,
                           Real a_BIalpha,
                           Real a_BIbeta);

  /// Destructor
  /**
   */
  virtual ~BarotropicInstabilityIBC();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysMappedIBC *new_physIBC();

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  virtual void setFortranCommon(const Real& a_gravity,
                                const Real& a_omega,
                                const Real& a_alpha);

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

  Real m_phi0;
  Real m_phi1;
  Real m_phi2;
  Real m_umax;
  Real m_hmean;
  Real m_hhat;
  Real m_BIalpha;
  Real m_BIbeta;
};

#include "NamespaceFooter.H"

#endif
