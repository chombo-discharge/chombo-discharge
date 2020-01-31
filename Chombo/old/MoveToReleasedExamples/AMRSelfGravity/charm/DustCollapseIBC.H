#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DUSTCOLLAPSEIBC_H_
#define _DUSTCOLLAPSEIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "PhysIBC.H"


///
/**
   General input for Hydro:
     a_gamma       - Gamma for polytropic, gamma-law gas
     a_artvisc     - Artificial viscosity coefficient
     a_rsTolerance - tolerance for the Riemann solver
     a_max_Rs_iter - max number of iteration in the Riemann solver
     a_maxMach     - max Mach number to switch from Etot to entropy (S) eq.

   Initial conditions parameters
     a_cloudRadius    - Initial radius of the dust sphere
     a_cloudDensity   - Initial density of the dust sphere

 */
class DustCollapseIBC : public PhysIBC
{
public:
  /// Null constructor
  /**
   */
  DustCollapseIBC();

  /// Constructor which defines parameters used by Fortran routines
  /**
   */
  DustCollapseIBC(const Real&  a_gamma,
               const Real&  a_radius,
               const Real&  a_density,
               const Real&  a_artvisc,
               const Real&  a_rsTolerance,
               const Real&  a_maxRsIter,
                  const Real&  a_maxMach);

  /// Destructor
  /**
   */
  ~DustCollapseIBC()
  {
  }

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon(const Real&  a_gamma,
                        const Real&  a_radius,
                        const Real&  a_density,
                        const Real&  a_artvisc,
                        const Real&  a_rsTolerance,
                        const Real&  a_maxRsIter,
                        const Real&  a_maxMach);

  /// Set the flag m_isFortranCommonSet to true
  /**
     Set the flag m_isFortranCommonSet to true so that new IBCs made with
     new_physIBC() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setFortranCommonSet();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysIBC *new_physIBC();

  /// Set boundary primitive values.
  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  /// Set boundary slopes
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  /// Set up initial conditions for gas
  void initialize(LevelData<FArrayBox>& a_U);

protected:
  // cloud density and radius (for spher. dust collapse and makes central
  // mass for the secondary infall)
  Real m_density, m_radius;

  // True if the Fortran common block has been set
  bool m_isFortranCommonSet;
};

#endif
