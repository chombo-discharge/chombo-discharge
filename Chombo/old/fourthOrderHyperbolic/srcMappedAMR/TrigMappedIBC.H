#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TRIGMAPPEDIBC_H_
#define _TRIGMAPPEDIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "PhysMappedIBC.H"

#include "UsingNamespace.H"

/// Example of a gaussian density distribution plus a constant density
/**
   Parameters:
     a_smallPressure  - Lower limit for pressure (returned)
     a_gamma          - Gamma for polytropic, gamma-law gas
     a_ambientDensity - Ambient density
     a_pressure       - If 0, use isentropic pressure
                        if 1, use the constant pressure of 1.0
     a_artvisc        - Artificial viscosity coefficient
 */
class TrigMappedIBC: public PhysMappedIBC
{
public:
  /// Null constructor
  /**
   */
  TrigMappedIBC();

  /// Constructor which defines parameters used by Fortran routines
  /**
   */
  TrigMappedIBC(Real&           a_smallPressure,
                const Real&     a_gamma,
                const Real&     a_ambientDensity,
                const Real&     a_artvisc);

  /// Destructor
  /**
   */
  virtual ~TrigMappedIBC();

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon(Real&           a_smallPressure,
                        const Real&     a_gamma,
                        const Real&     a_ambientDensity,
                        const Real&     a_artvisc);


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

  /// Set boundary primitive values.
  /**
   */
  void primBC(FArrayBox&             a_WGdnv,
              const FArrayBox&       a_Wextrap,
              const FArrayBox&       a_W,
              const FArrayBox *const a_unitNormalBasisPtr,
              const Interval&        a_velIntv,
              const int&             a_dir,
              const Side::LoHiSide&  a_side,
              const Real&            a_time);

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
                 const int                                    a_dir)
  {
    // No-op
  }

protected:
  // True if the Fortran common block has been set
  bool m_isFortranCommonSet;
};

#endif
