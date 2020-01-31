#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RAMPMAPPEDIBC_H_
#define _RAMPMAPPEDIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "PhysMappedIBC.H"

#include "UsingNamespace.H"

/// Example of a wedge.
/**
   Parameters:
     a_smallPressure - Lower limit for pressure (returned)
     a_gamma         - Gamma for polytropic, gamma-law gas
     a_ms            - Mach shock number
     a_X0Start       - Starting x-dir location of the shock in physical space
     a_artvisc       - Artificial viscosity coefficient
*/
class RampMappedIBC: public PhysMappedIBC
{
public:

  /// Null Constructor
  RampMappedIBC();

  /// Constructor which defines parameters used by Fortran routines
  RampMappedIBC(Real&      a_smallPressure,
                const Real a_gamma,
                const Real a_alpha,
                const Real a_ms,
                const Real a_X0Start,
                const Real a_artvisc);

  /// Destructor
  virtual ~RampMappedIBC();

  /// Sets parameters in a common block used by Fortran routines
  void setFortranCommon(Real&      a_smallPressure,
                        const Real a_gamma,
                        const Real a_alpha,
                        const Real a_ms,
                        const Real a_X0Start,
                        const Real a_artvisc);

  /// Factory method - this object is its own factory
  PhysMappedIBC* new_physIBC();

  /// Set up initial conditions
  virtual void initialize(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values.
  virtual void primBC(FArrayBox&             a_WGdnv,
                      const FArrayBox&       a_Wextrap,
                      const FArrayBox&       a_W,
                      const FArrayBox *const a_unitNormalBasisPtr,
                      const Interval&        a_velIntv,
                      const int&             a_dir,
                      const Side::LoHiSide&  a_side,
                      const Real&            a_time);

  /// Set boundary slopes
  virtual void setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
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
                 const int                                    a_dir);

protected:
  // True if the Fortran common block has been set
  static bool s_isFortranCommonSet;
  Real m_alpha;                       ///< Ramp angle in radians
  Real m_X0Start;                     ///< Starting x-Dir. location of the shock
};

#endif  // _RAMPMAPPEDIBC_H_
