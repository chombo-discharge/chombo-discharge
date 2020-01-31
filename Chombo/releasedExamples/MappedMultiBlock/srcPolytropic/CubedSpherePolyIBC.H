#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHEREPOLYIBC_H__
#define _CUBEDSPHEREPOLYIBC_H__

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "CubedSphereShellPanelCS.H"
#include "PhysMappedIBC.H"

///  solid body rotation with cosine bewll denisty perturbation
/**
   Parameters:
   a_smallPressure - Lower limit for pressure (returned)
   a_gamma         - Gamma for polytropic, gamma-law gas
   a_rho0          -base line density
   a_p0            -base line pressure
   a_hmax        - Maximum height of the cosine bells (set to 1)
   a_radius      - Radius of the cosine bells (set to 1/2)
   a_longitude   - Longitudes of centers of cosine bells (5*pi/6, 7*pi/6)
   a_latitude    - Latitudes of centers of cosine bells (0, 0)
   a_delta       - Difference in rho (set to 0.9)
*/
class CubedSpherePolyIBC: public PhysMappedIBC
{
public:

  /// Constructor which defines parameters used by Fortran routines
  CubedSpherePolyIBC(const Real a_gamma,
                     const Real a_r0,
                     const Real a_p0,
                     const Real&     a_radius,
                     const Vector<Real>&  a_longitude,
                     const Vector<Real>&  a_latitude,
                     const Real&     a_delta,
                     const int       a_power,
                     const Real& a_magnitude, 
                     const Real& a_angle, 
                     const Real& a_shift);

  /// Destructor
  virtual ~CubedSpherePolyIBC();

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
  /// Null Constructor
  CubedSpherePolyIBC();

  void getAdvVel(LevelData<FArrayBox>& a_advVel, Real a_time);

  void pointAdvVel(FArrayBox& a_velFab,
                   const FArrayBox& a_xiLoc, Real a_time, 
                   const CubedSphereShellPanelCS& coordSysBlock);
  
  /// Sets parameters in a common block used by Fortran routines
  void setFortranCommon(Real&      a_smallPressure,
                        const Real a_gamma);

  // True if the Fortran common block has been set
  static bool s_isFortranCommonSet;


  Real m_gamma;

  /// Background state
  Real m_r0;
  Real m_p0;
  /// Radius of the cosine bells (set to 1/2)
  Real m_radius;

  /// Longitudes of centers of cosine bells (5*pi/6, 7*pi/6)
  Vector<Real> m_longitude;

  /// Latitudes of centers of cosine bells (0, 0)
  Vector<Real> m_latitude;

  /// Difference in rho0
  Real m_delta;

  /// The power for the cosine bell (default 1)

  //stuff from solid body rotation class
  int  m_power;
  Real m_magnitude;
  Real m_angle;
  Real m_shift;

};

#endif  // _CUBEDSPHEREPOLYIBC_H_
