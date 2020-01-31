#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COSINEBELLSFLOWIBC_H_
#define _COSINEBELLSFLOWIBC_H_

#include <iostream>

#include "DefFlowIBC.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"

#include "NamespaceHeader.H"

/// Example of a cosine bell distribution
/**
   Parameters:
     a_hmax        - Maximum height of the cosine bells (set to 1)
     a_radius      - Radius of the cosine bells (set to 1/2)
     a_longitude   - Longitudes of centers of cosine bells (5*pi/6, 7*pi/6)
     a_latitude    - Latitudes of centers of cosine bells (0, 0)
     a_background  - Background phi (set to 0.1)
     a_delta       - Difference in phi (set to 0.9)
     a_period      - Period of flow, in days (set to 5)
     a_kappa       - Scalar multiple for velocity (set to 2)
     a_evalTime    - Initial time to evaluate profile (default 0)
 */
class CosineBellsFlowIBC: public DefFlowIBC
{
public:
  /// Null constructor
  /**
   */
  CosineBellsFlowIBC();

  /// Constructor which defines parameters
  /**
   */
  CosineBellsFlowIBC(
                     const Real&     a_hmax,
                     const Real&     a_radius,
                     const Vector<Real>&  a_longitude,
                     const Vector<Real>&  a_latitude,
                     const Real&     a_background,
                     const Real&     a_delta,
                     const Real&     a_period,
                     const Real&     a_kappa,
                     const Real&     a_evalTime = 0.);

  /// Destructor
  /**
   */
  virtual ~CosineBellsFlowIBC();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysMappedIBC *new_physIBC();

protected:

  void initializeUnified(LevelData<FArrayBox>& a_U,
                         bool a_includeJ);

  //  bool m_params_are_set;

  /// Maximum height of the cosine bells (set to 1)
  Real m_hmax;

  /// Radius of the cosine bells (set to 1/2)
  Real m_radius;

  /// Longitudes of centers of cosine bells (5*pi/6, 7*pi/6)
  Vector<Real> m_longitude;

  /// Latitudes of centers of cosine bells (0, 0)
  Vector<Real> m_latitude;

  /// Background phi (set to 0.1)
  Real m_background;

  /// Difference in phi (set to 0.9)
  Real m_delta;
};

#include "NamespaceFooter.H"

#endif
