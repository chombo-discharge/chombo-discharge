#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RTHETAPHICS_H_
#define _RTHETAPHICS_H_

#include "CONSTANTS.H"
#include "NewFourthOrderCoordSys.H"
#include "NamespaceHeader.H"

/// Mapping from spherical coordinates (Xi) to real coodinates.
/**
 *  <ul>
 *    <li> \f$\theta\f$ is the azimuthal angle (counter-clockwise from +x axis
 *         in xy plane)
 *    <li> \f$\varphi\f$ is the inclination angle measure from +z axis
 *  </ul>
 *
 *  \f$ x = r\sin\varphi\cos\theta \f$<br/>
 *  \f$ y = r\sin\varphi\sin\theta \f$<br/>
 *  \f$ z = r\cos\varphi \f$
 *
 *  <p> But to avoid singularities, we define<br/>
 *  \f$ r = \mbox(stretch)(\hat(r) + r_(\min)) \f$<br/>
 *  \f$ \theta = 2\pi\hat(\theta) \f$<br/>
 *  \f$ \varphi = \hat(\varphi)(\pi - 2\varphi_(\min)) + \varphi_(\min) \f$<br/>
 *  where \f$0 <= \hat(r), \hat(\theta), \hat(\varphi) <= 1 \f$.  stretch,
 *  \f$r_(\min)\f$, and \f$\varphi_(\min)\f$ are defined by input.  Note that
 *  \f$\varphi_(\min)\f$ is the minimum angle with the z-axis.
 *
 */

class RThetaPhiCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  RThetaPhiCS(const RealVect& a_dx,
              const Real      a_stretch = 1.,
              const Real      a_rMin    = 1.,
              const Real      a_phiMin  = Pi/4.);

  /// Destructor
  virtual ~RThetaPhiCS();

  /// Given coordinate in mapped space, return its location in real space
  virtual RealVect realCoord(const RealVect& a_hatXi) const;

  /// Given coordinate in real space, return its location in the mapped space
  virtual RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

private:
  /// Normalize the mapped coordinates \f$(0 <= \hat{\xi} <= 1)\f$
  RealVect normalize(const RealVect & a_Xi) const;

  /// Unnormalized the mapped coordinates
  RealVect unnormalize(const RealVect & a_hatXi) const;

  /// Return location in real space of an unnormalized Coord. from mapped space
  RealVect realCoordUN(const RealVect& a_Xi) const;

protected:
  Real m_stretch;                   ///< Stretching in the r direction
  Real m_rMin;                      ///< Min r (> 0 avoids singularity there)
  Real m_phiMin;                    ///< Min angle between r and z-axis. (>0
                                    ///< avoids singularity there)
};


/// Factory for RThetaPhiCS
/**
 */
class
RThetaPhiCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  RThetaPhiCSFactory(const Real a_stretch = 1.,
                     const Real a_rMin = 1.,
                     const Real a_phiMin = Pi/4.);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:

  ///
  Real m_stretch;
  Real m_rMin;
  Real m_phiMin;
};

#include "NamespaceFooter.H"

#endif
