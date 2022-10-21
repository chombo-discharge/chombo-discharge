/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   PerlinRodSdf.cpp
  @brief  Implementation of CD_PerlinRodSdf.H
  @author Robert Marskar
*/

// Chombo includes
#include <TransformIF.H>
#include <IntersectionIF.H>

// Our includes
#include <CD_CylinderSdf.H>
#include <CD_PerlinRodSdf.H>
#include <CD_NamespaceHeader.H>

PerlinRodSdf::PerlinRodSdf(const Real&     a_rad,
                           const RealVect& a_center1,
                           const RealVect& a_center2,
                           const bool&     a_inside,
                           const Real&     a_noiseAmp,
                           const RealVect& a_noiseFreq,
                           const Real&     a_persistence,
                           const int&      a_octaves,
                           const bool&     a_reseed)
{

  // Fix up center2
  const RealVect axis    = a_center2 - a_center1;
  const RealVect center2 = a_center2 - axis * a_rad / axis.vectorLength();

  // Cylinder and graded noise sphere
  BaseIF*       cyl = static_cast<BaseIF*>(new CylinderSdf(a_center1, center2, a_rad, a_inside));
  const BaseIF* sph = static_cast<BaseIF*>(new GradedPerlinSphereSdf(a_rad + a_noiseAmp,
                                                                     center2,
                                                                     a_inside,
                                                                     a_noiseAmp,
                                                                     a_noiseFreq,
                                                                     a_persistence,
                                                                     a_octaves,
                                                                     a_reseed));

  // Rotate sph so that it aligns with center2 - a_center1
#if CH_SPACEDIM == 2
  const Real   theta = atan2(axis[0], axis[1]);
  TransformIF* trans = new TransformIF(*sph);
  trans->rotate(theta, center2);
#elif CH_SPACEDIM == 3 // By default, the sphere noise points along the negative SpaceDim  direction. Rotate
  // that vector to center2-center1
  TransformIF*   trans = new TransformIF(*sph);
  const RealVect vec   = RealVect(1.E-8, 0.0, 1); // Need something small so that axis is not already parallel
  trans->rotate(vec, center2 - a_center1, center2);
#endif

  // By default, sphere noise is along the SpaceDim direction. Rotate that vector to center2 - center1
  Vector<BaseIF*> parts(2);
  parts[0] = cyl;
  parts[1] = trans;
  m_baseif = RefCountedPtr<BaseIF>(new IntersectionIF(parts));

  delete cyl;
  delete sph;
  delete trans;
}

PerlinRodSdf::PerlinRodSdf(const PerlinRodSdf& a_inputIF) : PerlinSphereSdf(a_inputIF)
{
  m_baseif = a_inputIF.m_baseif;
}

PerlinRodSdf::~PerlinRodSdf() {}

Real
PerlinRodSdf::value(const RealVect& a_pos) const
{
  return m_baseif->value(a_pos);
}

BaseIF*
PerlinRodSdf::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new PerlinRodSdf(*this));
}

#include <CD_NamespaceFooter.H>
