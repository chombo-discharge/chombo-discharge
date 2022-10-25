/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeedleIF.cpp
  @brief  Implementation of CD_NeedleIF.H
  @author Fanny Skirbekk
*/

// Chombo includes
#include <UnionIF.H>
#include <SmoothIntersection.H>

// Our includes
#include <CD_NeedleIF.H>
#include <CD_CylinderSdf.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

using Vec3 = EBGeometry::Vec3T<Real>;

NeedleIF::NeedleIF(const Real& a_length,
                   const Real& a_radius,
                   const bool& a_fluidInside,
                   const Real& a_tipRadius,
                   const Real& a_angle,
                   const Real& a_cornerCurve)
{
  m_tipRadius = a_tipRadius;

  constexpr Real pi = 3.14159265358979323846;
  // a_angle is entire opening angle, dividing by two to get half of the opening angle
  const Real tipLength = (a_radius - m_tipRadius) / std::tan((a_angle / 2) * pi / 180); //double

  RealVect centerFront(CH_SPACEDIM), centerBack(CH_SPACEDIM);
  centerFront[1] = tipLength;
  centerBack[1]  = a_length;

  // the center of the needle tip is set to origo in order for the rotation to work more easily
  const Vec3 centerT(0.0, 0.0, 0.0);

  // Build the needle-parts:
  Vector<BaseIF*> isects;
  isects.push_back(
    static_cast<BaseIF*>(new CylinderSdf(centerFront, centerBack, (a_radius - m_tipRadius), a_fluidInside)));

  // flipinside=true for cone since EBGeometry and Chombo has opposing sign conventions/logic regarding the flipinside..
  // the center of the needle tip is set to origo in order for the rotation to work more easily
  auto cone = std::make_shared<EBGeometry::ConeSDF<Real>>(centerT, tipLength, a_angle, true);

  // Cone rotation will only work as expected if the cone tip is placed in origo.
  cone->rotate(90, 0);

  isects.push_back(static_cast<BaseIF*>(new EBGeometryIF(cone, false)));

  // Build the needle
  m_baseif = RefCountedPtr<BaseIF>(new SmoothIntersection(isects, a_cornerCurve));

  // Delete everything we have allocated so far
  for (int i = 0; i < isects.size(); ++i) {
    delete isects[i];
  }
}

NeedleIF::NeedleIF(const NeedleIF& a_inputIF)
{
  this->m_baseif    = a_inputIF.m_baseif;
  this->m_tipRadius = a_inputIF.m_tipRadius;
}

Real
NeedleIF::value(const RealVect& a_point) const
{
  // round of the point of the needle with + m_tipRadius. Gives a curvature radius of m_tipRadius
  return m_baseif->value(a_point) + m_tipRadius;
}

BaseIF*
NeedleIF::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new NeedleIF(*this));
}

#include <CD_NamespaceFooter.H>
