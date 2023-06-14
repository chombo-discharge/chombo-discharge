/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Copyright © 2023 Robert Marskar
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeedleIF.cpp
  @brief  Implementation of CD_NeedleIF.H
  @author Fanny Skirbekk
  @author Robert Marskar
*/

// Our includes
#include <CD_Units.H>
#include <CD_NeedleIF.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

using Vec3 = EBGeometry::Vec3T<Real>;

NeedleIF::NeedleIF(const Real& a_length,
                   const Real& a_radius,
                   const Real& a_tipRadius,
                   const Real& a_angle,
                   const Real& a_cornerCurve,
                   const bool& a_flipInside)
{
  CH_TIME("NeedleIF::NeedleIF");

  CH_assert(a_radius > a_tipRadius);

  // a_angle is entire opening angle, dividing by two to get half of the opening angle
  const Real bodyRadius = a_radius - a_tipRadius;
  const Real tipLength  = bodyRadius / std::tan((0.5 * a_angle) * Units::pi / 180.0);

  // The center of the needle tip is set to origo in order for the rotation to work more easily
  const Vec3 bodyBackPosition(0.0, 0.0, -a_length + a_tipRadius);
  const Vec3 bodyFrontPosition(0.0, 0.0, -tipLength);
  const Vec3 needleTipCenter(0.0, 0.0, -a_tipRadius);

  std::shared_ptr<EBGeometry::ImplicitFunction<Real>> cone;
  std::shared_ptr<EBGeometry::ImplicitFunction<Real>> cylinder;
  std::shared_ptr<EBGeometry::ImplicitFunction<Real>> needle;

  cylinder = std::make_shared<EBGeometry::CapsuleSDF<Real>>(bodyBackPosition, bodyFrontPosition, bodyRadius);
  cone     = std::make_shared<EBGeometry::ConeSDF<Real>>(needleTipCenter, tipLength, a_angle);
  needle   = EBGeometry::SmoothUnion(cylinder, cone, a_cornerCurve);

  // Rotate so that geometry aligns with the y-axis.
  needle = EBGeometry::Rotate<Real>(needle, 90.0, 0);

  m_tipRadius        = a_tipRadius;
  m_implicitFunction = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(needle, !a_flipInside));
}

NeedleIF::NeedleIF(const NeedleIF& a_inputIF)
{
  CH_TIME("NeedleIF::NeedleIF(copy constructor)");

  this->m_implicitFunction = a_inputIF.m_implicitFunction;
  this->m_tipRadius        = a_inputIF.m_tipRadius;
}

Real
NeedleIF::value(const RealVect& a_point) const
{
  CH_TIME("NeedleIF::value");

  // round of the point of the needle with + m_tipRadius. Gives a curvature radius of m_tipRadius
  return m_implicitFunction->value(a_point) + m_tipRadius;
}

BaseIF*
NeedleIF::newImplicitFunction() const
{
  CH_TIME("NeedleIF::newImplicitFunction");

  return static_cast<BaseIF*>(new NeedleIF(*this));
}

#include <CD_NamespaceFooter.H>
