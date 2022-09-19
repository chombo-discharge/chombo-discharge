/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file  CD_CylinderSdf.cpp
  @brief CD_Implementation of CylinderSdf.H
  @author Robert Marskar
*/

// Chombo includes
#include <PolyGeom.H>

// Our includes
#include <CD_CylinderSdf.H>
#include <CD_NamespaceHeader.H>

CylinderSdf::CylinderSdf(const RealVect& a_center1,
                         const RealVect& a_center2,
                         const Real&     a_radius,
                         const bool&     a_fluidInside)
{
  m_endPoint1   = a_center1;
  m_endPoint2   = a_center2;
  m_center      = 0.5 * (m_endPoint1 + m_endPoint2);
  m_top         = m_endPoint2 - m_endPoint1;
  m_length      = m_top.vectorLength();
  m_axis        = m_top / m_length;
  m_radius      = a_radius;
  m_fluidInside = a_fluidInside;
}

CylinderSdf::CylinderSdf(const CylinderSdf& a_inputIF)
{
  m_endPoint1   = a_inputIF.m_endPoint1;
  m_endPoint2   = a_inputIF.m_endPoint2;
  m_center      = a_inputIF.m_center;
  m_top         = a_inputIF.m_top;
  m_length      = a_inputIF.m_length;
  m_axis        = a_inputIF.m_axis;
  m_radius      = a_inputIF.m_radius;
  m_fluidInside = a_inputIF.m_fluidInside;
}

Real
CylinderSdf::value(const RealVect& a_point) const
{

  // Translate cylinder center to origin and do all calculations for "positive half" of the cylinder.
  const RealVect newPoint  = a_point - m_center;
  const Real     paraComp  = PolyGeom::dot(newPoint, m_axis);
  const RealVect orthoVec  = newPoint - paraComp * m_axis;
  const Real     orthoComp = orthoVec.vectorLength();

  const Real f = orthoComp - m_radius; // Distance from cylinder wall.    Inside = Negative, Outside = Positive
  const Real g = abs(paraComp) -
                 0.5 * m_length; // Distance from cylinder end cap. Inside = Negative, Outside = Positive

  // This sets retval to be negative inside the cylinder, and positive outside.
  Real retval = 0.0;
  if (f <= 0. && g <= 0.) { // Point lies within the cylinder. Either short end or wall is closest point.
    retval = (abs(f) <= abs(g)) ? f : g;
  }
  else if (f <= 0. && g > 0.) { // Point lies inside radius but outside length. Short end is the closest point
    retval = g;
  }
  else if (f > 0. && g <= 0.) { // Point lies outside radius but inside length. Cylinder wall is the closest point.
    retval = f;
  }
  else if (f > 0. && g > 0.) { // Point lies outside both. The cylinder corner is the closest point.
    retval = sqrt(f * f + g * g);
  }
  else {
    MayDay::Error("CylinderSdf::value - logic bust, maybe center1 = center2?");
  }

  // The above value of retval is correct for m_fluidInside. Otherwise, flip it.
  if (!m_fluidInside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
CylinderSdf::newImplicitFunction() const
{
  return (BaseIF*)(new CylinderSdf(m_endPoint1, m_endPoint2, m_radius, m_fluidInside));
}

#include <CD_NamespaceFooter.H>
