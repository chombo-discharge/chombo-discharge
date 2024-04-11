/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TorusSdf.cpp
  @brief  Implementation of CD_TorusSdf.H
  @author Robert Marskar
*/

// Our includes
#include <CD_TorusSdf.H>
#include <CD_NamespaceHeader.H>

TorusSdf::TorusSdf(const RealVect a_center,
                   const Real     a_majorRadius,
                   const Real     a_minorRadius,
                   const bool     a_fluidInside)
{
  m_center      = a_center;
  m_majorRadius = a_majorRadius;
  m_minorRadius = a_minorRadius;
  m_fluidInside = a_fluidInside;
}

TorusSdf::TorusSdf(const TorusSdf& a_inputIF)
{
  m_center      = a_inputIF.m_center;
  m_majorRadius = a_inputIF.m_majorRadius;
  m_minorRadius = a_inputIF.m_minorRadius;
  m_fluidInside = a_inputIF.m_fluidInside;
}

TorusSdf::~TorusSdf()
{}

Real
TorusSdf::value(const RealVect& a_point) const
{
  const RealVect p = a_point - m_center;

  Real radius = 0.0;
  for (int dir = 0; dir < 2; dir++) {
    const Real cur = p[dir];
    radius += cur * cur;
  }
  radius = sqrt(radius) - m_majorRadius;

  Real retval = radius * radius;
#if CH_SPACEDIM == 3
  retval += p[2] * p[2];
#endif

  retval = sqrt(retval) - m_minorRadius; // Positive outside.

  if (!m_fluidInside) { // Make sure negative outside, if fluid is outside.
    retval = -retval;
  }

  return retval;
}

BaseIF*
TorusSdf::newImplicitFunction() const
{
  return (BaseIF*)new TorusSdf(m_center, m_majorRadius, m_minorRadius, m_fluidInside);
}

#include <CD_NamespaceFooter.H>
