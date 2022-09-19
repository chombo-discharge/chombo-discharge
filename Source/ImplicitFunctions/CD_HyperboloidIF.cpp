/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_HyperboloidIF.cpp
  @brief  Implementation CD_HyperboloidIF.H
  @author Robert Marskar
*/

// Our includes
#include <CD_HyperboloidIF.H>
#include <CD_NamespaceHeader.H>

HyperboloidIF::HyperboloidIF(const RealVect& a_radii, const RealVect& a_center, const bool& a_inside)
{
  m_radii  = a_radii;
  m_center = a_center;
  m_inside = a_inside;
  m_radii2 = m_radii * m_radii;
  m_sign   = RealVect::Unit;

  m_sign[SpaceDim - 1] = -1.;
}

HyperboloidIF::HyperboloidIF(const HyperboloidIF& a_inputIF)
{
  m_radii  = a_inputIF.m_radii;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;
  m_radii2 = a_inputIF.m_radii2;
  m_sign   = a_inputIF.m_sign;
}

Real
HyperboloidIF::value(const RealVect& a_point) const
{

  Real retval = 1.;

  for (int dir = 0; dir < SpaceDim - 1; dir++) {
    Real cur;
    cur = a_point[dir] - m_center[dir];
    retval += m_sign[dir] * cur * cur / m_radii2[dir];
  }

  retval = (a_point[SpaceDim - 1] - m_center[SpaceDim - 1] + m_radii[SpaceDim - 1]) -
           m_radii[SpaceDim - 1] * sqrt(retval);

  // Change sign if outside
  if (m_inside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
HyperboloidIF::newImplicitFunction() const
{
  HyperboloidIF* hyperboloidPtr = new HyperboloidIF(m_radii, m_center, m_inside);

  return static_cast<BaseIF*>(hyperboloidPtr);
}

#include <CD_NamespaceFooter.H>
