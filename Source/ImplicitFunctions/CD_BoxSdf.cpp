/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxSdf.cpp
  @brief  Implements CD_BoxSdf.H
  @author Robert Marskar
*/

// Chombo includes
#include <PlaneIF.H>
#include <SmoothUnion.H>
#include <UnionIF.H>
#include <IntersectionIF.H>

// our includes
#include <CD_BoxSdf.H>
#include <CD_NamespaceHeader.H>

BoxSdf::BoxSdf(const RealVect& a_loCorner, const RealVect& a_hiCorner, const bool& a_fluidInside)
{
  m_loCorner    = a_loCorner;
  m_hiCorner    = a_hiCorner;
  m_fluidInside = a_fluidInside;
}

BoxSdf::BoxSdf(const BoxSdf& a_inputIF)
{
  m_loCorner    = a_inputIF.m_loCorner;
  m_hiCorner    = a_inputIF.m_hiCorner;
  m_fluidInside = a_inputIF.m_fluidInside;
}

BoxSdf::~BoxSdf()
{}

Real
BoxSdf::value(const RealVect& a_pos) const
{
  // TLDR: Min(0.0, Max(dx, dy)) is the shortest distance from a point inside the rectangle to one of the edges. It becomes zero if dx > 0 or dy > 0. The
  //       second term max(RealVect::Zero, delta).vectorLength() is funky. In principle we should only take the distance to be sqrt(dx^2 + dy^2 + dz^2) if
  //       we are closest to a "corner", but if we are closest to one of the edges we have dx*dy < 0.0. We can just set the other component to zero, which is what
  //       max(RealVect::Zero, v) does. Then take the distance anyways.

  const RealVect delta = RealVect(D_DECL(Max(m_loCorner[0] - a_pos[0], a_pos[0] - m_hiCorner[0]),
                                         Max(m_loCorner[1] - a_pos[1], a_pos[1] - m_hiCorner[1]),
                                         Max(m_loCorner[2] - a_pos[2], a_pos[2] - m_hiCorner[2])));

  Real retval = Min((Real)0.0, delta[delta.maxDir(false)]) +
                max(RealVect::Zero, delta)
                  .vectorLength(); // Negative inside box. NOTE: This is RealVect::max and not std::max

  if (!m_fluidInside) { // Flip so lsf is positive inside box and negative outside box.
    retval = -retval;
  }

  return retval;
}

BaseIF*
BoxSdf::newImplicitFunction() const
{
  return (BaseIF*)(new BoxSdf(m_loCorner, m_hiCorner, m_fluidInside));
}

#include <CD_NamespaceFooter.H>
