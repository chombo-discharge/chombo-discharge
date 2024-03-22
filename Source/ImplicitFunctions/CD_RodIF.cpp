/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodIF.cpp
  @brief  Implementation of CD_RodIF.H
  @author Robert Marskar
*/

// Chombo includes
#include <SphereIF.H>
#include <IntersectionIF.H>
#include <PolyGeom.H>

// Our includes
#include <CD_RodIF.H>
#include <CD_SphereSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_NamespaceHeader.H>

RodIF::RodIF(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_fluidInside)
{
  const RealVect axis    = (a_center2 - a_center1);
  const RealVect axisVec = axis / axis.vectorLength();

  // Find two new centers where can place cylinder edges and spheres
  const RealVect c1 = a_center1 + axisVec * a_radius;
  const RealVect c2 = a_center2 - axisVec * a_radius;

  // Build the cylinder
  Vector<BaseIF*> isects;
  isects.push_back(static_cast<BaseIF*>(new CylinderSdf(c1, c2, a_radius, a_fluidInside)));
  isects.push_back(static_cast<BaseIF*>(new SphereSdf(c1, a_radius, a_fluidInside)));
  isects.push_back(static_cast<BaseIF*>(new SphereSdf(c2, a_radius, a_fluidInside)));

  // Build the rod
  m_baseif = RefCountedPtr<BaseIF>(new IntersectionIF(isects));

  // Delete everything we allocated so far
  for (int i = 0; i < isects.size(); i++) {
    delete isects[i];
  }
}

RodIF::RodIF(const RodIF& a_inputIF)
{
  this->m_baseif = a_inputIF.m_baseif;
}

Real
RodIF::value(const RealVect& a_point) const
{
  return m_baseif->value(a_point);
}

BaseIF*
RodIF::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new RodIF(*this));
}

#include <CD_NamespaceFooter.H>
