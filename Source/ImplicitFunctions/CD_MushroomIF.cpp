/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file    CD_MushroomIF.cpp
  @brief   Implementation of CD_MushroomIF.H
  @author  Robert Marskar
*/

// Chombo includes
#include <SmoothIntersection.H>

// Our includes
#include <CD_MushroomIF.H>
#include <CD_RoundedBoxIF.H>
#include <CD_BoxSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_NamespaceHeader.H>

MushroomIF::MushroomIF(const RealVect a_center,
                       const Real     a_R,
                       const Real     a_r,
                       const Real     a_L,
                       const Real     a_d,
                       const Real     a_curv,
                       const bool     a_fluidInside)
{
  const RealVect up = BASISREALV(SpaceDim - 1);

  // Make the two parts, the circular plate (head) and the support (stem).
  BaseIF* head = (BaseIF*)(new RoundedCylinderIF(a_center, a_center - a_d * up, a_R, a_curv, a_fluidInside));
  BaseIF* stem = (BaseIF*)(new CylinderSdf(a_center - 0.5 * a_d * up, a_center - (a_d + a_L) * up, a_r, a_fluidInside));

  // Smoosh them into a smooth intersection.
  Vector<BaseIF*> parts;
  parts.push_back(stem);
  parts.push_back(head);

  m_baseIF      = RefCountedPtr<BaseIF>(new SmoothIntersection(parts, a_curv));
  m_fluidInside = a_fluidInside;

  delete parts[0];
  delete parts[1];
}

MushroomIF::MushroomIF(const MushroomIF& a_inputIF)
{
  m_baseIF      = a_inputIF.m_baseIF;
  m_fluidInside = a_inputIF.m_fluidInside;
}

MushroomIF::~MushroomIF()
{}

Real
MushroomIF::value(const RealVect& a_pos) const
{
  Real retval = m_baseIF->value(a_pos);

  if (m_fluidInside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
MushroomIF::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new MushroomIF(*this));
}

#include <CD_NamespaceFooter.H>
