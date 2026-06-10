/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_HollowCylinderIF.cpp
  @brief  Implementation of CD_HollowCylinderIF.H
  @author Robert Marskar
*/

// Chombo includes
#include <SmoothUnion.H>
#include <UnionIF.H>
#include <IntersectionIF.H>

// Our includes
#include <CD_CylinderSdf.H>
#include <CD_HollowCylinderIF.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_NamespaceHeader.H>

HollowCylinderIF::HollowCylinderIF(const RealVect& a_center1,
                                   const RealVect& a_center2,
                                   const Real      a_majorRadius,
                                   const Real      a_minorRadius,
                                   const Real      a_outerCurvature,
                                   const Real      a_innerCurvature,
                                   const bool      a_fluidInside)
{

  // Make the SmoothUnion of this stuff.
  Vector<BaseIF*> parts;

  RealVect axis = (a_center2 - a_center1);
  axis          = axis / axis.vectorLength();

  auto*   bigCylinder   = (BaseIF*)(new RoundedCylinderIF(a_center1,
                                                      a_center2,
                                                      a_majorRadius,
                                                      a_outerCurvature,
                                                      a_fluidInside));
  BaseIF* smallCylinder = (BaseIF*)(new CylinderSdf(a_center1, a_center2, a_minorRadius, !a_fluidInside));

  parts.push_back(bigCylinder);
  parts.push_back(smallCylinder);

  // Make union
  m_baseIF = RefCountedPtr<BaseIF>(new SmoothUnion(parts, a_innerCurvature));
}

HollowCylinderIF::HollowCylinderIF(const HollowCylinderIF& a_inputIF) : m_baseIF(a_inputIF.m_baseIF)
{}

Real
HollowCylinderIF::value(const RealVect& a_point) const
{
  return m_baseIF->value(a_point);
}

BaseIF*
HollowCylinderIF::newImplicitFunction() const
{
  return (BaseIF*)(new HollowCylinderIF(*this));
}

#include <CD_NamespaceFooter.H>
