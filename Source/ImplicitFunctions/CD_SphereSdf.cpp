/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_SphereSdf.cpp
  @brief  Implementation of CD_SphereSdf.H
  @author Robert Marskar
*/

// Our includes
#include <CD_SphereSdf.H>
#include <CD_NamespaceHeader.H>

SphereSdf::SphereSdf(const RealVect& a_center, const Real& a_radius, const bool& a_fluidInside)
{
  m_center      = a_center;
  m_radius      = a_radius;
  m_fluidInside = a_fluidInside;
}

SphereSdf::SphereSdf(const SphereSdf& a_inputIF)
{
  m_center      = a_inputIF.m_center;
  m_radius      = a_inputIF.m_radius;
  m_fluidInside = a_inputIF.m_fluidInside;
}

SphereSdf::~SphereSdf()
{}

Real
SphereSdf::value(const RealVect& a_point) const
{

  const RealVect newPos = a_point - m_center;

  Real retval = newPos.vectorLength() - m_radius; // Yields negative inside

  if (!m_fluidInside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
SphereSdf::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new SphereSdf(*this));
}

#include <CD_NamespaceFooter.H>
