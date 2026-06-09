/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MFMultigridInterpolator.cpp
  @brief  Implementation of CD_MFMultigridInterpolator.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

MFMultigridInterpolator::MFMultigridInterpolator()
{}

MFMultigridInterpolator::MFMultigridInterpolator(const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators)
{
  this->define(a_interpolators);
}

MFMultigridInterpolator::~MFMultigridInterpolator()
{}

MFMultigridInterpolator&
MFMultigridInterpolator::operator=(const MFMultigridInterpolator& a_other)
{
  m_interpolators = a_other.m_interpolators;

  return *this;
}

void
MFMultigridInterpolator::define(const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators)
{
  m_interpolators = a_interpolators;
}

RefCountedPtr<EBMultigridInterpolator>&
MFMultigridInterpolator::getInterpolator(const int a_phase) const
{
  CH_assert(a_phase >= 0 && a_phase < m_interpolators.size());
  return m_interpolators[a_phase];
}

int
MFMultigridInterpolator::getGhostCF() const
{
  return m_interpolators.front()->getGhostCF();
}

#include <CD_NamespaceFooter.H>
