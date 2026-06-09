/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MFReflux.cpp
  @brief  Implementation of CD_MFReflux.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFReflux.H>
#include <CD_NamespaceHeader.H>

MFReflux::MFReflux()
{}

MFReflux::MFReflux(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs)
{
  this->define(a_fluxRegs);
}

MFReflux::~MFReflux()
{}

void
MFReflux::define(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs)
{
  m_fluxRegs = a_fluxRegs;
}

const RefCountedPtr<EBReflux>&
MFReflux::getFluxRegPointer(const int a_phase) const
{
  CH_assert(a_phase >= 0 && a_phase < m_fluxRegs.size());
  return m_fluxRegs[a_phase];
}

EBReflux&
MFReflux::getFluxReg(const int a_phase)
{
  CH_assert(a_phase >= 0 && a_phase < m_fluxRegs.size());
  CH_assert(!m_fluxRegs[a_phase].isNull());
  return *m_fluxRegs[a_phase];
}

const EBReflux&
MFReflux::getFluxReg(const int a_phase) const
{
  CH_assert(a_phase >= 0 && a_phase < m_fluxRegs.size());
  CH_assert(!m_fluxRegs[a_phase].isNull());
  return *m_fluxRegs[a_phase];
}

#include <CD_NamespaceFooter.H>
