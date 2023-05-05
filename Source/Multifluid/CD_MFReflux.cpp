/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFReflux.cpp
  @brief  Implementation of CD_MFReflux.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFReflux.H>
#include <CD_NamespaceHeader.H>

MFReflux::MFReflux() {}

MFReflux::MFReflux(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs) { this->define(a_fluxRegs); }

MFReflux::~MFReflux() {}

void
MFReflux::define(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs)
{
  m_fluxRegs = a_fluxRegs;
}

const RefCountedPtr<EBReflux>&
MFReflux::getFluxRegPointer(const int a_phase) const
{
  return m_fluxRegs[a_phase];
}

EBReflux&
MFReflux::getFluxReg(const int a_phase)
{
  return *m_fluxRegs[a_phase];
}

const EBReflux&
MFReflux::getFluxReg(const int a_phase) const
{
  return *m_fluxRegs[a_phase];
}

#include <CD_NamespaceFooter.H>
