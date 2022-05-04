/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFFluxReg.cpp
  @brief  Implementation of CD_MFFluxReg.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFFluxReg.H>
#include <CD_NamespaceHeader.H>

MFFluxReg::MFFluxReg() {}

MFFluxReg::MFFluxReg(const Vector<RefCountedPtr<EBFluxRegister>>& a_fluxRegs) { this->define(a_fluxRegs); }

MFFluxReg::~MFFluxReg() {}

void
MFFluxReg::define(const Vector<RefCountedPtr<EBFluxRegister>>& a_fluxRegs)
{
  m_fluxRegs = a_fluxRegs;
}

const RefCountedPtr<EBFluxRegister>&
MFFluxReg::getFluxRegPointer(const int a_phase) const
{
  return m_fluxRegs[a_phase];
}

EBFluxRegister&
MFFluxReg::getFluxReg(const int a_phase)
{
  return *m_fluxRegs[a_phase];
}

const EBFluxRegister&
MFFluxReg::getFluxReg(const int a_phase) const
{
  return *m_fluxRegs[a_phase];
}

#include <CD_NamespaceFooter.H>
