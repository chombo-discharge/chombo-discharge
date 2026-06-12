/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MFReflux.H
  @brief  Declaration of a class for refluxing in a multiphase context.
  @author Robert Marskar
*/

#ifndef CD_MFREFLUX_H
#define CD_MFREFLUX_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_EBReflux.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class which wraps EBRefluxs in multiphase.
*/
class MFReflux
{
public:
  /**
    @brief Weak constructor. Must subsequently call define.
  */
  MFReflux();

  /**
    @brief Full constructor. Calls define
    @param[in] a_fluxRegs Flux registers on each phase
  */
  MFReflux(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs);

  /**
    @brief Destructor (does nothing)
  */
  ~MFReflux();

  /**
    @brief Define function
    @param[in] a_fluxRegs Flux registers on each phase
  */
  void
  define(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs);

  /**
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
    @return Const reference-counted pointer to flux register for the given phase.
  */
  const RefCountedPtr<EBReflux>&
  getFluxRegPointer(const int a_phase) const;

  /**
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
    @return Reference to flux register for the given phase.
  */
  EBReflux&
  getFluxReg(const int a_phase);

  /**
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
    @return Const reference to flux register for the given phase.
  */
  const EBReflux&
  getFluxReg(const int a_phase) const;

protected:
  /**
    @brief Flux register for each phase.
  */
  Vector<RefCountedPtr<EBReflux>> m_fluxRegs;
};

#include <CD_NamespaceFooter.H>

#endif
