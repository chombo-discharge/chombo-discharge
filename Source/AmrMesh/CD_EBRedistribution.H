/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBRedistribution.H
 * @brief  Declaration of an abstract interface for the cut-cell stabilization of explicit finite-volume updates.
 * @author Robert Marskar
 */

#ifndef CD_EBREDISTRIBUTION_H
#define CD_EBREDISTRIBUTION_H

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Abstract interface for the cut-cell stabilization of explicit finite-volume updates.
 * @details An explicit update phi^{n+1} = phi^n - dt*D with the conservative divergence D = 1/(kappa*dx) * sum(fluxes)
 * is unstable in small cut cells at a time step set by the regular cells. An implementation of this class turns the
 * conservative divergence into the right-hand side R of a stabilized update,
 *
 *    phi^{n+1} = phi^n - dt * R,
 *
 * which is conservative. R is per unit volume, i.e. not weighted by the volume fraction.
 *
 * Depending on the implementation, R depends on phi^n and on dt. Callers must therefore pass the state that the update
 * starts from and the time step of the update, use R only with that time step, and not add together right-hand sides
 * that were computed separately. Several divergences that enter the same update are added together before they are
 * passed to redistribute().
 */
class EBRedistribution
{
public:
  /**
   * @brief Default constructor.
   */
  EBRedistribution() noexcept;

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other object
   */
  EBRedistribution(const EBRedistribution& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other object
   */
  EBRedistribution(const EBRedistribution&& a_other) = delete;

  /**
   * @brief Disallowed copy assignment.
   * @param[in] a_other Other object
   * @return This object
   */
  EBRedistribution&
  operator=(const EBRedistribution& a_other) = delete;

  /**
   * @brief Disallowed move assignment.
   * @param[in] a_other Other object
   * @return This object
   */
  EBRedistribution&
  operator=(const EBRedistribution&& a_other) = delete;

  /**
   * @brief Destructor.
   */
  virtual ~EBRedistribution() noexcept;

  /**
   * @brief Turn the conservative divergence into the right-hand side of the stabilized update.
   * @details On input, a_divG holds kappa*D, the conservative divergence multiplied by the volume fraction, with
   * coarse cells under finer levels averaged and ghost cells filled. On output, it holds R, such that the update is
   * phi^{n+1} = phi^n - a_dt * R.
   * @param[in,out] a_divG On input kappa*D, on output R.
   * @param[in]     a_phi  State that the update starts from.
   * @param[in]     a_dt   Time step of the update. It is signed: an update phi^n + dt*D passes -dt.
   */
  virtual void
  redistribute(EBAMRCellData& a_divG, const EBAMRCellData& a_phi, const Real a_dt) noexcept = 0;
};

#include <CD_NamespaceFooter.H>

#endif
