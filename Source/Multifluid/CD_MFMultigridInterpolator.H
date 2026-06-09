/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MFMultigridInterpolator.H
  @brief  Declaration of a wrapper class for holding multifluid EBMultigridInterpolators
  @author Robert Marskar
*/

#ifndef CD_MFMULTIGRIDINTERPOLATOR_H
#define CD_MFMULTIGRIDINTERPOLATOR_H

// Std includes
#include <map>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Wrapper class for holding multifluid EBMultigridInterpolators.
*/
class MFMultigridInterpolator
{
public:
  /**
    @brief Weak constructor. Must subsequently call define.
  */
  MFMultigridInterpolator();

  /**
    @brief Defining constructor.
    @param[in] a_interpolators Interpolators on each phase.
  */
  MFMultigridInterpolator(const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators);

  /**
    @brief Destructor. Does nothing.
  */
  ~MFMultigridInterpolator();

  /**
    @brief Copy assignment
    @param[in] a_other Other interpolator
    @return Reference to this object.
  */
  MFMultigridInterpolator&
  operator=(const MFMultigridInterpolator& a_other);

  /**
    @brief Define function
    @param[in] a_interpolators Interpolators on each phase.
  */
  void
  define(const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators);

  /**
    @brief Get interpolator for specified phase
    @param[in] a_phase Phase index.
    @return Reference to the interpolator for the given phase.
  */
  RefCountedPtr<EBMultigridInterpolator>&
  getInterpolator(const int a_phase) const;

  /**
    @brief Get number of filled ghost cells
    @return Number of ghost cells filled across coarse-fine interface.
  */
  int
  getGhostCF() const;

protected:
  /**
    @brief Interpolators on each phase.
  */
  mutable Vector<RefCountedPtr<EBMultigridInterpolator>> m_interpolators;
};

#include <CD_NamespaceFooter.H>

#endif
