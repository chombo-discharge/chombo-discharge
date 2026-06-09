/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MFCoarAve.H
  @brief  Wrapper class for holding multifluid EBCoarAves
  @author Robert Marskar
*/

#ifndef CD_MFCOARAVE_H
#define CD_MFCOARAVE_H

// Std includes
#include <map>

// Our includes
#include <CD_EBCoarAve.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class for coarsening data in a multifluid context.
*/
class MFCoarAve
{
public:
  /**
    @brief Constructor. Must subsequently call define.
  */
  MFCoarAve();

  /**
    @brief Destructor (does nothing).
  */
  ~MFCoarAve();

  /**
    @brief Full constructor. Calls define
    @param[in] a_aveOps Averaging operators on each phase.
  */
  MFCoarAve(const Vector<RefCountedPtr<EBCoarAve>>& a_aveOps);

  /**
    @brief Define function
    @param[in] a_aveOps Averaging operators on each phase.
  */
  void
  define(const Vector<RefCountedPtr<EBCoarAve>>& a_aveOps);

  /**
    @brief Get coarsening utility for specified phase.
    @param[in] a_phase Phase
    @return Const reference to the averaging operator for the given phase.
  */
  const RefCountedPtr<EBCoarAve>&
  getAveOp(const int a_phase) const;

protected:
  /**
    @brief Averaging operators for each phase.
  */
  Vector<RefCountedPtr<EBCoarAve>> m_aveOps;
};

#include <CD_NamespaceFooter.H>

#endif
