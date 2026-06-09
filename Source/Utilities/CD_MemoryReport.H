/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MemoryReport.H
  @brief  Headers for utility functions that report memory usage
  @author Robert Marskar
*/

#ifndef CD_MEMORYREPORT_H
#define CD_MEMORYREPORT_H

// Chombo includes
#include <REAL.H>
#include <Vector.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Namespace containing utility functions for reporting memory usage in various formats
*/
namespace MemoryReport {

  /**
    @brief Print the global maximum and minimum memory usage. This includes peak memory and unfreed memory. 
  */
  void
  getMaxMinMemoryUsage();

  /**
    @brief Print the global maximum and minimum memory usage. This includes peak memory and unfreed memory. 
    @param[out] a_maxPeak    Maximum peak memory (maximized over ranks)
    @param[out] a_minPeak    Minimum peak memory (minimuzed over ranks)
    @param[out] a_maxUnfreed Maximum unfreed memory (maximized over ranks)
    @param[out] a_minUnfreed Minimum unfreed memory (minimized over ranks)
  */
  void
  getMaxMinMemoryUsage(Real& a_maxPeak, Real& a_minPeak, Real& a_maxUnfreed, Real& a_minUnfreed);

  /**
    @brief Get peak and unfreed memory usage per rank. 
    @details If using MPI, this will do an Allgather so that every rank can see usage. 
    @param[out] a_peak    Peak memory usage
    @param[in]  a_unfreed Unfreed memory usage
  */
  void
  getMemoryUsage(Vector<Real>& a_peak, Vector<Real>& a_unfreed);
} // namespace MemoryReport

#include <CD_NamespaceFooter.H>

#endif
