/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MemoryReport.H
 * @brief  Headers for utility functions that report memory usage
 * @author Robert Marskar
 */

#ifndef CD_MEMORYREPORT_H
#define CD_MEMORYREPORT_H

// Chombo includes
#include <REAL.H>
#include <Vector.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief Namespace containing utility functions for reporting memory usage in various formats
 */
namespace MemoryReport {

/**
 * @brief Print the global maximum and minimum memory usage. This includes peak memory and unfreed memory.
 */
void
getMaxMinMemoryUsage();

/**
 * @brief Print the global maximum and minimum memory usage. This includes peak memory and unfreed memory.
 * @param[out] a_maxPeak    Maximum peak memory (maximized over ranks)
 * @param[out] a_minPeak    Minimum peak memory (minimuzed over ranks)
 * @param[out] a_maxUnfreed Maximum unfreed memory (maximized over ranks)
 * @param[out] a_minUnfreed Minimum unfreed memory (minimized over ranks)
 */
void
getMaxMinMemoryUsage(Real& a_maxPeak, Real& a_minPeak, Real& a_maxUnfreed, Real& a_minUnfreed);

/**
 * @brief Get peak and unfreed memory usage per rank.
 * @details If using MPI, this will do an Allgather so that every rank can see usage.
 * @param[out] a_peak    Peak memory usage
 * @param[in]  a_unfreed Unfreed memory usage
 */
void
getMemoryUsage(Vector<Real>& a_peak, Vector<Real>& a_unfreed);
/**
 * @brief Get the calling process' current and peak resident set size.
 * @details This is the memory the operating system has actually given the process, read from
 * /proc/self/status (VmRSS and VmHWM). Unlike every other figure in this namespace it counts the
 * whole process -- mesh arenas, particle arenas, every std::vector buffer, MPI internals and the
 * binary itself -- so it is the quantity a job is killed for exceeding, and the only one that can
 * tell real growth from a miscounting allocator. When a tracked figure climbs and this does not, the
 * growth is in the counter rather than in the process.
 *
 * Both outputs are zero on platforms without /proc, so a zero reading means "unavailable" rather
 * than "no memory in use".
 * @param[out] a_currentRSS Current resident set size (VmRSS), in bytes.
 * @param[out] a_peakRSS    Peak resident set size (VmHWM), in bytes.
 */
void
getResidentSetSize(Real& a_currentRSS, Real& a_peakRSS);

/**
 * @brief Get the global maximum and minimum resident set size over all ranks.
 * @details See getResidentSetSize(). The maximum is what decides an out-of-memory failure, because
 * the job dies when the worst rank exceeds its limit and not when the average does.
 * @param[out] a_maxCurrent Maximum current RSS, in bytes (maximized over ranks).
 * @param[out] a_minCurrent Minimum current RSS, in bytes (minimized over ranks).
 * @param[out] a_maxPeak    Maximum peak RSS, in bytes (maximized over ranks).
 * @param[out] a_minPeak    Minimum peak RSS, in bytes (minimized over ranks).
 */
void
getMaxMinResidentSetSize(Real& a_maxCurrent, Real& a_minCurrent, Real& a_maxPeak, Real& a_minPeak);
} // namespace MemoryReport

#include <CD_NamespaceFooter.H>

#endif
