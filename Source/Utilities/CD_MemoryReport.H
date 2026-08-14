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

// Std includes
#include <iosfwd>
#include <string>

// Chombo includes
#include <REAL.H>
#include <Vector.H>

// Our includes
#include <CD_ParticleMemory.H>
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

/**
 * @brief Format a byte count as megabytes for a report line.
 * @details Fixed to three decimals, i.e. a resolution of about a kilobyte. Megabytes are what a
 * reader wants for a figure in the gigabyte range, but rounding to whole megabytes destroys the
 * signal these reports exist to carry: a ratchet of tens of kilobytes per step reads as an
 * unchanging integer, and a counter that is quietly climbing looks flat. Formatted into its own
 * stream so the caller's formatting flags are left alone.
 * @param[in] a_bytes Byte count.
 * @return The value in MB, as a fixed-point string with three decimals.
 */
std::string
bytesAsMB(const long long a_bytes);

/**
 * @brief Format a byte count as megabytes for a report line.
 * @details Overload for the byte counts this namespace reports as Real. See bytesAsMB(long long).
 * @param[in] a_bytes Byte count.
 * @return The value in MB, as a fixed-point string with three decimals.
 */
std::string
bytesAsMB(const Real a_bytes);

/**
 * @brief Write the full memory report as a table.
 * @details Three sources that do not overlap, gathered here so that they cannot drift apart: the
 * resident set is the whole process, Chombo's counter sees only its own Arena (BaseFab-backed mesh
 * data), and ParticleMemory sees only particle storage. Neither counter can adjudicate itself -- a
 * tracked figure that climbs while the resident set does not is a bug in the counter rather than
 * memory the process is holding -- which is why all three are printed together.
 *
 * The table has two independent axes, and conflating them has cost real time. A "high-water" row is
 * over *time* and never falls; the columns are over *ranks* at this instant. The figure a job is
 * killed for exceeding is the high-water of the worst rank, i.e. the rightmost value of a high-water
 * row. A dash marks a quantity that does not exist rather than one that is zero: there is no useful
 * minimum over ranks of a high-water mark, because the ranks reach theirs at different moments.
 *
 * Prints nothing unless Chombo was compiled with memory tracking. The arena figures are the reason
 * the table exists, and a table reporting two of its three sources would invite exactly the confusion
 * between counter and process that the third column is there to prevent. The guard lives here rather
 * than at the call sites, so callers need no preprocessor conditionals of their own.
 *
 * Collective: the rank columns are reductions, so every rank must call this.
 * @param[in,out] a_out    Stream to print to.
 * @param[in]     a_indent Prefix written before every line, so the table can sit inside any report.
 */
void
printMemoryTable(std::ostream& a_out, const std::string& a_indent);
} // namespace MemoryReport

#include <CD_NamespaceFooter.H>

#endif
