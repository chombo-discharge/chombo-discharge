/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MemoryReport.cpp
 * @brief  Implementation of CD_MemoryReport.H
 * @author Robert Marskar
 */

// Std includes
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

// Chombo includes
#include <memtrack.H>
#include <parstream.H>
#include <memusage.H>
#include <SPMD.H>
#include <CH_Timer.H>

// Our includes
#include <CD_MemoryReport.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

void
MemoryReport::getMaxMinMemoryUsage()
{
  CH_TIME("MemoryReport::getMaxMinMemoryUsage");
#ifdef CH_USE_MEMORY_TRACKING
  Real maxPeak;
  Real minPeak;
  Real maxUnfreed;
  Real minUnfreed;

  MemoryReport::getMaxMinMemoryUsage(maxPeak, minPeak, maxUnfreed, minUnfreed);

  pout() << "MemoryReport::getMaxMinMemoryUsage:" << "\t Max peak = " << 1.0 * maxPeak
         << "\t Min peak = " << 1.0 * minPeak << "\t Max unfreed = " << 1.0 * maxUnfreed
         << "\t Min unfreed = " << 1.0 * minUnfreed << endl;
#endif
}

void
MemoryReport::getMaxMinMemoryUsage(Real& a_maxPeak, Real& a_minPeak, Real& a_maxUnfreed, Real& a_minUnfreed)
{
  CH_TIME("MemoryReport::getMaxMinMemoryUsage");

  constexpr int BytesPerMB = 1024 * 1024;

  // Gets usage in bytes.
  long long curMemLL  = 0LL;
  long long peakMemLL = 0LL;
#ifdef CH_USE_MEMORY_TRACKING
  overallMemoryUsage(curMemLL, peakMemLL);
#endif

  const int unfreedMem = int(curMemLL);
  const int peakMem    = int(peakMemLL);

  int maxPeakMem    = 0;
  int minPeakMem    = 0;
  int maxUnfreedMem = 0;
  int minUnfreedMem = 0;

  // Find maximum/minimum usage.
#ifdef CH_MPI
  MPI_Allreduce(&peakMem, &maxPeakMem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  MPI_Allreduce(&peakMem, &minPeakMem, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
  MPI_Allreduce(&unfreedMem, &maxUnfreedMem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  MPI_Allreduce(&unfreedMem, &minUnfreedMem, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
#else
  maxPeakMem    = peakMem;
  minPeakMem    = peakMem;
  maxUnfreedMem = unfreedMem;
  minUnfreedMem = unfreedMem;
#endif

  // Convert to real.
  a_maxPeak    = 1.0 * maxPeakMem / BytesPerMB;
  a_minPeak    = 1.0 * minPeakMem / BytesPerMB;
  a_maxUnfreed = 1.0 * maxUnfreedMem / BytesPerMB;
  a_minUnfreed = 1.0 * minUnfreedMem / BytesPerMB;
}

void
MemoryReport::getMemoryUsage(Vector<Real>& a_peak, Vector<Real>& a_unfreed)
{
  CH_TIME("MemoryReport::getMemoryUsage");

  constexpr int BytesPerMB = 1024 * 1024;

  long long curMemLL  = 0LL;
  long long peakMemLL = 0LL;
#ifdef CH_USE_MEMORY_TRACKING
  overallMemoryUsage(curMemLL, peakMemLL);
#endif

  const int unfreedMem = int(curMemLL);
  const int peakMem    = int(peakMemLL);

#ifdef CH_MPI
  int* unfreed = (int*)malloc(numProc() * sizeof(int)); // new int[numProc()];
  int* peak    = (int*)malloc(numProc() * sizeof(int)); // new int[numProc()];

  MPI_Allgather(&peakMem, 1, MPI_INT, peak, 1, MPI_INT, Chombo_MPI::comm);
  MPI_Allgather(&unfreedMem, 1, MPI_INT, unfreed, 1, MPI_INT, Chombo_MPI::comm);

  a_peak.resize(numProc());
  a_unfreed.resize(numProc());

  for (int i = 0; i < numProc(); i++) {
    a_peak[i]    = 1.0 * peak[i] / BytesPerMB;
    a_unfreed[i] = 1.0 * unfreed[i] / BytesPerMB;
  }

  free(unfreed);
  free(peak);
#else
  a_peak.resize(1);
  a_unfreed.resize(1);

  a_peak[0]    = peakMem / BytesPerMB;
  a_unfreed[0] = unfreedMem / BytesPerMB;
#endif
}

void
MemoryReport::getResidentSetSize(Real& a_currentRSS, Real& a_peakRSS)
{
  CH_TIME("MemoryReport::getResidentSetSize");

  a_currentRSS = 0.0;
  a_peakRSS    = 0.0;

  // /proc/self/status reports these as 'VmRSS:<whitespace><number> kB'. Absent on non-Linux
  // platforms, in which case both outputs are left at zero.
  std::ifstream status("/proc/self/status");

  if (status.is_open()) {
    std::string line;

    while (std::getline(status, line)) {
      const bool isCurrent = line.compare(0, 6, "VmRSS:") == 0;
      const bool isPeak    = line.compare(0, 6, "VmHWM:") == 0;

      if (isCurrent || isPeak) {
        std::istringstream parser(line.substr(6));

        double kiloBytes = 0.0;
        parser >> kiloBytes;

        if (isCurrent) {
          a_currentRSS = 1024.0 * kiloBytes;
        }
        else {
          a_peakRSS = 1024.0 * kiloBytes;
        }
      }
    }
  }
}

void
MemoryReport::getMaxMinResidentSetSize(Real& a_maxCurrent, Real& a_minCurrent, Real& a_maxPeak, Real& a_minPeak)
{
  CH_TIME("MemoryReport::getMaxMinResidentSetSize");

  Real currentRSS = 0.0;
  Real peakRSS    = 0.0;

  MemoryReport::getResidentSetSize(currentRSS, peakRSS);

  a_maxCurrent = ParallelOps::max(currentRSS);
  a_minCurrent = ParallelOps::min(currentRSS);
  a_maxPeak    = ParallelOps::max(peakRSS);
  a_minPeak    = ParallelOps::min(peakRSS);
}

std::string
MemoryReport::bytesAsMB(const long long a_bytes)
{
  constexpr double bytesPerMB = 1024.0 * 1024.0;

  std::ostringstream os;

  os << std::fixed << std::setprecision(3) << (static_cast<double>(a_bytes) / bytesPerMB);

  return os.str();
}

std::string
MemoryReport::bytesAsMB(const Real a_bytes)
{
  return MemoryReport::bytesAsMB(static_cast<long long>(a_bytes));
}

void
MemoryReport::printMemoryTable([[maybe_unused]] std::ostream& a_out, [[maybe_unused]] const std::string& a_indent)
{
  CH_TIME("MemoryReport::printMemoryTable");

#ifdef CH_USE_MEMORY_TRACKING
  constexpr int labelWidth  = 24;
  constexpr int columnWidth = 10;

  const std::string none = "-";

  auto rule = [&a_indent, &a_out]() -> void {
    a_out << a_indent << std::string(labelWidth + 3 * columnWidth, '-') << endl;
  };

  auto row = [&a_indent, &a_out](const std::string& a_label,
                                 const std::string& a_local,
                                 const std::string& a_max,
                                 const std::string& a_min) -> void {
    a_out << a_indent << std::left << std::setw(labelWidth) << a_label << std::right << std::setw(columnWidth)
          << a_local << std::setw(columnWidth) << a_max << std::setw(columnWidth) << a_min << endl;
  };

  rule();
  row("Memory (MB)", "this rank", "rank max", "rank min");
  rule();

  Real currentRSS = 0.0;
  Real peakRSS    = 0.0;

  getResidentSetSize(currentRSS, peakRSS);

  if (currentRSS > 0.0) {
    Real maxCurrent = currentRSS;
    Real minCurrent = currentRSS;
    Real maxPeak    = peakRSS;
    Real minPeak    = peakRSS;

#ifdef CH_MPI
    getMaxMinResidentSetSize(maxCurrent, minCurrent, maxPeak, minPeak);
#endif

    row("Resident set", bytesAsMB(currentRSS), bytesAsMB(maxCurrent), bytesAsMB(minCurrent));
    row("  high-water", bytesAsMB(peakRSS), bytesAsMB(maxPeak), none);
    rule();
  }

  long long meshBytes = 0LL;
  long long meshPeak  = 0LL;

  overallMemoryUsage(meshBytes, meshPeak);

#ifdef CH_MPI
  row("Mesh arenas",
      bytesAsMB(meshBytes),
      bytesAsMB(ParallelOps::max(meshBytes)),
      bytesAsMB(ParallelOps::min(meshBytes)));
  row("  high-water", bytesAsMB(meshPeak), bytesAsMB(ParallelOps::max(meshPeak)), none);
#else
  row("Mesh arenas", bytesAsMB(meshBytes), bytesAsMB(meshBytes), bytesAsMB(meshBytes));
  row("  high-water", bytesAsMB(meshPeak), bytesAsMB(meshPeak), none);
#endif
  rule();

  const long long arenaBytes  = static_cast<long long>(ParticleMemory::getBytes(ParticleMemory::Kind::Container));
  const long long arenaPeak   = static_cast<long long>(ParticleMemory::getPeak(ParticleMemory::Kind::Container));
  const long long bufferBytes = static_cast<long long>(ParticleMemory::getBytes(ParticleMemory::Kind::Buffer));
  const long long bufferPeak  = static_cast<long long>(ParticleMemory::getPeak(ParticleMemory::Kind::Buffer));
  const long long totalBytes  = static_cast<long long>(ParticleMemory::getAllocatedBytes());
  const long long totalPeak   = static_cast<long long>(ParticleMemory::getPeakBytes());

#ifdef CH_MPI
  row("Particle arenas",
      bytesAsMB(arenaBytes),
      bytesAsMB(ParallelOps::max(arenaBytes)),
      bytesAsMB(ParallelOps::min(arenaBytes)));
  row("  high-water", bytesAsMB(arenaPeak), bytesAsMB(ParallelOps::max(arenaPeak)), none);
  row("Particle buffers",
      bytesAsMB(bufferBytes),
      bytesAsMB(ParallelOps::max(bufferBytes)),
      bytesAsMB(ParallelOps::min(bufferBytes)));
  row("  high-water", bytesAsMB(bufferPeak), bytesAsMB(ParallelOps::max(bufferPeak)), none);
  row("Particle total",
      bytesAsMB(totalBytes),
      bytesAsMB(ParallelOps::max(totalBytes)),
      bytesAsMB(ParallelOps::min(totalBytes)));
  row("  high-water", bytesAsMB(totalPeak), bytesAsMB(ParallelOps::max(totalPeak)), none);
#else
  row("Particle arenas", bytesAsMB(arenaBytes), bytesAsMB(arenaBytes), bytesAsMB(arenaBytes));
  row("  high-water", bytesAsMB(arenaPeak), bytesAsMB(arenaPeak), none);
  row("Particle buffers", bytesAsMB(bufferBytes), bytesAsMB(bufferBytes), bytesAsMB(bufferBytes));
  row("  high-water", bytesAsMB(bufferPeak), bytesAsMB(bufferPeak), none);
  row("Particle total", bytesAsMB(totalBytes), bytesAsMB(totalBytes), bytesAsMB(totalBytes));
  row("  high-water", bytesAsMB(totalPeak), bytesAsMB(totalPeak), none);
#endif
  rule();
#endif
}

#include <CD_NamespaceFooter.H>
