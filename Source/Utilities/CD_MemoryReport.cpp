/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MemoryReport.cpp
  @brief  Implementation of CD_MemoryReport.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBLevelDataOps.H>
#include <memtrack.H>
#include <memusage.H>
#include <CH_Timer.H>

// Our includes
#include <CD_MemoryReport.H>
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

  pout() << "MemoryReport::getMaxMinMemoryUsage:"
         << "\t Max peak = " << 1.0 * maxPeak << "\t Min peak = " << 1.0 * minPeak
         << "\t Max unfreed = " << 1.0 * maxUnfreed << "\t Min unfreed = " << 1.0 * minUnfreed << endl;
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
  int* unfreed = (int*)malloc(numProc() * sizeof(int)); //new int[numProc()];
  int* peak    = (int*)malloc(numProc() * sizeof(int)); //new int[numProc()];

  MPI_Allgather(&peakMem, 1, MPI_INT, peak, 1, MPI_INT, Chombo_MPI::comm);
  MPI_Allgather(&unfreedMem, 1, MPI_INT, unfreed, 1, MPI_INT, Chombo_MPI::comm);

  a_peak.resize(numProc());
  a_unfreed.resize(numProc());

  for (int i = 0; i < numProc(); i++) {
    a_peak[i]    = 1.0 * peak[i] / BytesPerMB;
    a_unfreed[i] = 1.0 * unfreed[i] / BytesPerMB;
  }

  delete unfreed;
  delete peak;
#else
  a_peak.resize(1);
  a_unfreed.resize(1);

  a_peak[0]    = peakMem / BytesPerMB;
  a_unfreed[0] = unfreedMem / BytesPerMB;
#endif
}

#include <CD_NamespaceFooter.H>
