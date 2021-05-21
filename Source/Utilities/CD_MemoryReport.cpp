/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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

// Our includes
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

void MemoryReport::getMaxMinMemoryUsage(){
  Real max_peak, min_peak, min_unfreed, max_unfreed;

  MemoryReport::getMaxMinMemoryUsage(max_peak, min_peak, min_unfreed, max_unfreed);

  pout() << "MemoryReport::getMaxMinMemoryUsage:" 
	 << "\t max peak = "       << 1.0*max_peak
	 << "\t min peak = "    << 1.0*min_peak
	 << "\t max unfreed = " << 1.0*max_unfreed
	 << "\t min unfreed = " << 1.0*min_unfreed << endl;
}

void MemoryReport::getMaxMinMemoryUsage(Real& a_max_peak, Real& a_min_peak, Real& a_min_unfreed, Real& a_max_unfreed){
  int BytesPerMB = 1024*1024;
  long long curMem;
  long long peakMem;
  overallMemoryUsage(curMem, peakMem);

  int unfreed_mem = curMem;
  int peak_mem    = peakMem;

  int max_unfreed_mem;
  int max_peak_mem;
  int min_unfreed_mem;
  int min_peak_mem;

  int result1 = MPI_Allreduce(&unfreed_mem, &max_unfreed_mem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  int result2 = MPI_Allreduce(&peak_mem,    &max_peak_mem,    1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  int result3 = MPI_Allreduce(&unfreed_mem, &min_unfreed_mem, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
  int result4 = MPI_Allreduce(&peak_mem,    &min_peak_mem,    1, MPI_INT, MPI_MIN, Chombo_MPI::comm);

  a_max_peak    = 1.0*max_peak_mem/BytesPerMB;
  a_min_peak    = 1.0*min_peak_mem/BytesPerMB;
  a_max_unfreed = 1.0*max_unfreed_mem/BytesPerMB;
  a_min_unfreed = 1.0*min_unfreed_mem/BytesPerMB;
}

void MemoryReport::getMemoryUsage(Vector<Real>& a_peak, Vector<Real>& a_unfreed){
  const int BytesPerMB = 1024*1024;

  long long curMem, peakMem;
  overallMemoryUsage(curMem, peakMem);

  int unfreed_mem = curMem;
  int peak_mem    = peakMem;


  int* unfreed = (int*) malloc(numProc()*sizeof(int));//new int[numProc()];
  int* peak    = (int*) malloc(numProc()*sizeof(int));//new int[numProc()];

  MPI_Allgather(&peak_mem,    1, MPI_INT, peak,    1, MPI_INT, Chombo_MPI::comm);
  MPI_Allgather(&unfreed_mem, 1, MPI_INT, unfreed, 1, MPI_INT, Chombo_MPI::comm);

  a_peak.resize(numProc());
  a_unfreed.resize(numProc());
  for (int i = 0; i < numProc(); i++){
    a_peak[i]    = 1.0*peak[i]/BytesPerMB;
    a_unfreed[i] = 1.0*unfreed[i]/BytesPerMB;
  }

  delete unfreed;
  delete peak;
}

#include <CD_NamespaceFooter.H>
