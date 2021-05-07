/*!
  @file data_ops.H
  @brief Agglomeration of useful data operations
  @author Robert Marskar
  @date Nov 2017
*/

#include "memrep.H"

#include <EBLevelDataOps.H>
#include <memtrack.H>
#include <memusage.H>

namespace ChomboDischarge {

  void memrep::get_max_min_memory(){
    Real max_peak, min_peak, min_unfreed, max_unfreed;

    memrep::get_max_min_memory(max_peak, min_peak, min_unfreed, max_unfreed);

    pout() << "memrep::get_max_min_memory:" 
	   << "\t max peak = "       << 1.0*max_peak
	   << "\t min peak = "    << 1.0*min_peak
	   << "\t max unfreed = " << 1.0*max_unfreed
	   << "\t min unfreed = " << 1.0*min_unfreed << endl;
  }

  void memrep::get_max_min_memory(Real& a_max_peak, Real& a_min_peak, Real& a_min_unfreed, Real& a_max_unfreed){
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

  void memrep::get_memory(Vector<Real>& a_peak, Vector<Real>& a_unfreed){
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
}
