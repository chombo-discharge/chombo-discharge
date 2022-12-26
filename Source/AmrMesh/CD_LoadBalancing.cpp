/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file    CD_LoadBalancing.cpp
  @details Implementation of CD_LoadBalancing.H
  @author  Robert Marskar
*/

// Our includes
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

void
LoadBalancing::makeBalance(Vector<int>& a_ranks, const Vector<Box>& a_boxes)
{
  CH_TIME("LoadBalancing::makeBalance");

  LoadBalance(a_ranks, a_boxes);
}

void
LoadBalancing::roundRobin(Vector<int>& a_ranks, const Vector<Box>& a_boxes)
{
  CH_TIME("LoadBalancing::roundRobin");

  const int nProcs = numProc();
  const int nBoxes = a_boxes.size();

  a_ranks.resize(nBoxes);
  for (int ibox = 0; ibox < nBoxes; ibox++) {
    a_ranks[ibox] = ibox % nProcs;
  }
}

void
LoadBalancing::sort(Vector<Box>& a_boxes, const BoxSorting a_which)
{
  CH_TIME("LoadBalancing::sort");

  // The LoadBalancing::sort routines takes pairs of boxes/loads. Just use dummy loads here.
  if (a_which == BoxSorting::Morton) {
    mortonOrdering(a_boxes);
  }
  else { // Call chombo-discharge code
    Vector<int> dummy(a_boxes.size(), 0);

    LoadBalancing::sort(a_boxes, dummy, a_which);
  }
}

void
LoadBalancing::gatherBoxes(Vector<Box>& a_boxes)
{
  CH_TIME("LoadBalancing::gatherBoxes");

#ifdef CH_MPI
  // TLDR: This code does a gather operation on the loads and boxes. They are gather globally in this way:
  //       (BoxesForMPIRank=0, BoxesForMPIRank=1, BoxesForMPIRank=2, ....)

  // 1. Linearize local boxes
  int  send_size   = 2 * CH_SPACEDIM;            // Message size for one box
  int  send_count  = a_boxes.size() * send_size; // Number of elements sent from this rank
  int* send_buffer = new int[send_count];        // Send buffer for this MPI rank
  int* send_buf2   = send_buffer; // Backup address. Going to monkey with pointer increments on send buffer

  // Linearize a_boxes onto send_buffer
  for (int i = 0; i < a_boxes.size(); i++, send_buffer += send_size) {
    const Box& b = a_boxes[i];
    D_TERM6(send_buffer[0] = b.smallEnd(0); send_buffer[1] = b.bigEnd(0);, send_buffer[2] = b.smallEnd(1);
            send_buffer[3] = b.bigEnd(1);
            , send_buffer[4] = b.smallEnd(2);
            send_buffer[5] = b.bigEnd(2);
            , send_buffer[6] = b.smallEnd(3);
            send_buffer[7] = b.bigEnd(3);
            , send_buffer[8] = b.smallEnd(4);
            send_buffer[9] = b.bigEnd(4);
            , send_buffer[10] = b.smallEnd(5);
            send_buffer[11] = b.bigEnd(5););
  }
  send_buffer = send_buf2; // Revert point to start of array

  // 2. Get the number of elements sent from each rank
  int* send_counts = new int[numProc()];
  MPI_Allgather(&send_count, 1, MPI_INT, send_counts, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0]   = 0;
  for (int i = 0; i < numProc() - 1; i++) {
    offsets[i + 1] = offsets[i] + send_counts[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++) {
    total_count += send_counts[i];
  }
  int* recv_buffer = new int[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, send_count, MPI_INT, recv_buffer, send_counts, offsets, MPI_INT, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_boxes.resize(0);
  int* recv_buf2 = recv_buffer; // Going to monkey with pointer increments again
  for (int i = 0; i < total_count / send_size; i++, recv_buffer += send_size) {
    IntVect lo, hi;
    D_TERM6(lo[0] = recv_buffer[0]; hi[0] = recv_buffer[1];, lo[1] = recv_buffer[2]; hi[1] = recv_buffer[3];
            , lo[2]                                                                        = recv_buffer[4];
            hi[2] = recv_buffer[5];
            , lo[3] = recv_buffer[6];
            hi[3] = recv_buffer[7];
            , lo[4] = recv_buffer[8];
            hi[4] = recv_buffer[9];
            , lo[5] = recv_buffer[10];
            hi[5] = recv_buffer[11];);

    a_boxes.push_back(Box(lo, hi));
  }
  recv_buffer = recv_buf2;

  delete[] send_buffer;
  delete[] send_counts;
  delete[] offsets;
  delete[] recv_buffer;
#endif
}

void
LoadBalancing::gatherLoads(Vector<Real>& a_loads)
{
  CH_TIME("LoadBalancing::gatherLoads");

#ifdef CH_MPI
  // TLDR: This code does a gather operation on the loads. They are gather globally in this way:
  //       (LoadsForRank=0, LoadsForRank=1, LoadsForRank2=2, ....)

  // 1. Linearize local boxes onto appropriate memory structure
  int   send_size   = 1;                          // Message size for one element. Just one Real.
  int   mySendCount = a_loads.size() * send_size; // Number of elements sent from this rank
  Real* send_buffer = new Real[mySendCount];      // Send buffer for this MPI rank
  for (int i = 0; i < a_loads.size(); i++) {
    send_buffer[i] = a_loads[i];
  }

  // 2. Get the number of elements sent from each rank
  int* allSendCount = new int[numProc()];
  MPI_Allgather(&mySendCount, 1, MPI_INT, allSendCount, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0]   = 0;
  for (int i = 0; i < numProc() - 1; i++) {
    offsets[i + 1] = offsets[i] + allSendCount[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++) {
    total_count += allSendCount[i];
  }
  Real* recv_buffer = new Real[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer,
                 mySendCount,
                 MPI_CH_REAL,
                 recv_buffer,
                 allSendCount,
                 offsets,
                 MPI_CH_REAL,
                 Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_loads.resize(total_count);
  for (int i = 0; i < total_count; i++) {
    a_loads[i] = recv_buffer[i];
  }

  delete[] recv_buffer;
  delete[] offsets;
  delete[] allSendCount;
  delete[] send_buffer;
#endif
}

void
LoadBalancing::gatherLoads(Vector<long>& a_loads)
{
  CH_TIME("LoadBalancing::gatherLoads");

#ifdef CH_MPI
  // TLDR: This code does a gather operation on the loads. They are gather globally in this way:
  //       (LoadsForRank=0, LoadsForRank=1, LoadsForRank2=2, ....)

  // 1. Linearize local boxes onto appropriate memory structure
  int   send_size   = 1;                          // Message size for one element. Just one long.
  int   mySendCount = a_loads.size() * send_size; // Number of elements sent from this rank
  long* send_buffer = new long[mySendCount];      // Send buffer for this MPI rank
  for (int i = 0; i < a_loads.size(); i++) {
    send_buffer[i] = a_loads[i];
  }

  // 2. Get the number of elements sent from each rank
  int* allSendCount = new int[numProc()];
  MPI_Allgather(&mySendCount, 1, MPI_INT, allSendCount, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0]   = 0;
  for (int i = 0; i < numProc() - 1; i++) {
    offsets[i + 1] = offsets[i] + allSendCount[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++) {
    total_count += allSendCount[i];
  }
  long* recv_buffer = new long[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, mySendCount, MPI_LONG, recv_buffer, allSendCount, offsets, MPI_LONG, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_loads.resize(total_count);
  for (int i = 0; i < total_count; i++) {
    a_loads[i] = recv_buffer[i];
  }

  delete[] offsets;
  delete[] allSendCount;
  delete[] recv_buffer;
  delete[] send_buffer;
#endif
}

void
LoadBalancing::gatherLoads(Vector<int>& a_loads)
{
  CH_TIME("LoadBalancing::gatherLoads");

#ifdef CH_MPI
  // TLDR: This code does a gather operation on the loads. They are gather globally in this way:
  //       (LoadsForRank=0, LoadsForRank=1, LoadsForRank2=2, ....)

  // 1. Linearize local boxes onto appropriate memory structure
  int  send_size   = 1;                          // Message size for one element. Just one int.
  int  mySendCount = a_loads.size() * send_size; // Number of elements sent from this rank
  int* send_buffer = new int[mySendCount];       // Send buffer for this MPI rank
  for (int i = 0; i < a_loads.size(); i++) {
    send_buffer[i] = a_loads[i];
  }

  // 2. Get the number of elements sent from each rank
  int* allSendCount = new int[numProc()];
  MPI_Allgather(&mySendCount, 1, MPI_INT, allSendCount, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0]   = 0;
  for (int i = 0; i < numProc() - 1; i++) {
    offsets[i + 1] = offsets[i] + allSendCount[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++) {
    total_count += allSendCount[i];
  }
  int* recv_buffer = new int[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, mySendCount, MPI_INT, recv_buffer, allSendCount, offsets, MPI_INT, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_loads.resize(total_count);
  for (int i = 0; i < total_count; i++) {
    a_loads[i] = recv_buffer[i];
  }

  delete[] recv_buffer;
  delete[] send_buffer;
  delete[] offsets;
  delete[] allSendCount;
#endif
}

void
LoadBalancing::gatherBoxesAndLoads(Vector<Box>& a_boxes, Vector<int>& a_loads)
{
  CH_TIME("LoadBalancing::gatherBoxesAndLoads");
#ifdef CH_MPI
  LoadBalancing::gatherBoxes(a_boxes);
  LoadBalancing::gatherLoads(a_loads);
#endif
}

int
LoadBalancing::maxBits(std::vector<Box>::iterator a_first, std::vector<Box>::iterator a_last)
{
  CH_TIME("LoadBalancing::maxBits");

  int maxSize = 0;
  for (std::vector<Box>::iterator p = a_first; p < a_last; ++p) {
    IntVect small = p->smallEnd();
    D_EXPR6(maxSize = Max(maxSize, std::abs(small[0])),
            maxSize = Max(maxSize, std::abs(small[1])),
            maxSize = Max(maxSize, std::abs(small[2])),
            maxSize = Max(maxSize, std::abs(small[3])),
            maxSize = Max(maxSize, std::abs(small[4])),
            maxSize = Max(maxSize, std::abs(small[5])));
  }
  int bits;
  for (bits = 8 * sizeof(int) - 2; bits > 0; bits--) {
    const int N = (1 << bits);

    if (maxSize / N > 0) {
      break;
    }
  }
  bits++;
  return bits;
}

#include <CD_NamespaceFooter.H>
