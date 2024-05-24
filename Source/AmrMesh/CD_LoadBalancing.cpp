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
LoadBalancing::makeBalance2(Vector<int>&        a_ranks,
                            Loads&              a_rankLoads,
                            const Vector<Real>& a_boxLoads,
                            const Vector<Box>&  a_boxes)
{
  CH_TIME("LoadBalancing::makeBalance");

  // Minimum number of grid subsets is the number of boxe, and the maximum number of grid subsets
  // is the number of ranks.
  const int numBoxes   = a_boxes.size();
  const int numRanks   = numProc();
  const int numSubsets = std::min(numBoxes, numRanks);

  if (numSubsets > 0) {

    // Figure out the total and target load (load per subset) on this level.
    Real totalLoad        = 0.0;
    Real staticTargetLoad = 0.0;
    for (int ibox = 0; ibox < numBoxes; ibox++) {
      totalLoad += a_boxLoads[ibox];
    }

    staticTargetLoad = totalLoad / numSubsets;

    // Build the grid subsets. When we do this we iterate through the boxes and try to ensure that we partition
    // the subsets such that the subsetLoad is as close to the dynamic targetLoad as possible.
    //
    // The pair contains the starting index for the subset and the computational load for the subset.
    using Span   = std::pair<int, int>;
    using Subset = std::pair<Span, Real>;

    std::vector<Subset> subsets(numSubsets);

    pout() << endl;
    pout() << "//// start subset report" << endl;
    pout() << "num boxes = " << numBoxes << endl;
    pout() << "num subsets = " << numSubsets << endl;
    pout() << "target load = " << totalLoad / numSubsets << endl;
    pout() << "total load = " << totalLoad << endl;

    int firstSubsetBox = 0;

    Real remainingLoad = totalLoad;

    for (int curSubset = 0; curSubset < numSubsets; curSubset++) {

      // The firstSubsetBox is the index for the first box in this subset (always assigned).
      Real subsetLoad = a_boxLoads[firstSubsetBox];

      int lastSubsetBox = firstSubsetBox;

      const int subsetsLeft = numSubsets - curSubset - 1;
      const int boxesLeft   = numBoxes - firstSubsetBox - 1;

      if (boxesLeft > subsetsLeft) {
        for (int ibox = firstSubsetBox + 1; ibox < numBoxes; ibox++) {

          // Check if we should add this box - we do this by making sure that the dynamically moving target load stays as close
          // to the static load as possible. We have two scenarios; if we add the current box to the subset then the remaining load becomes
          // smaller but the number of subsets stay the same. If we don't add the current box to the subset the remaining load stays the same
          // but the number of subsets decreases by one.
          const int  remainingSubsets      = numSubsets - curSubset - 1;
          const Real dynamicLoadWithBox    = std::abs(remainingLoad - subsetLoad - a_boxLoads[ibox]) / remainingSubsets;
          const Real dynamicLoadWithoutBox = std::abs(remainingLoad - subsetLoad) / (remainingSubsets - 1);
          const Real loadDiffWithBox       = std::abs(dynamicLoadWithBox - std::abs(staticTargetLoad));
          const Real loadDiffWithoutBox    = std::abs(dynamicLoadWithoutBox - std::abs(staticTargetLoad));
          const bool addBoxToSubset        = std::abs(loadDiffWithBox) < std::abs(loadDiffWithoutBox);

          if (addBoxToSubset) {
            subsetLoad += a_boxLoads[ibox];

            lastSubsetBox = ibox;
          }
          else {
            lastSubsetBox = ibox - 1;

            break;
          }

          // Hook for catching case when we add too many boxes to this subset. Each remaining subset must have at least one box.
          if (numBoxes - lastSubsetBox - 1 <= subsetsLeft) {
            break;
          }
        }
      }

      // Create the subset
      pout() << "create: "
             << "\t" << firstSubsetBox << "\t" << lastSubsetBox << "\t" << subsetLoad << "\t" << staticTargetLoad
             << endl;
      subsets[curSubset] = std::make_pair(std::make_pair(firstSubsetBox, lastSubsetBox), subsetLoad);

      // Update the remaining load and starting box for next iteration
      remainingLoad  = remainingLoad - subsetLoad;
      firstSubsetBox = lastSubsetBox + 1;
      lastSubsetBox  = -1;
    }

#if 1 // Debug hook - remove later
    Real loadSum = 0.0;

    for (const auto& s : subsets) {
      loadSum += s.second;
    }

    pout() << "actual subsets = " << subsets.size() << endl;
    pout() << "actual load = " << loadSum << endl;

    for (int i = 0; i < subsets.size(); i++) {
      pout() << i << "\t" << subsets[i].first.first << "\t" << subsets[i].first.second << "\t" << subsets[i].second
             << endl;
    }
    pout() << "//// end subset report" << endl;
    pout() << endl;

    if (subsets.size() != numSubsets) {
      MayDay::Abort("subset size is wrong");
    }
    if (subsets.front().first.first != 0) {
      MayDay::Abort("first index is wrong");
    }
    if (subsets.back().first.second != numBoxes - 1) {
      MayDay::Abort("second index is wrong");
    }
    for (int i = 0; i < numSubsets - 1; i++) {
      if (subsets[i].first.second + 1 != subsets[i + 1].first.first) {
        MayDay::Abort("indices are wrong");
      }
    }

    if (loadSum != totalLoad) {
      std::cout << loadSum << "\t" << totalLoad << std::endl;
      MayDay::Abort("load sum is wrong");
    }

#endif

    // Sort the subsets from largest to smallest computational load.
    std::sort(subsets.begin(), subsets.end(), [](const Subset& A, const Subset& B) -> bool {
      return A.second >= B.second;
    });

    // Get the accumulated loads per rank, sorted by lowest-to-highest load.
    const std::vector<std::pair<int, Real>>& sortedRankLoads = a_rankLoads.getSortedLoads();

    // Assign the most expensive grid subset to the rank with the lowest accumulated load.
  }
  else {
    a_ranks.resize(0);
  }

#if 1 // Chombo fallback code - to be removed later.
  Vector<long int> proxyLoads(a_boxLoads.size(), 1L);

  LoadBalance(a_ranks, proxyLoads, a_boxes);
#endif
}

void
LoadBalancing::sort(Vector<Box>& a_boxes, const BoxSorting a_which)
{
  CH_TIME("LoadBalancing::sort");

  // The LoadBalancing::sort routines takes pairs of boxes/loads. Just use dummy loads here.
  Vector<int> dummy(a_boxes.size(), 0);

  LoadBalancing::sort(a_boxes, dummy, a_which);
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
